/** \class VertexFinder
 *
 *  Cluster vertices from tracks
 *
 *  \authors A. Hart, M. Selvaggi
 *
 */


#include "modules/VertexFinder.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMatrixT.h"
#include "TVector3.h"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <map>
#include <string>

using namespace std;

static const Double_t mm  = 1.;
static const Double_t m = 1000.*mm;
static const Double_t ns  = 1.;
static const Double_t s = 1.e+9 *ns;
static const Double_t c_light   = 2.99792458e+8 * m/s;

//------------------------------------------------------------------------------

VertexFinder::VertexFinder() :
  fSigma(0), fMinPT(0), fMaxEta(0), fSeedMinPT(0), fMinNDF(0), fGrowSeeds(0)
{
}

//------------------------------------------------------------------------------

VertexFinder::~VertexFinder()
{
}

//------------------------------------------------------------------------------

void VertexFinder::Init()
{
  fSigma = GetDouble("Sigma", 3.0);
  fMinPT = GetDouble("MinPT", 0.1);
  fMaxEta = GetDouble("MaxEta", 10.0);
  fSeedMinPT = GetDouble("SeedMinPT", 5.0);
  fMinNDF = GetInt("MinNDF", 4);
  fGrowSeeds = GetInt("GrowSeeds", 1);

  fInputArray = ImportArray(GetString("InputArray", "TrackSmearing/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
  fVertexOutputArray = ExportArray(GetString("VertexOutputArray", "vertices"));
}

//------------------------------------------------------------------------------

void VertexFinder::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

static Bool_t secondAscending (pair<UInt_t, Double_t> pair0, pair<UInt_t, Double_t> pair1)
{
  return (pair0.second < pair1.second);
}

static Bool_t secondDescending (pair<UInt_t, Double_t> pair0, pair<UInt_t, Double_t> pair1)
{
  return (pair0.second > pair1.second);
}

//------------------------------------------------------------------------------

void VertexFinder::Process()
{
  Candidate *candidate;

  // Clear the track and cluster maps before starting
  trackIDToDouble.clear ();
  trackIDToInt.clear ();
  trackIDToBool.clear ();
  clusterIDToDouble.clear ();
  clusterIDToInt.clear ();
  clusterIDToBool.clear ();
  trackPT.clear ();
  clusterSumPT2.clear ();

  // Create the initial cluster seeds
  createSeeds ();

  // In order of descending seed pt, grow each cluster. If a cluster ends up with
  // fewer than MinNDF tracks, release the tracks for other clusters to claim.
  sort (clusterSumPT2.begin (), clusterSumPT2.end (), secondDescending);
  for (vector<pair<UInt_t, Double_t> >::const_iterator cluster = clusterSumPT2.begin (); cluster != clusterSumPT2.end (); cluster++)
    {
      // Skip the cluster if it no longer has any tracks
      if (!clusterIDToInt.at (cluster->first).at ("ndf"))
        continue;
        
      // Grow the cluster if GrowSeeds is true
      if (fGrowSeeds)
        growCluster (cluster->first);

      // If the cluster still has fewer than MinNDF tracks, release the tracks;
      // otherwise, mark the seed track as claimed
      
      if ((Int_t) clusterIDToInt.at (cluster->first).at ("ndf") < fMinNDF)
        {
          for (map<UInt_t, map<string, Int_t> >::iterator track = trackIDToInt.begin (); track != trackIDToInt.end (); track++)
            {
              if (track->second.at ("clusterIndex") != (Int_t) cluster->first)
                continue;
              track->second["clusterIndex"] = -1;
              trackIDToBool[track->first]["claimed"] = false;
            }
        }
      else
        trackIDToBool[clusterIDToInt.at (cluster->first).at ("seed")]["claimed"] = true;
    }

  // Add tracks to the output array after updating their ClusterIndex.
  fItInputArray->Reset ();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
    {
      if (candidate->Momentum.Pt () < fMinPT || fabs (candidate->Momentum.Eta ()) > fMaxEta)
        continue;
      candidate->ClusterIndex = trackIDToInt.at (candidate->GetUniqueID ()).at ("clusterIndex");
      fOutputArray->Add(candidate);
    }

  // Add clusters with at least MinNDF tracks to the output array in order of
  // descending sum(pt**2).
  clusterSumPT2.clear ();
  for (map<UInt_t, map<string, Int_t> >::const_iterator cluster = clusterIDToInt.begin (); cluster != clusterIDToInt.end (); cluster++)
  {
    
    if (cluster->second.at ("ndf") < fMinNDF)
      continue;
    clusterSumPT2.push_back (make_pair (cluster->first, clusterIDToDouble.at (cluster->first).at ("sumPT2")));
  }
  sort (clusterSumPT2.begin (), clusterSumPT2.end (), secondDescending);
  
  for (vector<pair<UInt_t, Double_t> >::const_iterator cluster = clusterSumPT2.begin (); cluster != clusterSumPT2.end (); cluster++)
  {
    DelphesFactory *factory = GetFactory();
    candidate = factory->NewCandidate();

    candidate->ClusterIndex = cluster->first;
    candidate->ClusterNDF = clusterIDToInt.at (cluster->first).at ("ndf");
    candidate->ClusterSigma = fSigma;
    candidate->SumPT2 = cluster->second;
    candidate->Position.SetXYZT(0.0, 0.0, clusterIDToDouble.at (cluster->first).at ("z"), 0.0);
    candidate->PositionError.SetXYZT(0.0, 0.0, clusterIDToDouble.at (cluster->first).at ("ez"), 0.0);

    fVertexOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

void VertexFinder::createSeeds ()
{
  Candidate *candidate;
  UInt_t clusterIndex = 0, maxSeeds = 0;

  // Loop over all tracks, initializing some variables.
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
    {
      if (candidate->Momentum.Pt () < fMinPT || fabs (candidate->Momentum.Eta ()) > fMaxEta)
        continue;

      trackIDToDouble[candidate->GetUniqueID ()]["pt"] = candidate->Momentum.Pt ();
      trackIDToDouble[candidate->GetUniqueID ()]["ept"] = candidate->ErrorPT ? candidate->ErrorPT : 1.0e-15;;
      trackIDToDouble[candidate->GetUniqueID ()]["eta"] = candidate->Momentum.Eta ();

      trackIDToDouble[candidate->GetUniqueID ()]["z"] = candidate->DZ;
      trackIDToDouble[candidate->GetUniqueID ()]["ez"] = candidate->ErrorDZ ? candidate->ErrorDZ : 1.0e-15;

      trackIDToInt[candidate->GetUniqueID ()]["clusterIndex"] = -1;
      trackIDToInt[candidate->GetUniqueID ()]["interactionIndex"] = candidate->IsPU;

      trackIDToBool[candidate->GetUniqueID ()]["claimed"] = false;

      trackPT.push_back (make_pair (candidate->GetUniqueID (), candidate->Momentum.Pt ()));
    }

  // Sort tracks by pt and leave only the SeedMinPT highest pt ones in the
  // trackPT vector.
  sort (trackPT.begin (), trackPT.end (), secondDescending);
  for (vector<pair<UInt_t, Double_t> >::const_iterator track = trackPT.begin (); track != trackPT.end (); track++, maxSeeds++)
    {
      if (track->second < fSeedMinPT)
        break;
    }
  // If there are no tracks with pt above MinSeedPT, create just one seed from
  // the highest pt track.
  if (!maxSeeds)
    maxSeeds++;
  if (trackPT.size () > maxSeeds)
    {
      trackPT.erase (trackPT.begin () + maxSeeds, trackPT.end ());
    }

  // Create the seeds from the SeedMinPT highest pt tracks.
  for (vector<pair<UInt_t, Double_t> >::const_iterator track = trackPT.begin (); track != trackPT.end (); track++, clusterIndex++)
    {
      addTrackToCluster (track->first, clusterIndex);
      clusterSumPT2.push_back (make_pair (clusterIndex, track->second * track->second));
    }
}

//------------------------------------------------------------------------------

void VertexFinder::growCluster (const UInt_t clusterIndex)
{
  Bool_t done = false;
  UInt_t nearestID;
  Int_t oldClusterIndex;
  Double_t nearestDistance;
  vector<UInt_t> nearTracks;
  nearTracks.clear ();

  // Grow the cluster until there are no more tracks within Sigma standard
  // deviations of the cluster.
  while (!done)
    {
      done = true;
      nearestID = 0;
      nearestDistance = -1.0;

      // These two loops are for finding the nearest track to the cluster. The
      // first time, the ID of each track within 10*Sigma of the cluster is
      // saved in the nearTracks vector; subsequently, to save time, only the
      // tracks in this vector are checked.
      if (!nearTracks.size ())
        {
        
            for (map<UInt_t, map<string, Double_t> >::const_iterator track = trackIDToDouble.begin (); track != trackIDToDouble.end (); track++)
            {
              if (trackIDToBool.at (track->first).at ("claimed") || trackIDToInt.at (track->first).at ("clusterIndex") == (Int_t) clusterIndex)
                continue;
                
              Double_t distance = fabs (clusterIDToDouble.at (clusterIndex).at ("z") - track->second.at ("z")) / hypot (clusterIDToDouble.at (clusterIndex).at ("ez"), track->second.at ("ez"));
              if (nearestDistance < 0.0 || distance < nearestDistance)
                {
                  nearestID = track->first;
                  nearestDistance = distance;
                }
              if (distance < 10.0 * fSigma)
                nearTracks.push_back (track->first);
            }
        }
      
      else
        {
          for (vector<UInt_t>::const_iterator track = nearTracks.begin (); track != nearTracks.end (); track++)
            {
              if (trackIDToBool.at (*track).at ("claimed") || trackIDToInt.at (*track).at ("clusterIndex") == (Int_t) clusterIndex)
                continue;
              Double_t distance = fabs (clusterIDToDouble.at (clusterIndex).at ("z") - trackIDToDouble.at (*track).at ("z")) / hypot (clusterIDToDouble.at (clusterIndex).at ("ez"), trackIDToDouble.at (*track).at ("ez"));
              if (nearestDistance < 0.0 || distance < nearestDistance)
                {
                  nearestID = *track;
                  nearestDistance = distance;
                }
            }
        }
      
      // If no tracks within Sigma of the cluster were found, stop growing.
      done = nearestDistance > fSigma || nearestDistance < 0.0;
      if (done)
        {
          continue;
        }

      // Add the nearest track within Sigma to the cluster. If it already
      // belonged to another cluster, remove it from that cluster first.
      if (nearestDistance < fSigma)
        {
          oldClusterIndex = trackIDToInt.at (nearestID).at ("clusterIndex");
          if (oldClusterIndex >= 0)
            removeTrackFromCluster (nearestID, oldClusterIndex);

          trackIDToBool[nearestID]["claimed"] = true;
          addTrackToCluster (nearestID, clusterIndex);
        }
    }
}

//------------------------------------------------------------------------------

Double_t VertexFinder::weight (const UInt_t trackID)
{
  return ((trackIDToDouble.at (trackID).at ("pt") / (trackIDToDouble.at (trackID).at ("ept") * trackIDToDouble.at (trackID).at ("ez"))) * (trackIDToDouble.at (trackID).at ("pt") / (trackIDToDouble.at (trackID).at ("ept") * trackIDToDouble.at (trackID).at ("ez"))));
}

//------------------------------------------------------------------------------

void VertexFinder::removeTrackFromCluster (const UInt_t trackID, const UInt_t clusterID)
{
  Double_t wz = weight (trackID);

  trackIDToInt[trackID]["clusterIndex"] = -1;
  clusterIDToInt[clusterID]["ndf"]--;

  clusterIDToDouble[clusterID]["sumZ"] -= wz * trackIDToDouble.at (trackID).at ("z");
  clusterIDToDouble[clusterID]["errorSumZ"] -= wz * trackIDToDouble.at (trackID).at ("ez") * trackIDToDouble.at (trackID).at ("ez");
  clusterIDToDouble[clusterID]["sumOfWeightsZ"] -= wz;
  clusterIDToDouble[clusterID]["z"] = clusterIDToDouble.at (clusterID).at ("sumZ") / clusterIDToDouble.at (clusterID).at ("sumOfWeightsZ");
  clusterIDToDouble[clusterID]["ez"] = sqrt ((1.0 / clusterIDToInt.at (clusterID).at ("ndf")) * (clusterIDToDouble.at (clusterID).at ("errorSumZ") / clusterIDToDouble.at (clusterID).at ("sumOfWeightsZ")));
  clusterIDToDouble[clusterID]["sumPT2"] -= trackIDToDouble.at (trackID).at ("pt") * trackIDToDouble.at (trackID).at ("pt");
}

//------------------------------------------------------------------------------

void VertexFinder::addTrackToCluster (const UInt_t trackID, const UInt_t clusterID)
{
  Double_t wz = weight (trackID);

  if (!clusterIDToInt.count (clusterID))
    {
      clusterIDToInt[clusterID]["ndf"] = 0;
      clusterIDToInt[clusterID]["seed"] = trackID;
      clusterIDToDouble[clusterID]["sumZ"] = 0.0;
      clusterIDToDouble[clusterID]["errorSumZ"] = 0.0;
      clusterIDToDouble[clusterID]["sumOfWeightsZ"] = 0.0;
      clusterIDToDouble[clusterID]["sumPT2"] = 0.0;
    }

  trackIDToInt[trackID]["clusterIndex"] = clusterID;
  clusterIDToInt[clusterID]["ndf"]++;

  clusterIDToDouble[clusterID]["sumZ"] += wz * trackIDToDouble.at (trackID).at ("z");
  clusterIDToDouble[clusterID]["errorSumZ"] += wz * trackIDToDouble.at (trackID).at ("ez") * trackIDToDouble.at (trackID).at ("ez");
  clusterIDToDouble[clusterID]["sumOfWeightsZ"] += wz;
  clusterIDToDouble[clusterID]["z"] = clusterIDToDouble.at (clusterID).at ("sumZ") / clusterIDToDouble.at (clusterID).at ("sumOfWeightsZ");
  clusterIDToDouble[clusterID]["ez"] = sqrt ((1.0 / clusterIDToInt.at (clusterID).at ("ndf")) * (clusterIDToDouble.at (clusterID).at ("errorSumZ") / clusterIDToDouble.at (clusterID).at ("sumOfWeightsZ")));
  clusterIDToDouble[clusterID]["sumPT2"] += trackIDToDouble.at (trackID).at ("pt") * trackIDToDouble.at (trackID).at ("pt");
}

//------------------------------------------------------------------------------
