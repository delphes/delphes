/** \class VertexSorter
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include <unordered_map>
#include "modules/VertexSorter.h"

//#include "CLHEP/Units/GlobalSystemOfUnits.h"
//#include "CLHEP/Units/GlobalPhysicalConstants.h"

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

static const double mm  = 1.;
static const double m = 1000.*mm;
static const double ns  = 1.;
static const double s = 1.e+9 *ns;
static const double c_light   = 2.99792458e+8 * m/s;

//------------------------------------------------------------------------------

VertexSorter::VertexSorter() :
  fInputArray(NULL), fTrackInputArray(NULL), fItTrackInputArray(NULL), fJetInputArray(NULL), fItJetInputArray(NULL), fOutputArray(NULL)
{
}

//------------------------------------------------------------------------------

VertexSorter::~VertexSorter()
{
}

//------------------------------------------------------------------------------

void VertexSorter::Init()
{
  fInputArray = ImportArray(GetString("InputArray", "VertexFinder/clusters"));

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "VertexFinder/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  if (string (GetString("JetInputArray", "")) != "")
    {
      fJetInputArray = ImportArray(GetString("JetInputArray", ""));
      fItJetInputArray = fJetInputArray->MakeIterator();
    }

  if (string (GetString("BeamSpotInputArray", "")) != "")
    {
      fBeamSpotInputArray = ImportArray(GetString("BeamSpotInputArray", ""));
      fItBeamSpotInputArray = fBeamSpotInputArray->MakeIterator();
    }

  fOutputArray = ExportArray(GetString("OutputArray", "clusters"));

  fMethod = GetString ("Method", "BTV");
}

//------------------------------------------------------------------------------

void VertexSorter::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
}

//------------------------------------------------------------------------------
//
bool VertexSorter::secondDescending (pair<unsigned, double> pair0, pair<unsigned, double> pair1)
{
  return (pair0.second > pair1.second);
}
bool VertexSorter::secondAscending (pair<unsigned, double> pair0, pair<unsigned, double> pair1)
{
  return (pair0.second < pair1.second);
}

void VertexSorter::Process()
{
  Candidate *candidate, *jetCandidate, *beamSpotCandidate;
  unordered_map<int, unsigned> clusterIDToIndex;
  unordered_map<int, double> clusterIDToSumPT2;
  vector<pair<int, double> > sortedClusterIDs;

  for (int iCluster = 0; iCluster < fInputArray->GetEntries (); iCluster++)
    {
      const Candidate &cluster = *((Candidate *) fInputArray->At (iCluster));
      clusterIDToIndex[cluster.ClusterIndex] = iCluster;
      clusterIDToSumPT2[cluster.ClusterIndex] = 0.0;
    }

  if (fMethod == "BTV")
    {
      if (!fJetInputArray)
        {
          cout << "BTV PV sorting selected, but no jet collection given!" << endl;
          throw 0;
        }

      fItTrackInputArray->Reset();
      while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next())))
        {
          if (candidate->Momentum.Pt () < 1.0)
            continue;
          if (candidate->ClusterIndex < 0)
            continue;
          TLorentzVector p (candidate->Momentum.Px (), candidate->Momentum.Py (), candidate->Momentum.Pz (), candidate->Momentum.E ());
          bool isInJet = false;

          fItJetInputArray->Reset();
          while((jetCandidate = static_cast<Candidate*>(fItJetInputArray->Next())))
            {
              if (jetCandidate->Momentum.Pt () < 30.0)
                continue;
              TLorentzVector q (jetCandidate->Momentum.Px (), jetCandidate->Momentum.Py (), jetCandidate->Momentum.Pz (), jetCandidate->Momentum.E ());

              if (p.DeltaR (q) > 0.4)
                continue;
              isInJet = true;
              break;
            }
          if (!isInJet)
            continue;

          clusterIDToSumPT2.at (candidate->ClusterIndex) += candidate->Momentum.Pt () * candidate->Momentum.Pt ();
        }

      for (const auto &clusterID : clusterIDToSumPT2)
        sortedClusterIDs.push_back (make_pair (clusterID.first, clusterID.second));
      sort (sortedClusterIDs.begin (), sortedClusterIDs.end (), secondDescending);
    }
  else if (fMethod == "GenClosest")
    {
      if (!fBeamSpotInputArray)
        {
          cout << "GenClosest PV sorting selected, but no beamspot collection given!" << endl;
          throw 0;
        }
      if (!fBeamSpotInputArray->GetSize ())
        {
          cout << "Beamspot collection is empty!" << endl;
          throw 0;
        }

      beamSpotCandidate = (Candidate *) fBeamSpotInputArray->At (0);
      for (int iCluster = 0; iCluster < fInputArray->GetEntries (); iCluster++)
        {
          const Candidate &cluster = *((Candidate *) fInputArray->At (iCluster));
          sortedClusterIDs.push_back (make_pair (cluster.ClusterIndex, fabs (cluster.Position.Z () - beamSpotCandidate->Position.Z ())));
        }
      sort (sortedClusterIDs.begin (), sortedClusterIDs.end (), secondAscending);
    }
  else if (fMethod == "GenBest")
    {
      fItTrackInputArray->Reset();
      while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next())))
        {
          if (candidate->IsPU)
            continue;
          for (const auto &clusterID : clusterIDToIndex)
            {
              if (candidate->ClusterIndex != clusterID.first)
                continue;
              clusterIDToSumPT2.at (clusterID.first) += candidate->Momentum.Pt () * candidate->Momentum.Pt ();
            }
        }

      for (const auto &clusterID : clusterIDToSumPT2)
        sortedClusterIDs.push_back (make_pair (clusterID.first, clusterID.second));
      sort (sortedClusterIDs.begin (), sortedClusterIDs.end (), secondDescending);
    }
  else
    {
      cout << "\"" << fMethod << "\" is not a valid sorting method!" << endl;
      cout << "Valid methods are:" << endl;
      cout << "  BTV" << endl;
      cout << "  GenClosest" << endl;
      cout << "  GenBest" << endl;
      throw 0;
    }

  for (const auto &clusterID : sortedClusterIDs)
    {
      Candidate *cluster = (Candidate *) fInputArray->At (clusterIDToIndex.at (clusterID.first));
      if (fMethod == "BTV")
        cluster->BTVSumPT2 = clusterID.second;
      else if (fMethod == "GenClosest")
        cluster->GenDeltaZ = clusterID.second;
      else if (fMethod == "GenBest")
        cluster->GenSumPT2 = clusterID.second;
      fOutputArray->Add (cluster);
    }
}

//------------------------------------------------------------------------------
