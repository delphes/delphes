/** \class VertexSorter
 *
 *
 *  Sorts vertices according to different criteria
 *
 *  \authors A. Hart, M. Selvaggi
 *
 *
*/

#include "modules/VertexSorter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

static const Double_t mm = 1.;
static const Double_t m = 1000. * mm;
static const Double_t ns = 1.;
static const Double_t s = 1.e+9 * ns;
static const Double_t c_light = 2.99792458e+8 * m / s;

//------------------------------------------------------------------------------

void VertexSorter::Init()
{
  fMethod = GetString("Method", "BTV");

  // import input arrays
  GetFactory()->EventModel()->Attach(GetString("InputArray", "VertexFinder/vertices"), fInputArray);
  GetFactory()->EventModel()->Attach(GetString("TrackInputArray", "VertexFinder/tracks"), fTrackInputArray);

  if(const auto jet_input_array_label = string(GetString("JetInputArray", "")); !jet_input_array_label.empty())
    GetFactory()->EventModel()->Attach(jet_input_array_label, fJetInputArray);

  try
  { // import beamspot
    GetFactory()->EventModel()->Attach(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"), fBeamSpotInputArray);
  }
  catch(runtime_error &e)
  {
  }

  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "clusters"));
}

//------------------------------------------------------------------------------

void VertexSorter::Finish()
{
}

//------------------------------------------------------------------------------

static Bool_t secondDescending(pair<UInt_t, Double_t> pair0, pair<UInt_t, Double_t> pair1)
{
  return (pair0.second > pair1.second);
}

static Bool_t secondAscending(pair<UInt_t, Double_t> pair0, pair<UInt_t, Double_t> pair1)
{
  return (pair0.second < pair1.second);
}

//------------------------------------------------------------------------------

void VertexSorter::Process()
{
  map<Int_t, UInt_t> clusterIDToIndex;
  map<Int_t, UInt_t>::const_iterator itClusterIDToIndex;
  map<Int_t, Double_t> clusterIDToSumPT2;
  map<Int_t, Double_t>::const_iterator itClusterIDToSumPT2;
  vector<pair<Int_t, Double_t> > sortedClusterIDs;
  vector<pair<Int_t, Double_t> >::const_iterator itSortedClusterIDs;

  for(size_t iCluster = 0; iCluster < fInputArray->size(); iCluster++)
  {
    const auto &cluster = fInputArray->at(iCluster);
    clusterIDToIndex[cluster.ClusterIndex] = iCluster;
    clusterIDToSumPT2[cluster.ClusterIndex] = 0.0;
  }

  if(fMethod == "BTV")
  {
    if(!fJetInputArray)
    {
      cout << "BTV PV sorting selected, but no jet collection given!" << endl;
      throw 0;
    }

    for(const auto &candidate : *fTrackInputArray)
    {
      if(candidate.Momentum.Pt() < 1.0)
        continue;
      if(candidate.ClusterIndex < 0)
        continue;
      TLorentzVector p(candidate.Momentum.Px(), candidate.Momentum.Py(), candidate.Momentum.Pz(), candidate.Momentum.E());
      Bool_t isInJet = false;

      for(const auto &jetCandidate : *fJetInputArray)
      {
        if(jetCandidate.Momentum.Pt() < 30.0)
          continue;
        TLorentzVector q(jetCandidate.Momentum.Px(), jetCandidate.Momentum.Py(), jetCandidate.Momentum.Pz(), jetCandidate.Momentum.E());

        if(p.DeltaR(q) > 0.4)
          continue;
        isInJet = true;
        break;
      }
      if(!isInJet)
        continue;

      clusterIDToSumPT2.at(candidate.ClusterIndex) += candidate.Momentum.Pt() * candidate.Momentum.Pt();
    }

    for(itClusterIDToSumPT2 = clusterIDToSumPT2.begin(); itClusterIDToSumPT2 != clusterIDToSumPT2.end(); ++itClusterIDToSumPT2)
      sortedClusterIDs.push_back(make_pair(itClusterIDToSumPT2->first, itClusterIDToSumPT2->second));
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondDescending);
  }
  else if(fMethod == "GenClosest")
  {
    if(!fBeamSpotInputArray)
    {
      cout << "GenClosest PV sorting selected, but no beamspot collection given!" << endl;
      throw 0;
    }
    if(fBeamSpotInputArray->empty())
    {
      cout << "Beamspot collection is empty!" << endl;
      throw 0;
    }

    const auto &beamSpotCandidate = fBeamSpotInputArray->at(0);
    for(size_t iCluster = 0; iCluster < fInputArray->size(); iCluster++)
    {
      const auto &cluster = fInputArray->at(iCluster);
      sortedClusterIDs.push_back(make_pair(cluster.ClusterIndex, fabs(cluster.Position.Z() - beamSpotCandidate.Position.Z())));
    }
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondAscending);
  }
  else if(fMethod == "GenBest")
  {
    for(const auto &candidate : *fTrackInputArray)
    {
      if(candidate.IsPU)
        continue;
      for(itClusterIDToIndex = clusterIDToIndex.begin(); itClusterIDToIndex != clusterIDToIndex.end(); ++itClusterIDToIndex)
      {
        if(candidate.ClusterIndex != itClusterIDToIndex->first)
          continue;
        clusterIDToSumPT2.at(itClusterIDToIndex->first) += candidate.Momentum.Pt() * candidate.Momentum.Pt();
      }
    }

    for(itClusterIDToSumPT2 = clusterIDToSumPT2.begin(); itClusterIDToSumPT2 != clusterIDToSumPT2.end(); ++itClusterIDToSumPT2)
      sortedClusterIDs.push_back(make_pair(itClusterIDToSumPT2->first, itClusterIDToSumPT2->second));
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondDescending);
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
  for(itSortedClusterIDs = sortedClusterIDs.begin(); itSortedClusterIDs != sortedClusterIDs.end(); ++itSortedClusterIDs)
  {
    auto &cluster = fInputArray->at(clusterIDToIndex.at(itSortedClusterIDs->first)); //TODO: do we want modification of input cluster?
    if(fMethod == "BTV")
      cluster.BTVSumPT2 = itSortedClusterIDs->second;
    else if(fMethod == "GenClosest")
      cluster.GenDeltaZ = itSortedClusterIDs->second;
    else if(fMethod == "GenBest")
      cluster.GenSumPT2 = itSortedClusterIDs->second;
    fOutputArray->emplace_back(cluster);
  }
}

//------------------------------------------------------------------------------
