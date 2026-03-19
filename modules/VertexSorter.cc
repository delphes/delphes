/** \class VertexSorter
 *
 *
 *  Sorts vertices according to different criteria
 *
 *  \authors A. Hart, M. Selvaggi
 *
 *
*/

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

using namespace std;

class VertexSorter: public DelphesModule
{
public:
  VertexSorter() = default;

  void Init() override;
  void Process() override;

private:
  CandidatesCollection fInputArray;
  CandidatesCollection fTrackInputArray;
  CandidatesCollection fJetInputArray;
  CandidatesCollection fBeamSpotInputArray;

  CandidatesCollection fOutputArray;

  std::string fMethod;
};

//------------------------------------------------------------------------------

void VertexSorter::Init()
{
  fInputArray = ImportArray(GetString("InputArray", "VertexFinder/vertices"));
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "VertexFinder/tracks"));
  if(string(GetString("JetInputArray", "")) != "")
    fJetInputArray = ImportArray(GetString("JetInputArray", ""));

  // import beamspot
  try
  {
    fBeamSpotInputArray = ImportArray(GetString("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle"));
  }
  catch(runtime_error &e)
  {
  }

  fOutputArray = ExportArray(GetString("OutputArray", "clusters"));

  fMethod = GetString("Method", "BTV");
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
  fOutputArray->clear();

  Candidate *beamSpotCandidate = nullptr;
  map<Int_t, UInt_t> clusterIDToIndex;
  map<Int_t, UInt_t>::const_iterator itClusterIDToIndex;
  map<Int_t, Double_t> clusterIDToSumPT2;
  map<Int_t, Double_t>::const_iterator itClusterIDToSumPT2;
  vector<pair<Int_t, Double_t> > sortedClusterIDs;
  vector<pair<Int_t, Double_t> >::const_iterator itSortedClusterIDs;

  for(size_t iCluster = 0; iCluster < fInputArray->size(); iCluster++)
  {
    const Candidate &cluster = *((Candidate *)fInputArray->at(iCluster));
    clusterIDToIndex[cluster.ClusterIndex] = iCluster;
    clusterIDToSumPT2[cluster.ClusterIndex] = 0.0;
  }

  if(fMethod == "BTV")
  {
    if(fJetInputArray->empty())
    {
      cout << "BTV PV sorting selected, but no jet collection given!" << endl;
      throw 0;
    }

    for(const auto &candidate : *fTrackInputArray)
    {
      if(candidate->Momentum.Pt() < 1.0)
        continue;
      if(candidate->ClusterIndex < 0)
        continue;
      TLorentzVector p(candidate->Momentum.Px(), candidate->Momentum.Py(), candidate->Momentum.Pz(), candidate->Momentum.E());
      Bool_t isInJet = false;

      for(const auto &jetCandidate : *fJetInputArray)
      {
        if(jetCandidate->Momentum.Pt() < 30.0)
          continue;
        TLorentzVector q(jetCandidate->Momentum.Px(), jetCandidate->Momentum.Py(), jetCandidate->Momentum.Pz(), jetCandidate->Momentum.E());

        if(p.DeltaR(q) > 0.4)
          continue;
        isInJet = true;
        break;
      }
      if(!isInJet)
        continue;

      clusterIDToSumPT2.at(candidate->ClusterIndex) += candidate->Momentum.Pt() * candidate->Momentum.Pt();
    }

    for(itClusterIDToSumPT2 = clusterIDToSumPT2.begin(); itClusterIDToSumPT2 != clusterIDToSumPT2.end(); ++itClusterIDToSumPT2)
      sortedClusterIDs.push_back(make_pair(itClusterIDToSumPT2->first, itClusterIDToSumPT2->second));
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondDescending);
  }
  else if(fMethod == "GenClosest")
  {
    if(fBeamSpotInputArray->empty())
    {
      cout << "Beamspot collection is empty!" << endl;
      throw 0;
    }

    beamSpotCandidate = (Candidate *)fBeamSpotInputArray->at(0);
    for(size_t iCluster = 0; iCluster < fInputArray->size(); iCluster++)
    {
      const Candidate &cluster = *((Candidate *)fInputArray->at(iCluster));
      sortedClusterIDs.push_back(make_pair(cluster.ClusterIndex, fabs(cluster.Position.Z() - beamSpotCandidate->Position.Z())));
    }
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondAscending);
  }
  else if(fMethod == "GenBest")
  {
    for(const auto &candidate : *fTrackInputArray)
    {
      if(candidate->IsPU)
        continue;
      for(itClusterIDToIndex = clusterIDToIndex.begin(); itClusterIDToIndex != clusterIDToIndex.end(); ++itClusterIDToIndex)
      {
        if(candidate->ClusterIndex != itClusterIDToIndex->first)
          continue;
        clusterIDToSumPT2.at(itClusterIDToIndex->first) += candidate->Momentum.Pt() * candidate->Momentum.Pt();
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
    Candidate *cluster = (Candidate *)fInputArray->at(clusterIDToIndex.at(itSortedClusterIDs->first));
    if(fMethod == "BTV")
      cluster->BTVSumPT2 = itSortedClusterIDs->second;
    else if(fMethod == "GenClosest")
      cluster->GenDeltaZ = itSortedClusterIDs->second;
    else if(fMethod == "GenBest")
      cluster->GenSumPT2 = itSortedClusterIDs->second;
    fOutputArray->emplace_back(cluster);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("VertexSorter", VertexSorter);
