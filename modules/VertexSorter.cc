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
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>

class VertexSorter: public DelphesModule
{
public:
  explicit VertexSorter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fMethod(Steer<std::string>("Method", "BTV")) {}

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "VertexFinder/vertices"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "VertexFinder/tracks"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "clusters"));
    if(const std::string jetInputArray = Steer<std::string>("JetInputArray", ""); !jetInputArray.empty())
      fJetInputArray = ImportArray(jetInputArray);

    // import beamspot
    if(const std::string beamSpotLabel = Steer<std::string>("BeamSpotInputArray", "BeamSpotFilter/beamSpotParticle");
      !beamSpotLabel.empty() && GetFactory()->Has(beamSpotLabel))
      fBeamSpotInputArray = ImportArray(beamSpotLabel);
  }
  void Process() override;

private:
  const std::string fMethod;

  CandidatesCollection fInputArray;
  CandidatesCollection fTrackInputArray;
  CandidatesCollection fJetInputArray;
  CandidatesCollection fBeamSpotInputArray;
  CandidatesCollection fOutputArray;
};

//------------------------------------------------------------------------------

static bool secondDescending(std::pair<unsigned int, double> pair0, std::pair<unsigned int, double> pair1)
{
  return (pair0.second > pair1.second);
}

static bool secondAscending(std::pair<unsigned int, double> pair0, std::pair<unsigned int, double> pair1)
{
  return (pair0.second < pair1.second);
}

//------------------------------------------------------------------------------

void VertexSorter::Process()
{
  fOutputArray->clear();

  Candidate *beamSpotCandidate = nullptr;
  std::map<int, unsigned int> clusterIDToIndex;
  std::map<int, unsigned int>::const_iterator itClusterIDToIndex;
  std::map<int, double> clusterIDToSumPT2;
  std::map<int, double>::const_iterator itClusterIDToSumPT2;
  std::vector<std::pair<int, double> > sortedClusterIDs;
  std::vector<std::pair<int, double> >::const_iterator itSortedClusterIDs;

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
      std::cout << "BTV PV sorting selected, but no jet collection given!" << std::endl;
      throw 0;
    }

    for(Candidate *const &candidate : *fTrackInputArray)
    {
      if(candidate->Momentum.Pt() < 1.0)
        continue;
      if(candidate->ClusterIndex < 0)
        continue;
      TLorentzVector p(candidate->Momentum.Px(), candidate->Momentum.Py(), candidate->Momentum.Pz(), candidate->Momentum.E());
      bool isInJet = false;

      for(Candidate *const &jetCandidate : *fJetInputArray)
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
      sortedClusterIDs.push_back(std::make_pair(itClusterIDToSumPT2->first, itClusterIDToSumPT2->second));
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondDescending);
  }
  else if(fMethod == "GenClosest")
  {
    if(fBeamSpotInputArray->empty())
    {
      std::cout << "Beamspot collection is empty!" << std::endl;
      throw 0;
    }

    beamSpotCandidate = (Candidate *)fBeamSpotInputArray->at(0);
    for(size_t iCluster = 0; iCluster < fInputArray->size(); iCluster++)
    {
      const Candidate &cluster = *((Candidate *)fInputArray->at(iCluster));
      sortedClusterIDs.push_back(std::make_pair(cluster.ClusterIndex, fabs(cluster.Position.Z() - beamSpotCandidate->Position.Z())));
    }
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondAscending);
  }
  else if(fMethod == "GenBest")
  {
    for(Candidate *const &candidate : *fTrackInputArray)
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
      sortedClusterIDs.push_back(std::make_pair(itClusterIDToSumPT2->first, itClusterIDToSumPT2->second));
    sort(sortedClusterIDs.begin(), sortedClusterIDs.end(), secondDescending);
  }
  else
  {
    std::cout << "\"" << fMethod << "\" is not a valid sorting method!" << std::endl;
    std::cout << "Valid methods are:" << std::endl;
    std::cout << "  BTV" << std::endl;
    std::cout << "  GenClosest" << std::endl;
    std::cout << "  GenBest" << std::endl;
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
