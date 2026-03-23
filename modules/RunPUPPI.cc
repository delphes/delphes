#include <PUPPI/AlgoObj.hh>
#include <PUPPI/PuppiContainer.hh>
#include <PUPPI/RecoObj2.hh>

#include <fastjet/PseudoJet.hh>

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace fastjet;

class RunPUPPI: public DelphesModule
{
public:
  explicit RunPUPPI(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    // puppi parameters
    fApplyNoLep(Steer<bool>("UseNoLep", true)),
    fEtaMinBin(Steer<std::vector<float> >("EtaMinBin")), // read eta min ranges
    fEtaMaxBin(Steer<std::vector<float> >("EtaMaxBin")), // read eta max ranges
    fPtMinBin(Steer<std::vector<float> >("PtMinBin")), // read pt min value
    fConeSizeBin(Steer<std::vector<float> >("ConeSizeBin")), // read cone size
    fRMSPtMinBin(Steer<std::vector<float> >("RMSPtMinBin")), // read RMS min pt
    fRMSScaleFactorBin(Steer<std::vector<float> >("RMSScaleFactorBin")), // read RMS scale factor
    fNeutralMinEBin(Steer<std::vector<float> >("NeutralMinEBin")), // read neutral pt min cut
    fNeutralPtSlope(Steer<std::vector<float> >("NeutralPtSlope")), // read neutral pt min slope
    //fApplyCHS(Steer<std::vector<bool> >("ApplyCHS")), // read apply chs
    fUseCharged(Steer<std::vector<bool> >("UseCharged")), // read use charged
    fApplyLowPUCorr(Steer<std::vector<bool> >("ApplyLowPUCorr")), // read apply chs correction
    fMetricId(Steer<std::vector<int> >("MetricId")), // read metric id
    fCombId(Steer<std::vector<int> >("CombId")) // scheme for combining
  {
    // create output array
    // Create algorithm list for puppi
    std::vector<AlgoObj> puppiAlgo;
    if(puppiAlgo.empty())
    {
      if(!(fEtaMinBin.size() == fEtaMaxBin.size() and fEtaMinBin.size() == fPtMinBin.size() and fEtaMinBin.size() == fConeSizeBin.size() and fEtaMinBin.size() == fRMSPtMinBin.size()
           and fEtaMinBin.size() == fRMSScaleFactorBin.size() and fEtaMinBin.size() == fNeutralMinEBin.size() and fEtaMinBin.size() == fNeutralPtSlope.size()
           and fEtaMinBin.size() == fUseCharged.size()
           and fEtaMinBin.size() == fApplyLowPUCorr.size() and fEtaMinBin.size() == fMetricId.size()))
      {
        std::cerr << " Error in PUPPI configuration, algo info should have the same size --> exit from the code" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    for(size_t iAlgo = 0; iAlgo < fEtaMinBin.size(); iAlgo++)
    {
      AlgoObj algoTmp;
      algoTmp.etaMin = fEtaMinBin.at(iAlgo);
      algoTmp.etaMax = fEtaMaxBin.at(iAlgo);
      algoTmp.ptMin = fPtMinBin.at(iAlgo);
      algoTmp.minNeutralPt = fNeutralMinEBin.at(iAlgo);
      algoTmp.minNeutralPtSlope = fNeutralPtSlope.at(iAlgo);
      //Eta Extrapolation stuff is missing
      //Loop through file requiring algos for same bins to be adjacent
      while(iAlgo < fEtaMinBin.size() and algoTmp.etaMin == fEtaMinBin.at(iAlgo) and algoTmp.etaMax == fEtaMaxBin.at(iAlgo))
      {
        AlgoSubObj algoSubTmp;
        algoSubTmp.metricId = fMetricId.at(iAlgo);
        algoSubTmp.useCharged = fUseCharged.at(iAlgo);
        algoSubTmp.applyLowPUCorr = fApplyLowPUCorr.at(iAlgo);
        algoSubTmp.combId = fCombId.at(iAlgo);
        algoSubTmp.coneSize = fConeSizeBin.at(iAlgo);
        algoSubTmp.rmsPtMin = fRMSPtMinBin.at(iAlgo);
        algoSubTmp.rmsScaleFactor = fRMSScaleFactorBin.at(iAlgo);
        algoTmp.subAlgos.push_back(algoSubTmp);
        iAlgo++;
      }
      iAlgo--;
      //if(std::find(puppiAlgo.begin(),puppiAlgo.end(),algoTmp) != puppiAlgo.end()) continue;
      puppiAlgo.push_back(algoTmp);
    }
    fPuppi = std::make_unique<PuppiContainer>(true,
      Steer<bool>("UseExp", false),
      Steer<double>("MinPuppiWeight", 0.01),
      puppiAlgo);
  }

  void Init() override
  {
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "Calorimeter/towers"));
    fNeutralInputArray = ImportArray(Steer<std::string>("NeutralInputArray", "Calorimeter/towers"));
    fPVInputArray = ImportArray(Steer<std::string>("PVInputArray", "PV"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "puppiParticles"));
    fOutputTrackArray = ExportArray(Steer<std::string>("OutputArrayTracks", "puppiTracks"));
    fOutputNeutralArray = ExportArray(Steer<std::string>("OutputArrayNeutrals", "puppiNeutrals"));
  }
  void Process() override;

private:
  // puppi parameters
  const bool fApplyNoLep;

  const std::vector<float> fEtaMinBin;
  const std::vector<float> fEtaMaxBin;
  const std::vector<float> fPtMinBin;
  const std::vector<float> fConeSizeBin;
  const std::vector<float> fRMSPtMinBin;
  const std::vector<float> fRMSScaleFactorBin;
  const std::vector<float> fNeutralMinEBin;
  const std::vector<float> fNeutralPtSlope;
  //const std::vector<bool> fApplyCHS;
  const std::vector<bool> fUseCharged;
  const std::vector<bool> fApplyLowPUCorr;
  const std::vector<int> fMetricId;
  const std::vector<int> fCombId;

  std::unique_ptr<PuppiContainer> fPuppi;

  CandidatesCollection fTrackInputArray;
  CandidatesCollection fNeutralInputArray; //!
  CandidatesCollection fPVInputArray; //!

  CandidatesCollection fOutputArray;
  CandidatesCollection fOutputTrackArray;
  CandidatesCollection fOutputNeutralArray;
};

//------------------------------------------------------------------------------

void RunPUPPI::Process()
{
  fOutputArray->clear();
  fOutputTrackArray->clear();
  fOutputNeutralArray->clear();

  Candidate *candidate, *particle;
  TLorentzVector momentum;

  //DelphesFactory *factory = GetFactory();

  // loop over input objects
  std::vector<Candidate *> InputParticles;
  InputParticles.clear();

  // take the leading vertex
  float PVZ = 0.;
  Candidate *pv = static_cast<Candidate *>(fPVInputArray->at(0));
  if(pv) PVZ = pv->Position.Z();
  // Fill input particles for puppi
  std::vector<RecoObj> puppiInputVector;
  puppiInputVector.clear();
  // Loop on charge track candidate
  for(Candidate *const &candidate : *fTrackInputArray)
  {
    momentum = candidate->Momentum;
    RecoObj curRecoObj;
    curRecoObj.pt = momentum.Pt();
    curRecoObj.eta = momentum.Eta();
    curRecoObj.phi = momentum.Phi();
    curRecoObj.m = momentum.M();
    particle = static_cast<Candidate *>(candidate->GetCandidates().at(0)); //if(fApplyNoLep && std::abs(candidate->PID) == 11) continue; //Dumb cut to minimize the nolepton on electron
    //if(fApplyNoLep && std::abs(candidate->PID) == 13) continue;
    if(candidate->IsRecoPU and candidate->Charge != 0)
    { // if it comes fromPU vertexes after the resolution smearing and the dZ matching within resolution
      curRecoObj.id = 2;
      curRecoObj.vtxId = 0.7 * (fPVInputArray->size()); //Hack apply reco vtx efficiency of 70% for calibration
      if(std::abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(std::abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(std::abs(candidate->PID) == 22)
        curRecoObj.pfType = 4;
      else
        curRecoObj.pfType = 1;
      curRecoObj.dZ = particle->Position.Z() - PVZ;
    }
    else if(!candidate->IsRecoPU && candidate->Charge != 0)
    {
      curRecoObj.id = 1; // charge from LV
      curRecoObj.vtxId = 1; // from PV
      if(std::abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(std::abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(std::abs(candidate->PID) == 22)
        curRecoObj.pfType = 4;
      else
        curRecoObj.pfType = 1;
      curRecoObj.dZ = particle->Position.Z() - PVZ;
    }
    else
    {
      std::cerr << " RunPUPPI: problem with a charged track --> it has charge 0 " << std::endl;
      continue;
    }

    puppiInputVector.push_back(curRecoObj);
    InputParticles.push_back(candidate);
  }

  // Loop on neutral calo cells
  for(Candidate *const &candidate : *fNeutralInputArray)
  {
    momentum = candidate->Momentum;
    RecoObj curRecoObj;
    curRecoObj.pt = momentum.Pt();
    curRecoObj.eta = momentum.Eta();
    curRecoObj.phi = momentum.Phi();
    curRecoObj.m = momentum.M();
    curRecoObj.charge = 0;
    particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));
    if(candidate->Charge == 0)
    {
      curRecoObj.id = 0; // neutrals have id==0
      curRecoObj.vtxId = 0; // neutrals have vtxId==0
      if(std::abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(std::abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(std::abs(candidate->PID) == 22)
        curRecoObj.pfType = 4;
      else
        curRecoObj.pfType = 5;
      curRecoObj.dZ = particle->Position.Z() - PVZ;
    }
    else
    {
      std::cerr << " RunPUPPI: problem with a neutrals cells --> it has charge !=0 " << std::endl;
      continue;
    }
    puppiInputVector.push_back(curRecoObj);
    InputParticles.push_back(candidate);
  }
  // Create PUPPI container
  fPuppi->initialize(puppiInputVector);
  fPuppi->puppiWeights();
  std::vector<PseudoJet> puppiParticles = fPuppi->puppiParticles();

  // Loop on final particles
  for(std::vector<PseudoJet>::iterator it = puppiParticles.begin(); it != puppiParticles.end(); it++)
  {
    if(it->user_index() <= int(InputParticles.size()))
    {
      candidate = static_cast<Candidate *>(InputParticles.at(it->user_index())->Clone());
      candidate->Momentum.SetPxPyPzE(it->px(), it->py(), it->pz(), it->e());
      fOutputArray->emplace_back(candidate);
      if(puppiInputVector.at(it->user_index()).id == 1 or puppiInputVector.at(it->user_index()).id == 2)
        fOutputTrackArray->emplace_back(candidate);
      else if(puppiInputVector.at(it->user_index()).id == 0)
        fOutputNeutralArray->emplace_back(candidate);
    }
    else
    {
      std::cerr << " particle not found in the input Array --> skip " << std::endl;
      continue;
    }
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("RunPUPPI", RunPUPPI);
