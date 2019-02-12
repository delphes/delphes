#include "modules/RunPUPPI.h"

#include "PUPPI/AlgoObj.hh"
#include "PUPPI/PuppiContainer.hh"
#include "PUPPI/RecoObj2.hh"

#include "fastjet/PseudoJet.hh"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace fastjet;

//------------------------------------------------------------------------------
RunPUPPI::RunPUPPI() :
  fItTrackInputArray(0),
  fItNeutralInputArray(0)
{
}

//------------------------------------------------------------------------------
RunPUPPI::~RunPUPPI() {}

//------------------------------------------------------------------------------

void RunPUPPI::Init()
{
  // input collection
  fTrackInputArray = ImportArray(GetString("TrackInputArray", "Calorimeter/towers"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "Calorimeter/towers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();
  fPVInputArray = ImportArray(GetString("PVInputArray", "PV"));
  fPVItInputArray = fPVInputArray->MakeIterator();
  // puppi parameters
  fApplyNoLep = GetBool("UseNoLep", true);
  fMinPuppiWeight = GetDouble("MinPuppiWeight", 0.01);
  fUseExp = GetBool("UseExp", false);
  // read eta min ranges
  ExRootConfParam param = GetParam("EtaMinBin");
  fEtaMinBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fEtaMinBin.push_back(param[iMap].GetDouble());
  // read eta max ranges
  param = GetParam("EtaMaxBin");
  fEtaMaxBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fEtaMaxBin.push_back(param[iMap].GetDouble());
  // read pt min value
  param = GetParam("PtMinBin");
  fPtMinBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fPtMinBin.push_back(param[iMap].GetDouble());
  // read cone size
  param = GetParam("ConeSizeBin");
  fConeSizeBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fConeSizeBin.push_back(param[iMap].GetDouble());
  // read RMS min pt
  param = GetParam("RMSPtMinBin");
  fRMSPtMinBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fRMSPtMinBin.push_back(param[iMap].GetDouble());
  // read RMS scale factor
  param = GetParam("RMSScaleFactorBin");
  fRMSScaleFactorBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fRMSScaleFactorBin.push_back(param[iMap].GetDouble());
  // read neutral pt min cut
  param = GetParam("NeutralMinEBin");
  fNeutralMinEBin.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fNeutralMinEBin.push_back(param[iMap].GetDouble());
  // read neutral pt min slope
  param = GetParam("NeutralPtSlope");
  fNeutralPtSlope.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fNeutralPtSlope.push_back(param[iMap].GetDouble());
  // read apply chs
  //param = GetParam("ApplyCHS");
  //fApplyCHS.clear();
  //for(int iMap = 0; iMap < param.GetSize(); ++iMap) fApplyCHS.push_back(param[iMap].GetBool());
  // read use charged
  param = GetParam("UseCharged");
  fUseCharged.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fUseCharged.push_back(param[iMap].GetBool());
  // read apply chs correction
  param = GetParam("ApplyLowPUCorr");
  fApplyLowPUCorr.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fApplyLowPUCorr.push_back(param[iMap].GetBool());
  // read metric id
  param = GetParam("MetricId");
  fMetricId.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fMetricId.push_back(param[iMap].GetInt());
  // scheme for combining
  param = GetParam("CombId");
  fCombId.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fCombId.push_back(param[iMap].GetInt());
  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "puppiParticles"));
  fOutputTrackArray = ExportArray(GetString("OutputArrayTracks", "puppiTracks"));
  fOutputNeutralArray = ExportArray(GetString("OutputArrayNeutrals", "puppiNeutrals"));
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
  fPuppi = new PuppiContainer(true, fUseExp, fMinPuppiWeight, puppiAlgo);
}

//------------------------------------------------------------------------------

void RunPUPPI::Finish()
{
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(fPuppi) delete fPuppi;
}

//------------------------------------------------------------------------------

void RunPUPPI::Process()
{

  Candidate *candidate, *particle;
  TLorentzVector momentum;

  //DelphesFactory *factory = GetFactory();

  // loop over input objects
  fItTrackInputArray->Reset();
  fItNeutralInputArray->Reset();
  fPVItInputArray->Reset();

  std::vector<Candidate *> InputParticles;
  InputParticles.clear();

  // take the leading vertex
  float PVZ = 0.;
  Candidate *pv = static_cast<Candidate *>(fPVItInputArray->Next());
  if(pv) PVZ = pv->Position.Z();
  // Fill input particles for puppi
  std::vector<RecoObj> puppiInputVector;
  puppiInputVector.clear();
  int lNBad = 0;
  // Loop on charge track candidate
  while((candidate = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {
    momentum = candidate->Momentum;
    RecoObj curRecoObj;
    curRecoObj.pt = momentum.Pt();
    curRecoObj.eta = momentum.Eta();
    curRecoObj.phi = momentum.Phi();
    curRecoObj.m = momentum.M();
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0)); //if(fApplyNoLep && TMath::Abs(candidate->PID) == 11) continue; //Dumb cut to minimize the nolepton on electron
    //if(fApplyNoLep && TMath::Abs(candidate->PID) == 13) continue;
    if(candidate->IsRecoPU and candidate->Charge != 0)
    { // if it comes fromPU vertexes after the resolution smearing and the dZ matching within resolution
      lNBad++;
      curRecoObj.id = 2;
      curRecoObj.vtxId = 0.7 * (fPVInputArray->GetEntries()); //Hack apply reco vtx efficiency of 70% for calibration
      if(TMath::Abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(TMath::Abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(TMath::Abs(candidate->PID) == 22)
        curRecoObj.pfType = 4;
      else
        curRecoObj.pfType = 1;
      curRecoObj.dZ = particle->Position.Z() - PVZ;
    }
    else if(!candidate->IsRecoPU && candidate->Charge != 0)
    {
      curRecoObj.id = 1; // charge from LV
      curRecoObj.vtxId = 1; // from PV
      if(TMath::Abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(TMath::Abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(TMath::Abs(candidate->PID) == 22)
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
  while((candidate = static_cast<Candidate *>(fItNeutralInputArray->Next())))
  {
    momentum = candidate->Momentum;
    RecoObj curRecoObj;
    curRecoObj.pt = momentum.Pt();
    curRecoObj.eta = momentum.Eta();
    curRecoObj.phi = momentum.Phi();
    curRecoObj.m = momentum.M();
    curRecoObj.charge = 0;
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));
    if(candidate->Charge == 0)
    {
      curRecoObj.id = 0; // neutrals have id==0
      curRecoObj.vtxId = 0; // neutrals have vtxId==0
      if(TMath::Abs(candidate->PID) == 11)
        curRecoObj.pfType = 2;
      else if(TMath::Abs(candidate->PID) == 13)
        curRecoObj.pfType = 3;
      else if(TMath::Abs(candidate->PID) == 22)
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
      fOutputArray->Add(candidate);
      if(puppiInputVector.at(it->user_index()).id == 1 or puppiInputVector.at(it->user_index()).id == 2)
        fOutputTrackArray->Add(candidate);
      else if(puppiInputVector.at(it->user_index()).id == 0)
        fOutputNeutralArray->Add(candidate);
    }
    else
    {
      std::cerr << " particle not found in the input Array --> skip " << std::endl;
      continue;
    }
  }
}
