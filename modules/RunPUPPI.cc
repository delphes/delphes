#include "modules/RunPUPPI.h"

#include "PUPPI/puppiCleanContainer.hh"
#include "PUPPI/RecoObj.hh"
#include "PUPPI/puppiParticle.hh"
#include "PUPPI/puppiAlgoBin.hh"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

//------------------------------------------------------------------------------
RunPUPPI::RunPUPPI() :
  fItTrackInputArray(0), 
  fItNeutralInputArray(0)
{}

//------------------------------------------------------------------------------
RunPUPPI::~RunPUPPI(){}

//------------------------------------------------------------------------------

void RunPUPPI::Init(){

  // input collection
  fTrackInputArray     = ImportArray(GetString("TrackInputArray", "Calorimeter/towers"));
  fItTrackInputArray   = fTrackInputArray->MakeIterator();
  fNeutralInputArray   = ImportArray(GetString("NeutralInputArray", "Calorimeter/towers"));
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();
  fPVInputArray        = ImportArray(GetString("PVInputArray", "PV"));
  fPVItInputArray      = fPVInputArray->MakeIterator();


  // puppi parameters                                     
  fMinPuppiWeight = GetDouble("MinPuppiWeight", 0.01);
  fUseExp         = GetBool("UseExp", false);

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
  param = GetParam("ApplyCHS");
  fApplyCHS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fApplyCHS.push_back(param[iMap].GetBool());

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

  // create output array
  fOutputArray        = ExportArray(GetString("OutputArray", "puppiParticles"));
  fOutputTrackArray   = ExportArray(GetString("OutputArrayTracks", "puppiTracks"));
  fOutputNeutralArray = ExportArray(GetString("OutputArrayNeutrals", "puppiNeutrals"));
}

//------------------------------------------------------------------------------

void RunPUPPI::Finish(){
  if(fItTrackInputArray)   delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
}

//------------------------------------------------------------------------------

void RunPUPPI::Process(){

  Candidate *candidate, *particle;
  TLorentzVector momentum;

  DelphesFactory *factory = GetFactory();

  // loop over input objects
  fItTrackInputArray->Reset();
  fItNeutralInputArray->Reset();
  fPVItInputArray->Reset();

  std::vector<Candidate *> InputParticles;
  InputParticles.clear();

  // take the leading vertex 
  float PVZ = 0.;
  Candidate *pv = static_cast<Candidate*>(fPVItInputArray->Next());
  if (pv) PVZ = pv->Position.Z();

  // Fill input particles for puppi
  std::vector<RecoObj> puppiInputVector;
  puppiInputVector.clear();

  // Loop on charge track candidate
  while((candidate = static_cast<Candidate*>(fItTrackInputArray->Next()))){   

      momentum = candidate->Momentum;
      
      RecoObj curRecoObj;
      curRecoObj.pt  = momentum.Pt();
      curRecoObj.eta = momentum.Eta();
      curRecoObj.phi = momentum.Phi();
      curRecoObj.m   = momentum.M();  
      particle = static_cast<Candidate*>(candidate->GetCandidates()->Last());
      if (candidate->IsRecoPU and candidate->Charge !=0) { // if it comes fromPU vertexes after the resolution smearing and the dZ matching within resolution
	curRecoObj.id    = 2;
	curRecoObj.vtxId = candidate->IsPU;
	if(TMath::Abs(candidate->PID) == 11)      curRecoObj.pfType = 2;
	else if(TMath::Abs(candidate->PID) == 13) curRecoObj.pfType = 3;
	else if(TMath::Abs(candidate->PID) == 22) curRecoObj.pfType = 4;
        else curRecoObj.pfType = 1;
        curRecoObj.dZ = particle->Position.Z()-PVZ;
      } 
      else if(!candidate->IsRecoPU and candidate->Charge !=0) {
	curRecoObj.id    = 1;  // charge from LV
        curRecoObj.vtxId = 1; // from PV
	if(TMath::Abs(candidate->PID) == 11)      curRecoObj.pfType = 2;
	else if(TMath::Abs(candidate->PID) == 13) curRecoObj.pfType = 3;
	else if(TMath::Abs(candidate->PID) == 22) curRecoObj.pfType = 4;
        else curRecoObj.pfType = 1;
        curRecoObj.dZ = particle->Position.Z()-PVZ;
      }
      else {
	std::cerr<<" RunPUPPI: problem with a charged track --> it has charge 0 "<<std::endl;
        continue;
      }

      puppiInputVector.push_back(curRecoObj);
      InputParticles.push_back(candidate);
  }

  // Loop on neutral calo cells 
  while((candidate = static_cast<Candidate*>(fItNeutralInputArray->Next()))){

      momentum = candidate->Momentum;

      RecoObj curRecoObj;
      curRecoObj.pt  = momentum.Pt();
      curRecoObj.eta = momentum.Eta();
      curRecoObj.phi = momentum.Phi();
      curRecoObj.m   = momentum.M();
      particle = static_cast<Candidate*>(candidate->GetCandidates()->Last());


      if(candidate->Charge == 0){
	curRecoObj.id    = 0; // neutrals have id==0
	curRecoObj.vtxId = 0; // neutrals have vtxId==0 
	if(TMath::Abs(candidate->PID) == 11)      curRecoObj.pfType = 2;
	else if(TMath::Abs(candidate->PID) == 13) curRecoObj.pfType = 3;
	else if(TMath::Abs(candidate->PID) == 22) curRecoObj.pfType = 4;
        else curRecoObj.pfType = 5;
        curRecoObj.dZ = particle->Position.Z()-PVZ;
      }
      else{
	std::cerr<<" RunPUPPI: problem with a neutrals cells --> it has charge !=0 "<<std::endl;
        continue;
      }
      puppiInputVector.push_back(curRecoObj);
      InputParticles.push_back(candidate);
  }

  // Create algorithm list for puppi
  std::vector<puppiAlgoBin> puppiAlgo;
  if(puppiAlgo.empty()){
   if(!(fEtaMinBin.size() == fEtaMaxBin.size() and fEtaMinBin.size() == fPtMinBin.size() and fEtaMinBin.size() == fConeSizeBin.size() and fEtaMinBin.size() == fRMSPtMinBin.size()
       and fEtaMinBin.size() == fRMSScaleFactorBin.size() and fEtaMinBin.size() == fNeutralMinEBin.size() and  fEtaMinBin.size() == fNeutralPtSlope.size() 
       and fEtaMinBin.size() == fApplyCHS.size()  and fEtaMinBin.size() == fUseCharged.size()
       and fEtaMinBin.size() == fApplyLowPUCorr.size() and fEtaMinBin.size() == fMetricId.size())) {
    std::cerr<<" Error in PUPPI configuration, algo info should have the same size --> exit from the code"<<std::endl;
    std::exit(EXIT_FAILURE);
   } 

   for( size_t iAlgo =  0 ; iAlgo < fEtaMinBin.size() ; iAlgo++){
    puppiAlgoBin algoTmp ;
    algoTmp.fEtaMin_ = fEtaMinBin.at(iAlgo);
    algoTmp.fEtaMax_ = fEtaMaxBin.at(iAlgo);
    algoTmp.fPtMin_  = fPtMinBin.at(iAlgo);
    algoTmp.fConeSize_        = fConeSizeBin.at(iAlgo);
    algoTmp.fRMSPtMin_        = fRMSPtMinBin.at(iAlgo);
    algoTmp.fRMSScaleFactor_  = fRMSScaleFactorBin.at(iAlgo);
    algoTmp.fNeutralMinE_     = fNeutralMinEBin.at(iAlgo);
    algoTmp.fNeutralPtSlope_  = fNeutralPtSlope.at(iAlgo);
    algoTmp.fApplyCHS_        = fApplyCHS.at(iAlgo);
    algoTmp.fUseCharged_      = fUseCharged.at(iAlgo);
    algoTmp.fApplyLowPUCorr_  = fApplyLowPUCorr.at(iAlgo);
    algoTmp.fMetricId_        = fMetricId.at(iAlgo);
    if(std::find(puppiAlgo.begin(),puppiAlgo.end(),algoTmp) != puppiAlgo.end()) continue;    
    puppiAlgo.push_back(algoTmp);     
   }
  }  

  // Create PUPPI container
  puppiCleanContainer curEvent(puppiInputVector,puppiAlgo,fMinPuppiWeight,fUseExp);
  std::vector<fastjet::PseudoJet> puppiParticles = curEvent.puppiEvent();

  // Loop on final particles
  for (std::vector<fastjet::PseudoJet>::iterator it = puppiParticles.begin() ; it != puppiParticles.end() ; it++) {
    if(it->user_index() <= int(InputParticles.size())){      
      candidate = static_cast<Candidate *>(InputParticles.at(it->user_index())->Clone());
      candidate->Momentum.SetPxPyPzE(it->px(),it->py(),it->pz(),it->e());
      fOutputArray->Add(candidate);
      if( puppiInputVector.at(it->user_index()).id == 1 or puppiInputVector.at(it->user_index()).id == 2) fOutputTrackArray->Add(candidate);
      else if (puppiInputVector.at(it->user_index()).id == 0) fOutputNeutralArray->Add(candidate);
    }
    else{ 
      std::cerr<<" particle not found in the input Array --> skip "<<std::endl;
      continue;
    }     
  }

}
