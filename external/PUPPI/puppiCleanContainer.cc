#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "fastjet/Selector.hh"

#include <algorithm>

#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"

using namespace std;

// ------------- Constructor
puppiCleanContainer::puppiCleanContainer(std::vector<RecoObj> inParticles, 
                                         std::vector<puppiAlgoBin> puppiAlgo,
                                         float minPuppiWeight,
                                         bool  fUseExp){

    // take the input particles
    fRecoParticles_.clear();
    fRecoParticles_ = inParticles;

    // puppi algo 
    puppiAlgo_.clear();
    puppiAlgo_ = puppiAlgo; 

    // min puppi weight
    fMinPuppiWeight_ = minPuppiWeight;

    //Clear everything
    fPFParticles_.clear();
    fPFchsParticles_.clear();
    fChargedPV_.clear();
    fChargedNoPV_.clear();
    fPuppiWeights_.clear();

    fNPV_    = 1 ;
    fPVFrac_ = 0.;
    fUseExp_ = fUseExp;

    //Link to the RecoObjects --> loop on the input particles
    for (unsigned int i = 0; i < fRecoParticles_.size(); i++){
        fastjet::PseudoJet curPseudoJet;
        curPseudoJet.reset_PtYPhiM (fRecoParticles_[i].pt,fRecoParticles_[i].eta,fRecoParticles_[i].phi,fRecoParticles_[i].m);
        curPseudoJet.set_user_index(fRecoParticles_[i].id);  
        // fill vector of pseudojets for internal references
        fPFParticles_.push_back(curPseudoJet);
        if(fRecoParticles_[i].id <= 1) fPFchsParticles_.push_back(curPseudoJet);    //Remove Charged particles associated to other vertex
        if(fRecoParticles_[i].id == 1) fChargedPV_.push_back(curPseudoJet);         //Take Charged particles associated to PV
        if(fRecoParticles_[i].id == 2) fChargedNoPV_.push_back(curPseudoJet);
        if(fRecoParticles_[i].id >= 0) fPVFrac_++ ;
	if(fNPV_ < fRecoParticles_[i].vtxId) fNPV_ = fRecoParticles_[i].vtxId;

    }

    fPVFrac_ = double(fChargedPV_.size())/fPVFrac_;
}

// ------------- De-Constructor
puppiCleanContainer::~puppiCleanContainer(){}

// main function to compute puppi Event
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent(){

    // output particles
    std::vector<fastjet::PseudoJet> particles;
    particles.clear();

    std::vector<int> pPupId ; 
    std::vector<puppiParticle> partTmp ; // temp puppi particle vector; make a clone of the same particle for all the algo in which it is contained
  
    // calculate puppi metric, RMS and mean value for all the algorithms
    for(size_t iPuppiAlgo = 0; iPuppiAlgo < puppiAlgo_.size(); iPuppiAlgo++){
      getRMSAvg(iPuppiAlgo,fPFParticles_,fChargedPV_); // give all the particles in the event and the charged one
    }
  
    int npart = 0;  

    // Loop on all the incoming particles
    for(size_t iPart = 0; iPart < fPFParticles_.size(); iPart++) {

      float pWeight = 1; // default weight
      pPupId.clear();
      pPupId = getPuppiId(fPFParticles_[iPart].pt(),fPFParticles_[iPart].eta(),puppiAlgo_); // take into account only algo eta

      //////////////////////////////////////////      
      // acceptance check of the puppi algorithm
      //////////////////////////////////////////      

      if(pPupId.empty()) { // out acceptance... no algorithm found
        fPuppiWeights_.push_back(pWeight); // take the particle as it is
        fastjet::PseudoJet curjet(pWeight*fPFParticles_[iPart].px(),pWeight*fPFParticles_[iPart].py(),pWeight*fPFParticles_[iPart].pz(),pWeight*fPFParticles_[iPart].e());    
        curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
        particles.push_back(curjet); // fill the output collection                                                                                                            
	continue; //go to the next particle
      }
      
      ///////////////      
      //  PT check //
      ///////////////      
      for(size_t iAlgo = 0; iAlgo < pPupId.size(); iAlgo++){ // loop on all the available algo for that region       
	if(fPFParticles_.at(iPart).pt() < puppiAlgo_.at(pPupId.at(iAlgo)).fPtMin_){ // low momentum particles should be cut by puppi method
	pWeight = 0; // if this particle is under the pT threshold of one the algorithm, put the weight as zero      
	break;
	}
      }

      if(pWeight == 0){ 
        fPuppiWeights_.push_back(0); // puppi weight is zero
        continue;
      }
 
      /////////////////////////////////////      
      // fill the p-values for Z-vertex //
      /////////////////////////////////////      

      double pChi2 = 0;   
      if(fUseExp_){ // use vertex-z resolution
       //Compute an Experimental Puppi Weight with delta Z info (very simple example)
       if(iPart <= fRecoParticles_.size() and fRecoParticles_[iPart].id == fPFParticles_.at(iPart).user_index()){
	pChi2 = getChi2FromdZ(fRecoParticles_[iPart].dZ); // get the probability fiven the dZ of the particle wrt the leading vertex
        if(fRecoParticles_[iPart].pfType > 3) pChi2 = 0; // not use this info for neutrals
       }
      }
      

      /////////////////////////////////////      
      // found  the particle in all the algorithm
      /////////////////////////////////////      
      partTmp.clear();
      for(size_t iAlgo = 0; iAlgo < pPupId.size(); iAlgo++){ // loop on all the algo found

       int found = 0;  //  found index
       if(fabs(fPFParticles_[iPart].user_index()) <= 1 and puppiAlgo_.at(pPupId.at(iAlgo)).fUseCharged_){ // charged or neutral from PV
	 for(size_t puppiIt = 0 ; puppiIt < puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPV_.size(); puppiIt++){ // Loop on PV particles
	   if(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPV_.at(puppiIt).fPosition_ == int(iPart)){
	    partTmp.push_back(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPV_.at(puppiIt));  // take the puppi particle
            found = 1 ;
            break;
	  }
	}
       }
       else if ((fabs(fPFParticles_[iPart].user_index()) <= 1 and !puppiAlgo_.at(pPupId.at(iAlgo)).fUseCharged_) or fabs(fPFParticles_[iPart].user_index()) >= 2){
	for(size_t puppiIt = 0; puppiIt < puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPU_.size(); puppiIt++){
          if(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPU_.at(puppiIt).fPosition_ == int(iPart)){
            partTmp.push_back(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesPU_.at(puppiIt));
            found = 1;
            break;
	  }
	}
       }
      
      /////////////////////////////////////      
       // means that is inside the NULL vector for some reasons
      /////////////////////////////////////      

       if(found == 0){
	for(size_t puppiIt = 0; puppiIt < puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesNULL_.size(); puppiIt++){
          if(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesNULL_.at(puppiIt).fPosition_ == int(iPart)){
            partTmp.push_back(puppiAlgo_.at(pPupId.at(iAlgo)).fPuppiParticlesNULL_.at(puppiIt));
            found = 1 ;
            break;
	  }
	}
       }
      }

      if(partTmp.size() != pPupId.size()){ // not found the particle in one of the algorithms 
        pWeight = 1 ;
	fPuppiWeights_.push_back(pWeight);
        fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());    
        curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
        particles.push_back(curjet); // by default is one, so 4V is not chaged                                                                                                
        continue;
      }
      
      /////////////////////////////////////      
      //Check the Pval  
      /////////////////////////////////////      

      bool badPVal = false ;
      for(size_t iPuppi = 0; iPuppi < partTmp.size() ; iPuppi++){
	if(partTmp.at(iPuppi).fPval_ == -999){ // if the default is found as PVal, leave the particle as it is in the output
	 pWeight = 1 ;
         badPVal = true ;
	}
      }

      if(badPVal){
 	 fPuppiWeights_.push_back(pWeight);
         fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());    
         curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
         particles.push_back(curjet);                                                                                                                                                
         continue;
	
      }

      
      // compute combining the weight for all the algorithm
      pWeight = compute(pChi2,partTmp,puppiAlgo_,pPupId);      
      
      //Basic Weight Checks
      if( std::isinf(pWeight) || std::isnan(pWeight)){
	std::cerr << "====> Weight is nan : pt " << fPFParticles_[iPart].pt() << " -- eta : " << fPFParticles_[iPart].eta() << " -- id : " << fPFParticles_[iPart].user_index() << std::endl;
	pWeight = 1; // set the default to avoid problems
      }

      //Basic Cuts
      if(pWeight < fMinPuppiWeight_) pWeight = 0; //==> Elminate the low Weight stuff
      
      //threshold cut on the neutral Pt
      for(size_t iPuppi = 0; iPuppi < pPupId.size(); iPuppi++){
        if(fPFParticles_[iPart].user_index() == 1 && puppiAlgo_.at(pPupId.at(iPuppi)).fApplyCHS_ ) pWeight = 1; // charged from LV 
        if(fPFParticles_[iPart].user_index() == 2 && puppiAlgo_.at(pPupId.at(iPuppi)).fApplyCHS_ ) pWeight = 0; // charged from PU
	if(pWeight*fPFParticles_[iPart].pt() < getNeutralPtCut(puppiAlgo_.at(pPupId.at(iPuppi)).fNeutralMinE_,puppiAlgo_.at(pPupId.at(iPuppi)).fNeutralPtSlope_,fNPV_) && fPFParticles_[iPart].user_index() == 0 ) // if don't pass one of the algo neutral pt condition the particle is cut
       pWeight = 0; 
      }

      fPuppiWeights_.push_back(pWeight); // push back the weight

      //Now get rid of the thrown out weights for the particle collection
      if(pWeight == 0) continue; // if zero don't fill the particle in the output
      npart++;

      //Produce
      fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());           
      curjet.set_user_index(iPart);                                                                                                                            
      particles.push_back(curjet);
      
    }
    
    return particles;
      
}

// compute puppi metric, RMS and median for PU particle for each algo
void puppiCleanContainer::getRMSAvg(const int & iPuppiAlgo, std::vector<fastjet::PseudoJet> & particlesAll, std::vector<fastjet::PseudoJet> &chargedPV) { 

  std::vector<puppiParticle> puppiParticles; // puppi particles to be set for a specific algo
  puppiParticles.clear();

  // Loop on all the particles of the event  
    
  for(size_t iPart = 0; iPart < particlesAll.size(); iPart++ ) { 

    float pVal    = -999;
    bool  pPupId  = isGoodPuppiId(particlesAll[iPart].pt(),particlesAll[iPart].eta(),puppiAlgo_.at(iPuppiAlgo)); // get the puppi id algo asaf of eta and phi of the particle
    // does not exsist and algorithm for this particle, store -999 as pVal
    if(pPupId == false) continue;
    // apply CHS in puppi metric computation -> use only LV hadrons to compute the metric for each particle
    if(puppiAlgo_.at(iPuppiAlgo).fUseCharged_)  
       pVal = goodVar(particlesAll[iPart], chargedPV,    puppiAlgo_.at(iPuppiAlgo).fMetricId_,puppiAlgo_.at(iPuppiAlgo).fConeSize_);
    else if(!puppiAlgo_.at(iPuppiAlgo).fUseCharged_) 
       pVal = goodVar(particlesAll[iPart], particlesAll, puppiAlgo_.at(iPuppiAlgo).fMetricId_,puppiAlgo_.at(iPuppiAlgo).fConeSize_);

    // fill the value
    if(std::isnan(pVal) || std::isinf(pVal)) std::cout << "====>  Value is Nan " << pVal << " == " << particlesAll[iPart].pt() << " -- " << particlesAll[iPart].eta() << std::endl;
    if(std::isnan(pVal) || std::isinf(pVal)) continue;
    
    puppiParticles.push_back(puppiParticle(particlesAll.at(iPart).pt(),particlesAll.at(iPart).eta(),pVal,particlesAll.at(iPart).user_index(),iPart));
  }
  
  // set the puppi particles for the algorithm
  puppiAlgo_.at(iPuppiAlgo).setPuppiParticles(puppiParticles);
  // compute RMS, median and mean value  
  computeMedRMS(iPuppiAlgo);
  
}


float puppiCleanContainer::goodVar(const fastjet::PseudoJet & particle, const std::vector<fastjet::PseudoJet> & particleAll, const int & pPupId, const float & coneSize) {
  float lPup = 0;
  lPup = var_within_R(pPupId,particleAll,particle,coneSize);
  return lPup;
}

float puppiCleanContainer::var_within_R(const int & pPupId, const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R){

  if(pPupId == -1) return 1;
  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  float var = 0;

  for(size_t iPart = 0; iPart < near_particles.size(); iPart++){

    double pDEta = near_particles[iPart].eta()-centre.eta();
    double pDPhi = fabs(near_particles[iPart].phi()-centre.phi());
    if(pDPhi > 2.*3.14159265-pDPhi) pDPhi = 2.*3.14159265-pDPhi;
    double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);

    if(pDR < 0.0001) continue;
    if(pDR == 0)    continue;

    if(pPupId == 0) var += (near_particles[iPart].pt()/(pDR*pDR));
    if(pPupId == 1) var += near_particles[iPart].pt();
    if(pPupId == 2) var += (1./pDR)*(1./pDR);
    if(pPupId == 3) var += (1./pDR)*(1./pDR);
    if(pPupId == 4) var += near_particles[iPart].pt();
    if(pPupId == 5) var += (near_particles[iPart].pt()/pDR)*(near_particles[iPart].pt()/pDR);
  }

  if(pPupId == 0 && var != 0) var = log(var);
  if(pPupId == 3 && var != 0) var = log(var);
  if(pPupId == 5 && var != 0) var = log(var);
  return var;
 
}


float puppiCleanContainer::pt_within_R(const std::vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet & centre, const float & R){

  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  double answer = 0.0;
  for(size_t iPart = 0; iPart<near_particles.size(); iPart++){
    answer += near_particles[iPart].pt();
  }

  return answer;
}

fastjet::PseudoJet puppiCleanContainer::flow_within_R(const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R){

  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  fastjet::PseudoJet flow;
  for(unsigned int i=0; i<near_particles.size(); i++){
    flow += near_particles[i];
  }
  return flow;

}

// compute median, mean value and RMS for puppi
void puppiCleanContainer::computeMedRMS(const int & puppiAlgo) {

  if(puppiAlgo > int(puppiAlgo_.size())  ) return;
  if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size() == 0) return;

  // sort in pVal increasing order
  std::sort(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.begin(),puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.end(),puppiValSort());

  // if apply correction
  float lCorr = 1.;
  if(puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_) lCorr *= 1.-fPVFrac_;

  // count the position of the last particle with pval zero coming from PU
  int lNum0 = 0;
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size(); i0++) {
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_[i0].fPval_ == 0) lNum0 = i0;
  }

  // take the median value on PU particles
  int lNHalfway = lNum0 + int(float(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size()-lNum0)*0.50*lCorr);
  puppiAlgo_.at(puppiAlgo).fMedian_ = puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(lNHalfway).fPval_;
  float lMed = puppiAlgo_.at(puppiAlgo).fMedian_; //Just to make the readability easier

  // take the RMS
  int lNRMS = 0;
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size(); i0++) {
    puppiAlgo_.at(puppiAlgo).fMean_ += puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_;
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_ == 0) continue;
    if(!puppiAlgo_.at(puppiAlgo).fUseCharged_ && puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_ && puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_ > lMed) continue;
    lNRMS++;
    puppiAlgo_.at(puppiAlgo).fRMS_ += (puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_-lMed)*( puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_-lMed);
  }

  puppiAlgo_.at(puppiAlgo).fMean_ /= puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size();
  if(lNRMS > 0) puppiAlgo_.at(puppiAlgo).fRMS_/=lNRMS;
  if(puppiAlgo_.at(puppiAlgo).fRMS_ == 0) puppiAlgo_.at(puppiAlgo).fRMS_ = 1e-5;
  puppiAlgo_.at(puppiAlgo).fRMS_  = sqrt(puppiAlgo_.at(puppiAlgo).fRMS_);
  puppiAlgo_.at(puppiAlgo).fRMS_ *= puppiAlgo_.at(puppiAlgo).fRMSScaleFactor_;

  if(!puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_) return;

  //Adjust the p-value to correspond to the median
  std::sort(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.begin(),puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.end(),puppiValSort());
  int lNPV = 0; 
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.size(); i0++){ 
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_[i0].fPval_ <= lMed ) lNPV++;
  }

  // it helps in puppi the median value close to the mean one when a lot of pval 0 are present
  float lAdjust = 1.5*float(lNPV)/float(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.size()+puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size());
  if(lAdjust > 0) puppiAlgo_.at(puppiAlgo).fMedian_ -= sqrt(ROOT::Math::chisquared_quantile(lAdjust,1.)*puppiAlgo_.at(puppiAlgo).fRMS_);

}

float puppiCleanContainer::getNeutralPtCut(const float & fNeutralMinE, const float & fNeutralPtSlope, const int & fNPV) {
  return fNeutralMinE + fNPV * fNeutralPtSlope;
}

// take the type of algorithm : return a vector since more than one algo can be defined for the same eta region
std::vector<int> puppiCleanContainer::getPuppiId(const float & pt, const float & eta, const std::vector<puppiAlgoBin> & puppiAlgos){
  std::vector<int> PuppiId ;
  for(size_t iPuppiAlgo = 0; iPuppiAlgo < puppiAlgos.size() ; iPuppiAlgo++){
    if(fabs(eta) <= puppiAlgos[iPuppiAlgo].fEtaMin_) continue;
    if(fabs(eta) > puppiAlgos[iPuppiAlgo].fEtaMax_) continue;
    PuppiId.push_back(int(iPuppiAlgo));
  }
  return PuppiId;  
}

//check if a particle is good for an Algo definition
bool puppiCleanContainer::isGoodPuppiId(const float & pt, const float & eta, const puppiAlgoBin & puppiAlgo){
  if(fabs(eta) <= puppiAlgo.fEtaMin_) return false;
  if(fabs(eta) > puppiAlgo.fEtaMax_) return false;
  return true;

}



// ----------------------
float puppiCleanContainer::compute(const float & chi2, const std::vector<puppiParticle> & particles, const std::vector<puppiAlgoBin> & puppiAlgos, const std::vector<int> & pPupId) {

  if(particles.size() != pPupId.size() ) return 0; // default check

  float lVal  = 0.;
  float lPVal = 1.;
  int   lNDOF = 0;

  for( size_t iAlgo = 0; iAlgo < pPupId.size(); iAlgo++){

    if(puppiAlgos.at(pPupId.at(iAlgo)).fPuppiParticlesPU_.size() + puppiAlgos.at(pPupId.at(iAlgo)).fPuppiParticlesPV_.size()  == 0) return 1; 

    if(iAlgo > 0 ){
      float pPVal = ROOT::Math::chisquared_cdf(lVal,lNDOF); // take a chi2 value since the should be multiplied (multiply weight and not summing chi2)
      lPVal *= pPVal;
      lNDOF  = 0;
      lVal   = 0; 
    }

    if(puppiAlgos.at(pPupId.at(iAlgo)).fMetricId_ == -1) continue;

    float pVal = particles.at(iAlgo).fPval_ ;

    if(puppiAlgos.at(pPupId.at(iAlgo)).fMetricId_ == 0 && pVal == 0) pVal = puppiAlgos.at(pPupId.at(iAlgo)).fMedian_;
    if(puppiAlgos.at(pPupId.at(iAlgo)).fMetricId_ == 3 && pVal == 0) pVal = puppiAlgos.at(pPupId.at(iAlgo)).fMedian_;
    if(puppiAlgos.at(pPupId.at(iAlgo)).fMetricId_ == 5 && pVal == 0) pVal = puppiAlgos.at(pPupId.at(iAlgo)).fMedian_;
    
    lVal += (pVal-puppiAlgos.at(pPupId.at(iAlgo)).fMedian_)*(fabs(pVal-puppiAlgos.at(pPupId.at(iAlgo)).fMedian_))/puppiAlgos.at(pPupId.at(iAlgo)).fRMS_/puppiAlgos.at(pPupId.at(iAlgo)).fRMS_;
    lNDOF++;
    if(chi2 != 0) lNDOF++; 
    if(chi2 != 0) lVal+=chi2; //Add external Chi2 to first element
  }

  lPVal *= ROOT::Math::chisquared_cdf(lVal,lNDOF);
  return lPVal;

}

float puppiCleanContainer::getChi2FromdZ(float iDZ) {
   //We need to obtain prob of PU + (1-Prob of LV)
   // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm (its really more like 1mm)
   //double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.2)*2.; //*2 is to do it double sided
   //Take iDZ to be corrected by sigma already
   double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),1.)*2.; //*2 is to do it double sided
   double lProbPU = 1-lProbLV;
   if(lProbPU <= 0) lProbPU = 1e-16; //Quick Trick to through out infs
   if(lProbPU >= 0) lProbPU = 1-1e-16; //Ditto
   double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
   lChi2PU*=lChi2PU;
   return lChi2PU;
}
