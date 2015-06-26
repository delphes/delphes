#ifndef PUPPIALGOBIN_HH
#define PUPPIALGOBIN_HH

#include <iostream>
#include <vector>
#include <cmath>

//...................... class to identify the puppi Algorithm parameters to be used                                                                                                       
class puppiAlgoBin {

 public :

  puppiAlgoBin(){

    fEtaMin_   = 0.;
    fEtaMax_   = 0.;
    fPtMin_    = 0.;
    fConeSize_ = 0.;
    fRMSPtMin_ = 0.;
    fRMSScaleFactor_ = 1.0;
    fNeutralMinE_    = 0.;
    fNeutralPtSlope_ = 0.;
    fApplyCHS_   = true;
    fUseCharged_ = false ;
    fApplyLowPUCorr_ = false;
    fMetricId_ = 0;
    fRMS_      = 0.;
    fMean_     = 0.;
    fMedian_   = 0.;
  };

  puppiAlgoBin(const float & etaMin,         const float & etaMax, const float & ptMin, const float & coneSize, const float & RMSPtMin,
               const float & RMSScaleFactor, const float & neutralMinE, const float & neutralPtSlope,	                                     
               const bool  & applyCHS,       const bool & useCharged, const bool & applyLowPUCorr, const int & metricId):
    fEtaMin_(etaMin),
    fEtaMax_(etaMax),
    fPtMin_(ptMin),
    fConeSize_(coneSize),
    fRMSPtMin_(RMSPtMin),
    fRMSScaleFactor_(RMSScaleFactor),
    fNeutralMinE_(neutralMinE),
    fNeutralPtSlope_(neutralPtSlope),
    fApplyCHS_(applyCHS),
    fUseCharged_(useCharged),
    fApplyLowPUCorr_(applyLowPUCorr),
    fMetricId_(metricId)
  {
    fRMS_    = 0.;
    fMean_   = 0.;
    fMedian_ = 0.;

  };

  virtual ~puppiAlgoBin(){};

  bool operator == (const puppiAlgoBin & algo2) const {
    if( fEtaMin_ == algo2.fEtaMin_ and 
        fEtaMax_ == algo2.fEtaMax_ and 
        fApplyCHS_ == algo2.fApplyCHS_ and 
        fUseCharged_ == algo2.fUseCharged_ and 
        fApplyLowPUCorr_ == algo2.fApplyLowPUCorr_ and 
        fMetricId_ == algo2.fMetricId_ and 
        fPtMin_ == algo2.fPtMin_ and 
        fConeSize_ == algo2.fConeSize_ and 
        fRMSPtMin_ == algo2.fRMSPtMin_ and 
        fRMSScaleFactor_ == algo2.fRMSScaleFactor_ and 
        fNeutralMinE_ == algo2.fNeutralMinE_ and 
        fNeutralPtSlope_ == algo2.fNeutralPtSlope_ ) return true;
    else return false;
  };

  void setPuppiParticles(const std::vector<puppiParticle> & puppiParticles){ // set the particles used by the current algorithm                                             
    fPuppiParticlesPU_.clear();
    fPuppiParticlesPV_.clear();
    fPuppiParticlesNULL_.clear();
  
    for(size_t iPart = 0; iPart < puppiParticles.size(); iPart++){ // loop on the puppiParticles                            
      if(puppiParticles.at(iPart).fPval_ == -999){ // default PVal
        fPuppiParticlesNULL_.push_back(puppiParticles.at(iPart));        
        continue ; 
      }
                                                 
      if(puppiParticles.at(iPart).fPt_ <  fRMSPtMin_){  // should be discarded in the RMS computation
        fPuppiParticlesNULL_.push_back(puppiParticles.at(iPart)); 
        continue ; 
      }
     
      if(puppiParticles.at(iPart).fPt_ < fPtMin_){ // under the chosen pt cut
        fPuppiParticlesNULL_.push_back(puppiParticles.at(iPart)); 
        continue ; 
      }

      // zero are neutral particle, 1 LV charged, 2PU      
      if(fUseCharged_ && (fabs(puppiParticles.at(iPart).fParticleId_) >=0 && fabs(puppiParticles.at(iPart).fParticleId_) <=1)) fPuppiParticlesPV_.push_back(puppiParticles.at(iPart));
      if(fUseCharged_ && fabs(puppiParticles.at(iPart).fParticleId_) < 2) continue;
      fPuppiParticlesPU_.push_back(puppiParticles.at(iPart));
    }
  };

  float fEtaMin_;
  float fEtaMax_;
  float fPtMin_ ;
  float fConeSize_;
  float fRMSPtMin_;
  float fRMSScaleFactor_;
  float fNeutralMinE_;
  float fNeutralPtSlope_;

  bool  fApplyCHS_;
  bool  fUseCharged_;
  bool  fApplyLowPUCorr_;

  int   fMetricId_;

  float fRMS_;
  float fMean_;
  float fMedian_;

  std::vector<puppiParticle> fPuppiParticlesPU_;
  std::vector<puppiParticle> fPuppiParticlesPV_;
  std::vector<puppiParticle> fPuppiParticlesNULL_;

};

#endif
