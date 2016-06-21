#include "PuppiAlgo.hh"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TMath.h"

PuppiAlgo::PuppiAlgo(AlgoObj &iAlgo) { 
  fEtaMin             = iAlgo.etaMin;
  fEtaMax             = iAlgo.etaMax;
  fPtMin              = iAlgo.ptMin;
  fNeutralPtMin       = iAlgo.minNeutralPt;
  fNeutralPtSlope     = iAlgo.minNeutralPtSlope;
  fRMSEtaSF           = iAlgo.rmsEtaSF;
  fMedEtaSF           = iAlgo.medEtaSF;
  fEtaMaxExtrap       = iAlgo.etaMaxExtrap;
  std::vector<AlgoSubObj> lAlgos = iAlgo.subAlgos;
  fNAlgos = lAlgos.size();
  //std::cout << "==> "  << fEtaMin << " - " << fEtaMax << " - " << fPtMin  << " - " << fNeutralPtMin  << " - " << fNeutralPtSlope  << " - " << fRMSEtaSF  << " - " << std::endl;
  for(unsigned int i0 = 0; i0 < lAlgos.size(); i0++)  { 
    int    pAlgoId      = lAlgos[i0].metricId;
    bool   pCharged     = lAlgos[i0].useCharged;
    bool   pWeight0     = lAlgos[i0].applyLowPUCorr;
    int    pComb        = lAlgos[i0].combId;
    double pConeSize    = lAlgos[i0].coneSize;
    double pRMSPtMin    = lAlgos[i0].rmsPtMin;
    double pRMSSF       = lAlgos[i0].rmsScaleFactor;
    //std::cout << "Algo==> " << i0 << " - " << pAlgoId << " - " << pCharged << " - " << pWeight0 << " - " << pComb << " - " << pConeSize << " - " << pRMSPtMin << " - " << std::endl;
    fAlgoId        .push_back(pAlgoId);
    fCharged       .push_back(pCharged);
    fAdjust        .push_back(pWeight0);
    fCombId        .push_back(pComb);
    fConeSize      .push_back(pConeSize);
    fRMSPtMin      .push_back(pRMSPtMin);
    fRMSScaleFactor.push_back(pRMSSF);
    double pRMS  = 0; 
    double pMed  = 0; 
    double pMean = 0;
    int    pNCount = 0; 
    fRMS   .push_back(pRMS);
    fMedian.push_back(pMed);
    fMean  .push_back(pMean);
    fNCount.push_back(pNCount);
    /*
      tmprms .clear();
      tmpmed .clear();
      for (unsigned int j0 = 0; j0 < fEtaMin.size(); j0++){
      tmprms.push_back(pRMS);
      tmpmed.push_back(pMed);
      }
      fRMS_perEta.push_back(tmprms);
      fMedian_perEta.push_back(tmpmed);
    */
  }
}
PuppiAlgo::~PuppiAlgo() { 
  fPups  .clear();
  fPupsPV.clear();
}
void PuppiAlgo::reset() { 
  fPups  .clear();
  fPupsPV.clear();
  for(unsigned int i0 = 0; i0 < fNAlgos; i0++) { 
    fMedian[i0] =  0; 
    fRMS   [i0] =  0;
    fMean  [i0] =  0;
    fNCount[i0] =  0;
  }
}
void PuppiAlgo::add(const fastjet::PseudoJet &iParticle,const double &iVal,const unsigned int iAlgo) { 
  if(iParticle.pt() < fRMSPtMin[iAlgo]) return;
  if(fCharged[iAlgo] && fabs(iParticle.user_index())  < 1) return;
  if(fCharged[iAlgo] && (fabs(iParticle.user_index()) >=1 && fabs(iParticle.user_index()) <=2)) fPupsPV.push_back(iVal);
  if(fCharged[iAlgo] && fabs(iParticle.user_index()) < 3) return;
  fPups.push_back(iVal);
  fNCount[iAlgo]++;
  /* eta extrap (cmssw version)
     if ((std::abs(iParticle.eta()) < fEtaMaxExtrap) && (std::abs(puppi_register) >= 3)){
        fPups.push_back(iVal);
        // fPupsPV.push_back(iVal);        
        fNCount[iAlgo]++;
	}
    // for the low PU case, correction.  for checking that the PU-only median will be below the PV particles
    if(std::abs(iParticle.eta()) < fEtaMaxExtrap && (std::abs(puppi_register) >=1 && std::abs(puppi_register) <=2)) fPupsPV.push_back(iVal);
  */
}
void PuppiAlgo::computeMedRMS(const unsigned int &iAlgo,const double &iPVFrac) { 
  if(iAlgo >= fNAlgos   ) return;
  if(fNCount[iAlgo] == 0) return;
  int lNBefore = 0; 
  for(unsigned int i0 = 0; i0 < iAlgo; i0++) lNBefore += fNCount[i0];
  std::sort(fPups.begin()+lNBefore,fPups.begin()+lNBefore+fNCount[iAlgo]);
  int lNum0 = 0;
  for(int i0 = lNBefore; i0 < lNBefore+fNCount[iAlgo]; i0++) { 
    if(fPups[i0] == 0) lNum0 = i0-lNBefore; 
  }
  //lNum0 = 0; 
  int lNHalfway = lNBefore + lNum0 + int( double( fNCount[iAlgo]-lNum0 )*0.50);
  fMedian[iAlgo] = fPups[lNHalfway];
  double lMed = fMedian[iAlgo];  //Just to make the readability easier

  int lNRMS = 0; 
  for(int i0 = lNBefore; i0 < lNBefore+fNCount[iAlgo]; i0++) {
    fMean[iAlgo] += fPups[i0];
    if(fPups[i0] == 0) continue;
    if(fAdjust[iAlgo] && fPups[i0] > lMed) continue;
    lNRMS++;
    fRMS [iAlgo] += (fPups[i0]-lMed)*(fPups[i0]-lMed);
  }
  fMean[iAlgo]/=fNCount[iAlgo];
  if(lNRMS > 0) fRMS [iAlgo]/=lNRMS;
  if(fRMS[iAlgo] == 0) fRMS[iAlgo] = 1e-5;

  fRMS [iAlgo] = sqrt(fRMS[iAlgo]);
  fRMS [iAlgo] *= fRMSScaleFactor[iAlgo];
  if(fAdjust[iAlgo]){ 
    //Adjust the p-value to correspond to the median
    std::sort(fPupsPV.begin(),fPupsPV.end());
    int lNPV = 0; 
    for(unsigned int i0 = 0; i0 < fPupsPV.size(); i0++) if(fPupsPV[i0] <= lMed ) lNPV++;
    double lAdjust = double(lNPV)/double(lNPV+0.5*fNCount[iAlgo]);
    if(lAdjust > 0) {
      fMedian[iAlgo] -= sqrt(ROOT::Math::chisquared_quantile(lAdjust,1.)*fRMS[iAlgo]);
      fRMS[iAlgo]    -= sqrt(ROOT::Math::chisquared_quantile(lAdjust,1.)*fRMS[iAlgo]);
    }        
  }
  /* eta extrapolated version
  for (unsigned int j0 = 0; j0 < fEtaMin.size(); j0++){
    fRMS_perEta[iAlgo][j0]    = fRMS[iAlgo]*fRMSEtaSF[j0];
    fMedian_perEta[iAlgo][j0] = fMedian[iAlgo]*fMedEtaSF[j0];
  }    
  */
}
//This code is probably a bit confusing
double PuppiAlgo::compute(std::vector<double> &iVals,double iChi2) { 
  if(fAlgoId[0] == -1) return 1;
  double lVal  = 0.;
  double lPVal = 1.;
  int    lNDOF = 0; 
  for(unsigned int i0 = 0; i0 < fNAlgos; i0++) { 
    if(fNCount[i0] == 0) return 1.;   //in the NoPU case return 1.
    if(fCombId[i0] == 1 && i0 > 0) {  //Compute the previous p-value so that p-values can be multiplieed
      double pPVal = ROOT::Math::chisquared_cdf(lVal,lNDOF);
      lPVal *= pPVal;
      lNDOF = 0; 
      lVal  = 0; 
    }
    double pVal = iVals[i0];
    //Special Check for any algo with log(0) 
    if(fAlgoId[i0] == 0 && iVals[i0] == 0) pVal = fMedian[i0];
    if(fAlgoId[i0] == 3 && iVals[i0] == 0) pVal = fMedian[i0];
    if(fAlgoId[i0] == 5 && iVals[i0] == 0) pVal = fMedian[i0];
    lVal += (pVal-fMedian[i0])*(fabs(pVal-fMedian[i0]))/fRMS[i0]/fRMS[i0];
    lNDOF++;
    if(i0 == 0 && iChi2 != 0) lNDOF++;      //Add external Chi2 to first element
    if(i0 == 0 && iChi2 != 0) lVal+=iChi2;  //Add external Chi2 to first element
  }
  //Top it off with the last calc
  lPVal *= ROOT::Math::chisquared_cdf(lVal,lNDOF);
  return lPVal;
}
double PuppiAlgo::neutralPt(int iNPV) { 
  return fNeutralPtMin + iNPV * fNeutralPtSlope;
}
int PuppiAlgo::numAlgos() { 
  return fNAlgos;
}
double PuppiAlgo::ptMin() { 
  return fPtMin;
}
double PuppiAlgo::etaMin() { 
  return fEtaMin;
}
double PuppiAlgo::etaMax() { 
  return fEtaMax;
}
int PuppiAlgo::algoId(const unsigned int &iAlgo) { 
  assert(iAlgo < fNAlgos);
  return fAlgoId[iAlgo];
}
bool PuppiAlgo::isCharged(const unsigned int &iAlgo) { 
  assert(iAlgo < fNAlgos);
  return fCharged[iAlgo];
}
double PuppiAlgo::coneSize(const unsigned int &iAlgo) { 
  assert(iAlgo < fNAlgos);
  return fConeSize[iAlgo];
}
