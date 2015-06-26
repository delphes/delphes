#ifndef RunPUPPI_h
#define RunPUPPI_h

#include "classes/DelphesModule.h"
#include <vector>

class TObjArray;
class TIterator;


class RunPUPPI: public DelphesModule {

 public:

  RunPUPPI();
  ~RunPUPPI();

  void Init();
  void Process();
  void Finish();
  
 private:

  TIterator *fItTrackInputArray;
  TIterator *fItNeutralInputArray; //!
  TIterator *fPVItInputArray; //!                                                                                                                                                      
  
  const TObjArray *fTrackInputArray;
  const TObjArray *fNeutralInputArray; //!
  const TObjArray *fPVInputArray; //!                                                                                                                                                     
 
  // puppi parameters
  float fMinPuppiWeight;
  bool fUseExp;
  
  std::vector<float> fEtaMinBin ;
  std::vector<float> fEtaMaxBin ;
  std::vector<float> fPtMinBin ;
  std::vector<float> fConeSizeBin ;
  std::vector<float> fRMSPtMinBin ;
  std::vector<float> fRMSScaleFactorBin ;
  std::vector<float> fNeutralMinEBin;
  std::vector<float> fNeutralPtSlope;
  std::vector<bool>  fApplyCHS;
  std::vector<bool>  fUseCharged;
  std::vector<bool>  fApplyLowPUCorr;
  std::vector<int>   fMetricId;

  TObjArray *fOutputArray;
  TObjArray *fOutputTrackArray;
  TObjArray *fOutputNeutralArray;

  ClassDef(RunPUPPI, 1)
};

#endif
