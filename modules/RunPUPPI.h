#ifndef RunPUPPI_h
#define RunPUPPI_h

#include "classes/DelphesModule.h"
#include <vector>

class TObjArray;
class TIterator;
class PuppiContainer;

class RunPUPPI: public DelphesModule
{

public:
  RunPUPPI();
  ~RunPUPPI();

  void Init();
  void Process();
  void Finish();

private:
  TIterator *fItTrackInputArray = nullptr;
  TIterator *fItNeutralInputArray = nullptr; //!
  TIterator *fPVItInputArray = nullptr; //!

  const TObjArray *fTrackInputArray = nullptr;
  const TObjArray *fNeutralInputArray = nullptr; //!
  const TObjArray *fPVInputArray = nullptr; //!
  PuppiContainer *fPuppi = nullptr;
  // puppi parameters
  bool fApplyNoLep;
  double fMinPuppiWeight;
  bool fUseExp;

  std::vector<float> fEtaMinBin;
  std::vector<float> fEtaMaxBin;
  std::vector<float> fPtMinBin;
  std::vector<float> fConeSizeBin;
  std::vector<float> fRMSPtMinBin;
  std::vector<float> fRMSScaleFactorBin;
  std::vector<float> fNeutralMinEBin;
  std::vector<float> fNeutralPtSlope;
  std::vector<bool> fApplyCHS;
  std::vector<bool> fUseCharged;
  std::vector<bool> fApplyLowPUCorr;
  std::vector<int> fMetricId;
  std::vector<int> fCombId;

  TObjArray *fOutputArray = nullptr;
  TObjArray *fOutputTrackArray = nullptr;
  TObjArray *fOutputNeutralArray = nullptr;

  ClassDef(RunPUPPI, 1)
};

#endif
