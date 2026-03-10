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
  std::unique_ptr<PuppiContainer> fPuppi;
  // puppi parameters
  bool fApplyNoLep;
  double fMinPuppiWeight;
  bool fUseExp;

  const TObjArray *fTrackInputArray{nullptr};
  std::unique_ptr<TIterator> fItTrackInputArray;

  const TObjArray *fNeutralInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItNeutralInputArray; //!

  const TObjArray *fPVInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fPVItInputArray; //!

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

  TObjArray *fOutputArray{nullptr};
  TObjArray *fOutputTrackArray{nullptr};
  TObjArray *fOutputNeutralArray{nullptr};

  ClassDef(RunPUPPI, 1)
};

#endif
