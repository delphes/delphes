#ifndef RunPUPPI_h
#define RunPUPPI_h

#include "classes/DelphesModule.h"
#include <vector>

class Candidate;
class PuppiContainer;

class RunPUPPI : public DelphesModule
{

public:
  RunPUPPI() = default;

  void Init();
  void Process();
  void Finish();

private:
  PuppiContainer *fPuppi;
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

  InputHandle<std::vector<Candidate> > fTrackInputArray; //!
  InputHandle<std::vector<Candidate> > fNeutralInputArray; //!
  InputHandle<std::vector<Candidate> > fPVInputArray; //!
  OutputHandle<std::vector<Candidate> > fOutputArray; //!
  OutputHandle<std::vector<Candidate> > fOutputTrackArray; //!
  OutputHandle<std::vector<Candidate> > fOutputNeutralArray; //!

  ClassDef(RunPUPPI, 1)
};

#endif
