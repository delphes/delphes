#ifndef TrackSmearing_h
#define TrackSmearing_h

/** \class TrackSmearing
 *
 *  Performs d0, dZ, p, Theta, Phi smearing of tracks.
 *
 *
 *
 *  \author A. Hart, M. Selvaggi
 *
 */

#include "classes/DelphesModule.h"

class TIterator;
class TObjArray;
class DelphesFormula;

class TrackSmearing: public DelphesModule
{
public:
  TrackSmearing();
  ~TrackSmearing();

  void Init();
  void Process();
  void Finish();

private:
  Double_t ptError(const Double_t, const Double_t, const Double_t, const Double_t);

  Double_t fBz;

  DelphesFormula *fD0Formula = nullptr; //!
  std::string fD0ResolutionFile;
  std::string fD0ResolutionHist;
  Bool_t fUseD0Formula;

  DelphesFormula *fDZFormula = nullptr; //!
  std::string fDZResolutionFile;
  std::string fDZResolutionHist;
  Bool_t fUseDZFormula;

  DelphesFormula *fPFormula = nullptr; //!
  std::string fPResolutionFile;
  std::string fPResolutionHist;
  Bool_t fUsePFormula;

  DelphesFormula *fCtgThetaFormula = nullptr; //!
  std::string fCtgThetaResolutionFile;
  std::string fCtgThetaResolutionHist;
  Bool_t fUseCtgThetaFormula;

  DelphesFormula *fPhiFormula = nullptr; //!
  std::string fPhiResolutionFile;
  std::string fPhiResolutionHist;
  Bool_t fUsePhiFormula;

  Bool_t fApplyToPileUp;

  TIterator *fItInputArray = nullptr; //!

  const TObjArray *fInputArray = nullptr; //!
  const TObjArray *fBeamSpotInputArray = nullptr; //!

  TObjArray *fOutputArray = nullptr; //!

  ClassDef(TrackSmearing, 1)
};

#endif
