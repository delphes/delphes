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

  DelphesFormula *fD0Formula; //!
  std::string fD0ResolutionFile;
  std::string fD0ResolutionHist;
  Bool_t fUseD0Formula;

  DelphesFormula *fDZFormula; //!
  std::string fDZResolutionFile;
  std::string fDZResolutionHist;
  Bool_t fUseDZFormula;

  DelphesFormula *fPFormula; //!
  std::string fPResolutionFile;
  std::string fPResolutionHist;
  Bool_t fUsePFormula;

  DelphesFormula *fCtgThetaFormula; //!
  std::string fCtgThetaResolutionFile;
  std::string fCtgThetaResolutionHist;
  Bool_t fUseCtgThetaFormula;

  DelphesFormula *fPhiFormula; //!
  std::string fPhiResolutionFile;
  std::string fPhiResolutionHist;
  Bool_t fUsePhiFormula;

  Bool_t fApplyToPileUp;

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!
  const TObjArray *fBeamSpotInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(TrackSmearing, 1)
};

#endif
