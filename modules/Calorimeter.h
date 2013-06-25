#ifndef Calorimeter_h
#define Calorimeter_h

/** \class Calorimeter
 *
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  preselects towers hit by photons and creates energy flow objects.
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <map>
#include <set>
#include <vector>

class TObjArray;
class DelphesFormula;
class Candidate;

class Calorimeter: public DelphesModule
{
public:

  Calorimeter();
  ~Calorimeter();

  void Init();
  void Process();
  void Finish();
  Double_t Draw_lognormal(Double_t mu, Double_t sigma);
private:

  typedef std::map< Long64_t, std::pair< Double_t, Double_t > > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fTowerECalEnergy, fTowerHCalEnergy;
  Double_t fTowerECalNeutralEnergy, fTowerHCalNeutralEnergy;
  Int_t fTowerNeutralHits, fTowerPhotonHits, fTowerElectronHits, fTowerTrackHits, fTowerAllHits;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fECalFractions;
  std::vector < Double_t > fHCalFractions;

  DelphesFormula *fECalResolutionFormula; //!
  DelphesFormula *fHCalResolutionFormula; //!

  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fTowerOutputArray; //!
  TObjArray *fPhotonOutputArray; //!

  TObjArray *fEFlowTrackOutputArray; //!
  TObjArray *fEFlowTowerOutputArray; //!

  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  TObjArray *fTowerPhotonArray; //!
  TIterator *fItTowerPhotonArray; //!

  void FinalizeTower();
 
  ClassDef(Calorimeter, 1)
};

#endif
