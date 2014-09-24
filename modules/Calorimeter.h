#ifndef Calorimeter_h
#define Calorimeter_h

/** \class Calorimeter
 *
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
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

private:

  typedef std::map< Long64_t, std::pair< Double_t, Double_t > > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fTowerECalEnergy, fTowerHCalEnergy;
  Double_t fTrackECalEnergy, fTrackHCalEnergy;
  
  Double_t fTowerECalTime, fTowerHCalTime;
  Double_t fTrackECalTime, fTrackHCalTime;
   
  Double_t fTowerECalWeightTime, fTowerHCalWeightTime;
  Double_t fTrackECalWeightTime, fTrackHCalWeightTime;
  
  Int_t fTowerTrackHits, fTowerPhotonHits;

  Double_t fEcalEnergyMin;
  Double_t fHcalEnergyMin;
  
  Double_t fEcalSigmaMin;
  Double_t fHcalSigmaMin;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fTowerECalFractions;
  std::vector < Double_t > fTowerHCalFractions;

  std::vector < Double_t > fTrackECalFractions;
  std::vector < Double_t > fTrackHCalFractions;

  DelphesFormula *fECalResolutionFormula; //!
  DelphesFormula *fHCalResolutionFormula; //!

  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fTowerOutputArray; //!
  TObjArray *fPhotonOutputArray; //!
 
  TObjArray *fEFlowTrackOutputArray; //!
  TObjArray *fEFlowPhotonOutputArray; //!
  TObjArray *fEFlowNeutralHadronOutputArray; //!

  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);

  ClassDef(Calorimeter, 1)
};

#endif
