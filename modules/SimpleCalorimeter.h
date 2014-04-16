#ifndef SimpleCalorimeter_h
#define SimpleCalorimeter_h

/** \class SimpleCalorimeter
 *
 *  Fills SimpleCalorimeter towers, performs SimpleCalorimeter resolution smearing,
 *  and creates energy flow objects (tracks, photons, and neutral hadrons).
 *
 *  $Date: 2014-04-16 15:29:31 +0200 (Wed, 16 Apr 2014) $
 *  $Revision: 1364 $
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

class SimpleCalorimeter: public DelphesModule
{
public:

  SimpleCalorimeter();
  ~SimpleCalorimeter();

  void Init();
  void Process();
  void Finish();

private:

  typedef std::map< Long64_t, Double_t > TFractionMap; //!
  typedef std::map< Double_t, std::set< Double_t > > TBinMap; //!

  Candidate *fTower;
  Double_t fTowerEta, fTowerPhi, fTowerEdges[4];
  Double_t fTowerEnergy;
  Double_t fTrackEnergy;
  
  Double_t fTowerTime;
  Double_t fTrackTime;
   
  Double_t fTowerWeightTime;
  Double_t fTrackWeightTime;
  
  Int_t fTowerTrackHits, fTowerPhotonHits;

  TFractionMap fFractionMap; //!
  TBinMap fBinMap; //!

  std::vector < Double_t > fEtaBins;
  std::vector < std::vector < Double_t >* > fPhiBins;

  std::vector < Long64_t > fTowerHits;

  std::vector < Double_t > fTowerFractions;
  
  std::vector < Double_t > fTrackFractions;
 
  DelphesFormula *fResolutionFormula; //!
 
  TIterator *fItParticleInputArray; //!
  TIterator *fItTrackInputArray; //!

  const TObjArray *fParticleInputArray; //!
  const TObjArray *fTrackInputArray; //!

  TObjArray *fTowerOutputArray; //!
 
  TObjArray *fEFlowTowerOutputArray; //!
 
  TObjArray *fTowerTrackArray; //!
  TIterator *fItTowerTrackArray; //!

  void FinalizeTower();
  Double_t LogNormal(Double_t mean, Double_t sigma);

  ClassDef(SimpleCalorimeter, 1)
};

#endif
