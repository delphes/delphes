/** \class OldCalorimeter
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

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesModule.h"

#include <TLorentzVector.h>
#include <TRandom3.h>

#include <set>

using namespace std;

class OldCalorimeter: public DelphesModule
{
public:
  explicit OldCalorimeter(const DelphesParameters &moduleParams) :
    DelphesModule(moduleParams),
    fPhiBins(Steer<std::unordered_map<double, std::vector<double> > >("EtaPhiBins")),
    fFractionMap(Steer<TFractionMap>("EnergyFraction")), // read energy fractions for different particles
    fECalResolutionFormula(std::make_unique<DelphesFormula>()),
    fHCalResolutionFormula(std::make_unique<DelphesFormula>()),
    fTowerECalArray(std::make_shared<std::vector<Candidate *> >()),
    fTowerHCalArray(std::make_shared<std::vector<Candidate *> >()),
    fTowerTrackArray(std::make_shared<std::vector<Candidate *> >()),
    fTowerECalTrackArray(std::make_shared<std::vector<Candidate *> >()),
    fTowerHCalTrackArray(std::make_shared<std::vector<Candidate *> >())
  {
    for(const std::pair<double, std::vector<double> > etaPhiBins : fPhiBins) // auto would avoid a copy
      fEtaBins.emplace_back(etaPhiBins.first);
    std::sort(fEtaBins.begin(), fEtaBins.end());

    // set default energy fractions values
    fFractionMap[0] = std::make_pair(0., 1.);
    /*TFractionMap::iterator itFractionMap;
    for(itFractionMap = fFractionMap.begin(); itFractionMap != fFractionMap.end(); ++itFractionMap)
      cout << itFractionMap->first << "   " << itFractionMap->second.first << "   " << itFractionMap->second.second << endl;*/

    // read resolution formulas
    fECalResolutionFormula->Compile(Steer<std::string>("ECalResolutionFormula", "0"));
    fHCalResolutionFormula->Compile(Steer<std::string>("HCalResolutionFormula", "0"));
  }

  void Init() override
  {
    fParticleInputArray = ImportArray(Steer<std::string>("ParticleInputArray", "ParticlePropagator/particles"));
    fTrackInputArray = ImportArray(Steer<std::string>("TrackInputArray", "ParticlePropagator/tracks"));
    fTowerOutputArray = ExportArray(Steer<std::string>("TowerOutputArray", "towers"));
    fPhotonOutputArray = ExportArray(Steer<std::string>("PhotonOutputArray", "photons"));
    fEFlowTrackOutputArray = ExportArray(Steer<std::string>("EFlowTrackOutputArray", "eflowTracks"));
    fEFlowTowerOutputArray = ExportArray(Steer<std::string>("EFlowTowerOutputArray", "eflowTowers"));
  }
  void Process() override;

private:
  typedef std::map<unsigned long long, std::pair<double, double> > TFractionMap; //!
  typedef std::map<double, std::set<double> > TBinMap; //!

  void FinalizeTower();
  double LogNormal(double mean, double sigma);

  const std::unordered_map<double, std::vector<double> > fPhiBins;
  std::vector<double> fEtaBins;

  TFractionMap fFractionMap; //!

  const std::unique_ptr<DelphesFormula> fECalResolutionFormula; //!
  const std::unique_ptr<DelphesFormula> fHCalResolutionFormula; //!

  Candidate *fTower{nullptr};
  double fTowerEta, fTowerPhi, fTowerEdges[4];
  double fTowerECalEnergy, fTowerHCalEnergy;
  double fTowerECalNeutralEnergy, fTowerHCalNeutralEnergy;
  int fTowerPhotonHits, fTowerECalHits, fTowerHCalHits, fTowerAllHits;
  int fTowerECalTrackHits, fTowerHCalTrackHits, fTowerTrackAllHits;

  TBinMap fBinMap; //!

  std::vector<unsigned long long> fTowerHits;

  std::vector<double> fECalFractions;
  std::vector<double> fHCalFractions;

  // unsaved output collections
  const CandidatesCollection fTowerECalArray; //!
  const CandidatesCollection fTowerHCalArray; //!
  const CandidatesCollection fTowerTrackArray; //!
  const CandidatesCollection fTowerECalTrackArray; //!
  const CandidatesCollection fTowerHCalTrackArray; //!

  // input collections
  CandidatesCollection fParticleInputArray; //!
  CandidatesCollection fTrackInputArray; //!

  // saved output collections
  CandidatesCollection fTowerOutputArray; //!
  CandidatesCollection fPhotonOutputArray; //!
  CandidatesCollection fEFlowTrackOutputArray; //!
  CandidatesCollection fEFlowTowerOutputArray; //!
};

//------------------------------------------------------------------------------

void OldCalorimeter::Process()
{
  fTowerOutputArray->clear();
  fPhotonOutputArray->clear();
  fEFlowTrackOutputArray->clear();
  fEFlowTowerOutputArray->clear();

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalFractions.clear();
  fHCalFractions.clear();

  // loop over all particles
  size_t number = 0;
  for(Candidate *const &particle : *fParticleInputArray)
  {
    const TLorentzVector &particlePosition = particle->Position;

    const int pdgCode = std::abs(particle->PID);

    TFractionMap::iterator itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
      itFractionMap = fFractionMap.find(0);

    const double ecalFraction = itFractionMap->second.first;
    const double hcalFraction = itFractionMap->second.second;

    fECalFractions.push_back(ecalFraction);
    fHCalFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), particlePosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    const short phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    short flags = 0;
    flags |= (ecalFraction >= 1.0E-9) << 1;
    flags |= (hcalFraction >= 1.0E-9) << 2;
    flags |= (pdgCode == 11 || pdgCode == 22) << 3;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    const unsigned long long towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);
    ++number;
  }

  // loop over all tracks
  number = 0;
  for(Candidate *const &track : *fTrackInputArray)
  {
    const TLorentzVector &trackPosition = track->Position;

    const int pdgCode = std::abs(track->PID);

    TFractionMap::iterator itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
      itFractionMap = fFractionMap.find(0);

    const double ecalFraction = itFractionMap->second.first;
    const double hcalFraction = itFractionMap->second.second;

    // find eta bin [1, fEtaBins.size - 1]
    std::vector<double>::const_iterator itEtaBin = std::lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    const short etaBin = std::distance(fEtaBins.cbegin(), itEtaBin);

    // phi bins for given eta bin
    const std::vector<double> &phiBins = fPhiBins.at(*itEtaBin);

    // find phi bin [1, phiBins.size - 1]
    std::vector<double>::const_iterator itPhiBin = std::lower_bound(phiBins.begin(), phiBins.end(), trackPosition.Phi());
    if(itPhiBin == phiBins.begin() || itPhiBin == phiBins.end()) continue;
    const short phiBin = std::distance(phiBins.cbegin(), itPhiBin);

    short flags = 1;
    flags |= (ecalFraction >= 1.0E-9) << 1;
    flags |= (hcalFraction >= 1.0E-9) << 2;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    const unsigned long long towerHit = ((unsigned long long)(etaBin) << 48) | ((unsigned long long)(phiBin) << 32) | ((unsigned long long)(flags) << 24) | (unsigned long long)(number);

    fTowerHits.push_back(towerHit);
    ++number;
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  unsigned long long towerEtaPhi = 0;
  fTower = 0;
  for(const unsigned long long &towerHit : fTowerHits)
  {
    const short flags = (towerHit >> 24) & 0x00000000000000FFLL;
    const int number = (towerHit) & 0x0000000000FFFFFFLL;
    const unsigned long long hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();

      const short phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      const short etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      // phi bins for given eta bin
      const std::vector<double> &phiBins = fPhiBins.at(fEtaBins.at(etaBin));

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5 * (fEtaBins.at(etaBin - 1) + fEtaBins.at(etaBin));
      fTowerPhi = 0.5 * (phiBins.at(phiBin - 1) + phiBins.at(phiBin));

      fTowerEdges[0] = fEtaBins.at(etaBin - 1);
      fTowerEdges[1] = fEtaBins.at(etaBin);
      fTowerEdges[2] = phiBins.at(phiBin - 1);
      fTowerEdges[3] = phiBins.at(phiBin);

      fTowerECalEnergy = 0.0;
      fTowerHCalEnergy = 0.0;

      fTowerPhotonHits = 0;

      fTowerAllHits = 0;
      fTowerECalHits = 0;
      fTowerHCalHits = 0;

      fTowerTrackAllHits = 0;
      fTowerECalTrackHits = 0;
      fTowerHCalTrackHits = 0;

      fTowerECalArray->clear();
      fTowerHCalArray->clear();

      fTowerTrackArray->clear();
      fTowerECalTrackArray->clear();
      fTowerHCalTrackArray->clear();
    }

    // check for track hits
    if(flags & 1)
    {
      Candidate *track = static_cast<Candidate *>(fTrackInputArray->at(number));

      ++fTowerTrackAllHits;
      fTowerTrackArray->emplace_back(track);

      // check for track ECAL hits
      if(flags & 2)
      {
        ++fTowerECalTrackHits;
        fTowerECalTrackArray->emplace_back(track);
      }

      // check for track HCAL hits
      if(flags & 4)
      {
        ++fTowerHCalTrackHits;
        fTowerHCalTrackArray->emplace_back(track);
      }
      continue;
    }

    ++fTowerAllHits;

    // check for ECAL hits
    if(flags & 2)
    {
      Candidate *particle = nullptr; //FIXME
      ++fTowerECalHits;
      fTowerECalArray->emplace_back(particle);
    }

    // check for HCAL hits
    if(flags & 4)
    {
      Candidate *particle = nullptr; //FIXME
      ++fTowerHCalHits;
      fTowerHCalArray->emplace_back(particle);
    }

    // check for photon and electron hits in current tower
    if(flags & 8) ++fTowerPhotonHits;

    Candidate *particle = static_cast<Candidate *>(fParticleInputArray->at(number));
    TLorentzVector &momentum = particle->Momentum;

    // fill current tower
    const double ecalEnergy = momentum.E() * fECalFractions[number];
    const double hcalEnergy = momentum.E() * fHCalFractions[number];

    fTowerECalEnergy += ecalEnergy;
    fTowerHCalEnergy += hcalEnergy;

    fTower->AddCandidate(particle);
  }

  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void OldCalorimeter::FinalizeTower()
{
  const CandidatesCollection *towerTrackArray = nullptr;

  if(!fTower) return;

  //  ecalEnergy = gRandom->Gaus(fTowerECalEnergy, fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy));
  //  if(ecalEnergy < 0.0) ecalEnergy = 0.0;

  const double ecalEnergy = LogNormal(fTowerECalEnergy, fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy));

  //  hcalEnergy = gRandom->Gaus(fTowerHCalEnergy, fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy));
  //  if(hcalEnergy < 0.0) hcalEnergy = 0.0;

  const double hcalEnergy = LogNormal(fTowerHCalEnergy, fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy));

  const double energy = ecalEnergy + hcalEnergy;

  //  eta = fTowerEta;
  //  phi = fTowerPhi;

  const double eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
  const double phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);

  const double pt = energy / std::cosh(eta);

  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.0);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
  fTower->Ehad = hcalEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  // fill calorimeter towers and photon candidates
  if(energy > 0.0)
  {
    if(fTowerPhotonHits > 0 && fTowerTrackAllHits == 0)
    {
      fPhotonOutputArray->emplace_back(fTower);
    }

    fTowerOutputArray->emplace_back(fTower);
  }

  // fill energy flow candidates
  if(fTowerTrackAllHits == fTowerAllHits)
  {
    for(Candidate *const &track : *fTowerTrackArray)
      fEFlowTrackOutputArray->emplace_back(track);
  }
  else if(fTowerTrackAllHits > 0 && fTowerECalHits + fTowerHCalHits == fTowerAllHits)
  {
    if(fTowerECalHits == fTowerECalTrackHits && fTowerHCalHits == fTowerHCalTrackHits)
    {
      towerTrackArray = &fTowerTrackArray;
    }
    else if(fTowerECalHits == fTowerECalTrackHits)
    {
      towerTrackArray = &fTowerECalTrackArray;

      if(hcalEnergy > 0.0)
      {
        DelphesFactory *factory = GetFactory();

        // create new tower
        Candidate *tower = factory->NewCandidate();

        for(Candidate *const &particle : *fTowerHCalArray)
        {
          tower->AddCandidate(particle);
        }

        const double pt = hcalEnergy / std::cosh(eta);

        tower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.0);
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, hcalEnergy);
        tower->Eem = 0.0;
        tower->Ehad = hcalEnergy;

        tower->Edges[0] = fTowerEdges[0];
        tower->Edges[1] = fTowerEdges[1];
        tower->Edges[2] = fTowerEdges[2];
        tower->Edges[3] = fTowerEdges[3];

        fEFlowTowerOutputArray->emplace_back(tower);
      }
    }
    else if(fTowerHCalHits == fTowerHCalTrackHits)
    {
      towerTrackArray = &fTowerHCalTrackArray;

      if(ecalEnergy > 0.0)
      {
        DelphesFactory *factory = GetFactory();

        // create new tower
        Candidate *tower = factory->NewCandidate();

        for(Candidate *const &particle : *fTowerECalArray)
        {
          tower->AddCandidate(particle);
        }

        const double pt = ecalEnergy / std::cosh(eta);

        tower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.0);
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, ecalEnergy);
        tower->Eem = ecalEnergy;
        tower->Ehad = 0.0;

        tower->Edges[0] = fTowerEdges[0];
        tower->Edges[1] = fTowerEdges[1];
        tower->Edges[2] = fTowerEdges[2];
        tower->Edges[3] = fTowerEdges[3];

        fEFlowTowerOutputArray->emplace_back(tower);
      }
    }
    else
      fEFlowTowerOutputArray->emplace_back(fTower);

    if(towerTrackArray)
    {
      for(Candidate *const &track : **towerTrackArray)
        fEFlowTrackOutputArray->emplace_back(track);
    }
  }
  else if(energy > 0.0)
    fEFlowTowerOutputArray->emplace_back(fTower);
}

//------------------------------------------------------------------------------

double OldCalorimeter::LogNormal(double mean, double sigma)
{
  if(mean > 0.0)
  {
    const double b = std::sqrt(std::log((1.0 + (sigma * sigma) / (mean * mean)))),
                 a = std::log(mean) - 0.5 * b * b;
    return std::exp(a + b * gRandom->Gaus(0, 1));
  }
  else
    return 0.0;
}

//------------------------------------------------------------------------------

REGISTER_MODULE("OldCalorimeter", OldCalorimeter);
