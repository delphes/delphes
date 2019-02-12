
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

#include "modules/OldCalorimeter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

OldCalorimeter::OldCalorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerECalArray(0), fItTowerECalArray(0),
  fTowerHCalArray(0), fItTowerHCalArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0),
  fTowerECalTrackArray(0), fItTowerECalTrackArray(0),
  fTowerHCalTrackArray(0), fItTowerHCalTrackArray(0)
{
  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  fTowerECalArray = new TObjArray;
  fItTowerECalArray = fTowerECalArray->MakeIterator();
  fTowerHCalArray = new TObjArray;
  fItTowerHCalArray = fTowerHCalArray->MakeIterator();

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();
  fTowerECalTrackArray = new TObjArray;
  fItTowerECalTrackArray = fTowerECalTrackArray->MakeIterator();
  fTowerHCalTrackArray = new TObjArray;
  fItTowerHCalTrackArray = fTowerHCalTrackArray->MakeIterator();
}

//------------------------------------------------------------------------------

OldCalorimeter::~OldCalorimeter()
{
  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;

  if(fTowerECalArray) delete fTowerECalArray;
  if(fItTowerECalArray) delete fItTowerECalArray;
  if(fTowerHCalArray) delete fTowerHCalArray;
  if(fItTowerHCalArray) delete fItTowerHCalArray;

  if(fTowerTrackArray) delete fTowerTrackArray;
  if(fItTowerTrackArray) delete fItTowerTrackArray;
  if(fTowerECalTrackArray) delete fTowerECalTrackArray;
  if(fItTowerECalTrackArray) delete fItTowerECalTrackArray;
  if(fTowerHCalTrackArray) delete fTowerHCalTrackArray;
  if(fItTowerHCalTrackArray) delete fItTowerHCalTrackArray;
}

//------------------------------------------------------------------------------

void OldCalorimeter::Init()
{
  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins, sizeFractions;
  Double_t ecalFraction, hcalFraction;
  TBinMap::iterator itEtaBin;
  set<Double_t>::iterator itPhiBin;
  vector<Double_t> *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size / 2; ++i)
  {
    paramEtaBins = param[i * 2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i * 2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      for(k = 0; k < sizePhiBins; ++k)
      {
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin)
  {
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector<double>(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin)
    {
      phiBins->push_back(*itPhiBin);
    }
  }

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param.GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  fFractionMap[0] = make_pair(0.0, 1.0);

  for(i = 0; i < size / 2; ++i)
  {
    paramFractions = param[i * 2 + 1];
    sizeFractions = paramFractions.GetSize();

    ecalFraction = paramFractions[0].GetDouble();
    hcalFraction = paramFractions[1].GetDouble();

    fFractionMap[param[i * 2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }
  /*
  TFractionMap::iterator itFractionMap;
  for(itFractionMap = fFractionMap.begin(); itFractionMap != fFractionMap.end(); ++itFractionMap)
  {
    cout << itFractionMap->first << "   " << itFractionMap->second.first  << "   " << itFractionMap->second.second << endl;
  }
*/
  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  fTowerOutputArray = ExportArray(GetString("TowerOutputArray", "towers"));
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));

  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "eflowTowers"));
}

//------------------------------------------------------------------------------

void OldCalorimeter::Finish()
{
  vector<vector<Double_t> *>::iterator itPhiBin;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  for(itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin)
  {
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void OldCalorimeter::Process()
{
  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector<Double_t>::iterator itEtaBin;
  vector<Double_t>::iterator itPhiBin;
  vector<Double_t> *phiBins;

  vector<Long64_t>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fECalFractions.clear();
  fHCalFractions.clear();

  // loop over all particles
  fItParticleInputArray->Reset();
  number = -1;
  while((particle = static_cast<Candidate *>(fItParticleInputArray->Next())))
  {
    const TLorentzVector &particlePosition = particle->Position;
    ++number;

    pdgCode = TMath::Abs(particle->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fECalFractions.push_back(ecalFraction);
    fHCalFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), particlePosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 0;
    flags |= (ecalFraction >= 1.0E-9) << 1;
    flags |= (hcalFraction >= 1.0E-9) << 2;
    flags |= (pdgCode == 11 || pdgCode == 22) << 3;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // loop over all tracks
  fItTrackInputArray->Reset();
  number = -1;
  while((track = static_cast<Candidate *>(fItTrackInputArray->Next())))
  {
    const TLorentzVector &trackPosition = track->Position;
    ++number;

    pdgCode = TMath::Abs(track->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end())
    {
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1;
    flags |= (ecalFraction >= 1.0E-9) << 1;
    flags |= (hcalFraction >= 1.0E-9) << 2;

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end());

  // loop over all hits
  towerEtaPhi = 0;
  fTower = 0;
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits)
  {
    towerHit = (*itTowerHits);
    flags = (towerHit >> 24) & 0x00000000000000FFLL;
    number = (towerHit)&0x0000000000FFFFFFLL;
    hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi)
    {
      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower
      fTower = factory->NewCandidate();

      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5 * (fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5 * ((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fTowerECalEnergy = 0.0;
      fTowerHCalEnergy = 0.0;

      fTowerPhotonHits = 0;

      fTowerAllHits = 0;
      fTowerECalHits = 0;
      fTowerHCalHits = 0;

      fTowerTrackAllHits = 0;
      fTowerECalTrackHits = 0;
      fTowerHCalTrackHits = 0;

      fTowerECalArray->Clear();
      fTowerHCalArray->Clear();

      fTowerTrackArray->Clear();
      fTowerECalTrackArray->Clear();
      fTowerHCalTrackArray->Clear();
    }

    // check for track hits
    if(flags & 1)
    {
      track = static_cast<Candidate *>(fTrackInputArray->At(number));

      ++fTowerTrackAllHits;
      fTowerTrackArray->Add(track);

      // check for track ECAL hits
      if(flags & 2)
      {
        ++fTowerECalTrackHits;
        fTowerECalTrackArray->Add(track);
      }

      // check for track HCAL hits
      if(flags & 4)
      {
        ++fTowerHCalTrackHits;
        fTowerHCalTrackArray->Add(track);
      }
      continue;
    }

    ++fTowerAllHits;

    // check for ECAL hits
    if(flags & 2)
    {
      ++fTowerECalHits;
      fTowerECalArray->Add(particle);
    }

    // check for HCAL hits
    if(flags & 4)
    {
      ++fTowerHCalHits;
      fTowerHCalArray->Add(particle);
    }

    // check for photon and electron hits in current tower
    if(flags & 8) ++fTowerPhotonHits;

    particle = static_cast<Candidate *>(fParticleInputArray->At(number));
    momentum = particle->Momentum;

    // fill current tower
    ecalEnergy = momentum.E() * fECalFractions[number];
    hcalEnergy = momentum.E() * fHCalFractions[number];

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
  Candidate *particle, *track, *tower;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  TIterator *itTowerTrackArray;

  if(!fTower) return;

  //  ecalEnergy = gRandom->Gaus(fTowerECalEnergy, fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy));
  //  if(ecalEnergy < 0.0) ecalEnergy = 0.0;

  ecalEnergy = LogNormal(fTowerECalEnergy, fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy));

  //  hcalEnergy = gRandom->Gaus(fTowerHCalEnergy, fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy));
  //  if(hcalEnergy < 0.0) hcalEnergy = 0.0;

  hcalEnergy = LogNormal(fTowerHCalEnergy, fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy));

  energy = ecalEnergy + hcalEnergy;

  //  eta = fTowerEta;
  //  phi = fTowerPhi;

  eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
  phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);

  pt = energy / TMath::CosH(eta);

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
      fPhotonOutputArray->Add(fTower);
    }

    fTowerOutputArray->Add(fTower);
  }

  // fill energy flow candidates
  if(fTowerTrackAllHits == fTowerAllHits)
  {
    fItTowerTrackArray->Reset();
    while((track = static_cast<Candidate *>(fItTowerTrackArray->Next())))
    {
      fEFlowTrackOutputArray->Add(track);
    }
  }
  else if(fTowerTrackAllHits > 0 && fTowerECalHits + fTowerHCalHits == fTowerAllHits)
  {
    if(fTowerECalHits == fTowerECalTrackHits && fTowerHCalHits == fTowerHCalTrackHits)
    {
      itTowerTrackArray = fItTowerTrackArray;
    }
    else if(fTowerECalHits == fTowerECalTrackHits)
    {
      itTowerTrackArray = fItTowerECalTrackArray;

      if(hcalEnergy > 0.0)
      {
        DelphesFactory *factory = GetFactory();

        // create new tower
        tower = factory->NewCandidate();

        fItTowerHCalArray->Reset();
        while((particle = static_cast<Candidate *>(fItTowerHCalArray->Next())))
        {
          tower->AddCandidate(particle);
        }

        pt = hcalEnergy / TMath::CosH(eta);

        tower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.0);
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, hcalEnergy);
        tower->Eem = 0.0;
        tower->Ehad = hcalEnergy;

        tower->Edges[0] = fTowerEdges[0];
        tower->Edges[1] = fTowerEdges[1];
        tower->Edges[2] = fTowerEdges[2];
        tower->Edges[3] = fTowerEdges[3];

        fEFlowTowerOutputArray->Add(tower);
      }
    }
    else if(fTowerHCalHits == fTowerHCalTrackHits)
    {
      itTowerTrackArray = fItTowerHCalTrackArray;

      if(ecalEnergy > 0.0)
      {
        DelphesFactory *factory = GetFactory();

        // create new tower
        tower = factory->NewCandidate();

        fItTowerECalArray->Reset();
        while((particle = static_cast<Candidate *>(fItTowerECalArray->Next())))
        {
          tower->AddCandidate(particle);
        }

        pt = ecalEnergy / TMath::CosH(eta);

        tower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.0);
        tower->Momentum.SetPtEtaPhiE(pt, eta, phi, ecalEnergy);
        tower->Eem = ecalEnergy;
        tower->Ehad = 0.0;

        tower->Edges[0] = fTowerEdges[0];
        tower->Edges[1] = fTowerEdges[1];
        tower->Edges[2] = fTowerEdges[2];
        tower->Edges[3] = fTowerEdges[3];

        fEFlowTowerOutputArray->Add(tower);
      }
    }
    else
    {
      itTowerTrackArray = 0;
      fEFlowTowerOutputArray->Add(fTower);
    }

    if(itTowerTrackArray)
    {
      itTowerTrackArray->Reset();
      while((track = static_cast<Candidate *>(itTowerTrackArray->Next())))
      {
        fEFlowTrackOutputArray->Add(track);
      }
    }
  }
  else if(energy > 0.0)
  {
    fEFlowTowerOutputArray->Add(fTower);
  }
}

//------------------------------------------------------------------------------

Double_t OldCalorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma * sigma) / (mean * mean))));
    a = TMath::Log(mean) - 0.5 * b * b;

    return TMath::Exp(a + b * gRandom->Gaus(0, 1));
  }
  else
  {
    return 0.0;
  }
}
