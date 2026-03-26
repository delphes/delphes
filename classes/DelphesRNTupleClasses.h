/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DelphesRNTupleClasses_h
#define DelphesRNTupleClasses_h

/**
 *
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

// Dependencies (#includes)

#include <Math/LorentzVector.h>

#include "classes/SortableObject.h"

class Candidate;

//---------------------------------------------------------------------------

namespace DelphesRNTuple
{
class GenParticle: public ::GenParticle
{
public:
  explicit GenParticle(const Candidate &cand);
  ClassDef(GenParticle, 4)
};

class Photon: public SortableObject
{
public:
  Photon() = default;
  explicit Photon(const Candidate &);

  Float_t PT{0.f}; // photon transverse momentum
  Float_t Eta{0.f}; // photon pseudorapidity
  Float_t Phi{0.f}; // photon azimuthal angle

  Float_t E{0.f}; // photon energy

  Float_t T{0.f}; // particle arrival time of flight

  Float_t EhadOverEem{0.f}; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  //TRefArray Particles; // references to generated particles

  Float_t IsolationVar{0.f}; // isolation variable
  Float_t IsolationVarRhoCorr{0.f}; // isolation variable
  Float_t SumPtCharged{0.f}; // isolation variable
  Float_t SumPtNeutral{0.f}; // isolation variable
  Float_t SumPtChargedPU{0.f}; // isolation variable
  Float_t SumPt{0.f}; // isolation variable

  Int_t Status{0}; // 1: prompt, 2: non prompt, 3: fake

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(Photon, 4)
};

//---------------------------------------------------------------------------

class Electron: public SortableObject
{
public:
  Electron() = default;
  explicit Electron(const Candidate &);

  Float_t PT{0.f}; // electron transverse momentum
  Float_t Eta{0.f}; // electron pseudorapidity
  Float_t Phi{0.f}; // electron azimuthal angle

  Float_t T{0.f}; // particle arrival time of flight

  Int_t Charge{0}; // electron charge

  Float_t EhadOverEem{0.f}; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  //TRef Particle; // reference to generated particle

  Float_t IsolationVar{0.f}; // isolation variable
  Float_t IsolationVarRhoCorr{0.f}; // isolation variable
  Float_t SumPtCharged{0.f}; // isolation variable
  Float_t SumPtNeutral{0.f}; // isolation variable
  Float_t SumPtChargedPU{0.f}; // isolation variable
  Float_t SumPt{0.f}; // isolation variable

  Float_t D0{0.f}; // track transverse impact parameter
  Float_t DZ{0.f}; // track longitudinal impact parameter
  Float_t ErrorD0{0.f}; // track transverse impact parameter error
  Float_t ErrorDZ{0.f}; // track longitudinal impact parameter error

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(Electron, 4)
};

//---------------------------------------------------------------------------

class Muon: public SortableObject
{
public:
  Muon() = default;
  explicit Muon(const Candidate &);

  Float_t PT{0.f}; // muon transverse momentum
  Float_t Eta{0.f}; // muon pseudorapidity
  Float_t Phi{0.f}; // muon azimuthal angle

  Float_t T{0.f}; // particle arrival time of flight

  Int_t Charge; // muon charge

  //TRef Particle; // reference to generated particle

  Float_t IsolationVar{0.f}; // isolation variable
  Float_t IsolationVarRhoCorr{0.f}; // isolation variable
  Float_t SumPtCharged{0.f}; // isolation variable
  Float_t SumPtNeutral{0.f}; // isolation variable
  Float_t SumPtChargedPU{0.f}; // isolation variable
  Float_t SumPt{0.f}; // isolation variable

  Float_t D0{0.f}; // track transverse impact parameter
  Float_t DZ{0.f}; // track longitudinal impact parameter
  Float_t ErrorD0{0.f}; // track transverse impact parameter error
  Float_t ErrorDZ{0.f}; // track longitudinal impact parameter error

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(Muon, 4)
};

//---------------------------------------------------------------------------

class Jet: public SortableObject
{
public:
  Jet() = default;
  explicit Jet(const Candidate &);

  Float_t PT{0.f}; // jet transverse momentum
  Float_t Eta{0.f}; // jet pseudorapidity
  Float_t Phi{0.f}; // jet azimuthal angle

  Float_t T{0.f}; //particle arrival time of flight

  Float_t Mass{0.f}; // jet invariant mass

  Float_t DeltaEta{0.f}; // jet radius in pseudorapidity
  Float_t DeltaPhi{0.f}; // jet radius in azimuthal angle

  UInt_t Flavor; // jet flavor
  UInt_t FlavorAlgo; // jet flavor
  UInt_t FlavorPhys; // jet flavor
  UInt_t TauFlavor; // jet flavor according to Tau tagging module

  UInt_t BTag; // 0 or 1 for a jet that has been tagged as containing a heavy quark
  UInt_t BTagAlgo; // 0 or 1 for a jet that has been tagged as containing a heavy quark
  UInt_t BTagPhys; // 0 or 1 for a jet that has been tagged as containing a heavy quark

  UInt_t TauTag; // 0 or 1 for a jet that has been tagged as a tau
  Float_t TauWeight{0.f}; // probability for jet to be identified as tau

  Int_t Charge; // tau charge

  Float_t EhadOverEem{0.f}; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

  Int_t NCharged; // number of charged constituents
  Int_t NNeutrals; // number of neutral constituents

  Float_t NeutralEnergyFraction{0.f}; // charged energy fraction
  Float_t ChargedEnergyFraction{0.f}; // neutral energy fraction

  Float_t Beta{0.f}; // (sum pt of charged pile-up constituents)/(sum pt of charged constituents)
  Float_t BetaStar{0.f}; // (sum pt of charged constituents coming from hard interaction)/(sum pt of charged constituents)
  Float_t MeanSqDeltaR{0.f}; // average distance (squared) between constituent and jet weighted by pt (squared) of constituent
  Float_t PTD{0.f}; // average pt between constituent and jet weighted by pt of constituent
  Float_t FracPt[5]; // (sum pt of constituents within a ring 0.1*i < DeltaR < 0.1*(i+1))/(sum pt of constituents)

  Float_t Tau[5]; // N-subjettiness

  ROOT::Math::XYZTVector SoftDroppedJet;
  ROOT::Math::XYZTVector SoftDroppedSubJet1;
  ROOT::Math::XYZTVector SoftDroppedSubJet2;

  ROOT::Math::XYZTVector TrimmedP4[5]; // first entry (i = 0) is the total Trimmed Jet 4-momenta and from i = 1 to 4 are the trimmed subjets 4-momenta
  ROOT::Math::XYZTVector PrunedP4[5]; // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta
  ROOT::Math::XYZTVector SoftDroppedP4[5]; // first entry (i = 0) is the total SoftDropped Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta

  Int_t NSubJetsTrimmed{0}; // number of subjets trimmed
  Int_t NSubJetsPruned{0}; // number of subjets pruned
  Int_t NSubJetsSoftDropped{0}; // number of subjets soft-dropped

  Double_t ExclYmerge12;
  Double_t ExclYmerge23;
  Double_t ExclYmerge34;
  Double_t ExclYmerge45;
  Double_t ExclYmerge56;

  //std::vector<CandidatePtr> Constituents; // references to constituents
  //std::vector<CandidatePtr> Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;
  ROOT::Math::XYZTVector Area;

  ClassDef(Jet, 5)
};
} // namespace DelphesRNTuple

#endif // DelphesRNTupleClasses_h
