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

class Tower: public SortableObject
{
public:
  Tower() = default;
  explicit Tower(const Candidate &);

  Float_t ET{0.f}; // calorimeter tower transverse energy
  Float_t Eta{0.f}; // calorimeter tower pseudorapidity
  Float_t Phi{0.f}; // calorimeter tower azimuthal angle

  Float_t E{0.f}; // calorimeter tower energy

  Float_t T{0.f}; // calo deposit time, averaged by sqrt(EM energy) over all particles
  Float_t X{0.f}; // calo tower position
  Float_t Y{0.f}; // calo tower position
  Float_t Z{0.f}; // calo tower position

  Int_t NTimeHits{0}; // number of hits contributing to time measurement

  Float_t Eem{0.f}; // calorimeter tower electromagnetic energy
  Float_t Ehad{0.f}; // calorimeter tower hadronic energy
  Float_t Etrk{0.f}; // total charged energy hitting tower

  Float_t Edges[4]{0.f}; // calorimeter tower edges

  Int_t IsPU{0}; // 0 or 1 for particles from pile-up interactions
  Int_t IsRecoPU{0}; // 0 or 1 for reconstructed particles from pile-up
  Float_t HardEnergyFraction{0.f}; // fraction of hard scattering vs PU energy

  //TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(Tower, 5)
};

class Track: public SortableObject
{
public:
  Track() = default;
  explicit Track(const Candidate &);

  Int_t PID{0}; // HEP ID number

  Int_t Charge{0}; // track charge

  Int_t IsPU{0}; // 0 or 1 for particles from pile-up interactions
  Int_t IsRecoPU{0}; // 0 or 1 for reconstructed particles from pile-up
  Float_t HardEnergyFraction{0.f}; // fraction of hard scattering vs PU energy in the particle flow candidate

  Float_t P{0.f}; // track momentum
  Float_t PT{0.f}; // track transverse momentum
  Float_t Eta{0.f}; // track pseudorapidity
  Float_t Phi{0.f}; // track azimuthal angle
  Float_t CtgTheta{0.f}; // track cotangent of theta
  Float_t C{0.f}; // track curvature inverse
  Float_t Mass{0.f}; // particle mass

  Float_t EtaOuter{0.f}; // track pseudorapidity at the tracker edge
  Float_t PhiOuter{0.f}; // track azimuthal angle at the tracker edge

  Float_t T{0.f}; // track vertex position (t component)
  Float_t X{0.f}; // track vertex position (x component)
  Float_t Y{0.f}; // track vertex position (y component)
  Float_t Z{0.f}; // track vertex position (z component)

  Float_t TOuter{0.f}; // track position (t component) at the tracker edge
  Float_t XOuter{0.f}; // track position (x component) at the tracker edge
  Float_t YOuter{0.f}; // track position (y component) at the tracker edge
  Float_t ZOuter{0.f}; // track position (z component) at the tracker edge

  Float_t Xd{0.f}; // X coordinate of point of closest approach to vertex
  Float_t Yd{0.f}; // Y coordinate of point of closest approach to vertex
  Float_t Zd{0.f}; // Z coordinate of point of closest approach to vertex

  Float_t XFirstHit{0.f}; // X coordinate of point of closest approach to vertex
  Float_t YFirstHit{0.f}; // Y coordinate of point of closest approach to vertex
  Float_t ZFirstHit{0.f}; // Z coordinate of point of closest approach to vertex

  Float_t L{0.f}; // track path length
  Float_t D0{0.f}; // track transverse impact parameter
  Float_t DZ{0.f}; // track longitudinal impact parameter
  Float_t Nclusters{0.f}; // Number of ionization clusters
  Float_t dNdx{0.f}; // Number of ionization clusters

  Float_t ErrorP{0.f}; // track momentum error
  Float_t ErrorPT{0.f}; // track transverse momentum error
  Float_t ErrorPhi{0.f}; // track azimuthal angle error
  Float_t ErrorCtgTheta{0.f}; // track cotangent of theta error

  Float_t ErrorT{0.f}; // time measurement error
  Float_t ErrorD0{0.f}; // track transverse impact parameter error
  Float_t ErrorDZ{0.f}; // track longitudinal impact parameter error
  Float_t ErrorC{0.f}; // track curvature error

  // track covariance off-diagonal terms
  Float_t ErrorD0Phi{0.f};
  Float_t ErrorD0C{0.f};
  Float_t ErrorD0DZ{0.f};
  Float_t ErrorD0CtgTheta{0.f};
  Float_t ErrorPhiC{0.f};
  Float_t ErrorPhiDZ{0.f};
  Float_t ErrorPhiCtgTheta{0.f};
  Float_t ErrorCDZ{0.f};
  Float_t ErrorCCtgTheta{0.f};
  Float_t ErrorDZCtgTheta{0.f};

  //TRef Particle; // reference to generated particle

  Int_t VertexIndex{-1}; // reference to vertex

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;
  TMatrixDSym CovarianceMatrix() const;

  ClassDef(Track, 4)
};
} // namespace DelphesRNTuple

#endif // DelphesRNTupleClasses_h
