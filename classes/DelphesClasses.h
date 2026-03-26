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

#ifndef DelphesClasses_h
#define DelphesClasses_h

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

#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TObject.h"
#include "TRef.h"
#include "TRefArray.h"

#include "classes/SortableObject.h"

class DelphesFactory;

class Candidate;

//---------------------------------------------------------------------------

class Event: public TObject
{
public:
  Long64_t Number{0}; // event number

  Float_t ReadTime{0.f}; // read time
  Float_t ProcTime{0.f}; // processing time

  ClassDef(Event, 1)
};

//---------------------------------------------------------------------------

class LHCOEvent: public Event
{
public:
  Int_t Trigger{0}; // trigger word

  ClassDef(LHCOEvent, 1)
};

//---------------------------------------------------------------------------

class LHEFEvent: public Event
{
public:
  Int_t ProcessID{0}; // subprocess code for the event | hepup.IDPRUP

  Float_t Weight{0.f}; // weight for the event | hepup.XWGTUP
  Float_t CrossSection{0.f}; // cross-section (read from init, implemented only for Wizard evgen)
  Float_t ScalePDF{0.f}; // scale in GeV used in the calculation of the PDFs in the event | hepup.SCALUP
  Float_t AlphaQED{0.f}; // value of the QED coupling used in the event | hepup.AQEDUP
  Float_t AlphaQCD{0.f}; // value of the QCD coupling used in the event | hepup.AQCDUP

  ClassDef(LHEFEvent, 3)
};

//---------------------------------------------------------------------------

class LHEFWeight: public TObject
{
public:
  Int_t ID{0}; // weight ID
  Float_t Weight{0.f}; // weight value

  ClassDef(LHEFWeight, 1)
};

//---------------------------------------------------------------------------

class HepMCEvent: public Event
{
public:
  Int_t ProcessID; // unique signal process id | signal_process_id()
  Int_t MPI; // number of multi parton interactions | mpi()

  Float_t Weight{0.f}; // weight for the event
  Float_t CrossSection{0.f}; // cross-section in pb
  Float_t CrossSectionError{0.f}; // cross-section error in pb

  Float_t Scale{0.f}; // energy scale, see hep-ph/0109068 | event_scale()
  Float_t AlphaQED{0.f}; // QED coupling, see hep-ph/0109068 | alphaQED()
  Float_t AlphaQCD{0.f}; // QCD coupling, see hep-ph/0109068 | alphaQCD()

  Int_t ID1; // flavour code of first parton | pdf_info()->id1()
  Int_t ID2; // flavour code of second parton | pdf_info()->id2()

  Float_t X1{0.f}; // fraction of beam momentum carried by first parton ("beam side") | pdf_info()->x1()
  Float_t X2{0.f}; // fraction of beam momentum carried by second parton ("target side") | pdf_info()->x2()

  Float_t ScalePDF{0.f}; // Q-scale used in evaluation of PDF's (in GeV) | pdf_info()->scalePDF()

  Float_t PDF1{0.f}; // PDF (id1, x1, Q) | pdf_info()->pdf1()
  Float_t PDF2{0.f}; // PDF (id2, x2, Q) | pdf_info()->pdf2()

  ClassDef(HepMCEvent, 3)
};

//---------------------------------------------------------------------------

class GenParticle: public SortableObject
{
public:
  GenParticle() = default;
  explicit GenParticle(const Candidate &);

  Int_t PID; // particle HEP ID number | hepevt.idhep[number]

  Int_t Status; // particle status | hepevt.isthep[number]
  Int_t IsPU; // 0 or 1 for particles from pile-up interactions

  Int_t M1; // particle 1st mother | hepevt.jmohep[number][0] - 1
  Int_t M2; // particle 2nd mother | hepevt.jmohep[number][1] - 1

  Int_t D1; // particle 1st daughter | hepevt.jdahep[number][0] - 1
  Int_t D2; // particle last daughter | hepevt.jdahep[number][1] - 1

  Int_t Charge; // particle charge

  Float_t Mass{0.f}; // particle mass

  Float_t E{0.f}; // particle energy | hepevt.phep[number][3]
  Float_t Px{0.f}; // particle momentum vector (x component) | hepevt.phep[number][0]
  Float_t Py{0.f}; // particle momentum vector (y component) | hepevt.phep[number][1]
  Float_t Pz{0.f}; // particle momentum vector (z component) | hepevt.phep[number][2]

  Float_t P{0.f}; // particle momentum
  Float_t PT{0.f}; // particle transverse momentum
  Float_t Eta{0.f}; // particle pseudorapidity
  Float_t Phi{0.f}; // particle azimuthal angle
  Float_t Rapidity{0.f}; // particle rapidity

  Float_t T{0.f}; // particle vertex position (t component) | hepevt.vhep[number][3]
  Float_t X{0.f}; // particle vertex position (x component) | hepevt.vhep[number][0]
  Float_t Y{0.f}; // particle vertex position (y component) | hepevt.vhep[number][1]
  Float_t Z{0.f}; // particle vertex position (z component) | hepevt.vhep[number][2]

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(GenParticle, 2)
};

//---------------------------------------------------------------------------

class Vertex: public SortableObject
{
public:
  Vertex() = default;
  explicit Vertex(const Candidate &);

  Float_t T{0.f}; // vertex position (t component)
  Float_t X{0.f}; // vertex position (x component)
  Float_t Y{0.f}; // vertex position (y component)
  Float_t Z{0.f}; // vertex position (z component)

  Double_t ErrorT{0.}; // vertex position error (t component)
  Double_t ErrorX{0.}; // vertex position error (x component)
  Double_t ErrorY{0.}; // vertex position error (y component)
  Double_t ErrorZ{0.}; // vertex position error (z component)

  Int_t Index{-1}; // vertex index
  Int_t NDF{0}; // number of degrees of freedom

  Double_t Sigma{0.}; // vertex position (z component) error
  Double_t SumPT2{0.}; // sum pt^2 of tracks attached to the vertex
  Double_t GenSumPT2{0.}; // sum pt^2 of gen tracks attached to the vertex

  Double_t GenDeltaZ{0.}; // distance in z to closest generated vertex
  Double_t BTVSumPT2{0.}; // sum pt^2 of tracks attached to the secondary vertex

  TRefArray Constituents; // references to constituents

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  ClassDef(Vertex, 3)
};

//---------------------------------------------------------------------------

class MissingET: public TObject
{
public:
  MissingET() = default;
  explicit MissingET(const Candidate &);

  Float_t MET{0.f}; // mising transverse energy
  Float_t Eta{0.f}; // mising energy pseudorapidity
  Float_t Phi{0.f}; // mising energy azimuthal angle

  TLorentzVector P4() const;

  ClassDef(MissingET, 1)
};

//---------------------------------------------------------------------------

class ScalarHT: public TObject
{
public:
  ScalarHT() = default;
  explicit ScalarHT(const Candidate &);

  Float_t HT{0.f}; // scalar sum of transverse momenta

  ClassDef(ScalarHT, 1)
};

//---------------------------------------------------------------------------

class Rho: public TObject
{
public:
  Float_t Rho{0.f}; // rho energy density
  Float_t Edges[2]{0.f}; // pseudorapidity range edges

  ClassDef(Rho, 1)
};

//---------------------------------------------------------------------------

class Weight: public TObject
{
public:
  Float_t Weight{0.f}; // weight for the event

  ClassDef(Weight, 1)
};

//---------------------------------------------------------------------------

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

  TRefArray Particles; // references to generated particles

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

  TRef Particle; // reference to generated particle

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

  TRef Particle; // reference to generated particle

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

  TLorentzVector SoftDroppedJet;
  TLorentzVector SoftDroppedSubJet1;
  TLorentzVector SoftDroppedSubJet2;

  TLorentzVector TrimmedP4[5]; // first entry (i = 0) is the total Trimmed Jet 4-momenta and from i = 1 to 4 are the trimmed subjets 4-momenta
  TLorentzVector PrunedP4[5]; // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta
  TLorentzVector SoftDroppedP4[5]; // first entry (i = 0) is the total SoftDropped Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta

  Int_t NSubJetsTrimmed{0}; // number of subjets trimmed
  Int_t NSubJetsPruned{0}; // number of subjets pruned
  Int_t NSubJetsSoftDropped{0}; // number of subjets soft-dropped

  Double_t ExclYmerge12;
  Double_t ExclYmerge23;
  Double_t ExclYmerge34;
  Double_t ExclYmerge45;
  Double_t ExclYmerge56;

  TRefArray Constituents; // references to constituents
  TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;
  TLorentzVector Area;

  ClassDef(Jet, 5)
};

//---------------------------------------------------------------------------

class Track: public SortableObject
{
public:
  Track() = default;
  explicit Track(const Candidate &);

  Int_t PID{0}; // HEP ID number

  Int_t Charge{0}; // track charge

  Int_t IsPU; // 0 or 1 for particles from pile-up interactions
  Int_t IsRecoPU; // 0 or 1 for reconstructed particles from pile-up
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

  TRef Particle; // reference to generated particle

  Int_t VertexIndex{-1}; // reference to vertex

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;
  TMatrixDSym CovarianceMatrix() const;

  ClassDef(Track, 4)
};

//---------------------------------------------------------------------------

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

  Int_t IsPU; // 0 or 1 for particles from pile-up interactions
  Int_t IsRecoPU; // 0 or 1 for reconstructed particles from pile-up
  Float_t HardEnergyFraction{0.f}; // fraction of hard scattering vs PU energy

  TRefArray Particles; // references to generated particles

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;

  ClassDef(Tower, 5)
};

//---------------------------------------------------------------------------

class ParticleFlowCandidate: public SortableObject
{

public:
  ParticleFlowCandidate() = default;
  explicit ParticleFlowCandidate(const Candidate &);

  Int_t PID; // HEP ID number

  Int_t Charge; // track charge

  Int_t IsPU; // 0 or 1 for particles from pile-up interactions
  Int_t IsRecoPU; // 0 or 1 for reconstructed particles from pile-up
  Float_t HardEnergyFraction{0.f}; // fraction of hard scattering vs PU energy in the particle flow candidate

  Float_t E{0.f}; // reconstructed energy [GeV]
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

  Int_t VertexIndex{-1}; // reference to vertex

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  TLorentzVector P4() const;
  TMatrixDSym CovarianceMatrix() const;

  Int_t NTimeHits{0}; // number of hits contributing to time measurement

  Float_t Eem{0.f}; // calorimeter tower electromagnetic energy
  Float_t Ehad{0.f}; // calorimeter tower hadronic energy
  Float_t Etrk{0.f}; // total charged energy hitting tower

  Float_t Edges[4]{0.f}; // calorimeter tower edges

  TRefArray Particles; // references to generated particles

  ClassDef(ParticleFlowCandidate, 4)
};

//---------------------------------------------------------------------------

class HectorHit: public SortableObject
{
public:
  HectorHit() = default;
  explicit HectorHit(const Candidate &);

  Float_t E{0.f}; // reconstructed energy [GeV]

  Float_t Tx{0.f}; // angle of the momentum in the horizontal (x,z) plane [urad]
  Float_t Ty{0.f}; // angle of the momentum in the verical (y,z) plane [urad]

  Float_t T{0.f}; // time of flight to the detector [s]

  Float_t X{0.f}; // horizontal distance to the beam [um]
  Float_t Y{0.f}; // vertical distance to the beam [um]
  Float_t S{0.f}; // distance to the interaction point [m]

  TRef Particle; // reference to generated particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  ClassDef(HectorHit, 1)
};
//---------------------------------------------------------------------------

class CscCluster: public SortableObject
{
public:
  CscCluster() = default;
  explicit CscCluster(const Candidate &);

  Float_t Eta{0.f}; // eta of LLP
  Float_t Phi{0.f}; // phi of LLP
  Float_t PT{0.f}; // pt of LLP
  Float_t Px{0.f}; // px of LLP
  Float_t Py{0.f}; // py of LLP
  Float_t Pz{0.f}; // pz of LLP
  Float_t E{0.f}; // E of LLP
  Float_t Ehad{0.f}; // had energy of LLP
  Float_t Eem{0.f}; // em energy of LLP
  Float_t pid; // LLP pid
  Float_t T{0.f}; // LLP decay time-photon travel time
  Float_t X{0.f}; // LLP decay x
  Float_t Y{0.f}; // LLP decay y
  Float_t Z{0.f}; // LLP decay z
  Float_t R{0.f}; // LLP decay z
  Float_t beta{0.f}; // LLP beta
  Float_t ctau{0.f}; // LLP ctau

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  ClassDef(CscCluster, 5)
};

//---------------------------------------------------------------------------

class Candidate: public SortableObject
{
  friend class DelphesFactory;

public:
  Candidate();

  Int_t PID{0};

  Int_t Status{0};
  Int_t M1{-1}, M2{-1}, D1{-1}, D2{-1};

  Int_t Charge{0};

  Float_t Mass{0.f};

  Int_t IsPU{0};
  Int_t IsRecoPU{0};

  Int_t IsConstituent{0};
  Int_t IsFromConversion{0};

  UInt_t Flavor{0};
  UInt_t FlavorAlgo{0};
  UInt_t FlavorPhys{0};
  UInt_t TauFlavor{0};

  UInt_t BTag{0};
  UInt_t BTagAlgo{0};
  UInt_t BTagPhys{0};

  UInt_t TauTag{0};
  Float_t TauWeight{0.f};

  Float_t Eem{0.f};
  Float_t Ehad{0.f};
  Float_t Etrk{0.f};

  Float_t Edges[4]{0.f};
  Float_t DeltaEta{0.f};
  Float_t DeltaPhi{0.f};

  TLorentzVector Momentum, Position, InitialPosition, PositionError, DecayPosition, Area;

  Float_t L{0.f}; // path length
  Float_t DZ{0.f};
  Float_t ErrorDZ{0.f};
  Float_t ErrorT{0.f}; // path length
  Float_t D0{0.f};
  Float_t ErrorD0{0.f};
  Float_t C{0.f};
  Float_t ErrorC{0.f};
  Float_t P{0.f};
  Float_t ErrorP{0.f};
  Float_t PT{0.f};
  Float_t ErrorPT{0.f};
  Float_t CtgTheta{0.f};
  Float_t ErrorCtgTheta{0.f};
  Float_t Phi{0.f};
  Float_t ErrorPhi{0.f};

  Float_t Nclusters{0.f}; // Number of ionization clusters
  Float_t dNdx{0.f}; // Number of ionization clusters per unit length

  Float_t Xd{0.f};
  Float_t Yd{0.f};
  Float_t Zd{0.f};

  Float_t XFirstHit{0.f};
  Float_t YFirstHit{0.f};
  Float_t ZFirstHit{0.f};

  // tracking resolution

  Float_t TrackResolution{0.f};

  // PileUpJetID variables

  Int_t NCharged{0};
  Int_t NNeutrals{0};
  Float_t Beta{0.f};
  Float_t BetaStar{0.f};
  Float_t MeanSqDeltaR{0.f};
  Float_t PTD{0.f};
  Float_t FracPt[5]{0.f};
  Float_t NeutralEnergyFraction{0.f}; // charged energy fraction
  Float_t ChargedEnergyFraction{0.f}; // neutral energy fraction

  // Timing information

  Int_t NTimeHits{-1};
  std::vector<std::pair<Float_t, Float_t> > ECalEnergyTimePairs;

  // Isolation variables

  Float_t IsolationVar{-999.f};
  Float_t IsolationVarRhoCorr{-999.f};
  Float_t SumPtCharged{-999.f};
  Float_t SumPtNeutral{-999.f};
  Float_t SumPtChargedPU{-999.f};
  Float_t SumPt{-999.f};

  // ACTS compliant 6x6 track covariance (D0, phi, Curvature, dz, ctg(theta))

  TMatrixDSym TrackCovariance;

  // vertex variables

  Int_t ClusterIndex{-1};
  Int_t ClusterNDF{-99};
  Double_t ClusterSigma{0.};
  Double_t SumPT2{0.};
  Double_t BTVSumPT2{0.};
  Double_t GenDeltaZ{0.};
  Double_t GenSumPT2{0.};

  // N-subjettiness variables

  Float_t Tau[5]{0.f};

  // Other Substructure variables

  TLorentzVector SoftDroppedJet;
  TLorentzVector SoftDroppedSubJet1;
  TLorentzVector SoftDroppedSubJet2;

  TLorentzVector TrimmedP4[5]; // first entry (i = 0) is the total Trimmed Jet 4-momenta and from i = 1 to 4 are the trimmed subjets 4-momenta
  TLorentzVector PrunedP4[5]; // first entry (i = 0) is the total Pruned Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta
  TLorentzVector SoftDroppedP4[5]; // first entry (i = 0) is the total SoftDropped Jet 4-momenta and from i = 1 to 4 are the pruned subjets 4-momenta

  Int_t NSubJetsTrimmed{0}; // number of subjets trimmed
  Int_t NSubJetsPruned{0}; // number of subjets pruned
  Int_t NSubJetsSoftDropped{0}; // number of subjets soft-dropped

  // Exclusive clustering variables
  Double_t ExclYmerge12{0.};
  Double_t ExclYmerge23{0.};
  Double_t ExclYmerge34{0.};
  Double_t ExclYmerge45{0.};
  Double_t ExclYmerge56{0.};

  // event characteristics variables
  Double_t ParticleDensity{0.}; // particle multiplicity density in the proximity of the particle

  static CompBase *fgCompare; //!
  const CompBase *GetCompare() const { return fgCompare; }

  void AddCandidate(const Candidate *object);
  const std::vector<Candidate *> &GetCandidates() const;

  Bool_t Overlaps(const Candidate *object) const;

  virtual void Copy(TObject &object) const;
  virtual TObject *Clone(const char *newname = "") const;
  virtual void Clear(Option_t *option = "");

private:
  DelphesFactory *fFactory{nullptr}; //!
  std::vector<Candidate *> fArray{}; //!

  void SetFactory(DelphesFactory *factory) { fFactory = factory; }

  ClassDef(Candidate, 7)
};

using CandidatesCollection = std::shared_ptr<std::vector<Candidate *> >;
//using CandidatesCollection = std::vector<Candidate *>;

#endif // DelphesClasses_h
