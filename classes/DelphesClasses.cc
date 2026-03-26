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

/**
 *
 *  Definition of classes to be stored in the root tree.
 *  Function CompareXYZ sorts objects by the variable XYZ that MUST be
 *  present in the data members of the root tree class of the branch.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/SortableObject.h"

static constexpr Double_t c_light = 2.99792458E8;

CompBase *Candidate::fgCompare = CompMomentumPt<Candidate>::Instance();

Candidate::Candidate() : TrackCovariance(5) {}

//------------------------------------------------------------------------------

void Candidate::AddCandidate(const Candidate *object)
{
  fArray.push_back(const_cast<Candidate *>(object));
}

//------------------------------------------------------------------------------

const std::vector<Candidate *> &Candidate::GetCandidates() const
{
  return fArray;
}

//------------------------------------------------------------------------------

Bool_t Candidate::Overlaps(const Candidate *object) const
{
  if(object->GetUniqueID() == GetUniqueID()) return true;

  for(const Candidate *candidate : fArray)
  {
    if(candidate->Overlaps(object)) return true;
  }

  for(const Candidate *candidate : object->fArray)
  {
    if(candidate->Overlaps(this)) return true;
  }

  return false;
}

//------------------------------------------------------------------------------

TObject *Candidate::Clone(const char * /*newname*/) const
{
  Candidate *object = fFactory->NewCandidate();
  Copy(*object);
  return object;
}

//------------------------------------------------------------------------------

void Candidate::Copy(TObject &obj) const
{
  Candidate &object = static_cast<Candidate &>(obj);
  object = *this;
  object.fFactory = fFactory;
  object.fArray.clear();

  // copy cluster timing info
  std::copy(ECalEnergyTimePairs.begin(), ECalEnergyTimePairs.end(), std::back_inserter(object.ECalEnergyTimePairs));

  for(const Candidate *candidate : fArray)
    object.AddCandidate(candidate);
}

//------------------------------------------------------------------------------

void Candidate::Clear(Option_t * /*option*/)
{
  *this = Candidate();
  SetUniqueID(0);
  ResetBit(kIsReferenced);
  TrackCovariance.Zero();
  ECalEnergyTimePairs.clear();
  fArray.clear();
}

//------------------------------------------------------------------------------

CompBase *CscCluster::fgCompare = CompE<CscCluster>::Instance();

CscCluster::CscCluster(const Candidate &cand) :
  Eta(std::fabs(cand.Momentum.CosTheta()) == 1.0 ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta()),
  Phi(cand.Momentum.Phi()), PT(cand.Momentum.Pt()),
  Px(cand.Momentum.Px()), Py(cand.Momentum.Py()), Pz(cand.Momentum.Pz()), E(cand.Momentum.E()),
  Ehad(cand.Ehad), Eem(cand.Eem),
  pid(cand.PID),
  X(cand.DecayPosition.X()), Y(cand.DecayPosition.Y()), Z(cand.DecayPosition.Z()),
  R(std::hypot(X, Y)), // LLP distance in transverse plane
  beta(cand.Momentum.P() / cand.Momentum.E())
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());

  const double gamma = 1.0 / std::sqrt(1 - beta * beta);
  const double decayDistance = std::hypot(cand.DecayPosition.X(), cand.DecayPosition.Y(), cand.DecayPosition.Z()); // mm
  ctau = decayDistance / (beta * gamma); // LLP travel time in its rest frame
  T = decayDistance * (1. / beta - 1) * 1.0E-3 / c_light * 1e9; // ns
}

//------------------------------------------------------------------------------

CompBase *Electron::fgCompare = CompPT<Electron>::Instance();

Electron::Electron(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  T(cand.Position.T() * 1.e-3 / c_light),
  Charge(cand.Charge),
  Particle(cand.GetCandidates().at(0)),
  // Isolation variables
  IsolationVar(cand.IsolationVar), IsolationVarRhoCorr(cand.IsolationVarRhoCorr),
  SumPtCharged(cand.SumPtCharged), SumPtNeutral(cand.SumPtNeutral), SumPtChargedPU(cand.SumPtChargedPU), SumPt(cand.SumPt),
  D0(cand.D0), DZ(cand.DZ), ErrorD0(cand.ErrorD0), ErrorDZ(cand.ErrorDZ) // displacement
{
}

TLorentzVector Electron::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *GenParticle::fgCompare = 0;

GenParticle::GenParticle(const Candidate &cand) :
  PID(cand.PID), Status(cand.Status), IsPU(cand.IsPU),
  M1(cand.M1), M2(cand.M2), D1(cand.D1), D2(cand.D2),
  Charge(cand.Charge), Mass(cand.Mass),
  E(cand.Momentum.E()), Px(cand.Momentum.Px()), Py(cand.Momentum.Py()), Pz(cand.Momentum.Pz()),
  P(cand.Momentum.P()), PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  Rapidity((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Rapidity())),
  T(cand.Position.T() * 1.e-3 / c_light), X(cand.Position.X()), Y(cand.Position.Y()), Z(cand.Position.Z())
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());
}

TLorentzVector GenParticle::P4() const { return {Px, Py, Pz, E}; }

//------------------------------------------------------------------------------

CompBase *HectorHit::fgCompare = CompE<HectorHit>::Instance();

HectorHit::HectorHit(const Candidate &cand) :
  E(cand.Momentum.E()),
  Tx(cand.Momentum.Px()), Ty(cand.Momentum.Py()),
  T(cand.Position.T()), X(cand.Position.X()), Y(cand.Position.Y()), S(cand.Position.Z()),
  Particle(cand.GetCandidates().at(0)) {}

//------------------------------------------------------------------------------

CompBase *Jet::fgCompare = CompPT<Jet>::Instance();

Jet::Jet(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  T(cand.Position.T() * 1.e-3 / c_light),
  Mass(cand.Momentum.M()),
  DeltaEta(cand.DeltaEta), DeltaPhi(cand.DeltaPhi),
  Flavor(cand.Flavor), FlavorAlgo(cand.FlavorAlgo), FlavorPhys(cand.FlavorPhys), TauFlavor(cand.TauFlavor),
  BTag(cand.BTag), BTagAlgo(cand.BTagAlgo), BTagPhys(cand.BTagPhys),
  TauTag(cand.TauTag), TauWeight(cand.TauWeight),
  Charge(cand.Charge),
  // Pile-Up Jet ID variables
  NCharged(cand.NCharged), NNeutrals(cand.NNeutrals),
  NeutralEnergyFraction(cand.NeutralEnergyFraction), ChargedEnergyFraction(cand.ChargedEnergyFraction),
  Beta(cand.Beta), BetaStar(cand.BetaStar), MeanSqDeltaR(cand.MeanSqDeltaR), PTD(cand.PTD),
  SoftDroppedJet(cand.SoftDroppedJet), SoftDroppedSubJet1(cand.SoftDroppedSubJet1), SoftDroppedSubJet2(cand.SoftDroppedSubJet2),
  // Sub-structure variables
  NSubJetsTrimmed(cand.NSubJetsTrimmed), NSubJetsPruned(cand.NSubJetsPruned), NSubJetsSoftDropped(cand.NSubJetsSoftDropped),
  // exclusive clustering variables
  ExclYmerge12(cand.ExclYmerge12), ExclYmerge23(cand.ExclYmerge23), ExclYmerge34(cand.ExclYmerge34), ExclYmerge45(cand.ExclYmerge45), ExclYmerge56(cand.ExclYmerge56),
  Area(cand.Area)
{
  double ecalEnergy = 0., hcalEnergy = 0.;
  for(Candidate *const &constituent : cand.GetCandidates())
  {
    Constituents.Add(constituent);
    ecalEnergy += constituent->Eem;
    hcalEnergy += constituent->Ehad;
  }
  EhadOverEem = ecalEnergy > 0.0 ? hcalEnergy / ecalEnergy : 999.9;

  for(size_t i = 0; i < 5; ++i)
  {
    FracPt[i] = cand.FracPt[i];
    Tau[i] = cand.Tau[i];
    TrimmedP4[i] = cand.TrimmedP4[i];
    PrunedP4[i] = cand.PrunedP4[i];
    SoftDroppedP4[i] = cand.SoftDroppedP4[i];
  }
}

TLorentzVector Jet::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

MissingET::MissingET(const Candidate &cand) : MET(cand.Momentum.Pt()), Eta((-cand.Momentum).Eta()), Phi((-cand.Momentum).Phi()) {}

TLorentzVector MissingET::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(MET, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *Muon::fgCompare = CompPT<Muon>::Instance();

Muon::Muon(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  T(cand.Position.T() * 1.e-3 / c_light),
  Charge(cand.Charge), Particle(cand.GetCandidates().at(0)),
  // Isolation variables
  IsolationVar(cand.IsolationVar), IsolationVarRhoCorr(cand.IsolationVarRhoCorr),
  SumPtCharged(cand.SumPtCharged), SumPtNeutral(cand.SumPtNeutral), SumPtChargedPU(cand.SumPtChargedPU), SumPt(cand.SumPt),
  // Displacement variables
  D0(cand.D0), DZ(cand.DZ), ErrorD0(cand.ErrorD0), ErrorDZ(cand.ErrorDZ)
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());
}

TLorentzVector Muon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *ParticleFlowCandidate::fgCompare = CompE<ParticleFlowCandidate>::Instance();

ParticleFlowCandidate::ParticleFlowCandidate(const Candidate &cand) :
  PID(cand.PID), Charge(cand.Charge),
  IsPU(cand.IsPU), IsRecoPU(cand.IsRecoPU),
  HardEnergyFraction(std::fabs(cand.Charge) > 0. ? (cand.IsPU ? 0. : 1.) : cand.BetaStar),
  E(cand.Momentum.E()), P(cand.Momentum.P()), PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  CtgTheta((std::tan(cand.Momentum.Theta()) != 0) ? 1 / std::tan(cand.Momentum.Theta()) : 1e10),
  C(cand.C),
  Mass(cand.Momentum.M()),
  EtaOuter(Eta), PhiOuter(cand.Position.Phi()),
  T(cand.InitialPosition.T() * 1.e-3 / c_light), X(cand.InitialPosition.X()), Y(cand.InitialPosition.Y()), Z(cand.InitialPosition.Z()),
  TOuter(cand.Position.T() * 1.e-3 / c_light), XOuter(cand.Position.X()), YOuter(cand.Position.Y()), ZOuter(cand.Position.Z()),
  Xd(cand.Xd), Yd(cand.Yd), Zd(cand.Zd),
  XFirstHit(cand.XFirstHit), YFirstHit(cand.YFirstHit), ZFirstHit(cand.ZFirstHit),
  L(cand.L), D0(cand.D0), DZ(cand.DZ), Nclusters(cand.Nclusters), dNdx(cand.dNdx),
  ErrorP(cand.ErrorP), ErrorPT(cand.ErrorPT), ErrorPhi(cand.ErrorPhi), ErrorCtgTheta(cand.ErrorCtgTheta),
  // diagonal covariance matrix terms
  ErrorT(cand.ErrorT * 1.e-3 / c_light), ErrorD0(cand.ErrorD0), ErrorDZ(cand.ErrorDZ), ErrorC(cand.ErrorC),
  // add some offdiagonal covariance matrix elements
  ErrorD0Phi(cand.TrackCovariance(0, 1)), ErrorD0C(cand.TrackCovariance(0, 2)), ErrorD0DZ(cand.TrackCovariance(0, 3)), ErrorD0CtgTheta(cand.TrackCovariance(0, 4)),
  ErrorPhiC(cand.TrackCovariance(1, 2)), ErrorPhiDZ(cand.TrackCovariance(1, 3)), ErrorPhiCtgTheta(cand.TrackCovariance(1, 4)),
  ErrorCDZ(cand.TrackCovariance(2, 3)), ErrorCCtgTheta(cand.TrackCovariance(2, 4)),
  ErrorDZCtgTheta(cand.TrackCovariance(3, 4)),
  VertexIndex(cand.ClusterIndex),
  NTimeHits(cand.NTimeHits),
  Eem(cand.Eem), Ehad(cand.Ehad), Etrk(cand.Etrk)
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());
  for(size_t i = 0; i < 4; ++i)
    Edges[i] = cand.Edges[i];
}

TMatrixDSym ParticleFlowCandidate::CovarianceMatrix() const
{
  TMatrixDSym Cv;
  Cv.ResizeTo(5, 5);

  // convert diagonal term to original units
  Cv(0, 0) = TMath::Power(ErrorD0, 2.);
  Cv(1, 1) = TMath::Power(ErrorPhi, 2.);
  Cv(2, 2) = TMath::Power(ErrorC, 2.);
  Cv(3, 3) = TMath::Power(ErrorDZ, 2.);
  Cv(4, 4) = TMath::Power(ErrorCtgTheta, 2.);

  // off diagonal terms
  Cv(0, 1) = ErrorD0Phi;
  Cv(0, 2) = ErrorD0C;
  Cv(0, 3) = ErrorD0DZ;
  Cv(0, 4) = ErrorD0CtgTheta;
  Cv(1, 2) = ErrorPhiC;
  Cv(1, 3) = ErrorPhiDZ;
  Cv(1, 4) = ErrorPhiCtgTheta;
  Cv(2, 3) = ErrorCDZ;
  Cv(2, 4) = ErrorCCtgTheta;
  Cv(3, 4) = ErrorDZCtgTheta;

  Cv(1, 0) = Cv(0, 1);
  Cv(2, 0) = Cv(0, 2);
  Cv(3, 0) = Cv(0, 3);
  Cv(4, 0) = Cv(0, 4);
  Cv(2, 1) = Cv(1, 2);
  Cv(3, 1) = Cv(1, 3);
  Cv(4, 1) = Cv(1, 4);
  Cv(3, 2) = Cv(2, 3);
  Cv(4, 2) = Cv(2, 4);
  Cv(4, 3) = Cv(3, 4);

  return Cv;
}

TLorentzVector ParticleFlowCandidate::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *Photon::fgCompare = CompPT<Photon>::Instance();

Photon::Photon(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  E(cand.Momentum.E()), T(cand.Position.T() * 1.e-3 / c_light),
  EhadOverEem(cand.Eem > 0. ? cand.Ehad / cand.Eem : 999.9),
  // Isolation variables
  IsolationVar(cand.IsolationVar), IsolationVarRhoCorr(cand.IsolationVarRhoCorr),
  SumPtCharged(cand.SumPtCharged), SumPtNeutral(cand.SumPtNeutral), SumPtChargedPU(cand.SumPtChargedPU), SumPt(cand.SumPt),
  Status(cand.Status)
{
}

TLorentzVector Photon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

ScalarHT::ScalarHT(const Candidate &cand) : HT(cand.Momentum.Pt()) {}

//------------------------------------------------------------------------------

CompBase *Tower::fgCompare = CompE<Tower>::Instance();

Tower::Tower(const Candidate &cand) :
  ET(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? ((cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9) : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  E(cand.Momentum.E()),
  T(cand.Position.T() * 1.e-3 / c_light),
  X(cand.Position.X()), Y(cand.Position.Y()), Z(cand.Position.Z()),
  NTimeHits(cand.NTimeHits),
  Eem(cand.Eem), Ehad(cand.Ehad), Etrk(cand.Etrk),
  IsPU(cand.IsPU), IsRecoPU(cand.IsRecoPU),
  HardEnergyFraction(cand.BetaStar)
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());
  for(size_t i = 0; i < 4; ++i)
    Edges[i] = cand.Edges[i];
}

TLorentzVector Tower::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiE(ET, Eta, Phi, E);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *Track::fgCompare = CompPT<Track>::Instance();

Track::Track(const Candidate &cand) :
  PID(cand.PID), Charge(cand.Charge),
  IsPU(cand.IsPU), IsRecoPU(cand.IsRecoPU), HardEnergyFraction(cand.IsPU ? 0. : 1.),
  P(cand.Momentum.P()), PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  CtgTheta((std::tan(cand.Momentum.Theta()) != 0) ? 1 / std::tan(cand.Momentum.Theta()) : 1e10),
  C(cand.C), Mass(cand.Momentum.M()),
  EtaOuter(Eta), PhiOuter(cand.Position.Phi()),
  T(cand.InitialPosition.T() * 1.e-3 / c_light),
  X(cand.InitialPosition.X()), Y(cand.InitialPosition.Y()), Z(cand.InitialPosition.Z()),
  TOuter(cand.Position.T() * 1.0E-3 / c_light),
  XOuter(cand.Position.X()), YOuter(cand.Position.Y()), ZOuter(cand.Position.Z()),
  Xd(cand.Xd), Yd(cand.Yd), Zd(cand.Zd),
  XFirstHit(cand.XFirstHit), YFirstHit(cand.YFirstHit), ZFirstHit(cand.ZFirstHit),
  L(cand.L), D0(cand.D0), DZ(cand.DZ),
  Nclusters(cand.Nclusters), dNdx(cand.dNdx),
  // diagonal covariance matrix terms
  ErrorP(cand.ErrorP), ErrorPT(cand.ErrorPT), ErrorPhi(cand.ErrorPhi), ErrorCtgTheta(cand.ErrorCtgTheta),
  ErrorT(cand.ErrorT * 1.0E-3 / c_light), ErrorD0(cand.ErrorD0), ErrorDZ(cand.ErrorDZ), ErrorC(cand.ErrorC),
  // add some offdiagonal covariance matrix elements
  ErrorD0Phi(cand.TrackCovariance(0, 1) * 1.e3), ErrorD0C(cand.TrackCovariance(0, 2)), ErrorD0DZ(cand.TrackCovariance(0, 3) * 1.e6), ErrorD0CtgTheta(cand.TrackCovariance(0, 4) * 1.e3),
  ErrorPhiC(cand.TrackCovariance(1, 2) * 1.e-3), ErrorPhiDZ(cand.TrackCovariance(1, 3) * 1.e3), ErrorPhiCtgTheta(cand.TrackCovariance(1, 4)),
  ErrorCDZ(cand.TrackCovariance(2, 3)), ErrorCCtgTheta(cand.TrackCovariance(2, 4) * 1.e-3),
  ErrorDZCtgTheta(cand.TrackCovariance(3, 4) * 1.e3),
  Particle(cand.GetCandidates().at(0)),
  VertexIndex(cand.ClusterIndex)
{
  SetBit(kIsReferenced);
  SetUniqueID(cand.GetUniqueID());
}

TMatrixDSym Track::CovarianceMatrix() const
{
  TMatrixDSym Cv;
  Cv.ResizeTo(5, 5);

  // convert diagonal term to original units
  Cv(0, 0) = TMath::Power(ErrorD0, 2.);
  Cv(1, 1) = TMath::Power(ErrorPhi, 2.);
  Cv(2, 2) = TMath::Power(ErrorC, 2.);
  Cv(3, 3) = TMath::Power(ErrorDZ, 2.);
  Cv(4, 4) = TMath::Power(ErrorCtgTheta, 2.);

  // off diagonal terms
  Cv(0, 1) = ErrorD0Phi;
  Cv(0, 2) = ErrorD0C;
  Cv(0, 3) = ErrorD0DZ;
  Cv(0, 4) = ErrorD0CtgTheta;
  Cv(1, 2) = ErrorPhiC;
  Cv(1, 3) = ErrorPhiDZ;
  Cv(1, 4) = ErrorPhiCtgTheta;
  Cv(2, 3) = ErrorCDZ;
  Cv(2, 4) = ErrorCCtgTheta;
  Cv(3, 4) = ErrorDZCtgTheta;

  Cv(1, 0) = Cv(0, 1);
  Cv(2, 0) = Cv(0, 2);
  Cv(3, 0) = Cv(0, 3);
  Cv(4, 0) = Cv(0, 4);
  Cv(2, 1) = Cv(1, 2);
  Cv(3, 1) = Cv(1, 3);
  Cv(4, 1) = Cv(1, 4);
  Cv(3, 2) = Cv(2, 3);
  Cv(4, 2) = Cv(2, 4);
  Cv(4, 3) = Cv(3, 4);

  return Cv;
}

TLorentzVector Track::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

CompBase *Vertex::fgCompare = CompSumPT2<Vertex>::Instance();

Vertex::Vertex(const Candidate &cand) :
  T(cand.Position.T() * 1.e-3 / c_light), X(cand.Position.X()), Y(cand.Position.Y()), Z(cand.Position.Z()),
  ErrorT(cand.PositionError.T() * 1.e-3 / c_light), ErrorX(cand.PositionError.X()), ErrorY(cand.PositionError.Y()), ErrorZ(cand.PositionError.Z()),
  Index(cand.ClusterIndex), NDF(cand.ClusterNDF),
  Sigma(cand.ClusterSigma), SumPT2(cand.SumPT2), GenSumPT2(cand.GenSumPT2),
  GenDeltaZ(cand.GenDeltaZ), BTVSumPT2(cand.BTVSumPT2)
{
  for(Candidate *const &constituent : cand.GetCandidates())
    Constituents.Add(const_cast<Candidate *>(constituent));
}

//------------------------------------------------------------------------------
