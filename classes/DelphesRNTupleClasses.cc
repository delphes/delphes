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
#include "classes/DelphesRNTupleClasses.h"
#include "classes/SortableObject.h"

static constexpr Double_t c_light = 2.99792458E8;

namespace DelphesRNTuple
{

//------------------------------------------------------------------------------

CompBase *Electron::fgCompare = CompPT<Electron>::Instance();

Electron::Electron(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  T(cand.Position.T() * 1.e-3 / c_light),
  Charge(cand.Charge),
  //Particle(cand.GetCandidates().at(0)),
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

GenParticle::GenParticle(const Candidate &cand) : ::GenParticle(cand)
{
  SetBit(kIsReferenced, false);
}

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
    //Constituents.emplace_back(constituent);
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

CompBase *Muon::fgCompare = CompPT<Muon>::Instance();

Muon::Muon(const Candidate &cand) :
  PT(cand.Momentum.Pt()),
  Eta((std::fabs(cand.Momentum.CosTheta()) == 1. ? (cand.Momentum.Pz() >= 0. ? 1. : -1.) * 999.9 : cand.Momentum.Eta())),
  Phi(cand.Momentum.Phi()),
  T(cand.Position.T() * 1.e-3 / c_light),
  Charge(cand.Charge), //Particle(cand.GetCandidates().at(0)),
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

} // namespace DelphesRNTuple
