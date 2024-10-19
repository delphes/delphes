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

CompBase *GenParticle::fgCompare = 0;
CompBase *Photon::fgCompare = CompPT<Photon>::Instance();
CompBase *Electron::fgCompare = CompPT<Electron>::Instance();
CompBase *Muon::fgCompare = CompPT<Muon>::Instance();

CompBase *Jet::fgCompare = CompPT<Jet>::Instance();
CompBase *Track::fgCompare = CompPT<Track>::Instance();
CompBase *Tower::fgCompare = CompE<Tower>::Instance();
CompBase *ParticleFlowCandidate::fgCompare = CompE<ParticleFlowCandidate>::Instance();
CompBase *HectorHit::fgCompare = CompE<HectorHit>::Instance();
CompBase *Vertex::fgCompare = CompSumPT2<Vertex>::Instance();
CompBase *Candidate::fgCompare = CompMomentumPt<Candidate>::Instance();
CompBase *CscCluster::fgCompare = CompE<CscCluster>::Instance();

//------------------------------------------------------------------------------

TLorentzVector GenParticle::P4() const
{
  TLorentzVector vec;
  vec.SetPxPyPzE(Px, Py, Pz, E);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector MissingET::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(MET, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Photon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Electron::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Muon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, 0.0);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Jet::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector Track::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

TMatrixDSym Track::CovarianceMatrix() const
{
  TMatrixDSym Cv;
  Cv.ResizeTo(5, 5);

  // convert diagonal term to original units
  Cv(0, 0)=TMath::Power(ErrorD0, 2.);
  Cv(1, 1)=TMath::Power(ErrorPhi, 2.);
  Cv(2, 2)=TMath::Power(ErrorC, 2.);
  Cv(3, 3)=TMath::Power(ErrorDZ, 2.);
  Cv(4, 4)=TMath::Power(ErrorCtgTheta, 2.);

  // off diagonal terms
  Cv(0, 1)=ErrorD0Phi;
  Cv(0, 2)=ErrorD0C;
  Cv(0, 3)=ErrorD0DZ;
  Cv(0, 4)=ErrorD0CtgTheta;
  Cv(1, 2)=ErrorPhiC;
  Cv(1, 3)=ErrorPhiDZ;
  Cv(1, 4)=ErrorPhiCtgTheta;
  Cv(2, 3)=ErrorCDZ;
  Cv(2, 4)=ErrorCCtgTheta;
  Cv(3, 4)=ErrorDZCtgTheta;

  Cv(1, 0)=Cv(0, 1);
  Cv(2, 0)=Cv(0, 2);
  Cv(3, 0)=Cv(0, 3);
  Cv(4, 0)=Cv(0, 4);
  Cv(2, 1)=Cv(1, 2);
  Cv(3, 1)=Cv(1, 3);
  Cv(4, 1)=Cv(1, 4);
  Cv(3, 2)=Cv(2, 3);
  Cv(4, 2)=Cv(2, 4);
  Cv(4, 3)=Cv(3, 4);

  return Cv;
}


//------------------------------------------------------------------------------

TLorentzVector Tower::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiE(ET, Eta, Phi, E);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector ParticleFlowCandidate::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

TMatrixDSym ParticleFlowCandidate::CovarianceMatrix() const
{
  TMatrixDSym Cv;
  Cv.ResizeTo(5, 5);

  // convert diagonal term to original units
  Cv(0, 0)=TMath::Power(ErrorD0, 2.);
  Cv(1, 1)=TMath::Power(ErrorPhi, 2.);
  Cv(2, 2)=TMath::Power(ErrorC, 2.);
  Cv(3, 3)=TMath::Power(ErrorDZ, 2.);
  Cv(4, 4)=TMath::Power(ErrorCtgTheta, 2.);

  // off diagonal terms
  Cv(0, 1)=ErrorD0Phi;
  Cv(0, 2)=ErrorD0C;
  Cv(0, 3)=ErrorD0DZ;
  Cv(0, 4)=ErrorD0CtgTheta;
  Cv(1, 2)=ErrorPhiC;
  Cv(1, 3)=ErrorPhiDZ;
  Cv(1, 4)=ErrorPhiCtgTheta;
  Cv(2, 3)=ErrorCDZ;
  Cv(2, 4)=ErrorCCtgTheta;
  Cv(3, 4)=ErrorDZCtgTheta;

  Cv(1, 0)=Cv(0, 1);
  Cv(2, 0)=Cv(0, 2);
  Cv(3, 0)=Cv(0, 3);
  Cv(4, 0)=Cv(0, 4);
  Cv(2, 1)=Cv(1, 2);
  Cv(3, 1)=Cv(1, 3);
  Cv(4, 1)=Cv(1, 4);
  Cv(3, 2)=Cv(2, 3);
  Cv(4, 2)=Cv(2, 4);
  Cv(4, 3)=Cv(3, 4);

  return Cv;
}

//------------------------------------------------------------------------------

Candidate::Candidate() :
  PID(0), Status(0), M1(-1), M2(-1), D1(-1), D2(-1),
  Charge(0), Mass(0.0),
  IsPU(0), IsRecoPU(0), IsConstituent(0), IsFromConversion(0),
  Flavor(0), FlavorAlgo(0), FlavorPhys(0), TauFlavor(0),
  BTag(0), BTagAlgo(0), BTagPhys(0),
  TauTag(0), TauWeight(0.0), Eem(0.0), Ehad(0.0), Etrk(0.0),
  DeltaEta(0.0), DeltaPhi(0.0),
  Momentum(0.0, 0.0, 0.0, 0.0),
  Position(0.0, 0.0, 0.0, 0.0),
  InitialPosition(0.0, 0.0, 0.0, 0.0),
  PositionError(0.0, 0.0, 0.0, 0.0),
  DecayPosition(0.0, 0.0, 0.0, 0.0),
  Area(0.0, 0.0, 0.0, 0.0),
  L(0),
  DZ(0), ErrorDZ(0),
  D0(0), ErrorD0(0),
  C(0), ErrorC(0),
  P(0), ErrorP(0),
  PT(0), ErrorPT(0),
  CtgTheta(0), ErrorCtgTheta(0),
  Phi(0), ErrorPhi(0),
  Nclusters(0.0),
  dNdx(0.0),
  Xd(0), Yd(0), Zd(0),
  XFirstHit(0), YFirstHit(0), ZFirstHit(0),
  TrackResolution(0),
  NCharged(0),
  NNeutrals(0),
  Beta(0),
  BetaStar(0),
  MeanSqDeltaR(0),
  PTD(0),
  NeutralEnergyFraction(0),  // neutral energy fraction
  ChargedEnergyFraction(0),  // charged energy fraction
  NTimeHits(-1),
  IsolationVar(-999),
  IsolationVarRhoCorr(-999),
  SumPtCharged(-999),
  SumPtNeutral(-999),
  SumPtChargedPU(-999),
  SumPt(-999),
  TrackCovariance(5),
  ClusterIndex(-1), ClusterNDF(0), ClusterSigma(0), SumPT2(0), BTVSumPT2(0), GenDeltaZ(0), GenSumPT2(0),
  NSubJetsTrimmed(0),
  NSubJetsPruned(0),
  NSubJetsSoftDropped(0),
  ExclYmerge12(0),
  ExclYmerge23(0),
  ExclYmerge34(0),
  ExclYmerge45(0),
  ExclYmerge56(0),
  ParticleDensity(0),
  fFactory(0),
  fArray(0)
{
  int i;
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  FracPt[0] = 0.0;
  FracPt[1] = 0.0;
  FracPt[2] = 0.0;
  FracPt[3] = 0.0;
  FracPt[4] = 0.0;
  Tau[0] = 0.0;
  Tau[1] = 0.0;
  Tau[2] = 0.0;
  Tau[3] = 0.0;
  Tau[4] = 0.0;

  SoftDroppedJet.SetXYZT(0.0, 0.0, 0.0, 0.0);
  SoftDroppedSubJet1.SetXYZT(0.0, 0.0, 0.0, 0.0);
  SoftDroppedSubJet2.SetXYZT(0.0, 0.0, 0.0, 0.0);

  for(i = 0; i < 5; ++i)
  {
    TrimmedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
    PrunedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
    SoftDroppedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
  }
}

//------------------------------------------------------------------------------

void Candidate::AddCandidate(Candidate *object)
{
  if(!fArray) fArray = fFactory->NewArray();
  fArray->Add(object);
}

//------------------------------------------------------------------------------

TObjArray *Candidate::GetCandidates()
{
  if(!fArray) fArray = fFactory->NewArray();
  return fArray;
}

//------------------------------------------------------------------------------

Bool_t Candidate::Overlaps(const Candidate *object) const
{
  const Candidate *candidate;

  if(object->GetUniqueID() == GetUniqueID()) return kTRUE;

  if(fArray)
  {
    TIter it(fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(object)) return kTRUE;
    }
  }

  if(object->fArray)
  {
    TIter it(object->fArray);
    while((candidate = static_cast<Candidate *>(it.Next())))
    {
      if(candidate->Overlaps(this)) return kTRUE;
    }
  }

  return kFALSE;
}

//------------------------------------------------------------------------------

TObject *Candidate::Clone(const char *newname) const
{
  Candidate *object = fFactory->NewCandidate();
  Copy(*object);
  return object;
}

//------------------------------------------------------------------------------

void Candidate::Copy(TObject &obj) const
{
  Candidate &object = static_cast<Candidate &>(obj);
  Candidate *candidate;

  object.PID = PID;
  object.Status = Status;
  object.M1 = M1;
  object.M2 = M2;
  object.D1 = D1;
  object.D2 = D2;
  object.Charge = Charge;
  object.Mass = Mass;
  object.IsPU = IsPU;
  object.IsRecoPU = IsRecoPU;
  object.IsConstituent = IsConstituent;
  object.IsFromConversion = IsFromConversion;
  object.ClusterIndex = ClusterIndex;
  object.ClusterNDF = ClusterNDF;
  object.ClusterSigma = ClusterSigma;
  object.SumPT2 = SumPT2;
  object.BTVSumPT2 = BTVSumPT2;
  object.GenDeltaZ = GenDeltaZ;
  object.GenSumPT2 = GenSumPT2;
  object.Flavor = Flavor;
  object.FlavorAlgo = FlavorAlgo;
  object.FlavorPhys = FlavorPhys;
  object.TauFlavor = TauFlavor;
  object.BTag = BTag;
  object.BTagAlgo = BTagAlgo;
  object.BTagPhys = BTagPhys;
  object.TauTag = TauTag;
  object.TauWeight = TauWeight;
  object.Eem = Eem;
  object.Ehad = Ehad;
  object.Etrk = Etrk;
  object.Edges[0] = Edges[0];
  object.Edges[1] = Edges[1];
  object.Edges[2] = Edges[2];
  object.Edges[3] = Edges[3];
  object.DeltaEta = DeltaEta;
  object.DeltaPhi = DeltaPhi;
  object.Momentum = Momentum;
  object.Position = Position;
  object.InitialPosition = InitialPosition;
  object.DecayPosition = DecayPosition;
  object.PositionError = PositionError;
  object.Area = Area;
  object.L = L;
  object.ErrorT = ErrorT;
  object.D0 = D0;
  object.ErrorD0 = ErrorD0;
  object.DZ = DZ;
  object.ErrorDZ = ErrorDZ;
  object.P = P;
  object.ErrorP = ErrorP;
  object.C = C;
  object.ErrorC = ErrorC;
  object.PT = PT;
  object.ErrorPT = ErrorPT;
  object.CtgTheta = CtgTheta;
  object.ErrorCtgTheta = ErrorCtgTheta;
  object.Phi = Phi;
  object.ErrorPhi = ErrorPhi;
  object.Xd = Xd;
  object.Yd = Yd;
  object.Zd = Zd;
  object.XFirstHit = XFirstHit;
  object.YFirstHit = YFirstHit;
  object.ZFirstHit = ZFirstHit;
  object.Nclusters = Nclusters;
  object.dNdx = dNdx;
  object.TrackResolution = TrackResolution;
  object.NCharged = NCharged;
  object.NNeutrals = NNeutrals;
  object.NeutralEnergyFraction = NeutralEnergyFraction;
  object.ChargedEnergyFraction = ChargedEnergyFraction;
  object.Beta = Beta;
  object.BetaStar = BetaStar;
  object.MeanSqDeltaR = MeanSqDeltaR;
  object.PTD = PTD;
  object.NTimeHits = NTimeHits;
  object.IsolationVar = IsolationVar;
  object.IsolationVarRhoCorr = IsolationVarRhoCorr;
  object.SumPtCharged = SumPtCharged;
  object.SumPtNeutral = SumPtNeutral;
  object.SumPtChargedPU = SumPtChargedPU;
  object.SumPt = SumPt;
  object.ClusterIndex = ClusterIndex;
  object.ClusterNDF = ClusterNDF;
  object.ClusterSigma = ClusterSigma;
  object.SumPT2 = SumPT2;

  object.FracPt[0] = FracPt[0];
  object.FracPt[1] = FracPt[1];
  object.FracPt[2] = FracPt[2];
  object.FracPt[3] = FracPt[3];
  object.FracPt[4] = FracPt[4];
  object.Tau[0] = Tau[0];
  object.Tau[1] = Tau[1];
  object.Tau[2] = Tau[2];
  object.Tau[3] = Tau[3];
  object.Tau[4] = Tau[4];

  object.TrimmedP4[0] = TrimmedP4[0];
  object.TrimmedP4[1] = TrimmedP4[1];
  object.TrimmedP4[2] = TrimmedP4[2];
  object.TrimmedP4[3] = TrimmedP4[3];
  object.TrimmedP4[4] = TrimmedP4[4];
  object.PrunedP4[0] = PrunedP4[0];
  object.PrunedP4[1] = PrunedP4[1];
  object.PrunedP4[2] = PrunedP4[2];
  object.PrunedP4[3] = PrunedP4[3];
  object.PrunedP4[4] = PrunedP4[4];
  object.SoftDroppedP4[0] = SoftDroppedP4[0];
  object.SoftDroppedP4[1] = SoftDroppedP4[1];
  object.SoftDroppedP4[2] = SoftDroppedP4[2];
  object.SoftDroppedP4[3] = SoftDroppedP4[3];
  object.SoftDroppedP4[4] = SoftDroppedP4[4];

  object.NSubJetsTrimmed = NSubJetsTrimmed;
  object.NSubJetsPruned = NSubJetsPruned;
  object.NSubJetsSoftDropped = NSubJetsSoftDropped;
  object.ExclYmerge12 = ExclYmerge12;
  object.ExclYmerge23 = ExclYmerge23;
  object.ExclYmerge34 = ExclYmerge34;
  object.ExclYmerge45 = ExclYmerge45;
  object.ExclYmerge56 = ExclYmerge56;

  object.SoftDroppedJet = SoftDroppedJet;
  object.SoftDroppedSubJet1 = SoftDroppedSubJet1;
  object.SoftDroppedSubJet2 = SoftDroppedSubJet2;
  object.TrackCovariance = TrackCovariance;
  object.fFactory = fFactory;
  object.fArray = 0;

  // copy cluster timing info
  copy(ECalEnergyTimePairs.begin(), ECalEnergyTimePairs.end(), back_inserter(object.ECalEnergyTimePairs));

  if(fArray && fArray->GetEntriesFast() > 0)
  {
    TIter itArray(fArray);
    TObjArray *array = object.GetCandidates();
    while((candidate = static_cast<Candidate *>(itArray.Next())))
    {
      array->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------

void Candidate::Clear(Option_t *option)
{
  int i;
  SetUniqueID(0);
  ResetBit(kIsReferenced);
  PID = 0;
  Status = 0;
  M1 = -1;
  M2 = -1;
  D1 = -1;
  D2 = -1;
  Charge = 0;
  Mass = 0.0;
  IsPU = 0;
  IsRecoPU = 0;
  IsConstituent = 0;
  IsFromConversion = 0;
  Flavor = 0;
  FlavorAlgo = 0;
  FlavorPhys = 0;
  TauFlavor = 0;
  BTag = 0;
  BTagAlgo = 0;
  BTagPhys = 0;
  TauTag = 0;
  TauWeight = 0.0;
  Eem = 0.0;
  Ehad = 0.0;
  Etrk = 0.0;
  Edges[0] = 0.0;
  Edges[1] = 0.0;
  Edges[2] = 0.0;
  Edges[3] = 0.0;
  DeltaEta = 0.0;
  DeltaPhi = 0.0;
  Momentum.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Position.SetXYZT(0.0, 0.0, 0.0, 0.0);
  InitialPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  DecayPosition.SetXYZT(0.0, 0.0, 0.0, 0.0);
  Area.SetXYZT(0.0, 0.0, 0.0, 0.0);
  TrackCovariance.Zero();
  L = 0.0;
  ErrorT = 0.0;
  D0 = 0.0;
  ErrorD0 = 0.0;
  DZ = 0.0;
  ErrorDZ = 0.0;
  P = 0.0;
  ErrorP = 0.0;
  C = 0.0;
  ErrorC = 0.0;
  PT = 0.0;
  ErrorPT = 0.0;
  CtgTheta = 0.0;
  ErrorCtgTheta = 0.0;
  Phi = 0.0;
  ErrorPhi = 0.0;
  Xd = 0.0;
  Yd = 0.0;
  Zd = 0.0;
  XFirstHit = 0.0;
  YFirstHit = 0.0;
  ZFirstHit = 0.0;
  Nclusters = 0.0;
  dNdx = 0.0;
  TrackResolution = 0.0;
  NCharged = 0;
  NNeutrals = 0;
  Beta = 0.0;
  BetaStar = 0.0;
  MeanSqDeltaR = 0.0;
  PTD = 0.0;

  NTimeHits = 0;
  ECalEnergyTimePairs.clear();

  IsolationVar = -999;
  IsolationVarRhoCorr = -999;
  SumPtCharged = -999;
  SumPtNeutral = -999;
  SumPtChargedPU = -999;
  SumPt = -999;

  ClusterIndex = -1;
  ClusterNDF = -99;
  ClusterSigma = 0.0;
  SumPT2 = 0.0;
  BTVSumPT2 = 0.0;
  GenDeltaZ = 0.0;
  GenSumPT2 = 0.0;

  FracPt[0] = 0.0;
  FracPt[1] = 0.0;
  FracPt[2] = 0.0;
  FracPt[3] = 0.0;
  FracPt[4] = 0.0;
  Tau[0] = 0.0;
  Tau[1] = 0.0;
  Tau[2] = 0.0;
  Tau[3] = 0.0;
  Tau[4] = 0.0;

  SoftDroppedJet.SetXYZT(0.0, 0.0, 0.0, 0.0);
  SoftDroppedSubJet1.SetXYZT(0.0, 0.0, 0.0, 0.0);
  SoftDroppedSubJet2.SetXYZT(0.0, 0.0, 0.0, 0.0);

  ExclYmerge12 = 0.0;
  ExclYmerge23 = 0.0;
  ExclYmerge34 = 0.0;
  ExclYmerge45 = 0.0;
  ExclYmerge56 = 0.0;
  ParticleDensity = 0.0;


  for(i = 0; i < 5; ++i)
  {
    TrimmedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
    PrunedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
    SoftDroppedP4[i].SetXYZT(0.0, 0.0, 0.0, 0.0);
  }

  NSubJetsTrimmed = 0;
  NSubJetsPruned = 0;
  NSubJetsSoftDropped = 0;

  fArray = 0;
}
