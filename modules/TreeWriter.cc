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

/** \class TreeWriter
 *
 *  Fills ROOT tree branches.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>

#include <set>

using namespace std;

struct SortCandidates
{
  bool operator()(const Candidate *lhs, const Candidate *rhs) const { return lhs->Compare(rhs); }
};

class TreeWriter: public DelphesModule
{
public:
  TreeWriter() = default;

  void Init() override;
  void Process() override;

private:
  void FillParticles(Candidate *const &candidate, TRefArray *array);

  void ProcessParticles(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessVertices(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessTracks(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessTowers(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessParticleFlowCandidates(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessPhotons(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessElectrons(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessMuons(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessCscCluster(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessTauJets(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessJets(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessMissingET(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessScalarHT(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessRho(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessWeight(ExRootTreeBranch *branch, const CandidatesCollection &array);
  void ProcessHectorHit(ExRootTreeBranch *branch, const CandidatesCollection &array);

#if !defined(__CINT__) && !defined(__CLING__)
  typedef void (TreeWriter::*TProcessMethod)(ExRootTreeBranch *, const CandidatesCollection &); //!

  typedef std::map<ExRootTreeBranch *, std::pair<TProcessMethod, CandidatesCollection> > TBranchMap; //!

  TBranchMap fBranchMap; //!

  std::map<TClass *, TProcessMethod> fClassMap; //!
#endif
};

//------------------------------------------------------------------------------

void TreeWriter::Init()
{
  fClassMap[GenParticle::Class()] = &TreeWriter::ProcessParticles;
  fClassMap[Vertex::Class()] = &TreeWriter::ProcessVertices;
  fClassMap[Track::Class()] = &TreeWriter::ProcessTracks;
  fClassMap[Tower::Class()] = &TreeWriter::ProcessTowers;
  fClassMap[ParticleFlowCandidate::Class()] = &TreeWriter::ProcessParticleFlowCandidates;
  fClassMap[Photon::Class()] = &TreeWriter::ProcessPhotons;
  fClassMap[Electron::Class()] = &TreeWriter::ProcessElectrons;
  fClassMap[Muon::Class()] = &TreeWriter::ProcessMuons;
  fClassMap[CscCluster::Class()] = &TreeWriter::ProcessCscCluster;
  fClassMap[Jet::Class()] = &TreeWriter::ProcessJets;
  fClassMap[MissingET::Class()] = &TreeWriter::ProcessMissingET;
  fClassMap[ScalarHT::Class()] = &TreeWriter::ProcessScalarHT;
  fClassMap[Rho::Class()] = &TreeWriter::ProcessRho;
  fClassMap[Weight::Class()] = &TreeWriter::ProcessWeight;
  fClassMap[HectorHit::Class()] = &TreeWriter::ProcessHectorHit;

  TBranchMap::iterator itBranchMap;
  map<TClass *, TProcessMethod>::iterator itClassMap;

  // read branch configuration and
  // import array with output from filter/classifier/jetfinder modules

  ExRootConfParam param = GetParam("Branch");
  Long_t i, size;
  TString branchName, branchClassName, branchInputArray;
  TClass *branchClass;
  ExRootTreeBranch *branch;

  size = param.GetSize();
  for(i = 0; i < size / 3; ++i)
  {
    branchInputArray = param[i * 3].GetString();
    branchName = param[i * 3 + 1].GetString();
    branchClassName = param[i * 3 + 2].GetString();

    branchClass = gROOT->GetClass(branchClassName);

    if(!branchClass)
    {
      cout << "** ERROR: cannot find class '" << branchClassName << "'" << endl;
      continue;
    }

    itClassMap = fClassMap.find(branchClass);
    if(itClassMap == fClassMap.end())
    {
      cout << "** ERROR: cannot create branch for class '" << branchClassName << "'" << endl;
      continue;
    }

    CandidatesCollection array = ImportArray(branchInputArray);
    branch = NewBranch(branchName, branchClass);

    fBranchMap.insert(make_pair(branch, make_pair(itClassMap->second, array)));
  }

  param = GetParam("Info");
  TString infoName;
  Double_t infoValue;

  size = param.GetSize();
  for(i = 0; i < size / 2; ++i)
  {
    infoName = param[i * 2].GetString();
    infoValue = param[i * 2 + 1].GetDouble();

    AddInfo(infoName, infoValue);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::FillParticles(Candidate *const &candidate, TRefArray *array)
{
  std::set<Candidate *> unique_particles;
  for(Candidate *const &sub_candidate : candidate->GetCandidates())
  {
    // particle
    if(sub_candidate->GetCandidates().empty())
    {
      unique_particles.insert(sub_candidate);
      continue;
    }

    // track
    if(Candidate *const &sub_sub_candidate = static_cast<Candidate *>(sub_candidate->GetCandidates().at(0));
      sub_sub_candidate->GetCandidates().empty())
    {
      unique_particles.insert(sub_sub_candidate);
      continue;
    }

    // tower
    for(Candidate *const &sub_sub_candidate : sub_candidate->GetCandidates())
    {
      if(sub_sub_candidate->GetCandidates().empty()) continue;
      if(Candidate *sub_sub_sub_candidate = static_cast<Candidate *>(sub_sub_candidate->GetCandidates().at(0));
        sub_sub_sub_candidate && sub_sub_sub_candidate->GetCandidates().empty())
        unique_particles.insert(sub_sub_sub_candidate);
    }
  }

  array->Clear();
  for(Candidate *const &item : unique_particles)
    array->Add(item);
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessParticles(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta, rapidity;

  const Double_t c_light = 2.99792458E8;

  // loop over all particles
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    GenParticle *entry = static_cast<GenParticle *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());
    rapidity = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Rapidity());

    entry->PID = candidate->PID;

    entry->Status = candidate->Status;
    entry->IsPU = candidate->IsPU;

    entry->M1 = candidate->M1;
    entry->M2 = candidate->M2;

    entry->D1 = candidate->D1;
    entry->D2 = candidate->D2;

    entry->Charge = candidate->Charge;
    entry->Mass = candidate->Mass;

    entry->E = momentum.E();
    entry->Px = momentum.Px();
    entry->Py = momentum.Py();
    entry->Pz = momentum.Pz();
    entry->P = momentum.P();

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;

    entry->Rapidity = rapidity;

    entry->X = position.X();
    entry->Y = position.Y();
    entry->Z = position.Z();
    entry->T = position.T() * 1.0E-3 / c_light;
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessVertices(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  const Double_t c_light = 2.99792458E8;

  Double_t x, y, z, t, xError, yError, zError, tError, sigma, sumPT2, btvSumPT2, genDeltaZ, genSumPT2;
  UInt_t index, ndf;

  CompBase *compare = Candidate::fgCompare;
  Candidate::fgCompare = CompSumPT2<Candidate>::Instance();
  std::sort(array->begin(), array->end(), SortCandidates{});
  Candidate::fgCompare = compare;

  // loop over all vertices
  for(Candidate *const &candidate : *array)
  {
    index = candidate->ClusterIndex;
    ndf = candidate->ClusterNDF;
    sigma = candidate->ClusterSigma;
    sumPT2 = candidate->SumPT2;
    btvSumPT2 = candidate->BTVSumPT2;
    genDeltaZ = candidate->GenDeltaZ;
    genSumPT2 = candidate->GenSumPT2;

    x = candidate->Position.X();
    y = candidate->Position.Y();
    z = candidate->Position.Z();
    t = candidate->Position.T() * 1.0E-3 / c_light;

    xError = candidate->PositionError.X();
    yError = candidate->PositionError.Y();
    zError = candidate->PositionError.Z();
    tError = candidate->PositionError.T() * 1.0E-3 / c_light;

    Vertex *entry = static_cast<Vertex *>(branch->NewEntry());

    entry->Index = index;
    entry->NDF = ndf;
    entry->Sigma = sigma;
    entry->SumPT2 = sumPT2;
    entry->BTVSumPT2 = btvSumPT2;
    entry->GenDeltaZ = genDeltaZ;
    entry->GenSumPT2 = genSumPT2;

    entry->X = x;
    entry->Y = y;
    entry->Z = z;
    entry->T = t;

    entry->ErrorX = xError;
    entry->ErrorY = yError;
    entry->ErrorZ = zError;
    entry->ErrorT = tError;

    entry->Constituents.Clear();
    for(Candidate *const &constituent : candidate->GetCandidates())
    {
      entry->Constituents.Add(const_cast<Candidate *>(constituent));
    }
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessTracks(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signz, cosTheta, eta, p, ctgTheta, phi, m;
  const Double_t c_light = 2.99792458E8;

  // loop over all tracks
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &position = candidate->Position;

    cosTheta = TMath::Abs(position.CosTheta());
    signz = (position.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signz * 999.9 : position.Eta());

    Track *entry = static_cast<Track *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->PID = candidate->PID;

    entry->Charge = candidate->Charge;

    entry->EtaOuter = eta;
    entry->PhiOuter = position.Phi();

    entry->XOuter = position.X();
    entry->YOuter = position.Y();
    entry->ZOuter = position.Z();
    entry->TOuter = position.T() * 1.0E-3 / c_light;

    entry->L = candidate->L;

    entry->D0 = candidate->D0;
    entry->DZ = candidate->DZ;
    entry->Nclusters = candidate->Nclusters;
    entry->dNdx = candidate->dNdx;

    entry->ErrorP = candidate->ErrorP;
    entry->ErrorPT = candidate->ErrorPT;

    // diagonal covariance matrix terms
    entry->ErrorD0 = candidate->ErrorD0;
    entry->ErrorC = candidate->ErrorC;
    entry->ErrorPhi = candidate->ErrorPhi;
    entry->ErrorDZ = candidate->ErrorDZ;
    entry->ErrorCtgTheta = candidate->ErrorCtgTheta;

    // add some offdiagonal covariance matrix elements
    entry->ErrorD0Phi = candidate->TrackCovariance(0, 1) * 1.e3;
    entry->ErrorD0C = candidate->TrackCovariance(0, 2);
    entry->ErrorD0DZ = candidate->TrackCovariance(0, 3) * 1.e6;
    entry->ErrorD0CtgTheta = candidate->TrackCovariance(0, 4) * 1.e3;
    entry->ErrorPhiC = candidate->TrackCovariance(1, 2) * 1.e-3;
    entry->ErrorPhiDZ = candidate->TrackCovariance(1, 3) * 1.e3;
    entry->ErrorPhiCtgTheta = candidate->TrackCovariance(1, 4);
    entry->ErrorCDZ = candidate->TrackCovariance(2, 3);
    entry->ErrorCCtgTheta = candidate->TrackCovariance(2, 4) * 1.e-3;
    entry->ErrorDZCtgTheta = candidate->TrackCovariance(3, 4) * 1.e3;

    entry->Xd = candidate->Xd;
    entry->Yd = candidate->Yd;
    entry->Zd = candidate->Zd;

    entry->XFirstHit = candidate->XFirstHit;
    entry->YFirstHit = candidate->YFirstHit;
    entry->ZFirstHit = candidate->ZFirstHit;

    const TLorentzVector &momentum = candidate->Momentum;

    pt = momentum.Pt();
    p = momentum.P();
    phi = momentum.Phi();
    m = momentum.M();
    ctgTheta = (TMath::Tan(momentum.Theta()) != 0) ? 1 / TMath::Tan(momentum.Theta()) : 1e10;

    cosTheta = TMath::Abs(momentum.CosTheta());
    signz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signz * 999.9 : momentum.Eta());

    entry->P = p;
    entry->PT = pt;
    entry->Eta = eta;
    entry->Phi = phi;
    entry->CtgTheta = ctgTheta;
    entry->C = candidate->C;
    entry->Mass = m;

    Candidate *particle = static_cast<Candidate *>(candidate->GetCandidates().at(0));
    //const TLorentzVector &initialPosition = particle->Position;
    const TLorentzVector &initialPosition = candidate->InitialPosition;

    entry->X = initialPosition.X();
    entry->Y = initialPosition.Y();
    entry->Z = initialPosition.Z();
    entry->T = initialPosition.T() * 1.0E-3 / c_light;
    entry->ErrorT = candidate->ErrorT * 1.0E-3 / c_light;

    entry->Particle = particle;

    entry->VertexIndex = candidate->ClusterIndex;

    entry->IsPU = candidate->IsPU;
    entry->IsRecoPU = candidate->IsRecoPU;
    entry->HardEnergyFraction = candidate->IsPU ? 0.0 : 1.0;
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessTowers(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta;
  const Double_t c_light = 2.99792458E8;

  // loop over all towers
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    Tower *entry = static_cast<Tower *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->ET = pt;
    entry->E = momentum.E();
    entry->Eem = candidate->Eem;
    entry->Ehad = candidate->Ehad;
    entry->Etrk = candidate->Etrk;
    entry->Edges[0] = candidate->Edges[0];
    entry->Edges[1] = candidate->Edges[1];
    entry->Edges[2] = candidate->Edges[2];
    entry->Edges[3] = candidate->Edges[3];

    entry->T = position.T() * 1.0E-3 / c_light;
    entry->X = position.X();
    entry->Y = position.Y();
    entry->Z = position.Z();

    entry->NTimeHits = candidate->NTimeHits;

    entry->IsPU = candidate->IsPU;
    entry->IsRecoPU = candidate->IsRecoPU;

    entry->HardEnergyFraction = candidate->BetaStar;

    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessParticleFlowCandidates(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t e, pt, signz, cosTheta, eta, p, ctgTheta, phi, m;
  const Double_t c_light = 2.99792458E8;

  // loop over all tracks
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &position = candidate->Position;

    cosTheta = TMath::Abs(position.CosTheta());
    signz = (position.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signz * 999.9 : position.Eta());

    ParticleFlowCandidate *entry = static_cast<ParticleFlowCandidate *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->PID = candidate->PID;

    entry->IsPU = candidate->IsPU;
    entry->IsRecoPU = candidate->IsRecoPU;

    entry->Charge = candidate->Charge;

    if(TMath::Abs(entry->Charge) > 0.)
    {
      entry->HardEnergyFraction = entry->IsPU ? 0.0 : 1.0;
    }
    else
    {
      entry->HardEnergyFraction = candidate->BetaStar;
    }

    entry->EtaOuter = eta;
    entry->PhiOuter = position.Phi();

    entry->XOuter = position.X();
    entry->YOuter = position.Y();
    entry->ZOuter = position.Z();
    entry->TOuter = position.T() * 1.0E-3 / c_light;

    entry->L = candidate->L;

    entry->D0 = candidate->D0;
    entry->DZ = candidate->DZ;
    entry->Nclusters = candidate->Nclusters;
    entry->dNdx = candidate->dNdx;

    entry->ErrorP = candidate->ErrorP;
    entry->ErrorPT = candidate->ErrorPT;
    entry->ErrorCtgTheta = candidate->ErrorCtgTheta;

    // diagonal covariance matrix terms

    entry->ErrorD0 = candidate->ErrorD0;
    entry->ErrorC = candidate->ErrorC;
    entry->ErrorPhi = candidate->ErrorPhi;
    entry->ErrorDZ = candidate->ErrorDZ;
    entry->ErrorCtgTheta = candidate->ErrorCtgTheta;

    // add some offdiagonal covariance matrix elements
    entry->ErrorD0Phi = candidate->TrackCovariance(0, 1);
    entry->ErrorD0C = candidate->TrackCovariance(0, 2);
    entry->ErrorD0DZ = candidate->TrackCovariance(0, 3);
    entry->ErrorD0CtgTheta = candidate->TrackCovariance(0, 4);
    entry->ErrorPhiC = candidate->TrackCovariance(1, 2);
    entry->ErrorPhiDZ = candidate->TrackCovariance(1, 3);
    entry->ErrorPhiCtgTheta = candidate->TrackCovariance(1, 4);
    entry->ErrorCDZ = candidate->TrackCovariance(2, 3);
    entry->ErrorCCtgTheta = candidate->TrackCovariance(2, 4);
    entry->ErrorDZCtgTheta = candidate->TrackCovariance(3, 4);

    entry->Xd = candidate->Xd;
    entry->Yd = candidate->Yd;
    entry->Zd = candidate->Zd;

    entry->XFirstHit = candidate->XFirstHit;
    entry->YFirstHit = candidate->YFirstHit;
    entry->ZFirstHit = candidate->ZFirstHit;

    const TLorentzVector &momentum = candidate->Momentum;

    e = momentum.E();
    pt = momentum.Pt();
    p = momentum.P();
    phi = momentum.Phi();
    m = momentum.M();
    ctgTheta = (TMath::Tan(momentum.Theta()) != 0) ? 1 / TMath::Tan(momentum.Theta()) : 1e10;

    cosTheta = TMath::Abs(momentum.CosTheta());
    signz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signz * 999.9 : momentum.Eta());

    entry->E = e;
    entry->P = p;
    entry->PT = pt;
    entry->Eta = eta;
    entry->Phi = phi;
    entry->CtgTheta = ctgTheta;
    entry->C = candidate->C;
    entry->Mass = m;

    const TLorentzVector &initialPosition = candidate->InitialPosition;

    entry->X = initialPosition.X();
    entry->Y = initialPosition.Y();
    entry->Z = initialPosition.Z();
    entry->T = initialPosition.T() * 1.0E-3 / c_light;
    entry->ErrorT = candidate->ErrorT * 1.0E-3 / c_light;

    entry->VertexIndex = candidate->ClusterIndex;

    entry->Eem = candidate->Eem;
    entry->Ehad = candidate->Ehad;
    entry->Etrk = candidate->Etrk;
    entry->Edges[0] = candidate->Edges[0];
    entry->Edges[1] = candidate->Edges[1];
    entry->Edges[2] = candidate->Edges[2];
    entry->Edges[3] = candidate->Edges[3];

    //entry->T = position.T() * 1.0E-3 / c_light;
    entry->NTimeHits = candidate->NTimeHits;

    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessPhotons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta;
  const Double_t c_light = 2.99792458E8;

  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all photons
  for(Candidate *const &candidate : *array)
  {
    //TIter it1(candidate->GetCandidates());
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    Photon *entry = static_cast<Photon *>(branch->NewEntry());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;
    entry->E = momentum.E();
    entry->T = position.T() * 1.0E-3 / c_light;

    // Isolation variables

    entry->IsolationVar = candidate->IsolationVar;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->SumPtCharged = candidate->SumPtCharged;
    entry->SumPtNeutral = candidate->SumPtNeutral;
    entry->SumPtChargedPU = candidate->SumPtChargedPU;
    entry->SumPt = candidate->SumPt;

    entry->EhadOverEem = candidate->Eem > 0.0 ? candidate->Ehad / candidate->Eem : 999.9;

    // 1: prompt -- 2: non prompt -- 3: fake
    entry->Status = candidate->Status;

    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessElectrons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta;
  const Double_t c_light = 2.99792458E8;

  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all electrons
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    Electron *entry = static_cast<Electron *>(branch->NewEntry());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;

    entry->T = position.T() * 1.0E-3 / c_light;

    // displacement
    entry->D0 = candidate->D0;
    entry->ErrorD0 = candidate->ErrorD0;
    entry->DZ = candidate->DZ;
    entry->ErrorDZ = candidate->ErrorDZ;

    // Isolation variables
    entry->IsolationVar = candidate->IsolationVar;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->SumPtCharged = candidate->SumPtCharged;
    entry->SumPtNeutral = candidate->SumPtNeutral;
    entry->SumPtChargedPU = candidate->SumPtChargedPU;
    entry->SumPt = candidate->SumPt;

    entry->Charge = candidate->Charge;

    entry->EhadOverEem = 0.0;

    entry->Particle = candidate->GetCandidates().at(0);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMuons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta;

  const Double_t c_light = 2.99792458E8;

  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all muons
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    Muon *entry = static_cast<Muon *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;

    entry->T = position.T() * 1.0E-3 / c_light;

    // displacement
    entry->D0 = candidate->D0;
    entry->ErrorD0 = candidate->ErrorD0;
    entry->DZ = candidate->DZ;
    entry->ErrorDZ = candidate->ErrorDZ;

    // Isolation variables

    entry->IsolationVar = candidate->IsolationVar;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->SumPtCharged = candidate->SumPtCharged;
    entry->SumPtNeutral = candidate->SumPtNeutral;
    entry->SumPtChargedPU = candidate->SumPtChargedPU;
    entry->SumPt = candidate->SumPt;

    entry->Charge = candidate->Charge;

    entry->Particle = candidate->GetCandidates().at(0);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessJets(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t pt, signPz, cosTheta, eta;
  Double_t ecalEnergy, hcalEnergy;
  const Double_t c_light = 2.99792458E8;
  Int_t i;

  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all jets
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    Jet *entry = static_cast<Jet *>(branch->NewEntry());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;

    entry->T = position.T() * 1.0E-3 / c_light;

    entry->Mass = momentum.M();

    entry->Area = candidate->Area;

    entry->DeltaEta = candidate->DeltaEta;
    entry->DeltaPhi = candidate->DeltaPhi;

    entry->Flavor = candidate->Flavor;
    entry->FlavorAlgo = candidate->FlavorAlgo;
    entry->FlavorPhys = candidate->FlavorPhys;

    entry->BTag = candidate->BTag;

    entry->BTagAlgo = candidate->BTagAlgo;
    entry->BTagPhys = candidate->BTagPhys;

    entry->TauFlavor = candidate->TauFlavor;
    entry->TauTag = candidate->TauTag;
    entry->TauWeight = candidate->TauWeight;

    entry->Charge = candidate->Charge;

    entry->Constituents.Clear();
    ecalEnergy = 0.0;
    hcalEnergy = 0.0;
    for(Candidate *const &constituent : candidate->GetCandidates())
    {
      entry->Constituents.Add(constituent);
      ecalEnergy += constituent->Eem;
      hcalEnergy += constituent->Ehad;
    }

    entry->EhadOverEem = ecalEnergy > 0.0 ? hcalEnergy / ecalEnergy : 999.9;

    //---   Pile-Up Jet ID variables ----

    entry->NCharged = candidate->NCharged;
    entry->NNeutrals = candidate->NNeutrals;

    entry->NeutralEnergyFraction = candidate->NeutralEnergyFraction;
    entry->ChargedEnergyFraction = candidate->ChargedEnergyFraction;
    entry->Beta = candidate->Beta;
    entry->BetaStar = candidate->BetaStar;
    entry->MeanSqDeltaR = candidate->MeanSqDeltaR;
    entry->PTD = candidate->PTD;

    //--- Sub-structure variables ----

    entry->NSubJetsTrimmed = candidate->NSubJetsTrimmed;
    entry->NSubJetsPruned = candidate->NSubJetsPruned;
    entry->NSubJetsSoftDropped = candidate->NSubJetsSoftDropped;

    entry->SoftDroppedJet = candidate->SoftDroppedJet;
    entry->SoftDroppedSubJet1 = candidate->SoftDroppedSubJet1;
    entry->SoftDroppedSubJet2 = candidate->SoftDroppedSubJet2;

    for(i = 0; i < 5; i++)
    {
      entry->FracPt[i] = candidate->FracPt[i];
      entry->Tau[i] = candidate->Tau[i];
      entry->TrimmedP4[i] = candidate->TrimmedP4[i];
      entry->PrunedP4[i] = candidate->PrunedP4[i];
      entry->SoftDroppedP4[i] = candidate->SoftDroppedP4[i];
    }

    //--- exclusive clustering variables ---
    entry->ExclYmerge12 = candidate->ExclYmerge12;
    entry->ExclYmerge23 = candidate->ExclYmerge23;
    entry->ExclYmerge34 = candidate->ExclYmerge34;
    entry->ExclYmerge45 = candidate->ExclYmerge45;
    entry->ExclYmerge56 = candidate->ExclYmerge56;

    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMissingET(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // get the first entry
  if(!array->empty())
  {
    const TLorentzVector &momentum = array->at(0)->Momentum;

    MissingET *entry = static_cast<MissingET *>(branch->NewEntry());

    entry->Eta = (-momentum).Eta();
    entry->Phi = (-momentum).Phi();
    entry->MET = momentum.Pt();
  }
}
//------------------------------------------------------------------------------

void TreeWriter::ProcessCscCluster(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  Double_t signPz, cosTheta, eta;

  const Double_t c_light = 2.99792458E8; // in unit of m/s

  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all clusters
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->DecayPosition;

    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz * 999.9 : momentum.Eta());

    CscCluster *entry = static_cast<CscCluster *>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();

    entry->PT = momentum.Pt(); // pt of LLP
    entry->Px = momentum.Px(); // px of LLP
    entry->Py = momentum.Py(); // py of LLP
    entry->Pz = momentum.Pz(); // pz of LLP
    entry->E = momentum.E(); // E of LLP
    entry->pid = candidate->PID; // LLP pid
    entry->Eem = candidate->Eem; // LLP Eem
    entry->Ehad = candidate->Ehad; // LLP Ehad
    Double_t beta = momentum.P() / momentum.E();
    Double_t gamma = 1.0 / sqrt(1 - beta * beta);
    Double_t decayDistance = sqrt(pow(position.X(), 2) + pow(position.Y(), 2) + pow(position.Z(), 2)); // mm
    entry->beta = beta; // LLP pid
    entry->ctau = decayDistance / (beta * gamma); // LLP travel time in its rest frame
    entry->T = decayDistance * (1. / beta - 1) * 1.0E-3 / c_light * 1e9; // ns
    entry->X = position.X(); // LLP decay x
    entry->Y = position.Y(); //  LLP decay y
    entry->Z = position.Z(); //  LLP decay z
    entry->R = sqrt(pow(position.X(),2)+pow(position.Y(),2)); // LLP distance in transverse plane
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessScalarHT(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // get the first entry
  if(!array->empty())
  {
    ScalarHT *entry = static_cast<ScalarHT *>(branch->NewEntry());
    entry->HT = array->at(0)->Momentum.Pt();
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessRho(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all rho
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &momentum = candidate->Momentum;

    Rho *entry = static_cast<Rho *>(branch->NewEntry());

    entry->Rho = momentum.E();
    entry->Edges[0] = candidate->Edges[0];
    entry->Edges[1] = candidate->Edges[1];
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessWeight(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // get the first entry
  if(!array->empty())
  {
    Weight *entry = static_cast<Weight *>(branch->NewEntry());
    entry->Weight = array->at(0)->Momentum.E();
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessHectorHit(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all roman pot hits
  for(Candidate *const &candidate : *array)
  {
    const TLorentzVector &position = candidate->Position;
    const TLorentzVector &momentum = candidate->Momentum;

    HectorHit *entry = static_cast<HectorHit *>(branch->NewEntry());

    entry->E = momentum.E();

    entry->Tx = momentum.Px();
    entry->Ty = momentum.Py();

    entry->T = position.T();

    entry->X = position.X();
    entry->Y = position.Y();
    entry->S = position.Z();

    entry->Particle = candidate->GetCandidates().at(0);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::Process()
{
  for(const auto &[branch, method_array] : fBranchMap)
  {
    auto &[method, array] = method_array;
    (this->*method)(branch, array);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TreeWriter", TreeWriter);
