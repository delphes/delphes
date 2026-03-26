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
#include "classes/DelphesFactory.h"
#include "classes/DelphesWriter.h"

#include <ExRootAnalysis/ExRootTreeBranch.h>
#include <ExRootAnalysis/ExRootTreeWriter.h>

#include <TFile.h>
#include <TLorentzVector.h>

#include <set>

struct SortCandidates
{
  bool operator()(const Candidate *lhs, const Candidate *rhs) const { return lhs->Compare(rhs); }
};

class TreeWriter: public DelphesWriter
{
public:
  explicit TreeWriter(const DelphesParameters &moduleParams) :
    DelphesWriter(moduleParams)
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
    fClassMap[HectorHit::Class()] = &TreeWriter::ProcessHectorHit;

    // in case we are coping with vectors, we let ExRootAnalysis/ROOT know about the base object
    fExtraTypeInfos[&typeid(std::vector<Weight>)] = &typeid(Weight);
    fExtraTypeInfos[&typeid(std::vector<LHEFWeight>)] = &typeid(LHEFWeight);

    fExtraClassMap[&typeid(std::vector<Weight>)] = &TreeWriter::ProcessWeight;
    fExtraClassMap[&typeid(LHEFWeight)] = &TreeWriter::ProcessLHEFWeight;
    fExtraClassMap[&typeid(LHEFEvent)] = &TreeWriter::ProcessLHEFEvent;
    fExtraClassMap[&typeid(LHCOEvent)] = &TreeWriter::ProcessLHCOEvent;
    fExtraClassMap[&typeid(HepMCEvent)] = &TreeWriter::ProcessHepMCEvent;

    for(const std::pair<std::string, double> &info :
      Steer<std::vector<std::pair<std::string, double> > >("Info"))
      AddInfo(info.first, info.second);
  }
  ~TreeWriter()
  {
    if(fTreeWriter) fTreeWriter->Write();
  }

  void Init() override;
  void Process() override
  {
    if(!fTreeWriter) // safety check for output tree writer
      throw std::runtime_error("Output tree writer was not properly initialised.");
    for(const auto &[branch, method_array] : fBranchMap)
    {
      auto &[method, array] = method_array;
      (this->*method)(branch, array);
    }
    for(const auto &[branch, method_array] : fExtraBranchMap)
    {
      auto &[method, array] = method_array;
      (this->*method)(branch, array);
    }
    fTreeWriter->Fill();
    fTreeWriter->Clear();
  }
  void AddInfo(std::string_view name, double value) override
  {
    fExtraInfo[std::string{name}] = value;
  }

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
  void ProcessHectorHit(ExRootTreeBranch *branch, const CandidatesCollection &array);

  void ProcessWeight(ExRootTreeBranch *branch, void *const &);
  void ProcessLHEFWeight(ExRootTreeBranch *branch, void *const &);
  void ProcessLHEFEvent(ExRootTreeBranch *branch, void *const &);
  void ProcessLHCOEvent(ExRootTreeBranch *branch, void *const &);
  void ProcessHepMCEvent(ExRootTreeBranch *branch, void *const &);

#if !defined(__CINT__) && !defined(__CLING__)
  typedef void (TreeWriter::*TProcessMethod)(ExRootTreeBranch *, const CandidatesCollection &); //!
  typedef void (TreeWriter::*TVoidProcessMethod)(ExRootTreeBranch *, void *const &); //!

  typedef std::map<ExRootTreeBranch *, std::pair<TProcessMethod, CandidatesCollection> > TBranchMap; //!

  TBranchMap fBranchMap; //!

  std::unique_ptr<TFile> fOutputFile;
  std::unique_ptr<ExRootTreeWriter> fTreeWriter;

  std::map<TClass *, TProcessMethod> fClassMap; //!
#endif
  std::map<std::string, double> fExtraInfo;
  std::unordered_map<ExRootTreeBranch *, std::pair<TVoidProcessMethod, void *> > fExtraBranchMap;
  std::unordered_map<const std::type_info *, const std::type_info *> fExtraTypeInfos;
  std::unordered_map<const std::type_info *, TVoidProcessMethod> fExtraClassMap; //!
  std::unordered_map<std::string, void *> fMemorySlots;
};

//------------------------------------------------------------------------------

void TreeWriter::Init()
{
  if(fOutputFile = std::make_unique<TFile>(GetOutputFile().data(), "CREATE");
    !fOutputFile || !fOutputFile->IsOpen() || fOutputFile->IsZombie())
  {
    std::ostringstream message;
    message << "can't create output file '" << GetOutputFile() << "'.";
    throw std::runtime_error(message.str());
  }
  fTreeWriter = std::make_unique<ExRootTreeWriter>(fOutputFile.get(), "Delphes");
  for(const std::pair<std::string, double> extraInfo : fExtraInfo)
    fTreeWriter->AddInfo(extraInfo.first.data(), extraInfo.second);

  // import array with output from filter/classifier/jetfinder modules
  for(const std::array<std::string, 3> &branchInfo :
    Steer<std::vector<std::array<std::string, 3> > >("Branch"))
  {
    const std::string &branchInputArray = branchInfo.at(0), &branchName = branchInfo.at(1), &branchClassName = branchInfo.at(2);
    TClass *branchClass = TClass::GetClass(branchClassName.data());
    if(!branchClass)
      throw std::runtime_error("** ERROR: cannot find class '" + branchClassName + "'");

    if(fClassMap.count(branchClass) == 0)
      throw std::runtime_error("** ERROR: cannot create branch for class '" + branchClassName + "'");

    fBranchMap.insert(std::make_pair(
      fTreeWriter->NewBranch(branchName.data(), branchClass), // ExRootAnalysis tree branch
      std::make_pair(fClassMap.at(branchClass), ImportArray(branchInputArray))));
  }
  for(const std::pair<std::string, const std::type_info *> &branchInfo : GetFactory()->GetExportCollections())
  {
    const std::type_info *typeInfo = nullptr;
    if(fExtraClassMap.count(branchInfo.second) == 0)
    {
      if(fExtraTypeInfos.count(branchInfo.second) > 0)
        typeInfo = fExtraTypeInfos.at(branchInfo.second);
    }
    else
      typeInfo = branchInfo.second;
    if(typeInfo == nullptr)
    {
      std::cerr << "** WARNING: did not find a proper conversion rule for '" << branchInfo.first << "' object "
                << "with (mangled) type '" << branchInfo.second->name() << "'. Will not store it into the output tree.";
      continue;
    }
    TClass *branchClass = TClass::GetClass(*typeInfo);
    fMemorySlots[branchInfo.first] = GetFactory()->Attach(branchInfo.first);
    fExtraBranchMap.insert(std::make_pair(
      fTreeWriter->NewBranch(branchInfo.first.data(), branchClass), // ExRootAnalysis tree branch
      std::make_pair(fExtraClassMap.at(typeInfo), fMemorySlots.at(branchInfo.first))));
  }
  fTreeWriter->Clear();
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
  // loop over all particles
  for(Candidate *const &candidate : *array)
  {
    GenParticle *entry = static_cast<GenParticle *>(branch->NewEntry());
    *entry = GenParticle(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessVertices(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  CompBase *compare = Candidate::fgCompare;
  Candidate::fgCompare = CompSumPT2<Candidate>::Instance();
  std::sort(array->begin(), array->end(), SortCandidates{});
  Candidate::fgCompare = compare;

  // loop over all vertices
  for(Candidate *const &candidate : *array)
  {
    Vertex *entry = static_cast<Vertex *>(branch->NewEntry());
    *entry = Vertex(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessTracks(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all tracks
  for(Candidate *const &candidate : *array)
  {
    Track *entry = static_cast<Track *>(branch->NewEntry());
    *entry = Track(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessTowers(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all towers
  for(Candidate *const &candidate : *array)
  {
    Tower *entry = static_cast<Tower *>(branch->NewEntry());
    *entry = Tower(*candidate);
    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessParticleFlowCandidates(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all tracks
  for(Candidate *const &candidate : *array)
  {
    ParticleFlowCandidate *entry = static_cast<ParticleFlowCandidate *>(branch->NewEntry());
    *entry = ParticleFlowCandidate(*candidate);
    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessPhotons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all photons
  for(Candidate *const &candidate : *array)
  {
    Photon *entry = static_cast<Photon *>(branch->NewEntry());
    *entry = Photon(*candidate);
    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessElectrons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all electrons
  for(Candidate *const &candidate : *array)
  {
    Electron *entry = static_cast<Electron *>(branch->NewEntry());
    *entry = Electron(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMuons(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all muons
  for(Candidate *const &candidate : *array)
  {
    Muon *entry = static_cast<Muon *>(branch->NewEntry());
    *entry = Muon(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessJets(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all jets
  for(Candidate *const &candidate : *array)
  {
    Jet *entry = static_cast<Jet *>(branch->NewEntry());
    *entry = Jet(*candidate);
    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMissingET(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  if(!array->empty())
  {
    MissingET *entry = static_cast<MissingET *>(branch->NewEntry());
    *entry = MissingET(*array->at(0)); // get the first entry
  }
}
//------------------------------------------------------------------------------

void TreeWriter::ProcessCscCluster(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  std::sort(array->begin(), array->end(), SortCandidates{});

  // loop over all clusters
  for(Candidate *const &candidate : *array)
  {
    CscCluster *entry = static_cast<CscCluster *>(branch->NewEntry());
    *entry = CscCluster(*candidate);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessScalarHT(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  if(!array->empty())
  {
    ScalarHT *entry = static_cast<ScalarHT *>(branch->NewEntry());
    *entry = ScalarHT(*array->at(0)); // get the first entry
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessRho(ExRootTreeBranch *branch, const CandidatesCollection &array)
{ // has to stay like this w/o constructor from Candidate, as Rho.Rho exists...
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

void TreeWriter::ProcessWeight(ExRootTreeBranch *branch, void *const &array)
{
  std::vector<Weight> *weightsArray = reinterpret_cast<std::vector<Weight> *>(array);
  for(const Weight &weight : *weightsArray)
  {
    Weight *element = static_cast<Weight *>(branch->NewEntry());
    *element = weight;
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessLHEFWeight(ExRootTreeBranch *branch, void *const &array)
{
  std::vector<LHEFWeight> *weightsArray = reinterpret_cast<std::vector<LHEFWeight> *>(array);
  for(const LHEFWeight &weight : *weightsArray)
  {
    LHEFWeight *element = static_cast<LHEFWeight *>(branch->NewEntry());
    *element = weight;
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessHepMCEvent(ExRootTreeBranch *branch, void *const &array)
{
  HepMCEvent *entry = static_cast<HepMCEvent *>(branch->NewEntry());
  *entry = *reinterpret_cast<HepMCEvent *>(array);
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessLHEFEvent(ExRootTreeBranch *branch, void *const &array)
{
  LHEFEvent *entry = static_cast<LHEFEvent *>(branch->NewEntry());
  *entry = *reinterpret_cast<LHEFEvent *>(array);
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessLHCOEvent(ExRootTreeBranch *branch, void *const &array)
{
  LHCOEvent *entry = static_cast<LHCOEvent *>(branch->NewEntry());
  *entry = *reinterpret_cast<LHCOEvent *>(array);
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessHectorHit(ExRootTreeBranch *branch, const CandidatesCollection &array)
{
  // loop over all roman pot hits
  for(Candidate *const &candidate : *array)
  {
    HectorHit *entry = static_cast<HectorHit *>(branch->NewEntry());
    *entry = HectorHit(*candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("TreeWriter", TreeWriter);
