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

/** \class RNTupleWriter
 *
 *  Fills ROOT RNTuple entries
 *
 *  \author L. Forthomme - AGH, Krakow
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesRNTupleClasses.h"
#include "classes/DelphesWriter.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <functional>

struct SortCandidates
{
  bool operator()(const Candidate *lhs, const Candidate *rhs) const { return lhs->Compare(rhs); }
};

class RNTupleWriter: public DelphesWriter
{
public:
  explicit RNTupleWriter(const DelphesParameters &moduleParams) : DelphesWriter(moduleParams) {}

  void Init() override
  {
    if(fOutputFile = std::make_unique<TFile>(GetOutputFile().data(), "CREATE");
      !fOutputFile || !fOutputFile->IsOpen() || fOutputFile->IsZombie())
    {
      std::ostringstream message;
      message << "can't create output file '" << GetOutputFile() << "'.";
      throw std::runtime_error(message.str());
    }
    std::unique_ptr<ROOT::RNTupleModel> eventModel = ROOT::RNTupleModel::Create();
    for(const std::array<std::string, 3> &branchInfo : Steer<std::vector<std::array<std::string, 3> > >("Branch"))
    {
      const std::string &branchInputArray = branchInfo.at(0), &branchName = branchInfo.at(1), &branchClassName = branchInfo.at(2);
      if(branchClassName == "CscCluster")
        eventModel->MakeField<std::vector<CscCluster> >(branchName);
      else if(branchClassName == "Electron")
        eventModel->MakeField<std::vector<DelphesRNTuple::Electron> >(branchName);
      else if(branchClassName == "GenParticle")
        eventModel->MakeField<std::vector<DelphesRNTuple::GenParticle> >(branchName);
      else if(branchClassName == "HectorHit")
        eventModel->MakeField<std::vector<HectorHit> >(branchName);
      else if(branchClassName == "Jet")
        eventModel->MakeField<std::vector<DelphesRNTuple::Jet> >(branchName);
      else if(branchClassName == "MissingET")
        eventModel->MakeField<MissingET>(branchName);
      else if(branchClassName == "Muon")
        eventModel->MakeField<std::vector<DelphesRNTuple::Muon> >(branchName);
      else if(branchClassName == "ParticleFlowCandidate")
        eventModel->MakeField<std::vector<ParticleFlowCandidate> >(branchName);
      else if(branchClassName == "Photon")
        eventModel->MakeField<std::vector<DelphesRNTuple::Photon> >(branchName);
      else if(branchClassName == "Rho")
        eventModel->MakeField<Rho>(branchName);
      else if(branchClassName == "ScalarHT")
        eventModel->MakeField<ScalarHT>(branchName);
      //else if(branchClassName == "TauJet") eventModel->MakeField<std::vector<TauJet> >(branchName);
      else if(branchClassName == "Tower")
        eventModel->MakeField<std::vector<Tower> >(branchName);
      else if(branchClassName == "Track")
        eventModel->MakeField<std::vector<Track> >(branchName);
      else
      {
        std::cout << "** WARNING: branch class name '" << branchName << "' is not supported for output." << std::endl;
        continue;
      }
      fObjTypes[branchInputArray] = std::make_pair(branchClassName, branchName);
    }
    for(const std::pair<std::string, const std::type_info *> &branchInfo : GetFactory()->GetExportCollections())
    {
      if(*branchInfo.second == typeid(std::vector<Weight>))
        eventModel->MakeField<std::vector<Weight> >(branchInfo.first);
      else if(*branchInfo.second == typeid(std::vector<LHEFWeight>))
        eventModel->MakeField<std::vector<LHEFWeight> >(branchInfo.first);
      else if(*branchInfo.second == typeid(LHEFEvent))
        eventModel->MakeField<LHEFEvent>(branchInfo.first);
      else if(*branchInfo.second == typeid(LHCOEvent))
        eventModel->MakeField<LHCOEvent>(branchInfo.first);
      else if(*branchInfo.second == typeid(HepMCEvent))
        eventModel->MakeField<HepMCEvent>(branchInfo.first);
      else
      {
        std::cout << "** WARNING: branch name '" << branchInfo.first << "' with class type '"
                  << branchInfo.second->name() << "' is not supported for output." << std::endl;
        continue;
      }
    }
    fEventWriter = ROOT::RNTupleWriter::Append(std::move(eventModel), "Delphes", *fOutputFile);

    { // write user-defined information into a dedicated RNTuple
      std::unique_ptr<ROOT::RNTupleModel> extraInfoModel = ROOT::RNTupleModel::Create();
      for(const std::pair<std::string, double> &extraInfo : Steer<std::vector<std::pair<std::string, double> > >("Info"))
        (void)extraInfoModel->MakeField<double>(extraInfo.first);
      const std::unique_ptr<ROOT::RNTupleWriter> extraInfoWriter =
        ROOT::RNTupleWriter::Append(std::move(extraInfoModel), "Info", *fOutputFile);
      std::unique_ptr<ROOT::REntry> extraModelEntry = extraInfoWriter->CreateEntry();
      for(const std::pair<std::string, double> &extraInfo : Steer<std::vector<std::pair<std::string, double> > >("Info"))
        *extraModelEntry->GetPtr<double>(extraInfo.first) = extraInfo.second;
      extraInfoWriter->Fill();
    }
  }
  void Process() override
  {
    std::unique_ptr<ROOT::REntry> eventEntry = fEventWriter->CreateEntry();
    for(const std::pair<std::string, std::pair<std::string, std::string> > collType : fObjTypes)
    {
      if(fMethodMap.count(collType.second.first) > 0)
        fMethodMap.at(collType.second.first)(*eventEntry, collType.second.second, ImportArray(collType.first));
    }
    for(const std::pair<std::string, const std::type_info *> &branchInfo : GetFactory()->GetExportCollections())
    {
      if(*branchInfo.second == typeid(std::vector<Weight>))
        *eventEntry->GetPtr<std::vector<Weight> >(branchInfo.first) = *GetFactory()->Attach<std::vector<Weight> >(branchInfo.first);
      else if(*branchInfo.second == typeid(std::vector<LHEFWeight>))
        *eventEntry->GetPtr<std::vector<LHEFWeight> >(branchInfo.first) = *GetFactory()->Attach<std::vector<LHEFWeight> >(branchInfo.first);
      else if(*branchInfo.second == typeid(LHEFEvent))
        *eventEntry->GetPtr<LHEFEvent>(branchInfo.first) = *GetFactory()->Attach<LHEFEvent>(branchInfo.first);
      else if(*branchInfo.second == typeid(LHCOEvent))
        *eventEntry->GetPtr<LHCOEvent>(branchInfo.first) = *GetFactory()->Attach<LHCOEvent>(branchInfo.first);
      else if(*branchInfo.second == typeid(HepMCEvent))
        *eventEntry->GetPtr<HepMCEvent>(branchInfo.first) = *GetFactory()->Attach<HepMCEvent>(branchInfo.first);
    }
    fEventWriter->Fill(*eventEntry);
  }

private:
  std::unique_ptr<TFile> fOutputFile;
  std::unique_ptr<ROOT::RNTupleWriter> fEventWriter;
  std::unordered_map<std::string, std::pair<std::string, std::string> > fObjTypes;

  std::unordered_map<std::string, std::function<void(ROOT::REntry &, std::string_view, const CandidatesCollection &)> > fMethodMap{
    {"CscCluster", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<CscCluster> > cscClusters = entry.GetPtr<std::vector<CscCluster> >(collName);
       for(Candidate *const &candidate : *array)
         cscClusters->emplace_back(*candidate);
     }},
    {"Electron", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<DelphesRNTuple::Electron> > electrons =
         entry.GetPtr<std::vector<DelphesRNTuple::Electron> >(collName);
       for(Candidate *const &candidate : *array)
         electrons->emplace_back(*candidate);
     }},
    {"GenParticle", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<std::vector<DelphesRNTuple::GenParticle> > genParts =
         entry.GetPtr<std::vector<DelphesRNTuple::GenParticle> >(collName);
       for(Candidate *const &candidate : *array)
         genParts->emplace_back(*candidate);
     }},
    {"HectorHit", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<std::vector<HectorHit> > hectorHits = entry.GetPtr<std::vector<HectorHit> >(collName);
       for(Candidate *const &candidate : *array)
         hectorHits->emplace_back(*candidate);
     }},
    {"Jet", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<DelphesRNTuple::Jet> > jets = entry.GetPtr<std::vector<DelphesRNTuple::Jet> >(collName);
       for(Candidate *const &candidate : *array)
         jets->emplace_back(*candidate);
     }},
    {"MissingET", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<MissingET> missingET = entry.GetPtr<MissingET>(collName);
       *missingET = MissingET(*array->at(0));
     }},
    {"Muon", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<DelphesRNTuple::Muon> > muons =
         entry.GetPtr<std::vector<DelphesRNTuple::Muon> >(collName);
       for(Candidate *const &candidate : *array)
         muons->emplace_back(*candidate);
     }},
    {"ParticleFlowCandidate", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<std::vector<ParticleFlowCandidate> > pfCands =
         entry.GetPtr<std::vector<ParticleFlowCandidate> >(collName);
       for(Candidate *const &candidate : *array)
         pfCands->emplace_back(*candidate);
     }},
    {"Photon", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<DelphesRNTuple::Photon> > photons =
         entry.GetPtr<std::vector<DelphesRNTuple::Photon> >(collName);
       for(Candidate *const &candidate : *array)
         photons->emplace_back(*candidate);
     }},
    {"Rho", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<Rho> rho = entry.GetPtr<Rho>(collName);
       const Candidate &rhoValue = *array->at(0);
       rho->Rho = rhoValue.Momentum.E();
       for(size_t i = 0; i < 2; ++i) rho->Edges[i] = rhoValue.Edges[i];
     }},
    {"ScalarHT", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<ScalarHT> scalarHT = entry.GetPtr<ScalarHT>(collName);
       *scalarHT = ScalarHT(*array->at(0));
     }},
    {"TauJet", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       // empty, for now
     }},
    {"Tower", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<std::vector<Tower> > towers = entry.GetPtr<std::vector<Tower> >(collName);
       for(Candidate *const &candidate : *array)
         towers->emplace_back(*candidate);
     }},
    {"Track", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::shared_ptr<std::vector<Track> > tracks = entry.GetPtr<std::vector<Track> >(collName);
       for(Candidate *const &candidate : *array)
         tracks->emplace_back(*candidate);
     }},
    {"Vertices", [](ROOT::REntry &entry, std::string_view collName, const CandidatesCollection &array) {
       std::sort(array->begin(), array->end(), SortCandidates{});
       std::shared_ptr<std::vector<Vertex> > vertices = entry.GetPtr<std::vector<Vertex> >(collName);
       for(Candidate *const &candidate : *array)
         vertices->emplace_back(*candidate);
     }}};
};

//------------------------------------------------------------------------------

REGISTER_MODULE("RNTupleWriter", RNTupleWriter);
