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

/** \class DelphesPythia8Reader
 *
 *  Reads Pythia8 event blocks
 *
 *  \author L. Forthomme - AGH, Krakow
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"

#include <ExRootAnalysis/ExRootProgressBar.h>

#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/CombineMatchingInput.h>

class DelphesPythia8Reader: public DelphesReader
{
public:
  explicit DelphesPythia8Reader(const DelphesParameters &readerParams) :
    DelphesReader(readerParams),
    fPythia(std::make_unique<Pythia8::Pythia>())
  {
    for(const std::string &configCmd : Steer<std::vector<std::string> >("pythiaConfig"))
      fPythia->readString(configCmd);
  }

  void SetFactory(DelphesFactory *factory) override
  {
    DelphesReader::SetFactory(factory);
    fEventInfo = GetFactory()->Book<HepMCEvent>("Event", true);
    fWeightInfo = GetFactory()->Book<std::vector<Weight> >("Weights", true);
  }

  bool ReadEvent() override;
  void SetReadoutTime(double readoutTime) override { fEventInfo->ReadTime = readoutTime; }
  void SetProcessingTime(double procTime) override { fEventInfo->ProcTime = procTime; }

  void LoadInputFile(std::string_view inputFile) override { fPythia->readFile(std::string{inputFile}); }
  void Clear() override
  {
    fAllParticleOutputArray->clear();
    fStableParticleOutputArray->clear();
    fPartonOutputArray->clear();
    fWeightInfo->clear();
  }

private:
  const std::unique_ptr<Pythia8::Pythia> fPythia;
  bool fInitialised{false};

  std::shared_ptr<HepMCEvent> fEventInfo;
  std::shared_ptr<std::vector<Weight> > fWeightInfo;
};

//---------------------------------------------------------------------------

bool DelphesPythia8Reader::ReadEvent()
{
  if(!fInitialised)
  {
    // jet matching
#if PYTHIA_VERSION_INTEGER < 8300
    Pythia8::CombineMatchingInput *combined = 0;
    Pythia8::UserHooks *matching = combined->getHook(*pythia);
    if(!matching)
      throw std::runtime_error("Pythia8 cannot do matching");
    fPythia->setUserHooksPtr(matching);
#else
    Pythia8::CombineMatchingInput combined;
    combined.setHook(*fPythia);
#endif

    fPythia->init();
    fEventInfo->Number = 0;
    fEventInfo->CrossSection = fPythia->info.sigmaGen() * 1.e9; // mb -> pb
    fEventInfo->ProcessID = fPythia->info.code();
    fEventInfo->MPI = 1;
    fEventInfo->Weight = fPythia->info.weight();
    fEventInfo->Scale = fPythia->info.QRen();
    fEventInfo->ScalePDF = fPythia->info.scalup();
    fEventInfo->ID1 = fPythia->info.id1();
    fEventInfo->ID2 = fPythia->info.id2();
    fEventInfo->X1 = fPythia->info.x1();
    fEventInfo->X2 = fPythia->info.x2();
    fEventInfo->PDF1 = fPythia->info.pdf1();
    fEventInfo->PDF2 = fPythia->info.pdf2();
    fEventInfo->AlphaQED = fPythia->info.alphaEM();
    fEventInfo->AlphaQCD = fPythia->info.alphaS();

    fInitialised = true;
  }
  ResetTimer();
  if(!fPythia->next()) return false;

  DelphesFactory *factory = GetFactory();

  fEventInfo->Number += 1;
#if PYTHIA_VERSION_INTEGER > 8300
  for(const double &weight : fPythia->info.weightValueVector())
    fWeightInfo->emplace_back().Weight = weight;
#endif

  for(int i = 1 /*skip the two-beam system*/; i < fPythia->event.size(); ++i)
  {
    Pythia8::Particle &pyPart = fPythia->event[i];
    Candidate *candidate = factory->NewCandidate();
    candidate->PID = pyPart.id();
    candidate->Status = pyPart.statusHepMC();
    candidate->Charge = pyPart.charge();
    candidate->Mass = pyPart.mCalc();

    candidate->Momentum.SetPxPyPzE(pyPart.px(), pyPart.py(), pyPart.pz(), pyPart.e());
    candidate->Position.SetXYZT(pyPart.xProd(), pyPart.yProd(), pyPart.zProd(), pyPart.tProd());

    candidate->M1 = pyPart.mother1() - 1;
    candidate->M2 = pyPart.mother2() - 1;

    candidate->D1 = pyPart.daughter1() - 1;
    candidate->D2 = pyPart.daughter2() - 1;

    fAllParticleOutputArray->emplace_back(candidate);

    if(candidate->Status == 1) fStableParticleOutputArray->emplace_back(candidate);
    if(pyPart.isParton()) fPartonOutputArray->emplace_back(candidate);
  }
  SetReadoutTime(ElapsedTime());
  return true;
}

//---------------------------------------------------------------------------

REGISTER_READER("Pythia8", DelphesPythia8Reader);
