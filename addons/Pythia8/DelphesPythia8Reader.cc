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
#include "classes/DelphesLHEFReader.h"

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
    if(const std::string configFile = Steer<std::string>("pythiaConfigFile"); !configFile.empty())
      fPythia->readFile(configFile);
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
  void SetReadoutTime(double readoutTime) override
  {
    fEventInfo->ReadTime = readoutTime;
    if(fLHEReader) fLHEReader->SetReadoutTime(readoutTime);
  }
  void SetProcessingTime(double procTime) override
  {
    fEventInfo->ProcTime = procTime;
    if(fLHEReader) fLHEReader->SetProcessingTime(procTime);
  }

  void LoadInputFile(std::string_view inputFile) override { fPythia->readFile(std::string{inputFile}); }
  void Clear() override
  {
    fAllParticleOutputArray->clear();
    fStableParticleOutputArray->clear();
    fPartonOutputArray->clear();
    if(fAllParticleOutputArray) fAllParticleOutputArray->clear();
    if(fStableParticleOutputArray) fStableParticleOutputArray->clear();
    if(fPartonOutputArray) fPartonOutputArray->clear();
    fWeightInfo->clear();
  }

private:
  /** Single-particle gun. The particle must be a colour singlet.
   * \param[in] id particle flavour
   * \param[in] pMax maximum particle momentum
   * \itarm[in] etaMax maximum particle pseudo-rapidity
   * \note If theta < 0 then random choice over solid angle. From pythia8 example 21
   */
  void fillParticle(int id, double pMax, double etaMax, Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
  {
    event.reset(); // Reset event record to allow for new event.

    // Generate uniform pt and eta.
    // pMin = 0.1 GeV for single particles
    const double pp = std::pow(10, -1. + (std::log10(pMax) + 1.) * rndm.flat());
    const double eta = (2. * rndm.flat() - 1.) * etaMax;
    const double phi = 2. * M_PI * rndm.flat();
    const double mm = pdt.mSel(id);
    const double ee = Pythia8::sqrtpos(pp * pp + mm * mm);
    const double pt = pp / std::cosh(eta);

    // Store the particle in the event record.
    event.append(id, 1, 0, 0, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee, mm);
  }
  void fillPartons(int id, double pMax, double etaMax, Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
  {
    event.reset(); // Reset event record to allow for new event.

    // Generate uniform pt and eta.
    // pMin = 1 GeV for jets
    const double pp = std::pow(10, std::log10(pMax) * rndm.flat());
    const double eta = (2. * rndm.flat() - 1.) * etaMax;
    const double phi = 2. * M_PI * rndm.flat();
    const double mm = pdt.mSel(id);
    const double ee = Pythia8::sqrtpos(pp * pp + mm * mm);
    const double pt = pp / std::cosh(eta);

    if((id == 4 || id == 5) && pt < 10.) return;
    if(id == 21) // particular case for gluons
    {
      event.append(21, 23, 101, 102, +pt * std::cos(phi), +pt * std::sin(phi), +pt * std::sinh(eta), ee);
      event.append(21, 23, 102, 101, -pt * std::cos(phi), -pt * std::sin(phi), -pt * std::sinh(eta), ee);
      return;
    }
    event.append(+id, 23, 101, 0, +pt * std::cos(phi), +pt * std::sin(phi), +pt * std::sinh(eta), ee, mm);
    event.append(-id, 23, 0, 101, -pt * std::cos(phi), -pt * std::sin(phi), -pt * std::sinh(eta), ee, mm);
  }

  const std::unique_ptr<Pythia8::Pythia> fPythia;
  bool fInitialised{false};

  std::shared_ptr<HepMCEvent> fEventInfo;
  std::shared_ptr<std::vector<Weight> > fWeightInfo;

  // external LHEF input mode objects
  std::unique_ptr<DelphesLHEFReader> fLHEReader;
  std::shared_ptr<HepMCEvent> fEventInfoLHEF;
  std::shared_ptr<std::vector<Weight> > fWeightsInfoLHEF;
  CandidatesCollection fAllParticleOutputArrayLHEF;
  CandidatesCollection fStableParticleOutputArrayLHEF;
  CandidatesCollection fPartonOutputArrayLHEF;

  // particle gun mode flags
  bool fSpareFlag1{false};
  int fSpareMode1{-1};
  double fSpareParm1{0.};
  double fSpareParm2{0.};
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

    fSpareFlag1 = fPythia->flag("Main:spareFlag1");
    fSpareMode1 = fPythia->mode("Main:spareMode1");
    fSpareParm1 = fPythia->parm("Main:spareParm1");
    fSpareParm2 = fPythia->parm("Main:spareParm2");
    if(!fSpareFlag1)
      if(const auto inputFile = fPythia->word("Beams:LHEF"); !inputFile.empty())
      {
        fLHEReader = std::make_unique<DelphesLHEFReader>(DelphesParameters{});
        fLHEReader->LoadInputFile(inputFile);

        fEventInfoLHEF = GetFactory()->Book<HepMCEvent>("EventLHEF", true);
        fWeightsInfoLHEF = GetFactory()->Book<std::vector<Weight> >("WeightLHEF", true);

        fAllParticleOutputArrayLHEF = ImportArray("Delphes/allParticlesLHEF");
        fStableParticleOutputArrayLHEF = ImportArray("Delphes/stableParticlesLHEF");
        fPartonOutputArrayLHEF = ImportArray("Delphes/partonsLHEF");
      }

    fInitialised = true;
  }
  ResetTimer();

  if(fSpareFlag1)
  {
    if((fSpareMode1 >= 1 && fSpareMode1 <= 5) || fSpareMode1 == 21)
      fillPartons(fSpareMode1, fSpareParm1, fSpareParm2, fPythia->event, fPythia->particleData, fPythia->rndm);
    else
      fillParticle(fSpareMode1, fSpareParm1, fSpareParm2, fPythia->event, fPythia->particleData, fPythia->rndm);
  }
  else if(fLHEReader)
  {
    fLHEReader->Clear();
    if(!fLHEReader->ReadEvent()) return false; // failed to read the next event in parallel to Pythia
  }

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
