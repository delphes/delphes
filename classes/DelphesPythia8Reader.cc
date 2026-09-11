/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCLouvain), Belgium
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
 *  Reads Pythia8 data
 *
 *  \author P. Demin - UCLouvain, Louvain-la-Neuve
 *
 */

#include "classes/DelphesPythia8Reader.h"

#include <sstream>

#include <cmath>

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"

#include "Pythia.h"
#include "Pythia8Plugins/CombineMatchingInput.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"

using namespace std;

//---------------------------------------------------------------------------

DelphesPythia8Reader::DelphesPythia8Reader()
{
  fPythia = new Pythia8::Pythia;
  fCombine = new Pythia8::CombineMatchingInput;

  fCombine->setHook(*fPythia);

  fPDG = TDatabasePDG::Instance();
}

//---------------------------------------------------------------------------

DelphesPythia8Reader::~DelphesPythia8Reader()
{
  if(fCombine) delete fCombine;
  if(fPythia) delete fPythia;
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::OpenInputFile(const char *inputFileName)
{
  stringstream message;

  if(!fPythia->readFile(inputFileName))
  {
    message << "can't read Pythia8 configuration file " << inputFileName << endl;
    throw runtime_error(message.str());
  }

  fEventCounter = 0;
  fErrorCounter = 0;

  fNumberOfEvents = fPythia->mode("Main:numberOfEvents");
  fTimesAllowErrors = fPythia->mode("Main:timesAllowErrors");

  fFrameType = fPythia->mode("Beams:frameType");
  fFileNameLHEF = fPythia->word("Beams:LHEF");

  fSpareFlag1 = fPythia->flag("Main:spareFlag1");
  fSpareMode1 = fPythia->mode("Main:spareMode1");
  fSpareParm1 = fPythia->parm("Main:spareParm1");
  fSpareParm2 = fPythia->parm("Main:spareParm2");

  fPythia->init();
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::CloseInputFile()
{
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::SetInputFile(FILE *inputFile)
{
}

//---------------------------------------------------------------------------

bool DelphesPythia8Reader::ReadLHEF()
{
  return fFrameType == 4;
}

//---------------------------------------------------------------------------

string DelphesPythia8Reader::FileNameLHEF()
{
  return fFileNameLHEF;
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::Clear(Option_t * /*option*/)
{
  fEventReady = kFALSE;
}

//---------------------------------------------------------------------------

bool DelphesPythia8Reader::EventReady()
{
  return fEventReady;
}

//---------------------------------------------------------------------------

int DelphesPythia8Reader::EventNumber()
{
  return fPythia->info.nTried();
}

//---------------------------------------------------------------------------

/*
Single-particle gun. The particle must be a colour singlet.
Input: flavour, energy, direction (theta, phi).
If theta < 0 then random choice over solid angle.
Optional final argument to put particle at rest => E = m.
from pythia8 example 21
*/

static void fillParticle(int id, double pMax, double etaMax,
  Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
{
  Double_t pt, eta, phi, pp, ee, mm;

  // Reset event record to allow for new event.
  event.reset();

  // Generate uniform pt and eta.

  // pMin = 0.1 GeV for single particles
  pp = pow(10, -1.0 + (log10(pMax) + 1.0) * rndm.flat());
  eta = (2.0 * rndm.flat() - 1.0) * etaMax;
  phi = 2.0 * M_PI * rndm.flat();
  mm = pdt.mSel(id);
  ee = Pythia8::sqrtpos(pp * pp + mm * mm);
  pt = pp / cosh(eta);

  // Store the particle in the event record.
  event.append(id, 1, 0, 0, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee, mm);
}

//---------------------------------------------------------------------------

static void fillPartons(int id, double pMax, double etaMax,
  Pythia8::Event &event, Pythia8::ParticleData &pdt, Pythia8::Rndm &rndm)
{
  Double_t pt, eta, phi, pp, ee, mm;

  // Reset event record to allow for new event.
  event.reset();

  // Generate uniform pt and eta.

  // pMin = 1 GeV for jets
  pp = pow(10, log10(pMax) * rndm.flat());
  eta = (2.0 * rndm.flat() - 1.0) * etaMax;
  phi = 2.0 * M_PI * rndm.flat();
  mm = pdt.mSel(id);
  ee = Pythia8::sqrtpos(pp * pp + mm * mm);
  pt = pp / cosh(eta);

  if((id == 4 || id == 5) && pt < 10.0) return;

  if(id == 21)
  {
    event.append(21, 23, 101, 102, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee);
    event.append(21, 23, 102, 101, -pt * cos(phi), -pt * sin(phi), -pt * sinh(eta), ee);
  }
  else
  {
    event.append(id, 23, 101, 0, pt * cos(phi), pt * sin(phi), pt * sinh(eta), ee, mm);
    event.append(-id, 23, 0, 101, -pt * cos(phi), -pt * sin(phi), -pt * sinh(eta), ee, mm);
  }
}

//---------------------------------------------------------------------------

bool DelphesPythia8Reader::ReadEvent(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  ++fEventCounter;

  if(fEventCounter > fNumberOfEvents) return kFALSE;

  if(fSpareFlag1)
  {
    if((fSpareMode1 >= 1 && fSpareMode1 <= 5) || fSpareMode1 == 21)
    {
      fillPartons(fSpareMode1, fSpareParm1, fSpareParm2, fPythia->event, fPythia->particleData, fPythia->rndm);
    }
    else
    {
      fillParticle(fSpareMode1, fSpareParm1, fSpareParm2, fPythia->event, fPythia->particleData, fPythia->rndm);
    }
  }

  fEventReady = fPythia->next();

  if(!fEventReady)
  {
    ++fErrorCounter;

    if(fPythia->info.atEndOfFile()) return kFALSE;

    if(fErrorCounter > fTimesAllowErrors) return kFALSE;

    return kTRUE;
  }

  AnalyzeParticles(factory, allParticleOutputArray,
    stableParticleOutputArray, partonOutputArray);

  return kTRUE;
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::AnalyzeEvent(ExRootTreeBranch *branch, long long eventNumber,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  HepMCEvent *element;

  element = static_cast<HepMCEvent *>(branch->NewEntry());

  element->Number = eventNumber;

  element->ProcessID = fPythia->info.code();
  element->MPI = 1;
  element->Weight = fPythia->info.weight();

  element->Scale = fPythia->info.QRen();
  element->AlphaQED = fPythia->info.alphaEM();
  element->AlphaQCD = fPythia->info.alphaS();

  element->ID1 = fPythia->info.id1();
  element->ID2 = fPythia->info.id2();
  element->X1 = fPythia->info.x1();
  element->X2 = fPythia->info.x2();
  element->ScalePDF = fPythia->info.QFac();
  element->PDF1 = fPythia->info.pdf1();
  element->PDF2 = fPythia->info.pdf2();

  element->ReadTime = readStopWatch ? readStopWatch->RealTime() : 0;
  element->ProcTime = procStopWatch ? procStopWatch->RealTime() : 0;
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::AnalyzeWeight(ExRootTreeBranch *branch)
{
  Int_t i;
  Weight *element;

  for(i = 0; i < fPythia->info.numberOfWeights(); ++i)
  {
    element = static_cast<Weight *>(branch->NewEntry());

    element->Weight = fPythia->info.weightValueByIndex(i);
  }
}

//---------------------------------------------------------------------------

void DelphesPythia8Reader::AnalyzeParticles(DelphesFactory *factory,
  TObjArray *allParticleOutputArray,
  TObjArray *stableParticleOutputArray,
  TObjArray *partonOutputArray)
{
  Int_t i;
  Candidate *candidate;
  TParticlePDG *pdgParticle;

  Int_t pdgCode, pid, status;
  Double_t px, py, pz, e, mass;
  Double_t x, y, z, t;
  Double_t x_decay, y_decay, z_decay, t_decay;

  for(i = 1; i < fPythia->event.size(); ++i)
  {
    Pythia8::Particle &particle = fPythia->event[i];

    pid = particle.id();
    status = particle.statusHepMC();

    px = particle.px();
    py = particle.py();
    pz = particle.pz();
    e = particle.e();
    mass = particle.m();

    x = particle.xProd();
    y = particle.yProd();
    z = particle.zProd();
    t = particle.tProd();

    x_decay = particle.xDec();
    y_decay = particle.yDec();
    z_decay = particle.zDec();
    t_decay = particle.tDec();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;

    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;

    pdgParticle = fPDG->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = mass;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);
    candidate->DecayPosition.SetXYZT(x_decay, y_decay, z_decay, t_decay);

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1)
    {
      stableParticleOutputArray->Add(candidate);
    }
    else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
    {
      partonOutputArray->Add(candidate);
    }
  }
}

//---------------------------------------------------------------------------
