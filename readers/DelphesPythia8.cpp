#include <stdexcept>
#include <iostream>
#include <sstream>

#include <signal.h>

#include "Pythia.h"

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//---------------------------------------------------------------------------

void ConvertInput(Pythia8::Pythia *pythia, DelphesFactory *factory, TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray)
{
  int i;

  Candidate *candidate;
  TDatabasePDG *pdg;
  TParticlePDG *pdgParticle;
  Int_t pdgCode;

  Int_t pid, status;
  Double_t px, py, pz, e;
  Double_t x, y, z, t;

  pdg = TDatabasePDG::Instance();

  for(i = 0; i < pythia->event.size(); ++i)
  {
    Pythia8::Particle &particle = pythia->event[i];

    pid = particle.id();
    status = particle.status();
    px = particle.px(); py = particle.py(); pz = particle.pz(); e = particle.e();
    x = particle.xProd(); y = particle.yProd(); z = particle.zProd(); t = particle.tProd();

    candidate = factory->NewCandidate();

    candidate->PID = pid;
    pdgCode = TMath::Abs(candidate->PID);

    candidate->Status = status;

    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;

    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;

    pdgParticle = pdg->GetParticle(pid);
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
    candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

    candidate->Momentum.SetPxPyPzE(px, py, pz, e);

    candidate->Position.SetXYZT(x, y, z, t);

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

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesPythia8";
  stringstream message;
  TFile *outputFile = 0;
  TStopwatch readStopWatch, procStopWatch;
  ExRootTreeWriter *treeWriter = 0;
  ExRootTreeBranch *branchEvent = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  Long64_t eventCounter, errorCounter;
  Long64_t numberOfEvents, timesAllowErrors;

  Pythia8::Pythia *pythia = 0;

  if(argc != 4)
  {
    cout << " Usage: " << appName << " config_file" << " pythia_card" << " output_file" << endl;
    cout << " config_file - configuration file in Tcl format," << endl;
    cout << " pythia_card - Pythia8 configuration file," << endl;
    cout << " output_file - output file in ROOT format." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    outputFile = TFile::Open(argv[3], "CREATE");

    if(outputFile == NULL)
    {
      message << "can't create output file " << argv[3];
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

    branchEvent = treeWriter->NewBranch("Event", LHEFEvent::Class());

    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);

    factory = modularDelphes->GetFactory();
    allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    // Initialize pythia
    pythia = new Pythia8::Pythia;

    if(pythia == NULL)
    {
      throw runtime_error("can't create Pythia instance");
    }

    // Read in commands from configuration file
    pythia->readFile(argv[2]);

    // Extract settings to be used in the main program
    numberOfEvents = pythia->mode("Main:numberOfEvents");
    timesAllowErrors = pythia->mode("Main:timesAllowErrors");

    pythia->init();

    ExRootProgressBar progressBar(numberOfEvents - 1);

    // Loop over all events
    errorCounter = 0;
    treeWriter->Clear();
    modularDelphes->Clear();
    readStopWatch.Start();
    for(eventCounter = 0; eventCounter < numberOfEvents && !interrupted; ++eventCounter)
    {
      if(!pythia->next())
      {
        // If failure because reached end of file then exit event loop
        if (pythia->info.atEndOfFile())
        {
          cerr << "Aborted since reached end of Les Houches Event File" << endl;
          break;
        }

        // First few failures write off as "acceptable" errors, then quit
        if (++errorCounter < timesAllowErrors) continue;
        cerr << "Event generation aborted prematurely, owing to error!" << endl;
        break;
      }

      readStopWatch.Stop();

      procStopWatch.Start();
      ConvertInput(pythia, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray);
      modularDelphes->ProcessTask();
      procStopWatch.Stop();

      treeWriter->Fill();

      treeWriter->Clear();
      modularDelphes->Clear();

      readStopWatch.Start();
      progressBar.Update(eventCounter);
    }

    progressBar.Finish();

    pythia->statistics();

    modularDelphes->FinishTask();
    treeWriter->Write();

    cout << "** Exiting..." << endl;

    delete pythia;
    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputFile;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
