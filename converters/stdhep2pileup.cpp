#include <stdexcept>
#include <iostream>
#include <sstream>

#include <signal.h>

#include "TROOT.h"
#include "TApplication.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesSTDHEPReader.h"
#include "classes/DelphesPileUpWriter.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"

using namespace std;

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "stdhep2pileup";
  stringstream message;
  FILE *inputFile = 0;
  DelphesFactory *factory = 0;
  TObjArray *stableParticleOutputArray = 0, *allParticleOutputArray = 0, *partonOutputArray = 0;
  TIterator *itParticle = 0;
  Candidate *candidate = 0;
  DelphesPileUpWriter *writer = 0;
  DelphesSTDHEPReader *reader = 0;
  Int_t i;
  Long64_t length, eventCounter;

  if(argc < 2)
  {
    cout << " Usage: " << appName << " output_file" << " [input_file(s)]" << endl;
    cout << " output_file - output binary pile-up file," << endl;
    cout << " input_file(s) - input file(s) in STDHEP format," << endl;
    cout << " with no input_file, or when input_file is -, read standard input." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  try
  {
    writer = new DelphesPileUpWriter(argv[1]);

    factory = new DelphesFactory("ObjectFactory");
    allParticleOutputArray = factory->NewPermanentArray();
    stableParticleOutputArray = factory->NewPermanentArray();
    partonOutputArray = factory->NewPermanentArray();

    itParticle = stableParticleOutputArray->MakeIterator();

    reader = new DelphesSTDHEPReader;

    i = 2;
    do
    {
      if(interrupted) break;

      if(i == argc || strncmp(argv[i], "-", 2) == 0)
      {
        cout << "** Reading standard input" << endl;
        inputFile = stdin;
        length = -1;
      }
      else
      {
        cout << "** Reading " << argv[i] << endl;
        inputFile = fopen(argv[i], "r");

        if(inputFile == NULL)
        {
          message << "can't open " << argv[i];
          throw runtime_error(message.str());
        }

        fseek(inputFile, 0L, SEEK_END);
        length = ftello(inputFile);
        fseek(inputFile, 0L, SEEK_SET);

        if(length <= 0)
        {
          fclose(inputFile);
          ++i;
          continue;
        }
      }

      reader->SetInputFile(inputFile);

      ExRootProgressBar progressBar(length);

      // Loop over all objects
      eventCounter = 0;
      factory->Clear();
      reader->Clear();
      while(reader->ReadBlock(factory, allParticleOutputArray,
        stableParticleOutputArray, partonOutputArray) && !interrupted)
      {
        if(reader->EventReady())
        {
          ++eventCounter;

          itParticle->Reset();
          while((candidate = static_cast<Candidate*>(itParticle->Next())))
          {
            const TLorentzVector &position = candidate->Position;
            const TLorentzVector &momentum = candidate->Momentum;
            writer->WriteParticle(candidate->PID,
              position.X(), position.Y(), position.Z(), position.T(),
              momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
          }

          writer->WriteEntry();

          factory->Clear();
          reader->Clear();
        }
        progressBar.Update(ftello(inputFile), eventCounter);
      }

      fseek(inputFile, 0L, SEEK_END);
      progressBar.Update(ftello(inputFile), eventCounter, kTRUE);
      progressBar.Finish();

      if(inputFile != stdin) fclose(inputFile);

      ++i;
    }
    while(i < argc);

    writer->WriteIndex();

    cout << "** Exiting..." << endl;

    delete reader;
    delete factory;
    delete writer;

    return 0;
  }
  catch(runtime_error &e)
  {
    if(writer) delete writer;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
