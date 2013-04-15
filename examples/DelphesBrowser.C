/*
root -l examples/DelphesBrowser.C
*/

{
  #include <stdexcept>
  gSystem->Load("libDelphes");

  TFile *outputFile = TFile::Open("test.root", "RECREATE");
  ExRootTreeWriter *treeWriter = new ExRootTreeWriter(outputFile, "Delphes");

  ExRootConfReader *confReader = new ExRootConfReader;
  
  Delphes *modularDelphes = new Delphes("Delphes");

  TObjArray *allParticleOutputArray = modularDelphes->ExportArray("allParticles");
  TObjArray *stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
  TObjArray *partonOutputArray = modularDelphes->ExportArray("partons");

  confReader->ReadFile("examples/delphes_card_CMS.tcl");

  modularDelphes->SetConfReader(confReader);
  modularDelphes->SetTreeWriter(treeWriter);

  modularDelphes->InitTask();

  TBrowser browser;
}