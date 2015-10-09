/*
Prints complete input particle arborescence for the first 100 events. Useful for debugging purposes.
root -l examples/Example5.C'("delphes_output.root")'
*/

//------------------------------------------------------------------------------

void Example5(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    if(entry>100) break;
    
    cout<<"" <<endl;
    cout<<"--------- New Event ---------" <<endl;
    cout<<"" <<endl;
 
    // loop over all input particles in the event
    for(Int_t i=0; i < branchParticle->GetEntriesFast(); i++)
    {    
     GenParticle *gen = (GenParticle*) branchParticle->At(i);     
     cout<<"N: "<<i<<", St: "<<gen->Status<<", PID: "<<gen->PID<<", E: "<<gen->E<<", Px: "<<gen->Px<<", Py: "<<gen->Py<<", Pz: "<<gen->Pz<<", M1: "<<gen->M1<<", M2: "<<gen->M2<<", D1: "<<gen->D1<<", D2: "<<gen->D2<<endl;
    }
  }
}
