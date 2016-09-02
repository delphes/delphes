/*
Simple macro showing how to access vertex information, and tracks associated to vertices.
Computes vertex reconstruction efficiencies, merging, fake and duplicate rates.
The input file needs to be produced with cards/CMS_PhaseII/testVertexing.tcl

root -l examples/vertexAnalyzer.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void vertexAnalyzer(const char *inputFile)
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
  TClonesArray *branchGenVtx = treeReader->UseBranch("GenVertex");
  TClonesArray *branchVtx = treeReader->UseBranch("Vertex");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");

  Vertex *bestgenvtx, *bestvtx, *vtx, *recovtx;

  Track *track;
  GenParticle *particle;

  TObject *object;

  Double_t sumpt2;

  Int_t nGen = 0, nReco = 0, ndenReco = 0, nMerge = 0;
  Int_t nGenTot = 0, nRecoTot = 0, nFake = 0, nDuplicate = 0, nOverall = 0;
  Int_t nmatch;

  Double_t distance;
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    //cout<<"------------- New event --------------------  "<<endl;
    //cout<<"  "<<endl;
    //cout<<" Reconstructed vertices: "<<branchGenVtx->GetEntriesFast()<<endl;

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    sumpt2 = - 999;

    // If event contains at least 1 genvtx
    for(Int_t i = 0; i < branchGenVtx->GetEntriesFast(); ++i)
    {

      nGenTot++;
      // Take vtx
      vtx = (Vertex*) branchGenVtx->At(i);
      //cout<<"gen vertex, Sumpt2: "<<vtx->SumPT2<<", Z: "<<vtx->Z<<", T:"<<vtx->T<<endl;

      // Loop over gen vertex attached tracks
      for(Int_t j = 0; j < vtx->Constituents.GetEntriesFast(); ++j)
      {
        object = vtx->Constituents.At(j);

        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          //cout << "    GenPart pt: " << particle->PT << ", z: " << particle->Z << ", t: " << particle->T <<", PID: "<<particle->PID<<", Status: "<<particle->Status<<", PU: "<<particle->IsPU<<endl;

          if(particle->IsPU == 0) bestgenvtx = (Vertex*) branchGenVtx->At(i);
        }


      }

      nmatch = 0;
      // loop of reco vertices
      for(Int_t j = 0; j < branchVtx->GetEntriesFast(); ++j)
      {

        recovtx = (Vertex*) branchVtx->At(j);

        distance = TMath::Sqrt( (vtx->Z - recovtx->Z)*(vtx->Z - recovtx->Z)/(recovtx->ErrorZ*recovtx->ErrorZ) +
                            (vtx->T - recovtx->T)*(vtx->T - recovtx->T)/(recovtx->ErrorT*recovtx->ErrorT) );

        if(distance < 3.0)
        {
          nmatch++;
        }

      }

      if(nmatch > 1) nDuplicate++;
      if(nmatch > 0) nOverall++;

     }

     nGen++;

     Bool_t found = false;

     for(Int_t i = 0; i < branchVtx->GetEntriesFast(); ++i)
     {

       nRecoTot++;
       recovtx = (Vertex*) branchVtx->At(i);

       distance = TMath::Sqrt( (bestgenvtx->Z - recovtx->Z)*(bestgenvtx->Z - recovtx->Z)/(recovtx->ErrorZ*recovtx->ErrorZ) +
                            (bestgenvtx->T - recovtx->T)*(bestgenvtx->T - recovtx->T)/(recovtx->ErrorT*recovtx->ErrorT) );


       if(distance < 3.0 && !found)
       {
         nReco++;
         found = true;
       }


       nmatch = 0;
       // start loop over gen vertices to find out merging rate
       for(Int_t i = 0; i < branchGenVtx->GetEntriesFast(); ++i)
       {

         // Take vtx
         vtx = (Vertex*) branchGenVtx->At(i);
         // Loop over gen vertex attached tracks
         distance = TMath::Sqrt( (vtx->Z - recovtx->Z)*(vtx->Z - recovtx->Z)/(recovtx->ErrorZ*recovtx->ErrorZ) +
                            (vtx->T - recovtx->T)*(vtx->T - recovtx->T)/(recovtx->ErrorT*recovtx->ErrorT) );

         if(distance < 3.0)
         {
           nmatch++;
         }


       }

       if(nmatch > 1) nMerge++;
       if(nmatch == 0) nFake++;


     }

  } // end event loop

  cout<<"hard vertex reco eff: "<<double(nReco)/nGen<<endl;
  cout<<"all vertex reco eff: "<<double(nOverall)/nGenTot<<endl;
  cout<<"merge rate: "<<double(nMerge)/nRecoTot<<endl;
  cout<<"fake rate: "<<double(nFake)/nRecoTot<<endl;
  cout<<"duplicate rate: "<<double(nDuplicate)/nGenTot<<endl;
}

