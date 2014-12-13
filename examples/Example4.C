/*

This macro shows how to compute jet energy scale.
root -l examples/Example4.C'("delphes_output.root", "plots.root")'

The output ROOT file contains the pT(MC)/pT(Reco) distributions for various pT(Reco) and |eta| bins.
The peak value of such distribution is interpreted as the jet energy correction to be applied for that given pT(Reco), |eta| bin.

This can be done by modifying the "ScaleFormula" input parameter to the JetEnergyScale module in the delphes_card_XXX.tcl



e.g  a smooth function:


  set ScaleFormula { sqrt(3.0 - 0.1*(abs(eta)))^2 / pt + 1.0 ) }


or a binned function:


  set ScaleFormula {(abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 50.0)  * (1.10) +
                    (abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 100.0) * (1.05) +
                    (abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 100.0)               * (1.00) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 50.0)  * (1.10) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 100.0) * (1.05) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 100.0)               * (1.00)}


Be aware that a binned jet energy scale can produce "steps" in the corrected jet pt distribution ...



*/

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fJetPT;

  TH1 *fJetRes_Pt_20_50_Eta_0_25;
  TH1 *fJetRes_Pt_20_50_Eta_25_5;

  TH1 *fJetRes_Pt_50_100_Eta_0_25;
  TH1 *fJetRes_Pt_50_100_Eta_25_5;

  TH1 *fJetRes_Pt_100_200_Eta_0_25;
  TH1 *fJetRes_Pt_100_200_Eta_25_5;

  TH1 *fJetRes_Pt_200_500_Eta_0_25;
  TH1 *fJetRes_Pt_200_500_Eta_25_5;

  TH1 *fJetRes_Pt_500_inf_Eta_0_25;
  TH1 *fJetRes_Pt_500_inf_Eta_25_5;

};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fJetPT = result->AddHist1D(
    "jet_pt", "p_{T}^{jet}",
    "p_{T}^{jet}  GeV/c", "number of jets",
    100, 0.0, 1000.0);

  plots->fJetRes_Pt_20_50_Eta_0_25 = result->AddHist1D(
    "jet_delta_pt_20_50_cen", "p_{T}^{truth,parton}/p_{T}^{jet} , 20 < p_{T} < 50 , 0 < | #eta | < 2.5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_20_50_Eta_0_25->SetStats();

  plots->fJetRes_Pt_20_50_Eta_25_5 = result->AddHist1D(
    "jet_delta_pt_20_50_fwd", "p_{T}^{truth,parton}/p_{T}^{jet} , 20 < p_{T} < 50 , 2.5 < | #eta | < 5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_20_50_Eta_25_5->SetStats();

  plots->fJetRes_Pt_50_100_Eta_0_25 = result->AddHist1D(
    "jet_delta_pt_50_100_cen", "p_{T}^{truth,parton}/p_{T}^{jet} , 50 < p_{T} < 100 , 0 < | #eta | < 2.5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_50_100_Eta_0_25->SetStats();

  plots->fJetRes_Pt_50_100_Eta_25_5 = result->AddHist1D(
    "jet_delta_pt_50_100_fwd", "p_{T}^{truth,parton}/p_{T}^{jet} , 50 < p_{T} < 100 , 2.5 < | #eta | < 5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_50_100_Eta_25_5->SetStats();


  plots->fJetRes_Pt_100_200_Eta_0_25 = result->AddHist1D(
    "jet_delta_pt_100_200_cen", "p_{T}^{truth,parton}/p_{T}^{jet} , 100 < p_{T} < 200 , 0 < | #eta | < 2.5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_100_200_Eta_0_25->SetStats();

  plots->fJetRes_Pt_100_200_Eta_25_5 = result->AddHist1D(
    "jet_delta_pt_100_200_fwd", "p_{T}^{truth,parton}/p_{T}^{jet} , 100 < p_{T} < 200 , 2.5 < | #eta | < 5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_100_200_Eta_25_5->SetStats();

  plots->fJetRes_Pt_200_500_Eta_0_25 = result->AddHist1D(
    "jet_delta_pt_200_500_cen", "p_{T}^{truth,parton}/p_{T}^{jet} , 200 < p_{T} < 500 , 0 < | #eta | < 2.5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_200_500_Eta_0_25->SetStats();

  plots->fJetRes_Pt_200_500_Eta_25_5 = result->AddHist1D(
    "jet_delta_pt_200_500_fwd", "p_{T}^{truth,parton}/p_{T}^{jet} , 200 < p_{T} < 500 , 2.5 < | #eta | < 5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_200_500_Eta_25_5->SetStats();

  plots->fJetRes_Pt_500_inf_Eta_0_25 = result->AddHist1D(
    "jet_delta_pt_500_1000_cen", "p_{T}^{truth,parton}/p_{T}^{jet} , 500 < p_{T} < 1000, 0 < | #eta | < 2.5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_500_inf_Eta_0_25->SetStats();

  plots->fJetRes_Pt_500_inf_Eta_25_5 = result->AddHist1D(
    "jet_delta_pt_500_1000_fwd", "p_{T}^{truth,parton}/p_{T}^{jet} , 500 < p_{T} < 1000, 2.5 < | #eta | < 5 ",
    "p_{T}^{truth,parton}/p_{T}^{jet}", "number of jets",
    100, 0.0, 2.0);

  plots->fJetRes_Pt_500_inf_Eta_25_5->SetStats();


}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Jet *jet, *genjet;
  GenParticle *part;
  TObject *object;

  TLorentzVector JetMom, GenJetMom, BestGenJetMom;

  Float_t Dr;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    //  cout<<"--  New event -- "<<endl;

    if(entry%500 == 0) cout << "Event number: "<< entry <<endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {

      jet = (Jet*) branchJet->At(i);
      JetMom = jet-> P4();

      plots->fJetPT->Fill(JetMom.Pt());

      Dr = 999;

     // Loop over all hard partons in event
     for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
     {

        part = (GenParticle*) branchParticle->At(j);

        GenJetMom = part -> P4();

	//this is simply to avoid warnings from initial state particle having infite rapidity ...
	if(GenJetMom.Px() == 0 && GenJetMom.Py() == 0) continue;

        //take the closest parton candidate
        if( GenJetMom.DeltaR(JetMom) < Dr )
        {
           Dr = GenJetMom.DeltaR(JetMom);
           BestGenJetMom = GenJetMom;
        }

      }

     if(Dr < 0.3)
     {
       pt  = JetMom.Pt();
       eta = TMath::Abs(JetMom.Eta());


       if( pt > 20.0 && pt < 50.0 && eta > 0.0 && eta < 2.5 ) plots -> fJetRes_Pt_20_50_Eta_0_25->Fill(BestGenJetMom.Pt()/JetMom.Pt());
       if( pt > 20.0 && pt < 50.0 && eta > 2.5 && eta < 5.0 ) plots -> fJetRes_Pt_20_50_Eta_25_5->Fill(BestGenJetMom.Pt()/JetMom.Pt());

       if( pt > 50.0 && pt < 100.0 && eta > 0.0 && eta < 2.5 ) plots -> fJetRes_Pt_50_100_Eta_0_25->Fill(BestGenJetMom.Pt()/JetMom.Pt());
       if( pt > 50.0 && pt < 100.0 && eta > 2.5 && eta < 5.0 ) plots -> fJetRes_Pt_50_100_Eta_25_5->Fill(BestGenJetMom.Pt()/JetMom.Pt());

       if( pt > 100.0 && pt < 200.0 && eta > 0.0 && eta < 2.5 ) plots -> fJetRes_Pt_100_200_Eta_0_25->Fill(BestGenJetMom.Pt()/JetMom.Pt());
       if( pt > 100.0 && pt < 200.0 && eta > 2.5 && eta < 5.0 ) plots -> fJetRes_Pt_100_200_Eta_25_5->Fill(BestGenJetMom.Pt()/JetMom.Pt());

       if( pt > 200.0 && pt < 500.0 && eta > 0.0 && eta < 2.5 ) plots -> fJetRes_Pt_200_500_Eta_0_25->Fill(BestGenJetMom.Pt()/JetMom.Pt());
       if( pt > 200.0 && pt < 500.0 && eta > 2.5 && eta < 5.0 ) plots -> fJetRes_Pt_200_500_Eta_25_5->Fill(BestGenJetMom.Pt()/JetMom.Pt());

       if( pt > 500.0               && eta > 0.0 && eta < 2.5 ) plots -> fJetRes_Pt_500_inf_Eta_0_25->Fill(BestGenJetMom.Pt()/JetMom.Pt());
       if( pt > 500.0               && eta > 2.5 && eta < 5.0 ) plots -> fJetRes_Pt_500_inf_Eta_25_5->Fill(BestGenJetMom.Pt()/JetMom.Pt());


     }


    }
  }
}


//------------------------------------------------------------------------------

void Example4(const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  result->Write(outputFile);

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
