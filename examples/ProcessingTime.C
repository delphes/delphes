/*
root -l examples/ProcessingTime.C\(\"delphes_output.root\"\)
*/

//------------------------------------------------------------------------------

TGraphErrors gr;

//------------------------------------------------------------------------------

void ProcessingTime(const char *inputFile)
{
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  TH1F hist("time", "time", 50, 0, 0.01);
  Int_t i;

  for(i = 1; i < 8; ++i)
  {
    chain->Draw("Event.ProcTime >> time", TString::Format("Jet_size == %d", i+1));
    gr.SetPoint(i, i+1, hist.GetMean()*1000);
    gr.SetPointError(i, 0, hist.GetRMS()*1000);
  }

  gr.GetXaxis()->SetLimits(1.0, 9.0);
  gr.GetXaxis()->SetTitle("number of jets");
  gr.GetYaxis()->SetTitle("processing time per event, ms");
  gr.SetMarkerStyle(kFullDotMedium);
  gr.Draw("AP");
}

