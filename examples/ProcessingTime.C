/*
root -l examples/ProcessingTime.C\(\"delphes_output.root\"\)
*/

static const Font_t kExRootFont = 42;
static const Float_t kExRootFontSize = 0.04;
static const Color_t kExRootBackgroundColor = 10;

//------------------------------------------------------------------------------

TGraph gr;
TGraphErrors grerr;
TPaveText comment(0.20, 0.75, 0.50, 0.84, "brNDC");

TCanvas *canvas;

//------------------------------------------------------------------------------

void ProcessingTime(const char *inputFile)
{
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  TH1F hist("time", "time", 50, 0, 0.01);
  Int_t i;

  TDirectory *currentDirectory = gDirectory;

  // Graphics style parameters to avoid grey background on figures
  gStyle->SetCanvasColor(kExRootBackgroundColor);
  gStyle->SetStatColor(kExRootBackgroundColor);
  //  gStyle->SetTitleColor(kExRootBackgroundColor);
  gStyle->SetPadColor(kExRootBackgroundColor);

  gStyle->SetPadTopMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetStatFont(kExRootFont);
  gStyle->SetStatFontSize(kExRootFontSize);

  gStyle->SetTitleFont(kExRootFont, "");
  gStyle->SetTitleFont(kExRootFont, "X");
  gStyle->SetTitleFont(kExRootFont, "Y");
  gStyle->SetTitleFont(kExRootFont, "Z");
  gStyle->SetTitleSize(kExRootFontSize, "");
  gStyle->SetTitleSize(kExRootFontSize, "X");
  gStyle->SetTitleSize(kExRootFontSize, "Y");
  gStyle->SetTitleSize(kExRootFontSize, "Z");

  gStyle->SetLabelFont(kExRootFont, "X");
  gStyle->SetLabelFont(kExRootFont, "Y");
  gStyle->SetLabelFont(kExRootFont, "Z");
  gStyle->SetLabelSize(kExRootFontSize, "X");
  gStyle->SetLabelSize(kExRootFontSize, "Y");
  gStyle->SetLabelSize(kExRootFontSize, "Z");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTextFont(kExRootFont);
  gStyle->SetTextSize(kExRootFontSize);

  gStyle->SetOptStat(111110);
  // gStyle->SetOptFit(101);

  canvas = static_cast<TCanvas*>(gROOT->FindObject("c1"));
  if(canvas)
  {
    canvas->Clear();
    canvas->UseCurrentStyle();
    canvas->SetWindowSize(800, 650);
  }
  else
  {
    canvas = new TCanvas("c1", "c1", 800, 650);
  }
  canvas->SetGrid();
  canvas->SetHighLightColor(kExRootBackgroundColor);

  currentDirectory->cd();

  for(i = 0; i < 9; ++i)
  {
    chain->Draw("Event.ProcTime >> time", TString::Format("Jet_size == %d", i+2));
    gr.SetPoint(i, i+2, hist.GetMean()*1000);
    grerr.SetPoint(i, i+2, hist.GetMean()*1000);
    grerr.SetPointError(i, 0, hist.GetRMS()*1000);
  }

  grerr.GetXaxis()->SetLimits(1.0, 11.0);
  grerr.GetXaxis()->SetTitleOffset(1.5);
  grerr.GetYaxis()->SetTitleOffset(1.75);
  grerr.GetXaxis()->SetTitle("jet multiplicity");
  grerr.GetYaxis()->SetTitle("processing time per event, ms");
  gr.SetMarkerStyle(kFullCircle);
  gr.SetMarkerColor(kBlack);
  gr.SetMarkerSize(1);
  gr.SetLineColor(kBlack);
  gr.SetLineWidth(2);
  grerr.SetFillStyle(1001);
  grerr.SetFillColor(17);
  grerr.Draw("A3");
  gr.Draw("P");

  comment.SetTextSize(kExRootFontSize);
  comment.SetTextFont(kExRootFont);
  comment.SetTextAlign(22);
  comment.SetFillColor(kExRootBackgroundColor);
  comment.SetBorderSize(0);
  comment.AddText("ttbar + jets events");
  comment.Draw();
}

