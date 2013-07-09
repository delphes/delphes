/*
root -l examples/MemoryUsage.C\(\"ps_output.txt\"\)
*/

#include "Riostream.h"

static const Font_t kExRootFont = 42;
static const Float_t kExRootFontSize = 0.04;
static const Color_t kExRootBackgroundColor = 10;

//------------------------------------------------------------------------------

TGraph grvsz, grrss;
TLegend legend(0.41, 0.58, 0.77, 0.67);

TCanvas *canvas;

//------------------------------------------------------------------------------

void MemoryUsage(const char *inputFile)
{
  Int_t i, vsz, rss;
  ifstream in;

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

  grvsz.SetPoint(1, 0.0, 0.0);
  grrss.SetPoint(1, 0.0, 0.0);

  in.open(inputFile);
  i = 2;
  while(1)
  {
    in >> vsz  >> rss;
    if(!in.good()) break; 
    grvsz.SetPoint(i, (i-1)*0.1, vsz);
    grrss.SetPoint(i, (i-1)*0.1, rss);
    ++i;
  }

  grvsz.GetXaxis()->SetLimits(0, 30);
  grvsz.GetXaxis()->SetTitleOffset(1.5);
  grvsz.GetYaxis()->SetTitleOffset(1.75);
  grvsz.GetXaxis()->SetTitle("time, s");
  grvsz.GetYaxis()->SetTitle("memory usage, kB");
  grvsz.SetLineColor(kRed);
  grrss.SetLineColor(kMagenta);
  grvsz.SetLineWidth(2);
  grrss.SetLineWidth(2);
  grvsz.Draw("AL");
  grrss.Draw("L");

  legend.SetTextSize(kExRootFontSize);
  legend.SetTextFont(kExRootFont);
  legend.SetFillColor(kExRootBackgroundColor);
  legend.SetBorderSize(0);
  legend.AddEntry(&grvsz, "virtual memory size", "l");
  legend.AddEntry(&grrss, "resident set size", "l");
  legend.Draw();
}

