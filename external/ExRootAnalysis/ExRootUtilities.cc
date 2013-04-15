
/** \class ExRootUtilities
 *
 *  Functions simplifying ROOT tree analysis
 *
 *  $Date: 2008-06-04 13:57:57 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "ExRootAnalysis/ExRootUtilities.h"

#include "TROOT.h"
#include "TH1.h"
#include "TChain.h"

#include <iostream>
#include <fstream>

using namespace std;

static const Font_t kExRootFont = 42;
static const Float_t kExRootFontSize = 0.04;

void HistStyle(TH1 *hist, Bool_t stats)
{
  hist->SetLineWidth(2);
  hist->SetLineColor(kBlack);
  hist->SetMarkerStyle(kFullSquare);
  hist->SetMarkerColor(kBlack);

  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetTitleOffset(1.75);
  hist->GetZaxis()->SetTitleOffset(1.5);

  hist->GetXaxis()->SetTitleFont(kExRootFont);
  hist->GetYaxis()->SetTitleFont(kExRootFont);
  hist->GetZaxis()->SetTitleFont(kExRootFont);
  hist->GetXaxis()->SetTitleSize(kExRootFontSize);
  hist->GetYaxis()->SetTitleSize(kExRootFontSize);
  hist->GetZaxis()->SetTitleSize(kExRootFontSize);

  hist->GetXaxis()->SetLabelFont(kExRootFont);
  hist->GetYaxis()->SetLabelFont(kExRootFont);
  hist->GetZaxis()->SetLabelFont(kExRootFont);
  hist->GetXaxis()->SetLabelSize(kExRootFontSize);
  hist->GetYaxis()->SetLabelSize(kExRootFontSize);
  hist->GetZaxis()->SetLabelSize(kExRootFontSize);

  hist->SetStats(stats);
}

//------------------------------------------------------------------------------

Bool_t FillChain(TChain *chain, const char *inputFileList)
{
  ifstream infile(inputFileList);
  string buffer;

  if(!infile.is_open())
  {
    cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << endl;
    return kFALSE;
  }

  while(1)
  {
    infile >> buffer;
    if(!infile.good()) break;
    chain->Add(buffer.c_str());
  }

  return kTRUE;
}

//------------------------------------------------------------------------------

