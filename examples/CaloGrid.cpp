//calorimeter grid
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include "ExRootAnalysis/ExRootConfReader.h"
#include "classes/DelphesClasses.h"
#include "display/Delphes3DGeometry.h"

#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"

using namespace std;

Bool_t debug = false;
//Bool_t debug = true;

int main(int argc, char *argv[])
{

  if(argc != 3)
  {
    cout << " Usage: ./CaloGrid [detector card] [calo name]" << endl;
    cout << "Example: ./CaloGrid cards/delphes_card_CMS.tcl ECal" << endl;
    return 0;
  }

  TString card(argv[1]);

  ExRootConfReader *confReader = new ExRootConfReader;
  confReader->ReadFile(card);

  std::vector<std::string> calorimeters_;
  std::map<std::string, std::set<std::pair<Double_t, Int_t> > > caloBinning_;

  std::string s(argv[2]);
  std::replace(s.begin(), s.end(), ',', ' ');
  std::istringstream stream(s);
  std::string word;
  while(stream >> word) calorimeters_.push_back(word);

  caloBinning_.clear(); // calo binning

  TCanvas c("", "", 1600, 838);

  gPad->SetLeftMargin(0.16);
  gPad->SetTopMargin(0.16);
  gPad->SetBottomMargin(0.20);
  gStyle->SetOptStat(0000000);

  gStyle->SetTextFont(132);

  TH2F h2("h2", "", 1, -6, 6, 1, 0, 6.28);

  h2.GetXaxis()->SetTitle("#eta");
  h2.GetYaxis()->SetTitle("#phi");

  h2.GetXaxis()->SetTitleFont(132);
  h2.GetYaxis()->SetTitleFont(132);
  h2.GetZaxis()->SetTitleFont(132);
  h2.GetXaxis()->SetLabelFont(132);
  h2.GetYaxis()->SetLabelFont(132);

  h2.GetXaxis()->SetTitleOffset(1.4);
  h2.GetYaxis()->SetTitleOffset(1.1);
  h2.GetXaxis()->SetLabelOffset(0.02);
  h2.GetYaxis()->SetLabelOffset(0.02);
  h2.GetXaxis()->SetTitleSize(0.06);
  h2.GetYaxis()->SetTitleSize(0.06);
  h2.GetXaxis()->SetLabelSize(0.06);
  h2.GetYaxis()->SetLabelSize(0.06);

  h2.GetXaxis()->SetTickLength(0.0);
  h2.GetYaxis()->SetTickLength(0.0);

  h2.Draw();

  // fake loop just keeping it for convenience right now
  for(std::vector<std::string>::const_iterator calo = calorimeters_.begin(); calo != calorimeters_.end(); ++calo)
  {

    //first entry is eta bin, second is number of phi bins
    set<pair<Double_t, Int_t> > caloBinning;
    ExRootConfParam paramEtaBins, paramPhiBins;
    ExRootConfParam param = confReader->GetParam(Form("%s::EtaPhiBins", calo->c_str()));
    Int_t size = param.GetSize();

    for(int i = 0; i < size / 2; ++i)
    {
      paramEtaBins = param[i * 2];
      paramPhiBins = param[i * 2 + 1];
      assert(paramEtaBins.GetSize() == 1);

      caloBinning.insert(std::make_pair(paramEtaBins[0].GetDouble(), paramPhiBins.GetSize() - 1));
    }
    caloBinning_[*calo] = caloBinning;

    TLine *liney;
    TLine *linex;

    //loop over calo binning
    std::set<std::pair<Double_t, Int_t> >::iterator it;

    Int_t n = -1;
    for(it = caloBinning.begin(); it != caloBinning.end(); ++it)
    {
      n++;

      if(debug) cout << "-----------------------" << endl;
      if(debug) cout << it->first << "," << it->second << endl;
      liney = new TLine(it->first, 0, it->first, 6.28);
      liney->SetLineColor(kRed + 3);
      liney->Draw();

      set<std::pair<Double_t, Int_t> >::iterator it2 = it;
      it2--;

      for(int j = 0; j <= it->second; j++)
      {

        Double_t yval0 = 0 + 6.28 * j / it->second;

        if(debug) cout << it2->first << "," << yval0 << "," << it->first << "," << yval0 << endl;

        linex = new TLine(it->first, yval0, it2->first, yval0);
        linex->SetLineColor(kRed + 3);
        linex->Draw();
      }
    }
  }

  TString text = TString(s);
  TText *th1 = new TText(5.00, 6.45, text);
  th1->SetTextAlign(31);
  th1->SetTextFont(132);
  th1->SetTextSize(0.075);
  th1->Draw();

  TString output = TString(s);
  c.Print(output + ".png", "png");
  c.Print(output + ".pdf", "pdf");
}
