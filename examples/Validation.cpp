/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <utility>
#include <vector>
#include <typeinfo>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"

#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

double ptrangemin = 10;
double ptrangemax = 10000;
static const int Nbins = 20;

int objStyle = 1;
int trackStyle = 7;
int towerStyle = 3;

Color_t objColor = kBlack;
Color_t trackColor = kBlack;
Color_t towerColor = kBlack;

double effLegXmin = 0.22;
double effLegXmax = 0.7;
double effLegYmin = 0.22;
double effLegYmax = 0.5;

double resLegXmin = 0.62;
double resLegXmax = 0.9;
double resLegYmin = 0.52;
double resLegYmax = 0.85;

double topLeftLegXmin = 0.22;
double topLeftLegXmax = 0.7;
double topLeftLegYmin = 0.52;
double topLeftLegYmax = 0.85;


struct resolPlot
{
  TH1 *cenResolHist;
  TH1 *fwdResolHist;
  double ptmin;
  double ptmax;
  double xmin;
  double xmax;
  TString obj;

  resolPlot();
  resolPlot(double ptdown, double ptup, TString object);
  void set(double ptdown, double ptup, TString object, double xmin = 0, double xmax = 2);
  void print(){std::cout << ptmin << std::endl;}
};


resolPlot::resolPlot()
{
}

resolPlot::resolPlot(double ptdown, double ptup, TString object)
{
  this->set(ptdown,ptup,object);
}

void resolPlot::set(double ptdown, double ptup, TString object, double xmin, double xmax)
{
  ptmin = ptdown;
  ptmax = ptup;
  obj = object;

  cenResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", 200,  xmin, xmax);
  fwdResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", 200,  xmin, xmax);
}

void HistogramsCollection(std::vector<resolPlot> *histos, double ptmin, double ptmax, TString obj, double xmin = 0, double xmax = 2)
{
  double width;
  double ptdown;
  double ptup;
  resolPlot ptemp;

  for (int i = 0; i < Nbins; i++)
  {
    width = (ptmax - ptmin) / Nbins;
    ptdown = TMath::Power(10,ptmin + i * width );
    ptup = TMath::Power(10,ptmin + (i+1) * width );
    ptemp.set(ptdown, ptup, obj, xmin, xmax);
    histos->push_back(ptemp);
  }
}

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BinLogX(TH1*h)
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++)
  {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}


//------------------------------------------------------------------------------

template<typename T>
std::pair<TH1D*, TH1D*> GetEff(TClonesArray *branchReco, TClonesArray *branchParticle, TString name, int pdgID, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing "<< branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D *histGenPtcen = new TH1D(name+" gen spectra Pt",name+" gen spectra cen", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D *histRecoPtcen = new TH1D(name+" reco spectra Pt",name+" reco spectra cen", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D *histGenPtfwd  = new TH1D(name+" gen spectra Eta",name+" gen spectra fwd", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D *histRecoPtfwd = new TH1D(name+" reco spectra Eta",name+" reco spectra fwd", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));

  histGenPtcen->SetDirectory(0);
  histRecoPtcen->SetDirectory(0);
  histGenPtfwd->SetDirectory(0);
  histRecoPtfwd->SetDirectory(0);

  BinLogX(histGenPtcen);
  BinLogX(histRecoPtcen);
  BinLogX(histGenPtfwd);
  BinLogX(histRecoPtfwd);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all generated particle in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle*) branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;

      if (particle->PID == pdgID && genMomentum.Pt() > ptrangemin && genMomentum.Pt() < ptrangemax )
      {

        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T*)branchReco->At(j);
          recoMomentum = recoObj->P4();
          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
        //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;

          // take the closest parton candidate
          if(TMath::Abs(pdgID) == 5)
          {
            Jet *jet = (Jet *)recoObj;
            if(jet->BTag != 1) continue;
          }
          if(TMath::Abs(pdgID) == 15)
          {
            Jet *jet = (Jet *)recoObj;
            if(jet->TauTag != 1) continue;
          }
          if(genMomentum.DeltaR(recoMomentum) < deltaR)
          {
            deltaR = genMomentum.DeltaR(recoMomentum);
            bestRecoMomentum = recoMomentum;
          }
        }

        pt  = genMomentum.Pt();
        eta = genMomentum.Eta();

        if (TMath::Abs(eta) < 1.5)
        {
          histGenPtcen->Fill(pt);
          if(deltaR < 0.3) { histRecoPtcen->Fill(pt); }
        }
        else if (TMath::Abs(eta) < 2.5)
        {
          histGenPtfwd->Fill(pt);
          if(deltaR < 0.3) { histRecoPtfwd->Fill(pt); }

        }
      }
    }
  }


  std::pair<TH1D*,TH1D*> histos;

  histRecoPtcen->Divide(histGenPtcen);
  histRecoPtfwd->Divide(histGenPtfwd);

  histos.first = histRecoPtcen;
  histos.second = histRecoPtfwd;

  return histos;
}

template<typename T>
void GetEres(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T* recoObj;

  TLorentzVector recoMomentum, genMomentum, bestGenMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j, bin;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchReco->GetEntriesFast(); ++i)
    {
      recoObj = (T*) branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

     // Loop over all hard partons in event
     for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
     {
        particle = (GenParticle*) branchParticle->At(j);
        if (particle->PID == pdgID && particle->Status == 1)
        {
          genMomentum = particle->P4();

          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
          if(genMomentum.Px() == 0 && genMomentum.Py() == 0) continue;

          // take the closest parton candidate
          if(genMomentum.DeltaR(recoMomentum) < deltaR)
          {
             deltaR = genMomentum.DeltaR(recoMomentum);
             bestGenMomentum = genMomentum;
          }
        }
      }

      if(deltaR < 0.3)
      {
        pt  = bestGenMomentum.E();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
          if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta < 2.5)
          {
            if (eta < 1.5) {histos->at(bin).cenResolHist->Fill(recoMomentum.E()/bestGenMomentum.E());}
            else if (eta < 2.5) {histos->at(bin).fwdResolHist->Fill(recoMomentum.E()/bestGenMomentum.E());}
          }
        }
      }
    }
  }
}


template<typename T>
void GetPtres(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing pt resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T* recoObj;

  TLorentzVector recoMomentum, genMomentum, bestGenMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j, bin;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchReco->GetEntriesFast(); ++i)
    {
      recoObj = (T*) branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

     // Loop over all hard partons in event
     for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
     {
        particle = (GenParticle*) branchParticle->At(j);
        if (particle->PID == pdgID && particle->Status == 1)
        {
          genMomentum = particle->P4();

          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
          if(genMomentum.Px() == 0 && genMomentum.Py() == 0) continue;

          // take the closest parton candidate
          if(genMomentum.DeltaR(recoMomentum) < deltaR)
          {
             deltaR = genMomentum.DeltaR(recoMomentum);
             bestGenMomentum = genMomentum;
          }
        }
      }

      if(deltaR < 0.3)
      {
        pt  = bestGenMomentum.Pt();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
          if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta < 2.5)
          {
            if (eta < 1.5) {histos->at(bin).cenResolHist->Fill(recoMomentum.Pt()/bestGenMomentum.Pt());}
            else if (eta < 2.5) {histos->at(bin).fwdResolHist->Fill(recoMomentum.Pt()/bestGenMomentum.Pt());}
          }
        }
      }
    }
  }
}


void GetJetsEres(std::vector<resolPlot> *histos, TClonesArray *branchJet, TClonesArray *branchGenJet, ExRootTreeReader *treeReader)
{

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing resolution of " << branchJet->GetName() << " induced by " << branchGenJet->GetName() << endl;

  Jet *jet, *genjet;

  TLorentzVector jetMomentum, genJetMomentum, bestGenJetMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j, bin;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%10000 == 0) cout << "Event number: "<< entry <<endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < TMath::Min(2,branchJet->GetEntriesFast()); ++i) //branchJet->GetEntriesFast(); ++i)
    {

      jet = (Jet*) branchJet->At(i);
      jetMomentum = jet->P4();

      deltaR = 999;

     // Loop over all hard partons in event
     for(j = 0; j < TMath::Min(2,branchGenJet->GetEntriesFast()); ++j)
     {
        genjet = (Jet*) branchGenJet->At(j);

        genJetMomentum = genjet->P4();

        // this is simply to avoid warnings from initial state particle
        // having infite rapidity ...
        if(genJetMomentum.Px() == 0 && genJetMomentum.Py() == 0) continue;

        // take the closest parton candidate
        if(genJetMomentum.DeltaR(jetMomentum) < deltaR)
        {
           deltaR = genJetMomentum.DeltaR(jetMomentum);
           bestGenJetMomentum = genJetMomentum;
        }
      }

      if(deltaR < 0.25)
      {
        pt  = genJetMomentum.E();
        eta = TMath::Abs(genJetMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
            if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta < 1.5)
            {
                histos->at(bin).cenResolHist->Fill(jetMomentum.E()/bestGenJetMomentum.E());
            }
            else if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta < 2.5)
            {
                histos->at(bin).fwdResolHist->Fill(jetMomentum.E()/bestGenJetMomentum.E());
            }
        }
      }
    }
  }
}

void GetMetres(std::vector<resolPlot> *histos, TClonesArray *branchScalarHT, TClonesArray *branchMet, TClonesArray *branchJet, ExRootTreeReader *treeReader)
{

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing resolution of " << branchMet->GetName() << " vs " << branchScalarHT->GetName() << endl;

  MissingET *met;
  ScalarHT *scalarHT;

  Long64_t entry;

  Int_t bin;
  Double_t ht;

  Jet *jet;
  TLorentzVector p1, p2;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%10000 == 0) cout << "Event number: "<< entry <<endl;

    if (branchJet->GetEntriesFast() > 1)
    {

      jet = (Jet*) branchJet->At(0);
      p1 = jet->P4();
      jet = (Jet*) branchJet->At(1);
      p2 = jet->P4();

      met = (MissingET*) branchMet->At(0);
      scalarHT = (ScalarHT*) branchScalarHT->At(0);
      ht = scalarHT->HT;

      if(p1.Pt() < 0.75*ht/2) continue;
      if(p2.Pt() < 0.75*ht/2) continue;

      for (bin = 0; bin < Nbins; bin++)
      {
        if(ht > histos->at(bin).ptmin && ht < histos->at(bin).ptmax )
        {
          histos->at(bin).cenResolHist->Fill(met->P4().Px());
          histos->at(bin).fwdResolHist->Fill(met->P4().Py());
        }
      }
    }
  }
}


std::pair<Double_t, Double_t> GausFit(TH1* hist)
{
  TF1 *f1 = new TF1("f1", "gaus", hist->GetMean()-2*hist->GetRMS(), hist->GetMean()+2*hist->GetRMS());
  hist->Fit("f1","RQ");

  TF1 *f2 = new TF1("f2", "gaus", f1->GetParameter(1) - 2*f1->GetParameter(2), f1->GetParameter(1) + 2*f1->GetParameter(2));
  hist->Fit("f2","RQ");

  Double_t sig = f2->GetParameter(2);
  Double_t sigErr = f2->GetParError(2);

  delete f1;
  delete f2;
  return make_pair (sig, sigErr);
}


TGraphErrors EresGraph(std::vector<resolPlot> *histos, bool central, bool rms = false)
{
  Int_t bin;
  Int_t count = 0;
  TGraphErrors gr = TGraphErrors(Nbins/2);
  Double_t sig = 0;
  Double_t sigErr = 0;
  for (bin = 0; bin < Nbins; bin++)
  {
    if (central == true && histos->at(bin).cenResolHist->GetEntries() > 100)
    {
      std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).cenResolHist);
      if (rms == true)
      {
        gr.SetPoint(count,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sigvalues.second);
        gr.SetPointError(count,0, sigvalues.second); // to correct
      }
      else
      {
        gr.SetPoint(count,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sigvalues.first);
        gr.SetPointError(count,0, sigvalues.second);
      }
      count++;
    }

    else if (central == false && histos->at(bin).fwdResolHist->GetEntries() > 10)
    {
      std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).fwdResolHist);
      if (rms == true)
      {
        gr.SetPoint(count,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sigvalues.second);
        gr.SetPointError(count,0, sigvalues.second); // to correct
      }
      else
      {
        gr.SetPoint(count,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sigvalues.first);
        gr.SetPointError(count,0, sigvalues.second);
      }
      count++;
    }

  }
  return gr;
}


//------------------------------------------------------------------------------


// type 1 : object, 2 : track, 3 : tower

void addGraph(TMultiGraph *mg, TGraphErrors *gr, TLegend *leg, int type)
{

  gr->SetLineWidth(2);

  switch ( type )
  {
    case 1:
      gr->SetLineColor(objColor);
      gr->SetLineStyle(objStyle);
      std::cout << "Adding " << gr->GetName() << std::endl;
      mg->Add(gr);
      leg->AddEntry(gr,"Reco","l");
      break;

    case 2:
      gr->SetLineColor(trackColor);
      gr->SetLineStyle(trackStyle);
      mg->Add(gr);
      leg->AddEntry(gr,"Track","l");
      break;

    case 3:
      gr->SetLineColor(towerColor);
      gr->SetLineStyle(towerStyle);
      mg->Add(gr);
      leg->AddEntry(gr,"Tower","l");
      break;

    case 0:
      gr->SetLineColor(objColor);
      gr->SetLineStyle(objStyle);
      mg->Add(gr);
      break;

    default:
      std::cout << "wrong type, possibles choices are Object, Track and Tower" << std::endl;
      break;
  }
}

void addHist(TH1D *h, TLegend *leg, int type)
{
  h->SetLineWidth(2);

  switch ( type )
  {
    case 1:
      h->SetLineColor(objColor);
      h->SetLineStyle(objStyle);
      leg->AddEntry(h,"Reco","l");
      break;

    case 2:
      h->SetLineColor(trackColor);
      h->SetLineStyle(trackStyle);
      leg->AddEntry(h,"Track","l");
      break;

    case 3:
      h->SetLineColor(towerColor);
      h->SetLineStyle(towerStyle);
      leg->AddEntry(h,"Tower","l");
      break;

    case 0:
      h->SetLineColor(objColor);
      h->SetLineStyle(objStyle);
      break;

    default:
      std::cout << "wrong type, possibles choices are Object, Track and Tower" << std::endl;
      break;
  }
}

void DrawAxis(TMultiGraph *mg, TLegend *leg, double max, int type = 0)
{
  mg->SetMinimum(0.);
  mg->SetMaximum(max);
  mg->GetXaxis()->SetLimits(ptrangemin,ptrangemax);
  mg->GetYaxis()->SetTitle("resolution");
  if (type == 0) mg->GetXaxis()->SetTitle("E [GeV]");
  else mg->GetXaxis()->SetTitle("p_{T} [GeV]");
  mg->GetYaxis()->SetTitleSize(0.07);
  mg->GetXaxis()->SetTitleSize(0.07);
  mg->GetYaxis()->SetLabelSize(0.06);
  mg->GetXaxis()->SetLabelSize(0.06);
  mg->GetYaxis()->SetLabelOffset(0.03);
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetXaxis()->SetTitleOffset(1.4);

  mg->GetYaxis()->SetNdivisions(505);

  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  gStyle->SetOptTitle(0);
  gPad->SetLogx();
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->Modified();
  gPad->Update();

}

void DrawAxis(TH1D *h, TLegend *leg, int type = 0)
{

  h->GetYaxis()->SetRangeUser(0,1.0);
  if (type == 0) h->GetXaxis()->SetTitle("E [GeV]");
  else h->GetXaxis()->SetTitle("p_{T} [GeV]");
  h->GetYaxis()->SetTitle("efficiency");
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetLabelSize(0.06);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelOffset(0.03);
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetTitleOffset(1.4);

  h->GetYaxis()->SetNdivisions(505);

  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);

  gPad->Modified();
  gPad->Update();

}


void Validation(const char *inputFileElectron, const char *inputFileMuon, const char *inputFilePhoton, const char *inputFileJet, const char *inputFileBJet, const char *inputFileTauJet, const char *outputFile)
{
  TChain *chainElectron = new TChain("Delphes");
  chainElectron->Add(inputFileElectron);
  ExRootTreeReader *treeReaderElectron = new ExRootTreeReader(chainElectron);

  TChain *chainMuon = new TChain("Delphes");
  chainMuon->Add(inputFileMuon);
  ExRootTreeReader *treeReaderMuon = new ExRootTreeReader(chainMuon);

  TChain *chainPhoton = new TChain("Delphes");
  chainPhoton->Add(inputFilePhoton);
  ExRootTreeReader *treeReaderPhoton = new ExRootTreeReader(chainPhoton);

  TChain *chainJet = new TChain("Delphes");
  chainJet->Add(inputFileJet);
  ExRootTreeReader *treeReaderJet = new ExRootTreeReader(chainJet);

  TChain *chainBJet = new TChain("Delphes");
  chainBJet->Add(inputFileBJet);
  ExRootTreeReader *treeReaderBJet = new ExRootTreeReader(chainBJet);

  TChain *chainTauJet = new TChain("Delphes");
  chainTauJet->Add(inputFileTauJet);
  ExRootTreeReader *treeReaderTauJet = new ExRootTreeReader(chainTauJet);

  TClonesArray *branchParticleElectron = treeReaderElectron->UseBranch("Particle");
  TClonesArray *branchTrackElectron = treeReaderElectron->UseBranch("Track");
  TClonesArray *branchTowerElectron = treeReaderElectron->UseBranch("Tower");
  TClonesArray *branchElectron = treeReaderElectron->UseBranch("Electron");

  TClonesArray *branchParticleMuon = treeReaderMuon->UseBranch("Particle");
  TClonesArray *branchTrackMuon = treeReaderMuon->UseBranch("Track");
  TClonesArray *branchMuon = treeReaderMuon->UseBranch("Muon");

  TClonesArray *branchParticlePhoton = treeReaderPhoton->UseBranch("Particle");
  TClonesArray *branchTowerPhoton = treeReaderPhoton->UseBranch("Tower");
  TClonesArray *branchPhoton = treeReaderPhoton->UseBranch("Photon");

  TClonesArray *branchGenJet = treeReaderJet->UseBranch("GenJet");
  TClonesArray *branchPFJet = treeReaderJet->UseBranch("Jet");
  TClonesArray *branchCaloJet = treeReaderJet->UseBranch("CaloJet");

  TClonesArray *branchParticleBJet = treeReaderBJet->UseBranch("Particle");
  TClonesArray *branchPFBJet = treeReaderBJet->UseBranch("Jet");

  TClonesArray *branchParticleTauJet = treeReaderTauJet->UseBranch("Particle");
  TClonesArray *branchPFTauJet = treeReaderTauJet->UseBranch("Jet");

  TClonesArray *branchScalarHT = treeReaderJet->UseBranch("ScalarHT");
  TClonesArray *branchMet = treeReaderJet->UseBranch("MissingET");

  ///////////////
  // Electrons //
  ///////////////

  // Reconstruction efficiency
  TString elecs = "Electron";
  int elID = 11;
  std::pair<TH1D*,TH1D*> histos_el = GetEff<Electron>(branchElectron, branchParticleElectron, "Electron", elID, treeReaderElectron);

  // tracking reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_eltrack = GetEff<Track>(branchTrackElectron, branchParticleElectron, "electronTrack", elID, treeReaderElectron);

  // Tower reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_eltower = GetEff<Tower>(branchTowerElectron, branchParticleElectron, "electronTower", elID, treeReaderElectron);

  // Electron Energy Resolution
  std::vector<resolPlot> plots_el;
  HistogramsCollection(&plots_el, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electrons");
  GetEres<Electron>(&plots_el, branchElectron, branchParticleElectron, elID, treeReaderElectron);
  TGraphErrors gr_el = EresGraph(&plots_el, true);
  TGraphErrors gr_elFwd = EresGraph(&plots_el, false);
  gr_el.SetName("Electron");
  gr_elFwd.SetName("ElectronFwd");

  // Electron Track Energy Resolution
  std::vector<resolPlot> plots_eltrack;
  HistogramsCollection(&plots_eltrack, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTracks");
  GetEres<Track>(&plots_eltrack, branchTrackElectron, branchParticleElectron, elID, treeReaderElectron);
  TGraphErrors gr_eltrack = EresGraph(&plots_eltrack, true);
  TGraphErrors gr_eltrackFwd = EresGraph(&plots_eltrack, false);
  gr_eltrack.SetName("ElectronTracks");
  gr_eltrackFwd.SetName("ElectronTracksFwd");

  // Electron Tower Energy Resolution
  std::vector<resolPlot> plots_eltower;
  HistogramsCollection(&plots_eltower, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTower");
  GetEres<Tower>(&plots_eltower, branchTowerElectron, branchParticleElectron, elID, treeReaderElectron);
  TGraphErrors gr_eltower = EresGraph(&plots_eltower, true);
  TGraphErrors gr_eltowerFwd = EresGraph(&plots_eltower, false);
  gr_eltower.SetName("ElectronTower");
  gr_eltrackFwd.SetName("ElectronTracksFwd");

  // Canvases
  TString elEff = "electronEff";
  TCanvas *C_el1 = new TCanvas(elEff,elEff, 1600, 600);
  C_el1->Divide(2);
  C_el1->cd(1);
  TLegend *leg_el1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_el1->SetHeader("#splitline{electrons}{|#eta| < 1.5}");
  leg_el1->AddEntry("","","");

  gPad->SetLogx();
  histos_eltrack.first->Draw("][");
  addHist(histos_eltrack.first, leg_el1, 2);
  histos_el.first->Draw("same ][");
  addHist(histos_el.first, leg_el1, 1);
  DrawAxis(histos_eltrack.first, leg_el1,1);

  leg_el1->Draw();

  C_el1->cd(2);
  TLegend *leg_el2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_el2->SetHeader("#splitline{electrons}{1.5 < |#eta| < 2.5}");
  leg_el2->AddEntry("","","");

  gPad->SetLogx();
  histos_eltrack.second->Draw("][");
  addHist(histos_eltrack.second, leg_el2, 2);
  histos_el.second->Draw("same ][");
  addHist(histos_el.second, leg_el2, 1);

  DrawAxis(histos_eltrack.second, leg_el2, 1);
  leg_el2->Draw();

  TString elRes = "electronERes";
  TString elResFwd = "electronEResForward";
  TCanvas *C_el2 = new TCanvas(elRes,elRes, 1600, 600);
  C_el2->Divide(2);
  C_el2->cd(1);
  TMultiGraph *mg_el = new TMultiGraph(elRes,elRes);
  TLegend *leg_el = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_el->SetHeader("#splitline{electrons}{|#eta| < 1.5}");
  leg_el->AddEntry("","","");

  addGraph(mg_el, &gr_eltower, leg_el, 3);
  addGraph(mg_el, &gr_eltrack, leg_el, 2);
  addGraph(mg_el, &gr_el, leg_el, 1);

  mg_el->Draw("ACX");
  leg_el->Draw();

  DrawAxis(mg_el, leg_el, 0.1);

  C_el2->cd(2);
  TMultiGraph *mg_elFwd = new TMultiGraph(elResFwd,elResFwd);
  TLegend *leg_elFwd = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_elFwd->SetHeader("#splitline{electrons}{1.5 < |#eta| < 2.5}");
  leg_elFwd->AddEntry("","","");

  addGraph(mg_elFwd, &gr_eltowerFwd, leg_elFwd, 3);
  addGraph(mg_elFwd, &gr_eltrackFwd, leg_elFwd, 2);
  addGraph(mg_elFwd, &gr_elFwd, leg_elFwd, 1);

  mg_elFwd->Draw("ACX");
  leg_elFwd->Draw();

  DrawAxis(mg_elFwd, leg_elFwd, 0.2);

  TString pdfOutput(outputFile);
  pdfOutput.ReplaceAll(".root", ".pdf");

  C_el1->Print(pdfOutput+"(","pdf");
  C_el2->Print(pdfOutput,"pdf");

  gDirectory->cd(0);

  ///////////
  // Muons //
  ///////////

  // Reconstruction efficiency
  int muID = 13;
  std::pair<TH1D*,TH1D*> histos_mu = GetEff<Muon>(branchMuon, branchParticleMuon,"Muon", muID, treeReaderMuon);

  // muon tracking reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_mutrack = GetEff<Track>(branchTrackMuon, branchParticleMuon, "muonTrack", muID, treeReaderMuon);

  // Muon Pt Resolution
  std::vector<resolPlot> plots_mu;
  HistogramsCollection(&plots_mu, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "muons");
  GetPtres<Muon>(&plots_mu, branchMuon, branchParticleMuon, muID, treeReaderMuon);
  TGraphErrors gr_mu = EresGraph(&plots_mu, true);
  TGraphErrors gr_muFwd = EresGraph(&plots_mu, false);
  gr_mu.SetName("Muon");
  gr_muFwd.SetName("MuonFwd");

  // Muon Track Energy Resolution
  std::vector<resolPlot> plots_mutrack;
  HistogramsCollection(&plots_mutrack, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "muonsTracks");
  GetPtres<Track>(&plots_mutrack, branchTrackMuon, branchParticleMuon, muID, treeReaderMuon);
  TGraphErrors gr_mutrack = EresGraph(&plots_mutrack, true);
  TGraphErrors gr_mutrackFwd = EresGraph(&plots_mutrack, false);
  gr_mutrackFwd.SetName("MuonTracksFwd");

  // Canvas

  TString muEff = "muonEff";
  TCanvas *C_mu1 = new TCanvas(muEff,muEff, 1600, 600);
  C_mu1->Divide(2);
  C_mu1->cd(1);
  TLegend *leg_mu1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_mu1->SetHeader("#splitline{muons}{|#eta| < 1.5}");
  leg_mu1->AddEntry("","","");


  gPad->SetLogx();
  histos_mutrack.first->Draw("][");
  addHist(histos_mutrack.first, leg_mu1, 2);
  histos_mu.first->Draw("same ][");
  addHist(histos_mu.first, leg_mu1, 1);

  DrawAxis(histos_mutrack.first, leg_mu1, 1);

  leg_mu1->Draw();

  C_mu1->cd(2);
  TLegend *leg_mu2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_mu2->SetHeader("#splitline{muons}{1.5 < |#eta| < 2.5}");
  leg_mu2->AddEntry("","","");

  gPad->SetLogx();
  histos_mutrack.second->Draw("][");
  addHist(histos_mutrack.second, leg_mu2, 2);
  histos_mu.second->Draw("same ][");
  addHist(histos_mu.second, leg_mu2, 1);

  DrawAxis(histos_mutrack.second, leg_mu2, 1);
  leg_mu2->Draw();

  TString muRes = "muonERes";
  TString muResFwd = "muonEResFwd";

  TCanvas *C_mu = new TCanvas(muRes,muRes, 1600, 600);
  C_mu->Divide(2);
  C_mu->cd(1);
  TMultiGraph *mg_mu = new TMultiGraph(muRes,muRes);
  TLegend *leg_mu = new TLegend(topLeftLegXmin,topLeftLegYmin,topLeftLegXmax,topLeftLegYmax);
  leg_mu->SetHeader("#splitline{muons}{|#eta| < 1.5}");
  leg_mu->AddEntry("","","");

  addGraph(mg_mu, &gr_mutrack, leg_mu, 2);
  addGraph(mg_mu, &gr_mu, leg_mu, 1);

  mg_mu->Draw("ACX");
  leg_mu->Draw();

  DrawAxis(mg_mu, leg_mu, 0.3, 1);

  C_mu->cd(2);
  TMultiGraph *mg_muFwd = new TMultiGraph(muResFwd,muResFwd);
  TLegend *leg_muFwd = new TLegend(topLeftLegXmin,topLeftLegYmin,topLeftLegXmax,topLeftLegYmax);
  leg_muFwd->SetHeader("#splitline{muons}{1.5 < |#eta| < 2.5}");
  leg_muFwd->AddEntry("","","");

  addGraph(mg_muFwd, &gr_mutrackFwd, leg_muFwd, 2);
  addGraph(mg_muFwd, &gr_muFwd, leg_muFwd, 1);

  mg_muFwd->Draw("ACX");
  leg_muFwd->Draw();

  DrawAxis(mg_muFwd, leg_muFwd, 0.3, 1);

  C_mu1->Print(pdfOutput,"pdf");
  C_mu->Print(pdfOutput,"pdf");

  gDirectory->cd(0);

  /////////////
  // Photons //
  /////////////

  // Reconstruction efficiency
  int phID = 22;
  std::pair<TH1D*,TH1D*> histos_ph = GetEff<Electron>(branchPhoton, branchParticlePhoton, "Photon", phID, treeReaderPhoton);
  std::pair<TH1D*,TH1D*> histos_phtower = GetEff<Electron>(branchTowerPhoton, branchParticlePhoton, "Photon", phID, treeReaderPhoton);

  // Photon Energy Resolution
  std::vector<resolPlot> plots_ph;
  HistogramsCollection(&plots_ph, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "photons");
  GetEres<Photon>(&plots_ph, branchPhoton, branchParticlePhoton, phID, treeReaderPhoton);
  TGraphErrors gr_ph = EresGraph(&plots_ph, true);
  TGraphErrors gr_phFwd = EresGraph(&plots_ph, false);
  gr_ph.SetName("Photon");
  gr_phFwd.SetName("PhotonFwd");


  // Photon Tower Energy Resolution
  std::vector<resolPlot> plots_phtower;
  HistogramsCollection(&plots_phtower, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "photonsTower");
  GetEres<Tower>(&plots_phtower, branchTowerPhoton, branchParticlePhoton, phID, treeReaderPhoton);
  TGraphErrors gr_phtower = EresGraph(&plots_phtower, true);
  TGraphErrors gr_phtowerFwd = EresGraph(&plots_phtower, false);
  gr_phtower.SetName("PhotonTower");
  gr_phtowerFwd.SetName("PhotonTowerFwd");

  // Canvas

  TString phEff = "photonEff";
  TCanvas *C_ph1 = new TCanvas(phEff,phEff, 1600, 600);
  C_ph1->Divide(2);
  C_ph1->cd(1);
  TLegend *leg_ph1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_ph1->SetHeader("#splitline{photons}{|#eta| < 1.5}");
  leg_ph1->AddEntry("","","");


  gPad->SetLogx();
  histos_phtower.first->Draw("][");
  addHist(histos_phtower.first, leg_ph1, 3);
  histos_ph.first->Draw("same ][");
  addHist(histos_ph.first, leg_ph1, 1);

  DrawAxis(histos_phtower.first, leg_ph1, 1);
  leg_ph1->Draw();

  C_ph1->cd(2);
  TLegend *leg_ph2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);
  leg_ph2->SetHeader("#splitline{photons}{1.5 < |#eta| < 2.5}");
  leg_ph2->AddEntry("","","");


  gPad->SetLogx();
  histos_phtower.second->Draw("][");
  addHist(histos_phtower.second, leg_ph2, 3);
  histos_ph.second->Draw("same ][");
  addHist(histos_ph.second, leg_ph2, 1);

  DrawAxis(histos_phtower.second, leg_ph2, 1);
  leg_ph2->Draw();

  TString phRes = "phERes";
  TString phResFwd = "phEResFwd";

  TCanvas *C_ph = new TCanvas(phRes,phRes, 1600, 600);
  C_ph->Divide(2);
  C_ph->cd(1);
  TMultiGraph *mg_ph = new TMultiGraph(phRes,phRes);
  TLegend *leg_ph = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_ph->SetHeader("#splitline{photons}{|#eta| < 1.5}");
  leg_ph->AddEntry("","","");

  addGraph(mg_ph, &gr_phtower, leg_ph, 3);
  addGraph(mg_ph, &gr_ph, leg_ph, 1);

  mg_ph->Draw("ACX");
  leg_ph->Draw();

  DrawAxis(mg_ph, leg_ph, 0.1);

  C_ph->cd(2);
  TMultiGraph *mg_phFwd = new TMultiGraph(phResFwd,phResFwd);
  TLegend *leg_phFwd = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_phFwd->SetHeader("#splitline{photons}{1.5 < |#eta| < 2.5}");
  leg_phFwd->AddEntry("","","");

  addGraph(mg_phFwd, &gr_phtowerFwd, leg_phFwd, 3);
  addGraph(mg_phFwd, &gr_phFwd, leg_phFwd, 1);

  mg_phFwd->Draw("ACX");
  leg_phFwd->Draw();

  DrawAxis(mg_phFwd, leg_phFwd, 0.1);

  C_ph1->Print(pdfOutput,"pdf");
  C_ph->Print(pdfOutput,"pdf");

  gDirectory->cd(0);

  //////////
  // Jets //
  //////////

  // BJets Reconstruction efficiency
  int bID = 5;
  std::pair<TH1D*,TH1D*> histos_btag = GetEff<Jet>(branchPFBJet, branchParticleBJet,"BTag", bID, treeReaderBJet);

  // TauJets Reconstruction efficiency
  int tauID = 15;
  std::pair<TH1D*,TH1D*> histos_tautag = GetEff<Jet>(branchPFTauJet, branchParticleTauJet,"TauTag", tauID, treeReaderTauJet);

  // PFJets Energy Resolution
  std::vector<resolPlot> plots_pfjets;
  HistogramsCollection(&plots_pfjets, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "PFJet");
  GetJetsEres(&plots_pfjets, branchPFJet, branchGenJet, treeReaderJet);
  TGraphErrors gr_pfjets = EresGraph(&plots_pfjets, true);
  TGraphErrors gr_pfjetsFwd = EresGraph(&plots_pfjets, false);
  gr_pfjets.SetName("pfJet");
  gr_pfjetsFwd.SetName("pfJetFwd");

  // CaloJets Energy Resolution
  std::vector<resolPlot> plots_calojets;
  HistogramsCollection(&plots_calojets, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "CaloJet");
  GetJetsEres(&plots_calojets, branchCaloJet, branchGenJet, treeReaderJet);
  TGraphErrors gr_calojets = EresGraph(&plots_calojets, true);
  TGraphErrors gr_calojetsFwd = EresGraph(&plots_calojets, false);
  gr_calojets.SetName("caloJet");
  gr_calojetsFwd.SetName("caloJetFwd");

  // MET Resolution vs HT
  std::vector<resolPlot> plots_met;
  HistogramsCollection(&plots_met, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "MET", -500, 500);
  GetMetres(&plots_met, branchScalarHT, branchMet, branchPFJet, treeReaderJet);
  TGraphErrors gr_met = EresGraph(&plots_met, true);
  gr_calojets.SetName("MET");

  // Canvas
  TString btagEff = "btagEff";
  TCanvas *C_btag1 = new TCanvas(btagEff,btagEff, 1600, 600);
  C_btag1->Divide(2);
  C_btag1->cd(1);
  TLegend *leg_btag1 = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_btag1->SetHeader("#splitline{B-tagging}{|#eta| < 1.5}");
  leg_btag1->AddEntry("","","");

  gPad->SetLogx();
  histos_btag.first->Draw();
  addHist(histos_btag.first, leg_btag1, 0);

  DrawAxis(histos_btag.first, leg_btag1, 1);
  leg_btag1->Draw();

  C_btag1->cd(2);
  TLegend *leg_btag2 = new TLegend(resLegXmin,resLegYmin,resLegXmax+0.05,resLegYmax);
  leg_btag2->SetHeader("#splitline{B-tagging}{1.5 < |#eta| < 2.5}");
  leg_btag2->AddEntry("","","");

  gPad->SetLogx();
  histos_btag.second->Draw();
  addHist(histos_btag.second, leg_btag2, 0);

  DrawAxis(histos_btag.second, leg_btag2, 1);
  leg_btag2->Draw();
  C_btag1->cd(0);

  TString tautagEff = "tautagEff";
  TCanvas *C_tautag1 = new TCanvas(tautagEff,tautagEff, 1600, 600);
  C_tautag1->Divide(2);
  C_tautag1->cd(1);
  TLegend *leg_tautag1 = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_tautag1->SetHeader("#splitline{#tau-tagging}{|#eta| < 1.5}");
  leg_tautag1->AddEntry("","","");

  gPad->SetLogx();
  histos_tautag.first->Draw();
  addHist(histos_tautag.first, leg_tautag1, 0);

  DrawAxis(histos_tautag.first, leg_tautag1, 1);
  leg_tautag1->Draw();

  C_tautag1->cd(2);
  TLegend *leg_tautag2 = new TLegend(resLegXmin,resLegYmin,resLegXmax+0.05,resLegYmax);
  leg_tautag2->SetHeader("#splitline{#tau-tagging}{1.5 < |#eta| < 2.5}");
  leg_tautag2->AddEntry("","","");

  gPad->SetLogx();
  histos_tautag.second->Draw();
  addHist(histos_tautag.second, leg_tautag2, 0);

  DrawAxis(histos_tautag.second, leg_tautag2, 1);
  leg_tautag2->Draw();
  C_tautag1->cd(0);

  TString jetRes = "jetERes";
  TString jetResFwd = "jetEResFwd";
  TCanvas *C_jet = new TCanvas(jetRes,jetRes, 1600, 600);
  C_jet->Divide(2);

  C_jet->cd(1);
  TMultiGraph *mg_jet = new TMultiGraph(jetRes,jetRes);
  TLegend *leg_jet = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_jet->SetHeader("#splitline{jets}{|#eta| < 1.5}");
  leg_jet->AddEntry("","","");

  addGraph(mg_jet, &gr_calojets, leg_jet, 3);
  addGraph(mg_jet, &gr_pfjets, leg_jet, 1);

  mg_jet->Draw("ALX");
  leg_jet->Draw();

  DrawAxis(mg_jet, leg_jet, 0.5);

  C_jet->cd(2);
  TMultiGraph *mg_jetFwd = new TMultiGraph(jetResFwd,jetResFwd);
  TLegend *leg_jetFwd = new TLegend(resLegXmin,resLegYmin,resLegXmax,resLegYmax);
  leg_jetFwd->SetHeader("#splitline{jets}{1.5 < |#eta| < 2.5}");
  leg_jetFwd->AddEntry("","","");

  addGraph(mg_jetFwd, &gr_calojetsFwd, leg_jetFwd, 3);
  addGraph(mg_jetFwd, &gr_pfjetsFwd, leg_jetFwd, 1);

  mg_jetFwd->Draw("ALX");
  leg_jetFwd->Draw();

  DrawAxis(mg_jetFwd, leg_jetFwd, 0.5);

  TString metRes = "MetRes";
  TCanvas *C_met = new TCanvas(metRes,metRes, 800, 600);

  TMultiGraph *mg_met = new TMultiGraph(metRes,metRes);
  TLegend *leg_met = new TLegend(topLeftLegXmin,topLeftLegYmin+0.2,topLeftLegXmax-0.2,topLeftLegYmax);
  leg_met->SetHeader("E_{T}^{miss}");
  leg_met->AddEntry("","","");


  addGraph(mg_met, &gr_met, leg_met, 0);

  mg_met->Draw("ACX");
  leg_met->Draw();

  DrawAxis(mg_met, leg_met, 300);

  mg_met->GetXaxis()->SetTitle("H_{T} [GeV]");
  mg_met->GetYaxis()->SetTitle("#sigma(ME_{x}) [GeV]");

  C_jet->Print(pdfOutput,"pdf");
  C_btag1->Print(pdfOutput,"pdf");
  C_tautag1->Print(pdfOutput,"pdf");
  C_met->Print(pdfOutput+")","pdf");

  TFile *fout = new TFile(outputFile,"recreate");

  for (int bin = 0; bin < Nbins; bin++)
  {
    plots_el.at(bin).cenResolHist->Write();
    plots_eltrack.at(bin).cenResolHist->Write();
    plots_eltower.at(bin).cenResolHist->Write();
    plots_el.at(bin).fwdResolHist->Write();
    plots_eltrack.at(bin).fwdResolHist->Write();
    plots_eltower.at(bin).fwdResolHist->Write();
  }

  histos_el.first->Write();
  histos_el.second->Write();
  histos_eltrack.first->Write();
  histos_eltrack.second->Write();
  histos_eltower.first->Write();
  histos_eltower.second->Write();
  C_el1->Write();
  C_el2->Write();

  for (int bin = 0; bin < Nbins; bin++)
  {
    plots_mu.at(bin).cenResolHist->Write();
    plots_mutrack.at(bin).cenResolHist->Write();
    plots_mu.at(bin).fwdResolHist->Write();
    plots_mutrack.at(bin).fwdResolHist->Write();
  }

  histos_mu.first->Write();
  histos_mu.second->Write();
  histos_mutrack.first->Write();
  histos_mutrack.second->Write();
  C_mu1->Write();
  C_mu->Write();

  histos_ph.first->Write();
  histos_ph.second->Write();
  C_ph1->Write();
  C_ph->Write();

  for (int bin = 0; bin < Nbins; bin++)
  {
    plots_pfjets.at(bin).cenResolHist->Write();
    plots_pfjets.at(bin).fwdResolHist->Write();
    plots_calojets.at(bin).cenResolHist->Write();
    plots_calojets.at(bin).fwdResolHist->Write();
    plots_met.at(bin).cenResolHist->Write();
  }
  histos_btag.first->Write();
  histos_btag.second->Write();
  histos_tautag.first->Write();
  histos_tautag.second->Write();
  C_btag1->Write();
  C_tautag1->Write();
  C_jet->Write();
  C_met->Write();

  fout->Write();

  cout << "** Exiting..." << endl;

  delete treeReaderElectron;
  delete treeReaderMuon;
  delete treeReaderPhoton;
  delete treeReaderJet;
  delete treeReaderBJet;
  delete treeReaderTauJet;
  delete chainElectron;
  delete chainMuon;
  delete chainPhoton;
  delete chainJet;
  delete chainBJet;
  delete chainTauJet;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char *appName = "Validation";

  if(argc != 8)
  {
    cout << " Usage: " << appName << " input_file_electron input_file_muon input_file_photon input_file_jet input_file_bjet input_file_taujet output_file" << endl;
    cout << " input_file_electron  - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_muon - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_photon - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_jet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_bjet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_taujet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " output_file - output file in ROOT format" << endl;
    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  Validation(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
}


