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
#include <typeinfo>
#include <utility>
#include <vector>

#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TStyle.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

static const int Nbins = 50;

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

unsigned int k;

struct resolPlot
{
  TH1 *resolHist;
  double ptmin;
  double ptmax;
  double etamin;
  double etamax;
  double xmin;
  double xmax;
  TString obj;

  resolPlot();
  resolPlot(double ptdown, double ptup, TString object);
  resolPlot(double etadown, double etaup, double ptdown, double ptup, TString object);
  void set(double ptdown, double ptup, TString object, double xmin = 0, double xmax = 2);
  void set(double etadown, double etaup, double ptdown, double ptup, TString object, double xmin = 0, double xmax = 2);
  void print() { std::cout << ptmin << std::endl; }
};

resolPlot::resolPlot()
{
}

resolPlot::resolPlot(double ptdown, double ptup, TString object)
{
  this->set(ptdown, ptup, object);
}

resolPlot::resolPlot(double etadown, double etaup, double ptdown, double ptup, TString object)
{
  this->set(etadown, etaup, ptdown, ptup, object);
}

void resolPlot::set(double ptdown, double ptup, TString object, double xmin, double xmax)
{
  ptmin = ptdown;
  ptmax = ptup;
  obj = object;

  resolHist = new TH1D(obj + "_delta_pt_" + Form("%4.2f", ptmin) + "_" + Form("%4.2f", ptmax), obj + "_delta_pt_" + Form("%4.2f", ptmin) + "_" + Form("%4.2f", ptmax), 1000, xmin, xmax);
}

void resolPlot::set(double etadown, double etaup, double ptdown, double ptup, TString object, double xmin, double xmax)
{
  etamin = etadown;
  etamax = etaup;
  ptmin = ptdown;
  ptmax = ptup;
  obj = object;

  resolHist = new TH1D(obj + "_delta_pt_" + Form("%4.2f", ptmin) + "_" + Form("%4.2f", ptmax) + "_" + Form("%4.2f", etamin) + "_" + Form("%4.2f", etamax), obj + "_delta_pt_" + Form("%4.2f", ptmin) + "_" + Form("%4.2f", ptmax) + "_" + Form("%4.2f", etamin) + "_" + Form("%4.2f", etamax), 1000, xmin, xmax);
}

void HistogramsCollection(std::vector<resolPlot> *histos, double ptmin, double ptmax, TString obj, double xmin = 0, double xmax = 2)
{
  double width;
  double ptdown;
  double ptup;
  resolPlot ptemp;

  for(int i = 0; i < Nbins; i++)
  {
    width = (ptmax - ptmin) / Nbins;
    ptdown = TMath::Power(10, ptmin + i * width);
    ptup = TMath::Power(10, ptmin + (i + 1) * width);
    ptemp.set(ptdown, ptup, obj, xmin, xmax);
    histos->push_back(ptemp);
  }
}

void HistogramsCollectionVsEta(std::vector<resolPlot> *histos, double etamin, double etamax, double ptmin, double ptmax, TString obj, double xmin = 0, double xmax = 2)
{
  resolPlot ptemp;
  double width;
  double etadown;
  double etaup;

  for(int i = 0; i < Nbins; i++)
  {
    width = (etamax - etamin) / Nbins;
    etadown = etamin + i * width;
    etaup = etamin + (i + 1) * width;

    ptemp.set(etadown, etaup, ptmin, ptmax, obj, xmin, xmax);
    histos->push_back(ptemp);
  }
}

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BinLogX(TH1 *h)
{
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for(int i = 0; i <= bins; i++)
  {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete new_bins;
}

//------------------------------------------------------------------------------

template <typename T>
TH1D *GetEffPt(TClonesArray *branchReco, TClonesArray *branchParticle, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D *histGenPt = new TH1D(name + " gen spectra Pt", name + " gen spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));
  TH1D *histRecoPt = new TH1D(name + " reco spectra Pt", name + " reco spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));

  histGenPt->SetDirectory(0);
  histRecoPt->SetDirectory(0);

  BinLogX(histGenPt);
  BinLogX(histRecoPt);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all generated particle in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle *)branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;

      pt = genMomentum.Pt();
      eta = TMath::Abs(genMomentum.Eta());

      if(eta > etamax || eta < etamin) continue;

      if(particle->PID == pdgID && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax)
      //if (TMath::Abs(particle->PID) == pdgID && (particle->Status>20 && particle->Status <30) && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax )
      {
        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T *)branchReco->At(j);
          recoMomentum = recoObj->P4();
          //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;

          // take the closest parton candidate
          if(TMath::Abs(pdgID) == 5)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;

            //if(jet->BTag != ) continue;
          }

          if(TMath::Abs(pdgID) == 4)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;
          }

          if(TMath::Abs(pdgID) == 1)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;
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
        histGenPt->Fill(pt);
        if(deltaR < 0.3 && bestRecoMomentum.Pt() > 0.20 * pt)
        {
          histRecoPt->Fill(pt);
        }
      }
    }
  }

  histRecoPt->Sumw2();
  histGenPt->Sumw2();

  histRecoPt->Divide(histGenPt);
  histRecoPt->Scale(100.);

  return histRecoPt;
}

template <typename T>
TH1D *GetEffEta(TClonesArray *branchReco, TClonesArray *branchParticle, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D *histGenEta = new TH1D(name + " gen spectra Eta", name + " gen spectra", Nbins, etamin, etamax);
  TH1D *histRecoEta = new TH1D(name + " reco spectra Eta", name + " reco spectra", Nbins, etamin, etamax);

  histGenEta->SetDirectory(0);
  histRecoEta->SetDirectory(0);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all generated particle in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle *)branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;

      pt = genMomentum.Pt();
      eta = genMomentum.Eta();

      if(pt > ptmax || pt < ptmin) continue;

      if(particle->PID == pdgID && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax)
      //if (TMath::Abs(particle->PID) == pdgID && (particle->Status>20 && particle->Status <30) && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax )
      {
        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T *)branchReco->At(j);
          recoMomentum = recoObj->P4();
          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
          //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;

          // take the closest parton candidate
          if(TMath::Abs(pdgID) == 5)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;
          }

          if(TMath::Abs(pdgID) == 4)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;
          }

          if(TMath::Abs(pdgID) == 1)
          {
            Jet *jet = (Jet *)recoObj;
            if(!(jet->BTag & (1 << 0))) continue;
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

        histGenEta->Fill(eta);
        if(deltaR < 0.3)
        {
          histRecoEta->Fill(eta);
        }
      }
    }
  }

  histRecoEta->Sumw2();
  histGenEta->Sumw2();

  histRecoEta->Divide(histGenEta);
  histRecoEta->Scale(100.);

  return histRecoEta;
}

//------------------------------------------------------------------------------

template <typename T>
TH1D *GetJetEffPt(TClonesArray *branchJet, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchJet->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t pt, eta;
  Long64_t entry;

  Int_t j;

  TH1D *histGenPt = new TH1D(name + " gen spectra Pt", name + " gen spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));
  TH1D *histRecoPt = new TH1D(name + " reco spectra Pt", name + " reco spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));

  histGenPt->SetDirectory(0);
  histRecoPt->SetDirectory(0);

  BinLogX(histGenPt);
  BinLogX(histRecoPt);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all reco object in event
    for(j = 0; j < branchJet->GetEntriesFast(); ++j)
    {
      recoObj = (T *)branchJet->At(j);
      recoMomentum = recoObj->P4();
      pt = recoMomentum.Pt();
      eta = TMath::Abs(recoMomentum.Eta());
      Jet *jet = (Jet *)recoObj;

      if(eta > etamax || eta < etamin) continue;
      if(pt < ptmin || pt > ptmax) continue;

      Int_t flavor = jet->Flavor;
      if(flavor == 21) flavor = 0;

      if(TMath::Abs(pdgID) == 1)
      {

        if(flavor < 4)
        {
          histGenPt->Fill(pt);
          if(jet->BTag & (1 << 0)) histRecoPt->Fill(pt);
        }
      }
      if(TMath::Abs(pdgID) == 4)
      {
        if(flavor == 4)
        {
          histGenPt->Fill(pt);
          if(jet->BTag & (1 << 0)) histRecoPt->Fill(pt);
        }
      }
      if(TMath::Abs(pdgID) == 5)
      {
        if(flavor == 5)
        {
          histGenPt->Fill(pt);
          if(jet->BTag & (1 << 0)) histRecoPt->Fill(pt);
        }
      }
    }
  }

  histRecoPt->Sumw2();
  histGenPt->Sumw2();

  histRecoPt->Divide(histGenPt);
  histRecoPt->Scale(100.);

  return histRecoPt;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TH1D *GetJetEffEta(TClonesArray *branchJet, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchJet->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t pt, eta;
  Long64_t entry;

  Int_t j;

  TH1D *histGenEta = new TH1D(name + " gen spectra Eta", name + " gen spectra", Nbins, etamin, etamax);
  TH1D *histRecoEta = new TH1D(name + " reco spectra Eta", name + " reco spectra", Nbins, etamin, etamax);

  histGenEta->SetDirectory(0);
  histRecoEta->SetDirectory(0);
  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all reco object in event
    for(j = 0; j < branchJet->GetEntriesFast(); ++j)
    {
      recoObj = (T *)branchJet->At(j);
      recoMomentum = recoObj->P4();
      pt = recoMomentum.Pt();
      eta = recoMomentum.Eta();
      Jet *jet = (Jet *)recoObj;

      if(eta > etamax || eta < etamin) continue;
      if(pt < ptmin || pt > ptmax) continue;

      Int_t flavor = jet->Flavor;
      if(flavor == 21) flavor = 0;

      if(TMath::Abs(pdgID) == 1)
      {
        if(flavor == 1 || flavor == 21)
        {
          histGenEta->Fill(eta);
          if(jet->BTag & (1 << 0)) histRecoEta->Fill(eta);
        }
      }
      if(TMath::Abs(pdgID) == 4)
      {
        if(flavor == 4)
        {
          histGenEta->Fill(eta);
          if(jet->BTag & (1 << 0)) histRecoEta->Fill(eta);
        }
      }
      if(TMath::Abs(pdgID) == 5)
      {
        if(flavor == 5)
        {
          histGenEta->Fill(eta);
          if(jet->BTag & (1 << 0)) histRecoEta->Fill(eta);
        }
      }
    }
  }

  histRecoEta->Sumw2();
  histGenEta->Sumw2();

  histRecoEta->Divide(histGenEta);
  histRecoEta->Scale(100.);

  return histRecoEta;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TH1D *GetTauEffPt(TClonesArray *branchReco, TClonesArray *branchParticle, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D *histGenPt = new TH1D(name + " gen spectra Pt", name + " gen spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));
  TH1D *histRecoPt = new TH1D(name + " reco spectra Pt", name + " reco spectra cen", Nbins, TMath::Log10(ptmin), TMath::Log10(ptmax));

  histGenPt->SetDirectory(0);
  histRecoPt->SetDirectory(0);

  BinLogX(histGenPt);
  BinLogX(histRecoPt);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all generated particle in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle *)branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;

      pt = genMomentum.Pt();
      eta = TMath::Abs(genMomentum.Eta());

      if(eta > etamax || eta < etamin) continue;

      if(particle->PID == pdgID && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax)
      {
        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T *)branchReco->At(j);
          recoMomentum = recoObj->P4();
          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
          //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;

          if(TMath::Abs(pdgID) == 1)
          {
            Jet *jet = (Jet *)recoObj;
            if(jet->TauTag != 1) continue;
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

        histGenPt->Fill(pt);
        if(deltaR < 0.3)
        {
          histRecoPt->Fill(pt);
        }
      }
    }
  }

  histRecoPt->Sumw2();
  histGenPt->Sumw2();

  histRecoPt->Divide(histGenPt);
  histRecoPt->Scale(100.);
  if(TMath::Abs(pdgID) == 15) histRecoPt->Scale(1 / 0.648);

  return histRecoPt;
}

template <typename T>
TH1D *GetTauEffEta(TClonesArray *branchReco, TClonesArray *branchParticle, TString name, int pdgID, double ptmin, double ptmax, double etamin, double etamax, ExRootTreeReader *treeReader)
{

  cout << "** Computing Efficiency of reconstructing " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D *histGenEta = new TH1D(name + " gen spectra Eta", name + " gen spectra", Nbins, etamin, etamax);
  TH1D *histRecoEta = new TH1D(name + " reco spectra Eta", name + " reco spectra", Nbins, etamin, etamax);

  histGenEta->SetDirectory(0);
  histRecoEta->SetDirectory(0);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all generated particle in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle *)branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;

      pt = genMomentum.Pt();
      eta = genMomentum.Eta();

      if(pt > ptmax || pt < ptmin) continue;

      if(particle->PID == pdgID && genMomentum.Pt() > ptmin && genMomentum.Pt() < ptmax)
      {
        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T *)branchReco->At(j);
          recoMomentum = recoObj->P4();
          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
          //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;

          if(TMath::Abs(pdgID) == 1)
          {
            Jet *jet = (Jet *)recoObj;
            if(jet->TauTag != 1) continue;
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

        histGenEta->Fill(eta);
        if(deltaR < 0.3)
        {
          histRecoEta->Fill(eta);
        }
      }
    }
  }

  histRecoEta->Sumw2();
  histGenEta->Sumw2();

  histRecoEta->Divide(histGenEta);
  histRecoEta->Scale(100.);
  if(TMath::Abs(pdgID) == 15) histRecoEta->Scale(1 / 0.648);

  return histRecoEta;
}

template <typename T>
void GetPtres(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, Double_t etaMin, Double_t etaMax, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing pt resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T *recoObj;

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
      recoObj = (T *)branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
        particle = (GenParticle *)branchParticle->At(j);
        if(particle->PID == pdgID && particle->Status == 1)
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
        pt = bestGenMomentum.Pt();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for(bin = 0; bin < Nbins; bin++)
        {
          if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta > etaMin && eta < etaMax)
          {
            histos->at(bin).resolHist->Fill(recoMomentum.Pt() / bestGenMomentum.Pt());
          }
        }
      }
    }
  }
}

template <typename T>
void GetEres(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, Double_t etaMin, Double_t etaMax, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing e resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestGenMomentum;

  Float_t deltaR;
  Float_t e, eta;
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
      recoObj = (T *)branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
        particle = (GenParticle *)branchParticle->At(j);
        if(particle->PID == pdgID && particle->Status == 1)
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
        e = bestGenMomentum.E();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for(bin = 0; bin < Nbins; bin++)
        {
          if(e > histos->at(bin).ptmin && e < histos->at(bin).ptmax && eta > etaMin && eta < etaMax)
          {
            histos->at(bin).resolHist->Fill(recoMomentum.E() / bestGenMomentum.E());
          }
        }
      }
    }
  }
}

template <typename T>
void GetPtresVsEta(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, Double_t ptMin, Double_t ptMax, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing pt resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T *recoObj;

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
      recoObj = (T *)branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
        particle = (GenParticle *)branchParticle->At(j);
        if(particle->PID == pdgID && particle->Status == 1)
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
        pt = bestGenMomentum.Pt();
        eta = bestGenMomentum.Eta();

        for(bin = 0; bin < Nbins; bin++)
        {
          if(eta > histos->at(bin).etamin && eta < histos->at(bin).etamax && pt > ptMin && pt < ptMax)
          {
            histos->at(bin).resolHist->Fill(recoMomentum.Pt() / bestGenMomentum.Pt());
          }
        }
      }
    }
  }
}

template <typename T>
void GetEresVsEta(std::vector<resolPlot> *histos, TClonesArray *branchReco, TClonesArray *branchParticle, int pdgID, Double_t eMin, Double_t eMax, ExRootTreeReader *treeReader)
{
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Computing E resolution of " << branchReco->GetName() << " induced by " << branchParticle->GetName() << " with PID " << pdgID << endl;

  GenParticle *particle;
  T *recoObj;

  TLorentzVector recoMomentum, genMomentum, bestGenMomentum;

  Float_t deltaR;
  Float_t e, eta;
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
      recoObj = (T *)branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
      {
        particle = (GenParticle *)branchParticle->At(j);
        if(particle->PID == pdgID && particle->Status == 1)
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
        e = bestGenMomentum.E();
        eta = bestGenMomentum.Eta();

        for(bin = 0; bin < Nbins; bin++)
        {
          if(eta > histos->at(bin).etamin && eta < histos->at(bin).etamax && e > eMin && e < eMax)
          {
            histos->at(bin).resolHist->Fill(recoMomentum.E() / bestGenMomentum.E());
          }
        }
      }
    }
  }
}

void GetJetsEres(std::vector<resolPlot> *histos, TClonesArray *branchJet, TClonesArray *branchGenJet, ExRootTreeReader *treeReader, Double_t etaMin, Double_t etaMax)
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

    if(entry % 10000 == 0) cout << "Event number: " << entry << endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < TMath::Min(2, branchJet->GetEntriesFast()); ++i) //branchJet->GetEntriesFast(); ++i)
    {

      jet = (Jet *)branchJet->At(i);
      jetMomentum = jet->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < TMath::Min(2, branchGenJet->GetEntriesFast()); ++j)
      {
        genjet = (Jet *)branchGenJet->At(j);

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

      if(deltaR < 0.3)
      {
        pt = genJetMomentum.E();
        eta = genJetMomentum.Eta();

        for(bin = 0; bin < Nbins; bin++)
        {
          if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta < etaMax && eta > etaMin)
          {
            histos->at(bin).resolHist->Fill(jetMomentum.E() / bestGenJetMomentum.E());
          }
        }
      }
    }
  }
}

void GetJetsEresVsEta(std::vector<resolPlot> *histos, TClonesArray *branchJet, TClonesArray *branchGenJet, ExRootTreeReader *treeReader, Double_t eMin, Double_t eMax)
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

    if(entry % 10000 == 0) cout << "Event number: " << entry << endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < TMath::Min(2, branchJet->GetEntriesFast()); ++i) //branchJet->GetEntriesFast(); ++i)
    {

      jet = (Jet *)branchJet->At(i);
      jetMomentum = jet->P4();

      deltaR = 999;

      // Loop over all hard partons in event
      for(j = 0; j < TMath::Min(2, branchGenJet->GetEntriesFast()); ++j)
      {
        genjet = (Jet *)branchGenJet->At(j);

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

      if(deltaR < 0.3)
      {

        pt = genJetMomentum.E();
        eta = genJetMomentum.Eta();

        for(bin = 0; bin < Nbins; bin++)
        {
          if(eta > histos->at(bin).etamin && eta < histos->at(bin).etamax && pt < eMax && pt > eMin)
          {
            histos->at(bin).resolHist->Fill(jetMomentum.E() / bestGenJetMomentum.E());
          }
        }
      }
    }
  }
}

std::pair<Double_t, Double_t> GausFit(TH1 *hist)
{
  TF1 *f1 = new TF1("f1", "gaus", hist->GetMean() - 2 * hist->GetRMS(), hist->GetMean() + 2 * hist->GetRMS());
  hist->Fit("f1", "RQ");

  TF1 *f2 = new TF1("f2", "gaus", f1->GetParameter(1) - 2 * f1->GetParameter(2), f1->GetParameter(1) + 2 * f1->GetParameter(2));
  hist->Fit("f2", "RQ");

  Double_t sig = f2->GetParameter(2);
  Double_t sigErr = f2->GetParError(2);

  delete f1;
  delete f2;
  return make_pair(sig, sigErr);
}

TGraphErrors EresGraph(std::vector<resolPlot> *histos, bool rms = false)
{
  Int_t bin;
  Int_t count = 0;
  TGraphErrors gr = TGraphErrors(Nbins / 2);
  double val, error;
  for(bin = 0; bin < Nbins; bin++)
  {
    std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).resolHist);
    if(rms == true)
    {
      gr.SetPoint(count, (histos->at(bin).ptmin + histos->at(bin).ptmax) / 2.0, 100 * histos->at(bin).resolHist->GetRMS());
      //gr.SetPointError(count,0, 100*sigvalues.second); // to correct
      error = 100 * histos->at(bin).resolHist->GetRMSError();
      val = 100 * histos->at(bin).resolHist->GetRMS();
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error); // to correct
    }
    else
    {

      gr.SetPoint(count, (histos->at(bin).ptmin + histos->at(bin).ptmax) / 2.0, 100 * sigvalues.first);
      error = 100 * sigvalues.second;
      val = 100 * sigvalues.first;
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error); // to correct
      //gr.SetPointError(count,0, 100*sigvalues.second);
    }
    count++;
  }

  return gr;
}

TGraphErrors MetResGraph(std::vector<resolPlot> *histos, bool rms = false)
{
  Int_t bin;
  Int_t count = 0;
  TGraphErrors gr = TGraphErrors(Nbins / 2);
  double val, error;
  for(bin = 0; bin < Nbins; bin++)
  {
    std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).resolHist);
    if(rms == true)
    {
      gr.SetPoint(count, (histos->at(bin).ptmin + histos->at(bin).ptmax) / 2.0, histos->at(bin).resolHist->GetRMS());
      error = histos->at(bin).resolHist->GetRMSError();
      val = histos->at(bin).resolHist->GetRMS();
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error); // to correct
    }
    else
    {

      gr.SetPoint(count, (histos->at(bin).ptmin + histos->at(bin).ptmax) / 2.0, sigvalues.first);
      val = sigvalues.first;
      error = sigvalues.second;
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error);
      //gr.SetPointError(count,0, 100*sigvalues.second);
    }
    count++;
  }

  return gr;
}

TGraphErrors EresGraphVsEta(std::vector<resolPlot> *histos, bool rms = false)
{
  Int_t bin;
  Int_t count = 0;
  TGraphErrors gr = TGraphErrors(Nbins / 2);
  double val, error;
  for(bin = 0; bin < Nbins; bin++)
  {

    std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).resolHist);
    if(rms == true)
    {
      gr.SetPoint(count, (histos->at(bin).etamin + histos->at(bin).etamax) / 2.0, histos->at(bin).resolHist->GetRMS());
      error = 100 * histos->at(bin).resolHist->GetRMSError();
      val = 100 * histos->at(bin).resolHist->GetRMS();
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error); // to correct
    }
    else
    {
      gr.SetPoint(count, (histos->at(bin).etamin + histos->at(bin).etamax) / 2.0, 100 * sigvalues.first);
      val = 100 * sigvalues.first;
      error = 100 * sigvalues.second;
      if(error > 0.2 * val) error = 0.2 * val;
      gr.SetPointError(count, 0, error);
      //gr.SetPointError(count,0, 100*sigvalues.second);
    }
    count++;
  }

  return gr;
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

    if(entry % 10000 == 0) cout << "Event number: " << entry << endl;

    if(branchJet->GetEntriesFast() > 1)
    {

      jet = (Jet *)branchJet->At(0);
      p1 = jet->P4();
      jet = (Jet *)branchJet->At(1);
      p2 = jet->P4();

      met = (MissingET *)branchMet->At(0);
      scalarHT = (ScalarHT *)branchScalarHT->At(0);
      ht = scalarHT->HT;

      if(p1.Pt() < 0.75 * ht / 2) continue;
      if(p2.Pt() < 0.75 * ht / 2) continue;

      for(bin = 0; bin < Nbins; bin++)
      {
        if(ht > histos->at(bin).ptmin && ht < histos->at(bin).ptmax)
        {
          histos->at(bin).resolHist->Fill(met->P4().Px());
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

void addResoGraph(TMultiGraph *mg, TGraphErrors *gr, TLegend *leg, int style, Color_t color, TString text)
{

  gr->SetLineWidth(2);
  gr->SetLineColor(color);
  gr->SetMarkerStyle(style);
  gr->SetMarkerColor(color);
  gr->SetMarkerSize(1.5);

  std::cout << "Adding " << gr->GetName() << std::endl;
  mg->Add(gr);
  leg->AddEntry(gr, text, "p");
}

void DrawAxis(TMultiGraph *mg, TLegend *leg, double xmin, double xmax, double ymin, double ymax, TString tx, TString ty, bool logx = 0, bool logy = 0)
{
  mg->SetMinimum(ymin);
  mg->SetMaximum(ymax);
  mg->GetXaxis()->SetLimits(xmin, xmax);

  mg->GetXaxis()->SetTitle(tx);
  mg->GetYaxis()->SetTitle(ty);

  mg->GetYaxis()->SetTitleSize(0.07);
  mg->GetXaxis()->SetTitleSize(0.07);
  mg->GetYaxis()->SetLabelSize(0.06);
  mg->GetXaxis()->SetLabelSize(0.06);
  mg->GetYaxis()->SetLabelOffset(0.03);
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetXaxis()->SetTitleOffset(1.4);

  mg->GetXaxis()->SetTitleFont(132);
  mg->GetYaxis()->SetTitleFont(132);
  mg->GetXaxis()->SetLabelFont(132);
  mg->GetYaxis()->SetLabelFont(132);

  mg->GetYaxis()->SetNdivisions(505);

  leg->SetTextFont(132);
  leg->SetBorderSize(0);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  gStyle->SetOptTitle(0);

  if(logx) gPad->SetLogx();
  if(logy) gPad->SetLogy();

  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->Modified();
  gPad->Update();
}

void DelphesValidation(
  const char *inputFilePion,
  const char *inputFileElectron,
  const char *inputFileMuon,
  const char *inputFilePhoton,
  const char *inputFileNeutralHadron,
  const char *inputFileJet,
  const char *inputFileBJet,
  const char *inputFileCJet,
  const char *inputFileTauJet,
  const char *outputFile,
  const char *version)
{

  TChain *chainPion = new TChain("Delphes");
  chainPion->Add(inputFilePion);
  ExRootTreeReader *treeReaderPion = new ExRootTreeReader(chainPion);

  TChain *chainElectron = new TChain("Delphes");
  chainElectron->Add(inputFileElectron);
  ExRootTreeReader *treeReaderElectron = new ExRootTreeReader(chainElectron);

  TChain *chainMuon = new TChain("Delphes");
  chainMuon->Add(inputFileMuon);
  ExRootTreeReader *treeReaderMuon = new ExRootTreeReader(chainMuon);

  TChain *chainPhoton = new TChain("Delphes");
  chainPhoton->Add(inputFilePhoton);
  ExRootTreeReader *treeReaderPhoton = new ExRootTreeReader(chainPhoton);

  TChain *chainNeutralHadron = new TChain("Delphes");
  chainNeutralHadron->Add(inputFileNeutralHadron);
  ExRootTreeReader *treeReaderNeutralHadron = new ExRootTreeReader(chainNeutralHadron);

  TChain *chainJet = new TChain("Delphes");
  chainJet->Add(inputFileJet);
  ExRootTreeReader *treeReaderJet = new ExRootTreeReader(chainJet);

  TChain *chainBJet = new TChain("Delphes");
  chainBJet->Add(inputFileBJet);
  ExRootTreeReader *treeReaderBJet = new ExRootTreeReader(chainBJet);

  TChain *chainCJet = new TChain("Delphes");
  chainCJet->Add(inputFileCJet);
  ExRootTreeReader *treeReaderCJet = new ExRootTreeReader(chainCJet);

  TChain *chainTauJet = new TChain("Delphes");
  chainTauJet->Add(inputFileTauJet);
  ExRootTreeReader *treeReaderTauJet = new ExRootTreeReader(chainTauJet);

  TClonesArray *branchParticleElectron = treeReaderElectron->UseBranch("Particle");
  TClonesArray *branchTrackElectron = treeReaderElectron->UseBranch("Track");
  TClonesArray *branchElectron = treeReaderElectron->UseBranch("Electron");
  TClonesArray *branchElectronPF = treeReaderElectron->UseBranch("ElectronPF");

  TClonesArray *branchParticleMuon = treeReaderMuon->UseBranch("Particle");
  TClonesArray *branchTrackMuon = treeReaderMuon->UseBranch("Track");
  TClonesArray *branchMuon = treeReaderMuon->UseBranch("Muon");

  TClonesArray *branchParticlePion = treeReaderPion->UseBranch("Particle");
  TClonesArray *branchTrackPion = treeReaderPion->UseBranch("Track");
  TClonesArray *branchPion = treeReaderPion->UseBranch("Pion");

  TClonesArray *branchParticlePhoton = treeReaderPhoton->UseBranch("Particle");
  TClonesArray *branchTowerPhoton = treeReaderPhoton->UseBranch("Tower");
  TClonesArray *branchPhoton = treeReaderPhoton->UseBranch("Photon");

  TClonesArray *branchParticleNeutralHadron = treeReaderNeutralHadron->UseBranch("Particle");
  TClonesArray *branchTowerNeutralHadron = treeReaderNeutralHadron->UseBranch("Tower");

  TClonesArray *branchGenJet = treeReaderJet->UseBranch("GenJet");
  TClonesArray *branchParticleJet = treeReaderJet->UseBranch("Particle");
  TClonesArray *branchPFJet = treeReaderJet->UseBranch("PFJet");
  TClonesArray *branchCaloJet = treeReaderJet->UseBranch("CaloJet");
  TClonesArray *branchJet = treeReaderJet->UseBranch("Jet");

  TClonesArray *branchParticleBJet = treeReaderBJet->UseBranch("Particle");
  TClonesArray *branchPFBJet = treeReaderBJet->UseBranch("Jet");

  TClonesArray *branchParticleCJet = treeReaderCJet->UseBranch("Particle");
  TClonesArray *branchPFCJet = treeReaderCJet->UseBranch("Jet");

  TClonesArray *branchParticleTauJet = treeReaderTauJet->UseBranch("Particle");
  TClonesArray *branchPFTauJet = treeReaderTauJet->UseBranch("Jet");

  TClonesArray *branchGenScalarHT = treeReaderJet->UseBranch("GenScalarHT");
  TClonesArray *branchMet = treeReaderJet->UseBranch("PFMissingET");
  TClonesArray *branchCaloMet = treeReaderJet->UseBranch("CaloMissingET");

  std::vector<Color_t> colors;

  colors.push_back(kBlack);
  colors.push_back(kBlue);
  colors.push_back(kRed);
  colors.push_back(kGreen + 1);
  colors.push_back(kMagenta + 1);
  colors.push_back(kOrange);

  std::vector<Int_t> markerStyles;

  markerStyles.push_back(20);
  markerStyles.push_back(21);
  markerStyles.push_back(22);
  markerStyles.push_back(23);
  markerStyles.push_back(33);
  markerStyles.push_back(34);

  TString pdfOutput(outputFile);
  pdfOutput.ReplaceAll(".root", ".pdf");

  TString figPath = inputFilePion;
  figPath.ReplaceAll(".root", "");
  figPath.ReplaceAll("root", "www/fig");
  Int_t lastSlash = figPath.Last('/');
  Int_t sizePath = figPath.Length();
  figPath.Remove(lastSlash + 1, sizePath);

  TString header = pdfOutput;
  header.ReplaceAll(".pdf", "");
  header.ReplaceAll("validation_", "");
  lastSlash = header.Last('/');
  sizePath = header.Length();
  header.Remove(0, lastSlash + 1);

  TString vrs(version);

  TPaveText *pave = new TPaveText(0.0, 0.89, 0.94, 0.94, "NDC");
  pave->SetTextAlign(30);
  pave->SetTextFont(132);
  pave->SetBorderSize(0);
  pave->SetShadowColor(0);
  pave->SetFillColor(0);
  pave->SetFillStyle(0);
  pave->AddText("Delphes " + vrs + " - " + header);

  TString s_etaMin, s_etaMax, s_eta, s_pt, s_e;

  Double_t ptMin = 1.;
  Double_t ptMax = 50000.;

  Double_t etaMin = -6.;
  Double_t etaMax = 6.;

  std::vector<Double_t> ptVals;

  ptVals.push_back(1.);
  ptVals.push_back(10.);
  ptVals.push_back(100.);
  ptVals.push_back(1000.);
  ptVals.push_back(10000.);

  std::vector<Double_t> etaVals;

  etaVals.push_back(0.);
  etaVals.push_back(1.5);
  etaVals.push_back(2.5);
  etaVals.push_back(4.0);
  etaVals.push_back(6.0);

  const int n_etabins = etaVals.size() - 1;
  const int n_ptbins = ptVals.size();

  //////////////////////////
  // Tracking performance //
  //////////////////////////

  // --------- Pion Tracks  --------- //

  TMultiGraph *mg_trkpi_res_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkpi_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkpi_res_eta = new TMultiGraph("", "");
  TMultiGraph *mg_trkpi_eff_eta = new TMultiGraph("", "");

  TLegend *leg_trkpi_res_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_trkpi_eff_pt = (TLegend *)leg_trkpi_res_pt->Clone();
  TLegend *leg_trkpi_res_eta = (TLegend *)leg_trkpi_res_pt->Clone();
  TLegend *leg_trkpi_eff_eta = (TLegend *)leg_trkpi_res_eta->Clone();

  TGraphErrors *gr_trkpi_res_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkpi_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkpi_res_eta = new TGraphErrors[n_ptbins];
  TGraphErrors *gr_trkpi_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_trkpi_eff_pt, *h_trkpi_eff_eta;

  std::vector<resolPlot> *plots_trkpi_res_pt = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_trkpi_res_eta = new std::vector<resolPlot>[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    HistogramsCollection(&plots_trkpi_res_pt[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "trkpi");
    GetPtres<Track>(&plots_trkpi_res_pt[k], branchTrackPion, branchParticlePion, 211, etaVals.at(k), etaVals.at(k + 1), treeReaderPion);
    gr_trkpi_res_pt[k] = EresGraph(&plots_trkpi_res_pt[k]);

    h_trkpi_eff_pt = GetEffPt<Track>(branchTrackPion, branchParticlePion, "Pion", 211, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderPion);
    gr_trkpi_eff_pt[k] = TGraphErrors(h_trkpi_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#pi^{ #pm} , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_trkpi_res_pt[k].SetName("trkRes_" + s_etaMin + "_" + s_etaMax);
    gr_trkpi_eff_pt[k].SetName("trkEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_trkpi_res_pt, &gr_trkpi_res_pt[k], leg_trkpi_res_pt, markerStyles.at(k), colors.at(k), s_eta);
    addResoGraph(mg_trkpi_eff_pt, &gr_trkpi_eff_pt[k], leg_trkpi_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    HistogramsCollectionVsEta(&plots_trkpi_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "trkpi", 0.0, 2.0);
    GetPtresVsEta<Track>(&plots_trkpi_res_eta[k], branchTrackPion, branchParticlePion, 211, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderPion);
    gr_trkpi_res_eta[k] = EresGraphVsEta(&plots_trkpi_res_eta[k]);

    h_trkpi_eff_eta = GetEffEta<Track>(branchTrackPion, branchParticlePion, "Pion", 211, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderPion);
    gr_trkpi_eff_eta[k] = TGraphErrors(h_trkpi_eff_eta);

    s_pt = Form("#pi^{ #pm} , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("#pi^{ #pm} , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_trkpi_res_eta, &gr_trkpi_res_eta[k], leg_trkpi_res_eta, markerStyles.at(k), colors.at(k), s_pt);
    addResoGraph(mg_trkpi_eff_eta, &gr_trkpi_eff_eta[k], leg_trkpi_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_trkpi_res_pt = new TCanvas("", "", 800, 600);

  mg_trkpi_res_pt->Draw("APE");
  DrawAxis(mg_trkpi_res_pt, leg_trkpi_res_pt, ptMin, ptMax, 0.01, 100, "p_{T} [GeV]", "(track resolution in p_{T})/p_{T} (%)", true, true);
  leg_trkpi_res_pt->Draw();
  pave->Draw();

  c_trkpi_res_pt->Print(pdfOutput + "(", "pdf");
  c_trkpi_res_pt->Print(figPath + "img_trkpi_res_pt.pdf", "pdf");
  c_trkpi_res_pt->Print(figPath + "img_trkpi_res_pt.png", "png");

  TCanvas *c_trkpi_res_eta = new TCanvas("", "", 800, 600);

  mg_trkpi_res_eta->Draw("APE");
  DrawAxis(mg_trkpi_res_eta, leg_trkpi_res_eta, etaMin, etaMax, 0.01, 100, " #eta ", "(track resolution in p_{T})/p_{T} (%)", false, true);
  leg_trkpi_res_eta->Draw();
  pave->Draw();

  c_trkpi_res_eta->Print(pdfOutput, "pdf");
  c_trkpi_res_eta->Print(figPath + "img_trkpi_res_eta.pdf", "pdf");
  c_trkpi_res_eta->Print(figPath + "img_trkpi_res_eta.png", "png");

  TCanvas *c_trkpi_eff_pt = new TCanvas("", "", 800, 600);

  mg_trkpi_eff_pt->Draw("APE");
  DrawAxis(mg_trkpi_eff_pt, leg_trkpi_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "tracking efficiency (%)", true, false);
  leg_trkpi_eff_pt->Draw();
  pave->Draw();

  c_trkpi_eff_pt->Print(pdfOutput, "pdf");
  c_trkpi_eff_pt->Print(figPath + "img_trkpi_eff_pt.pdf", "pdf");
  c_trkpi_eff_pt->Print(figPath + "img_trkpi_eff_pt.png", "png");

  TCanvas *c_trkpi_eff_eta = new TCanvas("", "", 800, 600);

  mg_trkpi_eff_eta->Draw("APE");
  DrawAxis(mg_trkpi_eff_eta, leg_trkpi_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "tracking efficiency (%)", false, false);
  leg_trkpi_eff_eta->Draw();
  pave->Draw();

  c_trkpi_eff_eta->Print(pdfOutput, "pdf");
  c_trkpi_eff_eta->Print(figPath + "img_trkpi_eff_eta.pdf", "pdf");
  c_trkpi_eff_eta->Print(figPath + "img_trkpi_eff_eta.png", "png");

  // --------- Electron Tracks  --------- //

  TMultiGraph *mg_trkele_res_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkele_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkele_res_eta = new TMultiGraph("", "");
  TMultiGraph *mg_trkele_eff_eta = new TMultiGraph("", "");

  TLegend *leg_trkele_res_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_trkele_eff_pt = (TLegend *)leg_trkele_res_pt->Clone();
  TLegend *leg_trkele_res_eta = (TLegend *)leg_trkele_res_pt->Clone();
  TLegend *leg_trkele_eff_eta = (TLegend *)leg_trkele_res_eta->Clone();

  TGraphErrors *gr_trkele_res_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkele_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkele_res_eta = new TGraphErrors[n_ptbins];
  TGraphErrors *gr_trkele_eff_eta = new TGraphErrors[n_ptbins];

  TH1D *h_trkele_eff_pt, *h_trkele_eff_eta;

  std::vector<resolPlot> *plots_trkele_res_pt = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_trkele_res_eta = new std::vector<resolPlot>[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    HistogramsCollection(&plots_trkele_res_pt[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "trkele");
    GetPtres<Track>(&plots_trkele_res_pt[k], branchTrackElectron, branchParticleElectron, 11, etaVals.at(k), etaVals.at(k + 1), treeReaderElectron);
    gr_trkele_res_pt[k] = EresGraph(&plots_trkele_res_pt[k]);

    h_trkele_eff_pt = GetEffPt<Track>(branchTrackElectron, branchParticleElectron, "Electron", 11, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderElectron);
    gr_trkele_eff_pt[k] = TGraphErrors(h_trkele_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "e^{ #pm} , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_trkele_res_pt[k].SetName("trkRes_" + s_etaMin + "_" + s_etaMax);
    gr_trkele_eff_pt[k].SetName("trkEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_trkele_res_pt, &gr_trkele_res_pt[k], leg_trkele_res_pt, markerStyles.at(k), colors.at(k), s_eta);
    addResoGraph(mg_trkele_eff_pt, &gr_trkele_eff_pt[k], leg_trkele_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    HistogramsCollectionVsEta(&plots_trkele_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "trkele", 0.0, 2.0);
    GetPtresVsEta<Track>(&plots_trkele_res_eta[k], branchTrackElectron, branchParticleElectron, 11, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderElectron);
    gr_trkele_res_eta[k] = EresGraphVsEta(&plots_trkele_res_eta[k]);

    h_trkele_eff_eta = GetEffEta<Track>(branchTrackElectron, branchParticleElectron, "Electron", 11, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderElectron);
    gr_trkele_eff_eta[k] = TGraphErrors(h_trkele_eff_eta);

    s_pt = Form("e^{ #pm} , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("e^{ #pm} , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_trkele_res_eta, &gr_trkele_res_eta[k], leg_trkele_res_eta, markerStyles.at(k), colors.at(k), s_pt);
    addResoGraph(mg_trkele_eff_eta, &gr_trkele_eff_eta[k], leg_trkele_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_trkele_res_pt = new TCanvas("", "", 800, 600);

  mg_trkele_res_pt->Draw("APE");
  DrawAxis(mg_trkele_res_pt, leg_trkele_res_pt, ptMin, ptMax, 0.01, 100, "p_{T} [GeV]", "(track resolution in p_{T})/p_{T} (%)", true, true);
  leg_trkele_res_pt->Draw();
  pave->Draw();

  c_trkele_res_pt->Print(pdfOutput, "pdf");
  c_trkele_res_pt->Print(figPath + "img_trkele_res_pt.pdf", "pdf");
  c_trkele_res_pt->Print(figPath + "img_trkele_res_pt.png", "png");

  TCanvas *c_trkele_res_eta = new TCanvas("", "", 800, 600);

  mg_trkele_res_eta->Draw("APE");
  DrawAxis(mg_trkele_res_eta, leg_trkele_res_eta, etaMin, etaMax, 0.01, 100, " #eta ", "(track resolution in p_{T})/p_{T} (%)", false, true);
  leg_trkele_res_eta->Draw();
  pave->Draw();

  c_trkele_res_eta->Print(pdfOutput, "pdf");
  c_trkele_res_eta->Print(figPath + "img_trkele_res_eta.pdf", "pdf");
  c_trkele_res_eta->Print(figPath + "img_trkele_res_eta.png", "png");

  TCanvas *c_trkele_eff_pt = new TCanvas("", "", 800, 600);

  mg_trkele_eff_pt->Draw("APE");
  DrawAxis(mg_trkele_eff_pt, leg_trkele_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "tracking efficiency (%)", true, false);
  leg_trkele_eff_pt->Draw();
  pave->Draw();

  c_trkele_eff_pt->Print(pdfOutput, "pdf");
  c_trkele_eff_pt->Print(figPath + "img_trkele_eff_pt.pdf", "pdf");
  c_trkele_eff_pt->Print(figPath + "img_trkele_eff_pt.png", "png");

  TCanvas *c_trkele_eff_eta = new TCanvas("", "", 800, 600);

  mg_trkele_eff_eta->Draw("APE");
  DrawAxis(mg_trkele_eff_eta, leg_trkele_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "tracking efficiency (%)", false, false);
  leg_trkele_eff_eta->Draw();
  pave->Draw();

  c_trkele_eff_eta->Print(pdfOutput, "pdf");
  c_trkele_eff_eta->Print(figPath + "img_trkele_eff_eta.pdf", "pdf");
  c_trkele_eff_eta->Print(figPath + "img_trkele_eff_eta.png", "png");

  // --------- Muon Tracks  --------- //

  TMultiGraph *mg_trkmu_res_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkmu_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_trkmu_res_eta = new TMultiGraph("", "");
  TMultiGraph *mg_trkmu_eff_eta = new TMultiGraph("", "");

  TLegend *leg_trkmu_res_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_trkmu_eff_pt = (TLegend *)leg_trkmu_res_pt->Clone();

  TLegend *leg_trkmu_res_eta = (TLegend *)leg_trkmu_res_pt->Clone();
  TLegend *leg_trkmu_eff_eta = (TLegend *)leg_trkmu_res_eta->Clone();

  TGraphErrors *gr_trkmu_res_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkmu_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkmu_res_eta = new TGraphErrors[n_ptbins];
  TGraphErrors *gr_trkmu_eff_eta = new TGraphErrors[n_ptbins];

  TH1D *h_trkmu_eff_pt, *h_trkmu_eff_eta;

  std::vector<resolPlot> *plots_trkmu_res_pt = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_trkmu_res_eta = new std::vector<resolPlot>[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    HistogramsCollection(&plots_trkmu_res_pt[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "trkmu");
    GetPtres<Track>(&plots_trkmu_res_pt[k], branchTrackMuon, branchParticleMuon, 13, etaVals.at(k), etaVals.at(k + 1), treeReaderMuon);
    gr_trkmu_res_pt[k] = EresGraph(&plots_trkmu_res_pt[k]);

    h_trkmu_eff_pt = GetEffPt<Track>(branchTrackMuon, branchParticleMuon, "Muon", 13, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderMuon);
    gr_trkmu_eff_pt[k] = TGraphErrors(h_trkmu_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#mu^{ #pm} , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_trkmu_res_pt[k].SetName("trkRes_" + s_etaMin + "_" + s_etaMax);
    gr_trkmu_eff_pt[k].SetName("trkEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_trkmu_res_pt, &gr_trkmu_res_pt[k], leg_trkmu_res_pt, markerStyles.at(k), colors.at(k), s_eta);
    addResoGraph(mg_trkmu_eff_pt, &gr_trkmu_eff_pt[k], leg_trkmu_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    HistogramsCollectionVsEta(&plots_trkmu_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "trkmu", 0.0, 2.0);
    GetPtresVsEta<Track>(&plots_trkmu_res_eta[k], branchTrackMuon, branchParticleMuon, 13, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderMuon);
    gr_trkmu_res_eta[k] = EresGraphVsEta(&plots_trkmu_res_eta[k]);

    h_trkmu_eff_eta = GetEffEta<Track>(branchTrackMuon, branchParticleMuon, "Muon", 13, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderMuon);
    gr_trkmu_eff_eta[k] = TGraphErrors(h_trkmu_eff_eta);

    s_pt = Form("#mu^{ #pm} , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("#mu^{ #pm} , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_trkmu_res_eta, &gr_trkmu_res_eta[k], leg_trkmu_res_eta, markerStyles.at(k), colors.at(k), s_pt);
    addResoGraph(mg_trkmu_eff_eta, &gr_trkmu_eff_eta[k], leg_trkmu_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_trkmu_res_pt = new TCanvas("", "", 800, 600);

  mg_trkmu_res_pt->Draw("APE");
  DrawAxis(mg_trkmu_res_pt, leg_trkmu_res_pt, ptMin, ptMax, 0.01, 100, "p_{T} [GeV]", "(track resolution in p_{T})/p_{T} (%)", true, true);
  leg_trkmu_res_pt->Draw();
  pave->Draw();

  c_trkmu_res_pt->Print(pdfOutput, "pdf");
  c_trkmu_res_pt->Print(figPath + "img_trkmu_res_pt.pdf", "pdf");
  c_trkmu_res_pt->Print(figPath + "img_trkmu_res_pt.png", "png");

  TCanvas *c_trkmu_res_eta = new TCanvas("", "", 800, 600);

  mg_trkmu_res_eta->Draw("APE");
  DrawAxis(mg_trkmu_res_eta, leg_trkmu_res_eta, etaMin, etaMax, 0.01, 100, " #eta ", "(track resolution in p_{T})/p_{T} (%)", false, true);
  leg_trkmu_res_eta->Draw();
  pave->Draw();

  c_trkmu_res_eta->Print(pdfOutput, "pdf");
  c_trkmu_res_eta->Print(figPath + "img_trkmu_res_eta.pdf", "pdf");
  c_trkmu_res_eta->Print(figPath + "img_trkmu_res_eta.png", "png");

  TCanvas *c_trkmu_eff_pt = new TCanvas("", "", 800, 600);

  mg_trkmu_eff_pt->Draw("APE");
  DrawAxis(mg_trkmu_eff_pt, leg_trkmu_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "tracking efficiency (%)", true, false);
  leg_trkmu_eff_pt->Draw();
  pave->Draw();

  c_trkmu_eff_pt->Print(pdfOutput, "pdf");
  c_trkmu_eff_pt->Print(figPath + "img_trkmu_eff_pt.pdf", "pdf");
  c_trkmu_eff_pt->Print(figPath + "img_trkmu_eff_pt.png", "png");

  TCanvas *c_trkmu_eff_eta = new TCanvas("", "", 800, 600);

  mg_trkmu_eff_eta->Draw("APE");
  DrawAxis(mg_trkmu_eff_eta, leg_trkmu_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "tracking efficiency (%)", false, false);
  leg_trkmu_eff_eta->Draw();
  pave->Draw();

  c_trkmu_eff_eta->Print(pdfOutput, "pdf");
  c_trkmu_eff_eta->Print(figPath + "img_trkmu_eff_eta.pdf", "pdf");
  c_trkmu_eff_eta->Print(figPath + "img_trkmu_eff_eta.png", "png");

  //////////////////////
  // Ecal performance //
  //////////////////////

  TMultiGraph *mg_ecal_res_e = new TMultiGraph("", "");
  TMultiGraph *mg_ecal_res_eta = new TMultiGraph("", "");

  TLegend *leg_ecal_res_e = new TLegend(0.55, 0.64, 0.90, 0.90);
  TLegend *leg_ecal_res_eta = new TLegend(0.60, 0.59, 0.95, 0.90);

  TGraphErrors *gr_ecal_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_ecal_res_eta = new TGraphErrors[n_ptbins];

  std::vector<resolPlot> *plots_ecal_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_ecal_res_eta = new std::vector<resolPlot>[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    HistogramsCollection(&plots_ecal_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "ecal");
    GetEres<Tower>(&plots_ecal_res_e[k], branchTowerPhoton, branchParticlePhoton, 22, etaVals.at(k), etaVals.at(k + 1), treeReaderPhoton);
    gr_ecal_res_e[k] = EresGraph(&plots_ecal_res_e[k]);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#gamma , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_ecal_res_e[k].SetName("trkRes_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_ecal_res_e, &gr_ecal_res_e[k], leg_ecal_res_e, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    HistogramsCollectionVsEta(&plots_ecal_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "ecal", 0.0, 2.0);
    GetEresVsEta<Tower>(&plots_ecal_res_eta[k], branchTowerPhoton, branchParticlePhoton, 22, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderPhoton);
    gr_ecal_res_eta[k] = EresGraphVsEta(&plots_ecal_res_eta[k]);

    s_e = Form("#gamma , E = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_e = Form("#gamma , E = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_ecal_res_eta, &gr_ecal_res_eta[k], leg_ecal_res_eta, markerStyles.at(k), colors.at(k), s_e);
  }

  TCanvas *c_ecal_res_e = new TCanvas("", "", 800, 600);

  mg_ecal_res_e->Draw("APE");
  // DrawAxis(mg_ecal_res_e, leg_ecal_res_e, ptMin, ptMax, 0.5, 100, "E [GeV]", "(ECAL resolution in E)/E (%)", true, true);
  DrawAxis(mg_ecal_res_e, leg_ecal_res_e, ptMin, ptMax, 0.0, 20, "E [GeV]", "(ECAL resolution in E)/E (%)", true, false);
  leg_ecal_res_e->Draw();
  pave->Draw();

  c_ecal_res_e->Print(pdfOutput, "pdf");
  c_ecal_res_e->Print(figPath + "img_ecal_res_e.pdf", "pdf");
  c_ecal_res_e->Print(figPath + "img_ecal_res_e.png", "png");

  TCanvas *c_ecal_res_eta = new TCanvas("", "", 800, 600);

  mg_ecal_res_eta->Draw("APE");
  //DrawAxis(mg_ecal_res_eta, leg_ecal_res_eta, etaMin, etaMax, 0.5, 100, " #eta ", "(ECAL resolution in E)/E (%)", false, true);
  DrawAxis(mg_ecal_res_eta, leg_ecal_res_eta, etaMin, etaMax, 0.0, 20, " #eta ", "(ECAL resolution in E)/E (%)", false, false);
  leg_ecal_res_eta->Draw();
  pave->Draw();

  c_ecal_res_eta->Print(pdfOutput, "pdf");
  c_ecal_res_eta->Print(figPath + "img_ecal_res_eta.pdf", "pdf");
  c_ecal_res_eta->Print(figPath + "img_ecal_res_eta.png", "png");

  //////////////////////
  // Hcal performance //
  //////////////////////

  TMultiGraph *mg_hcal_res_e = new TMultiGraph("", "");
  TMultiGraph *mg_hcal_res_eta = new TMultiGraph("", "");

  TLegend *leg_hcal_res_e = new TLegend(0.55, 0.64, 0.90, 0.90);
  TLegend *leg_hcal_res_eta = new TLegend(0.60, 0.59, 0.95, 0.90);

  TGraphErrors *gr_hcal_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_hcal_res_eta = new TGraphErrors[n_ptbins];

  std::vector<resolPlot> *plots_hcal_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_hcal_res_eta = new std::vector<resolPlot>[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    HistogramsCollection(&plots_hcal_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "hcal");
    GetEres<Tower>(&plots_hcal_res_e[k], branchTowerNeutralHadron, branchParticleNeutralHadron, 2112, etaVals.at(k), etaVals.at(k + 1), treeReaderNeutralHadron);

    gr_hcal_res_e[k] = EresGraph(&plots_hcal_res_e[k]);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "n , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_hcal_res_e[k].SetName("trkRes_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_hcal_res_e, &gr_hcal_res_e[k], leg_hcal_res_e, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    HistogramsCollectionVsEta(&plots_hcal_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "hcal", 0.0, 2.0);
    GetEresVsEta<Tower>(&plots_hcal_res_eta[k], branchTowerNeutralHadron, branchParticleNeutralHadron, 2112, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderNeutralHadron);
    gr_hcal_res_eta[k] = EresGraphVsEta(&plots_hcal_res_eta[k]);

    s_e = Form("n , E = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_e = Form("n , E = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_hcal_res_eta, &gr_hcal_res_eta[k], leg_hcal_res_eta, markerStyles.at(k), colors.at(k), s_e);
  }

  TCanvas *c_hcal_res_e = new TCanvas("", "", 800, 600);

  mg_hcal_res_e->Draw("APE");
  //DrawAxis(mg_hcal_res_e, leg_hcal_res_e, ptMin, ptMax, 1, 100, "E [GeV]", "(HCAL resolution in E)/E (%)", true, true);
  DrawAxis(mg_hcal_res_e, leg_hcal_res_e, ptMin, ptMax, 0.0, 50, "E [GeV]", "(HCAL resolution in E)/E (%)", true, false);
  leg_hcal_res_e->Draw();
  pave->Draw();

  c_hcal_res_e->Print(pdfOutput, "pdf");
  c_hcal_res_e->Print(figPath + "img_hcal_res_e.pdf", "pdf");
  c_hcal_res_e->Print(figPath + "img_hcal_res_e.png", "png");

  TCanvas *c_hcal_res_eta = new TCanvas("", "", 800, 600);

  mg_hcal_res_eta->Draw("APE");
  //DrawAxis(mg_hcal_res_eta, leg_hcal_res_eta, etaMin, etaMax, 1, 100, " #eta ", "(HCAL resolution in E)/E (%)", false, true);
  DrawAxis(mg_hcal_res_eta, leg_hcal_res_eta, etaMin, etaMax, 0.0, 50, " #eta ", "(HCAL resolution in E)/E (%)", false, false);
  leg_hcal_res_eta->Draw();
  pave->Draw();

  c_hcal_res_eta->Print(pdfOutput, "pdf");
  c_hcal_res_eta->Print(figPath + "img_hcal_res_eta.pdf", "pdf");
  c_hcal_res_eta->Print(figPath + "img_hcal_res_eta.png", "png");

  ////////////////////
  // PF - electrons //
  ////////////////////

  TMultiGraph *mg_pfele_res_e[n_etabins];
  TMultiGraph *mg_pfele_res_eta[n_ptbins];

  TLegend *leg_pfele_res_e[n_etabins];
  TLegend *leg_pfele_res_eta[n_ptbins];

  TGraphErrors *gr_pfele_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_pfele_res_eta = new TGraphErrors[n_ptbins];
  TGraphErrors *gr_trkele_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkele_res_eeta = new TGraphErrors[n_ptbins];

  std::vector<resolPlot> *plots_pfele_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_pfele_res_eta = new std::vector<resolPlot>[n_ptbins];
  std::vector<resolPlot> *plots_trkele_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_trkele_res_eeta = new std::vector<resolPlot>[n_ptbins];

  TCanvas *c_pfele_res_e[n_etabins];
  TCanvas *c_pfele_res_eta[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    mg_pfele_res_e[k] = new TMultiGraph("", "");
    leg_pfele_res_e[k] = new TLegend(0.40, 0.60, 0.75, 0.90);

    HistogramsCollection(&plots_pfele_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "pfele");
    GetEres<Electron>(&plots_pfele_res_e[k], branchElectronPF, branchParticleElectron, 11, etaVals.at(k), etaVals.at(k + 1), treeReaderElectron);
    gr_pfele_res_e[k] = EresGraph(&plots_pfele_res_e[k]);

    HistogramsCollection(&plots_trkele_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "trkele");
    GetEres<Track>(&plots_trkele_res_e[k], branchTrackElectron, branchParticleElectron, 11, etaVals.at(k), etaVals.at(k + 1), treeReaderElectron);
    gr_trkele_res_e[k] = EresGraph(&plots_trkele_res_e[k]);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));
    s_eta = "e^{ #pm}, " + s_etaMin + " < | #eta | < " + s_etaMax;

    leg_pfele_res_e[k]->SetTextFont(132);
    leg_pfele_res_e[k]->SetHeader(s_eta);

    addResoGraph(mg_pfele_res_e[k], &gr_ecal_res_e[k], leg_pfele_res_e[k], markerStyles.at(0), colors.at(0), "ECAL");
    addResoGraph(mg_pfele_res_e[k], &gr_trkele_res_e[k], leg_pfele_res_e[k], markerStyles.at(1), colors.at(1), "Track");
    addResoGraph(mg_pfele_res_e[k], &gr_pfele_res_e[k], leg_pfele_res_e[k], markerStyles.at(2), colors.at(2), "Particle-flow");

    c_pfele_res_e[k] = new TCanvas("", "", 800, 600);

    mg_pfele_res_e[k]->Draw("APE");
    //DrawAxis(mg_pfele_res_e[k], leg_pfele_res_e[k], ptMin, ptMax, 0.1, 100, "E [GeV]", "(resolution in E)/E (%)", true, true);
    DrawAxis(mg_pfele_res_e[k], leg_pfele_res_e[k], ptMin, ptMax, 0.0, 20, "E [GeV]", "(resolution in E)/E (%)", true, false);
    leg_pfele_res_e[k]->Draw();
    pave->Draw();

    TString s_etarange = "eta_" + s_etaMin + "_" + s_etaMax + "_";

    c_pfele_res_e[k]->Print(pdfOutput, "pdf");
    c_pfele_res_e[k]->Print(figPath + "img_pfele_res_" + s_etarange + "e.pdf", "pdf");
    c_pfele_res_e[k]->Print(figPath + "img_pfele_res_" + s_etarange + "e.png", "png");
  }

  // loop over eta bins
  for(k = 0; k < ptVals.size(); k++)
  {

    mg_pfele_res_eta[k] = new TMultiGraph("", "");
    leg_pfele_res_eta[k] = new TLegend(0.40, 0.60, 0.75, 0.90);

    HistogramsCollectionVsEta(&plots_pfele_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "pfele", 0.0, 2.0);
    GetEresVsEta<Electron>(&plots_pfele_res_eta[k], branchElectronPF, branchParticleElectron, 11, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderElectron);
    gr_pfele_res_eta[k] = EresGraphVsEta(&plots_pfele_res_eta[k]);

    HistogramsCollectionVsEta(&plots_trkele_res_eeta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "trkele", 0.0, 2.0);
    GetEresVsEta<Track>(&plots_trkele_res_eeta[k], branchTrackElectron, branchParticleElectron, 11, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderElectron);
    gr_trkele_res_eeta[k] = EresGraphVsEta(&plots_trkele_res_eeta[k]);

    s_e = Form("e^{ #pm}, E = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_e = Form("e^{ #pm}, E = %.0f TeV", ptVals.at(k) / 1000.);

    leg_pfele_res_eta[k]->SetTextFont(132);
    leg_pfele_res_eta[k]->SetHeader(s_e);

    addResoGraph(mg_pfele_res_eta[k], &gr_ecal_res_eta[k], leg_pfele_res_eta[k], markerStyles.at(0), colors.at(0), "ECAL");
    addResoGraph(mg_pfele_res_eta[k], &gr_trkele_res_eeta[k], leg_pfele_res_eta[k], markerStyles.at(1), colors.at(1), "Track");
    addResoGraph(mg_pfele_res_eta[k], &gr_pfele_res_eta[k], leg_pfele_res_eta[k], markerStyles.at(2), colors.at(2), "Particle-flow");

    c_pfele_res_eta[k] = new TCanvas("", "", 800, 600);

    mg_pfele_res_eta[k]->Draw("APE");
    //DrawAxis(mg_pfele_res_eta[k], leg_pfele_res_eta[k], etaMin, etaMax, 0.1, 1000, "#eta", "(resolution in E)/E (%)", false, true);
    DrawAxis(mg_pfele_res_eta[k], leg_pfele_res_eta[k], etaMin, etaMax, 0.0, 50, "#eta", "(resolution in E)/E (%)", false, false);
    leg_pfele_res_eta[k]->Draw();
    pave->Draw();

    TString s_ptrange = Form("pt_%.0f_", ptVals.at(k));

    c_pfele_res_eta[k]->Print(pdfOutput, "pdf");
    c_pfele_res_eta[k]->Print(figPath + "img_pfele_res_" + s_ptrange + "eta.pdf", "pdf");
    c_pfele_res_eta[k]->Print(figPath + "img_pfele_res_" + s_ptrange + "eta.png", "png");
  }

  /////////////////
  // PF - Pions  //
  /////////////////

  TMultiGraph *mg_pfpi_res_e[n_etabins];
  TMultiGraph *mg_pfpi_res_eta[n_ptbins];

  TLegend *leg_pfpi_res_e[n_etabins];
  TLegend *leg_pfpi_res_eta[n_ptbins];

  TGraphErrors *gr_pfpi_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_pfpi_res_eta = new TGraphErrors[n_ptbins];

  TGraphErrors *gr_trkpi_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_trkpi_res_eeta = new TGraphErrors[n_ptbins];

  std::vector<resolPlot> *plots_pfpi_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_pfpi_res_eta = new std::vector<resolPlot>[n_ptbins];
  std::vector<resolPlot> *plots_trkpi_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_trkpi_res_eeta = new std::vector<resolPlot>[n_ptbins];

  TCanvas *c_pfpi_res_e[n_etabins];
  TCanvas *c_pfpi_res_eta[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {
    mg_pfpi_res_e[k] = new TMultiGraph("", "");
    leg_pfpi_res_e[k] = new TLegend(0.40, 0.60, 0.75, 0.90);

    HistogramsCollection(&plots_pfpi_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "pfpi");
    GetEres<Track>(&plots_pfpi_res_e[k], branchPion, branchParticlePion, 211, etaVals.at(k), etaVals.at(k + 1), treeReaderPion);
    gr_pfpi_res_e[k] = EresGraph(&plots_pfpi_res_e[k]);

    HistogramsCollection(&plots_trkpi_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "trkpi");
    GetEres<Track>(&plots_trkpi_res_e[k], branchTrackPion, branchParticlePion, 211, etaVals.at(k), etaVals.at(k + 1), treeReaderPion);
    gr_trkpi_res_e[k] = EresGraph(&plots_trkpi_res_e[k]);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));
    s_eta = "#pi^{ #pm}, " + s_etaMin + " < | #eta | < " + s_etaMax;

    leg_pfpi_res_e[k]->SetTextFont(132);
    leg_pfpi_res_e[k]->SetHeader(s_eta);

    addResoGraph(mg_pfpi_res_e[k], &gr_hcal_res_e[k], leg_pfpi_res_e[k], markerStyles.at(0), colors.at(0), "HCAL");
    addResoGraph(mg_pfpi_res_e[k], &gr_trkpi_res_e[k], leg_pfpi_res_e[k], markerStyles.at(1), colors.at(1), "Track");
    addResoGraph(mg_pfpi_res_e[k], &gr_pfpi_res_e[k], leg_pfpi_res_e[k], markerStyles.at(2), colors.at(2), "Particle-flow");

    c_pfpi_res_e[k] = new TCanvas("", "", 800, 600);

    mg_pfpi_res_e[k]->Draw("APE");
    //DrawAxis(mg_pfpi_res_e[k], leg_pfpi_res_e[k], ptMin, ptMax, 0.1, 100, "E [GeV]", "(resolution in E)/E (%)", true, true);
    DrawAxis(mg_pfpi_res_e[k], leg_pfpi_res_e[k], ptMin, ptMax, 0.1, 50, "E [GeV]", "(resolution in E)/E (%)", true, false);
    leg_pfpi_res_e[k]->Draw();
    pave->Draw();

    TString s_etarange = "eta_" + s_etaMin + "_" + s_etaMax + "_";

    c_pfpi_res_e[k]->Print(pdfOutput, "pdf");
    c_pfpi_res_e[k]->Print(figPath + "img_pfpi_res_" + s_etarange + "e.pdf", "pdf");
    c_pfpi_res_e[k]->Print(figPath + "img_pfpi_res_" + s_etarange + "e.png", "png");
  }

  // loop over eta bins
  for(k = 0; k < ptVals.size(); k++)
  {

    mg_pfpi_res_eta[k] = new TMultiGraph("", "");
    leg_pfpi_res_eta[k] = new TLegend(0.40, 0.60, 0.75, 0.90);

    HistogramsCollectionVsEta(&plots_pfpi_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "pfpi", 0.0, 2.0);
    GetEresVsEta<Track>(&plots_pfpi_res_eta[k], branchPion, branchParticlePion, 211, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderPion);
    gr_pfpi_res_eta[k] = EresGraphVsEta(&plots_pfpi_res_eta[k]);

    HistogramsCollectionVsEta(&plots_trkpi_res_eeta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "trkpi", 0.0, 2.0);
    GetEresVsEta<Track>(&plots_trkpi_res_eeta[k], branchPion, branchParticlePion, 211, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), treeReaderPion);
    gr_trkpi_res_eeta[k] = EresGraphVsEta(&plots_trkpi_res_eeta[k]);

    s_e = Form("#pi^{ #pm}, E = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_e = Form("#pi^{ #pm}, E = %.0f TeV", ptVals.at(k) / 1000.);

    leg_pfpi_res_eta[k]->SetTextFont(132);
    leg_pfpi_res_eta[k]->SetHeader(s_e);

    addResoGraph(mg_pfpi_res_eta[k], &gr_hcal_res_eta[k], leg_pfpi_res_eta[k], markerStyles.at(0), colors.at(0), "HCAL");
    addResoGraph(mg_pfpi_res_eta[k], &gr_trkpi_res_eeta[k], leg_pfpi_res_eta[k], markerStyles.at(1), colors.at(1), "Track");
    addResoGraph(mg_pfpi_res_eta[k], &gr_pfpi_res_eta[k], leg_pfpi_res_eta[k], markerStyles.at(2), colors.at(2), "Particle-flow");

    c_pfpi_res_eta[k] = new TCanvas("", "", 800, 600);

    mg_pfpi_res_eta[k]->Draw("APE");
    //DrawAxis(mg_pfpi_res_eta[k], leg_pfpi_res_eta[k], etaMin, etaMax, 0.1, 1000, "#eta", "(resolution in E)/E (%)", false, true);
    DrawAxis(mg_pfpi_res_eta[k], leg_pfpi_res_eta[k], etaMin, etaMax, 0.0, 50, "#eta", "(resolution in E)/E (%)", false, false);
    leg_pfpi_res_eta[k]->Draw();
    pave->Draw();

    TString s_ptrange = Form("pt_%.0f_", ptVals.at(k));

    c_pfpi_res_eta[k]->Print(pdfOutput, "pdf");
    c_pfpi_res_eta[k]->Print(figPath + "img_pfpi_res_" + s_ptrange + "eta.pdf", "pdf");
    c_pfpi_res_eta[k]->Print(figPath + "img_pfpi_res_" + s_ptrange + "eta.png", "png");
  }

  /////////////////
  // PF - jets   //
  /////////////////

  TMultiGraph *mg_pfjet_res_e[n_etabins];
  TMultiGraph *mg_pfjet_res_eta[n_ptbins];

  TLegend *leg_pfjet_res_e[n_etabins];
  TLegend *leg_pfjet_res_eta[n_ptbins];

  TGraphErrors *gr_pfjet_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_pfjet_res_eta = new TGraphErrors[n_ptbins];

  TGraphErrors *gr_cajet_res_e = new TGraphErrors[n_etabins];
  TGraphErrors *gr_cajet_res_eta = new TGraphErrors[n_ptbins];

  std::vector<resolPlot> *plots_pfjet_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_pfjet_res_eta = new std::vector<resolPlot>[n_ptbins];
  std::vector<resolPlot> *plots_cajet_res_e = new std::vector<resolPlot>[n_etabins];
  std::vector<resolPlot> *plots_cajet_res_eta = new std::vector<resolPlot>[n_ptbins];

  TCanvas *c_pfjet_res_e[n_etabins];
  TCanvas *c_pfjet_res_eta[n_ptbins];

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    mg_pfjet_res_e[k] = new TMultiGraph("", "");
    leg_pfjet_res_e[k] = new TLegend(0.40, 0.70, 0.90, 0.90);

    HistogramsCollection(&plots_pfjet_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "pfjet");
    HistogramsCollection(&plots_cajet_res_e[k], TMath::Log10(ptMin), TMath::Log10(ptMax), "cajet");

    GetJetsEres(&plots_pfjet_res_e[k], branchPFJet, branchGenJet, treeReaderJet, etaVals.at(k), etaVals.at(k + 1));
    GetJetsEres(&plots_cajet_res_e[k], branchCaloJet, branchGenJet, treeReaderJet, etaVals.at(k), etaVals.at(k + 1));

    gr_pfjet_res_e[k] = EresGraph(&plots_pfjet_res_e[k]);
    gr_cajet_res_e[k] = EresGraph(&plots_cajet_res_e[k]);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));
    s_eta = "anti-k_{T},  R = 0.4,  " + s_etaMin + " < | #eta | < " + s_etaMax;

    leg_pfjet_res_e[k]->SetTextFont(132);
    leg_pfjet_res_e[k]->SetHeader(s_eta);

    addResoGraph(mg_pfjet_res_e[k], &gr_cajet_res_e[k], leg_pfjet_res_e[k], markerStyles.at(0), colors.at(0), "Calorimeter Jets");
    addResoGraph(mg_pfjet_res_e[k], &gr_pfjet_res_e[k], leg_pfjet_res_e[k], markerStyles.at(1), colors.at(1), "Particle-flow Jets");

    c_pfjet_res_e[k] = new TCanvas("", "", 800, 600);

    mg_pfjet_res_e[k]->Draw("APE");
    //DrawAxis(mg_pfjet_res_e[k], leg_pfjet_res_e[k], 10, ptMax, 0.5, 100, "E [GeV]", "(resolution in E)/E (%)", true, true);
    DrawAxis(mg_pfjet_res_e[k], leg_pfjet_res_e[k], 10, ptMax, 0.0, 30, "E [GeV]", "(resolution in E)/E (%)", true, false);
    leg_pfjet_res_e[k]->Draw();
    pave->Draw();

    TString s_etarange = "eta_" + s_etaMin + "_" + s_etaMax + "_";

    c_pfjet_res_e[k]->Print(pdfOutput, "pdf");
    c_pfjet_res_e[k]->Print(figPath + "img_pfjet_res_" + s_etarange + "e.pdf", "pdf");
    c_pfjet_res_e[k]->Print(figPath + "img_pfjet_res_" + s_etarange + "e.png", "png");
  }

  // loop over eta bins
  for(k = 0; k < ptVals.size(); k++)
  {

    mg_pfjet_res_eta[k] = new TMultiGraph("", "");
    leg_pfjet_res_eta[k] = new TLegend(0.30, 0.70, 0.85, 0.90);

    HistogramsCollectionVsEta(&plots_pfjet_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "pfjet", 0.0, 2.0);
    HistogramsCollectionVsEta(&plots_cajet_res_eta[k], etaMin, etaMax, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), "cajet", 0.0, 2.0);

    GetJetsEresVsEta(&plots_pfjet_res_eta[k], branchPFJet, branchGenJet, treeReaderJet, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k));
    GetJetsEresVsEta(&plots_cajet_res_eta[k], branchCaloJet, branchGenJet, treeReaderJet, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k));

    gr_pfjet_res_eta[k] = EresGraphVsEta(&plots_pfjet_res_eta[k]);
    gr_cajet_res_eta[k] = EresGraphVsEta(&plots_cajet_res_eta[k]);

    s_e = Form("anti-k_{T},  R = 0.4,  jets, E = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_e = Form("anti-k_{T},  R = 0.4,  E = %.0f TeV", ptVals.at(k) / 1000.);

    leg_pfjet_res_eta[k]->SetTextFont(132);
    leg_pfjet_res_eta[k]->SetHeader(s_e);

    addResoGraph(mg_pfjet_res_eta[k], &gr_cajet_res_eta[k], leg_pfjet_res_eta[k], markerStyles.at(0), colors.at(0), "Calorimeter Jets");
    addResoGraph(mg_pfjet_res_eta[k], &gr_pfjet_res_eta[k], leg_pfjet_res_eta[k], markerStyles.at(1), colors.at(1), "Particle-flow Jets");

    c_pfjet_res_eta[k] = new TCanvas("", "", 800, 600);

    mg_pfjet_res_eta[k]->Draw("APE");
    //DrawAxis(mg_pfjet_res_eta[k], leg_pfjet_res_eta[k], etaMin, etaMax, 0.1, 1000, "#eta", "(resolution in E)/E (%)", false, true);
    DrawAxis(mg_pfjet_res_eta[k], leg_pfjet_res_eta[k], etaMin, etaMax, 0.0, 50, "#eta", "(resolution in E)/E (%)", false, false);
    leg_pfjet_res_eta[k]->Draw();
    pave->Draw();

    TString s_ptrange = Form("pt_%.0f_", ptVals.at(k));

    c_pfjet_res_eta[k]->Print(pdfOutput, "pdf");
    c_pfjet_res_eta[k]->Print(figPath + "img_pfjet_res_" + s_ptrange + "eta.pdf", "pdf");
    c_pfjet_res_eta[k]->Print(figPath + "img_pfjet_res_" + s_ptrange + "eta.png", "png");
  }

  /////////////////////
  // PF Missing ET  ///
  /////////////////////

  TMultiGraph *mg_met_res_ht = new TMultiGraph("", "");
  TLegend *leg_met_res_ht = new TLegend(0.60, 0.22, 0.90, 0.42);

  std::vector<resolPlot> plots_pfmet, plots_camet;

  HistogramsCollection(&plots_pfmet, TMath::Log10(ptMin), TMath::Log10(ptMax), "pfMET", -500, 500);
  HistogramsCollection(&plots_camet, TMath::Log10(ptMin), TMath::Log10(ptMax), "caMET", -500, 500);

  GetMetres(&plots_pfmet, branchGenScalarHT, branchMet, branchPFJet, treeReaderJet);
  GetMetres(&plots_camet, branchGenScalarHT, branchCaloMet, branchCaloJet, treeReaderJet);

  TGraphErrors gr_pfmet_res_ht = MetResGraph(&plots_pfmet, true);
  TGraphErrors gr_camet_res_ht = MetResGraph(&plots_camet, true);

  addResoGraph(mg_met_res_ht, &gr_camet_res_ht, leg_met_res_ht, markerStyles.at(0), colors.at(0), "Calorimeter E_{T}^{miss}");
  addResoGraph(mg_met_res_ht, &gr_pfmet_res_ht, leg_met_res_ht, markerStyles.at(1), colors.at(1), "Particle-flow E_{T}^{miss}");

  TCanvas *c_met_res_ht = new TCanvas("", "", 800, 600);

  mg_met_res_ht->Draw("APE");
  DrawAxis(mg_met_res_ht, leg_met_res_ht, 1000, 100000, 0.1, 1000, " #sum p_{T} [GeV]", "resolution in E_{x,y}^{miss} [GeV]", true, true);

  leg_met_res_ht->Draw();
  pave->Draw();
  c_met_res_ht->Print(pdfOutput, "pdf");
  c_met_res_ht->Print(figPath + "img_met_res_ht.pdf", "pdf");
  c_met_res_ht->Print(figPath + "img_met_res_ht.png", "png");

  /////////////////////////////////////////
  // Electron Reconstruction Efficiency ///
  /////////////////////////////////////////

  TMultiGraph *mg_recele_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recele_eff_eta = new TMultiGraph("", "");

  TLegend *leg_recele_eff_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_recele_eff_eta = new TLegend(0.55, 0.22, 0.90, 0.48);

  TGraphErrors *gr_recele_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recele_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recele_eff_pt, *h_recele_eff_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recele_eff_pt = GetEffPt<Electron>(branchElectron, branchParticleElectron, "Electron", 11, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderElectron);
    gr_recele_eff_pt[k] = TGraphErrors(h_recele_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "e^{ #pm} , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recele_eff_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recele_eff_pt, &gr_recele_eff_pt[k], leg_recele_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recele_eff_eta = GetEffEta<Electron>(branchElectron, branchParticleElectron, "Electron", 11, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderElectron);
    gr_recele_eff_eta[k] = TGraphErrors(h_recele_eff_eta);

    s_pt = Form("e^{ #pm} , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("e^{ #pm} , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recele_eff_eta, &gr_recele_eff_eta[k], leg_recele_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recele_eff_pt = new TCanvas("", "", 800, 600);

  mg_recele_eff_pt->Draw("APE");
  DrawAxis(mg_recele_eff_pt, leg_recele_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "reconstruction efficiency (%)", true, false);
  leg_recele_eff_pt->Draw();
  pave->Draw();

  c_recele_eff_pt->Print(pdfOutput, "pdf");
  c_recele_eff_pt->Print(figPath + "img_recele_eff_pt.pdf", "pdf");
  c_recele_eff_pt->Print(figPath + "img_recele_eff_pt.png", "png");

  TCanvas *c_recele_eff_eta = new TCanvas("", "", 800, 600);

  mg_recele_eff_eta->Draw("APE");
  DrawAxis(mg_recele_eff_eta, leg_recele_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "reconstruction efficiency (%)", false, false);
  leg_recele_eff_eta->Draw();
  pave->Draw();

  c_recele_eff_eta->Print(pdfOutput, "pdf");
  c_recele_eff_eta->Print(figPath + "img_recele_eff_eta.pdf", "pdf");
  c_recele_eff_eta->Print(figPath + "img_recele_eff_eta.png", "png");

  /////////////////////////////////////////
  // Muon Reconstruction Efficiency ///
  /////////////////////////////////////////

  TMultiGraph *mg_recmu_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recmu_eff_eta = new TMultiGraph("", "");

  TLegend *leg_recmu_eff_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_recmu_eff_eta = new TLegend(0.55, 0.22, 0.90, 0.48);

  TGraphErrors *gr_recmu_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recmu_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recmu_eff_pt, *h_recmu_eff_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recmu_eff_pt = GetEffPt<Muon>(branchMuon, branchParticleMuon, "muon", 13, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderMuon);
    gr_recmu_eff_pt[k] = TGraphErrors(h_recmu_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#mu^{ #pm} , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recmu_eff_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recmu_eff_pt, &gr_recmu_eff_pt[k], leg_recmu_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recmu_eff_eta = GetEffEta<Muon>(branchMuon, branchParticleMuon, "muon", 13, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderMuon);
    gr_recmu_eff_eta[k] = TGraphErrors(h_recmu_eff_eta);

    s_pt = Form("#mu^{ #pm} , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("#mu^{ #pm} , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recmu_eff_eta, &gr_recmu_eff_eta[k], leg_recmu_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recmu_eff_pt = new TCanvas("", "", 800, 600);

  mg_recmu_eff_pt->Draw("APE");
  DrawAxis(mg_recmu_eff_pt, leg_recmu_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "reconstruction efficiency (%)", true, false);
  leg_recmu_eff_pt->Draw();
  pave->Draw();

  c_recmu_eff_pt->Print(pdfOutput, "pdf");
  c_recmu_eff_pt->Print(figPath + "img_recmu_eff_pt.pdf", "pdf");
  c_recmu_eff_pt->Print(figPath + "img_recmu_eff_pt.png", "png");

  TCanvas *c_recmu_eff_eta = new TCanvas("", "", 800, 600);

  mg_recmu_eff_eta->Draw("APE");
  DrawAxis(mg_recmu_eff_eta, leg_recmu_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "reconstruction efficiency (%)", false, false);
  leg_recmu_eff_eta->Draw();
  pave->Draw();

  c_recmu_eff_eta->Print(pdfOutput, "pdf");
  c_recmu_eff_eta->Print(figPath + "img_recmu_eff_eta.pdf", "pdf");
  c_recmu_eff_eta->Print(figPath + "img_recmu_eff_eta.png", "png");

  /////////////////////////////////////////
  // Photon Reconstruction Efficiency   ///
  /////////////////////////////////////////

  TMultiGraph *mg_recpho_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recpho_eff_eta = new TMultiGraph("", "");

  TLegend *leg_recpho_eff_pt = new TLegend(0.55, 0.22, 0.90, 0.48);
  TLegend *leg_recpho_eff_eta = new TLegend(0.55, 0.22, 0.90, 0.48);

  TGraphErrors *gr_recpho_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recpho_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recpho_eff_pt, *h_recpho_eff_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recpho_eff_pt = GetEffPt<Photon>(branchPhoton, branchParticlePhoton, "Photon", 22, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderPhoton);
    gr_recpho_eff_pt[k] = TGraphErrors(h_recpho_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#gamma , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recpho_eff_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recpho_eff_pt, &gr_recpho_eff_pt[k], leg_recpho_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recpho_eff_eta = GetEffEta<Photon>(branchPhoton, branchParticlePhoton, "Photon", 22, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderPhoton);
    gr_recpho_eff_eta[k] = TGraphErrors(h_recpho_eff_eta);

    s_pt = Form("#gamma , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("#gamma , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recpho_eff_eta, &gr_recpho_eff_eta[k], leg_recpho_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recpho_eff_pt = new TCanvas("", "", 800, 600);

  mg_recpho_eff_pt->Draw("APE");
  DrawAxis(mg_recpho_eff_pt, leg_recpho_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "reconstruction efficiency (%)", true, false);
  leg_recpho_eff_pt->Draw();
  pave->Draw();

  c_recpho_eff_pt->Print(pdfOutput, "pdf");
  c_recpho_eff_pt->Print(figPath + "img_recpho_eff_pt.pdf", "pdf");
  c_recpho_eff_pt->Print(figPath + "img_recpho_eff_pt.png", "png");

  TCanvas *c_recpho_eff_eta = new TCanvas("", "", 800, 600);

  mg_recpho_eff_eta->Draw("APE");
  DrawAxis(mg_recpho_eff_eta, leg_recpho_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "reconstruction efficiency (%)", false, false);
  leg_recpho_eff_eta->Draw();
  pave->Draw();

  c_recpho_eff_eta->Print(pdfOutput, "pdf");
  c_recpho_eff_eta->Print(figPath + "img_recpho_eff_eta.pdf", "pdf");
  c_recpho_eff_eta->Print(figPath + "img_recpho_eff_eta.png", "png");

  /////////////////////////////////////////
  // B-jets  Efficiency/ mistag rates   ///
  /////////////////////////////////////////

  TMultiGraph *mg_recbjet_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recbjet_eff_eta = new TMultiGraph("", "");

  TLegend *leg_recbjet_eff_pt = new TLegend(0.50, 0.22, 0.90, 0.48);
  TLegend *leg_recbjet_eff_eta = new TLegend(0.50, 0.22, 0.90, 0.48);

  TGraphErrors *gr_recbjet_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recbjet_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recbjet_eff_pt, *h_recbjet_eff_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recbjet_eff_pt = GetJetEffPt<Jet>(branchPFBJet, "BJet", 5, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderBJet);
    //h_recbjet_eff_pt = GetEffPt<Jet>(branchPFBJet, branchParticleBJet, "BJet", 5, ptMin, ptMax, etaVals.at(k), etaVals.at(k+1), treeReaderBJet);
    gr_recbjet_eff_pt[k] = TGraphErrors(h_recbjet_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "b-jet , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recbjet_eff_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recbjet_eff_pt, &gr_recbjet_eff_pt[k], leg_recbjet_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recbjet_eff_eta = GetJetEffEta<Jet>(branchPFBJet, "BJet", 5, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderBJet);
    //h_recbjet_eff_eta = GetEffEta<Jet>(branchPFBJet, branchParticleBJet, "BJet", 5, 0.5*ptVals.at(k), 2.0*ptVals.at(k) ,etaMin, etaMax , treeReaderBJet);
    gr_recbjet_eff_eta[k] = TGraphErrors(h_recbjet_eff_eta);

    s_pt = Form("b-jet , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("b-jet , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recbjet_eff_eta, &gr_recbjet_eff_eta[k], leg_recbjet_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recbjet_eff_pt = new TCanvas("", "", 800, 600);

  mg_recbjet_eff_pt->Draw("APE");
  DrawAxis(mg_recbjet_eff_pt, leg_recbjet_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "b - tag efficiency (%)", true, false);
  leg_recbjet_eff_pt->Draw();
  pave->Draw();

  c_recbjet_eff_pt->Print(pdfOutput, "pdf");
  c_recbjet_eff_pt->Print(figPath + "img_recbjet_eff_pt.pdf", "pdf");
  c_recbjet_eff_pt->Print(figPath + "img_recbjet_eff_pt.png", "png");

  TCanvas *c_recbjet_eff_eta = new TCanvas("", "", 800, 600);

  mg_recbjet_eff_eta->Draw("APE");
  DrawAxis(mg_recbjet_eff_eta, leg_recbjet_eff_eta, etaMin, etaMax, 0.0, 100, " #eta ", "b - tag efficiency (%)", false, false);
  leg_recbjet_eff_eta->Draw();
  pave->Draw();

  c_recbjet_eff_eta->Print(pdfOutput, "pdf");
  c_recbjet_eff_eta->Print(figPath + "img_recbjet_eff_eta.pdf", "pdf");
  c_recbjet_eff_eta->Print(figPath + "img_recbjet_eff_eta.png", "png");

  // ------ c - mistag  ------

  TMultiGraph *mg_recbjet_cmis_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recbjet_cmis_eta = new TMultiGraph("", "");

  TLegend *leg_recbjet_cmis_pt = new TLegend(0.50, 0.64, 0.90, 0.90);
  TLegend *leg_recbjet_cmis_eta = new TLegend(0.50, 0.64, 0.90, 0.90);

  TGraphErrors *gr_recbjet_cmis_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recbjet_cmis_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recbjet_cmis_pt, *h_recbjet_cmis_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recbjet_cmis_pt = GetJetEffPt<Jet>(branchPFCJet, "CJet", 4, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderCJet);
    //h_recbjet_cmis_pt = GetEffPt<Jet>(branchPFCJet, branchParticleCJet, "CJet", 4, ptMin, ptMax, etaVals.at(k), etaVals.at(k+1), treeReaderCJet);
    gr_recbjet_cmis_pt[k] = TGraphErrors(h_recbjet_cmis_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "c-jet , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recbjet_cmis_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recbjet_cmis_pt, &gr_recbjet_cmis_pt[k], leg_recbjet_cmis_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recbjet_cmis_eta = GetJetEffEta<Jet>(branchPFCJet, "CJet", 4, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderCJet);
    //h_recbjet_cmis_eta = GetEffEta<Jet>(branchPFCJet, branchParticleCJet, "CJet", 4, 0.5*ptVals.at(k), 2.0*ptVals.at(k) ,etaMin, etaMax , treeReaderCJet);
    gr_recbjet_cmis_eta[k] = TGraphErrors(h_recbjet_cmis_eta);

    s_pt = Form("c-jet , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("c-jet , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recbjet_cmis_eta, &gr_recbjet_cmis_eta[k], leg_recbjet_cmis_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recbjet_cmis_pt = new TCanvas("", "", 800, 600);

  mg_recbjet_cmis_pt->Draw("APE");
  DrawAxis(mg_recbjet_cmis_pt, leg_recbjet_cmis_pt, ptMin, ptMax, 0.0, 20, "p_{T} [GeV]", "c - mistag rate (%)", true, false);
  leg_recbjet_cmis_pt->Draw();
  pave->Draw();

  c_recbjet_cmis_pt->Print(pdfOutput, "pdf");
  c_recbjet_cmis_pt->Print(figPath + "img_recbjet_cmis_pt.pdf", "pdf");
  c_recbjet_cmis_pt->Print(figPath + "img_recbjet_cmis_pt.png", "png");

  TCanvas *c_recbjet_cmis_eta = new TCanvas("", "", 800, 600);

  mg_recbjet_cmis_eta->Draw("APE");
  DrawAxis(mg_recbjet_cmis_eta, leg_recbjet_cmis_eta, etaMin, etaMax, 0.0, 20, " #eta ", "c - mistag rate (%)", false, false);
  leg_recbjet_cmis_eta->Draw();
  pave->Draw();

  c_recbjet_cmis_eta->Print(pdfOutput, "pdf");
  c_recbjet_cmis_eta->Print(figPath + "img_recbjet_cmis_eta.pdf", "pdf");
  c_recbjet_cmis_eta->Print(figPath + "img_recbjet_cmis_eta.png", "png");

  // ------ light - mistag  ------

  TMultiGraph *mg_recbjet_lmis_pt = new TMultiGraph("", "");
  TMultiGraph *mg_recbjet_lmis_eta = new TMultiGraph("", "");

  TLegend *leg_recbjet_lmis_pt = new TLegend(0.50, 0.64, 0.90, 0.90);
  TLegend *leg_recbjet_lmis_eta = new TLegend(0.50, 0.64, 0.90, 0.90);

  TGraphErrors *gr_recbjet_lmis_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_recbjet_lmis_eta = new TGraphErrors[n_ptbins];
  TH1D *h_recbjet_lmis_pt, *h_recbjet_lmis_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_recbjet_lmis_pt = GetJetEffPt<Jet>(branchJet, "Jet", 1, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderJet);
    //h_recbjet_lmis_pt = GetEffPt<Jet>(branchJet, branchParticleJet, "Jet", 1, ptMin, ptMax, etaVals.at(k), etaVals.at(k+1), treeReaderJet);
    gr_recbjet_lmis_pt[k] = TGraphErrors(h_recbjet_lmis_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "uds-jet , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_recbjet_lmis_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_recbjet_lmis_pt, &gr_recbjet_lmis_pt[k], leg_recbjet_lmis_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_recbjet_lmis_eta = GetJetEffEta<Jet>(branchJet, "Jet", 1, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderJet);
    //h_recbjet_lmis_eta = GetEffEta<Jet>(branchJet, branchParticleJet, "Jet", 1, 0.5*ptVals.at(k), 2.0*ptVals.at(k) ,etaMin, etaMax , treeReaderJet);
    gr_recbjet_lmis_eta[k] = TGraphErrors(h_recbjet_lmis_eta);

    s_pt = Form("uds-jet , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("uds-jet , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_recbjet_lmis_eta, &gr_recbjet_lmis_eta[k], leg_recbjet_lmis_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_recbjet_lmis_pt = new TCanvas("", "", 800, 600);

  mg_recbjet_lmis_pt->Draw("APE");

  DrawAxis(mg_recbjet_lmis_pt, leg_recbjet_lmis_pt, ptMin, ptMax, 0.0, 1.0, "p_{T} [GeV]", "light - mistag rate (%)", true, false);

  leg_recbjet_lmis_pt->Draw();
  pave->Draw();

  c_recbjet_lmis_pt->Print(pdfOutput, "pdf");
  c_recbjet_lmis_pt->Print(figPath + "img_recbjet_lmis_pt.pdf", "pdf");
  c_recbjet_lmis_pt->Print(figPath + "img_recbjet_lmis_pt.png", "png");

  TCanvas *c_recbjet_lmis_eta = new TCanvas("", "", 800, 600);

  mg_recbjet_lmis_eta->Draw("APE");
  DrawAxis(mg_recbjet_lmis_eta, leg_recbjet_lmis_eta, etaMin, etaMax, 0.0, 1.0, " #eta ", "light - mistag rate (%)", false, false);
  leg_recbjet_lmis_eta->Draw();
  pave->Draw();

  c_recbjet_lmis_eta->Print(pdfOutput, "pdf");
  c_recbjet_lmis_eta->Print(figPath + "img_recbjet_lmis_eta.pdf", "pdf");
  c_recbjet_lmis_eta->Print(figPath + "img_recbjet_lmis_eta.png", "png");

  ///////////////////////////////////////////
  // tau-jets  Efficiency/ mistag rates   ///
  ///////////////////////////////////////////

  TMultiGraph *mg_rectaujet_eff_pt = new TMultiGraph("", "");
  TMultiGraph *mg_rectaujet_eff_eta = new TMultiGraph("", "");

  TLegend *leg_rectaujet_eff_pt = new TLegend(0.50, 0.22, 0.90, 0.48);
  TLegend *leg_rectaujet_eff_eta = new TLegend(0.50, 0.22, 0.90, 0.48);

  TGraphErrors *gr_rectaujet_eff_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_rectaujet_eff_eta = new TGraphErrors[n_ptbins];
  TH1D *h_rectaujet_eff_pt, *h_rectaujet_eff_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_rectaujet_eff_pt = GetTauEffPt<Jet>(branchPFTauJet, branchParticleTauJet, "TauJet", 15, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderTauJet);
    gr_rectaujet_eff_pt[k] = TGraphErrors(h_rectaujet_eff_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "#tau-jet , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_rectaujet_eff_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_rectaujet_eff_pt, &gr_rectaujet_eff_pt[k], leg_rectaujet_eff_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_rectaujet_eff_eta = GetTauEffEta<Jet>(branchPFTauJet, branchParticleTauJet, "TauJet", 15, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderTauJet);
    gr_rectaujet_eff_eta[k] = TGraphErrors(h_rectaujet_eff_eta);

    s_pt = Form("#tau-jet , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("#tau-jet , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_rectaujet_eff_eta, &gr_rectaujet_eff_eta[k], leg_rectaujet_eff_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_rectaujet_eff_pt = new TCanvas("", "", 800, 600);

  mg_rectaujet_eff_pt->Draw("APE");
  DrawAxis(mg_rectaujet_eff_pt, leg_rectaujet_eff_pt, ptMin, ptMax, 0.0, 100, "p_{T} [GeV]", "#tau - tag efficiency (%)", true, false);
  leg_rectaujet_eff_pt->Draw();
  pave->Draw();

  c_rectaujet_eff_pt->Print(pdfOutput, "pdf");
  c_rectaujet_eff_pt->Print(figPath + "img_rectaujet_eff_pt.pdf", "pdf");
  c_rectaujet_eff_pt->Print(figPath + "img_rectaujet_eff_pt.png", "png");

  TCanvas *c_rectaujet_eff_eta = new TCanvas("", "", 800, 600);

  mg_rectaujet_eff_eta->Draw("APE");
  DrawAxis(mg_rectaujet_eff_eta, leg_rectaujet_eff_eta, etaMin, etaMax, 0.0, 100., " #eta ", "#tau - tag efficiency (%)", false, false);
  leg_rectaujet_eff_eta->Draw();
  pave->Draw();

  c_rectaujet_eff_eta->Print(pdfOutput, "pdf");
  c_rectaujet_eff_eta->Print(figPath + "img_rectaujet_eff_eta.pdf", "pdf");
  c_rectaujet_eff_eta->Print(figPath + "img_rectaujet_eff_eta.png", "png");

  //--------------- tau mistag rate ----------

  TMultiGraph *mg_rectaujet_mis_pt = new TMultiGraph("", "");
  TMultiGraph *mg_rectaujet_mis_eta = new TMultiGraph("", "");

  TLegend *leg_rectaujet_mis_pt = new TLegend(0.50, 0.64, 0.90, 0.90);
  TLegend *leg_rectaujet_mis_eta = new TLegend(0.50, 0.64, 0.90, 0.90);

  TGraphErrors *gr_rectaujet_mis_pt = new TGraphErrors[n_etabins];
  TGraphErrors *gr_rectaujet_mis_eta = new TGraphErrors[n_ptbins];
  TH1D *h_rectaujet_mis_pt, *h_rectaujet_mis_eta;

  // loop over eta bins
  for(k = 0; k < etaVals.size() - 1; k++)
  {

    h_rectaujet_mis_pt = GetTauEffPt<Jet>(branchJet, branchParticleJet, "TauJet", 1, ptMin, ptMax, etaVals.at(k), etaVals.at(k + 1), treeReaderJet);
    gr_rectaujet_mis_pt[k] = TGraphErrors(h_rectaujet_mis_pt);

    s_etaMin = Form("%.1f", etaVals.at(k));
    s_etaMax = Form("%.1f", etaVals.at(k + 1));

    s_eta = "uds-jet , " + s_etaMin + " < | #eta | < " + s_etaMax;

    gr_rectaujet_mis_pt[k].SetName("recEff_" + s_etaMin + "_" + s_etaMax);

    addResoGraph(mg_rectaujet_mis_pt, &gr_rectaujet_mis_pt[k], leg_rectaujet_mis_pt, markerStyles.at(k), colors.at(k), s_eta);
  }

  // loop over pt
  for(k = 0; k < ptVals.size(); k++)
  {
    h_rectaujet_mis_eta = GetTauEffEta<Jet>(branchJet, branchParticleJet, "TauJet", 1, 0.5 * ptVals.at(k), 2.0 * ptVals.at(k), etaMin, etaMax, treeReaderJet);
    gr_rectaujet_mis_eta[k] = TGraphErrors(h_rectaujet_mis_eta);

    s_pt = Form("uds-jet , p_{T} = %.0f GeV", ptVals.at(k));
    if(ptVals.at(k) >= 1000.) s_pt = Form("uds-jet , p_{T} = %.0f TeV", ptVals.at(k) / 1000.);

    addResoGraph(mg_rectaujet_mis_eta, &gr_rectaujet_mis_eta[k], leg_rectaujet_mis_eta, markerStyles.at(k), colors.at(k), s_pt);
  }

  TCanvas *c_rectaujet_mis_pt = new TCanvas("", "", 800, 600);

  mg_rectaujet_mis_pt->Draw("APE");
  DrawAxis(mg_rectaujet_mis_pt, leg_rectaujet_mis_pt, ptMin, ptMax, 0.0, 5., "p_{T} [GeV]", "#tau - mistag(%)", true, false);
  leg_rectaujet_mis_pt->Draw();
  pave->Draw();

  c_rectaujet_mis_pt->Print(pdfOutput, "pdf");
  c_rectaujet_mis_pt->Print(figPath + "img_rectaujet_mis_pt.pdf", "pdf");
  c_rectaujet_mis_pt->Print(figPath + "img_rectaujet_mis_pt.png", "png");

  TCanvas *c_rectaujet_mis_eta = new TCanvas("", "", 800, 600);

  mg_rectaujet_mis_eta->Draw("APE");
  DrawAxis(mg_rectaujet_mis_eta, leg_rectaujet_mis_eta, etaMin, etaMax, 0.0, 5., " #eta ", "#tau - mistag (%)", false, false);
  leg_rectaujet_mis_eta->Draw();
  pave->Draw();

  c_rectaujet_mis_eta->Print(pdfOutput + ")", "pdf");
  c_rectaujet_mis_eta->Print(figPath + "img_rectaujet_mis_eta.pdf", "pdf");
  c_rectaujet_mis_eta->Print(figPath + "img_rectaujet_mis_eta.png", "png");

  //   ----   store resolution histograms in the output (for leave efficiencies out) ---

  TFile *fout = new TFile(outputFile, "recreate");

  for(int bin = 0; bin < Nbins; bin++)
  {

    for(k = 0; k < etaVals.size() - 1; k++)
    {
      plots_trkpi_res_pt[k].at(bin).resolHist->Write();
      plots_trkele_res_pt[k].at(bin).resolHist->Write();
      plots_trkmu_res_pt[k].at(bin).resolHist->Write();
      plots_ecal_res_e[k].at(bin).resolHist->Write();
      plots_hcal_res_e[k].at(bin).resolHist->Write();
      plots_pfele_res_e[k].at(bin).resolHist->Write();
      plots_pfpi_res_e[k].at(bin).resolHist->Write();
      plots_pfjet_res_e[k].at(bin).resolHist->Write();
      plots_cajet_res_e[k].at(bin).resolHist->Write();
    }

    for(k = 0; k < ptVals.size(); k++)
    {
      plots_trkpi_res_eta[k].at(bin).resolHist->Write();
      plots_trkele_res_eta[k].at(bin).resolHist->Write();
      plots_trkmu_res_eta[k].at(bin).resolHist->Write();
      plots_ecal_res_eta[k].at(bin).resolHist->Write();
      plots_hcal_res_eta[k].at(bin).resolHist->Write();
      plots_pfele_res_eta[k].at(bin).resolHist->Write();
      plots_pfpi_res_eta[k].at(bin).resolHist->Write();
      plots_pfjet_res_eta[k].at(bin).resolHist->Write();
      plots_cajet_res_eta[k].at(bin).resolHist->Write();
    }

    plots_pfmet.at(bin).resolHist->Write();
    plots_camet.at(bin).resolHist->Write();
  }

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
  char *appName = "DelphesValidation";

  if(argc != 12)
  {
    cout << " Usage: " << appName << " input_file_electron input_file_muon input_file_photon input_file_jet input_file_bjet input_file_taujet output_file version" << endl;
    cout << " input_file_pion - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_electron - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_muon - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_photon - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_neutralhadron - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_jet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_bjet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_cjet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " input_file_taujet - input file in ROOT format ('Delphes' tree)," << endl;
    cout << " output_file - output file in ROOT format," << endl;
    cout << " version - Delphes version" << endl;

    return 1;
  }

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  DelphesValidation(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11]);
}
