#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include <TH1.h>
#include "TString.h"
#include "vector"
#include <TMath.h>
#include <iostream>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <typeinfo>
#include "TLorentzVector.h"

//------------------------------------------------------------------------------

double ptrangemin = 10;
double ptrangemax = 10000;
static const int Nbins = 20;

int objStyle = 1;
int trackStyle = 7;
int towerStyle = 5;

Color_t objColor = kBlack;
Color_t trackColor = kBlack;
Color_t towerColor = kBlack;

double effLegXmin = 0.22;
double effLegXmax = 0.5;
double effLegYmin = 0.22;
double effLegYmax = 0.4;

struct resolPlot
{
    TH1 *cenResolHist;
    TH1 *fwdResolHist;
    double ptmin;
    double ptmax;
    TString obj;

    resolPlot();
    resolPlot(double ptdown, double ptup, TString object);
    void set(double ptdown, double ptup, TString object);
    void print(){std::cout << ptmin << std::endl;}
};


resolPlot::resolPlot()
{
}

resolPlot::resolPlot(double ptdown, double ptup, TString object)
{
    this->set(ptdown,ptup,object);
}

void resolPlot::set(double ptdown, double ptup, TString object){
    ptmin = ptdown;
    ptmax = ptup;
    obj = object;

    cenResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", 500,  0, 2);
    fwdResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", 500, 0.4, 0.4);

}

void HistogramsCollection(std::vector<resolPlot> *histos, double ptmin, double ptmax, TString obj)
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
        ptemp.set(ptdown, ptup, obj);
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

   for (int i = 0; i <= bins; i++) {
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

  TH1D *histGenPt = new TH1D(name+" gen spectra Pt",name+" gen spectra Pt", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D *histRecoPt = new TH1D(name+" reco spectra Pt",name+" reco spectra Pt", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D *histGenEta  = new TH1D(name+" gen spectra Eta",name+" gen spectra Eta", 12, -3, 3);
  TH1D *histRecoEta = new TH1D(name+" reco spectra Eta",name+" reco spectra Eta", 12, -3, 3);


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

        histGenPt->Fill(pt);
        histGenEta->Fill(eta);

        if(deltaR < 0.3)
        {
          histRecoPt->Fill(pt);
          histRecoEta->Fill(eta); 
        }
      }
    }
  }


  std::pair<TH1D*,TH1D*> histos; 

  histRecoPt->Divide(histGenPt);
  histRecoEta->Divide(histGenEta);

  histos.first = histRecoPt;
  histos.second = histRecoEta;

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
        pt  = bestGenMomentum.Pt();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
            if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta > 0.0 && eta < 2.5) 
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
        pt  = genJetMomentum.Pt();
        eta = TMath::Abs(genJetMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
            if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta > 0.0 && eta < 1.5) 
            {
                histos->at(bin).cenResolHist->Fill(jetMomentum.Pt()/bestGenJetMomentum.Pt());
            }
        }
      }
    }
  }
}
std::pair<Double_t, Double_t> GausFit(TH1* hist)
{
    TF1 *f1 = new TF1("f1", "gaus", hist->GetMean()-2*hist->GetRMS(), hist->GetMean()+2*hist->GetRMS());
    hist->Fit("f1","RQ");
    Double_t sig = f1->GetParameter(2);
    Double_t sigErr = f1->GetParError(2);
    delete f1;
    return make_pair (sig, sigErr);
    //return make_pair (hist->GetRMS(), hist->GetRMSError());
}


TGraphErrors EresGraph(std::vector<resolPlot> *histos, bool central)
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
            std::cout << " pt : " << (histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0;
            std::cout << " mean : " << histos->at(bin).cenResolHist->GetMean() << " RMS : " << histos->at(bin).cenResolHist->GetRMS(); 
            std::cout << " entries : " << histos->at(bin).cenResolHist->GetEntries() << std::endl;
            std::pair<Double_t, Double_t> sigvalues = GausFit(histos->at(bin).cenResolHist);
            gr.SetPoint(count,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sigvalues.first);
            gr.SetPointError(count,0, sigvalues.second);
            count++;
        }
        /*
        else if (histos->at(bin).cenResolHist->GetEntries() > 10) 
        {
            histos->at(bin).fwdResolHist->Fit("gaus","","", -2*histos->at(bin).fwdResolHist->GetRMS(), 2*histos->at(bin).fwdResolHist->GetRMS());
            TF1 *f = histos->at(bin).fwdResolHist->GetFunction("gaus");
            Double_t sig = f->GetParameter(2);
            Double_t sigErr = f->GetParError(2);
            gr.SetPoint(bin,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, sig);
            gr.SetPointError(bin,0, sigErr);
        }
        */
    }
    return gr;
}


//------------------------------------------------------------------------------


// type 1 : object, 2 : track, 3 : tower

void addGraph(TMultiGraph *mg, TGraphErrors *gr, TLegend *leg, int type)
{

  gr->SetLineWidth(3);

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

    default:
      std::cout << "wrong type, possibles choices are Object, Track and Tower" << std::endl;
      break;
  }
}

void addHist(TH1D *h, TLegend *leg, int type)
{
  h->SetLineWidth(3);

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

    default:
      std::cout << "wrong type, possibles choices are Object, Track and Tower" << std::endl;
      break;
  }
}

void DrawAxis(TMultiGraph *mg, TLegend *leg, double max)
{
  mg->SetMinimum(0.);
  mg->SetMaximum(max);
  mg->GetXaxis()->SetLimits(ptrangemin,ptrangemax);
  mg->GetYaxis()->SetTitle("#sigma (E) / E_{gen}");
  mg->GetXaxis()->SetTitle("E_{gen}");
  mg->GetYaxis()->SetTitleSize(0.07);
  mg->GetXaxis()->SetTitleSize(0.07);
  mg->GetYaxis()->SetLabelSize(0.07);
  mg->GetXaxis()->SetLabelSize(0.07);

  leg->SetLineStyle(0);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);

  gStyle->SetOptTitle(0); 
  gPad->SetLogx();
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);
  gPad->Modified();
  gPad->Update();
  
}

void DrawAxis(TH1D *h, TLegend *leg, int type)
{

  h->GetYaxis()->SetRangeUser(0,1.2);
  if (type == 0) h->GetXaxis()->SetTitle("E_{gen}");
  else h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("#epsilon");
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->SetLabelSize(0.07);
  leg->SetLineStyle(0);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.2);
  gPad->SetLeftMargin(0.2);

  gPad->Modified();
  gPad->Update();
  
}


void Validation(const char *inputFile, const char *outputFile)
{
  //gSystem->Load("libDelphes");

  std::cout << "input file : " << inputFile << " " << " , output file : " << outputFile << std::endl;

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchPFJet = treeReader->UseBranch("Jet");
  TClonesArray *branchCaloJet = treeReader->UseBranch("CaloJet");

#ifdef ELECTRON

  ///////////////
  // Electrons //
  ///////////////

  // Reconstruction efficiency
  TString elecs = "Electron";
  int elID = 11;
  std::pair<TH1D*,TH1D*> histos_el = GetEff<Electron>(branchElectron, branchParticle, "Electron", elID, treeReader);

  // tracking reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_eltrack = GetEff<Track>(branchTrack, branchParticle, "electronTrack", elID, treeReader);

  // Tower reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_eltower = GetEff<Tower>(branchTower, branchParticle, "electronTower", elID, treeReader);

  // Electron Energy Resolution
  std::vector<resolPlot> plots_el;
  HistogramsCollection(&plots_el, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electrons");
  GetEres<Electron>( &plots_el, branchElectron, branchParticle, elID, treeReader);
  TGraphErrors gr_el = EresGraph(&plots_el, true);
  gr_el.SetName("Electron");

  // Electron Track Energy Resolution
  std::vector<resolPlot> plots_eltrack;
  HistogramsCollection(&plots_eltrack, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTracks");
  GetEres<Track>( &plots_eltrack, branchTrack, branchParticle, elID, treeReader);
  TGraphErrors gr_eltrack = EresGraph(&plots_eltrack, true);
  gr_eltrack.SetName("ElectronTracks");

  // Electron Tower Energy Resolution
  std::vector<resolPlot> plots_eltower;
  HistogramsCollection(&plots_eltower, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTower");
  GetEres<Tower>( &plots_eltower, branchTower, branchParticle, elID, treeReader);
  TGraphErrors gr_eltower = EresGraph(&plots_eltower, true);
  gr_eltower.SetName("ElectronTower");

  // Canvases
  TString elEff = "electronEff";
  TCanvas *C_el1 = new TCanvas(elEff,elEff, 1000, 500);
  C_el1->Divide(2);
  C_el1->cd(1);
  TLegend *leg_el1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  gPad->SetLogx();
  histos_eltrack.first->Draw();
  addHist(histos_eltrack.first, leg_el1, 2); 
  histos_eltower.first->Draw("same");
  addHist(histos_eltower.first, leg_el1, 3);
  histos_el.first->Draw("same");
  addHist(histos_el.first, leg_el1, 1);

  DrawAxis(histos_eltrack.first, leg_el1, 0);
  leg_el1->Draw();

  C_el1->cd(2);
  TLegend *leg_el2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  histos_eltrack.second->Draw();
  addHist(histos_eltrack.second, leg_el2, 2);
  histos_eltower.second->Draw("same");
  addHist(histos_eltower.second, leg_el2, 3);
  histos_el.second->Draw("same");
  addHist(histos_el.second, leg_el2, 1);

  DrawAxis(histos_eltrack.second, leg_el2, 1);
  delete(leg_el2);
  C_el1->cd(0);
 
  TString elRes = "electronERes";
  TCanvas *C_el2 = new TCanvas(elRes,elRes, 1000, 500);
  TMultiGraph *mg_el = new TMultiGraph(elRes,elRes);
  TLegend *leg_el = new TLegend(0.52,0.7,0.9,0.9);

  addGraph(mg_el, &gr_eltower, leg_el, 3);
  addGraph(mg_el, &gr_eltrack, leg_el, 2);
  addGraph(mg_el, &gr_el, leg_el, 1);

  mg_el->Draw("ACX");
  leg_el->Draw();

  DrawAxis(mg_el, leg_el, 0.1);
  
  C_el1->SaveAs(elEff+".eps");
  C_el2->SaveAs(elRes+".eps");

#endif

  gDirectory->cd(0);

#ifdef MUON

  ///////////
  // Muons //
  ///////////

  // Reconstruction efficiency
  int muID = 13;
  std::pair<TH1D*,TH1D*> histos_mu = GetEff<Muon>(branchMuon, branchParticle,"Muon", muID, treeReader);

  // muon tracking reconstruction efficiency
  std::pair <TH1D*,TH1D*> histos_mutrack = GetEff<Track>(branchTrack, branchParticle, "muonTrack", muID, treeReader);

  // Muon Energy Resolution
  std::vector<resolPlot> plots_mu;
  HistogramsCollection(&plots_mu, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "muons");
  GetEres<Muon>( &plots_mu, branchMuon, branchParticle, muID, treeReader);
  TGraphErrors gr_mu = EresGraph(&plots_mu, true);
  gr_mu.SetName("Muon");

  // Muon Track Energy Resolution
  std::vector<resolPlot> plots_mutrack;
  HistogramsCollection(&plots_mutrack, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "muonsTracks");
  GetEres<Track>( &plots_mutrack, branchTrack, branchParticle, muID, treeReader);
  TGraphErrors gr_mutrack = EresGraph(&plots_mutrack, true);
  gr_mutrack.SetName("MuonTracks");

  // Canvas

  TString muEff = "muonEff";
  TCanvas *C_mu1 = new TCanvas(muEff,muEff, 1000, 500);
  C_mu1->Divide(2);
  C_mu1->cd(1);
  TLegend *leg_mu1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  gPad->SetLogx();
  histos_mutrack.first->Draw();
  addHist(histos_mutrack.first, leg_mu1, 2);
  histos_mu.first->Draw("same");
  addHist(histos_mu.first, leg_mu1, 1);

  DrawAxis(histos_mutrack.first, leg_mu1, 0);
  leg_mu1->Draw();

  C_mu1->cd(2);
  TLegend *leg_mu2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  histos_mutrack.second->Draw();
  addHist(histos_mutrack.second, leg_mu2, 2);
  histos_mu.second->Draw("same");
  addHist(histos_mu.second, leg_mu2, 1);

  DrawAxis(histos_mutrack.second, leg_mu2, 1);

  TString muRes = "muonERes";
  TCanvas *C_mu = new TCanvas(muRes,muRes, 1000, 500);
  TMultiGraph *mg_mu = new TMultiGraph(muRes,muRes);
  TLegend *leg_mu = new TLegend(0.52,0.7,0.9,0.9);

  addGraph(mg_mu, &gr_mutrack, leg_mu, 2);
  addGraph(mg_mu, &gr_mu, leg_mu, 1);

  mg_mu->Draw("ACX");
  leg_mu->Draw();

  DrawAxis(mg_mu, leg_mu, 0.3);

  C_mu1->SaveAs(muEff+".eps");
  C_mu->SaveAs(muRes+".eps");

#endif

  gDirectory->cd(0);

#ifdef PHOTON

  /////////////
  // Photons //
  /////////////

  // Reconstruction efficiency
  int phID = 22;
  std::pair<TH1D*,TH1D*> histos_ph = GetEff<Electron>(branchPhoton, branchParticle, "Photon", phID, treeReader);
  std::pair<TH1D*,TH1D*> histos_phtower = GetEff<Electron>(branchTower, branchParticle, "Photon", phID, treeReader);

  // Photon Energy Resolution
  std::vector<resolPlot> plots_ph;
  HistogramsCollection(&plots_ph, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "photons");
  GetEres<Photon>( &plots_ph, branchPhoton, branchParticle, phID, treeReader);
  TGraphErrors gr_ph = EresGraph(&plots_ph, true);
  gr_ph.SetName("Photon");

  // Photon Tower Energy Resolution
  std::vector<resolPlot> plots_phtower;
  HistogramsCollection(&plots_phtower, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "photonsTower");
  GetEres<Tower>( &plots_phtower, branchTower, branchParticle, phID, treeReader);
  TGraphErrors gr_phtower = EresGraph(&plots_phtower, true);
  gr_phtower.SetName("PhotonTower");

  // Canvas

  TString phEff = "photonEff";
  TCanvas *C_ph1 = new TCanvas(phEff,phEff, 1000, 500);
  C_ph1->Divide(2);
  C_ph1->cd(1);
  TLegend *leg_ph1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  gPad->SetLogx();
  histos_phtower.first->Draw();
  addHist(histos_phtower.first, leg_ph1, 3);
  histos_ph.first->Draw("same");
  addHist(histos_ph.first, leg_ph1, 1);

  DrawAxis(histos_phtower.first, leg_ph1, 0);
  leg_ph1->Draw();

  C_ph1->cd(2);
  TLegend *leg_ph2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  histos_phtower.second->Draw("same");
  addHist(histos_phtower.second, leg_ph2, 3);
  histos_ph.second->Draw("same");
  addHist(histos_ph.second, leg_ph2, 1);

  DrawAxis(histos_phtower.second, leg_ph2, 1);

  C_ph1->SaveAs(phEff+".eps");

  TString phRes = "phERes";
  TCanvas *C_ph = new TCanvas(phRes,phRes, 1000, 500);
  TMultiGraph *mg_ph = new TMultiGraph(phRes,phRes);
  TLegend *leg_ph = new TLegend(0.52,0.7,0.9,0.9);

  addGraph(mg_ph, &gr_phtower, leg_ph, 3);
  addGraph(mg_ph, &gr_ph, leg_ph, 1);

  mg_ph->Draw("ACX");
  leg_ph->Draw();

  DrawAxis(mg_ph, leg_ph, 0.3);

  C_ph->SaveAs(phRes+".eps");

#endif

  gDirectory->cd(0);

#ifdef JET

  //////////
  // Jets //
  //////////

  // BJets Reconstruction efficiency
  int bID = 5;
  std::pair<TH1D*,TH1D*> histos_btag = GetEff<Jet>(branchPFJet, branchParticle,"BTag", bID, treeReader);

  // TauJets Reconstruction efficiency
  int tauID = 15;
  std::pair<TH1D*,TH1D*> histos_tautag = GetEff<Jet>(branchPFJet, branchParticle,"TauTag", tauID, treeReader);

  // PFJets Energy Resolution
  std::vector<resolPlot> plots_pfjets;
  HistogramsCollection(&plots_pfjets, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "PFJet");
  GetJetsEres( &plots_pfjets, branchPFJet, branchGenJet, treeReader);
  TGraphErrors gr_pfjets = EresGraph(&plots_pfjets, true);
  gr_pfjets.SetName("pfJet");


  // PFJets Energy Resolution
  std::vector<resolPlot> plots_calojets;
  HistogramsCollection(&plots_calojets, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "CaloJet");
  GetJetsEres( &plots_calojets, branchCaloJet, branchGenJet, treeReader);
  TGraphErrors gr_calojets = EresGraph(&plots_calojets, true);
  gr_calojets.SetName("caloJet");

  // Canvas

  TString btagEff = "btagEff";
  TCanvas *C_btag1 = new TCanvas(btagEff,btagEff, 1000, 500);
  C_btag1->Divide(2);
  C_btag1->cd(1);
  TLegend *leg_btag1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  gPad->SetLogx();
  histos_btag.first->Draw();
  addHist(histos_btag.first, leg_btag1, 1);

  DrawAxis(histos_btag.first, leg_btag1, 0);
  leg_btag1->Draw();

  C_btag1->cd(2);
  TLegend *leg_btag2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  histos_btag.second->Draw();
  addHist(histos_btag.second, leg_btag2, 1);

  DrawAxis(histos_btag.second, leg_btag2, 1);
  delete(leg_btag2);
  C_btag1->cd(0);

  TString tautagEff = "tautagEff";
  TCanvas *C_tautag1 = new TCanvas(tautagEff,tautagEff, 1000, 500);
  C_tautag1->Divide(2);
  C_tautag1->cd(1);
  TLegend *leg_tautag1 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  gPad->SetLogx();
  histos_tautag.first->Draw();
  addHist(histos_tautag.first, leg_tautag1, 1);

  DrawAxis(histos_tautag.first, leg_tautag1, 0);
  leg_tautag1->Draw();

  C_tautag1->cd(2);
  TLegend *leg_tautag2 = new TLegend(effLegXmin,effLegYmin,effLegXmax,effLegYmax);

  histos_tautag.second->Draw();
  addHist(histos_tautag.second, leg_tautag2, 1);

  DrawAxis(histos_tautag.second, leg_tautag2, 1);
  delete(leg_tautag2);
  C_tautag1->cd(0);

  TString jetRes = "jetERes";
  TCanvas *C_jet = new TCanvas(jetRes,jetRes, 1000, 500);
  TMultiGraph *mg_jet = new TMultiGraph(jetRes,jetRes);
  TLegend *leg_jet = new TLegend(0.52,0.7,0.9,0.9);

  addGraph(mg_jet, &gr_calojets, leg_jet, 3);
  addGraph(mg_jet, &gr_pfjets, leg_jet, 1);

  mg_jet->Draw("ACX");
  leg_jet->Draw();

  DrawAxis(mg_jet, leg_jet, 0.25);

  C_btag1->SaveAs(btagEff+".eps");
  C_jet->SaveAs(jetRes+".eps");

#endif
  
  /*
  // CaloJets Energy Resolution
  std::vector<resolPlot> plots_calojets;
  HistogramsCollection(&plots_calojets, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "caloJet");
  GetJetsEres( &plots_calojets, branchCaloJet, branchGenJet, treeReader);
  TGraphErrors gr_calojets = EresGraph(&plots_calojets, true);
  gr_calojets.SetName("caloJet");
  */

  TFile *fout = new TFile(outputFile,"recreate");

#ifdef ELECTRON
  for (int bin = 0; bin < Nbins; bin++)
  {
      plots_el.at(bin).cenResolHist->Write();
      plots_eltrack.at(bin).cenResolHist->Write();
      plots_eltower.at(bin).cenResolHist->Write();
  }
  //  gr.Write();
  histos_el.first->Write();
  //histos_el.second->Write();
  histos_eltrack.first->Write();
  //histos_eltrack.second->Write();
  histos_eltower.first->Write();
  C_el1->Write();
  C_el2->Write();
#endif

#ifdef MUON
  histos_mu.first->Write();
  histos_mu.second->Write();
  histos_mutrack.first->Write();
  histos_mutrack.second->Write();
  C_mu1->Write();
  C_mu->Write();
#endif

#ifdef PHOTON
  histos_ph.first->Write();
  histos_ph.second->Write();
  C_ph1->Write();
  C_ph->Write();
#endif

#ifdef JET
  for (int bin = 0; bin < Nbins; bin++)
  {
      plots_pfjets.at(bin).cenResolHist->Write();
      plots_calojets.at(bin).cenResolHist->Write();
  }
  histos_btag.first->Write();
  histos_btag.second->Write();
  histos_tautag.first->Write();
  histos_tautag.second->Write();
  C_btag1->Write();
  C_tautag1->Write();
  C_jet->Write();
#endif

  fout->Write();

  cout << "** Exiting..." << endl;

  //delete plots;
  //delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
