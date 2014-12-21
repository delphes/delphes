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

#include "display/DelphesPlotSummary.h"
#include "TRootEmbeddedCanvas.h"
#include <algorithm>

bool vecsorter (TLorentzVector i,TLorentzVector j) { return (i.Pt()>j.Pt()); }

DelphesPlotSummary::DelphesPlotSummary(TEveWindowTab* tab):tab_(tab) {}

DelphesPlotSummary::~DelphesPlotSummary() {}

void DelphesPlotSummary::Progress(Int_t p)
{
  Emit("Progress(Int_t)",p);
}

void DelphesPlotSummary::Init(std::vector<DelphesBranchBase*>& elements) {
  elements_ = &elements;
  // loop on the elements, and create tabs
  for(std::vector<DelphesBranchBase*>::iterator data=elements.begin();data<elements.end();++data) {
    // the canvas
    TEveWindowSlot* slot = tab_->NewSlot();
    TRootEmbeddedCanvas* trec = new TRootEmbeddedCanvas();
    TCanvas* canvas = trec->GetCanvas();
    TEveWindowFrame * wf = slot->MakeFrame(trec);
    wf->SetElementName((*data)->GetName());
    canvas->Divide(3,3);
    canvases_[(*data)->GetName()] = canvas;
    // the histograms
    TH1F* h;
    std::vector<TH1F*> histograms;
    histograms.reserve(9);
    h = new TH1F(Form("%sPt",(*data)->GetName()),Form("%s Pt",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("%sEta",(*data)->GetName()),Form("%s Eta",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("%sPhi",(*data)->GetName()),Form("%s Phi",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("l%sPt",(*data)->GetName()),Form("leading %s Pt",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("l%sEta",(*data)->GetName()),Form("leading %s Eta",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("l%sPhi",(*data)->GetName()),Form("leading %s Phi",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("sl%sPt",(*data)->GetName()),Form("subleading %s Pt",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("sl%sEta",(*data)->GetName()),Form("subleading %s Eta",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    h = new TH1F(Form("sl%sPhi",(*data)->GetName()),Form("subleading %s Phi",(*data)->GetName()),100,0,-1);
    histograms.push_back(h);
    histograms_[(*data)->GetName()] = histograms;
    // the event histograms
    TH1F* h1 = (TH1F*)histograms[0]->Clone(); h1->Reset(); h1->SetLineColor(kBlue);
    TH1F* h2 = (TH1F*)histograms[1]->Clone(); h2->Reset(); h2->SetLineColor(kBlue);
    TH1F* h3 = (TH1F*)histograms[2]->Clone(); h3->Reset(); h3->SetLineColor(kBlue);
    std::vector<TH1F*> hv;
    hv.push_back(h1);
    hv.push_back(h2);
    hv.push_back(h3);
    eventProfiles_[(*data)->GetName()] = hv;
    // the event markers
    TMarker *m;
    std::vector<TMarker*> mv;
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    m = new TMarker(0,0,29); m->SetMarkerColor(kBlue); m->SetMarkerSize(3);
    mv.push_back(m);
    eventMarkers_[(*data)->GetName()] = mv;
  }
}

void DelphesPlotSummary::FillSample(ExRootTreeReader* treeReader, Int_t event_id) {
  Int_t entries = treeReader->GetEntries();
  for(Int_t i=0;i<entries;++i) {
    treeReader->ReadEntry(i);
    for(std::vector<DelphesBranchBase*>::iterator element = elements_->begin();element<elements_->end();++element) {
      std::vector<TLorentzVector> vectors = (*element)->GetVectors();
      std::sort(vectors.begin(), vectors.end(), vecsorter); 
      std::vector<TH1F*> histograms = histograms_[(*element)->GetName()];
      for(std::vector<TLorentzVector>::iterator it=vectors.begin(); it<vectors.end();++it) {
        histograms[0]->Fill(it->Pt());
        histograms[1]->Fill(it->Eta());
        histograms[2]->Fill(it->Phi());
        if(it==vectors.begin()) {
          histograms[3]->Fill(it->Pt());
          histograms[4]->Fill(it->Eta());
          histograms[5]->Fill(it->Phi());
        }
        if(it==vectors.begin()+1) {
          histograms[6]->Fill(it->Pt());
          histograms[7]->Fill(it->Eta());
          histograms[8]->Fill(it->Phi());
        }
      }
    }
    Progress(int(100*i/entries));
  }
  treeReader->ReadEntry(event_id);
  Progress(100);
}

void DelphesPlotSummary::Draw() {
  for(std::map< TString, TCanvas* >::iterator it=canvases_.begin(); it!=canvases_.end(); ++it) {
    TCanvas* c = it->second;
    std::vector<TH1F*> histograms = histograms_[it->first];
    std::vector<TH1F*> eventProfiles = eventProfiles_[it->first];
    std::vector<TMarker*> eventMarkers = eventMarkers_[it->first];
    for(Int_t i=0;i<9;++i) {
      c->cd(i+1);
      if(histograms[i]->GetEntries()==0) continue;
      histograms[i]->Draw();
      if(i<3) {
        eventProfiles[i]->Draw("same");
      } else {
        eventMarkers[i-3]->Draw("same");
      }
    }
    c->Update();
  } 
}

void DelphesPlotSummary::FillEvent() {
  // clear event histograms and markers
  for(std::map< TString, std::vector<TH1F*> >::iterator hv = eventProfiles_.begin(); hv!=eventProfiles_.end();++hv) {
    for(std::vector<TH1F*>::iterator h = hv->second.begin(); h<hv->second.end();++h) {
      (*h)->Reset();
    }
  }
  for(std::map< TString, std::vector<TMarker*> >::iterator mv = eventMarkers_.begin(); mv!=eventMarkers_.end();++mv) {
    for(std::vector<TMarker*>::iterator m = mv->second.begin(); m<mv->second.end();++m) {
      (*m)->SetMarkerSize(0);
    }
  }
  // loop over the elements and fill markers with event data
  for(std::vector<DelphesBranchBase*>::iterator element = elements_->begin();element<elements_->end();++element) {
    std::vector<TLorentzVector> vectors = (*element)->GetVectors();
    std::sort(vectors.begin(), vectors.end(), vecsorter); 
    std::vector<TH1F*> hv = eventProfiles_[(*element)->GetName()];
    TH1F* h1 = hv[0]; h1->Reset();
    TH1F* h2 = hv[1]; h1->Reset();
    TH1F* h3 = hv[2]; h1->Reset();
    std::vector<TMarker*> mv = eventMarkers_[(*element)->GetName()];
    for(std::vector<TLorentzVector>::iterator it=vectors.begin(); it<vectors.end();++it) {
      h1->Fill(it->Pt());
      h2->Fill(it->Eta());
      h3->Fill(it->Phi());
      if(it==vectors.begin()) {
        mv[0]->SetX(it->Pt()); mv[0]->SetMarkerSize(3);
        mv[1]->SetX(it->Eta()); mv[1]->SetMarkerSize(3);
        mv[2]->SetX(it->Phi()); mv[2]->SetMarkerSize(3);
      }
      if(it==vectors.begin()+1) {
        mv[3]->SetX(it->Pt()); mv[3]->SetMarkerSize(3);
        mv[4]->SetX(it->Eta()); mv[4]->SetMarkerSize(3);
        mv[5]->SetX(it->Phi()); mv[5]->SetMarkerSize(3);
      }
    }
  }
}
