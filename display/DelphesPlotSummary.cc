#include "display/DelphesPlotSummary.h"
#include "TRootEmbeddedCanvas.h"


DelphesPlotSummary::DelphesPlotSummary(TEveWindowTab* tab):tab_(tab) {}

DelphesPlotSummary::~DelphesPlotSummary() {}

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
  }
}

void DelphesPlotSummary::FillSample(ExRootTreeReader* treeReader, Int_t event_id) {
  for(Int_t i=0;i<treeReader->GetEntries();++i) {
    treeReader->ReadEntry(i);
    for(std::vector<DelphesBranchBase*>::iterator element = elements_->begin();element<elements_->end();++element) {
      std::vector<TLorentzVector> vectors = (*element)->GetVectors();
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
  }
  treeReader->ReadEntry(event_id);
}

void DelphesPlotSummary::Draw() {
  for(std::map< TString, TCanvas* >::iterator it=canvases_.begin(); it!=canvases_.end(); ++it) {
    TCanvas* c = it->second;
    std::vector<TH1F*> histograms = histograms_[it->first];
    std::vector<TH1F*> eventProfiles = eventProfiles_[it->first];
    std::vector<TMarker*> eventMarkers = eventMarkers_[it->first];
    for(Int_t i=0;i<9;++i) {
      c->cd(i+1);
      histograms[i]->Draw();
      if(i<3) 
        eventProfiles[i]->Draw("same");
      else
        eventMarkers[i-3]->Draw("same");
    }
  } 
}

void DelphesPlotSummary::FillEvent() {
  // clear previous markers
  for(std::map< TString, std::vector<TMarker*> >::iterator mv = eventMarkers_.begin(); mv!=eventMarkers_.end();++mv) {
    for(std::vector<TMarker*>::iterator m = mv->second.begin(); m<mv->second.end();++m) {
      delete *m;
    }
  }
  eventMarkers_.clear();
  for(std::map< TString, std::vector<TH1F*> >::iterator hv = eventProfiles_.begin(); hv!=eventProfiles_.end();++hv) {
    for(std::vector<TH1F*>::iterator h = hv->second.begin(); h<hv->second.end();++h) {
      delete *h;
    }
  }
  eventProfiles_.clear();
  // loop over the elements and fill markers with event data
  TMarker *m;
  std::vector<TMarker*> mv;
  for(std::vector<DelphesBranchBase*>::iterator element = elements_->begin();element<elements_->end();++element) {
    std::vector<TLorentzVector> vectors = (*element)->GetVectors();
    std::vector<TH1F*> histograms = histograms_[(*element)->GetName()];
    TH1F* h1 = (TH1F*)histograms[0]->Clone(); h1->Reset(); h1->SetLineColor(kBlue);
    TH1F* h2 = (TH1F*)histograms[1]->Clone(); h2->Reset(); h2->SetLineColor(kBlue);
    TH1F* h3 = (TH1F*)histograms[2]->Clone(); h3->Reset(); h3->SetLineColor(kBlue);
    for(std::vector<TLorentzVector>::iterator it=vectors.begin(); it<vectors.end();++it) {
      h1->Fill(it->Pt());
      h2->Fill(it->Eta());
      h3->Fill(it->Phi());
      if(it==vectors.begin()) {
        m = new TMarker(it->Pt(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
        m = new TMarker(it->Eta(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
        m = new TMarker(it->Phi(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
      }
      if(it==vectors.begin()+1) {
        m = new TMarker(it->Pt(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
        m = new TMarker(it->Eta(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
        m = new TMarker(it->Phi(),0,29); m->SetMarkerColor(kBlue);
        mv.push_back(m);
      }
    }
    std::vector<TH1F*> hv;
    hv.push_back(h1);
    hv.push_back(h2);
    hv.push_back(h3);
    eventProfiles_[(*element)->GetName()] = hv;
    eventMarkers_[(*element)->GetName()] = mv;
  }
}
