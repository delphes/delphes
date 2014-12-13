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

#include <cassert>
#include <iostream>
#include <utility>
#include <algorithm>
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "external/ExRootAnalysis/ExRootConfReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "display/DelphesCaloData.h"
#include "display/DelphesBranchElement.h"
#include "display/Delphes3DGeometry.h"
#include "display/DelphesEventDisplay.h"
#include "classes/DelphesClasses.h"
#include "TEveElement.h"
#include "TEveJetCone.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveCalo.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveTrans.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"
#include "TEveArrow.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRootBrowser.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGProgressBar.h"
#include "TGNumberEntry.h"
#include "TRootEmbeddedCanvas.h"
#include "TClonesArray.h"
#include "TEveEventManager.h"
#include "TCanvas.h"
#include "TH1F.h"

DelphesEventDisplay::DelphesEventDisplay()
{
   event_id_ = 0;
   tkRadius_ = 1.29;
   totRadius_ = 2.0;
   tkHalfLength_ = 3.0;
   muHalfLength_ = 6.0;
   bz_ = 3.8;
   chain_ = new TChain("Delphes");
   treeReader_ = 0;
   delphesDisplay_ = 0;
   etaAxis_ = 0;
   phiAxis_ = 0;
}

DelphesEventDisplay::~DelphesEventDisplay()
{
   delete chain_;
}

void DelphesEventDisplay::EventChanged(Int_t e)
{
  if( e!= event_id_ ) {
    event_id_ = e;
    Emit("EventChanged(Int_t)",e);
    load_event();
  }
}


DelphesEventDisplay::DelphesEventDisplay(const char *configFile, const char *inputFile, Delphes3DGeometry& det3D)
{
   event_id_ = 0;
   tkRadius_ = 1.29;
   totRadius_ = 2.0;
   tkHalfLength_ = 3.0;
   bz_ = 3.8;
   chain_ = new TChain("Delphes");
   treeReader_ = 0;
   delphesDisplay_ = 0;

   // initialize the application
   TEveManager::Create(kTRUE, "IV");
   fStatusBar_ = gEve->GetBrowser()->GetStatusBar();
   TGeoManager* geom = gGeoManager;

   // build the detector
   tkRadius_ = det3D.getTrackerRadius();
   totRadius_ = det3D.getDetectorRadius();
   tkHalfLength_ = det3D.getTrackerHalfLength();
   muHalfLength_ = det3D.getDetectorHalfLength();
   bz_ = det3D.getBField();
   etaAxis_ = det3D.getCaloAxes().first;
   phiAxis_ = det3D.getCaloAxes().second;
   TGeoVolume* top = det3D.getDetector(false);
   geom->SetTopVolume(top);
   TEveElementList *geometry = new TEveElementList("Geometry");
   TObjArray* nodes = top->GetNodes();
   TIter itNodes(nodes);
   TGeoNode* nodeobj;
   TEveGeoTopNode* node;
   while((nodeobj = (TGeoNode*)itNodes.Next())) {
     node = new TEveGeoTopNode(gGeoManager,nodeobj);
     node->UseNodeTrans();
     geometry->AddElement(node);
   }

   // Create chain of root trees
   chain_->Add(inputFile);

   // Create object of class ExRootTreeReader
   fStatusBar_->SetText("Opening Delphes data file", 1);
   gSystem->ProcessEvents();
   treeReader_ = new ExRootTreeReader(chain_);

   // prepare data collections
   readConfig(configFile, elements_);
   for(std::vector<DelphesBranchBase*>::iterator element = elements_.begin(); element<elements_.end(); ++element) {
     DelphesBranchElement<TEveTrackList>*   item_v1 = dynamic_cast<DelphesBranchElement<TEveTrackList>*>(*element);
     DelphesBranchElement<TEveElementList>* item_v2 = dynamic_cast<DelphesBranchElement<TEveElementList>*>(*element);
     if(item_v1) gEve->AddElement(item_v1->GetContainer());
     if(item_v2) gEve->AddElement(item_v2->GetContainer());
   }

   // viewers and scenes
   delphesDisplay_ = new DelphesDisplay;
   gEve->AddGlobalElement(geometry);
   delphesDisplay_->ImportGeomRPhi(geometry);
   delphesDisplay_->ImportGeomRhoZ(geometry);
   // find the first calo data and use that to initialize the calo display
   for(std::vector<DelphesBranchBase*>::iterator data=elements_.begin();data<elements_.end();++data) {
     if(TString((*data)->GetType())=="Tower") { // we could also use GetClassName()=="DelphesCaloData"
       DelphesCaloData* container = dynamic_cast<DelphesBranchElement<DelphesCaloData>*>((*data))->GetContainer();
       assert(container);
       TEveCalo3D *calo3d = new TEveCalo3D(container);
       calo3d->SetBarrelRadius(tkRadius_);
       calo3d->SetEndCapPos(tkHalfLength_);
       gEve->AddGlobalElement(calo3d);
       delphesDisplay_->ImportCaloRPhi(calo3d);
       delphesDisplay_->ImportCaloRhoZ(calo3d);
       TEveCaloLego *lego = new TEveCaloLego(container);
       lego->InitMainTrans();
 //    lego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
       lego->RefMainTrans().SetScale(100, 100, TMath::Pi());
       lego->SetAutoRebin(kFALSE);
       lego->Set2DMode(TEveCaloLego::kValSizeOutline);
       delphesDisplay_->ImportCaloLego(lego);
       break;
     }
   }

   // the GUI: control panel, summary tab
   make_gui();

   //ready...
   fStatusBar_->SetText("Ready.", 1);
   gSystem->ProcessEvents();
   load_event();
   gEve->Redraw3D(kTRUE);   

}

// function that parses the config to extract the branches of interest and prepare containers
void DelphesEventDisplay::readConfig(const char *configFile, std::vector<DelphesBranchBase*>& elements) {
   ExRootConfReader *confReader = new ExRootConfReader;
   confReader->ReadFile(configFile);
   ExRootConfParam branches = confReader->GetParam("TreeWriter::Branch");
   Int_t nBranches = branches.GetSize()/3;
   DelphesBranchElement<TEveTrackList>* tlist;
   DelphesBranchElement<DelphesCaloData>* clist;
   DelphesBranchElement<TEveElementList>* elist;
   // first loop with all but tracks
   for(Int_t b = 0; b<nBranches; ++b) {
     TString input = branches[b*3].GetString();
     TString name = branches[b*3+1].GetString();
     TString className = branches[b*3+2].GetString();
     if(className=="Tower") {
       if(input.Contains("eflow",TString::kIgnoreCase) || name.Contains("eflow",TString::kIgnoreCase)) continue; //no eflow
       clist = new DelphesBranchElement<DelphesCaloData>(name,treeReader_->UseBranch(name),kBlack);
       clist->GetContainer()->SetEtaBins(etaAxis_);
       clist->GetContainer()->SetPhiBins(phiAxis_);
       elements.push_back(clist);
     } else if(className=="Jet") {
       if(input.Contains("GenJetFinder")) {
         elist = new DelphesBranchElement<TEveElementList>(name,treeReader_->UseBranch(name),kCyan);
         elist->GetContainer()->SetRnrSelf(false);
         elist->GetContainer()->SetRnrChildren(false);
         elist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
         elements.push_back(elist);
       } else {
         elist = new DelphesBranchElement<TEveElementList>(name,treeReader_->UseBranch(name),kYellow);
         elist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
         elements.push_back(elist);
       }
     } else if(className=="Electron") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,treeReader_->UseBranch(name),kRed);
       tlist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
       elements.push_back(tlist);
     } else if(className=="Photon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,treeReader_->UseBranch(name),kYellow);
       tlist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
       elements.push_back(tlist);
     } else if(className=="Muon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,treeReader_->UseBranch(name),kGreen);
       tlist->SetTrackingVolume(totRadius_, muHalfLength_, bz_);
       elements.push_back(tlist);
     } else if(className=="MissingET") {
       elist = new DelphesBranchElement<TEveElementList>(name,treeReader_->UseBranch(name),kViolet);
       elist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
       elements.push_back(elist);
     } else if(className=="GenParticle") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,treeReader_->UseBranch(name),kCyan);
       tlist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
       tlist->GetContainer()->SetRnrSelf(false);
       tlist->GetContainer()->SetRnrChildren(false);
       elements.push_back(tlist);
     } else {
       continue;
     }
   }
   // second loop for tracks
   for(Int_t b = 0; b<nBranches; ++b) {
     TString input = branches[b*3].GetString();
     TString name = branches[b*3+1].GetString();
     TString className = branches[b*3+2].GetString();
     if(className=="Track") {
       if(input.Contains("eflow",TString::kIgnoreCase) || name.Contains("eflow",TString::kIgnoreCase)) continue; //no eflow
       tlist = new DelphesBranchElement<TEveTrackList>(name,treeReader_->UseBranch(name),kBlue);
       tlist->SetTrackingVolume(tkRadius_, tkHalfLength_, bz_);
       elements.push_back(tlist);
     }
   }
}

void DelphesEventDisplay::load_event()
{
   // Load event specified in global event_id_.
   // The contents of previous event are removed.

   // safety
   if(event_id_ >= treeReader_->GetEntries() || event_id_<0 ) return;

   // message
   fStatusBar_->SetText(Form("Loading event %d.", event_id_), 1);
   gSystem->ProcessEvents();

   // clear the previous event
   gEve->GetViewers()->DeleteAnnotations();
   for(std::vector<DelphesBranchBase*>::iterator data=elements_.begin();data<elements_.end();++data) {
     (*data)->Reset();
   }

   // Load selected branches with data from specified event
   treeReader_->ReadEntry(event_id_);
   for(std::vector<DelphesBranchBase*>::iterator data=elements_.begin();data<elements_.end();++data) {
     (*data)->ReadBranch();
   }

   // update display
   TEveElement* top = (TEveElement*)gEve->GetCurrentEvent();
   delphesDisplay_->DestroyEventRPhi();
   delphesDisplay_->ImportEventRPhi(top);
   delphesDisplay_->DestroyEventRhoZ();
   delphesDisplay_->ImportEventRhoZ(top);
   update_html_summary();
   plotSummary_->FillEvent();
   plotSummary_->Draw();

   gEve->Redraw3D(kFALSE, kTRUE);
   fStatusBar_->SetText(Form("Loaded event %d.", event_id_), 1);
   gSystem->ProcessEvents();
}

void DelphesEventDisplay::update_html_summary()
{
   // Update summary of current event.

   TEveElement::List_i i;
   TEveElement::List_i j;
   Int_t k;
   TEveElement *el;
   DelphesHtmlObjTable *table;
   TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
   if (mgr) {
      htmlSummary_->Clear("D");
      for (i=mgr->BeginChildren(); i!=mgr->EndChildren(); ++i) {
         el = ((TEveElement*)(*i));
         if (el->IsA() == TEvePointSet::Class()) {
            TEvePointSet *ps = (TEvePointSet *)el;
            TString ename  = ps->GetElementName();
            TString etitle = ps->GetElementTitle();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            etitle.Remove(0, 2);
            Int_t nel = atoi(etitle.Data());
            table = htmlSummary_->AddTable(ename, 0, nel);
         }
         else if (el->IsA() == TEveTrackList::Class()) {
            TEveTrackList *tracks = (TEveTrackList *)el;
            TString ename  = tracks->GetElementName();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            table = htmlSummary_->AddTable(ename.Data(), 5, 
                     tracks->NumChildren(), kTRUE, "first");
            table->SetLabel(0, "Momentum");
            table->SetLabel(1, "P_t");
            table->SetLabel(2, "Phi");
            table->SetLabel(3, "Theta");
            table->SetLabel(4, "Eta");
            k=0;
            for (j=tracks->BeginChildren(); j!=tracks->EndChildren(); ++j) {
               Float_t p     = ((TEveTrack*)(*j))->GetMomentum().Mag();
               table->SetValue(0, k, p);
               Float_t pt    = ((TEveTrack*)(*j))->GetMomentum().Perp();
               table->SetValue(1, k, pt);
               Float_t phi   = ((TEveTrack*)(*j))->GetMomentum().Phi();
               table->SetValue(2, k, phi);
               Float_t theta = ((TEveTrack*)(*j))->GetMomentum().Theta();
               table->SetValue(3, k, theta);
               Float_t eta   = theta>0.0005 && theta<3.1413 ? ((TEveTrack*)(*j))->GetMomentum().Eta() : 1e10;
               table->SetValue(4, k, eta);
               ++k;
            }
         }
      }
      htmlSummary_->Build();
      gHtml_->Clear();
      gHtml_->ParseText((char*)htmlSummary_->Html().Data());
      gHtml_->Layout();
   }
  
}

/******************************************************************************/
// GUI
/******************************************************************************/

void DelphesEventDisplay::make_gui()
{
   // Create minimal GUI for event navigation.

   // add a tab on the left
   TEveBrowser* browser = gEve->GetBrowser();
   browser->SetWindowName("Delphes Event Display");
   browser->StartEmbedding(TRootBrowser::kLeft);

   // set the main title
   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("Delphes Event Display");
   frmMain->SetCleanup(kDeepCleanup);

   // build the navigation menu
   TString icondir;
   if(gSystem->Getenv("ROOTSYS"))
     icondir = Form("%s/icons/", gSystem->Getenv("ROOTSYS"));
   if(!gSystem->OpenDirectory(icondir)) 
     icondir = Form("%s/icons/", (const char*)gSystem->GetFromPipe("root-config --etcdir") );
   TGGroupFrame* vf = new TGGroupFrame(frmMain,"Event navigation",kVerticalFrame | kFitWidth );
   {
     TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
     {
        TGPictureButton* b = 0;

        b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
        hf->AddFrame(b, new TGLayoutHints(kLHintsLeft | kLHintsCenterY , 10, 2, 10, 10));
        b->Connect("Clicked()", "DelphesEventDisplay", this, "Bck()");

        TGNumberEntry* numberEntry = new TGNumberEntry(hf,0,9,-1,TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, treeReader_->GetEntries());
        hf->AddFrame(numberEntry, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY , 2, 0, 10, 10));
        this->Connect("EventChanged(Int_t)","TGNumberEntry",numberEntry,"SetIntNumber(Long_t)");
        numberEntry->GetNumberEntry()->Connect("TextChanged(char*)", "DelphesEventDisplay", this, "PreSetEv(char*)");
        numberEntry->GetNumberEntry()->Connect("ReturnPressed()", "DelphesEventDisplay", this, "GoTo()");
        
        b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
        hf->AddFrame(b, new TGLayoutHints(kLHintsRight | kLHintsCenterY , 2, 10, 10, 10));
        b->Connect("Clicked()", "DelphesEventDisplay", this, "Fwd()");
  
     }
     vf->AddFrame(hf, new TGLayoutHints(kLHintsExpandX , 2, 2, 2, 2));

     TGHProgressBar* progress = new TGHProgressBar(frmMain, TGProgressBar::kFancy, 100);
     progress->SetMax( treeReader_->GetEntries());
     progress->ShowPosition(kTRUE, kFALSE, "Event %.0f");
     progress->SetBarColor("green");
     vf->AddFrame(progress, new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 5));
     this->Connect("EventChanged(Int_t)","TGHProgressBar",progress,"SetPosition(Float_t)");
   }
   frmMain->AddFrame(vf, new TGLayoutHints(kLHintsExpandX , 5, 5, 5, 5));
   vf = new TGGroupFrame(frmMain,"Batch operations",kVerticalFrame | kFitWidth );
   {
        TGTextButton *b = new TGTextButton(vf, "Initialize Summary Plots");
        vf->AddFrame(b, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX, 10, 10, 10, 10));
        b->Connect("Clicked()", "DelphesEventDisplay", this, "InitSummaryPlots()");
   }
   frmMain->AddFrame(vf, new TGLayoutHints(kLHintsExpandX , 5, 5, 5, 5));

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();
   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);

   // the summary tab
   htmlSummary_ = new DelphesHtmlSummary("Delphes Event Display Summary Table");
   TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
   gHtml_ = new TGHtml(0, 100, 100);
   TEveWindowFrame *wf = slot->MakeFrame(gHtml_);
   gHtml_->MapSubwindows();
   wf->SetElementName("Summary tables");

   // plot tab
   slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
   TEveWindowTab* tab = slot->MakeTab();
   tab->SetElementName("Summary plots");
   tab->SetShowTitleBar(kFALSE);
   plotSummary_ = new DelphesPlotSummary(tab);
   plotSummary_->Init(elements_);
   plotSummary_->Connect("Progress(Int_t)", "DelphesEventDisplay", this, "DisplayProgress(Int_t)");

}


