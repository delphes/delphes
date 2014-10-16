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
#include "TClonesArray.h"

DelphesEventDisplay::DelphesEventDisplay()
{
   event_id_ = 0;
   tkRadius_ = 1.29;
   totRadius_ = 2.0;
   tkHalfLength_ = 3.0;
   bz_ = 3.8;
   chain_ = new TChain("Delphes");
   treeReader_ = 0;
   delphesDisplay_ = 0;
}

DelphesEventDisplay::~DelphesEventDisplay()
{
   delete chain_;
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
   TGeoManager* geom = gGeoManager;

   // build the detector
   tkRadius_ = det3D.getTrackerRadius();
   totRadius_ = det3D.getDetectorRadius();
   tkHalfLength_ = det3D.getTrackerHalfLength();
   bz_ = det3D.getBField();
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
   printf("*** Opening Delphes data file ***\n");
   treeReader_ = new ExRootTreeReader(chain_);

   // prepare data collections
   readConfig(configFile, det3D, elements_, arrays_);
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
     if(TString((*data)->GetType())=="tower") { // we could also use GetClassName()=="DelphesCaloData"
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
       lego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
       lego->SetAutoRebin(kFALSE);
       lego->Set2DMode(TEveCaloLego::kValSizeOutline);
       delphesDisplay_->ImportCaloLego(lego);
       break;
     }
   }

   make_gui();
   load_event();
   gEve->Redraw3D(kTRUE);   

}

// function that parses the config to extract the branches of interest and prepare containers
void DelphesEventDisplay::readConfig(const char *configFile, Delphes3DGeometry& det3D, std::vector<DelphesBranchBase*>& elements, std::vector<TClonesArray*>& arrays) {
   ExRootConfReader *confReader = new ExRootConfReader;
   confReader->ReadFile(configFile);
   Double_t tk_radius = det3D.getTrackerRadius();
   Double_t tk_length = det3D.getTrackerHalfLength();
   Double_t tk_Bz     = det3D.getBField();
   Double_t mu_radius = det3D.getDetectorRadius();
   Double_t mu_length = det3D.getDetectorHalfLength();
   TAxis*   etaAxis   = det3D.getCaloAxes().first;
   TAxis*   phiAxis   = det3D.getCaloAxes().second;
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
       clist = new DelphesBranchElement<DelphesCaloData>(name,"tower",kBlack);
       clist->GetContainer()->SetEtaBins(etaAxis);
       clist->GetContainer()->SetPhiBins(phiAxis);
       elements.push_back(clist);
     } else if(className=="Jet") {
       if(input.Contains("GenJetFinder")) {
         elist = new DelphesBranchElement<TEveElementList>(name,"jet",kCyan);
         elist->GetContainer()->SetRnrSelf(false);
         elist->GetContainer()->SetRnrChildren(false);
         elements.push_back(elist);
       } else {
         elements.push_back(new DelphesBranchElement<TEveElementList>(name,"jet",kYellow));
       }
     } else if(className=="Electron") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"electron",kRed);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else if(className=="Photon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"photon",kYellow);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., 0.);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else if(className=="Muon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"muon",kGreen);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(mu_radius);
       trkProp->SetMaxZ(mu_length);
     } else if(className=="MissingET") {
       elements.push_back(new DelphesBranchElement<TEveElementList>(name,"vector",kViolet));
     } else if(className=="GenParticle") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"genparticle",kCyan);
       elements.push_back(tlist);
       tlist->GetContainer()->SetRnrSelf(false);
       tlist->GetContainer()->SetRnrChildren(false);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else {
       continue;
     }
     arrays.push_back(treeReader_->UseBranch(name));
   }
   // second loop for tracks
   for(Int_t b = 0; b<nBranches; ++b) {
     TString input = branches[b*3].GetString();
     TString name = branches[b*3+1].GetString();
     TString className = branches[b*3+2].GetString();
     if(className=="Track") {
       if(input.Contains("eflow",TString::kIgnoreCase) || name.Contains("eflow",TString::kIgnoreCase)) continue; //no eflow
       tlist = new DelphesBranchElement<TEveTrackList>(name,"track",kBlue);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
       arrays.push_back(treeReader_->UseBranch(name));
     }
   }
}

void DelphesEventDisplay::load_event()
{
   // Load event specified in global event_id_.
   // The contents of previous event are removed.

   // safety
   if(event_id_ >= treeReader_->GetEntries() || event_id_<0 ) return;

   //TODO move this to the status bar ???
   printf("Loading event %d.\n", event_id_);

   // clear the previous event
   gEve->GetViewers()->DeleteAnnotations();
   for(std::vector<DelphesBranchBase*>::iterator data=elements_.begin();data<elements_.end();++data) {
     (*data)->Reset();
   }

   // Load selected branches with data from specified event
   treeReader_->ReadEntry(event_id_);

   // loop over selected branches, and apply the proper recipe to fill the collections.
   // this is basically to loop on arrays_ to fill elements_.
   std::vector<TClonesArray*>::iterator data = arrays_.begin();
   std::vector<DelphesBranchBase*>::iterator element = elements_.begin();
   for(; data<arrays_.end() && element<elements_.end(); ++data, ++element) {
     TString type = (*element)->GetType();
     // branch on the element type
     if(type=="tower") delphes_read_towers(*data,*element);
     else if(type=="track" || type=="photon" || type=="electron" || type=="muon" || type=="genparticle") delphes_read_tracks(*data,*element);
     else if(type=="jet") delphes_read_jets(*data,*element);
     else if(type=="vector") delphes_read_vectors(*data,*element);
   }

   // update display
   TEveElement* top = (TEveElement*)gEve->GetCurrentEvent();
   delphesDisplay_->DestroyEventRPhi();
   delphesDisplay_->ImportEventRPhi(top);
   delphesDisplay_->DestroyEventRhoZ();
   delphesDisplay_->ImportEventRhoZ(top);
   //update_html_summary();

   gEve->Redraw3D(kFALSE, kTRUE);
}

void DelphesEventDisplay::delphes_read_towers(TClonesArray* data, DelphesBranchBase* element) {
  DelphesCaloData* container = dynamic_cast<DelphesBranchElement<DelphesCaloData>*>(element)->GetContainer();
  assert(container);
  // Loop over all towers
  TIter itTower(data);
  Tower *tower;
  while((tower = (Tower *) itTower.Next()))
  {
    container->AddTower(tower->Edges[0], tower->Edges[1], tower->Edges[2], tower->Edges[3]);
    container->FillSlice(0, tower->Eem);
    container->FillSlice(1, tower->Ehad);
  }
  container->DataChanged();
}

void DelphesEventDisplay::delphes_read_tracks(TClonesArray* data, DelphesBranchBase* element) {
  TEveTrackList* container = dynamic_cast<DelphesBranchElement<TEveTrackList>*>(element)->GetContainer();
  assert(container);
  TString type = element->GetType();
  TIter itTrack(data);
  Int_t counter = 0;
  TEveTrack *eveTrack;
  TEveTrackPropagator *trkProp = container->GetPropagator();
  if(type=="track") {
    // Loop over all tracks
    Track *track;
    while((track = (Track *) itTrack.Next())) {
      TParticle pb(track->PID, 1, 0, 0, 0, 0,
                   track->P4().Px(), track->P4().Py(),
                   track->P4().Pz(), track->P4().E(),
                   track->X, track->Y, track->Z, 0.0);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="electron") {
    // Loop over all electrons
    Electron *electron;
    while((electron = (Electron *) itTrack.Next())) {
      TParticle pb(electron->Charge<0?11:-11, 1, 0, 0, 0, 0,
                   electron->P4().Px(), electron->P4().Py(),
                   electron->P4().Pz(), electron->P4().E(),
                   0., 0., 0., 0.);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="muon") {
    // Loop over all muons
    Muon *muon;
    while((muon = (Muon *) itTrack.Next())) {
      TParticle pb(muon->Charge<0?13:-13, 1, 0, 0, 0, 0,
                   muon->P4().Px(), muon->P4().Py(),
                   muon->P4().Pz(), muon->P4().E(),
                   0., 0., 0., 0.);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="photon") {
    // Loop over all photons
    Photon *photon;
    while((photon = (Photon *) itTrack.Next())) {
      TParticle pb(22, 1, 0, 0, 0, 0,
                   photon->P4().Px(), photon->P4().Py(),
                   photon->P4().Pz(), photon->P4().E(),
                   0., 0., 0., 0.);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      eveTrack->SetLineStyle(7);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="genparticle") {
    // Loop over all particles
    GenParticle *particle;
    while((particle = (GenParticle *) itTrack.Next())) {
      TParticle pb(particle->PID, particle->Status, particle->M1, particle->M2, particle->D1, particle->D2,
                   particle->P4().Px(), particle->P4().Py(),
                   particle->P4().Pz(), particle->P4().E(),
                   particle->X, particle->Y, particle->Z, particle->T);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      if(particle->Charge==0) eveTrack->SetLineStyle(7);
      eveTrack->MakeTrack();
    }
  }
}

void DelphesEventDisplay::delphes_read_jets(TClonesArray* data, DelphesBranchBase* element) {
  TEveElementList* container = dynamic_cast<DelphesBranchElement<TEveElementList>*>(element)->GetContainer();
  assert(container);
  TIter itJet(data);
  Jet *jet;
  TEveJetCone *eveJetCone;
  // Loop over all jets
  Int_t counter = 0;
  while((jet = (Jet *) itJet.Next()))
  {
    eveJetCone = new TEveJetCone();
    eveJetCone->SetTitle(Form("jet [%d]: Pt=%f, Eta=%f, \nPhi=%f, M=%f",counter,jet->PT, jet->Eta, jet->Phi, jet->Mass));
    eveJetCone->SetName(Form("jet [%d]", counter++));
    eveJetCone->SetMainTransparency(60);
    eveJetCone->SetLineColor(element->GetColor());
    eveJetCone->SetFillColor(element->GetColor());
    eveJetCone->SetCylinder(tkRadius_ - 10, tkHalfLength_ - 10);
    eveJetCone->SetPickable(kTRUE);
    eveJetCone->AddEllipticCone(jet->Eta, jet->Phi, jet->DeltaEta, jet->DeltaPhi);
    container->AddElement(eveJetCone);
  }
}

void DelphesEventDisplay::delphes_read_vectors(TClonesArray* data, DelphesBranchBase* element) {
  TEveElementList* container = dynamic_cast<DelphesBranchElement<TEveElementList>*>(element)->GetContainer();
  assert(container);
  TIter itMet(data);
  MissingET *MET;
  TEveArrow *eveMet;
  // Missing Et
  Double_t maxPt = 50.;
  // TODO to be changed as we don't have access to maxPt anymore. MET scale could be a general parameter set in GUI
  while((MET = (MissingET*) itMet.Next())) {
    eveMet = new TEveArrow((tkRadius_ * MET->MET/maxPt)*cos(MET->Phi), (tkRadius_ * MET->MET/maxPt)*sin(MET->Phi), 0., 0., 0., 0.);
    eveMet->SetMainColor(element->GetColor());
    eveMet->SetTubeR(0.04);
    eveMet->SetConeR(0.08);
    eveMet->SetConeL(0.10);
    eveMet->SetPickable(kTRUE);
    eveMet->SetName("Missing Et");
    eveMet->SetTitle(Form("Missing Et (%.1f GeV)",MET->MET));
    container->AddElement(eveMet);
  }
}

/******************************************************************************/
// GUI
/******************************************************************************/

void DelphesEventDisplay::make_gui()
{
   // Create minimal GUI for event navigation.
   // TODO: better GUI could be made based on the ch15 of the manual (Writing a GUI)

   // add a tab on the left
   TEveBrowser* browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   // set the main title
   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("Delphes Event Display");
   frmMain->SetCleanup(kDeepCleanup);

   // build the navigation menu
   TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
   {
      TString icondir;
      if(gSystem->Getenv("ROOTSYS"))
        icondir = Form("%s/icons/", gSystem->Getenv("ROOTSYS"));
      if(!gSystem->OpenDirectory(icondir)) 
        icondir = Form("%s/icons/", (const char*)gSystem->GetFromPipe("root-config --etcdir") );
      TGPictureButton* b = 0;

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "DelphesEventDisplay", this, "Bck()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "DelphesEventDisplay", this, "Fwd()");
   }
   frmMain->AddFrame(hf);
   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();
   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);
}

