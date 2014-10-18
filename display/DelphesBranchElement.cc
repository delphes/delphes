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

#include "display/DelphesBranchElement.h"
#include "classes/DelphesClasses.h"
#include "TEveJetCone.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveArrow.h"
#include "TEveVector.h"
#include <iostream>

// special case for calo towers
template<> DelphesBranchElement<DelphesCaloData>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color):DelphesBranchBase(name, branch, color) {
    data_ = new DelphesCaloData(2);
    data_->RefSliceInfo(0).Setup("ECAL", 0.1, kRed);
    data_->RefSliceInfo(1).Setup("HCAL", 0.1, kBlue);
    data_->IncDenyDestroy();
}
template<> void DelphesBranchElement<DelphesCaloData>::Reset() { data_->ClearTowers(); }
template<> void DelphesBranchElement<DelphesCaloData>::ReadBranch() {
  if(TString(GetType())=="Tower") {
    // Loop over all towers
    TIter itTower(branch_);
    Tower *tower;
    while((tower = (Tower *) itTower.Next())) {
      data_->AddTower(tower->Edges[0], tower->Edges[1], tower->Edges[2], tower->Edges[3]);
      data_->FillSlice(0, tower->Eem);
      data_->FillSlice(1, tower->Ehad);
    }
    data_->DataChanged();
  }
}

// special case for element lists
template<> DelphesBranchElement<TEveElementList>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color):DelphesBranchBase(name, branch, color) {
    data_ = new TEveElementList(name);
    data_->SetMainColor(color_);
}
template<> void DelphesBranchElement<TEveElementList>::Reset() { data_->DestroyElements(); }
template<> void DelphesBranchElement<TEveElementList>::ReadBranch() {
  if(TString(GetType())=="Jet") {
    TIter itJet(branch_);
    Jet *jet;
    TEveJetCone *eveJetCone;
    // Loop over all jets
    Int_t counter = 0;
    while((jet = (Jet *) itJet.Next())) {
      eveJetCone = new TEveJetCone();
      eveJetCone->SetTitle(Form("jet [%d]: Pt=%f, Eta=%f, \nPhi=%f, M=%f",counter,jet->PT, jet->Eta, jet->Phi, jet->Mass));
      eveJetCone->SetName(Form("jet [%d]", counter++));
      eveJetCone->SetMainTransparency(60);
      eveJetCone->SetLineColor(GetColor());
      eveJetCone->SetFillColor(GetColor());
      eveJetCone->SetCylinder(tkRadius_ - 10, tkHalfLength_ - 10);
      eveJetCone->SetPickable(kTRUE);
      eveJetCone->AddEllipticCone(jet->Eta, jet->Phi, jet->DeltaEta, jet->DeltaPhi);
      data_->AddElement(eveJetCone);
    }
  } else if(TString(GetType())=="MissingET") {
    TIter itMet(branch_);
    MissingET *MET;
    TEveArrow *eveMet;
    // Missing Et
    Double_t maxPt = 50.;
    // TODO to be changed as we don't have access to maxPt anymore. MET scale could be a general parameter set in GUI
    while((MET = (MissingET*) itMet.Next())) {
      eveMet = new TEveArrow((tkRadius_ * MET->MET/maxPt)*cos(MET->Phi), (tkRadius_ * MET->MET/maxPt)*sin(MET->Phi), 0., 0., 0., 0.);
      eveMet->SetMainColor(GetColor());
      eveMet->SetTubeR(0.04);
      eveMet->SetConeR(0.08);
      eveMet->SetConeL(0.10);
      eveMet->SetPickable(kTRUE);
      eveMet->SetName("Missing Et");
      eveMet->SetTitle(Form("Missing Et (%.1f GeV)",MET->MET));
      data_->AddElement(eveMet);
    }
  }
}

// special case for track lists
template<> DelphesBranchElement<TEveTrackList>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color):DelphesBranchBase(name, branch, color) {
  data_ = new TEveTrackList(name);
  data_->SetMainColor(color_);
  data_->SetMarkerColor(color_);
  data_->SetMarkerStyle(kCircle);
  data_->SetMarkerSize(0.5);
}
template<> void DelphesBranchElement<TEveTrackList>::SetTrackingVolume(Float_t r, Float_t l, Float_t Bz) {
  tkRadius_ = r; 
  tkHalfLength_ = l;
  tk_Bz_ = Bz;
  TEveTrackPropagator *trkProp = data_->GetPropagator();
  trkProp->SetMagField(0., 0., -tk_Bz_);
  trkProp->SetMaxR(tkRadius_);
  trkProp->SetMaxZ(tkHalfLength_);
}
template<> void DelphesBranchElement<TEveTrackList>::Reset() { data_->DestroyElements(); }
template<> void DelphesBranchElement<TEveTrackList>::ReadBranch() {
  TString type = GetType();
  TIter itTrack(branch_);
  Int_t counter = 0;
  TEveTrack *eveTrack;
  TEveTrackPropagator *trkProp = data_->GetPropagator();
  trkProp->SetMagField(0., 0., -tk_Bz_);
  trkProp->SetMaxR(tkRadius_);
  trkProp->SetMaxZ(tkHalfLength_);
  if(type=="Track") { // CASE 1: TRACKS
    Track *track;
    while((track = (Track *) itTrack.Next())) {
      TParticle pb(track->PID, 1, 0, 0, 0, 0,
                   track->P4().Px(), track->P4().Py(),
                   track->P4().Pz(), track->P4().E(),
                   track->X, track->Y, track->Z, 0.0);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(data_);
      data_->AddElement(eveTrack);
      eveTrack->SetLineColor(GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="Electron") { // CASE 2: ELECTRONS
    Electron *electron;
    while((electron = (Electron *) itTrack.Next())) {
      TParticle pb(electron->Charge<0?11:-11, 1, 0, 0, 0, 0,
                   electron->P4().Px(), electron->P4().Py(),
                   electron->P4().Pz(), electron->P4().E(),
                   0., 0., 0., 0.);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(data_);
      data_->AddElement(eveTrack);
      eveTrack->SetLineColor(GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="Muon") { // CASE 3: MUONS
    Muon *muon;
    while((muon = (Muon *) itTrack.Next())) {
      TParticle pb(muon->Charge<0?13:-13, 1, 0, 0, 0, 0,
                   muon->P4().Px(), muon->P4().Py(),
                   muon->P4().Pz(), muon->P4().E(),
                   0., 0., 0., 0.);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(data_);
      data_->AddElement(eveTrack);
      eveTrack->SetLineColor(GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="Photon") { // CASE 4: PHOTONS
    Photon *photon;
    while((photon = (Photon *) itTrack.Next())) {
      TParticle pb(22, 1, 0, 0, 0, 0,
                   photon->P4().Px(), photon->P4().Py(),
                   photon->P4().Pz(), photon->P4().E(),
                   0., 0., 0., 0.);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(data_);
      eveTrack->SetLineStyle(7);
      data_->AddElement(eveTrack);
      eveTrack->SetLineColor(GetColor());
      eveTrack->MakeTrack();
    }
  } else if(type=="GenParticle") { // CASE 5: GENPARTICLES
    GenParticle *particle;
    while((particle = (GenParticle *) itTrack.Next())) {
      TParticle pb(particle->PID, particle->Status, particle->M1, particle->M2, particle->D1, particle->D2,
                   particle->P4().Px(), particle->P4().Py(),
                   particle->P4().Pz(), particle->P4().E(),
                   particle->X, particle->Y, particle->Z, particle->T);
      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(data_);
      data_->AddElement(eveTrack);
      eveTrack->SetLineColor(GetColor());
      if(particle->Charge==0) eveTrack->SetLineStyle(7);
      eveTrack->MakeTrack();
    }
  }
}
