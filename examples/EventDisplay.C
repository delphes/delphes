/*
root -l examples/EventDisplay.C\(\"examples/delphes_card_CMS.tcl\",\"delphes_output.root\"\)
ShowEvent(1);
ShowEvent(2);
*/

//------------------------------------------------------------------------------

// radius of the barrel, in m
Double_t gRadius = 1.29;

// half-length of the barrel, in m
Double_t gHalfLength = 3.0;

// magnetic field
Double_t gBz = 3.8;

TAxis *gEtaAxis = 0;
TAxis *gPhiAxis = 0;

//------------------------------------------------------------------------------

#include <set>
#include <vector>

using namespace std;

class ExRootTreeReader;
class DelphesCaloData;
class DelphesDisplay;

TChain gChain("Delphes");

ExRootTreeReader *gTreeReader = 0;

TClonesArray *gBranchTower = 0;
TClonesArray *gBranchTrack = 0;
TClonesArray *gBranchJet = 0;

DelphesCaloData *gCaloData = 0;
TEveElementList *gJetList = 0;
TEveTrackList *gTrackList = 0;

DelphesDisplay *gDelphesDisplay = 0;

//------------------------------------------------------------------------------

void EventDisplay(const char *configFile, const char *inputFile)
{
  gSystem->Load("libDelphesDisplay");

  TEveManager::Create(kTRUE, "IV");

  ExRootConfParam param, paramEtaBins;
  Long_t i, j, size, sizeEtaBins;
  set< Double_t > etaSet;
  set< Double_t >::iterator itEtaSet;

  Double_t *etaBins;

  ExRootConfReader *confReader = new ExRootConfReader;
  confReader->ReadFile(configFile);

  gRadius = confReader->GetDouble("ParticlePropagator::Radius", 1.0);
  gHalfLength = confReader->GetDouble("ParticlePropagator::HalfLength", 3.0);
  gBz = confReader->GetDouble("ParticlePropagator::Bz", 0.0);

  // read eta and phi bins
  param = confReader->GetParam("Calorimeter::EtaPhiBins");
  size = param.GetSize();
  etaSet.clear();
  for(i = 0; i < size/2; ++i)
  {
    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j)
    {
      etaSet.insert(paramEtaBins[j].GetDouble());
    }
  }

  delete confReader;

  etaBins = new Double_t[etaSet.size()];
  i = 0;

  for(itEtaSet = etaSet.begin(); itEtaSet != etaSet.end(); ++itEtaSet)
  {
    etaBins[i] = *itEtaSet;
    ++i;
  }

  gEtaAxis = new TAxis(etaSet.size() - 1, etaBins);
  gPhiAxis = new TAxis(72, -TMath::Pi(), TMath::Pi());

  // Create chain of root trees
  gChain.Add(inputFile);

  // Create object of class ExRootTreeReader
  gTreeReader = new ExRootTreeReader(&gChain);

  // Get pointers to branches
  gBranchTower = gTreeReader->UseBranch("Tower");
  gBranchTrack = gTreeReader->UseBranch("Track");
  gBranchJet = gTreeReader->UseBranch("Jet");

  // data
  gCaloData = new DelphesCaloData(2);
  gCaloData->RefSliceInfo(0).Setup("ECAL", 0.1, kRed);
  gCaloData->RefSliceInfo(1).Setup("HCAL", 0.1, kBlue);
  gCaloData->SetEtaBins(gEtaAxis);
  gCaloData->SetPhiBins(gPhiAxis);
  gCaloData->IncDenyDestroy();

  gJetList = new TEveElementList("Jets");
  gEve->AddElement(gJetList);

  gTrackList = new TEveTrackList("Tracks");
  gTrackList->SetMainColor(kBlue);
  gTrackList->SetMarkerColor(kRed);
  gTrackList->SetMarkerStyle(kCircle);
  gTrackList->SetMarkerSize(0.5);
  gEve->AddElement(gTrackList);

  TEveTrackPropagator *trkProp = gTrackList->GetPropagator();
  trkProp->SetMagField(0.0, 0.0, -gBz);
  trkProp->SetMaxR(gRadius*100.0);
  trkProp->SetMaxZ(gHalfLength*100.0);

  // viewers and scenes

  TEveElementList *geometry = new TEveElementList("Geometry");

  TEveGeoShape *barell = new TEveGeoShape("Barell");
  barell->SetShape(new TGeoTube(gRadius*100.0 - 1, gRadius*100.0, gHalfLength*100.0));
  barell->SetMainColor(kCyan);
  barell->SetMainTransparency(80);
  geometry->AddElement(barell);

  TEveCalo3D *calo = new TEveCalo3D(gCaloData);
  calo->SetBarrelRadius(gRadius*100.0);
  calo->SetEndCapPos(gHalfLength*100.0);

  gStyle->SetPalette(1, 0);
  TEveCaloLego *lego = new TEveCaloLego(gCaloData);
  lego->InitMainTrans();
  lego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
  lego->SetAutoRebin(kFALSE);
  lego->Set2DMode(TEveCaloLego::kValSizeOutline);

  gDelphesDisplay = new DelphesDisplay;

  gEve->AddGlobalElement(geometry);
  gEve->AddGlobalElement(calo);

  gDelphesDisplay->ImportGeomRPhi(geometry);
  gDelphesDisplay->ImportCaloRPhi(calo);

  gDelphesDisplay->ImportGeomRhoZ(geometry);
  gDelphesDisplay->ImportCaloRhoZ(calo);

  gDelphesDisplay->ImportCaloLego(lego);

  gEve->Redraw3D(kTRUE);
}

//------------------------------------------------------------------------------

void ShowEvent(Long64_t event)
{
  TIter itTower(gBranchTower);
  TIter itTrack(gBranchTrack);
  TIter itJet(gBranchJet);

  Tower *tower;
  Track *track;
  Jet *jet;

  TEveJetCone *eveJetCone;
  TEveTrack *eveTrack;

  Int_t counter;

  TEveElement *currentEvent = gEve->GetCurrentEvent();

  TEveTrackPropagator *trkProp = gTrackList->GetPropagator();

  if(event >= gTreeReader->GetEntries()) return;

  // Load selected branches with data from specified event
  gTreeReader->ReadEntry(event);

  gCaloData->ClearTowers();
  gJetList->DestroyElements();
  gTrackList->DestroyElements();

  // Loop over all towers
  itTower.Reset();
  while((tower = (Tower *) itTower.Next()))
  {
    gCaloData->AddTower(tower->Edges[0], tower->Edges[1], tower->Edges[2], tower->Edges[3]);
    gCaloData->FillSlice(0, tower->Eem);
    gCaloData->FillSlice(1, tower->Ehad);
  }
  gCaloData->DataChanged();

  // Loop over all tracks
  itTrack.Reset();
  counter = 0;
  while((track = (Track *) itTrack.Next()))
  {
    TParticle pb(track->PID, 1, 0, 0, 0, 0,
                 track->P4().Px(), track->P4().Py(),
                 track->P4().Pz(), track->P4().E(),
                 track->X, track->Y, track->Z, 0.0);

    eveTrack = new TEveTrack(&pb, counter, trkProp);
    eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
    eveTrack->SetStdTitle();
    eveTrack->SetAttLineAttMarker(gTrackList);

    switch(TMath::Abs(track->PID))
    {
      case 11:
        eveTrack->SetLineColor(kRed);
        break;
      case 13:
        eveTrack->SetLineColor(kGreen);
        break;
      default:
        eveTrack->SetLineColor(kBlue);
    }
    gTrackList->AddElement(eveTrack);
    eveTrack->MakeTrack();
  }

  // Loop over all jets
  itJet.Reset();
  counter = 0;
  while((jet = (Jet *) itJet.Next()))
  {
    eveJetCone = new TEveJetCone();
    eveJetCone->SetName(Form("jet [%d]", counter++));
    eveJetCone->SetMainTransparency(60);
    eveJetCone->SetLineColor(kYellow);
    eveJetCone->SetCylinder(gRadius*100.0 - 10, gHalfLength*100.0 - 10);
    eveJetCone->SetPickable(kTRUE);
    eveJetCone->AddEllipticCone(jet->Eta, jet->Phi, jet->DeltaEta, jet->DeltaPhi);
    gJetList->AddElement(eveJetCone);
  }

  gDelphesDisplay->DestroyEventRPhi();
  gDelphesDisplay->ImportEventRPhi(currentEvent);

  gDelphesDisplay->DestroyEventRhoZ();
  gDelphesDisplay->ImportEventRhoZ(currentEvent);

  gEve->Redraw3D(kTRUE);
}
