
#include "TEveManager.h"
#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TGLWidget.h"
#include "TEveScene.h"
#include "TEveProjectionManager.h"
#include "TEveProjectionAxes.h"
#include "TEveBrowser.h"
#include "TEveWindow.h"
#include "TGTab.h"

#include "TEveCalo2DGL.h"
#include "TEveCalo3DGL.h"  
#include "TEveCaloLegoGL.h"
#include "TEveCaloLegoOverlay.h"
#include "TEveLegoEventHandler.h"

#include "display/DelphesDisplay.h"

//------------------------------------------------------------------------------

DelphesDisplay::DelphesDisplay()
{
  TEveProjectionAxes *axes;
  TEveWindowSlot *slot;
  TEveWindowPack *packH, *pack0, *pack1;

  fRPhiMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
  fRhoZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);

	fRPhiGeomScene = gEve->SpawnNewScene("RPhi Geometry");
	fRhoZGeomScene = gEve->SpawnNewScene("RhoZ Geometry");

	fRPhiCaloScene = gEve->SpawnNewScene("RPhi Calorimeter");
	fRhoZCaloScene = gEve->SpawnNewScene("RhoZ Calorimeter");
	fLegoCaloScene = gEve->SpawnNewScene("Lego Calorimeter");

	fRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data");
	fRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data");
  
  axes = new TEveProjectionAxes(fRPhiMgr);
  fRPhiGeomScene->AddElement(axes);

  axes = new TEveProjectionAxes(fRhoZMgr);
  fRhoZGeomScene->AddElement(axes);

  // close default tab
  gEve->GetBrowser()->GetTabRight()->CloseTab(0);

  // frames
  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  packH = slot->MakePack();

  packH->SetElementName("Delphes Display");
  packH->SetHorizontal();
  packH->SetShowTitleBar(kFALSE);

  pack0 = packH->NewSlot()->MakePack();
  pack0->SetShowTitleBar(kFALSE);

  pack1 = packH->NewSlot()->MakePack();
  pack1->SetShowTitleBar(kFALSE);

  pack1->NewSlot()->MakeCurrent();
  f3DimView = gEve->SpawnNewViewer("3D View", "");
  f3DimView->AddScene(gEve->GetGlobalScene());
  f3DimView->AddScene(gEve->GetEventScene());
  
  pack1->NewSlot()->MakeCurrent();
  fLegoView = gEve->SpawnNewViewer("Lego View", "");
  fLegoView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fLegoView->AddScene(fLegoCaloScene);
  
  pack0->NewSlot()->MakeCurrent();
  fRPhiView = gEve->SpawnNewViewer("RPhi View", "");
  fRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fRPhiView->AddScene(fRPhiGeomScene);
  fRPhiView->AddScene(fRPhiCaloScene);
  fRPhiView->AddScene(fRPhiEventScene);
        
  pack0->NewSlot()->MakeCurrent();
  fRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
  fRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fRhoZView->AddScene(fRhoZGeomScene);
  fRhoZView->AddScene(fRhoZCaloScene);
  fRhoZView->AddScene(fRhoZEventScene);
}

//------------------------------------------------------------------------------

DelphesDisplay::~DelphesDisplay()
{
}

//------------------------------------------------------------------------------

void DelphesDisplay::ImportGeomRPhi(TEveElement* el)
{ 
  fRPhiMgr->ImportElements(el, fRPhiGeomScene);
}

void DelphesDisplay::ImportGeomRhoZ(TEveElement* el)
{ 
  fRhoZMgr->ImportElements(el, fRhoZGeomScene);
}

void DelphesDisplay::ImportCaloRPhi(TEveCalo3D *calo)
{
  fRPhiMgr->ImportElements(calo, fRPhiCaloScene);
}

void DelphesDisplay::ImportCaloRhoZ(TEveCalo3D *calo)
{
  fRhoZMgr->ImportElements(calo, fRhoZCaloScene);
}

void DelphesDisplay::ImportCaloLego(TEveCaloLego *calo)
{
  TEveCaloLegoOverlay *overlay = new TEveCaloLegoOverlay();
 
  overlay->SetCaloLego(calo);
  fLegoView->GetGLViewer()->AddOverlayElement(overlay);

  fLegoCaloScene->AddElement(calo);
}

void DelphesDisplay::ImportEventRPhi(TEveElement* el)
{ 
  fRPhiMgr->ImportElements(el, fRPhiEventScene);
}

void DelphesDisplay::ImportEventRhoZ(TEveElement* el)
{ 
  fRhoZMgr->ImportElements(el, fRhoZEventScene);
}

//---------------------------------------------------------------------------

void DelphesDisplay::DestroyEventRPhi()
{
  fRPhiEventScene->DestroyElements();
}

void DelphesDisplay::DestroyEventRhoZ()
{
  fRhoZEventScene->DestroyElements();
}
//------------------------------------------------------------------------------

