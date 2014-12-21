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

