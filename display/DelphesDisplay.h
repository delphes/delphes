#ifndef DelphesDisplay_h
#define DelphesDisplay_h

class TEveProjectionManager;
class TEveElement;
class TEveCalo3D;
class TEveCaloLego;
class TEveViewer;
class TEveScene;

class DelphesDisplay
{
public:
  
  DelphesDisplay();

  virtual ~DelphesDisplay();

  void ImportGeomRPhi(TEveElement* el);
  void ImportGeomRhoZ(TEveElement* el);

  void ImportCaloRPhi(TEveCalo3D *calo);
  void ImportCaloRhoZ(TEveCalo3D *calo);
  void ImportCaloLego(TEveCaloLego *calo);

  void ImportEventRPhi(TEveElement* el);
  void ImportEventRhoZ(TEveElement* el);

  void DestroyEventRPhi();
  void DestroyEventRhoZ();

private:

  TEveProjectionManager *fRPhiMgr;
  TEveProjectionManager *fRhoZMgr;

  TEveViewer *fRPhiView;
  TEveViewer *fRhoZView;
  TEveViewer *f3DimView;
  TEveViewer *fLegoView;

  TEveScene *fRPhiGeomScene;
  TEveScene *fRhoZGeomScene;

  TEveScene *fRPhiCaloScene;
  TEveScene *fRhoZCaloScene;
  TEveScene *fLegoCaloScene;

  TEveScene *fRPhiEventScene;
  TEveScene *fRhoZEventScene;
};

#endif /* DelphesDisplay_h */

