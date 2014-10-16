
void evdisplay(const char* filename = "delphes_card_CMS.tcl", const char* ParticlePropagator="ParticlePropagator",
                                                              const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                                              const char* MuonEfficiency="MuonEfficiency",
                                                              const char* Calorimeters="Calorimeter") 
{

   // load the libraries
   gSystem->Load("libGeom");
   gSystem->Load("../libDelphesDisplay");

   // create the detector representation
   Delphes3DGeometry det3D(new TGeoManager("delphes", "Delphes geometry"));
   det3D.readFile(filename, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

   // create the application items
   DelphesEventDisplay* display = new DelphesEventDisplay("delphes_card_CMS.tcl", "../delphes_output.root", det3D);

/*
   // EClipType not exported to CINT (see TGLUtil.h):
   // 0 - no clip, 1 - clip plane, 2 - clip box
   TGLViewer *v = gEve->GetDefaultGLViewer();
   //Double_t plane[4] = { 0., 1., 0., 0. };
   //v->GetClipSet()->SetClipState(1,plane);
   //v->GetClipSet()->SetClipType(1);
   //v->ColorSet().Background().SetColor(kMagenta+4);
   //v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
   v->RefreshPadEditor(v);
   v->CurrentCamera().RotateRad(-1.2, 0.5);
   v->DoDraw();
*/
}

