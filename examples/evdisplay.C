
void evdisplay(const char* filename = "delphes_card_CMS.tcl", const char* ParticlePropagator="ParticlePropagator",
                                                              const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                                              const char* MuonEfficiency="MuonEfficiency",
                                                              const char* Calorimeters="Calorimeter") 
{

   // load the libraries
   gSystem->Load("libGeom");
   //gSystem->Load("../libDelphes");
   gSystem->Load("../libDelphesDisplay");

   // create the detector representation
   Delphes3DGeometry det3D(new TGeoManager("delphes", "Delphes geometry"));
   det3D.readFile(filename, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

   // create the application items
   DelphesEventDisplay display("delphes_card_CMS.tcl", "../delphes_output.root", det3D); //TODO root file as input cfg
/*
   make_gui();
   load_event();
   gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.

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

