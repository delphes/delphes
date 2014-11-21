/* Example:
 * root -l 'EventDisplay.C("delphes_card_CMS.tcl","../delphes_output.root")'
 */

void EventDisplay(const char* configfile = "delphes_card_CMS.tcl", const char* datafile = "delphes_output.root",
                  const char* ParticlePropagator="ParticlePropagator",
                  const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                  const char* MuonEfficiency="MuonEfficiency",
                  const char* Calorimeters="Calorimeter") 
{
   // load the libraries
   gSystem->Load("libGeom");
   gSystem->Load("libGuiHtml");
   gSystem->Load("../libDelphesDisplay");

   // create the detector representation
   Delphes3DGeometry det3D(new TGeoManager("delphes", "Delphes geometry"), true);
   det3D.readFile(configfile, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

   // create the application
   DelphesEventDisplay* display = new DelphesEventDisplay(configfile, datafile, det3D);

   // another view of the geometry, in another window and without transparency
   Delphes3DGeometry det3D_geom(new TGeoManager("delphes", "Delphes geometry"), false);
   det3D_geom.readFile(configfile, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);
   new TCanvas;
   det3D_geom.getDetector()->Draw();
}

