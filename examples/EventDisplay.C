/* Example:
 * root -l examples/EventDisplay.C'("cards/delphes_card_CMS.tcl","delphes_output.root")'
 * root -l examples/EventDisplay.C'("cards/delphes_card_FCC_basic.tcl","delphes_output.root","ParticlePropagator","ChargedHadronTrackingEfficiency","MuonTrackingEfficiency","Ecal,Hcal")'
 */

void EventDisplay(const char *configfile = "delphes_card_CMS.tcl",
                  const char *datafile = "delphes_output.root",
                  const char *ParticlePropagator = "ParticlePropagator",
                  const char *TrackingEfficiency = "ChargedHadronTrackingEfficiency",
                  const char *MuonEfficiency = "MuonEfficiency",
                  const char *Calorimeters = "Calorimeter",
                  bool displayGeometryOnly = false)
{
  // load the libraries
  gSystem->Load("libGeom");
  gSystem->Load("libGuiHtml");
  gSystem->Load("libDelphesDisplay");

  if(displayGeometryOnly)
  {
    // create the detector representation without transparency
    Delphes3DGeometry det3D_geom(new TGeoManager("delphes", "Delphes geometry"), false);
    det3D_geom.readFile(configfile, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

    // display
    det3D_geom.getDetector()->Draw("ogl");
  } 
  else
  {
    // create the detector representation
    Delphes3DGeometry det3D(new TGeoManager("delphes", "Delphes geometry"), true);
    det3D.readFile(configfile, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

    // create the application
    DelphesEventDisplay* display = new DelphesEventDisplay(configfile, datafile, det3D);
  }
}

