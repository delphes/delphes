#include <set>
#include <map>
#include <utility>
#include <vector>
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoArb8.h"
//#include "../external/ExRootAnalysis/ExRootConfReader.h"
#include "TF2.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"

using namespace std;

// TODO: asymmetric detector
// TODO: generalize for FCC-like config: >1 calorimeter & flexibility in module names
class Delphes3DGeometry {
   public:
     Delphes3DGeometry(TGeoManager *geom = NULL);
     ~Delphes3DGeometry() {}

     void readFile(const char* filename, const char* ParticlePropagator="ParticlePropagator",
                                         const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                         const char* MuonEfficiency="MuonEfficiency",
                                         const char* Calorimeters="Calorimeter");

     void setContingency(Double_t contingency) { contingency_ = contingency; }
     void setCaloBarrelThickness(Double_t thickness) { calo_barrel_thickness_ = thickness; }
     void setCaloEndcapThickness(Double_t thickness) { calo_endcap_thickness_ = thickness; }
     void setMuonSystemThickness(Double_t thickness) { muonSystem_thickn_ = thickness; }

     TGeoVolume* getDetector(bool withTowers = true);

   private:
     void addTracker(TGeoVolume *top);
     void addCalorimeters(TGeoVolume *top);
     void addMuonDets(TGeoVolume *top);
     void addCaloTowers(TGeoVolume *top);

   private:

     TGeoManager *geom_;

     TGeoMedium *vacuum_;
     TGeoMedium *tkmed_;
     TGeoMedium *calomed_;
     TGeoMedium *mudetmed_;

     Double_t contingency_;
     Double_t calo_barrel_thickness_;
     Double_t calo_endcap_thickness_;
     Double_t muonSystem_thickn_;
     Double_t tk_radius_;
     Double_t tk_length_;
     Double_t tk_etamax_;
     Double_t calo_endcap_etamax_;
     Double_t muonSystem_etamax_;
     Double_t calo_barrel_innerRadius_;
     Double_t calo_endcap_etamin_;
     Double_t calo_endcap_innerRadius1_;
     Double_t calo_endcap_innerRadius2_;
     Double_t calo_endcap_outerRadius1_;
     Double_t calo_endcap_outerRadius2_;
     Double_t calo_endcap_coneThickness_;
     Double_t calo_endcap_diskThickness_;

     set< pair<Double_t, Int_t> > caloBinning_;
     
};

Delphes3DGeometry::Delphes3DGeometry(TGeoManager *geom) {

   //--- the geometry manager
   geom_ = geom==NULL? gGeoManager : geom;

   //--- define some materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7); // placeholder

   //--- define some media
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
   TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
   vacuum_ = Vacuum;
   tkmed_ = Vacuum; // placeholder
   calomed_ = Al;   // placeholder
   mudetmed_ = Al;  // placeholder

   // custom parameters
   contingency_ = 10.;
   calo_barrel_thickness_ = 50.;
   calo_endcap_thickness_ = 75.;
   muonSystem_thickn_ = contingency_;

   // read these parameters from the Delphes Card (with default values)
   tk_radius_ = 120.;
   tk_length_ = 150.;
   tk_etamax_ = 3.0;
   calo_endcap_etamax_ = 2.6;
   muonSystem_etamax_  = 2.4;
}

void Delphes3DGeometry::readFile(const char *configFile,
                                 const char* ParticlePropagator, const char* TrackingEfficiency,
                                 const char* MuonEfficiency, const char* Calorimeters) {

   ExRootConfReader *confReader = new ExRootConfReader;
   confReader->ReadFile(configFile);

   tk_radius_ = confReader->GetDouble(Form("%s::Radius",ParticlePropagator), 1.0)*100;		// tk_radius
   tk_length_ = confReader->GetDouble(Form("%s::HalfLength",ParticlePropagator), 3.0)*100; 	// tk_length

   {
   TString tkEffFormula = confReader->GetString(Form("%s::EfficiencyFormula",TrackingEfficiency),"abs(eta)<3.0");
   tkEffFormula.ReplaceAll("pt","x");
   tkEffFormula.ReplaceAll("eta","y");
   tkEffFormula.ReplaceAll("phi","0.");
   TF2* tkEffFunction = new TF2("tkEff",tkEffFormula,0,1000,-10,10);
   TH1F etaHisto("eta","eta",100,5.,-5.);
   Double_t pt,eta;
   for(int i=0;i<1000;++i) {
     tkEffFunction->GetRandom2(pt,eta);
     etaHisto.Fill(eta);
   }
   Int_t bin = -1;
   bin = etaHisto.FindFirstBinAbove(0.5);
   Double_t etamin = (bin>-1) ? etaHisto.GetBinLowEdge(bin) : -10.;
   bin = etaHisto.FindLastBinAbove(0.5);
   Double_t etamax = (bin>-1) ? etaHisto.GetBinLowEdge(bin+1) : -10.;
   tk_etamax_ = TMath::Max(fabs(etamin),fabs(etamax)); 						// tk_etamax
   delete tkEffFunction;
   }

   {
   TString muonEffFormula = confReader->GetString(Form("%s::EfficiencyFormula",MuonEfficiency),"abs(eta)<2.0");
   muonEffFormula.ReplaceAll("pt","x");
   muonEffFormula.ReplaceAll("eta","y");
   muonEffFormula.ReplaceAll("phi","0.");
   TF2* muEffFunction = new TF2("muEff",muonEffFormula,0,1000,-10,10);
   TH1F etaHisto("eta2","eta2",100,5.,-5.);
   Double_t pt,eta;
   for(int i=0;i<1000;++i) {
     muEffFunction->GetRandom2(pt,eta);
     etaHisto.Fill(eta);
   }
   Int_t bin = -1;
   bin = etaHisto.FindFirstBinAbove(0.5);
   Double_t etamin = (bin>-1) ? etaHisto.GetBinLowEdge(bin) : -10.;
   bin = etaHisto.FindLastBinAbove(0.5);
   Double_t etamax = (bin>-1) ? etaHisto.GetBinLowEdge(bin+1) : -10.;
   muonSystem_etamax_ = TMath::Max(fabs(etamin),fabs(etamax));				// muonSystem_etamax
   delete muEffFunction;
   }
       
   caloBinning_.clear();								// calo binning
   ExRootConfParam paramEtaBins, paramPhiBins;
   ExRootConfParam param = confReader->GetParam(Form("%s::EtaPhiBins",Calorimeters));
   Int_t size = param.GetSize();
   for(int i = 0; i < size/2; ++i) {
     paramEtaBins = param[i*2];
     paramPhiBins = param[i*2+1];
     assert(paramEtaBins.GetSize()==1);
     caloBinning_.insert(std::make_pair(paramEtaBins[0].GetDouble(),paramPhiBins.GetSize()-1));
   }

   if (size>0) calo_endcap_etamax_ = TMath::Max(fabs(caloBinning_.begin()->first),fabs(caloBinning_.rbegin()->first)); // calo_endcap_etamax_

   delete confReader;

   calo_barrel_innerRadius_   = tk_radius_+contingency_;
   calo_endcap_etamin_        = -log(tk_radius_/(2*tk_length_));
   calo_endcap_innerRadius1_  = tk_length_*2.*exp(-calo_endcap_etamax_)/(1-exp(-2.*calo_endcap_etamax_));
   calo_endcap_innerRadius2_  = (tk_length_+calo_endcap_thickness_)*2.*exp(-calo_endcap_etamax_)/(1-exp(-2.*calo_endcap_etamax_));
   calo_endcap_outerRadius1_  = tk_radius_;
   calo_endcap_outerRadius2_  = tk_radius_+calo_barrel_thickness_;
   calo_endcap_coneThickness_ = calo_barrel_thickness_ * (1-exp(-2.*calo_endcap_etamin_)) / (2.*exp(-calo_endcap_etamin_));
   calo_endcap_diskThickness_ = TMath::Max(0.,calo_endcap_thickness_-calo_endcap_coneThickness_);
}

TGeoVolume* Delphes3DGeometry::getDetector(bool withTowers) {
   TGeoVolume *top = geom_->MakeBox("Delphes3DGeometry", vacuum_, 1500, 1500, 2300); // determine the size from what we know about the detector TODO
   addTracker(top);
   addCalorimeters(top); // TODO: allow for more than one calo
   addMuonDets(top);
   if (withTowers) {
     addCaloTowers(top);
   }
   return top;
}

//TODO: there should be a cut by two cones to limit the acceptance in eta
//typically from ChargedHadronTrackingEfficiency 
void Delphes3DGeometry::addTracker(TGeoVolume *top) {
   // tracker: a cylinder
   TGeoVolume *tracker = geom_->MakeTube("tracker", tkmed_, 0., tk_radius_, tk_length_);
   tracker->SetLineColor(kYellow);
   top->AddNode(tracker,1);
}

void Delphes3DGeometry::addCalorimeters(TGeoVolume *top) {
   // calorimeters: tube truncated in eta + cones
   
   /*TGeoTube *calo_barrel_cylinder =*/ new TGeoTube("calo_barrel_cylinder",calo_barrel_innerRadius_,tk_radius_+calo_barrel_thickness_+contingency_,tk_length_+calo_barrel_thickness_);
   /*TGeoCone *calo_endcap_cone =*/ new TGeoCone("calo_endcap_cone",calo_endcap_coneThickness_/2.,calo_endcap_innerRadius1_,calo_endcap_outerRadius1_,calo_endcap_innerRadius2_,calo_endcap_outerRadius2_);
   /*TGeoTube *calo_endcap_disk =*/ new TGeoTube("calo_endcap_disk",calo_endcap_innerRadius2_,tk_radius_+calo_barrel_thickness_,calo_endcap_diskThickness_/2.);
   TGeoTranslation *tr1 = new TGeoTranslation("tr1",0., 0., (calo_endcap_coneThickness_+calo_endcap_diskThickness_)/2.);
   tr1->RegisterYourself();
   TGeoCompositeShape *calo_endcap_cs = new TGeoCompositeShape("calo_endcap_cs","calo_endcap_cone+calo_endcap_disk:tr1");
   TGeoTranslation *trc1 = new TGeoTranslation("calo_endcap1_position",0.,0., tk_length_+calo_endcap_coneThickness_/2.);
   trc1->RegisterYourself();
   TGeoRotation *negz = new TGeoRotation("negz",0,180,0);
   TGeoCombiTrans  *trc2 = new TGeoCombiTrans("calo_endcap2_position",0.,0.,-(tk_length_+calo_endcap_coneThickness_/2.),negz);
   trc2->RegisterYourself();
   TGeoTranslation *trc1c = new TGeoTranslation("calo_endcap1_position_cont",0.,0., tk_length_+calo_endcap_coneThickness_/2.+contingency_);
   trc1c->RegisterYourself();
   TGeoCombiTrans  *trc2c = new TGeoCombiTrans("calo_endcap2_position_cont",0.,0.,-(tk_length_+calo_endcap_coneThickness_/2.)-contingency_,negz);
   trc2c->RegisterYourself();
   TGeoVolume *calo_endcap = new TGeoVolume("calo_endcap",calo_endcap_cs,calomed_);
   TGeoCompositeShape *calo_barrel_cs = new TGeoCompositeShape("calo_barrel_cs","calo_barrel_cylinder-calo_endcap_cs:calo_endcap1_position-calo_endcap_cs:calo_endcap2_position");
   TGeoVolume *calo_barrel = new TGeoVolume("calo_barrel",calo_barrel_cs,calomed_);
   calo_endcap->SetLineColor(kViolet);
   calo_endcap->SetFillColor(kViolet);
   calo_barrel->SetLineColor(kRed);
   top->AddNode(calo_endcap,1,trc1c);
   top->AddNode(calo_endcap,2,trc2c);
   top->AddNode(calo_barrel,1); 
}

void Delphes3DGeometry::addMuonDets(TGeoVolume *top) {

   // muon system: tube + disks
   Double_t muonSystem_radius = tk_radius_+calo_barrel_thickness_+2*contingency_;
   Double_t muonSystem_length = tk_length_+TMath::Max(calo_endcap_coneThickness_,calo_endcap_thickness_)+2*contingency_;
   Double_t muonSystem_rmin   = muonSystem_length*2.*exp(-muonSystem_etamax_)/(1-exp(-2.*muonSystem_etamax_));
   TGeoVolume *muon_barrel = geom_->MakeTube("muon_barrel",mudetmed_,muonSystem_radius,muonSystem_radius+muonSystem_thickn_,muonSystem_length);
   muon_barrel->SetLineColor(kBlue);
   top->AddNode(muon_barrel,1);
   TGeoVolume *muon_endcap = geom_->MakeTube("muon_endcap",mudetmed_,muonSystem_rmin,muonSystem_radius+muonSystem_thickn_,muonSystem_thickn_/2.);
   muon_endcap->SetLineColor(kBlue);
   TGeoTranslation *trm1 = new TGeoTranslation("muonEndcap1_position",0.,0.,muonSystem_length);
   trm1->RegisterYourself();
   TGeoTranslation *trm2 = new TGeoTranslation("muonEndcap2_position",0.,0.,-muonSystem_length);
   trm1->RegisterYourself();
   top->AddNode(muon_endcap,1,trm1);
   top->AddNode(muon_endcap,1,trm2);
}

void Delphes3DGeometry::addCaloTowers(TGeoVolume *top) {

   TGeoVolume* calo_endcap = top->GetNode("calo_endcap_1")->GetVolume();
   TGeoVolume* calo_barrel = top->GetNode("calo_barrel_1")->GetVolume();

   // calo towers in the barrel
   Double_t vertices[16] = {0.,0.,0.,0.,0.,0.,0.,0.}; // summit of the pyramid
   Double_t R  = tk_radius_+calo_barrel_thickness_+2*contingency_; // radius of the muons system = height of the pyramid
   Int_t nEtaBins = caloBinning_.size();
   // this rotation is to make the tower point "up"
   TGeoRotation* initTowerRot = new TGeoRotation("initTowerRot",0.,90.,0.);
   TGeoCombiTrans* initTower  = new TGeoCombiTrans("initTower",0.,-R/2.,0.,initTowerRot);
   initTower->RegisterYourself();
   // eta bins... we build one pyramid per eta slice and then translate it nphi times.
   // phi bins represented by rotations around z
   Double_t *y = new Double_t[nEtaBins];
   Double_t *dx = new Double_t[nEtaBins];
   Int_t *nphi = new Int_t[nEtaBins];
   Int_t etaslice = 0;
   std::map<std::pair<int,int>, TGeoRotation*> phirotations;
   for(set< pair<Double_t, Int_t> >::const_iterator bin=caloBinning_.begin(); bin!=caloBinning_.end();++bin) {
     if(abs(bin->first)>calo_endcap_etamin_) continue; // only in the barrel
     nphi[etaslice] = bin->second;
     y[etaslice] = 0.5*R*(1-exp(-2*bin->first))/exp(-bin->first);
     Double_t phiRotationAngle = 360./nphi[etaslice];
     dx[etaslice] = R*tan(TMath::Pi()*phiRotationAngle/360.);
     for(int phislice=0;phislice<nphi[etaslice];++phislice) {
       phirotations[make_pair(etaslice,phislice)] = new TGeoRotation(Form("phi%d_%d",etaslice,phislice),phiRotationAngle*phislice,0.,0.);
       phirotations[make_pair(etaslice,phislice)]->RegisterYourself();
     }
     ++etaslice;
   }
   nEtaBins = ++etaslice;
   for(int i=0;i<nEtaBins-1;++i) { // loop on the eta slices
     vertices[8]  = -dx[i]; vertices[9]  = y[i];
     vertices[10] = -dx[i]; vertices[11] = y[i+1];
     vertices[12] =  dx[i]; vertices[13] = y[i+1];
     vertices[14] =  dx[i]; vertices[15] = y[i];
     /*TGeoArb8 *tower =*/ new TGeoArb8(Form("tower%d",i),R/2., vertices); // tower in the proper eta slice, at phi=0
     // intersection between the tower and the calo_barrel
     TGeoCompositeShape *finaltower_cs = new TGeoCompositeShape(Form("ftower%d_cs",i),Form("tower%d:initTower*calo_barrel_cs",i));
     TGeoVolume *finaltower = new TGeoVolume(Form("ftower%d",i),finaltower_cs,calomed_);
     finaltower->SetLineColor(kRed);
     for(int j=0;j<nphi[i];++j) { // loop on the phi slices
       calo_barrel->AddNode(finaltower,j,phirotations[make_pair(i,j)]);
     }
   }
   delete[] y;
   delete[] dx;
   delete[] nphi;
   //the towers in the forward region
   R  = tk_length_+calo_endcap_thickness_+3*contingency_; // Z of the muons system = height of the pyramid
   nEtaBins = caloBinning_.size();
   // translation to bring the origin of the tower to (0,0,0)
   TGeoTranslation* towerdz = new TGeoTranslation("towerdz",0.,0.,R/2.-(tk_length_+calo_endcap_coneThickness_/2.));
   towerdz->RegisterYourself();
   // eta bins... we build one pyramid per eta slice and then translate it nphi times.
   Double_t *r = new Double_t[nEtaBins];
   nphi = new Int_t[nEtaBins];
   etaslice = 0;
   phirotations.clear();
   for(set< pair<Double_t, Int_t> >::const_iterator bin=caloBinning_.begin(); bin!=caloBinning_.end();++bin) {
     if(bin->first<calo_endcap_etamin_) continue; // only in the + endcap
     r[etaslice] = R*2*exp(-bin->first)/(1-exp(-2*bin->first)); 
     nphi[etaslice] = bin->second;
     Double_t phiRotationAngle = 360./nphi[etaslice];
     for(int phislice=0;phislice<nphi[etaslice];++phislice) {
       phirotations[make_pair(etaslice,phislice)] = new TGeoRotation(Form("forward_phi%d_%d",etaslice,phislice),phiRotationAngle*phislice,0.,0.);
       phirotations[make_pair(etaslice,phislice)]->RegisterYourself();
     }
     ++etaslice;
   }
   nEtaBins = etaslice;
   for(int i=0;i<nEtaBins;++i) { // loop on the eta slices
     vertices[8]  = -r[i+1]*sin(TMath::Pi()/20.); vertices[9]  = r[i+1]*cos(TMath::Pi()/20.);
     vertices[10] = -r[i]*sin(TMath::Pi()/20.);   vertices[11] = r[i]*cos(TMath::Pi()/20.);
     vertices[12] =  r[i]*sin(TMath::Pi()/20.);   vertices[13] = r[i]*cos(TMath::Pi()/20.);
     vertices[14] =  r[i+1]*sin(TMath::Pi()/20.); vertices[15] = r[i+1]*cos(TMath::Pi()/20.);
     /*TGeoArb8 *fwdtower =*/ new TGeoArb8(Form("fwdtower%d",i),R/2., vertices); // tower in the proper eta slice, at phi=0
     // intersection between the tower and the calo_endcap
     TGeoCompositeShape *finalfwdtower_cs = new TGeoCompositeShape(Form("ffwdtower%d_cs",i),Form("fwdtower%d:towerdz*calo_endcap_cs",i));
     TGeoVolume *finalfwdtower = new TGeoVolume(Form("ffwdtower%d",i),finalfwdtower_cs,calomed_);
     finalfwdtower->SetLineColor(kViolet); 
     for(int j=0;j<nphi[i];++j) { // loop on the phi slices
       calo_endcap->AddNode(finalfwdtower,j,phirotations[make_pair(i,j)]);
     }
   }
   delete[] r;
   delete[] nphi;
}

void geometry(const char* filename = "delphes_card_CMS.tcl", const char* ParticlePropagator="ParticlePropagator",
                                                             const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                                             const char* MuonEfficiency="MuonEfficiency",
                                                             const char* Calorimeters="Calorimeter") 
{
   gSystem->Load("libGeom");
   gSystem->Load("../libDelphes");
   TGeoManager *geom = new TGeoManager("delphes", "Delphes geometry");

   // make the top container volume -> designed to contain a "big" detector (ATLAS)
   TGeoVolume *top = geom->MakeBox("TOP", 0, 1500, 1500, 2300);
   geom->SetTopVolume(top);

   // build the detector
   Delphes3DGeometry det3D;
   det3D.readFile(filename,ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);
   top->AddNode(det3D.getDetector(true),1);

   // draw it
   geom->CloseGeometry();
   top->Draw();
}

