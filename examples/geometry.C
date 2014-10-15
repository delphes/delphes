#include <set>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>
#include <exception>
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoArb8.h"
//#include "external/ExRootAnalysis/ExRootConfReader.h"
//#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#include "display/DelphesCaloData.h"
//#include "display/DelphesDisplay.h"
#include "../display/DelphesBranchElement.h"
//#include "classes/DelphesClasses.h"
#include "TF2.h"
#include "TH1F.h"
#include "TChain.h"
#include "TEveElement.h"
#include "TEveJetCone.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveCalo.h"
#include "TMath.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveTrans.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"
#include "TRootBrowser.h"
#include "TGLViewer.h"
#include "TGButton.h"
#include "TCollection.h"
#include "TClonesArray.h"
#include "TGLClip.h"
#include "TEveArrow.h"

/*
 * assembly.C: sauvegarde as shape-extract -> implement in the geometry class (read/write)
 * histobrowser.C: intégration d'histogrammes dans le display (on pourrait avoir Pt, eta, phi pour les principales collections)
 * also from alice_esd: summary html table
 * 
 */
using namespace std;

// Forward declarations.
class Delphes3DGeometry;
class ExRootTreeReader;
class DelphesCaloData;
class DelphesDisplay;
class DelphesBranchBase;
template<typename EveContainer> class DelphesBranchElement;
void make_gui();
void load_event();
void delphes_read();
void delphes_read_towers(TClonesArray* data, DelphesBranchBase* element);
void delphes_read_tracks(TClonesArray* data, DelphesBranchBase* element);
void delphes_read_jets(TClonesArray* data, DelphesBranchBase* element);
void delphes_read_vectors(TClonesArray* data, DelphesBranchBase* element);
void readConfig(const char *configFile, const Delphes3DGeometry& det3D, std::vector<DelphesBranchBase*>& elements, std::vector<TClonesArray*>& arrays);

// Configuration and global variables.
Int_t event_id           = 0; // Current event id.
Double_t gRadius = 1.29;
Double_t gTotRadius = 2.0;
Double_t gHalfLength = 3.0;
Double_t gBz = 3.8;

TChain gChain("Delphes");

ExRootTreeReader *gTreeReader = 0;

std::vector<DelphesBranchBase*> gElements;
std::vector<TClonesArray*> gArrays;

DelphesDisplay *gDelphesDisplay = 0;

/******************************************************************************/
// Construction of the geometry
/******************************************************************************/

// TODO: asymmetric detector

class Delphes3DGeometry {
   public:
     Delphes3DGeometry(TGeoManager *geom = NULL);
     ~Delphes3DGeometry() {}

     void readFile(const char* filename, const char* ParticlePropagator="ParticlePropagator",
                                         const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                         const char* MuonEfficiency="MuonEfficiency",
                                         const char* Calorimeters="Calorimeter");

     void loadFromFile(const char* filename, const char* name="DelphesGeometry");
     void save(const char* filename, const char* name="DelphesGeometry");

     void setContingency(Double_t contingency) { contingency_ = contingency; }
     void setCaloBarrelThickness(Double_t thickness) { calo_barrel_thickness_ = thickness; }
     void setCaloEndcapThickness(Double_t thickness) { calo_endcap_thickness_ = thickness; }
     void setMuonSystemThickness(Double_t thickness) { muonSystem_thickn_ = thickness; }

     TGeoVolume* getDetector(bool withTowers = true);

     Double_t getTrackerRadius() const { return tk_radius_; }
     Double_t getDetectorRadius() const { return muonSystem_radius_; }
     Double_t getTrackerHalfLength() const { return tk_length_; }
     Double_t getDetectorHalfLength() const { return muonSystem_length_; }
     Double_t getBField() const { return tk_Bz_; }
     std::pair<TAxis*, TAxis*> getCaloAxes() { return std::make_pair(etaAxis_,phiAxis_); }

   private:
     std::pair<Double_t, Double_t> addTracker(TGeoVolume *top);
     std::pair<Double_t, Double_t> addCalorimeter(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength, set< pair<Double_t, Int_t> >& caloBinning);
     std::pair<Double_t, Double_t> addMuonDets(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength);
     void addCaloTowers(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength, set< pair<Double_t, Int_t> >& caloBinning);

   private:

     TGeoManager *geom_;

     TGeoMedium *vacuum_;
     TGeoMedium *tkmed_;
     TGeoMedium *calomed_;
     TGeoMedium *mudetmed_;

     TAxis* etaAxis_;
     TAxis* phiAxis_;

     Double_t contingency_;
     Double_t calo_barrel_thickness_;
     Double_t calo_endcap_thickness_;
     Double_t muonSystem_thickn_;
     Double_t muonSystem_radius_;
     Double_t muonSystem_length_;
     Double_t tk_radius_;
     Double_t tk_length_;
     Double_t tk_etamax_;
     Double_t tk_Bz_;

     std::vector<std::string> calorimeters_;
     std::vector<std::string> muondets_;

     std::map<std::string, Double_t> muonSystem_etamax_;
     std::map<std::string, set< pair<Double_t, Int_t> > > caloBinning_;
     
};

Delphes3DGeometry::Delphes3DGeometry(TGeoManager *geom) {

   //--- the geometry manager
   geom_ = geom==NULL? gGeoManager : geom;
   //gGeoManager->DefaultColors();

   //--- define some materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7); // placeholder
   //TODO: create different materials for different subdetectors???
   matVacuum->SetTransparency(85); //TODO: tune
   matAl->SetTransparency(85); //TODO: tune

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
   muonSystem_thickn_ = 10.;

   // read these parameters from the Delphes Card (with default values)
   etaAxis_   = NULL;
   phiAxis_   = NULL;
   tk_radius_ = 120.;
   tk_length_ = 150.;
   tk_etamax_ = 3.0;
   tk_Bz_     = 1.;
   muonSystem_radius_ = 200.;
}

void Delphes3DGeometry::readFile(const char *configFile,
                                 const char* ParticlePropagator, const char* TrackingEfficiency,
                                 const char* MuonEfficiency, const char* Calorimeters) {

   ExRootConfReader *confReader = new ExRootConfReader;
   confReader->ReadFile(configFile);

   tk_radius_ = confReader->GetDouble(Form("%s::Radius",ParticlePropagator), 1.0)*100.;		// tk_radius
   tk_length_ = confReader->GetDouble(Form("%s::HalfLength",ParticlePropagator), 3.0)*100.; 	// tk_length
   tk_Bz_     = confReader->GetDouble("ParticlePropagator::Bz", 0.0);                           // tk_Bz

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
   muondets_.push_back("muons");
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
   muonSystem_etamax_["muons"] = TMath::Max(fabs(etamin),fabs(etamax));			// muonSystem_etamax
   delete muEffFunction;
   }

   std::string s(Calorimeters);
   std::replace( s.begin(), s.end(), ',', ' ' );
   std::istringstream stream( s );
   std::string word;
   while (stream >> word) calorimeters_.push_back(word);

   caloBinning_.clear();								// calo binning
   for(std::vector<std::string>::const_iterator calo=calorimeters_.begin();calo!=calorimeters_.end(); ++calo) {
     set< pair<Double_t, Int_t> > caloBinning;
     ExRootConfParam paramEtaBins, paramPhiBins;
     ExRootConfParam param = confReader->GetParam(Form("%s::EtaPhiBins",calo->c_str()));
     Int_t size = param.GetSize();
     for(int i = 0; i < size/2; ++i) {
       paramEtaBins = param[i*2];
       paramPhiBins = param[i*2+1];
       assert(paramEtaBins.GetSize()==1);
       caloBinning.insert(std::make_pair(paramEtaBins[0].GetDouble(),paramPhiBins.GetSize()-1));
     }
     caloBinning_[*calo] = caloBinning;
   }

   set< pair<Double_t, Int_t> > caloBinning = caloBinning_[*calorimeters_.begin()];
   Double_t *etaBins = new Double_t[caloBinning.size()]; // note that this is the eta binning of the first calo
   unsigned int ii = 0;
   for(set< pair<Double_t, Int_t> >::const_iterator itEtaSet = caloBinning.begin(); itEtaSet != caloBinning.end(); ++itEtaSet) {
     etaBins[ii++] = itEtaSet->first;
   }
   etaAxis_ = new TAxis(caloBinning.size() - 1, etaBins);
   phiAxis_ = new TAxis(72, -TMath::Pi(), TMath::Pi()); // note that this is fixed while #phibins could vary, also with eta, which doesn't seem possible in ROOT

   muonSystem_radius_ = tk_radius_ + contingency_ + (contingency_+calo_barrel_thickness_)*calorimeters_.size() + muonSystem_thickn_;
   muonSystem_length_ = tk_length_ + contingency_ + (contingency_+calo_endcap_thickness_)*calorimeters_.size() + muonSystem_thickn_;

   delete confReader;

}

TGeoVolume* Delphes3DGeometry::getDetector(bool withTowers) {
   // compute the envelope
   Double_t system_radius = tk_radius_+calo_barrel_thickness_+3*contingency_;
   Double_t system_length = tk_length_+contingency_+(contingency_+calo_endcap_thickness_)*calorimeters_.size()+contingency_;
   // the detector volume
   TGeoVolume *top = geom_->MakeBox("Delphes3DGeometry", vacuum_, system_radius, system_radius, system_length);
   // build the detector
   std::pair<Double_t, Double_t> limits = addTracker(top);
   Double_t radius = limits.first;
   Double_t length = limits.second;
   for(std::vector<std::string>::const_iterator calo = calorimeters_.begin(); calo != calorimeters_.end(); ++calo) {
     limits = addCalorimeter(top,calo->c_str(),radius,length,caloBinning_[*calo]);
     if (withTowers) {
       addCaloTowers(top,calo->c_str(),radius,length,caloBinning_[*calo]);
     }
     radius = limits.first;
     length = limits.second;
   }
   for(std::vector<std::string>::const_iterator muon = muondets_.begin(); muon != muondets_.end(); ++muon) {
     limits = addMuonDets(top, muon->c_str(), radius, length);
     radius = limits.first;
     length = limits.second;
   }
   // return the result
   return top;
}

std::pair<Double_t, Double_t> Delphes3DGeometry::addTracker(TGeoVolume *top) {
   // tracker: a cylinder with two cones substracted
   new TGeoCone("forwardTkAcceptance",(tk_length_/2.+0.05),0.,tk_radius_,(tk_length_)*2.*exp(-tk_etamax_)/(1-exp(-2.*tk_etamax_)),tk_radius_);
   TGeoTranslation *tr1  = new TGeoTranslation("tkacc1",0., 0., tk_length_/2.);
   tr1->RegisterYourself();
   TGeoRotation *negz    = new TGeoRotation("tknegz",0,180,0);
   negz->RegisterYourself();
   TGeoCombiTrans  *tr2  = new TGeoCombiTrans("tkacc2",0.,0.,-tk_length_/2.,negz);
   tr2->RegisterYourself();
   TGeoCompositeShape* tracker_cs = new TGeoCompositeShape("tracker_cs","forwardTkAcceptance:tkacc1+forwardTkAcceptance:tkacc2");
   TGeoVolume *tracker = new TGeoVolume("tracker",tracker_cs,tkmed_);   
   tracker->SetLineColor(kYellow);
   top->AddNode(tracker,1);
   return std::make_pair(tk_radius_,tk_length_);
}

std::pair<Double_t, Double_t> Delphes3DGeometry::addCalorimeter(TGeoVolume *top, const char* name, 
                                                                Double_t innerBarrelRadius, Double_t innerBarrelLength, set< pair<Double_t, Int_t> >& caloBinning) {
   // parameters derived from the inputs
   Double_t calo_endcap_etamax        = TMath::Max(fabs(caloBinning.begin()->first),fabs(caloBinning.rbegin()->first));
   Double_t calo_barrel_innerRadius   = innerBarrelRadius+contingency_;
   Double_t calo_barrel_length        = innerBarrelLength + calo_barrel_thickness_;
   Double_t calo_endcap_etamin        = -log(innerBarrelRadius/(2*innerBarrelLength));
   Double_t calo_endcap_innerRadius1  = innerBarrelLength*2.*exp(-calo_endcap_etamax)/(1-exp(-2.*calo_endcap_etamax));
   Double_t calo_endcap_innerRadius2  = (innerBarrelLength+calo_endcap_thickness_)*2.*exp(-calo_endcap_etamax)/(1-exp(-2.*calo_endcap_etamax));
   Double_t calo_endcap_outerRadius1  = innerBarrelRadius;
   Double_t calo_endcap_outerRadius2  = innerBarrelRadius+calo_barrel_thickness_;
   Double_t calo_endcap_coneThickness = TMath::Min(calo_barrel_thickness_ * (1-exp(-2.*calo_endcap_etamin)) / (2.*exp(-calo_endcap_etamin)), calo_endcap_thickness_);
   Double_t calo_endcap_diskThickness = TMath::Max(0.,calo_endcap_thickness_-calo_endcap_coneThickness);

   // calorimeters: tube truncated in eta + cones
   new TGeoTube(Form("%s_barrel_cylinder",name),calo_barrel_innerRadius,calo_barrel_innerRadius+calo_barrel_thickness_,calo_barrel_length);
   new TGeoCone(Form("%s_endcap_cone",name),calo_endcap_coneThickness/2.,calo_endcap_innerRadius1,calo_endcap_outerRadius1,calo_endcap_innerRadius2,calo_endcap_outerRadius2);
   new TGeoTube(Form("%s_endcap_disk",name),calo_endcap_innerRadius2,tk_radius_+calo_barrel_thickness_,calo_endcap_diskThickness/2.);
   TGeoTranslation *tr1 = new TGeoTranslation(Form("%s_tr1",name),0., 0., (calo_endcap_coneThickness+calo_endcap_diskThickness)/2.);
   tr1->RegisterYourself();
   TGeoCompositeShape *calo_endcap_cs = new TGeoCompositeShape(Form("%s_endcap_cs",name),Form("%s_endcap_cone+%s_endcap_disk:%s_tr1",name,name,name));
   TGeoTranslation *trc1 = new TGeoTranslation(Form("%s_endcap1_position",name),0.,0., innerBarrelLength+calo_endcap_coneThickness/2.);
   trc1->RegisterYourself();
   TGeoRotation *negz = new TGeoRotation(Form("%s_negz",name),0,180,0);
   TGeoCombiTrans  *trc2 = new TGeoCombiTrans(Form("%s_endcap2_position",name),0.,0.,-(innerBarrelLength+calo_endcap_coneThickness/2.),negz);
   trc2->RegisterYourself();
   TGeoTranslation *trc1c = new TGeoTranslation(Form("%s_endcap1_position_cont",name),0.,0., innerBarrelLength+calo_endcap_coneThickness/2.+contingency_);
   trc1c->RegisterYourself();
   TGeoCombiTrans  *trc2c = new TGeoCombiTrans(Form("%s_endcap2_position_cont",name),0.,0.,-(innerBarrelLength+calo_endcap_coneThickness/2.)-contingency_,negz);
   trc2c->RegisterYourself();
   TGeoVolume *calo_endcap = new TGeoVolume(Form("%s_endcap",name),calo_endcap_cs,calomed_);
   TGeoCompositeShape *calo_barrel_cs = new TGeoCompositeShape(Form("%s_barrel_cs",name),
                                                               Form("%s_barrel_cylinder-%s_endcap_cs:%s_endcap1_position-%s_endcap_cs:%s_endcap2_position",name,name,name,name,name));
   TGeoVolume *calo_barrel = new TGeoVolume(Form("%s_barrel",name),calo_barrel_cs,calomed_);
   calo_endcap->SetLineColor(kViolet);
   calo_endcap->SetFillColor(kViolet);
   calo_barrel->SetLineColor(kRed);
   top->AddNode(calo_endcap,1,trc1c);
   top->AddNode(calo_endcap,2,trc2c);
   top->AddNode(calo_barrel,1); 
   return std::make_pair(calo_barrel_innerRadius+calo_barrel_thickness_,innerBarrelLength+calo_endcap_thickness_+contingency_);
}

std::pair<Double_t, Double_t> Delphes3DGeometry::addMuonDets(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength) {
   // muon system: tube + disks
   Double_t muonSystem_radius = innerBarrelRadius + contingency_;
   Double_t muonSystem_length = innerBarrelLength + contingency_;
   Double_t muonSystem_rmin   = muonSystem_length*2.*exp(-muonSystem_etamax_[name])/(1-exp(-2.*muonSystem_etamax_[name]));
   TGeoVolume *muon_barrel = geom_->MakeTube(Form("%s_barrel",name),mudetmed_,muonSystem_radius,muonSystem_radius+muonSystem_thickn_,muonSystem_length);
   muon_barrel->SetLineColor(kBlue);
   top->AddNode(muon_barrel,1);
   TGeoVolume *muon_endcap = geom_->MakeTube(Form("%s_endcap",name),mudetmed_,muonSystem_rmin,muonSystem_radius+muonSystem_thickn_,muonSystem_thickn_/2.);
   muon_endcap->SetLineColor(kBlue);
   TGeoTranslation *trm1 = new TGeoTranslation(Form("%sEndcap1_position",name),0.,0.,muonSystem_length);
   trm1->RegisterYourself();
   TGeoTranslation *trm2 = new TGeoTranslation(Form("%sEndcap2_position",name),0.,0.,-muonSystem_length);
   trm1->RegisterYourself();
   top->AddNode(muon_endcap,1,trm1);
   top->AddNode(muon_endcap,2,trm2);
   return std::make_pair(muonSystem_radius,muonSystem_length);
}

void Delphes3DGeometry::addCaloTowers(TGeoVolume *top, const char* name,
                                                       Double_t innerBarrelRadius, Double_t innerBarrelLength, set< pair<Double_t, Int_t> >& caloBinning) {

   TGeoVolume* calo_endcap = top->GetNode(Form("%s_endcap_1",name))->GetVolume();
   TGeoVolume* calo_barrel = top->GetNode(Form("%s_barrel_1",name))->GetVolume();
   Double_t calo_endcap_etamin = -log(innerBarrelRadius/(2*innerBarrelLength));
   Double_t calo_endcap_coneThickness = TMath::Min(calo_barrel_thickness_ * (1-exp(-2.*calo_endcap_etamin)) / (2.*exp(-calo_endcap_etamin)), calo_endcap_thickness_);

   // calo towers in the barrel
   Double_t vertices[16] = {0.,0.,0.,0.,0.,0.,0.,0.}; // summit of the pyramid
   Double_t R  = tk_radius_ + contingency_+(contingency_+calo_barrel_thickness_)*calorimeters_.size(); // radius of the muons system = height of the pyramid
   Int_t nEtaBins = caloBinning.size();
   // this rotation is to make the tower point "up"
   TGeoRotation* initTowerRot = new TGeoRotation(Form("%s_initTowerRot",name),0.,90.,0.);
   TGeoCombiTrans* initTower  = new TGeoCombiTrans(Form("%s_initTower",name),0.,-R/2.,0.,initTowerRot);
   initTower->RegisterYourself();
   // eta bins... we build one pyramid per eta slice and then translate it nphi times.
   // phi bins represented by rotations around z
   Double_t *y = new Double_t[nEtaBins];
   Double_t *dx = new Double_t[nEtaBins];
   Int_t *nphi = new Int_t[nEtaBins];
   Int_t etaslice = 0;
   std::map<std::pair<int,int>, TGeoRotation*> phirotations;
   for(set< pair<Double_t, Int_t> >::const_iterator bin=caloBinning.begin(); bin!=caloBinning.end();++bin) {
     if(abs(bin->first)>calo_endcap_etamin) continue; // only in the barrel
     nphi[etaslice] = bin->second;
     y[etaslice] = 0.5*R*(1-exp(-2*bin->first))/exp(-bin->first);
     Double_t phiRotationAngle = 360./nphi[etaslice];
     dx[etaslice] = R*tan(TMath::Pi()*phiRotationAngle/360.);
     for(int phislice=0;phislice<nphi[etaslice];++phislice) {
       phirotations[make_pair(etaslice,phislice)] = new TGeoRotation(Form("%s_phi%d_%d",name,etaslice,phislice),phiRotationAngle*phislice,0.,0.);
       phirotations[make_pair(etaslice,phislice)]->RegisterYourself();
     }
     ++etaslice;
   }
   nEtaBins = etaslice;
   for(int i=0;i<nEtaBins-1;++i) { // loop on the eta slices
     vertices[8]  = -dx[i]; vertices[9]  = y[i];
     vertices[10] = -dx[i]; vertices[11] = y[i+1];
     vertices[12] =  dx[i]; vertices[13] = y[i+1];
     vertices[14] =  dx[i]; vertices[15] = y[i];
     new TGeoArb8(Form("%s_tower%d",name,i),R/2., vertices); // tower in the proper eta slice, at phi=0
     // intersection between the tower and the calo_barrel
     TGeoCompositeShape *finaltower_cs = new TGeoCompositeShape(Form("%s_ftower%d_cs",name,i),Form("%s_tower%d:%s_initTower*%s_barrel_cs",name,i,name,name));
     TGeoVolume *finaltower = new TGeoVolume(Form("%s_ftower%d",name,i),finaltower_cs,calomed_);
     finaltower->SetLineColor(kRed);
     for(int j=0;j<nphi[i];++j) { // loop on the phi slices
       calo_barrel->AddNode(finaltower,j,phirotations[make_pair(i,j)]);
     }
   }
   delete[] y;
   delete[] dx;
   delete[] nphi;
   //the towers in the forward region
   R  = tk_length_+contingency_+(contingency_+calo_endcap_thickness_)*calorimeters_.size(); // Z of the muons system = height of the pyramid
   nEtaBins = caloBinning.size();
   // translation to bring the origin of the tower to (0,0,0) (well, not really as the endcap is not yet in place)
   TGeoTranslation* towerdz = new TGeoTranslation(Form("%s_towerdz",name),0.,0.,R/2.-(innerBarrelLength+calo_endcap_coneThickness/2.));
   towerdz->RegisterYourself();
   // eta bins... we build one pyramid per eta slice and then translate it nphi times.
   Double_t *r = new Double_t[nEtaBins];
   nphi = new Int_t[nEtaBins];
   etaslice = 0;
   phirotations.clear();
   for(set< pair<Double_t, Int_t> >::const_iterator bin=caloBinning.begin(); bin!=caloBinning.end();++bin) {
     if(bin->first<calo_endcap_etamin) continue; // only in the + endcap
     r[etaslice] = R*2*exp(-bin->first)/(1-exp(-2*bin->first)); 
     nphi[etaslice] = bin->second;
     Double_t phiRotationAngle = 360./nphi[etaslice];
     for(int phislice=0;phislice<nphi[etaslice];++phislice) {
       phirotations[make_pair(etaslice,phislice)] = new TGeoRotation(Form("%s_forward_phi%d_%d",name,etaslice,phislice),phiRotationAngle*phislice,0.,0.);
       phirotations[make_pair(etaslice,phislice)]->RegisterYourself();
     }
     ++etaslice;
   }
   nEtaBins = etaslice;
   for(int i=0;i<nEtaBins-1;++i) { // loop on the eta slices
     vertices[8]  = -r[i+1]*sin(TMath::Pi()/nphi[i]); vertices[9]  = r[i+1]*cos(TMath::Pi()/nphi[i]);
     vertices[10] = -r[i]*sin(TMath::Pi()/nphi[i]);   vertices[11] = r[i]*cos(TMath::Pi()/nphi[i]);
     vertices[12] =  r[i]*sin(TMath::Pi()/nphi[i]);   vertices[13] = r[i]*cos(TMath::Pi()/nphi[i]);
     vertices[14] =  r[i+1]*sin(TMath::Pi()/nphi[i]); vertices[15] = r[i+1]*cos(TMath::Pi()/nphi[i]);
     new TGeoArb8(Form("%sfwdtower%d",name,i),R/2., vertices); // tower in the proper eta slice, at phi=0
     // intersection between the tower and the calo_endcap
     TGeoCompositeShape *finalfwdtower_cs = new TGeoCompositeShape(Form("%sffwdtower%d_cs",name,i),Form("%sfwdtower%d:%s_towerdz*%s_endcap_cs",name,i,name,name));
     TGeoVolume *finalfwdtower = new TGeoVolume(Form("%sffwdtower%d",name,i),finalfwdtower_cs,calomed_);
     finalfwdtower->SetLineColor(kViolet); 
     for(int j=0;j<nphi[i];++j) { // loop on the phi slices
       calo_endcap->AddNode(finalfwdtower,j,phirotations[make_pair(i,j)]);
     }
   }
   delete[] r;
   delete[] nphi;
}

/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

// function that parses the config to extract the branches of interest and prepare containers
void readConfig(const char *configFile, Delphes3DGeometry& det3D, std::vector<DelphesBranchBase*>& elements, std::vector<TClonesArray*>& arrays) {
   ExRootConfReader *confReader = new ExRootConfReader;
   confReader->ReadFile(configFile);
   Double_t tk_radius = det3D.getTrackerRadius();
   Double_t tk_length = det3D.getTrackerHalfLength();
   Double_t tk_Bz     = det3D.getBField();
   Double_t mu_radius = det3D.getDetectorRadius();
   Double_t mu_length = det3D.getDetectorHalfLength();
   TAxis*   etaAxis   = det3D.getCaloAxes().first;
   TAxis*   phiAxis   = det3D.getCaloAxes().second;
   ExRootConfParam branches = confReader->GetParam("TreeWriter::Branch");
   Int_t nBranches = branches.GetSize()/3;
   DelphesBranchElement<TEveTrackList>* tlist;
   DelphesBranchElement<DelphesCaloData>* clist;
   DelphesBranchElement<TEveElementList>* elist;
   for(Int_t b = 0; b<nBranches; ++b) {
     TString input = branches[b*3].GetString();
     TString name = branches[b*3+1].GetString();
     TString className = branches[b*3+2].GetString();
     if(className=="Track") {
       if(input.Contains("eflow",TString::kIgnoreCase) || name.Contains("eflow",TString::kIgnoreCase)) continue; //no eflow
       tlist = new DelphesBranchElement<TEveTrackList>(name,"track",kBlue);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else if(className=="Tower") {
       if(input.Contains("eflow",TString::kIgnoreCase) || name.Contains("eflow",TString::kIgnoreCase)) continue; //no eflow
       clist = new DelphesBranchElement<DelphesCaloData>(name,"tower",kBlack);
       clist->GetContainer()->SetEtaBins(etaAxis);
       clist->GetContainer()->SetPhiBins(phiAxis);
       elements.push_back(clist);
     } else if(className=="Jet") {
       if(input.Contains("GenJetFinder")) {
         elist = new DelphesBranchElement<TEveElementList>(name,"jet",kCyan);
         elist->GetContainer()->SetRnrSelf(false);
         elist->GetContainer()->SetRnrChildren(false);
         elements.push_back(elist);
       } else {
         elements.push_back(new DelphesBranchElement<TEveElementList>(name,"jet",kYellow));
       }
     } else if(className=="Electron") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"track",kRed);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else if(className=="Photon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"photon",kYellow);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., 0.);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     } else if(className=="Muon") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"track",kGreen);
       elements.push_back(tlist);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(mu_radius);
       trkProp->SetMaxZ(mu_length);
     } else if(className=="MissingET") {
       elements.push_back(new DelphesBranchElement<TEveElementList>(name,"vector",kViolet));
     } else if(className=="GenParticle") {
       tlist = new DelphesBranchElement<TEveTrackList>(name,"track",kCyan);
       elements.push_back(tlist);
       tlist->GetContainer()->SetRnrSelf(false);
       tlist->GetContainer()->SetRnrChildren(false);
       TEveTrackPropagator *trkProp = tlist->GetContainer()->GetPropagator();
       trkProp->SetMagField(0., 0., -tk_Bz);
       trkProp->SetMaxR(tk_radius);
       trkProp->SetMaxZ(tk_length);
     }
//TODO one possible simplification could be to add the array to the element class.
     arrays.push_back(gTreeReader->UseBranch(name));
   }
}

void delphes_event_display(const char *configFile, const char *inputFile, Delphes3DGeometry& det3D)
{

   // initialize the application
   TEveManager::Create(kTRUE, "IV");
   TGeoManager* geom = gGeoManager;

   // build the detector
   gRadius = det3D.getTrackerRadius();
   gTotRadius = det3D.getDetectorRadius();
   gHalfLength = det3D.getTrackerHalfLength();
   gBz = det3D.getBField();

   //TODO specific to some classical detector... could use better the det3D
   TGeoVolume* top = det3D.getDetector(false);
   geom->SetTopVolume(top);
   TEveElementList *geometry = new TEveElementList("Geometry");
   TEveGeoTopNode* trk = new TEveGeoTopNode(gGeoManager, top->FindNode("tracker_1"));
   trk->SetVisLevel(6);
   geometry->AddElement(trk);
   TEveGeoTopNode* calo = new TEveGeoTopNode(gGeoManager, top->FindNode("Calorimeter_barrel_1"));
   calo->SetVisLevel(3);
   geometry->AddElement(calo);
   calo = new TEveGeoTopNode(gGeoManager, top->FindNode("Calorimeter_endcap_1"));
   calo->SetVisLevel(3);
   calo->UseNodeTrans();
   geometry->AddElement(calo);
   calo = new TEveGeoTopNode(gGeoManager, top->FindNode("Calorimeter_endcap_2"));
   calo->SetVisLevel(3);
   calo->UseNodeTrans();
   geometry->AddElement(calo);
   TEveGeoTopNode* muon = new TEveGeoTopNode(gGeoManager, top->FindNode("muons_barrel_1"));
   muon->SetVisLevel(4);
   geometry->AddElement(muon);
   muon = new TEveGeoTopNode(gGeoManager, top->FindNode("muons_endcap_1"));
   muon->SetVisLevel(4);
   muon->UseNodeTrans();
   geometry->AddElement(muon);
   muon = new TEveGeoTopNode(gGeoManager, top->FindNode("muons_endcap_2"));
   muon->SetVisLevel(4);
   muon->UseNodeTrans();
   geometry->AddElement(muon);
   //gGeoManager->DefaultColors();

   // Create chain of root trees
   gChain.Add(inputFile);

   // Create object of class ExRootTreeReader
   printf("*** Opening Delphes data file ***\n");
   gTreeReader = new ExRootTreeReader(&gChain);

   // prepare data collections
   readConfig(configFile, det3D, gElements, gArrays);

   // viewers and scenes
   gDelphesDisplay = new DelphesDisplay;
   gEve->AddGlobalElement(geometry);
   gDelphesDisplay->ImportGeomRPhi(geometry);
   gDelphesDisplay->ImportGeomRhoZ(geometry);
   // find the first calo data and use that to initialize the calo display
   for(std::vector<DelphesBranchBase*>::iterator data=gElements.begin();data<gElements.end();++data) {
     if(TString((*data)->GetType())=="tower") {
//TODO: why do I have to split this in two lines??? seems a cint bug?
       DelphesBranchElement<DelphesCaloData>* data_tmp = dynamic_cast<DelphesBranchElement<DelphesCaloData>*>(*data);
       DelphesCaloData* container =  data_tmp->GetContainer();
//       DelphesCaloData* container = dynamic_cast<DelphesBranchElement<DelphesCaloData>*>((*data))->GetContainer();
       assert(container);
       TEveCalo3D *calo3d = new TEveCalo3D(container);
       calo3d->SetBarrelRadius(gRadius);
       calo3d->SetEndCapPos(gHalfLength);
       gEve->AddGlobalElement(calo3d);
       gDelphesDisplay->ImportCaloRPhi(calo3d);
       gDelphesDisplay->ImportCaloRhoZ(calo3d);
       TEveCaloLego *lego = new TEveCaloLego(container);
       lego->InitMainTrans();
       lego->RefMainTrans().SetScale(TMath::TwoPi(), TMath::TwoPi(), TMath::Pi());
       lego->SetAutoRebin(kFALSE);
       lego->Set2DMode(TEveCaloLego::kValSizeOutline);
       gDelphesDisplay->ImportCaloLego(lego);
       break;
     }
   }
   gEve->Redraw3D(kTRUE);
}

//______________________________________________________________________________
void load_event()
{
   // Load event specified in global event_id.
   // The contents of previous event are removed.

   //TODO move this to the status bar ???
   printf("Loading event %d.\n", event_id);

   // clear the previous event
   gEve->GetViewers()->DeleteAnnotations();
   for(std::vector<DelphesBranchBase*>::iterator data=gElements.begin();data<gElements.end();++data) {
     (*data)->Reset();
   }

   // read the new event
   delphes_read();

   // update display
   TEveElement* top = (TEveElement*)gEve->GetCurrentEvent();
   gDelphesDisplay->DestroyEventRPhi();
   gDelphesDisplay->ImportEventRPhi(top);
   gDelphesDisplay->DestroyEventRhoZ();
   gDelphesDisplay->ImportEventRhoZ(top);
   //update_html_summary();
   gEve->Redraw3D(kFALSE, kTRUE);
}

void delphes_read()
{

  // safety
  if(event_id >= gTreeReader->GetEntries() || event_id<0 ) return;

  // Load selected branches with data from specified event
  gTreeReader->ReadEntry(event_id);

  // loop over selected branches, and apply the proper recipe to fill the collections.
  // this is basically to loop on gArrays to fill gElements.

//TODO: one option would be to have templated methods in the element classes. We could simply call "element.fill()"
  std::vector<TClonesArray*>::iterator data = gArrays.begin();
  std::vector<DelphesBranchBase*>::iterator element = gElements.begin();
  std::vector<TClonesArray*>::iterator data_tracks = gArrays.begin();
  std::vector<DelphesBranchBase*>::iterator element_tracks = gElements.begin();
  Int_t nTracks = 0;
  for(; data<gArrays.end() && element<gElements.end(); ++data, ++element) {
    TString type = (*element)->GetType();
    // keep the most generic track collection for the end
    if(type=="track" && TString((*element)->GetClassName())=="Track" && nTracks==0) {
      data_tracks = data;
      element_tracks = element;
      nTracks = (*data_tracks)->GetEntries();
      continue;
    }
    // branch on the element type
    if(type=="tower") delphes_read_towers(*data,*element);
    else if(type=="track" || type=="photon") delphes_read_tracks(*data,*element);
    else if(type=="jet") delphes_read_jets(*data,*element);
    else if(type=="vector") delphes_read_vectors(*data,*element);
  }
  // finish whith what we consider to be the main track collection
  if(nTracks>0) delphes_read_tracks(*data,*element);
}

void delphes_read_towers(TClonesArray* data, DelphesBranchBase* element) {
  DelphesCaloData* container = dynamic_cast<DelphesBranchElement<DelphesCaloData>*>(element)->GetContainer();
  assert(container);
  // Loop over all towers
  TIter itTower(data);
  Tower *tower;
  while((tower = (Tower *) itTower.Next()))
  {
    container->AddTower(tower->Edges[0], tower->Edges[1], tower->Edges[2], tower->Edges[3]);
    container->FillSlice(0, tower->Eem);
    container->FillSlice(1, tower->Ehad);
  }
  container->DataChanged();
}

void delphes_read_tracks(TClonesArray* data, DelphesBranchBase* element) {
  TEveTrackList* container = dynamic_cast<DelphesBranchElement<TEveTrackList>*>(element)->GetContainer();
  assert(container);
  TString className = element->GetClassName();
  TIter itTrack(data);
  Int_t counter = 0;
  TEveTrack *eveTrack;
  TEveTrackPropagator *trkProp = container->GetPropagator();
  if(className=="Track") {
    // Loop over all tracks
    Track *track;
    while((track = (Track *) itTrack.Next())) {
      TParticle pb(track->PID, 1, 0, 0, 0, 0,
                   track->P4().Px(), track->P4().Py(),
                   track->P4().Pz(), track->P4().E(),
                   track->X, track->Y, track->Z, 0.0);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(className=="Electron") {
    // Loop over all electrons
    Electron *electron;
    while((electron = (Electron *) itTrack.Next())) {
      TParticle pb(electron->Charge<0?11:-11, 1, 0, 0, 0, 0,
                   electron->P4().Px(), electron->P4().Py(),
                   electron->P4().Pz(), electron->P4().E(),
                   0., 0., 0., 0.);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
  }
  } else if(className=="Muon") {
    // Loop over all muons
    Muon *muon;
    while((muon = (Muon *) itTrack.Next())) {
      TParticle pb(muon->Charge<0?13:-13, 1, 0, 0, 0, 0,
                   muon->P4().Px(), muon->P4().Py(),
                   muon->P4().Pz(), muon->P4().E(),
                   0., 0., 0., 0.);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  } else if(className=="Photon") {
    // Loop over all photons
    Photon *photon;
    while((photon = (Photon *) itTrack.Next())) {
      TParticle pb(22, 1, 0, 0, 0, 0,
                   photon->P4().Px(), photon->P4().Py(),
                   photon->P4().Pz(), photon->P4().E(),
                   0., 0., 0., 0.);

      eveTrack = new TEveTrack(&pb, counter, trkProp);
      eveTrack->SetName(Form("%s [%d]", pb.GetName(), counter++));
      eveTrack->SetStdTitle();
      eveTrack->SetAttLineAttMarker(container);
      container->AddElement(eveTrack);
      eveTrack->SetLineColor(element->GetColor());
      eveTrack->MakeTrack();
    }
  }
}

void delphes_read_jets(TClonesArray* data, DelphesBranchBase* element) {
  TEveElementList* container = dynamic_cast<DelphesBranchElement<TEveElementList>*>(element)->GetContainer();
  assert(container);
  TIter itJet(data);
  Jet *jet;
  TEveJetCone *eveJetCone;
  // Loop over all jets
  Int_t counter = 0;
  while((jet = (Jet *) itJet.Next()))
  {
    eveJetCone = new TEveJetCone();
    eveJetCone->SetTitle(Form("jet [%d]: Pt=%f, Eta=%f, \nPhi=%f, M=%f",counter,jet->PT, jet->Eta, jet->Phi, jet->Mass));
    eveJetCone->SetName(Form("jet [%d]", counter++));
    eveJetCone->SetMainTransparency(60);
    eveJetCone->SetLineColor(element->GetColor());
    eveJetCone->SetFillColor(element->GetColor());
    eveJetCone->SetCylinder(gRadius - 10, gHalfLength - 10);
    eveJetCone->SetPickable(kTRUE);
    eveJetCone->AddEllipticCone(jet->Eta, jet->Phi, jet->DeltaEta, jet->DeltaPhi);
    container->AddElement(eveJetCone);
  }
}

void delphes_read_vectors(TClonesArray* data, DelphesBranchBase* element) {
  TEveElementList* container = dynamic_cast<DelphesBranchElement<TEveElementList>*>(element)->GetContainer();
  assert(container);
  TIter itMet(data);
  MissingET *MET;
  TEveArrow *eveMet;
  // Missing Et
  Double_t maxPt = 50.;
  // TODO to be changed as we don't have access to maxPt anymore. MET scale could be a general parameter set in GUI
  while((MET = (MissingET*) itMet.Next())) {
    eveMet = new TEveArrow((gRadius * MET->MET/maxPt)*cos(MET->Phi), (gRadius * MET->MET/maxPt)*sin(MET->Phi), 0., 0., 0., 0.);
    eveMet->SetMainColor(element->GetColor());
    eveMet->SetTubeR(0.04);
    eveMet->SetConeR(0.08);
    eveMet->SetConeL(0.10);
    eveMet->SetPickable(kTRUE);
    eveMet->SetName("Missing Et");
    eveMet->SetTitle(Form("Missing Et (%.1f GeV)",MET->MET));
    container->AddElement(eveMet);
  }
}

/******************************************************************************/
// GUI
/******************************************************************************/

//______________________________________________________________________________
// 
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{     
public:
   void Fwd()
   {  
      if (event_id < gTreeReader->GetEntries() - 1) {
         ++event_id;
         load_event();
      } else {
         printf("Already at last event.\n");
      }
   }
   void Bck()
   {
      if (event_id > 0) {
         --event_id;
         load_event();
      } else {
         printf("Already at first event.\n");
      }
   }
};

//______________________________________________________________________________
void make_gui()
{
   // Create minimal GUI for event navigation.
   // TODO: better GUI could be made based on the ch15 of the manual (Writing a GUI)

   // add a tab on the left
   TEveBrowser* browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   // set the main title
   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("Delphes Event Display");
   frmMain->SetCleanup(kDeepCleanup);

   // build the navigation menu
   TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
   {
      TString icondir;
      if(gSystem->Getenv("ROOTSYS"))
        icondir = Form("%s/icons/", gSystem->Getenv("ROOTSYS"));
      if(!gSystem->OpenDirectory(icondir)) 
        icondir = Form("%s/icons/", (const char*)gSystem->GetFromPipe("root-config --etcdir") );
      TGPictureButton* b = 0;
      EvNavHandler    *fh = new EvNavHandler;

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
   }
   frmMain->AddFrame(hf);
   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();
   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);
}

/******************************************************************************/
// MAIN 
/******************************************************************************/

void geometry(const char* filename = "delphes_card_CMS.tcl", const char* ParticlePropagator="ParticlePropagator",
                                                             const char* TrackingEfficiency="ChargedHadronTrackingEfficiency",
                                                             const char* MuonEfficiency="MuonEfficiency",
                                                             const char* Calorimeters="Calorimeter") 
{

   // load the libraries
   gSystem->Load("libGeom");
   gSystem->Load("../libDelphes");
   gSystem->Load("../libDelphesDisplay");

   // create the detector representation
   Delphes3DGeometry det3D(new TGeoManager("delphes", "Delphes geometry"));
   det3D.readFile(filename, ParticlePropagator, TrackingEfficiency, MuonEfficiency, Calorimeters);

   // create the application items
   delphes_event_display("delphes_card_CMS.tcl", "../delphes_output.root", det3D); //TODO root file as input cfg
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

}

