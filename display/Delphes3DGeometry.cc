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

#include "display/Delphes3DGeometry.h"
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cassert>
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoCone.h"
#include "TGeoArb8.h"
#include "external/ExRootAnalysis/ExRootConfReader.h"
#include "classes/DelphesClasses.h"
#include "TF2.h"
#include "TH1F.h"
#include "TMath.h"

using namespace std;

Delphes3DGeometry::Delphes3DGeometry(TGeoManager *geom, bool transp) {

   //--- the geometry manager
   geom_ = geom==NULL? gGeoManager : geom;
   //gGeoManager->DefaultColors();

   //--- define some materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7); // placeholder
   if(transp) {
     matVacuum->SetTransparency(85);
     matAl->SetTransparency(85);
   }

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

