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

#ifndef Delphes3DGeometry_h
#define Delphes3DGeometry_h

#include <set>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>
#include <sstream>
#include "TAxis.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"

// TODO: asymmetric detector

class Delphes3DGeometry {
   public:
     Delphes3DGeometry(TGeoManager *geom = NULL, bool transp = false);
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
     std::pair<Double_t, Double_t> addCalorimeter(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength, std::set< std::pair<Double_t, Int_t> >& caloBinning);
     std::pair<Double_t, Double_t> addMuonDets(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength);
     void addCaloTowers(TGeoVolume *top, const char* name, Double_t innerBarrelRadius, Double_t innerBarrelLength, std::set< std::pair<Double_t, Int_t> >& caloBinning);

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
     std::map<std::string, std::set< std::pair<Double_t, Int_t> > > caloBinning_;
     
};

#endif
