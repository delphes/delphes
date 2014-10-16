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

#include "display/DelphesBranchElement.h"

// special case for calo towers
template<> DelphesBranchElement<DelphesCaloData>::DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
      if(TString(type)=="tower") {
        data_ = new DelphesCaloData(2);
        data_->RefSliceInfo(0).Setup("ECAL", 0.1, kRed);
        data_->RefSliceInfo(1).Setup("HCAL", 0.1, kBlue);
        data_->IncDenyDestroy();
      } else {
        throw std::exception();
      }
    }
template<> void DelphesBranchElement<DelphesCaloData>::Reset() { data_->ClearTowers(); }

// special case for element lists
template<> DelphesBranchElement<TEveElementList>::DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
      if(TString(type)=="vector" || TString(type)=="jet") {
        data_ = new TEveElementList(name);
        data_->SetMainColor(color_);
      } else {
        throw std::exception();
      }
    }
template<> void DelphesBranchElement<TEveElementList>::Reset() { data_->DestroyElements(); }

//TODO: does the type really make sense?
// special case for track lists
template<> DelphesBranchElement<TEveTrackList>::DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
      if(TString(type)=="track") {
        data_ = new TEveTrackList(name);
        data_->SetMainColor(color_);
        data_->SetMarkerColor(color_);
        data_->SetMarkerStyle(kCircle);
        data_->SetMarkerSize(0.5);
      } else if(TString(type)=="photon") {
        data_ = new TEveTrackList(name);
        data_->SetMainColor(color_);
        data_->SetMarkerColor(color_);
        data_->SetMarkerStyle(kCircle);
        data_->SetMarkerSize(0.5);
      } else {
        data_ = new TEveTrackList(name);
        data_->SetMainColor(color_);
        data_->SetMarkerColor(color_);
        data_->SetMarkerStyle(kCircle);
        data_->SetMarkerSize(0.5);
        //throw std::exception();
      }
    }
template<> void DelphesBranchElement<TEveTrackList>::Reset() { data_->DestroyElements(); }

