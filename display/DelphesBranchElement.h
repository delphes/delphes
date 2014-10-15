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

#ifndef DelphesBranchElement_h
#define DelphesBranchElement_h

#include "TColor.h"
#include "TString.h"
#include <exception>
#include "display/DelphesCaloData.h"
#include "TEveElement.h"
#include "TEveTrack.h"

// virtual class to represent objects from a Delphes-tree branch
class DelphesBranchBase
{
  public:
    DelphesBranchBase():color_(kBlack) {}
    DelphesBranchBase(const char* name, const char*type, const enum EColor color):name_(name),type_(type),color_(color) {}
    virtual ~DelphesBranchBase() {};
    const char* GetName() const { return (const char*)name_; }
    const char* GetType() const { return (const char*)type_; }
    enum EColor GetColor() const { return color_; }
    virtual const char* GetClassName() = 0;
    virtual void Reset() = 0;

  protected:
    TString name_;
    TString type_; // needed for parsing the branch later on
    const enum EColor color_;
};

// concrete implementations. EveContainer can be a TrackList, ElementList or CaloData.
template<typename EveContainer> class DelphesBranchElement: public DelphesBranchBase
{
  public:
    // constructor
    DelphesBranchElement():DelphesBranchBase() {}
    DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {}

    // destructor
    virtual ~DelphesBranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    EveContainer* GetContainer() { return data_; }

    // resets the collection (before moving to the next event)
    virtual void Reset() = 0;

    // template class name
    virtual const char* GetClassName() { return data_->ClassName(); }

  private:
    EveContainer* data_;
};

// special case for calo towers
template<> class DelphesBranchElement<DelphesCaloData>: public DelphesBranchBase
{
  public:
    // constructor
    DelphesBranchElement():DelphesBranchBase() {}
    DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
      if(TString(type)=="tower") {
        data_ = new DelphesCaloData(2);
        data_->RefSliceInfo(0).Setup("ECAL", 0.1, kRed);
        data_->RefSliceInfo(1).Setup("HCAL", 0.1, kBlue);
        data_->IncDenyDestroy();
      } else {
        throw std::exception();
      }
    }

    // destructor
    virtual ~DelphesBranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    DelphesCaloData* GetContainer() { return data_; }

    // resets the collection (before moving to the next event)
    virtual void Reset() { data_->ClearTowers(); }

    // template class name
    virtual const char* GetClassName() { return data_->ClassName(); }

  private:
    DelphesCaloData* data_;
};
//template<> DelphesBranchElement<DelphesCaloData>::DelphesBranchElement(const char* name, const char*type, const enum EColor color);
//template<> void DelphesBranchElement<DelphesCaloData>::Reset();

// special case for element lists
template<> class DelphesBranchElement<TEveElementList>: public DelphesBranchBase
{
  public:
    // constructor
    DelphesBranchElement():DelphesBranchBase() {}
    DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
      if(TString(type)=="vector" || TString(type)=="jet") {
        data_ = new TEveElementList(name);
        data_->SetMainColor(color_);
      } else {
        throw std::exception();
      }
    }

    // destructor
    virtual ~DelphesBranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    TEveElementList* GetContainer() { return data_; }

    // resets the collection (before moving to the next event)
    virtual void Reset() { data_->DestroyElements(); }

    // template class name
    virtual const char* GetClassName() { return data_->ClassName(); }

  private:
    TEveElementList* data_;
};
//template<> DelphesBranchElement<TEveElementList>::DelphesBranchElement(const char* name, const char*type, const enum EColor color);
//template<> void DelphesBranchElement<TEveElementList>::Reset();

// special case for track lists
template<> class DelphesBranchElement<TEveTrackList>: public DelphesBranchBase
{
  public:
    // constructor
    DelphesBranchElement():DelphesBranchBase() {}
    DelphesBranchElement(const char* name, const char*type, const enum EColor color):DelphesBranchBase(name, type, color) {
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
        throw std::exception();
      }
    }

    // destructor
    virtual ~DelphesBranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    TEveTrackList* GetContainer() { return data_; }

    // resets the collection (before moving to the next event)
    virtual void Reset() { data_->DestroyElements(); }

    // template class name
    virtual const char* GetClassName() { return data_->ClassName(); }

  private:
    TEveTrackList* data_;
};
//template<> DelphesBranchElement<TEveTrackList>::DelphesBranchElement(const char* name, const char*type, const enum EColor color);
//template<> void DelphesBranchElement<TEveTrackList>::Reset();

#endif //DelphesBranchElement_h
