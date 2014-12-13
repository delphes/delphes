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
#include "TClonesArray.h"
#include "TClass.h"
#include <exception>
#include <iostream>
#include "display/DelphesCaloData.h"
#include "TEveElement.h"
#include "TEveTrack.h"

// virtual class to represent objects from a Delphes-tree branch
class DelphesBranchBase
{
  public:
    DelphesBranchBase(const char* name="", TClonesArray* branch=NULL, const enum EColor color=kBlack, Float_t maxPt=50.):name_(name),branch_(branch),color_(color) {}
    virtual ~DelphesBranchBase() {}
    const char* GetName() const { return (const char*)name_; }
    const char* GetType() const { return branch_ ? branch_->GetClass()->GetName() : "None"; }
    virtual const char* GetClassName() = 0;
    enum EColor GetColor() const { return color_; }
    virtual void Reset() = 0;
    virtual void SetTrackingVolume(Float_t r, Float_t l, Float_t Bz=0.) { tkRadius_ = r; tkHalfLength_ = l; tk_Bz_ = Bz; }
    virtual void ReadBranch() = 0;
    virtual std::vector<TLorentzVector> GetVectors() = 0;

  protected:
    TString name_;
    Float_t maxPt_;
    TClonesArray* branch_;
    const enum EColor color_;
    Float_t tkRadius_,tkHalfLength_, tk_Bz_;
};

// concrete implementations. EveContainer can be a TrackList, ElementList or CaloData.
template<typename EveContainer> class DelphesBranchElement: public DelphesBranchBase
{
  public:
    // constructor
    DelphesBranchElement(const char* name="", TClonesArray* branch=NULL, const enum EColor color=kBlack, Float_t maxPt=50.):DelphesBranchBase(name, branch, color, maxPt) {
      throw std::exception();
    }

    // destructor
    virtual ~DelphesBranchElement() { delete data_; }

    // get the container (ElementList, TrackList, or CaloData)
    EveContainer* GetContainer() { return data_; }

    // tracking volume
    virtual void SetTrackingVolume(Float_t r, Float_t l, Float_t Bz=0.) { tkRadius_ = r; tkHalfLength_ = l; tk_Bz_ = Bz; }

    // resets the collection (before moving to the next event)
    virtual void Reset() {};

    // template class name
    virtual const char* GetClassName() { return data_->ClassName(); }

    // read the branch and fill elements for display
    virtual void ReadBranch() {}

    // return the vector for all elements
    virtual std::vector<TLorentzVector> GetVectors() { std::vector<TLorentzVector> v; return v; }

  private:
    EveContainer* data_;
};

#if !defined(__CINT__) && !defined(__CLING__)

// special case for calo towers
template<> DelphesBranchElement<DelphesCaloData>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color, Float_t maxPt);
template<> void DelphesBranchElement<DelphesCaloData>::Reset();
template<> void DelphesBranchElement<DelphesCaloData>::ReadBranch();
template<> std::vector<TLorentzVector> DelphesBranchElement<DelphesCaloData>::GetVectors();

// special case for element lists
template<> DelphesBranchElement<TEveElementList>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color, Float_t maxPt);
template<> void DelphesBranchElement<TEveElementList>::Reset();
template<> void DelphesBranchElement<TEveElementList>::ReadBranch();
template<> std::vector<TLorentzVector> DelphesBranchElement<TEveElementList>::GetVectors();

// special case for track lists
template<> DelphesBranchElement<TEveTrackList>::DelphesBranchElement(const char* name, TClonesArray* branch, const enum EColor color, Float_t maxPt);
template<> void DelphesBranchElement<TEveTrackList>::SetTrackingVolume(Float_t r, Float_t l, Float_t Bz);
template<> void DelphesBranchElement<TEveTrackList>::Reset();
template<> void DelphesBranchElement<TEveTrackList>::ReadBranch();
template<> std::vector<TLorentzVector> DelphesBranchElement<TEveTrackList>::GetVectors();

#endif // CINT, CLING

#endif //DelphesBranchElement_h
