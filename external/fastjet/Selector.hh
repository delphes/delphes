#ifndef __FASTJET_SELECTOR_HH__
#define __FASTJET_SELECTOR_HH__

//FJSTARTHEADER
// $Id: Selector.hh 3711 2014-09-29 13:54:51Z salam $
//
// Copyright (c) 2009-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#include "fastjet/PseudoJet.hh"
#ifndef __FJCORE__
#include "fastjet/RangeDefinition.hh"  // for initialisation from a RangeDefinition
#endif  // __FJCORE__
#include <limits>
#include <cmath>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// @ingroup selectors
/// \class Selector
/// Class that encodes information about cuts and other selection
/// criteria that can be applied to PseudoJet(s).
///
class Selector;
//----------------------------------------------------------------------

/// @ingroup selectors
/// \class SelectorWorker
/// default selector worker is an abstract virtual base class
///
/// The Selector class is only an interface, it is the SelectorWorker
/// that really does the work. To implement various selectors, one
/// thus has to overload this class.
class SelectorWorker {
public:
  //----------------------------------------------------------
  // fundamental info
  //----------------------------------------------------------
  /// default dtor
  virtual ~SelectorWorker() {}

  //----------------------------------------------------------
  // basic operations for checking what gets selected
  //----------------------------------------------------------

  /// returns true if a given object passes the selection criterion,
  /// and is the main function that needs to be overloaded by derived
  /// workers. 
  ///
  /// NB: this function is used only if applies_jet_by_jet() returns
  /// true. If it does not, then derived classes are expected to
  /// (re)implement the terminator function()
  virtual bool pass(const PseudoJet & jet) const = 0;

  /// For each jet that does not pass the cuts, this routine sets the 
  /// pointer to 0.
  ///
  /// It does not assume that the PseudoJet* passed as argument are not NULL
  virtual void terminator(std::vector<const PseudoJet *> & jets) const {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (jets[i] && !pass(*jets[i])) jets[i] = NULL;
    }
  }

  /// returns true if this can be applied jet by jet
  virtual bool applies_jet_by_jet() const {return true;}

  /// returns a description of the worker
  virtual std::string description() const {return "missing description";}


  //----------------------------------------------------------
  // operations for dealing with reference jets
  //----------------------------------------------------------

  /// returns true if the worker is defined with respect to a reference jet
  virtual bool takes_reference() const { return false;}

  /// sets the reference jet for the selector
  /// NB: "reference" is commented to avoid unused-variable compiler warnings
  virtual void set_reference(const PseudoJet & /*reference*/){
    throw Error("set_reference(...) cannot be used for a selector worker that does not take a reference");
  }

  /// return a copy of the current object.
  ///
  /// This function is only called for objects that take a reference and need
  /// not be reimplemented otherwise.
  virtual SelectorWorker* copy(){ 
    throw Error("this SelectorWorker has nothing to copy");
  }

  //----------------------------------------------------------
  // operations for area and extent
  //----------------------------------------------------------

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    rapmax = std::numeric_limits<double>::infinity();
    rapmin = -rapmax; 
  }

  /// check if it is a geometric selector (i.e. only puts constraints
  /// on rapidity and azimuthal angle)
  virtual bool is_geometric() const { return false;}

  /// check if it has a finite area
  virtual bool has_finite_area() const;

  /// check if it has an analytically computable area
  virtual bool has_known_area() const { return false;}

  /// if it has a computable area, return it
  virtual double known_area() const{
    throw Error("this selector has no computable area");
  }

};

//----------------------------------------------------------------------
// class Selector
//
// Class that encodes information about cuts that 
class Selector{
public:
  /// default constructor produces a Selector whose action is undefined
  /// (any attempt to use it will lead to an error)
  Selector() {}

  /// constructor that causes the Selector to use the supplied worker
  ///
  /// Note that the Selector takes ownership of the pointer to the
  /// worker (and so will delete automatically when appropriate).
  Selector(SelectorWorker * worker_in) {_worker.reset(worker_in);}

#ifndef __FJCORE__
  /// ctor from a RangeDefinition
  ///
  /// This is provided for backward compatibility and will be removed in
  /// a future major release of FastJet
  ///
  /// Watch out that the Selector will only hold a pointer to the
  /// range so the selector will crash if one tries to use it after
  /// the range has gone out of scope. We thus strongly advise against
  /// the direct use of this constructor.
  Selector(const RangeDefinition &range);
#endif  // __FJCORE__

  /// dummy virtual dtor
  virtual ~Selector(){}

  /// return true if the jet passes the selection
  bool pass(const PseudoJet & jet) const {
    if (!validated_worker()->applies_jet_by_jet()) {
      throw Error("Cannot apply this selector to an individual jet");
    }
    return _worker->pass(jet);
  }

  /// an operator way of knowing whether a given jet passes the selection or not
  bool operator()(const PseudoJet & jet) const {
    return pass(jet);
  }

  /// Return a count of the objects that pass the selection.
  ///
  /// This will often be more efficient that getting the vector of objects that
  /// passes and then evaluating the size of the vector
  unsigned int count(const std::vector<PseudoJet> & jets) const;

  /// Return the 4-vector sum of the objects that pass the selection.
  ///
  /// This will often be more efficient that getting the vector of objects that
  /// passes and then evaluating the size of the vector
  PseudoJet sum(const std::vector<PseudoJet> & jets) const;

  /// Return the scalar pt sum of the objects that pass the selection.
  ///
  /// This will often be more efficient that getting the vector of objects that
  /// passes and then evaluating the size of the vector
  double scalar_pt_sum(const std::vector<PseudoJet> & jets) const;

  /// sift the input jets into two vectors -- those that pass the selector
  /// and those that do not
  void sift(const std::vector<PseudoJet> & jets,
		  std::vector<PseudoJet> & jets_that_pass,
		  std::vector<PseudoJet> & jets_that_fail) const;

  /// returns true if this can be applied jet by jet
  bool applies_jet_by_jet() const {
    return validated_worker()->applies_jet_by_jet();
  }

  /// returns a vector with the jets that pass the selection
  std::vector<PseudoJet> operator()(const std::vector<PseudoJet> & jets) const;

  /// For each jet that does not pass the cuts, this routine sets the 
  /// pointer to 0. 
  ///
  /// It is legitimate for some (or all) of the pointers that are
  /// passed to already be NULL.
  virtual void nullify_non_selected(std::vector<const PseudoJet *> & jets) const {
    validated_worker()->terminator(jets);
  }

  /// returns the rapidity range for which it may return "true"
  void get_rapidity_extent(double &rapmin, double &rapmax) const {
    return validated_worker()->get_rapidity_extent(rapmin, rapmax);
  }

  /// returns a textual description of the selector
  std::string description() const {
    return validated_worker()->description();
  }

  /// returns true if it is a geometric selector (i.e. one that only puts
  /// constraints on rapidities and azimuthal angles)
  bool is_geometric() const{
    return validated_worker()->is_geometric();
  }

  /// returns true if it has a meaningful and finite area (i.e. the
  /// Selector has the property that is_geometric() returns true and
  /// the rapidity extent is finite).
  bool has_finite_area() const{
    return validated_worker()->has_finite_area();
  }

#ifndef __FJCORE__
  /// returns the rapidity-phi area associated with the Selector
  /// (throws InvalidArea if the area does not make sense).
  ///
  /// If the result is not known analytically, the area will be
  /// estimated using a pseudo Monte Carlo method (as for jet areas),
  /// using the default ghost area from the GhostedAreaSpec class
  /// (0.01). The Monte Carlo estimate involves a time penalty
  /// proportional to the ratio of the rapidity extent of the Selector
  /// divided by the ghost area.
  double area() const;

  /// returns the rapidity-phi area associated with the Selector
  /// (throws InvalidArea if the area does not make sense).
  ///
  /// The behaviour is the as with the area() call, but with the
  /// ability to additionally specify the ghost area to be used in the
  /// case of a Monte Carlo area evaluation.
  ///
  double area(double ghost_area) const;
#endif  // __FJCORE__

  /// returns a (reference to) the underlying worker's shared pointer
  const SharedPtr<SelectorWorker> & worker() const {return _worker;}

  /// returns a worker if there is a valid one, otherwise throws an InvalidWorker error
  const SelectorWorker* validated_worker() const {
    const SelectorWorker* worker_ptr = _worker.get();
    if (worker_ptr == 0) throw InvalidWorker();
    return worker_ptr;
  }

  /// returns true if this can be applied jet by jet
  bool takes_reference() const {
    return validated_worker()->takes_reference();
  }

  /// set the reference jet for this Selector
  const Selector & set_reference(const PseudoJet &reference){

    // if the worker does not take a reference jet, do nothing 
    if (! validated_worker()->takes_reference()){
      return *this;
    }
    
    // since this is a non-const operation, make sure we have a
    // correct behaviour with respect to shared workers
    _copy_worker_if_needed();

    _worker->set_reference(reference);
    return *this;
  }

  /// class that gets thrown when a Selector is applied despite it not
  /// having a valid underlying worker.
  class InvalidWorker : public Error {
  public:
    InvalidWorker() : Error("Attempt to use Selector with no valid underlying worker") {}
  };

  /// class that gets thrown when the area is requested from a Selector for which
  /// the area is not meaningful
  class InvalidArea : public Error {
  public:
    InvalidArea() : Error("Attempt to obtain area from Selector for which this is not meaningful") {}
  };

  // some operators (applying directly on a Selector)
  //----------------------------------------------------------------------
  /// For 2 Selectors a and b, a &= b is eauivalent to a = a && b;
  Selector & operator &=(const Selector & b);

  /// For 2 Selectors a and b, a |= b is eauivalent to a = a || b;
  Selector & operator |=(const Selector & b);


protected:
  /// Helper for copying selector workers if needed
  ///
  /// The following is needed if we want to modify a selectors that
  /// shares a worker with another selector. In that case, we need to
  /// get another copy of the worker to avoid interferences
  ///
  /// Note that any non-const operation has to call this to behave
  /// correctly w.r.t shared workers!
  void _copy_worker_if_needed(){
    // do nothing if there's a sinlge user of the worker
    if (_worker.unique()) return;

    // call the worker's copy
    //std::cout << "will make a copy of " << description() << std::endl;
    _worker.reset(_worker->copy());
  }

private:
  SharedPtr<SelectorWorker> _worker; ///< the underlying worker
};


//----------------------------------------------------------------------
// a list of specific selectors
//----------------------------------------------------------------------

/// \addtogroup selectors
/// @{


// fundamental selectors
//----------------------------------------------------------------------

// "identity" selector that lets everything pass
Selector SelectorIdentity();

// logical operations
//----------------------------------------------------------------------

/// logical not applied on a selector
///
/// This will keep objects that do not pass the 's' selector
Selector operator!(const Selector & s);

/// logical or between two selectors
///
/// this will keep the objects that are selected by s1 or s2
Selector operator ||(const Selector & s1, const Selector & s2);


/// logical and between two selectors
///
/// this will keep the objects that are selected by both s1 and s2
/// 
/// watch out: for both s1 and s2, the selection is applied on the
///   original list of objects. For successive applications of two
///   selectors (convolution/multiplication) see the operator *
Selector operator&&(const Selector & s1, const Selector & s2);

/// successive application of 2 selectors
///
/// Apply the selector s2, then the selector s1.
///
/// watch out: the operator * acts like an operator product i.e. does
///   not commute. The order of its arguments is therefore important.
///   Whenever they commute (in particluar, when they apply jet by
///   jet), this would have the same effect as the logical &&.
Selector operator*(const Selector & s1, const Selector & s2);


// selection with kinematic cuts
//----------------------------------------------------------------------
Selector SelectorPtMin(double ptmin);                    ///< select objects with pt >= ptmin
Selector SelectorPtMax(double ptmax);                    ///< select objects with pt <= ptmax
Selector SelectorPtRange(double ptmin, double ptmax);    ///< select objects with ptmin <= pt <= ptmax

Selector SelectorEtMin(double Etmin);                    ///< select objects with Et >= Etmin
Selector SelectorEtMax(double Etmax);                    ///< select objects with Et <= Etmax
Selector SelectorEtRange(double Etmin, double Etmax);    ///< select objects with Etmin <= Et <= Etmax

Selector SelectorEMin(double Emin);                      ///< select objects with E >= Emin
Selector SelectorEMax(double Emax);                      ///< select objects with E <= Emax
Selector SelectorERange(double Emin, double Emax);       ///< select objects with Emin <= E <= Emax

Selector SelectorMassMin(double Mmin);                      ///< select objects with Mass >= Mmin
Selector SelectorMassMax(double Mmax);                      ///< select objects with Mass <= Mmax
Selector SelectorMassRange(double Mmin, double Mmax);       ///< select objects with Mmin <= Mass <= Mmax

Selector SelectorRapMin(double rapmin);                  ///< select objects with rap >= rapmin
Selector SelectorRapMax(double rapmax);                  ///< select objects with rap <= rapmax
Selector SelectorRapRange(double rapmin, double rapmax); ///< select objects with rapmin <= rap <= rapmax

Selector SelectorAbsRapMin(double absrapmin);                     ///< select objects with |rap| >= absrapmin
Selector SelectorAbsRapMax(double absrapmax);                     ///< select objects with |rap| <= absrapmax
Selector SelectorAbsRapRange(double absrapmin, double absrapmax); ///< select objects with absrapmin <= |rap| <= absrapmax

Selector SelectorEtaMin(double etamin);                  ///< select objects with eta >= etamin
Selector SelectorEtaMax(double etamax);                  ///< select objects with eta <= etamax
Selector SelectorEtaRange(double etamin, double etamax); ///< select objects with etamin <= eta <= etamax

Selector SelectorAbsEtaMin(double absetamin);                     ///< select objects with |eta| >= absetamin
Selector SelectorAbsEtaMax(double absetamax);                     ///< select objects with |eta| <= absetamax
Selector SelectorAbsEtaRange(double absetamin, double absetamax); ///< select objects with absetamin <= |eta| <= absetamax

Selector SelectorPhiRange(double phimin, double phimax); ///< select objects with phimin <= phi <= phimax

/// select objects with rapmin <= rap <= rapmax  &&  phimin <= phi <= phimax
///
/// Note that this is essentially a combination of SelectorRapRange
/// and SelectorPhiRange. We provide it as a Selector on its own in
/// order to use the known area (which would otherwise be lost by the &&
/// operator)
Selector SelectorRapPhiRange(double rapmin, double rapmax, double phimin, double phimax);

/// select the n hardest objects 
Selector SelectorNHardest(unsigned int n); 


// Selectors that take (require) a reference jet.
//----------------------------------------------------------------------

/// select objets within a distance 'radius' from the location of the
/// reference jet, set by Selector::set_reference(...)
Selector SelectorCircle(const double radius); 

/// select objets with distance from the reference jet is between 'radius_in'
/// and 'radius_out'; the reference jet is set by Selector::set_reference(...)
Selector SelectorDoughnut(const double radius_in, const double radius_out); 

/// select objets within a rapidity distance 'half_width' from the
/// location of the reference jet, set by Selector::set_reference(...)
Selector SelectorStrip(const double half_width);

/// select objets within rapidity distance 'half_rap_width' from the
/// reference jet and azimuthal-angle distance within 'half_phi_width'; the
/// reference jet is set by Selector::set_reference(...)
Selector SelectorRectangle(const double half_rap_width, const double half_phi_width);


/// select objects that carry at least a fraction "fraction" of the
/// reference jet. The reference jet must have been set with
/// Selector::set_reference(...)
Selector SelectorPtFractionMin(double fraction);


// additional (mostly helper) selectors
//----------------------------------------------------------------------

/// select PseudoJet with 0 momentum
Selector SelectorIsZero();

#ifndef __FJCORE__
/// select objects that are (or are only made of) ghosts.
/// PseudoJets for which has_area() are considered non-pure-ghost.
Selector SelectorIsPureGhost();
#endif  // __FJCORE__

/// @}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_SELECTOR_HH__

