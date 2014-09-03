//FJSTARTHEADER
// $Id: Selector.cc 3504 2014-08-01 06:07:54Z soyez $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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


#include <sstream>
#include <algorithm>
#include "fastjet/Selector.hh"
#ifndef __FJCORE__
#include "fastjet/GhostedAreaSpec.hh"  // for area support
#endif  // __FJCORE__

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
// implementations of some of the more complex bits of Selector
//----------------------------------------------------------------------

// implementation of the operator() acting on a vector of jets
std::vector<PseudoJet> Selector::operator()(const std::vector<PseudoJet> & jets) const {
  std::vector<PseudoJet> result;
  const SelectorWorker * worker_local = validated_worker();
  if (worker_local->applies_jet_by_jet()) {
    //if (false) {
    // for workers that apply jet by jet, this is more efficient
    for (std::vector<PseudoJet>::const_iterator jet = jets.begin(); 
         jet != jets.end(); jet++) {
      if (worker_local->pass(*jet)) result.push_back(*jet);
    }
  } else {
    // for workers that can only be applied to entire vectors,
    // go through the following
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) result.push_back(jets[i]);
    }
  }
  return result;
}


//----------------------------------------------------------------------
// count the number of jets that pass the cuts
unsigned int Selector::count(const std::vector<PseudoJet> & jets) const {
  unsigned n = 0;
  const SelectorWorker * worker_local = validated_worker();
  
  // separate strategies according to whether the worker applies jet by jet
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) n++;
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) n++;
    }
  }

  return n;
}

//----------------------------------------------------------------------
// sum the momenta of the jets that pass the cuts
PseudoJet Selector::sum(const std::vector<PseudoJet> & jets) const {
  PseudoJet this_sum(0,0,0,0);
  const SelectorWorker * worker_local = validated_worker();
  
  // separate strategies according to whether the worker applies jet by jet
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) this_sum += jets[i];
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) this_sum += jets[i];
    }
  }

  return this_sum;
}

//----------------------------------------------------------------------
// sum the (scalar) pt of the jets that pass the cuts
double Selector::scalar_pt_sum(const std::vector<PseudoJet> & jets) const {
  double this_sum = 0.0;
  const SelectorWorker * worker_local = validated_worker();
  
  // separate strategies according to whether the worker applies jet by jet
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) this_sum += jets[i].pt();
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) this_sum += jets[i].pt();
    }
  }

  return this_sum;
}


//----------------------------------------------------------------------
// sift the input jets into two vectors -- those that pass the selector
// and those that do not
void Selector::sift(const std::vector<PseudoJet> & jets,
		    std::vector<PseudoJet> & jets_that_pass,
		    std::vector<PseudoJet> & jets_that_fail
		    ) const {
  const SelectorWorker * worker_local = validated_worker();
  
  jets_that_pass.clear();
  jets_that_fail.clear();
  
  // separate strategies according to whether the worker applies jet by jet
  if (worker_local->applies_jet_by_jet()) {
    for (unsigned i = 0; i < jets.size(); i++) {
      if (worker_local->pass(jets[i])) {
	jets_that_pass.push_back(jets[i]);
      } else {
	jets_that_fail.push_back(jets[i]);
      }
    }
  } else {
    std::vector<const PseudoJet *> jetptrs(jets.size());
    for (unsigned i = 0; i < jets.size(); i++) {
      jetptrs[i] = & jets[i];
    }
    worker_local->terminator(jetptrs);
    for (unsigned i = 0; i < jetptrs.size(); i++) {
      if (jetptrs[i]) {
	jets_that_pass.push_back(jets[i]);
      } else {
	jets_that_fail.push_back(jets[i]);
      }
    }
  }
}

#ifndef __FJCORE__
// area using default ghost area
double Selector::area() const{
  return area(gas::def_ghost_area);
}

// implementation of the Selector's area function
double Selector::area(double ghost_area) const{
  if (! is_geometric()) throw InvalidArea();
  
  // has area will already check we've got a valid worker
  if (_worker->has_known_area()) return _worker->known_area();
  
  // generate a set of "ghosts"
  double rapmin, rapmax;
  get_rapidity_extent(rapmin, rapmax);
  GhostedAreaSpec ghost_spec(rapmin, rapmax, 1, ghost_area);
  std::vector<PseudoJet> ghosts;
  ghost_spec.add_ghosts(ghosts);
  
  // check what passes the selection
  return ghost_spec.ghost_area() * ((*this)(ghosts)).size();
}
#endif  // __FJCORE__


//----------------------------------------------------------------------
// implementations of some of the more complex bits of SelectorWorker
//----------------------------------------------------------------------
// check if it has a finite area
bool SelectorWorker::has_finite_area() const { 
  if (! is_geometric()) return false;
  double rapmin, rapmax;
  get_rapidity_extent(rapmin, rapmax);
  return (rapmax != std::numeric_limits<double>::infinity())
    &&  (-rapmin != std::numeric_limits<double>::infinity());
}



//----------------------------------------------------------------------
// very basic set of selectors (at the moment just the identity!)
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// helper for selecting the n hardest jets
class SW_Identity : public SelectorWorker {
public:
  /// ctor with specification of the number of objects to keep
  SW_Identity(){}

  /// just let everything pass
  virtual bool pass(const PseudoJet &) const {
    return true;
  }

  /// For each jet that does not pass the cuts, this routine sets the 
  /// pointer to 0. 
  virtual void terminator(vector<const PseudoJet *> &) const {
    // everything passes, hence nothing to nullify
    return;
  }
  
  /// returns a description of the worker
  virtual string description() const { return "Identity";}

  /// strictly speaking, this is geometric
  virtual bool is_geometric() const { return true;}
};


// returns an "identity" selector that lets everything pass
Selector SelectorIdentity() {
  return Selector(new SW_Identity);
}


//----------------------------------------------------------------------
// selector and workers for operators
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// helper for combining selectors with a logical not
class SW_Not : public SelectorWorker {
public:
  /// ctor
  SW_Not(const Selector & s) : _s(s) {}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_Not(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure that the "pass" can be applied on both selectors
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    
    return ! _s.pass(jet);
  } 

  /// returns true if this can be applied jet by jet
  virtual bool applies_jet_by_jet() const {return _s.applies_jet_by_jet();}

  /// select the jets in the list that pass both selectors
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    // if we can apply the selector jet-by-jet, call the base selector
    // that does exactly that
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }

    // check the effect of the selector we want to negate
    vector<const PseudoJet *> s_jets = jets;
    _s.worker()->terminator(s_jets);

    // now apply the negation: all the jets that pass the base
    // selector (i.e. are not NULL) have to be set to NULL
    for (unsigned int i=0; i<s_jets.size(); i++){
      if (s_jets[i]) jets[i] = NULL;
    }
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "!(" << _s.description() << ")";
    return ostr.str();
  }

  /// is geometric if the underlying selector is
  virtual bool is_geometric() const { return _s.is_geometric();}

  /// returns true if the worker can be set_referenced
  virtual bool takes_reference() const { return _s.takes_reference();}

  /// set the reference jet for this selector
  virtual void set_reference(const PseudoJet &ref) { _s.set_reference(ref);}

protected:
  Selector _s;
};


// logical not applied on a selector
Selector operator!(const Selector & s) {
  return Selector(new SW_Not(s));
}


//----------------------------------------------------------------------
/// Base class for binary operators
class SW_BinaryOperator: public SelectorWorker {
public:
  /// ctor
  SW_BinaryOperator(const Selector & s1, const Selector & s2) : _s1(s1), _s2(s2) {
    // stores info for more efficient access to the selector's properties

    // we can apply jet by jet only if this is the case for both sub-selectors
    _applies_jet_by_jet = _s1.applies_jet_by_jet() && _s2.applies_jet_by_jet();

    // the selector takes a reference if either of the sub-selectors does
    _takes_reference = _s1.takes_reference() || _s2.takes_reference();

    // we have a well-defined area provided the two objects have one
    _is_geometric = _s1.is_geometric() && _s2.is_geometric();
  }

  /// returns true if this can be applied jet by jet
  virtual bool applies_jet_by_jet() const {return _applies_jet_by_jet;}

  /// returns true if this takes a reference jet
  virtual bool takes_reference() const{ 
    return _takes_reference;
  }

  /// sets the reference jet
  virtual void set_reference(const PseudoJet &centre){
    _s1.set_reference(centre);
    _s2.set_reference(centre);
  }

  /// check if it has a finite area
  virtual bool is_geometric() const { return _is_geometric;} 

protected:
  Selector _s1, _s2;
  bool _applies_jet_by_jet;
  bool _takes_reference;
  bool _is_geometric;
};



//----------------------------------------------------------------------
/// helper for combining selectors with a logical and
class SW_And: public SW_BinaryOperator {
public:
  /// ctor
  SW_And(const Selector & s1, const Selector & s2) : SW_BinaryOperator(s1,s2){}

  /// return a copy of this
  virtual SelectorWorker* copy(){ return new SW_And(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure that the "pass" can be applied on both selectors
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    
    return _s1.pass(jet) && _s2.pass(jet);
  }

  /// select the jets in the list that pass both selectors
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    // if we can apply the selector jet-by-jet, call the base selector
    // that does exactly that
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }

    // check the effect of the first selector
    vector<const PseudoJet *> s1_jets = jets;
    _s1.worker()->terminator(s1_jets);

    // apply the second
    _s2.worker()->terminator(jets);

    // terminate the jets that wiuld be terminated by _s1
    for (unsigned int i=0; i<jets.size(); i++){
      if (! s1_jets[i]) jets[i] = NULL;
    }
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    double s1min, s1max, s2min, s2max;
    _s1.get_rapidity_extent(s1min, s1max);
    _s2.get_rapidity_extent(s2min, s2max);
    rapmax = min(s1max, s2max);
    rapmin = max(s1min, s2min);
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " && " << _s2.description() << ")";
    return ostr.str();
  }
};


// logical and between two selectors
Selector operator&&(const Selector & s1, const Selector & s2) {
  return Selector(new SW_And(s1,s2));
}



//----------------------------------------------------------------------
/// helper for combining selectors with a logical or
class SW_Or: public SW_BinaryOperator {
public:
  /// ctor
  SW_Or(const Selector & s1, const Selector & s2) : SW_BinaryOperator(s1,s2) {}

  /// return a copy of this
  virtual SelectorWorker* copy(){ return new SW_Or(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure that the "pass" can be applied on both selectors
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    
    return _s1.pass(jet) || _s2.pass(jet);
  }

  /// returns true if this can be applied jet by jet
  virtual bool applies_jet_by_jet() const {
    // watch out, even though it's the "OR" selector, to be applied jet
    // by jet, both the base selectors need to be jet-by-jet-applicable,
    // so the use of a && in the line below
    return _s1.applies_jet_by_jet() && _s2.applies_jet_by_jet();
  }

  /// select the jets in the list that pass both selectors
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    // if we can apply the selector jet-by-jet, call the base selector
    // that does exactly that
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }

    // check the effect of the first selector
    vector<const PseudoJet *> s1_jets = jets;
    _s1.worker()->terminator(s1_jets);

    // apply the second
    _s2.worker()->terminator(jets);

    // resurrect any jet that has been terminated by the second one
    // and not by the first one
    for (unsigned int i=0; i<jets.size(); i++){
      if (s1_jets[i]) jets[i] = s1_jets[i];
    }
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " || " << _s2.description() << ")";
    return ostr.str();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const {
    double s1min, s1max, s2min, s2max;
    _s1.get_rapidity_extent(s1min, s1max);
    _s2.get_rapidity_extent(s2min, s2max);
    rapmax = max(s1max, s2max);
    rapmin = min(s1min, s2min);
  }
};


// logical or between two selectors
Selector operator ||(const Selector & s1, const Selector & s2) {
  return Selector(new SW_Or(s1,s2));
}

//----------------------------------------------------------------------
/// helper for multiplying two selectors (in an operator-like way)
class SW_Mult: public SW_And {
public:
  /// ctor
  SW_Mult(const Selector & s1, const Selector & s2) : SW_And(s1,s2) {}

  /// return a copy of this
  virtual SelectorWorker* copy(){ return new SW_Mult(*this);}

  /// select the jets in the list that pass both selectors
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    // if we can apply the selector jet-by-jet, call the base selector
    // that does exactly that
    if (applies_jet_by_jet()){
      SelectorWorker::terminator(jets);
      return;
    }

    // first apply _s2
    _s2.worker()->terminator(jets);

    // then apply _s1
    _s1.worker()->terminator(jets);
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "(" << _s1.description() << " * " << _s2.description() << ")";
    return ostr.str();
  }
};


// logical and between two selectors
Selector operator*(const Selector & s1, const Selector & s2) {
  return Selector(new SW_Mult(s1,s2));
}


//----------------------------------------------------------------------
// selector and workers for kinematic cuts
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// a series of basic classes that allow easy implementations of
// min, max and ranges on a quantity-to-be-defined

// generic holder for a quantity
class QuantityBase{
public:
  QuantityBase(double q) : _q(q){}
  virtual ~QuantityBase(){}
  virtual double operator()(const PseudoJet & jet ) const =0;
  virtual string description() const =0;
  virtual bool is_geometric() const { return false;}
  virtual double comparison_value() const {return _q;}
  virtual double description_value() const {return comparison_value();}
protected:
  double _q;
};  

// generic holder for a squared quantity
class QuantitySquareBase : public QuantityBase{
public:
  QuantitySquareBase(double sqrtq) : QuantityBase(sqrtq*sqrtq), _sqrtq(sqrtq){}
  virtual double description_value() const {return _sqrtq;}
protected:
  double _sqrtq;
};  

// generic_quantity >= minimum
template<typename QuantityType>
class SW_QuantityMin : public SelectorWorker{
public:
  /// detfault ctor (initialises the pt cut)
  SW_QuantityMin(double qmin) : _qmin(qmin) {}

  /// returns true is the given object passes the selection pt cut
  virtual bool pass(const PseudoJet & jet) const {return _qmin(jet) >= _qmin.comparison_value();}

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmin.description() << " >= " << _qmin.description_value();
    return ostr.str();
  }

  virtual bool is_geometric() const { return _qmin.is_geometric();}

protected:
  QuantityType _qmin;     ///< the cut
};


// generic_quantity <= maximum
template<typename QuantityType>
class SW_QuantityMax : public SelectorWorker {
public:
  /// detfault ctor (initialises the pt cut)
  SW_QuantityMax(double qmax) : _qmax(qmax) {}

  /// returns true is the given object passes the selection pt cut
  virtual bool pass(const PseudoJet & jet) const {return _qmax(jet) <= _qmax.comparison_value();}

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmax.description() << " <= " << _qmax.description_value();
    return ostr.str();
  }

  virtual bool is_geometric() const { return _qmax.is_geometric();}

protected:
  QuantityType _qmax;   ///< the cut
};


// generic quantity in [minimum:maximum]
template<typename QuantityType>
class SW_QuantityRange : public SelectorWorker {
public:
  /// detfault ctor (initialises the pt cut)
  SW_QuantityRange(double qmin, double qmax) : _qmin(qmin), _qmax(qmax) {}

  /// returns true is the given object passes the selection pt cut
  virtual bool pass(const PseudoJet & jet) const {
    double q = _qmin(jet); // we could identically use _qmax
    return (q >= _qmin.comparison_value()) && (q <= _qmax.comparison_value());
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << _qmin.description_value() << " <= " << _qmin.description() << " <= " << _qmax.description_value();
    return ostr.str();
  }

  virtual bool is_geometric() const { return _qmin.is_geometric();}

protected:
  QuantityType _qmin;   // the lower cut 
  QuantityType _qmax;   // the upper cut
};


//----------------------------------------------------------------------
/// helper class for selecting on pt
class QuantityPt2 : public QuantitySquareBase{
public:
  QuantityPt2(double pt) : QuantitySquareBase(pt){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.perp2();}
  virtual string description() const {return "pt";}
};  

// returns a selector for a minimum pt
Selector SelectorPtMin(double ptmin) {
  return Selector(new SW_QuantityMin<QuantityPt2>(ptmin));
}

// returns a selector for a maximum pt
Selector SelectorPtMax(double ptmax) {
  return Selector(new SW_QuantityMax<QuantityPt2>(ptmax));
}

// returns a selector for a pt range
Selector SelectorPtRange(double ptmin, double ptmax) {
  return Selector(new SW_QuantityRange<QuantityPt2>(ptmin, ptmax));
}


//----------------------------------------------------------------------
/// helper class for selecting on transverse energy
class QuantityEt2 : public QuantitySquareBase{
public:
  QuantityEt2(double Et) : QuantitySquareBase(Et){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.Et2();}
  virtual string description() const {return "Et";}
};  

// returns a selector for a minimum Et
Selector SelectorEtMin(double Etmin) {
  return Selector(new SW_QuantityMin<QuantityEt2>(Etmin));
}

// returns a selector for a maximum Et
Selector SelectorEtMax(double Etmax) {
  return Selector(new SW_QuantityMax<QuantityEt2>(Etmax));
}

// returns a selector for a Et range
Selector SelectorEtRange(double Etmin, double Etmax) {
  return Selector(new SW_QuantityRange<QuantityEt2>(Etmin, Etmax));
}


//----------------------------------------------------------------------
/// helper class for selecting on energy
class QuantityE : public QuantityBase{
public:
  QuantityE(double E) : QuantityBase(E){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.E();}
  virtual string description() const {return "E";}
};  

// returns a selector for a minimum E
Selector SelectorEMin(double Emin) {
  return Selector(new SW_QuantityMin<QuantityE>(Emin));
}

// returns a selector for a maximum E
Selector SelectorEMax(double Emax) {
  return Selector(new SW_QuantityMax<QuantityE>(Emax));
}

// returns a selector for a E range
Selector SelectorERange(double Emin, double Emax) {
  return Selector(new SW_QuantityRange<QuantityE>(Emin, Emax));
}


//----------------------------------------------------------------------
/// helper class for selecting on mass
class QuantityM2 : public QuantitySquareBase{
public:
  QuantityM2(double m) : QuantitySquareBase(m){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.m2();}
  virtual string description() const {return "mass";}
};  

// returns a selector for a minimum mass
Selector SelectorMassMin(double mmin) {
  return Selector(new SW_QuantityMin<QuantityM2>(mmin));
}

// returns a selector for a maximum mass
Selector SelectorMassMax(double mmax) {
  return Selector(new SW_QuantityMax<QuantityM2>(mmax));
}

// returns a selector for a mass range
Selector SelectorMassRange(double mmin, double mmax) {
  return Selector(new SW_QuantityRange<QuantityM2>(mmin, mmax));
}



//----------------------------------------------------------------------
/// helper for selecting on rapidities: quantity
class QuantityRap : public QuantityBase{
public:
  QuantityRap(double rap) : QuantityBase(rap){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.rap();}
  virtual string description() const {return "rap";}
  virtual bool is_geometric() const { return true;}
};  


/// helper for selecting on rapidities: min
class SW_RapMin : public SW_QuantityMin<QuantityRap>{
public:
  SW_RapMin(double rapmin) : SW_QuantityMin<QuantityRap>(rapmin){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = std::numeric_limits<double>::max();     
    rapmin = _qmin.comparison_value();
  }
};

/// helper for selecting on rapidities: max
class SW_RapMax : public SW_QuantityMax<QuantityRap>{
public:
  SW_RapMax(double rapmax) : SW_QuantityMax<QuantityRap>(rapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = _qmax.comparison_value(); 
    rapmin = -std::numeric_limits<double>::max();
  }
};

/// helper for selecting on rapidities: range
class SW_RapRange : public SW_QuantityRange<QuantityRap>{
public:
  SW_RapRange(double rapmin, double rapmax) : SW_QuantityRange<QuantityRap>(rapmin, rapmax){
    assert(rapmin<=rapmax);
  }
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax = _qmax.comparison_value();      
    rapmin = _qmin.comparison_value(); 
  }
  virtual bool has_known_area() const { return true;} ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * (_qmax.comparison_value()-_qmin.comparison_value());
  }
};

// returns a selector for a minimum rapidity
Selector SelectorRapMin(double rapmin) {
  return Selector(new SW_RapMin(rapmin));
}

// returns a selector for a maximum rapidity
Selector SelectorRapMax(double rapmax) {
  return Selector(new SW_RapMax(rapmax));
}

// returns a selector for a rapidity range
Selector SelectorRapRange(double rapmin, double rapmax) {
  return Selector(new SW_RapRange(rapmin, rapmax));
}


//----------------------------------------------------------------------
/// helper for selecting on |rapidities|
class QuantityAbsRap : public QuantityBase{
public:
  QuantityAbsRap(double absrap) : QuantityBase(absrap){}
  virtual double operator()(const PseudoJet & jet ) const { return abs(jet.rap());}
  virtual string description() const {return "|rap|";}
  virtual bool is_geometric() const { return true;}
};  


/// helper for selecting on |rapidities|: max
class SW_AbsRapMax : public SW_QuantityMax<QuantityAbsRap>{
public:
  SW_AbsRapMax(double absrapmax) : SW_QuantityMax<QuantityAbsRap>(absrapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax =  _qmax.comparison_value(); 
    rapmin = -_qmax.comparison_value();
  }
  virtual bool has_known_area() const { return true;}   ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * _qmax.comparison_value();
  }
};

/// helper for selecting on |rapidities|: max
class SW_AbsRapRange : public SW_QuantityRange<QuantityAbsRap>{
public:
  SW_AbsRapRange(double absrapmin, double absrapmax) : SW_QuantityRange<QuantityAbsRap>(absrapmin, absrapmax){}
  virtual void get_rapidity_extent(double &rapmin, double & rapmax) const{
    rapmax =  _qmax.comparison_value(); 
    rapmin = -_qmax.comparison_value();
  }
  virtual bool has_known_area() const { return true;} ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * (_qmax.comparison_value()-max(_qmin.comparison_value(),0.0)); // this should handle properly absrapmin<0
  }
};

// returns a selector for a minimum |rapidity|
Selector SelectorAbsRapMin(double absrapmin) {
  return Selector(new SW_QuantityMin<QuantityAbsRap>(absrapmin));
}

// returns a selector for a maximum |rapidity|
Selector SelectorAbsRapMax(double absrapmax) {
  return Selector(new SW_AbsRapMax(absrapmax));
}

// returns a selector for a |rapidity| range
Selector SelectorAbsRapRange(double rapmin, double rapmax) {
  return Selector(new SW_AbsRapRange(rapmin, rapmax));
}


//----------------------------------------------------------------------
/// helper for selecting on pseudo-rapidities
class QuantityEta : public QuantityBase{
public:
  QuantityEta(double eta) : QuantityBase(eta){}
  virtual double operator()(const PseudoJet & jet ) const { return jet.eta();}
  virtual string description() const {return "eta";}
  // virtual bool is_geometric() const { return true;} // not strictly only y and phi-dependent
};  

// returns a selector for a pseudo-minimum rapidity
Selector SelectorEtaMin(double etamin) {
  return Selector(new SW_QuantityMin<QuantityEta>(etamin));
}

// returns a selector for a pseudo-maximum rapidity
Selector SelectorEtaMax(double etamax) {
  return Selector(new SW_QuantityMax<QuantityEta>(etamax));
}

// returns a selector for a pseudo-rapidity range
Selector SelectorEtaRange(double etamin, double etamax) {
  return Selector(new SW_QuantityRange<QuantityEta>(etamin, etamax));
}


//----------------------------------------------------------------------
/// helper for selecting on |pseudo-rapidities|
class QuantityAbsEta : public QuantityBase{
public:
  QuantityAbsEta(double abseta) : QuantityBase(abseta){}
  virtual double operator()(const PseudoJet & jet ) const { return abs(jet.eta());}
  virtual string description() const {return "|eta|";}
  virtual bool is_geometric() const { return true;}
};  

// returns a selector for a minimum |pseudo-rapidity|
Selector SelectorAbsEtaMin(double absetamin) {
  return Selector(new SW_QuantityMin<QuantityAbsEta>(absetamin));
}

// returns a selector for a maximum |pseudo-rapidity|
Selector SelectorAbsEtaMax(double absetamax) {
  return Selector(new SW_QuantityMax<QuantityAbsEta>(absetamax));
}

// returns a selector for a |pseudo-rapidity| range
Selector SelectorAbsEtaRange(double absetamin, double absetamax) {
  return Selector(new SW_QuantityRange<QuantityAbsEta>(absetamin, absetamax));
}


//----------------------------------------------------------------------
/// helper for selecting on azimuthal angle
///
/// Note that the bounds have to be specified as min<max
/// phimin has to be > -2pi
/// phimax has to be <  4pi
class SW_PhiRange : public SelectorWorker {
public:
  /// detfault ctor (initialises the pt cut)
  SW_PhiRange(double phimin, double phimax) : _phimin(phimin), _phimax(phimax){
    assert(_phimin<_phimax);
    assert(_phimin>-twopi);
    assert(_phimax<2*twopi);

    _phispan = _phimax - _phimin;
  }

  /// returns true is the given object passes the selection pt cut
  virtual bool pass(const PseudoJet & jet) const {
    double dphi=jet.phi()-_phimin;
    if (dphi >= twopi) dphi -= twopi;
    if (dphi < 0)      dphi += twopi;
    return (dphi <= _phispan);
  }

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << _phimin << " <= phi <= " << _phimax;
    return ostr.str();
  }

  virtual bool is_geometric() const { return true;}

protected:
  double _phimin;   // the lower cut 
  double _phimax;   // the upper cut
  double _phispan;  // the span of the range
};


// returns a selector for a phi range
Selector SelectorPhiRange(double phimin, double phimax) {
  return Selector(new SW_PhiRange(phimin, phimax));
}

//----------------------------------------------------------------------
/// helper for selecting on both rapidity and azimuthal angle
class SW_RapPhiRange : public SW_And{
public:
  SW_RapPhiRange(double rapmin, double rapmax, double phimin, double phimax)
    : SW_And(SelectorRapRange(rapmin, rapmax), SelectorPhiRange(phimin, phimax)){
    _known_area = ((phimax-phimin > twopi) ? twopi : phimax-phimin) * (rapmax-rapmin);
  }

  /// if it has a computable area, return it
  virtual double known_area() const{
    return _known_area;
  }

protected:
  double _known_area;
};

Selector SelectorRapPhiRange(double rapmin, double rapmax, double phimin, double phimax) {
  return Selector(new SW_RapPhiRange(rapmin, rapmax, phimin, phimax));
}


//----------------------------------------------------------------------
/// helper for selecting the n hardest jets
class SW_NHardest : public SelectorWorker {
public:
  /// ctor with specification of the number of objects to keep
  SW_NHardest(unsigned int n) : _n(n) {};

  /// pass makes no sense here normally the parent selector will throw
  /// an error but for internal use in the SW, we'll throw one from
  /// here by security
  virtual bool pass(const PseudoJet &) const {
    if (!applies_jet_by_jet())
      throw Error("Cannot apply this selector worker to an individual jet");
    return false;
  }

  /// For each jet that does not pass the cuts, this routine sets the 
  /// pointer to 0. 
  virtual void terminator(vector<const PseudoJet *> & jets) const {
    // nothing to do if the size is too small
    if (jets.size() < _n) return;

    // do we want to first chech if things are already ordered before
    // going through the ordering process? For now, no. Maybe carry
    // out timing tests at some point to establish the optimal
    // strategy.

    vector<double> minus_pt2(jets.size());
    vector<unsigned int> indices(jets.size());

    for (unsigned int i=0; i<jets.size(); i++){
      indices[i] = i;

      // we need to make sure that the object has not already been
      // nullified.  Note that if we have less than _n jets, this
      // whole n-hardest selection will not have any effect.
      minus_pt2[i] = jets[i] ? -jets[i]->perp2() : 0.0;
    }
    
    IndexedSortHelper sort_helper(& minus_pt2);
    
    partial_sort(indices.begin(), indices.begin()+_n, indices.end(), sort_helper);
    
    for (unsigned int i=_n; i<jets.size(); i++)
      jets[indices[i]] = NULL;
  }
  
  /// returns true if this can be applied jet by jet
  virtual bool applies_jet_by_jet() const {return false;}
  
  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << _n << " hardest";
    return ostr.str();
  }
  
protected:
  unsigned int _n;
};


// returns a selector for the n hardest jets
Selector SelectorNHardest(unsigned int n) {
  return Selector(new SW_NHardest(n));
}



//----------------------------------------------------------------------
// selector and workers for geometric ranges
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// a generic class for objects that contain a position
class SW_WithReference : public SelectorWorker{
public:
  /// ctor
  SW_WithReference() : _is_initialised(false){};

  /// returns true if the worker takes a reference jet
  virtual bool takes_reference() const { return true;}

  /// sets the reference jet
  virtual void set_reference(const PseudoJet &centre){
    _is_initialised = true;
    _reference = centre;
  }

protected:
  PseudoJet _reference;
  bool _is_initialised;
};

//----------------------------------------------------------------------
/// helper for selecting on objects within a distance 'radius' of a reference
class SW_Circle : public SW_WithReference {
public:
  SW_Circle(const double radius) : _radius2(radius*radius) {}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_Circle(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorCircle (or any selector that requires a reference), you first have to call set_reference(...)");
    
    return jet.squared_distance(_reference) <= _radius2;
  } 

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "distance from the centre <= " << sqrt(_radius2);
    return ostr.str();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorCircle (or any selector that requires a reference), you first have to call set_reference(...)");
    
    rapmax = _reference.rap()+sqrt(_radius2);
    rapmin = _reference.rap()-sqrt(_radius2);
  }

  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return pi * _radius2;
  }

protected:
  double _radius2;
};


// select on objets within a distance 'radius' of a variable location
Selector SelectorCircle(const double radius) {
  return Selector(new SW_Circle(radius));
}


//----------------------------------------------------------------------
/// helper for selecting on objects with a distance to a reference
/// betwene 'radius_in' and 'radius_out'
class SW_Doughnut : public SW_WithReference {
public:
  SW_Doughnut(const double radius_in, const double radius_out)
    : _radius_in2(radius_in*radius_in), _radius_out2(radius_out*radius_out) {}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_Doughnut(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorDoughnut (or any selector that requires a reference), you first have to call set_reference(...)");

    double distance2 = jet.squared_distance(_reference);

    return (distance2 <= _radius_out2) && (distance2 >= _radius_in2);
  } 

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << sqrt(_radius_in2) << " <= distance from the centre <= " << sqrt(_radius_out2);
    return ostr.str();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorDoughnut (or any selector that requires a reference), you first have to call set_reference(...)");

    rapmax = _reference.rap()+sqrt(_radius_out2);
    rapmin = _reference.rap()-sqrt(_radius_out2);
  }

  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return pi * (_radius_out2-_radius_in2);
  }

protected:
  double _radius_in2, _radius_out2;
};



// select on objets with distance from the centre is between 'radius_in' and 'radius_out' 
Selector SelectorDoughnut(const double radius_in, const double radius_out) {
  return Selector(new SW_Doughnut(radius_in, radius_out));
}


//----------------------------------------------------------------------
/// helper for selecting on objects with rapidity within a distance 'delta' of a reference
class SW_Strip : public SW_WithReference {
public:
  SW_Strip(const double delta) : _delta(delta) {}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_Strip(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorStrip (or any selector that requires a reference), you first have to call set_reference(...)");
    
    return abs(jet.rap()-_reference.rap()) <= _delta;
  } 

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "|rap - rap_reference| <= " << _delta;
    return ostr.str();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorStrip (or any selector that requires a reference), you first have to call set_reference(...)");
    
    rapmax = _reference.rap()+_delta;
    rapmin = _reference.rap()-_delta;
  }

  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return twopi * 2 * _delta;
  }

protected:
  double _delta;
};


// select on objets within a distance 'radius' of a variable location
Selector SelectorStrip(const double half_width) {
  return Selector(new SW_Strip(half_width));
}


//----------------------------------------------------------------------
/// helper for selecting on objects with rapidity within a distance
/// 'delta_rap' of a reference and phi within a distanve delta_phi of
/// a reference
class SW_Rectangle : public SW_WithReference {
public:
  SW_Rectangle(const double delta_rap, const double delta_phi)
    : _delta_rap(delta_rap),  _delta_phi(delta_phi) {}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_Rectangle(*this);}

  /// returns true if a given object passes the selection criterium
  /// this has to be overloaded by derived workers
  virtual bool pass(const PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorRectangle (or any selector that requires a reference), you first have to call set_reference(...)");

    return (abs(jet.rap()-_reference.rap()) <= _delta_rap) && (abs(jet.delta_phi_to(_reference)) <= _delta_phi);
  } 

  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "|rap - rap_reference| <= " << _delta_rap << " && |phi - phi_reference| <= " << _delta_phi ;
    return ostr.str();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorRectangle (or any selector that requires a reference), you first have to call set_reference(...)");

    rapmax = _reference.rap()+_delta_rap;
    rapmin = _reference.rap()-_delta_rap;
  }

  virtual bool is_geometric() const { return true;}    ///< implies a finite area
  virtual bool has_finite_area() const { return true;} ///< regardless of the reference 
  virtual bool has_known_area() const { return true;}  ///< the area is analytically known
  virtual double known_area() const { 
    return 4 * _delta_rap * _delta_phi;
  }

protected:
  double _delta_rap, _delta_phi;
};


// select on objets within a distance 'radius' of a variable location
Selector SelectorRectangle(const double half_rap_width, const double half_phi_width) {
  return Selector(new SW_Rectangle(half_rap_width, half_phi_width));
}


//----------------------------------------------------------------------
/// helper for selecting the jets that carry at least a given fraction
/// of the reference jet
class SW_PtFractionMin : public SW_WithReference {
public:
  /// ctor with specification of the number of objects to keep
  SW_PtFractionMin(double fraction) : _fraction2(fraction*fraction){}

  /// return a copy of the current object
  virtual SelectorWorker* copy(){ return new SW_PtFractionMin(*this);}

  /// return true if the jet carries a large enough fraction of the reference.
  /// Throw an error if the reference is not initialised.
  virtual bool pass(const PseudoJet & jet) const {
    // make sure the centre is initialised
    if (! _is_initialised)
      throw Error("To use a SelectorPtFractionMin (or any selector that requires a reference), you first have to call set_reference(...)");

    // otherwise, just call that method on the jet
    return (jet.perp2() >= _fraction2*_reference.perp2());
  }
  
  /// returns a description of the worker
  virtual string description() const {
    ostringstream ostr;
    ostr << "pt >= " << sqrt(_fraction2) << "* pt_ref";
    return ostr.str();
  }

protected:
  double _fraction2;
};


// select objects that carry at least a fraction "fraction" of the reference jet
// (Note that this selectir takes a reference)
Selector SelectorPtFractionMin(double fraction){
  return Selector(new SW_PtFractionMin(fraction));
}


//----------------------------------------------------------------------
// additional (mostly helper) selectors
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// helper for selecting the 0-momentum jets
class SW_IsZero : public SelectorWorker {
public:
  /// ctor
  SW_IsZero(){}

  /// return true if the jet has zero momentum
  virtual bool pass(const PseudoJet & jet) const {
    return jet==0;
  }
  
  /// rereturns a description of the worker
  virtual string description() const { return "zero";}
};


// select objects with zero momentum
Selector SelectorIsZero(){
  return Selector(new SW_IsZero());
}


//----------------------------------------------------------------------
#ifndef __FJCORE__
/// helper for selecting the pure ghost
class SW_IsPureGhost : public SelectorWorker {
public:
  /// ctor
  SW_IsPureGhost(){}

  /// return true if the jet is a pure-ghost jet
  virtual bool pass(const PseudoJet & jet) const {
    // if the jet has no area support then it's certainly not a ghost
    if (!jet.has_area()) return false;

    // otherwise, just call that method on the jet
    return jet.is_pure_ghost();
  }
  
  /// rereturns a description of the worker
  virtual string description() const { return "pure ghost";}
};


// select objects that are (or are only made of) ghosts
Selector SelectorIsPureGhost(){
  return Selector(new SW_IsPureGhost());
}

//----------------------------------------------------------------------
// Selector and workers for obtaining a Selector from an old
// RangeDefinition
//
// This is mostly intended for backward compatibility and is likely to
// be removed in a future major release of FastJet
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// helper for selecting on both rapidity and azimuthal angle
class SW_RangeDefinition : public SelectorWorker{
public:
  /// ctor from a RangeDefinition
  SW_RangeDefinition(const RangeDefinition &range) : _range(&range){}

  /// transfer the selection creterium to the underlying RangeDefinition
  virtual bool pass(const PseudoJet & jet) const {
    return _range->is_in_range(jet);
  } 

  /// returns a description of the worker
  virtual string description() const {
    return _range->description();
  }

  /// returns the rapidity range for which it may return "true"
  virtual void get_rapidity_extent(double & rapmin, double & rapmax) const{
    _range->get_rap_limits(rapmin, rapmax);
  }

  /// check if it has a finite area
  virtual bool is_geometric() const { return true;}

  /// check if it has an analytically computable area
  virtual bool has_known_area() const { return true;}
  
  /// if it has a computable area, return it
  virtual double known_area() const{
    return _range->area();
  }

protected:
  const RangeDefinition *_range;
};


// ctor from a RangeDefinition
//----------------------------------------------------------------------
//
// This is provided for backward compatibility and will be removed in
// a future major release of FastJet
Selector::Selector(const RangeDefinition &range) {
  _worker.reset(new SW_RangeDefinition(range));
}
#endif  // __FJCORE__


// operators applying directly on a Selector
//----------------------------------------------------------------------

// operator &=
// For 2 Selectors a and b, a &= b is eauivalent to a = a & b;
Selector & Selector::operator &=(const Selector & b){
  _worker.reset(new SW_And(*this, b));
  return *this;
}

// operator &=
// For 2 Selectors a and b, a &= b is eauivalent to a = a & b;
Selector & Selector::operator |=(const Selector & b){
  _worker.reset(new SW_Or(*this, b));
  return *this;
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
