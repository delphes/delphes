//FJSTARTHEADER
// $Id: RangeDefinition.hh 3433 2014-07-23 08:17:03Z salam $
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

#ifndef __FASTJET_RANGEDEFINITION_HH__
#define __FASTJET_RANGEDEFINITION_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/Error.hh"
#include "fastjet/LimitedWarning.hh"
#include<sstream>
#include<iostream>
#include<string>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup area_classes
/// \class RangeDefinition
/// class for holding a range definition specification, given by limits
/// on rapidity and azimuth.
///
class RangeDefinition {
public:
  /// default constructor
  RangeDefinition() { _warn_deprecated(); }

  /// constructor for a range definition given by |y|<rapmax
  RangeDefinition(double rapmax) {  _warn_deprecated(); 
                     assert ( rapmax > 0.0 );
                     _rapmax = rapmax;
		     _rapmin = -rapmax;
		     _phimin = 0.0;
		     _phimax = twopi;
		     _total_area = 2.0*rapmax*twopi;
                     _phispan = _phimax-_phimin; }
  
  /// destructor does nothing
  virtual ~RangeDefinition() {}
     
  /// constructor for a range definition given by 
  /// rapmin <= y <= rapmax, phimin <= phi <= phimax
  RangeDefinition(double rapmin, double rapmax, 
                  double phimin = 0.0, double phimax = twopi) {
                     _warn_deprecated(); 
                     assert ( rapmin < rapmax);
                     assert ( phimin < phimax);
                     assert ( phimin > -twopi );
                     assert ( phimax < 2*twopi);
                     _rapmax = rapmax;
		     _rapmin = rapmin;
		     _phimin = phimin;
		     _phimax = phimax;
		     if (_phimax-_phimin > twopi)
		       _total_area = (_rapmax - _rapmin)*twopi;
		     else
		       _total_area = (_rapmax - _rapmin)*(_phimax - _phimin);
                     _phispan = _phimax-_phimin; }

  /// returns true if the range is localizable (i.e. set_position is
  /// meant to do something meaningful). 
  ///
  /// This version of the class is not localizable and so it returns
  /// false.
  ///
  /// For localizable classes override this function with a function
  /// that returns true
  virtual inline bool is_localizable() const { return false; }


  /// place the range on the rap-phi position
  ///
  /// THIS DOES NOT DO ANYTHING FOR THIS CLASS AND IS ONLY THERE
  /// TO FACILITATE DERIVED CLASSES
  ///
  /// DON'T NECESSARILY COUNT ON IT IN THE FUTURE EITHER???
  inline void set_position(const double & rap, const double & phi) {
     if (! is_localizable() ) {
       std::ostringstream err;
       err << description() << 
         "\nThis range is not localizable. set_position() should not be used on it.";         
       throw Error(err.str()); 
     } else {
       _rapjet = rap;
       _phijet = phi;
     }
  }

  /// place the range on the jet position
  inline void set_position(const PseudoJet & jet) {
     set_position(jet.rap(),jet.phi());
  }

  /// return bool according to whether the jet is within the given range
  inline bool is_in_range(const PseudoJet & jet) const {
    double rap = jet.rap();
    double phi = jet.phi();
    return is_in_range(rap,phi);
  }
  
  /// return bool according to whether a (rap,phi) point is in range
  virtual inline bool is_in_range(double rap, double phi) const {
    double dphi=phi-_phimin;
    if (dphi >= twopi) dphi -= twopi;
    if (dphi < 0)      dphi += twopi;
    return  ( rap  >= _rapmin && 
	      rap  <= _rapmax &&
	      dphi <= _phispan );
  }

  /// return the minimal and maximal rapidity of this range; remember to 
  /// replace this if you write a derived class with more complex ranges;
  virtual inline void get_rap_limits(double & rapmin, double & rapmax) const {
    rapmin = _rapmin;
    rapmax = _rapmax;
  }
  
  /// area of the range region
  virtual inline double area() const { return _total_area; }
  
  /// textual description of range
  virtual inline std::string description() const {
    std::ostringstream ostr;
    ostr << "Range: " << _rapmin << " <= y <= "   << _rapmax << ", "
                      << _phimin << " <= phi <= " << _phimax ;
    return ostr.str();
}

protected:
  double _total_area;  // total area of specified range

  /// calculate, and set  _total_area, by calculating which of points on 
  /// a grid (npoints * npoints from -rapmax..rapmax,0..2pi) are contained
  /// in the range; it takes a reasonable time with rapmax = 10,
  /// npoints = 100.
  void _numerical_total_area(double rapmax, int npoints) ;
  double _rapjet,_phijet; // jet position. only used in localizable derived classes
  
private:
  double _rapmin,_rapmax,_phimin,_phimax,_phispan;

  static LimitedWarning _warnings_deprecated;

  /// the use of RangeDefinition is deprecated since FastJet version
  /// 3.0 onwards. Please use Selector instead.  
  /// RangeDefinition is only provided for backward compatibility
  /// reasons and is not guaranteed to work in future releases of
  /// FastJet.
  void _warn_deprecated() const; 
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __FASTJET_RANGEDEFINITION_HH__
