#ifndef __FASTJET_FUNCTION_OF_PSEUDOJET_HH__
#define __FASTJET_FUNCTION_OF_PSEUDOJET_HH__

//FJSTARTHEADER
// $Id: FunctionOfPseudoJet.hh 3433 2014-07-23 08:17:03Z salam $
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

#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>

FASTJET_BEGIN_NAMESPACE

/// \class FunctionOfPseudoJet
/// base class providing interface for a generic function of a PseudoJet
///
/// This class serves as a base class to provide a standard interface
/// for a function that returns an object of a given (templated) type
/// that depends on a PseudoJet argument. The rationale for using a
/// class (rather than a pointer to a function) is that a class can be
/// constructed with (and store) additional arguments.
template<typename TOut>
class FunctionOfPseudoJet{
public:
  /// default ctor
  FunctionOfPseudoJet(){}

  // ctor that creates a constant function
  //----------
  // this declaration was present in versions of FJ from 3.0.0 to 3.0.6,
  // but never implemented. It is being removed from 3.0.7 upwards
  // to avoid misleading users
  //FunctionOfPseudoJet(const TOut &constant_value);

  /// default dtor (virtual to allow safe polymorphism)
  virtual ~FunctionOfPseudoJet(){}

  /// returns a description of the function (an empty string by
  /// default)
  virtual std::string description() const{ return "";}

  /// the action of the function
  /// this _has_ to be overloaded in derived classes
  ///  \param pj   the PseudoJet input to the function
  virtual TOut result(const PseudoJet &pj) const = 0;

  /// apply the function using the "traditional" () operator.
  /// By default, this just calls the apply(...) method above.
  ///  \param pj   the PseudoJet input to the function
  TOut operator()(const PseudoJet &pj) const { return result(pj);}

  /// apply the function on a vector of PseudoJet, returning a vector
  /// of the results.
  /// This just calls apply on every PseudoJet in the vector.
  ///  \param pjs  the vector of PseudoJet inputs to the function
  std::vector<TOut> operator()(const std::vector<PseudoJet> &pjs) const {
    std::vector<TOut> res(pjs.size());
    for (unsigned int i=0; i<pjs.size(); i++)
      res[i] = result(pjs[i]);
    return res;
  }
};

// The following functions will not be for FJ3.0, because passing a
// reference does not work when the argument is a temporary, which can
// lead to hard-to-diagnose run-time errors. A workaround is to to
// have a pointer rather than a reference as argument, since this
// provides a clearer signal to the user that the object must remain
// in scope.
//
//
// // Selectors created from the ordering between a FunctionOfPseudoJet
// // and a constant
// //----------------------------------------------------------------------
// 
// /// 'larger than' operator
// ///
// /// Select jets for which the given function returns a result larger
// /// than the specified constant
// Selector operator >(const FunctionOfPseudoJet<double> & fn, const double & cut);
// 
// /// 'smaller than' operator
// ///
// /// Select jets for which the given function returns a result smaller
// /// than the specified constant
// Selector operator <(const FunctionOfPseudoJet<double> & fn, const double & cut);
// 
// /// 'larger or equal' operator
// ///
// /// Select jets for which the given function returns a result larger or equal 
// /// to the specified constant
// Selector operator >=(const FunctionOfPseudoJet<double> & fn, const double & cut);
// 
// /// 'smaller or equal' operator
// ///
// /// Select jets for which the given function returns a result smaller or equal
// /// to the specified constant
// Selector operator <=(const FunctionOfPseudoJet<double> & fn, const double & cut);


FASTJET_END_NAMESPACE

#endif  // __FASTJET_FUNCTION_OF_PSEUDOJET_HH__
