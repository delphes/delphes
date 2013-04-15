#ifndef __FASTJET_BACKGROUND_ESTIMATOR_BASE_HH__
#define __FASTJET_BACKGROUND_ESTIMATOR_BASE_HH__

//STARTHEADER
// $Id: BackgroundEstimatorBase.hh 2689 2011-11-14 14:51:06Z soyez $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include <fastjet/ClusterSequenceAreaBase.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/Error.hh>
#include <iostream>

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh


/// @ingroup tools_background
/// \class BackgroundEstimatorBase
///
/// Abstract base class that provides the basic interface for classes
/// that estimate levels of background radiation in hadrion and
/// heavy-ion collider events.
///
///
class BackgroundEstimatorBase {
public:
  /// @name  constructors and destructors
  //\{
  //----------------------------------------------------------------
  BackgroundEstimatorBase() : _rescaling_class(0) {}
  //\}

  /// a default virtual destructor that does nothing
  virtual ~BackgroundEstimatorBase() {}


  /// @name setting a new event
  //\{
  //----------------------------------------------------------------

  /// tell the background estimator that it has a new event, composed
  /// of the specified particles.
  virtual void set_particles(const std::vector<PseudoJet> & particles) = 0;

  //\}

  /// @name  retrieving fundamental information
  //\{
  //----------------------------------------------------------------

  /// get rho, the background density per unit area
  virtual double rho() const = 0;

  /// get sigma, the background fluctuations per unit area; must be
  /// multipled by sqrt(area) to get fluctuations for a region of a
  /// given area.
  virtual double sigma() const { 
    throw Error("sigma() not supported for this Background Estimator");
  }

  /// get rho, the background density per unit area, locally at the
  /// position of a given jet. Note that this is not const, because a
  /// user may then wish to query other aspects of the background that
  /// could depend on the position of the jet last used for a rho(jet)
  /// determination.
  virtual double rho(const PseudoJet & jet) = 0;

  /// get sigma, the background fluctuations per unit area, locally at
  /// the position of a given jet. As for rho(jet), it is non-const.
  virtual double sigma(const PseudoJet & /*jet*/) { 
    throw Error("sigma(jet) not supported for this Background Estimator");
  }

  /// returns true if this background estimator has support for
  /// determination of sigma
  virtual bool has_sigma() {return false;}
  //\}
  

  /// @name configuring the behaviour
  //\{
  //----------------------------------------------------------------

  /// Set a pointer to a class that calculates the rescaling factor as
  /// a function of the jet (position). Note that the rescaling factor
  /// is used both in the determination of the "global" rho (the pt/A
  /// of each jet is divided by this factor) and when asking for a
  /// local rho (the result is multiplied by this factor).
  ///
  /// The BackgroundRescalingYPolynomial class can be used to get a
  /// rescaling that depends just on rapidity.
  virtual void set_rescaling_class(const FunctionOfPseudoJet<double> * rescaling_class_in) { _rescaling_class = rescaling_class_in; }

  /// return the pointer to the jet density class
  const FunctionOfPseudoJet<double> *  rescaling_class() const{
    return _rescaling_class;
  }

  //\}

  /// @name description
  //\{
  //----------------------------------------------------------------

  /// returns a textual description of the background estimator
  virtual std::string description() const = 0;

  //\}

protected:
  /// @name helpers for derived classes
  ///
  /// Note that these helpers are related to median-based estimation
  /// of the background, so there is no guarantee that they will
  /// remain in this base class in the long term
  //\{
  //----------------------------------------------------------------

  /// given a quantity in a vector (e.g. pt_over_area) and knowledge
  /// about the number of empty jets, calculate the median and
  /// stand_dev_if_gaussian (roughly from the 16th percentile)
  ///
  /// If do_fj2_calculation is set to true then this performs FastJet
  /// 2.X estimation of the standard deviation, which has a spurious
  /// offset in the limit of a small number of jets.
  void _median_and_stddev(const std::vector<double> & quantity_vector, 
                          double n_empty_jets, 
                          double & median, 
                          double & stand_dev_if_gaussian,
                          bool do_fj2_calculation = false
                          ) const;

  /// computes a percentile of a given _sorted_ vector
  ///  \param sorted_quantity_vector   the vector contains the data sample
  ///  \param percentile               the percentile (defined between 0 and 1) to compute
  ///  \param nempty                   an additional number of 0's
  ///                                  (considered at the beginning of 
  ///                                  the quantity vector)
  ///  \param do_fj2_calculation       carry out the calculation as it
  ///                                  was done in fj2 (suffers from "edge effects")
  double _percentile(const std::vector<double> & sorted_quantity_vector, 
                     const double percentile, 
                     const double nempty=0.0,
                     const bool do_fj2_calculation = false) const;

  //\}

  const FunctionOfPseudoJet<double> * _rescaling_class;
  static LimitedWarning _warnings_empty_area;
};



//----------------------------------------------------------------------
/// @ingroup tools_background
/// A background rescaling that is a simple polynomial in y
class BackgroundRescalingYPolynomial : public FunctionOfPseudoJet<double> {
public:
  /// construct a background rescaling polynomial of the form
  /// a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4
  ///
  /// The following values give a reasonable reproduction of the
  /// Pythia8 tune 4C background shape for pp collisions at
  /// sqrt(s)=7TeV:
  ///
  /// - a0 =  1.157
  /// - a1 =  0
  /// - a2 = -0.0266
  /// - a3 =  0
  /// - a4 =  0.000048
  ///
  BackgroundRescalingYPolynomial(double a0=1, 
                                 double a1=0, 
                                 double a2=0, 
                                 double a3=0, 
                                 double a4=0) : _a0(a0), _a1(a1), _a2(a2), _a3(a3), _a4(a4) {}

  /// return the rescaling factor associated with this jet
  virtual double result(const PseudoJet & jet) const;
private:
  double _a0, _a1, _a2, _a3, _a4;
};





FASTJET_END_NAMESPACE

#endif  // __BACKGROUND_ESTIMATOR_BASE_HH__

