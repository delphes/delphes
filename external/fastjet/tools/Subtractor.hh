//FJSTARTHEADER
// $Id: Subtractor.hh 3670 2014-09-08 14:17:59Z soyez $
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

#ifndef __FASTJET_TOOLS_SUBTRACTOR_HH__
#define __FASTJET_TOOLS_SUBTRACTOR_HH__

#include "fastjet/internal/base.hh"     // namespace macros (include explicitly to help Doxygen)
#include "fastjet/tools/Transformer.hh" // to derive Subtractor from Transformer
#include "fastjet/tools/BackgroundEstimatorBase.hh" // used as a ctor argument

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
/// @ingroup tools_background
/// \class Subtractor
/// Class that helps perform jet background subtraction.
///
/// This class derives from Transformer and makes use of a pointer to
/// a BackgroundEstimatorBase object in order to determine the background
/// in the vicinity of a given jet and then subtract area*background from
/// the jet. It can also be initialised with a specific fixed value for the 
/// background pt density.
///
/// \section input Input conditions
/// 
/// The original jet must have area support (4-vector)
///
/// \section output Output/interface
/// 
/// The underlying structure of the returned, subtracted jet
/// (i.e. constituents, pieces, etc.) is identical to that of the
/// original jet.
///
class Subtractor : public Transformer{
public:
  /// define a subtractor based on a BackgroundEstimator
  Subtractor(BackgroundEstimatorBase * bge) : 
    _bge(bge), _rho(-1.0) { set_defaults(); }

  /// define a subtractor that uses a fixed value of rho, the background
  /// pt density per unit area (which must be positive)
  Subtractor(double rho);

  /// define a subtractor that uses a fixed value of rho and rho_m;
  /// both must be >= 0;
  Subtractor(double rho, double rho_m);

  /// default constructor
  Subtractor() : _bge(0), _rho(_invalid_rho) { set_defaults(); }

  /// default dtor
  virtual ~Subtractor(){};

  /// @name configuring the behaviour
  //\{
  //----------------------------------------------------------------

  /// reset all parameters to default values
  ///
  /// Note: by default, the rho_m term is not included and the safety
  /// test for the mass is not done. This is mostly for backwards
  /// compatibility with FastJet 3.0 and is highly likely to change in
  /// a future release of FastJet
  void set_defaults();

  /// when 'use_rho_m' is true, include in the subtraction the
  /// correction from rho_m, the purely longitudinal,
  /// particle-mass-induced component of the background density per
  /// unit area
  ///
  /// Note: this will be switched off by default (for backwards
  /// compatibility with FastJet 3.0) but is highly likely to change
  /// in a future release of FastJet
  void set_use_rho_m(bool use_rho_m_in = true){
    if (_bge == 0  && _rho_m < 0) {
      throw Error("Subtractor: rho_m support works only for Subtractors constructed with a background estimator or an explicit rho_m value");
    }
    _use_rho_m=use_rho_m_in;
  }
  
  /// returns whether or not the rho_m component is used
  bool use_rho_m() const{ return _use_rho_m;}

  /// when 'safe_mass' is true, ensure that the mass of the subtracted
  /// 4-vector remain positive
  ///
  /// when true, if the subtracted mass is negative, we return a
  /// 4-vector with 0 mass, pt and phi from the subtracted 4-vector
  /// and the rapidity of the original, unsubtracted jet.
  ///
  /// Note: this will be switched off by default (for backwards
  /// compatibility with FastJet 3.0) but is highly likely to change
  /// in a future release of FastJet
  void set_safe_mass(bool safe_mass_in=true){ _safe_mass=safe_mass_in;}

  /// returns whether or not safety tests on the mass are included
  bool safe_mass() const{ return _safe_mass;}

  /// This is mostly intended for cherge-hadron-subtracted type of
  /// events where we wich to use vertex information to improve the
  /// subtraction.
  ///
  /// Given the following parameters:
  ///   \param sel_known_vertex    selects the particles with a
  ///                              known vertex origin
  ///   \param sel_leading_vertex  amongst the particles with a
  ///                              known vertex origin, select those
  ///                              coming from the leading vertex
  /// Momentum identified as coming from the leading vertex will be
  /// kept, momentum identified as coming from a non-leading vertex
  /// will be eliminated and a regular area-median subtraction will be
  /// applied on the 4-vector sum of the particles with unknown vertex
  /// origin.
  ///
  /// When this is set, we shall ensure that the pt of the subtracted
  /// 4-vector is at least the pt of the particles that are known to
  /// come from the leading vertex (if it fails, subtraction returns
  /// the component that is known to come from the leading vertex ---
  /// or, the original unsubtracted jet if it contains no particles
  /// from the leading vertex).  Furthermore, when safe_mass() is on, we
  /// also impose a similar constraint on the mass of the subtracted
  /// 4-vector (if the test fails, the longitudinal part of the
  /// subtracted 4-vector is taken from the component that is known to
  /// come from the leading vertex).
  void set_known_selectors(const Selector &sel_known_vertex,
			   const Selector &sel_leading_vertex){
    _sel_known_vertex   = sel_known_vertex;
    _sel_leading_vertex = sel_leading_vertex;
  }

  //\}

  /// @name description and action
  //\{
  //----------------------------------------------------------------

  /// returns a jet that's subtracted
  ///
  /// \param jet    the jet that is to be subtracted
  /// \return       the subtracted jet
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// class description
  virtual std::string description() const;

  //\}
protected:
  /// compute the 4-vector that should be subtracted from the given
  /// jet
  PseudoJet _amount_to_subtract(const PseudoJet &jet) const;

  /// the tool used to estimate the background
  /// if has to be mutable in case its underlying selector takes a reference jet
  mutable BackgroundEstimatorBase * _bge;
  /// the fixed value of rho and/or rho_m to use if the user has selected that option
  double _rho, _rho_m;

  // configuration parameters/flags
  bool _use_rho_m;   ///< include the rho_m correction
  bool _safe_mass;   ///< ensures that the subtracted mass is +ve

  Selector _sel_known_vertex;   ///< selects the particles with a
				///< known vertex origin
  Selector _sel_leading_vertex; ///< amongst the particles with a
				///< known vertex origin, select those
				///< coming from the leading vertex

  /// a value of rho that is used as a default to label that the stored
  /// rho is not valid for subtraction. 
  //
  // NB: there are two reasons for not having the value written here:
  // 1) that it caused problems on karnak with g++ 4.0.1 and 2) that
  // we anyway like -infinity as a default, and since that's a function,
  // that's not allowed in an include file.
  static const double _invalid_rho;

  mutable LimitedWarning _unused_rho_m_warning;
};

FASTJET_END_NAMESPACE

#endif  // __FASTJET_TOOLS_SUBTRACTOR_HH__

