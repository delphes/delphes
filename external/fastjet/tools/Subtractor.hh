//STARTHEADER
// $Id: Subtractor.hh 2577 2011-09-13 15:11:38Z salam $
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

#ifndef __FASTJET_TOOLS_SUBTRACTOR_HH__
#define __FASTJET_TOOLS_SUBTRACTOR_HH__

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
    _bge(bge), _rho(-1.0) {}

  /// define a subtractor that uses a fixed value of rho, the background
  /// pt density per unit area (which must be positive)
  Subtractor(double rho);

  /// default constructor
  Subtractor() : _bge(0), _rho(_invalid_rho) {}

  /// default dtor
  virtual ~Subtractor(){};

  /// returns a jet that's subtracted
  ///
  /// \param jet    the jet that is to be subtracted
  /// \return       the subtracted jet
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// class description
  virtual std::string description() const;

protected:

  /// the tool used to estimate the background
  /// if has to be mutable in case its underlying selector takes a reference jet
  mutable BackgroundEstimatorBase * _bge;
  /// the fixed value of rho to use if the user has selected that option
  double _rho;

  /// a value of rho that is used as a default to label that the stored
  /// rho is not valid for subtraction. 
  //
  // NB: there are two reasons for not having the value written here:
  // 1) that it caused problems on karnak with g++ 4.0.1 and 2) that
  // we anyway like -infinity as a default, and since that's a function,
  // that's not allowed in an include file.
  static const double _invalid_rho;
};

FASTJET_END_NAMESPACE

#endif  // __FASTJET_TOOLS_SUBTRACTOR_HH__

