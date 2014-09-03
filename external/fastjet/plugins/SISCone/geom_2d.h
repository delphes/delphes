// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: geom_2d.h                                                           //
// Description: header file for two-dimensional geometry tools               //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 268                                                          $//
// $Date:: 2009-03-12 21:24:16 +0100 (Thu, 12 Mar 2009)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __GEOM_2D_H__
#define __GEOM_2D_H__

#include <iostream>
#include <math.h>
#include "defines.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197
#endif

namespace siscone{

/// return a result that corresponds to phi, but in the
/// range (-pi..pi]; the result is only correct if -3pi < phi <= 3pi
inline double phi_in_range(double phi) {
  if      (phi <= -M_PI) phi += twopi;
  else if (phi >   M_PI) phi -= twopi;
  return phi;
}

/// return the difference between the two phi values, 
/// placed in the correct range (-pi..pi], , assuming that phi1,phi2
/// are already in the correct range.
inline double dphi(double phi1, double phi2) {
  return phi_in_range(phi1-phi2);
}


/// return the absolute difference between the two phi values, 
/// placed in the correct range, assuming that phi1,phi2 are already
/// in the correct range.
inline double abs_dphi(double phi1, double phi2) {
  double delta = fabs(phi1-phi2);
  return delta > M_PI ? twopi-delta : delta;
}

/// return the square of the argument
inline double pow2(double x) {return x*x;}


/** 
 * \class Ctwovect
 * \brief class for holding a two-vector
 */
class Ctwovect {
public:
  /// default ctor
  Ctwovect() : x(0.0), y(0.0) {}

  /// ctor with initialisation
  /// \param _x   first coordinate
  /// \param _y   second coordinate
  Ctwovect(double _x, double _y) : x(_x), y(_y) {}

  /// vector coordinates
  double x, y;

  /// norm (modulud square) of the vector
  inline double mod2() const {return pow2(x)+pow2(y);}

  /// modulus of the vector
  inline double modulus() const {return sqrt(mod2());}
};


/// dot product of two 2-vectors
/// \param a   first 2-vect
/// \param b   second 2-vect
/// \return a.b is returned
inline double dot_product(const Ctwovect & a, const Ctwovect & b) {
  return a.x*b.x + a.y*b.y;
}


/// cross product of two 2-vectors
/// \param a   first 2-vect
/// \param b   second 2-vect
/// \return a x b is returned
inline double cross_product(const Ctwovect & a, const Ctwovect & b) {
  return a.x*b.y - a.y*b.x;
}


/** 
 * \class Ceta_phi_range
 * \brief class for holding a covering range in eta-phi
 *
 * This class deals with ranges in the eta-phi plane. It
 * implements methods to test if two ranges overlap and
 * to take the union of two overlapping intervals.
 */
class Ceta_phi_range{
public:
  /// default ctor
  Ceta_phi_range();

  /// ctor with initialisation
  /// we initialise with a centre (in eta,phi) and a radius
  /// \param c_eta   eta coordinate of the centre
  /// \param c_phi   phi coordinate of the centre
  /// \param R       radius
  Ceta_phi_range(double c_eta, double c_phi, double R);

  /// assignment of range
  /// \param r   range to assign to current one
  Ceta_phi_range& operator = (const Ceta_phi_range &r);

  /// add a particle to the range
  /// \param eta  eta coordinate of the particle
  /// \param phi  phi coordinate of the particle
  /// \return 0 on success, 1 on error
  int add_particle(const double eta, const double phi);

  /// eta range as a binary coding of covered cells
  unsigned int eta_range;     

  /// phi range as a binary coding of covered cells
  unsigned int phi_range;     

  // extremal value for eta
  static double eta_min;  ///< minimal value for eta
  static double eta_max;  ///< maximal value for eta

private:
  /// return the cell index corrsponding to an eta value
  inline unsigned int get_eta_cell(double eta){
    return (unsigned int) (1 << ((int) (32*((eta-eta_min)/(eta_max-eta_min)))));
  }

  /// return the cell index corrsponding to a phi value
  inline unsigned int get_phi_cell(double phi){
    return (unsigned int) (1 << ((int) (32*phi/twopi+16)%32));
  }
};

/// test overlap
/// \param  r1  first range
/// \param  r2  second range
/// \return true if overlap, false otherwise.
bool is_range_overlap(const Ceta_phi_range &r1, const Ceta_phi_range &r2);

/// compute union
/// Note: we assume that the two intervals overlap
/// \param  r1  first range
/// \param  r2  second range
/// \return union of the two ranges
const Ceta_phi_range range_union(const Ceta_phi_range &r1, const Ceta_phi_range &r2);

}

#endif
