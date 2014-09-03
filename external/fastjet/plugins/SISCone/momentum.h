// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: momentum.h                                                          //
// Description: header file for 4-momentum class Cmomentum                   //
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
// $Revision:: 163                                                          $//
// $Date:: 2007-04-26 22:31:02 +0200 (Thu, 26 Apr 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <vector>
#include <math.h>
#include "reference.h"
#include "geom_2d.h"
#include "defines.h"

namespace siscone{

/**
 * \class Cmomentum
 * \brief base class for dynamic coordinates management
 *
 * This class contains the information for particle or group of
 * particles management.
 * It includes all Lorentz properties as well as tools for summing them.
 * Note: 'sums' over phi angles are indeed averages. This allows to
 *       deal with periodicity at each step
 */
class Cmomentum{
 public:
  /// default ctor
  Cmomentum();

  /// ctor with initialisation
  Cmomentum(double _px, double _py, double _pz, double _E);

  /// ctor with detailed initialisation
  Cmomentum(double _eta, double _phi, Creference _ref);

  /// default dtor
  ~Cmomentum();

  /// computes pT
  inline double perp() const {return sqrt(perp2());}

  /// computes pT^2
  inline double perp2() const {return px*px+py*py;}

  /// computes m
  inline double mass() const {return sqrt(mass2());}

  /// computes m^2
  inline double mass2() const {return perpmass2()-perp2();}

  /// transverse mass, mt = sqrt(pt^2+m^2) = sqrt(E^2 - pz^2)
  inline double perpmass() const {return sqrt((E-pz)*(E+pz));}

  /// transverse mass squared, mt^2 = pt^2+m^2 = E^2 - pz^2
  inline double perpmass2() const {return (E-pz)*(E+pz);}

  /// computes transverse energy
  inline double Et() const {return E/sqrt(1.0+pz*pz/perp2());}

  /// computes transverse energy (squared)
  inline double Et2() const {return E*E/(1.0+pz*pz/perp2());}

  /// assignment of vectors
  Cmomentum& operator = (const Cmomentum &v);

  /// addition of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  const Cmomentum operator + (const Cmomentum &v);

  /// incrementation of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  Cmomentum& operator += (const Cmomentum &v);

  /// decrementation of vectors
  /// !!! WARNING !!! no updating of eta and phi !!!
  Cmomentum& operator -= (const Cmomentum &v);

  /// build eta-phi from 4-momentum info
  /// !!!                WARNING                   !!!
  /// !!! computing eta and phi is time-consuming  !!!
  /// !!! use this whenever you need eta or phi    !!!
  /// !!! automatically called for single-particle !!!
  void build_etaphi();

  double px;        ///< x-momentum
  double py;        ///< y-momentum
  double pz;        ///< z-momentum
  double E;         ///< energy

  double eta;       ///< particle pseudo-rapidity 
  double phi;       ///< particle azimuthal angle
  int parent_index; ///< particle number in the parent list
  int index;        ///< internal particle number

  //////////////////////////////////////////////
  // the following part is used for checksums //
  //////////////////////////////////////////////
  Creference ref;   ///< reference number for the vector
};

/// ordering of two vectors
/// this is by default done w.r.t. their references
bool operator < (const Cmomentum &v1, const Cmomentum &v2);

/// ordering of vectors in eta (e.g. used in collinear tests)
bool momentum_eta_less(const Cmomentum &v1, const Cmomentum &v2);

/// ordering of vectors in pt
bool momentum_pt_less(const Cmomentum &v1, const Cmomentum &v2);


//////////////////////////
// some handy utilities //
//////////////////////////

/// get distance between to eta-phi points
/// \param eta  eta coordinate of first point
/// \param phi  phi coordinate of first point
/// \param v    vector defining the second point
inline double get_distance(double eta, double phi, Cmomentum *v){
  double dx, dy;

  dx = eta - v->eta;
  dy = fabs(phi - v->phi);
  if (dy>M_PI) 
    dy -= twopi;

  return dx*dx+dy*dy;
}

}

#endif
