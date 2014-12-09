// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: vicinity.h                                                          //
// Description: header file for particle vicinity (Cvicinity class)          //
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
// $Revision:: 123                                                          $//
// $Date:: 2007-03-01 02:52:16 +0100 (Thu, 01 Mar 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __VICINITY_H__
#define __VICINITY_H__

#include <vector>
#include <list>
#include "momentum.h"
#include "defines.h"
#include "quadtree.h"

namespace siscone{

  

/**
 * \class Cvicinity_inclusion
 * \brief a class to keep track of inclusion status in cone and in cocircular region
 *        while using minimal resources
 */
class Cvicinity_inclusion {
public:
  /// default ctor
  Cvicinity_inclusion() : cone(false), cocirc(false) {}

  bool cone;    ///< flag for particle inclusion in the cone
  bool cocirc;  ///< flag for particle inclusion in the border
};


/**
 * \class Cvicinity_elm
 * \brief element in the vicinity of a parent.
 *
 * class used to manage one points in the vicinity 
 * of a parent point.
 */
class Cvicinity_elm{
 public:
  /// pointer to the second borderline particle
  Cmomentum *v;

  /// variable to tell if the particle is inside or outside the cone 
  Cvicinity_inclusion *is_inside;   

  // centre variables
  double eta;              ///< eta coordinate of the center
  double phi;              ///< phi coordinate of the center
  double angle;            ///< angle with parent
  bool side;               ///< true if angle on the positive side, false otherwise
  double cocircular_range; ///< amount by which the angle can be varied while
                           ///< maintaining this point within co-circularity margin

  /// list of elements co-circular with this one
  /// NB: empty list uses less mem than vector
  std::list<Cvicinity_elm * > cocircular;                                          
};

/// ordering pointers to Cvicinity_elm
bool ve_less(Cvicinity_elm *ve1, Cvicinity_elm *ve2);


/**
 * \class Cvicinity
 * \brief list of element in the vicinity of a parent.
 *
 * class used to manage the points which are in the vicinity 
 * of a parent point.
 */
class Cvicinity{
 public:
  /// default constructor
  Cvicinity();

  /// constructor with initialisation (see set_particle_list)
  Cvicinity(std::vector<Cmomentum> &_particle_list);

  /// default destructor
  ~Cvicinity();

  /**
   * set the particle_list
   * \param _particle_list   list of particles (type Cmomentum)
   */ 
  void set_particle_list(std::vector<Cmomentum> &_particle_list);

  /**
   * build the vicinity list from the list of points.
   * \param _parent    reference particle
   * \param _VR        vicinity radius
   */
  void build(Cmomentum *_parent, double _VR);

  // cone kinematical information
  Cmomentum *parent;         ///< parent vector
  double VR;                 ///< radius of the vicinity
  double VR2;                ///< squared radius of the vicinity
  double R;                  ///< normal radius
  double R2;                 ///< squared normal radius
  double inv_R_EPS_COCIRC;   ///< R / EPSILON_COCIRCULAR
  double inv_R_2EPS_COCIRC;  ///< R / (2*EPSILON_COCIRCULAR)

  // particle list information
  int n_part;                                 ///< number of particles
  std::vector<Cmomentum> plist;               ///< the list of particles
  std::vector<Cvicinity_inclusion> pincluded; ///< the inclusion state of particles
  Cvicinity_elm *ve_list;                     ///< list of vicinity elements built from particle list (size=2*n)
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  Cquadtree *quadtree;                        ///< quadtree used for final stability tests
#endif

  // vicinity information
  std::vector<Cvicinity_elm*> vicinity;       ///< list of points in parent's vicinity
  unsigned int vicinity_size;                 ///< number of elements in vicinity

 protected:
  /**
   * append a particle to the 'vicinity' list after
   * having tested it and computed the angular-ordering quantities
   * \param v   vector to test
   */
  void append_to_vicinity(Cmomentum *v);

  // internal variables
  double pcx;    ///< parent centre (eta)
  double pcy;    ///< parent centre (phi)
};

}

#endif
