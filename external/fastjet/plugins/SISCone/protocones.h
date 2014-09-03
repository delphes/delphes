// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: protocones.h                                                        //
// Description: header file for stable cones determination (Cstable_cones)   //
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
// $Revision:: 224                                                          $//
// $Date:: 2008-05-16 19:58:30 +0200 (Fri, 16 May 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __PROTOCONES_H__
#define __PROTOCONES_H__

#include "momentum.h"
#include "vicinity.h"
#include <stdio.h>
#include <vector>
#include <list>
#include "hash.h"

#include "defines.h"

namespace siscone{

/**
 * \class Cborder_store
 * 
 * class for storing a border momentum (in context of co-circularity
 * checks).

 * This class essentially calculates angle of border point w.r.t.
 * circle center (eta & phi), and provides a store of information
 * about whether we are currently including this point in the
 * candidate 
 */
class Cborder_store{
public:
  /// default ctor
  Cborder_store(Cmomentum * momentum, double centre_eta, double centre_phi) : 
    mom(momentum),  is_in(false) {
    angle = atan2(mom->phi - centre_phi, mom->eta - centre_eta);
  }

  Cmomentum * mom;  ///< particle momentum
  double angle;     ///< angle w.r.t. circle centre
  bool   is_in;     ///< inclusion status of the particle
};


/// allows easy sorting of Cborder_store objects (which need to be
/// ordered in angle).
inline bool operator<(const Cborder_store & a, const Cborder_store & b) {
  return a.angle < b.angle;
}


/**
 * \class Cstable_cones
 * \brief Computes the list of stable comes from a particle list.
 *
 * This class does the first fundamental task of te cone algorithm:
 * it is used to compute the list of stable cones given a list
 * of particles.
 */
class Cstable_cones : public Cvicinity{
 public:
  /// default ctor
  Cstable_cones();

  /// ctor with initialisation (sse init for details)
  Cstable_cones(std::vector<Cmomentum> &_particle_list);

  /// default dtor
  ~Cstable_cones();

  /**
   * initialisation
   * \param _particle_list  list of particles
   */
  void init(std::vector<Cmomentum> &_particle_list);

  /**
   * compute stable cones.
   * This function really does the job i.e. computes
   * the list of stable cones (in a seedless way)
   * \param _radius   radius of the cones
   * \return The number of stable cones found is returned
   */
  int get_stable_cones(double _radius);

  /// list of stable cones
  std::vector<Cmomentum> protocones;

  /// list of candidates
  hash_cones *hc;

  /// total number of tested cones
  int nb_tot;
#ifdef DEBUG_STABLE_CONES
  int nb_hash_cones, nb_hash_occupied;
#endif

 protected:
  /// cone radius
  double R;

  /// cone radius SQUARED
  double R2;

 private:
  /// cone with a given particle as parent
  /// this reduction to a single vector assumes we trust the checksums
  Cmomentum cone;

  /// child particle, taken in the 'vicinity' list
  Cmomentum *child;

  /// centre of the tested cone 
  Cvicinity_elm *centre;

  /// index in the particle list;
  unsigned int centre_idx;

  /// first cone used in the vicinity list
  unsigned int first_cone;

  /**
   * initialise the cone.
   * We take the first particle in the angular ordering to compute this one
   * \return 0 on success, 1 on error
   */
  int init_cone();

  /**
   * test cones.
   * We check if the cone(s) build with the present parent and child 
   * are stable
   * \return 0 on success 1 on error
   */
  int test_cone();

  /**
   * update the cone
   * go to the next child for that parent and update 'cone' appropriately
   * \return 0 if update candidate found, 1 otherwise
   */
  int update_cone();

  /*
   * run through the vicinity of the current parent and for each child
   * indicate which members are cocircular...
   */
  void prepare_cocircular_lists();

  /**
   * check if we are in a situation of cocircularity.
   * if it is the case, update and test in the corresponding way
   * \return 'false' if no cocircularity detected, 'true' otherwise
   * Note that if cocircularity is detected, we need to 
   * recall 'update' from 'update' !!!
   */
  bool cocircular_check();

  /**
   * Routine for testing cocircular configurations in p^3 time,
   * rather than 2^p time;
   */
  void test_cone_cocircular(Cmomentum & borderless_cone, 
			    std::list<Cmomentum *> & border_list);

  /**
   * carry out the computations needed for the stability check of the
   * candidate, using the border_vect to indicate which particles
   * should / should not be in the stable cone; if the cone is stable
   * insert it into the hash.
   */
  void test_stability(Cmomentum & candidate, 
                      const std::vector<Cborder_store> & border_vect);

  /**
   * compute the cone contents by going once around the full set of
   * circles and tracking the entry/exit status each time -- this sets
   * up the inclusion information, which can then be directly used to
   * calculate the cone momentum.
   */
  void compute_cone_contents();

  /**
   * compute the cone momentum from particle list.
   * in this version, we use the 'pincluded' information
   * from the Cviinity class
   */
  void recompute_cone_contents();

  /*
   * if we have gone beyond the acceptable threshold of change, compute
   * the cone momentum from particle list.  in this version, we use the
   * 'pincluded' information from the Cvicinity class, but we don't
   * change the member cone, only the locally supplied one
   */
  void recompute_cone_contents_if_needed(Cmomentum & this_cone, double & this_dpt);

  /**
   * compute stability of all enumerated candidates.
   * For all candidate cones which are stable w.r.t. their border particles,
   * pass the last test: stability with quadtree intersection
   */
  int proceed_with_stability();

  /*
   * circle intersection.
   * computes the intersection with a circle of given centre and radius.
   * The output takes the form of a checkxor of the intersection's particles
   *  - cx    circle centre x coordinate
   *  - cy    circle centre y coordinate
   * return the checkxor for the intersection
   ******************************************************************/
  Creference circle_intersect(double cx, double cy);

  /// present candidate cone
  Cmomentum cone_candidate;

  /// in case of co-circular points, vector for them
  std::vector<Cmomentum*> child_list;

  /// list of cocircular enclusures already studied
  /// first element if cone contents, second is cone border
  std::vector< std::pair<Creference,Creference> > multiple_centre_done;

  // information for updating cone contents to avoid rounding errors
  double dpt;          ///< sums of Delta P_t

  /**
   * test if a particle is inside a cone of given centre.
   * check if the particle of coordinates 'v' is inside the circle of radius R 
   * centered at 'centre'.
   * \param centre   centre of the circle
   * \param v        particle to test
   * \return true if inside, false if outside
   */
  inline bool is_inside(Cmomentum *centre, Cmomentum *v);
};

/*
 * compute the absolute value of the difference between 2 angles.
 * We take care of the 2pi periodicity
 * \param angle1   first angle
 * \param angle2   second angle
 * \return the absolute value of the difference between the angles
 *****************************************************************/
inline double abs_dangle(double &angle1, double &angle2);

}
#endif
