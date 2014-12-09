///////////////////////////////////////////////////////////////////////////////
// File: vicinity.cpp                                                        //
// Description: source file for particle vicinity (Cvicinity class)          //
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

#include "vicinity.h"
#include <math.h>
#include <algorithm>
#include <iostream>

namespace siscone{

using namespace std;

/*************************************************************
 * Cvicinity_elm implementation                              *
 * element in the vicinity of a parent.                      *
 * class used to manage one points in the vicinity           *
 * of a parent point.                                        *
 *************************************************************/

// ordering pointers to Cvicinity_elm
//------------------------------------
bool ve_less(Cvicinity_elm *ve1, Cvicinity_elm *ve2){
  return ve1->angle < ve2->angle;
}


/*************************************************************
 * Cvicinity implementation                                  *
 * list of element in the vicinity of a parent.              *
 * class used to manage the points which are in the vicinity *
 * of a parent point. The construction of the list can be    *
 * made from a list of points or from a quadtree.            *
 *************************************************************/

// default constructor
//---------------------
Cvicinity::Cvicinity(){
  n_part = 0;

  ve_list = NULL;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  quadtree = NULL;
#endif

  parent = NULL;
  VR2 = VR = 0.0;

}

// constructor with initialisation
//---------------------------------
Cvicinity::Cvicinity(vector<Cmomentum> &_particle_list){
  parent = NULL;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  quadtree = NULL;
#endif
  VR2 = VR = 0.0;

  set_particle_list(_particle_list);
}

// default destructor
//--------------------
Cvicinity::~Cvicinity(){
  if (ve_list!=NULL)
    delete[] ve_list;

#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  if (quadtree!=NULL)
    delete quadtree;
#endif
}

/*
 * set the particle_list
 *  - particle_list   list of particles (type Cmomentum)
 *  - n               number of particles in the list
 ************************************************************/ 
void Cvicinity::set_particle_list(vector<Cmomentum> &_particle_list){
  int i,j;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  double eta_max=0.0;
#endif
  
  // if the particle list is not empty, destroy it !
  if (ve_list!=NULL){
    delete[] ve_list;
  }
  vicinity.clear();
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  if (quadtree!=NULL)
    delete quadtree;
#endif

  // allocate memory array for particles
  // Note: - we compute max for |eta|
  //       - we allocate indices to particles
  n_part = 0;
  plist.clear();
  pincluded.clear();
  for (i=0;i<(int) _particle_list.size();i++){
    // if a particle is colinear with the beam (infinite rapidity)
    // we do not take it into account
    if (fabs(_particle_list[i].pz)!=_particle_list[i].E){
      plist.push_back(_particle_list[i]);
      pincluded.push_back(Cvicinity_inclusion()); // zero inclusion status

      // the parent_index is handled in the split_merge because 
      // of our multiple-pass procedure.
      // Hence, it is not required here any longer.
      // plist[n_part].parent_index = i;
      plist[n_part].index = n_part;

      // make sure the reference is randomly created
      plist[n_part].ref.randomize();

#ifdef USE_QUADTREE_FOR_STABILITY_TEST
      if (fabs(plist[n_part].eta)>eta_max) eta_max=fabs(plist[n_part].eta);
#endif

      n_part++;
    }
  }

  // allocate quadtree and vicinity_elm list
  // note: we set phi in [-pi:pi] as it is the natural range for atan2!
  ve_list = new Cvicinity_elm[2*n_part];
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  eta_max+=0.1;
  quadtree = new Cquadtree(0.0, 0.0, eta_max, M_PI);
#endif

  // append particle to the vicinity_elm list
  j = 0;
  for (i=0;i<n_part;i++){
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
    quadtree->add(&plist[i]);
#endif
    ve_list[j].v = ve_list[j+1].v = &plist[i];
    ve_list[j].is_inside = ve_list[j+1].is_inside = &(pincluded[i]);
    j+=2;
  }

}


/*
 * build the vicinity list from a list of points.
 *  - _parent   reference particle
 *  - _VR       vicinity radius
 ************************************************************/
void Cvicinity::build(Cmomentum *_parent, double _VR){
  int i;

  // set parent and radius
  parent = _parent;
  VR  = _VR;
  VR2 = VR*VR;
  R2  = 0.25*VR2;
  R   = 0.5*VR;
  inv_R_EPS_COCIRC  = 1.0 / R / EPSILON_COCIRCULAR;
  inv_R_2EPS_COCIRC = 0.5 / R / EPSILON_COCIRCULAR;

  // clear vicinity
  vicinity.clear();

  // init parent variables
  pcx = parent->eta;
  pcy = parent->phi;

  // really browse the particle list
  for (i=0;i<n_part;i++){
    append_to_vicinity(&plist[i]);
  }

  // sort the vicinity
  sort(vicinity.begin(), vicinity.end(), ve_less);

  vicinity_size = vicinity.size();
}


/// strictly increasing function of the angle 
inline double sort_angle(double s, double c){
  if (s==0) return (c>0) ? 0.0 : 2.0;
  double t=c/s;
  return (s>0) ? 1-t/(1+fabs(t)) : 3-t/(1+fabs(t));
}


/*
 * append a particle to the 'vicinity' list after
 * having computed the angular-ordering quantities
 *  - v   vector to test
 **********************************************************/
void Cvicinity::append_to_vicinity(Cmomentum *v){
  double dx, dy, d2;

  // skip the particle itself)
  if (v==parent)
    return;

  int i=2*(v->index);

  // compute the distance of the i-th particle with the parent
  dx = v->eta - pcx;
  dy = v->phi - pcy;
    
  // pay attention to the periodicity in phi !
  if (dy>M_PI) 
    dy -= twopi;
  else if (dy<-M_PI) 
    dy += twopi;

  d2 = dx*dx+dy*dy;
    
  // really check if the distance is less than VR
  if (d2<VR2){
    double s,c,tmp;
    
    // compute the angles used for future ordering ...
    //  - build temporary variables used for the computation
    //d   = sqrt(d2);
    tmp = sqrt(VR2/d2-1);

    // first angle (+)
    c = 0.5*(dx-dy*tmp);  // cosine of (parent,child) pair w.r.t. horizontal
    s = 0.5*(dy+dx*tmp);  // sine   of (parent,child) pair w.r.t. horizontal
    ve_list[i].angle = sort_angle(s,c);
    ve_list[i].eta = pcx+c;
    ve_list[i].phi = phi_in_range(pcy+s);
    ve_list[i].side = true;
    ve_list[i].cocircular.clear();
    vicinity.push_back(&(ve_list[i]));

    // second angle (-)    
    c = 0.5*(dx+dy*tmp);  // cosine of (parent,child) pair w.r.t. horizontal
    s = 0.5*(dy-dx*tmp);  // sine   of (parent,child) pair w.r.t. horizontal
    ve_list[i+1].angle = sort_angle(s,c);
    ve_list[i+1].eta = pcx+c;
    ve_list[i+1].phi = phi_in_range(pcy+s);
    ve_list[i+1].side = false;
    ve_list[i+1].cocircular.clear();
    vicinity.push_back(&(ve_list[i+1]));

    // now work out the cocircularity range for the two points (range
    // of angle within which the points stay within a distance
    // EPSILON_COCIRCULAR of circule
    // P = parent; C = child; O = Origin (center of circle)
    Ctwovect OP(pcx - ve_list[i+1].eta, phi_in_range(pcy-ve_list[i+1].phi));
    Ctwovect OC(v->eta - ve_list[i+1].eta, 
                phi_in_range(v->phi-ve_list[i+1].phi));

    // two sources of error are (GPS CCN29-19) epsilon/(R sin theta)
    // and sqrt(2*epsilon/(R (1-cos theta))) and the way things work
    // out, it is the _smaller_ of the two that is relevant [NB have
    // changed definition of theta here relative to that used in
    // CCN29] [NB2: write things so as to avoid zero denominators and
    // to minimize the multiplications, divisions and above all sqrts
    // -- that means that c & s are defined including a factor of VR2]
    c = dot_product(OP,OC);
    s = fabs(cross_product(OP,OC));
    double inv_err1 = s * inv_R_EPS_COCIRC;
    double inv_err2_sq = (R2-c) * inv_R_2EPS_COCIRC;
    ve_list[i].cocircular_range = pow2(inv_err1) > inv_err2_sq ? 
                                                     1.0/inv_err1 : 
                                                     sqrt(1.0/inv_err2_sq);
    ve_list[i+1].cocircular_range = ve_list[i].cocircular_range;
  }
}

}
