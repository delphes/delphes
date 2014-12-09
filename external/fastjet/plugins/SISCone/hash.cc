///////////////////////////////////////////////////////////////////////////////
// File: hash.cpp                                                            //
// Description: source file for classes hash_element and hash_cones          //
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
// $Revision:: 225                                                          $//
// $Date:: 2008-05-20 16:59:47 +0200 (Tue, 20 May 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include "hash.h"
#include <iostream>

namespace siscone{

using namespace std;

/**************************************************************
 * implementation of hash_cones                               *
 * list of cones candidates.                                  *
 * We store in this class all the hash_elements and give      *
 * functions to manipulate them.                              *
 **************************************************************/

// constructor with initialisation
//  - _Np  number of particles
//  - _R2  cone radius (squared)
//-----------------------------------
hash_cones::hash_cones(int _Np, double _R2){
  int i;

  n_cones = 0;
#ifdef DEBUG_STABLE_CONES
  n_occupied_cells = 0;
#endif

  // determine hash size
  // for a ymax=5 and R=0.7, we observed an occupancy around 1/8 N^2 ~ N2 R2/4
  //mask = 1 << (int) (2*log(double(_Np))/log(2.0));
  //if (mask<=1) mask=2;
  int nbits = (int) (log(_Np*_R2*_Np/4.0)/log(2.0));
  if (nbits<1) nbits=1;
  mask = 1 << nbits;

  // create hash
  hash_array = new hash_element*[mask];
  mask--;

  // set the array to 0
  //? needed ?
  for (i=0;i<mask+1;i++)
    hash_array[i] = NULL;

  R2 = _R2;
}

// destructor
//------------
hash_cones::~hash_cones(){
  int i;
  hash_element *elm;

  for (i=0;i<mask+1;i++){
    while (hash_array[i]!=NULL){
      elm = hash_array[i];
      hash_array[i] = hash_array[i]->next;
      delete elm;
    }
  }

  delete[] hash_array;
}


/*
 * insert a new candidate into the hash.
 *  - v       4-momentum of the cone to add
 *  - parent  parent particle defining the cone
 *  - child   child particle defining the cone
 *  - p_io    whether the parent has to belong to the cone or not
 *  - c_io    whether the child has to belong to the cone or not
 * return 0 on success, 1 on error
 ***********************************************************************/
int hash_cones::insert(Cmomentum *v, Cmomentum *parent, Cmomentum *child, bool p_io, bool c_io){
  hash_element *elm;
  int index = (v->ref.ref[0]) & mask;

  // check the array cell corresponding to our reference
  elm = hash_array[index];

#ifdef DEBUG_STABLE_CONES
  if (elm==NULL)
    n_occupied_cells++;
#endif

  do{
    // if it is not present, add it
    if (elm==NULL){
      // create element
      elm = new hash_element;

      // set its varibles
      // Note: at this level, eta and phi have already been computed
      //       through Cmomentum::build_etaphi.
      elm->ref = v->ref;
      
      //compute vectors centre
      v->build_etaphi();
      elm->eta = v->eta;
      elm->phi = v->phi;
      // if at least one of the two is_inside tests gives a result != from the expected,
      // the || will be true hence !(...) false as wanted
      elm->is_stable = !((is_inside(v, parent)^p_io)||(is_inside(v, child)^c_io));
      //cout << "-- new status of " <<  v->ref[0] << ":" << elm->is_stable << endl;

      // update hash
      elm->next = hash_array[index];
      hash_array[index] = elm;
      
      n_cones++;
      return 0;
    }

    // if the cone is already there, simply update stability status
    if (v->ref == elm->ref){
      // there is only an update to perform to see if the cone is still stable
      if (elm->is_stable){
	v->build_etaphi();
	elm->is_stable = !((is_inside(v, parent)^p_io)||(is_inside(v, child)^c_io));
        //cout << " parent/child: " 
        //     << parent->ref[0] << ":" << is_inside(v, parent) << ":" << p_io << " "
        //     << child->ref[0] << ":" << is_inside(v, child) << ":" << c_io << endl;
        //cout << "-- rep status of " <<  v->ref[0] << ":" << elm->is_stable << endl;
        //cout << v->eta << " " << v->phi << endl;
        //cout << (child->eta) << " " << child->phi << endl;
      }
      return 0;
    }

    elm = elm->next;
  } while (1);

  return 1;
}

/*
 * insert a new candidate into the hash.
 *  - v       4-momentum of te cone to add
 * Note, in this case, we assume stability. We also assume
 * that eta and phi are computed for v
 * return 0 on success, 1 on error
 ***********************************************************************/
int hash_cones::insert(Cmomentum *v){
  hash_element *elm;
  int index = (v->ref.ref[0]) & mask;
  //cout << "-- stable candidate: " << v->ref[0] << ":" << endl;

  // check the array cell corresponding to our reference
  elm = hash_array[index];
  do{
    // if it is not present, add it
    if (elm==NULL){
      // create element
      elm = new hash_element;

      // set its varibles
      // Note: at this level, eta and phi have already been computed
      //       through Cmomentum::build_etaphi.
      elm->ref = v->ref;
      elm->eta = v->eta;
      elm->phi = v->phi;
      elm->is_stable = true;

      // update hash
      elm->next = hash_array[index];
      hash_array[index] = elm;
      
      n_cones++;
      return 0;
    }

    // if the cone is already there, we have nothing to do
    if (v->ref == elm->ref){
      return 0;
    }

    elm = elm->next;
  } while (1);

  return 1;
}

/*
 * test if a particle is inside a cone of given centre.
 * check if the particle of coordinates 'v' is inside the circle of radius R 
 * centered at 'centre'.
 *  - centre   centre of the circle
 *  - v        particle to test
 * return true if inside, false if outside
 ******************************************************************************/
inline bool hash_cones::is_inside(Cmomentum *centre, Cmomentum *v){
  double dx, dy;

  dx = centre->eta - v->eta;
  dy = fabs(centre->phi - v->phi);
  if (dy>M_PI) 
    dy -= 2.0*M_PI;
      
  return dx*dx+dy*dy<R2;
}

}
