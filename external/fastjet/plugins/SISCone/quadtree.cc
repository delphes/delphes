///////////////////////////////////////////////////////////////////////////////
// File: quadtree.cpp                                                        //
// Description: source file for quadtree management (Cquadtree class)        //
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
// $Revision:: 320                                                          $//
// $Date:: 2011-11-15 09:54:50 +0100 (Tue, 15 Nov 2011)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "quadtree.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

namespace siscone{

using namespace std;

/*******************************************************************
 * Cquadtree implementation                                        *
 * Implementation of a 2D quadtree.                                *
 * This class implements the traditional two-dimensional quadtree. *
 * The elements at each node are of 'Cmomentum' type.              *
 *******************************************************************/

// default ctor
//--------------
Cquadtree::Cquadtree(){
  v = NULL;

  children[0][0] = children[0][1] = children[1][0] = children[1][1] = NULL;
  has_child = false;
}


// ctor with initialisation (see init for details)
//--------------------------
Cquadtree::Cquadtree(double _x, double _y, double _half_size_x, double _half_size_y){
  v = NULL;

  children[0][0] = children[0][1] = children[1][0] = children[1][1] = NULL;
  has_child = false;

  init(_x, _y, _half_size_x, _half_size_y);
}


// default destructor
// at destruction, everything is destroyed except 
// physical values at the leaves
//------------------------------------------------
Cquadtree::~Cquadtree(){
  if (has_child){
    if (v!=NULL) delete v;
    delete children[0][0];
    delete children[0][1];
    delete children[1][0];
    delete children[1][1];
  }
}


/*
 * init the tree.
 * By initializing the tree, we mean setting the cell parameters
 * and preparing the object to act as a seed for a new tree.
 *  - _x           x-position of the center
 *  - _y           y-position of the center
 *  - half_size_x  half x-size of the cell
 *  - half_size_y  half y-size of the cell
 * return 0 on success, 1 on error. Note that if the cell
 *        is already filled, we return an error.
 ******************************************************************/
int Cquadtree::init(double _x, double _y, double _half_size_x, double _half_size_y){
  if (v!=NULL)
    return 1;

  centre_x = _x;
  centre_y = _y;
  half_size_x = _half_size_x;
  half_size_y = _half_size_y;

  return 0;
}


/*
 * adding a particle to the tree.
 * This method adds one vector to the quadtree structure which 
 * is updated consequently.
 *  - v   vector to add
 * return 0 on success 1 on error
 ******************************************************************/
int Cquadtree::add(Cmomentum *v_add){
  // Description of the method:
  // --------------------------
  // the addition process goes as follows:
  //  1. check if the cell is empty, in which case, add the particle 
  //     here and leave.
  //  2. If there is a unique particle already inside,
  //      (a) create children
  //      (b) forward the existing particle to the appropriate child
  //  3. Add current particle to this cell and forward to the 
  //     adequate child
  // NOTE: we assume in the whole procedure that the particle is 
  //       indeed inside the cell !

  // step 1: the case of empty cells
  if (v==NULL){
    v = v_add;
    return 0;
  }

  // step 2: additional work if 1! particle already present
  //         we use the fact that only 1-particle systems have no child
  if (!has_child){
    double new_half_size_x = 0.5*half_size_x;
    double new_half_size_y = 0.5*half_size_y;
    // create children
    children[0][0] = new Cquadtree(centre_x-new_half_size_x, centre_y-new_half_size_y,
				   new_half_size_x, new_half_size_y);
    children[0][1] = new Cquadtree(centre_x-new_half_size_x, centre_y+new_half_size_y,
				   new_half_size_x, new_half_size_y);
    children[1][0] = new Cquadtree(centre_x+new_half_size_x, centre_y-new_half_size_y,
				   new_half_size_x, new_half_size_y);
    children[1][1] = new Cquadtree(centre_x+new_half_size_x, centre_y+new_half_size_y,
				   new_half_size_x, new_half_size_y);

    has_child = true;

    // forward to child
    //? The following line assumes 'true'==1 and 'false'==0
    // Note: v being a single particle, eta and phi are correct
    children[v->eta>centre_x][v->phi>centre_y]->add(v);

    // copy physical params
    v = new Cmomentum(*v);
  }

  // step 3: add new particle
  // Note: v_add being a single particle, eta and phi are correct
  children[v_add->eta>centre_x][v_add->phi>centre_y]->add(v_add);
  *v+=*v_add;

  return 0;
}


/*
 * circle intersection.
 * computes the intersection with a circle of given centre and radius.
 * The output takes the form of a quadtree with all squares included 
 * in the circle.
 *  - cx    circle centre x coordinate
 *  - cy    circle centre y coordinate
 *  - cR2   circle radius SQUARED
 * return the checksum for the intersection
 ******************************************************************/
Creference Cquadtree::circle_intersect(double cx, double cy, double cR2){
  // Description of the method:
  // --------------------------
  // 1. check if cell is empty => no intersection
  // 2. if cell has 1! particle, check if it is inside the circle.
  //    If yes, add it and return, if not simply return.
  // 3. check if the circle intersects the square. If not, return.
  // 4. check if the square is inside the circle. 
  //    If yes, add it to qt and return.
  // 5. check intersections with children.

  // step 1: if there is no particle inside te square, no reason to go further
  if (v==NULL)
    return Creference();

  double dx, dy;

  // step 2: if there is only one particle inside the square, test if it is in
  //         the circle, in which case return associated reference
  if (!has_child){
    // compute the distance
    // Note: v has only one particle => eta and phi are defined
    dx = cx - v->eta;
    dy = fabs(cy - v->phi);
    if (dy>M_PI) 
      dy -= 2.0*M_PI;

    // test distance
    if (dx*dx+dy*dy<cR2){
      return v->ref;
    }

    return Creference();
  }

  // step 3: check if there is an intersection
  //double ryp, rym;
  double dx_c, dy_c;

  // store distance with the centre of the square
  dx_c = fabs(cx-centre_x);
  dy_c = fabs(cy-centre_y);
  if (dy_c>M_PI) dy_c = 2.0*M_PI-dy_c;

  // compute (minimal) the distance (pay attention to the periodicity in phi).
  dx = dx_c-half_size_x;
  if (dx<0) dx=0;
  dy = dy_c-half_size_y;
  if (dy<0) dy=0;

  // check the distance 
  if (dx*dx+dy*dy>=cR2){
    // no intersection
    return Creference();
  }

  // step 4: check if included

  // compute the (maximal) distance
  dx = dx_c+half_size_x;
  dy = dy_c+half_size_y;
  if (dy>M_PI) dy = M_PI;

  // compute the distance
  if (dx*dx+dy*dy<cR2){
    return v->ref;
  }

  // step 5: the square is not fully in. Recurse to children
  return children[0][0]->circle_intersect(cx, cy, cR2)
    + children[0][1]->circle_intersect(cx, cy, cR2)
    + children[1][0]->circle_intersect(cx, cy, cR2)
    + children[1][1]->circle_intersect(cx, cy, cR2);
}


/*
 * output a data file for drawing the grid.
 * This can be used to output a data file containing all the
 * grid subdivisions. The file contents is as follows:
 * first and second columns give center of the cell, the third 
 * gives the size.
 *  - flux  opened stream to write to
 * return 0 on success, 1 on error
 ******************************************************************/
int Cquadtree::save(FILE *flux){

  if (flux==NULL)
    return 1;

  if (has_child){
    fprintf(flux, "%e\t%e\t%e\t%e\n", centre_x, centre_y, half_size_x, half_size_y);
    children[0][0]->save(flux);
    children[0][1]->save(flux);
    children[1][0]->save(flux);
    children[1][1]->save(flux);
  }

  return 0;
}


/*
 * output a data file for drawing the tree leaves.
 * This can be used to output a data file containing all the
 * tree leaves. The file contents is as follows:
 * first and second columns give center of the cell, the third 
 * gives the size.
 *  - flux  opened stream to write to
 * return 0 on success, 1 on error
 ******************************************************************/
int Cquadtree::save_leaves(FILE *flux){

  if (flux==NULL)
    return 1;

  if (has_child){
    if (children[0][0]!=NULL) children[0][0]->save_leaves(flux);
    if (children[0][1]!=NULL) children[0][1]->save_leaves(flux);
    if (children[1][0]!=NULL) children[1][0]->save_leaves(flux);
    if (children[1][1]!=NULL) children[1][1]->save_leaves(flux);
  } else {
    fprintf(flux, "%e\t%e\t%e\t%e\n", centre_x, centre_y, half_size_x, half_size_y);
  }

  return 0;
}

}
