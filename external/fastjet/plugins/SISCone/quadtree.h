// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: quadtree.h                                                          //
// Description: header file for quadtree management (Cquadtree class)        //
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

#ifndef __QUADTREE_H__
#define __QUADTREE_H__

#include "momentum.h"
#include <stdio.h>

namespace siscone{

/**
 * \class Cquadtree
 * \brief Implementation of a 2D quadtree.
 *
 * This class implements the traditional two-dimensional quadtree.
 * The elements at each node are of 'Cmomentum' type.
 */
class Cquadtree{
 public:
  /// default ctor
  Cquadtree();

  /// ctor with initialisation (see init for details)
  Cquadtree(double _x, double _y, double _half_size_x, double _half_size_y);

  /// default destructor
  /// at destruction, everything is destroyed except 
  /// physical values at the leaves
  ~Cquadtree();

  /**
   * init the tree.
   * By initializing the tree, we mean setting the cell parameters
   * and preparing the object to act as a seed for a new tree.
   * \param _x            x-position of the center
   * \param _y            y-position of the center
   * \param _half_size_x  x-size of the cell
   * \param _half_size_y  y-size of the cell
   * \return 0 on success, 1 on error. Note that if the cell or its 
   *         parent is already filled, we return an error.
   */
  int init(double _x, double _y, double _half_size_x, double _half_size_y);

  /**
   * adding a particle to the tree.
   * This method adds one vector to the quadtree structure which 
   * is updated consequently.
   * \param v_add   vector to add
   * \return 0 on success 1 on error
   */
  int add(Cmomentum *v_add);

  /**
   * circle intersection.
   * computes the intersection with a circle of given centre and radius.
   * The output takes the form of a quadtree with all squares included 
   * in the circle.
   * \param cx    circle centre x coordinate
   * \param cy    circle centre y coordinate
   * \param cR2   circle radius SQUARED
   * \return the checksum for that intersection
   */
  Creference circle_intersect(double cx, double cy, double cR2);

  /**
   * output a data file for drawing the grid.
   * This can be used to output a data file containing all the
   * grid subdivisions. The file contents is as follows:
   * first and second columns give center of the cell, the third 
   * gives the size.
   * \param flux  opened stream to write to
   * \return 0 on success, 1 on error
   */
  int save(FILE *flux);

  /**
   * output a data file for drawing the tree leaves.
   * This can be used to output a data file containing all the
   * tree leaves. The file contents is as follows:
   * first and second columns give center of the cell, the third 
   * gives the size.
   * \param flux  opened stream to write to
   * \return 0 on success, 1 on error
   */
  int save_leaves(FILE *flux);

  double centre_x;           ///< x-position of the centre of the cell
  double centre_y;           ///< y-position of the centre of the cell
  double half_size_x;        ///< HALF size of the cell
  double half_size_y;        ///< HALF size of the cell

  Cmomentum *v;              ///< physical contents

  Cquadtree* children[2][2]; ///< sub-cells ( 0,1->left-right; 0,1->bottom,top)
  bool has_child;            ///< true if not a leaf
};

}
#endif
