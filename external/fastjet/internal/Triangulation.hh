#ifndef DROP_CGAL // in case we do not have the code for CGAL
#ifndef __FASTJET_TRIANGULATION__
#define __FASTJET_TRIANGULATION__

//STARTHEADER
// $Id: Triangulation.hh 2595 2011-09-23 09:05:04Z salam $
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


// file: examples/Triangulation_2/Voronoi.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// \if internal_doc
/// @ingroup internal
/// \struct K
/// the basic geometrical kernel that lies at the base of all CGAL
/// operations
/// \endif
#ifdef CGAL_SIMPLE_KERNEL
struct K : CGAL::Simple_cartesian<double> {};
#else
struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
#endif // CGAL_SIMPLE_KERNEL

// our extras to help us navigate, find distance, etc.
const int INFINITE_VERTEX=-1;
const int NEW_VERTEX=-2;
const double HUGE_DOUBLE=1e300;

/// \if internal_doc
/// @ingroup internal
/// \struct InitialisedInt
/// A class to provide an "int" with an initial value.
/// \endif
class InitialisedInt {
 private:
  int _val;
 public:
  inline InitialisedInt () {_val=NEW_VERTEX;};
  inline InitialisedInt& operator= (int value) {_val = value; return *this;};
  inline int val() const {return _val;};
};


// We can have triangulations with and without hierarchies -- those with 
// are able to guarantee N ln N time for the construction of a large
// triangulation, whereas those without go as N^{3/2} for points
// sufficiently uniformly distributed in a plane.
//
//#define NOHIERARCHY
#ifdef NOHIERARCHY
typedef CGAL::Triangulation_vertex_base_with_info_2<InitialisedInt,K> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>  Triangulation;
#else
typedef CGAL::Triangulation_vertex_base_with_info_2<InitialisedInt,K> Vbb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb> Vb;
typedef CGAL::Triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds>  Dt;
typedef CGAL::Triangulation_hierarchy_2<Dt> Triangulation;
#endif

typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Point          Point; /// CGAL Point structure
typedef Triangulation::Vertex_circulator Vertex_circulator;
typedef Triangulation::Face_circulator Face_circulator;
typedef Triangulation::Face_handle Face_handle;



FASTJET_END_NAMESPACE

#endif // __FASTJET_TRIANGULATION__
#endif //  DROP_CGAL 
