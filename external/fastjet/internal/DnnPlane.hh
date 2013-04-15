//STARTHEADER
// $Id: DnnPlane.hh 2577 2011-09-13 15:11:38Z salam $
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


#ifndef DROP_CGAL // in case we do not have the code for CGAL

#ifndef __FASTJET_DNNPLANE_HH__
#define __FASTJET_DNNPLANE_HH__

#include "fastjet/internal/Triangulation.hh"
#include "fastjet/internal/DynamicNearestNeighbours.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// \if internal_doc
/// @ingroup internal
/// \class DnnPlane
/// class derived from DynamicNearestNeighbours that provides an
/// implementation for the Euclidean plane
/// \endif
class DnnPlane : public DynamicNearestNeighbours {
 public:
  /// empty initaliser
  DnnPlane() {}

  /// Initialiser from a set of points on an Eta-Phi plane, where both
  /// eta and phi can have arbitrary ranges
  DnnPlane(const std::vector<EtaPhi> &, const bool & verbose = false );


  /// Returns the index of  the nearest neighbour of point labelled
  /// by ii (assumes ii is valid)
  int NearestNeighbourIndex(const int & ii) const ;

  /// Returns the distance to the nearest neighbour of point labelled
  /// by index ii (assumes ii is valid)
  double NearestNeighbourDistance(const int & ii) const ;

  /// Returns true iff the given index corresponds to a point that
  /// exists in the DNN structure (meaning that it has been added, and
  /// not removed in the meantime)
  bool Valid(const int & index) const;

  void RemoveAndAddPoints(const std::vector<int> & indices_to_remove,
			  const std::vector<EtaPhi> & points_to_add,
			  std::vector<int> & indices_added,
			  std::vector<int> & indices_of_updated_neighbours);

  /// returns the EtaPhi of point with index i.
  EtaPhi etaphi(const int i) const;
  /// returns the eta point with index i.
  double eta(const int i) const;
  /// returns the phi point with index i.
  double phi(const int i) const;

 private:

  /// Structure containing a vertex_handle and cached information on
  /// the nearest neighbour.
  struct SuperVertex {
    Vertex_handle vertex; // NULL indicates inexistence...
    double NNdistance;
    int NNindex;
    // later on for cylinder put a second vertex?
  };

  std::vector<SuperVertex> _supervertex;
  //set<Vertex_handle> _vertex_set;
  bool _verbose;

  static const bool _crash_on_coincidence = true;
  //static const bool _crash_on_coincidence = false;

  Triangulation _TR; /// CGAL object for dealing with triangulations

  /// calculates and returns the euclidean distance between points p1
  /// and p2
  inline double _euclid_distance(const Point& p1, const Point& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    return distx*distx+disty*disty;
  }

  //---------------------------------------------------------------------- 
  /// Determines the index and distance of the nearest neighbour to 
  /// point j and puts the information into the _supervertex entry for j
  void _SetNearest(const int & j);

  //----------------------------------------------------------------------
  /// Determines and stores the nearest neighbour of j.
  ///
  /// For each voronoi neighbour D of j if the distance between j and D
  /// is less than D's own nearest neighbour, then update the
  /// nearest-neighbour info in D; push D's index onto 
  /// indices_of_updated_neighbours
  ///
  /// Note that j is NOT pushed onto indices_of_updated_neighbours --
  /// if you want it there, put it there yourself.
  void _SetAndUpdateNearest(const int & j, 
			    std::vector<int> & indices_of_updated_neighbours);

  /// given a vertex_handle returned by CGAL on insertion of a new
  /// points, crash if it turns out that it corresponds to a vertex
  /// that we already knew about (usually because two points coincide)
  void _CrashIfVertexPresent(const Vertex_handle & vertex, 
			     const int & its_index);

};


// here follow some inline implementations of the simpler of the
// functions defined above

inline int DnnPlane::NearestNeighbourIndex(const int & ii) const {
  return _supervertex[ii].NNindex;}

inline double DnnPlane::NearestNeighbourDistance(const int & ii) const {
  return _supervertex[ii].NNdistance;}

inline bool DnnPlane::Valid(const int & index) const {
  if (index >= 0 && index < static_cast<int>(_supervertex.size())) {
    return (_supervertex[index].vertex != NULL);} else {return false;} }

inline EtaPhi DnnPlane::etaphi(const int i) const {
  Point * p = & (_supervertex[i].vertex->point());
  return EtaPhi(p->x(),p->y()); }

inline double DnnPlane::eta(const int i) const {
  return _supervertex[i].vertex->point().x(); }

inline double DnnPlane::phi(const int i) const {
  return _supervertex[i].vertex->point().y(); }


FASTJET_END_NAMESPACE

#endif //  __FASTJET_DNNPLANE_HH__

#endif // DROP_CGAL
