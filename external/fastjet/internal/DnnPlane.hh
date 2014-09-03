//FJSTARTHEADER
// $Id: DnnPlane.hh 3442 2014-07-24 07:20:49Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER


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
///
/// This class that uses CGAL Delaunay triangulation for most of the
/// work (it allows for easy and efficient removal and addition of
/// points and circulation over a point's neighbours). The treatment
/// of coincident points is not supported by CGAL and is implemented
/// according to the method specified in
/// issue-tracker/2012-02-CGAL-coincident/METHOD
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
  int NearestNeighbourIndex(const int ii) const ;

  /// Returns the distance to the nearest neighbour of point labelled
  /// by index ii (assumes ii is valid)
  double NearestNeighbourDistance(const int ii) const ;

  /// Returns true iff the given index corresponds to a point that
  /// exists in the DNN structure (meaning that it has been added, and
  /// not removed in the meantime)
  bool Valid(const int index) const;

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
    int coincidence;  // ==vertex->info.val() if no coincidence
                      // points to the coinciding SV in case of coincidence
    // later on for cylinder put a second vertex?
  };

  std::vector<SuperVertex> _supervertex;
  //set<Vertex_handle> _vertex_set;
  bool _verbose;

  //static const bool _crash_on_coincidence = true;
  static const bool _crash_on_coincidence = false;

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
  void _SetNearest(const int j);

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
  void _SetAndUpdateNearest(const int j, 
			    std::vector<int> & indices_of_updated_neighbours);

  /// given a vertex_handle returned by CGAL on insertion of a new
  /// points, returns the coinciding vertex's value if it turns out
  /// that it corresponds to a vertex that we already knew about
  /// (usually because two points coincide)
  int _CheckIfVertexPresent(const Vertex_handle & vertex, 
			    const int its_index);

  //----------------------------------------------------------------------
  /// if the distance between 'pref' and 'candidate' is smaller (or
  /// equal) than the one between 'pref' and 'near', return true and
  /// set 'mindist' to that distance. Note that it is assumed that
  /// 'mindist' is the euclidian distance between 'pref' and 'near'
  ///
  /// Note that the 'near' point is passed through its vertex rather
  /// than as a point. This allows us to handle cases where we have no min
  /// yet (near is the infinite vertex)
  inline bool _is_closer_to(const Point &pref, 
			    const Point &candidate,
			    const Vertex_handle &near,
			    double & dist,
			    double & mindist){
    dist = _euclid_distance(pref, candidate);
    return _is_closer_to_with_hint(pref, candidate, near, dist, mindist);
  }

  /// same as '_is_closer_to' except that 'dist' already contains the
  /// distance between 'pref' and 'candidate'
  inline bool _is_closer_to_with_hint(const Point &pref, 
				      const Point &candidate,
				      const Vertex_handle &near,
				      const double & dist,
				      double & mindist){
    
    // check if 'dist', the pre-computed distance between 'candidate'
    // and 'pref' is smaller than the distance between 'pref' and its
    // currently registered nearest neighbour 'near' (and update
    // things if it is)
    //
    // Interestingly enough, it has to be pointed out that the use of
    // 'abs' instead of 'std::abs' returns wrong results (apparently
    // ints without any compiler warning)
    //
    // The (near != NULL) test is there for one single reason: when
    // checking that a newly inserted point is not closer than a
    // previous NN, if that distance comparison involves a "nearly
    // degenerate" distance we need to access near->point. But
    // sometimes, in the course of RemoveAndAddPoints, its previous NN
    // has been deleted and its vertex (corresponding to 'near') set
    // to NULL. This is not a problem as all points having a deleted
    // point as NN will have their NN explicitly recomputed at the end
    // of RemoveAndAddPoints so here we should just make sure there is
    // no crash... that's done by checking (near != NULL)
    if ((std::abs(dist-mindist)<DISTANCE_FOR_CGAL_CHECKS) &&
	(near != NULL) &&
	(_euclid_distance(candidate, near->point())<DISTANCE_FOR_CGAL_CHECKS)){
      // we're in a situation where there might be a rounding issue,
      // use CGAL's distance computation to get it right
      //
      // Note that in the test right above,
      // (abs(dist-mindist)<1e-12) guarantees that the current
      // nearest point is not the infinite vertex and thus
      // nearest->point() is not ill-defined
      if (_verbose) std::cout << "using CGAL's distance ordering" << std::endl;
      if (CGAL::compare_distance_to_point(pref, candidate, near->point())!=CGAL::LARGER){
	mindist = dist;
	return true;
      }
    } else if (dist <= mindist) {
      // Note that the use of a <= in the above expression (instead of
      // a strict ordering <) is important in one case: when checking
      // if a new point is the new NN of one of the points in its
      // neighbourhood, in case of distances being ==, we are sure
      // that 'candidate' is in a cell adjacent to 'pref' while it may
      // no longer be the case for 'near'
      mindist = dist;
      return true;
    } 
    
    return false;
  }

  /// if a distance between a point and 2 others is smaller than this
  /// and the distance between the two points is also smaller than this
  /// then use CGAL to compare the distances. 
  static const double DISTANCE_FOR_CGAL_CHECKS;  
  
};


// here follow some inline implementations of the simpler of the
// functions defined above

inline int DnnPlane::NearestNeighbourIndex(const int ii) const {
  return _supervertex[ii].NNindex;}

inline double DnnPlane::NearestNeighbourDistance(const int ii) const {
  return _supervertex[ii].NNdistance;}

inline bool DnnPlane::Valid(const int index) const {
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
