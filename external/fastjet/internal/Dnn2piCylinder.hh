//STARTHEADER
// $Id: Dnn2piCylinder.hh 2577 2011-09-13 15:11:38Z salam $
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
#ifndef __FASTJET_DNN2PICYLINDER_HH__
#define __FASTJET_DNN2PICYLINDER_HH__

#include "fastjet/internal/DynamicNearestNeighbours.hh"
#include "fastjet/internal/DnnPlane.hh"
#include "fastjet/internal/numconsts.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// \if internal_doc
/// @ingroup internal
/// \class Dnn2piCylinder
/// class derived from DynamicNearestNeighbours that provides an
/// implementation for the surface of cylinder (using one 
/// DnnPlane object spanning 0--2pi).
/// \endif
class Dnn2piCylinder : public DynamicNearestNeighbours {
 public:
  /// empty initaliser
  Dnn2piCylinder() {}

  /// Initialiser from a set of points on an Eta-Phi plane, where
  /// eta can have an arbitrary ranges and phi must be in range
  /// 0 <= phi < 2pi;
  /// 
  /// NB: this class is more efficient than the plain Dnn4piCylinder
  /// class, but can give wrong answers when the nearest neighbour is
  /// further away than 2pi (in this case a point's nearest neighbour
  /// becomes itself, because it is considered to be a distance 2pi
  /// away). For the kt-algorithm (e.g.) this is actually not a
  /// problem (the distance need only be accurate when it is less than
  /// R, assuming R<2pi [not necessarily always the case as of
  /// 2010-11-19, when we've removed the requirement R<pi/2 in the
  /// JetDefinition constructor]), so we can tell the routine to
  /// ignore this problem -- alternatively the routine will crash if
  /// it detects it occurring (only when finding the nearest neighbour
  /// index, not its distance).
  Dnn2piCylinder(const std::vector<EtaPhi> &,
		 const bool & ignore_nearest_is_mirror = false,
		 const bool & verbose = false );

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

  ~Dnn2piCylinder();

 private:

  // our extras to help us navigate, find distance, etc.
  const static int INEXISTENT_VERTEX=-3;

  bool _verbose;

  bool _ignore_nearest_is_mirror;

  /// Picture of how things will work... Copy 0--pi part of the 0--2pi
  /// cylinder into a region 2pi--2pi+ a bit of a Euclidean plane. Below we
  /// show points labelled by + that have a mirror image in this
  /// manner, while points labelled by * do not have a mirror image.
  ///      	 			  
  ///      |           .     |		  
  ///      |	       .     |		  
  ///      |           .     |		  
  ///      |           .     |		  
  ///      |        2  .     |		  
  ///      |        *  .     |		  
  ///      | +         . +   |		  
  ///      | 0         . 1   |
  ///      |	       .     |
  ///      0          2pi   2pi + a bit
  ///   	     
  /// Each "true" point has its true "cylinder" index (the index that
  /// is known externally to this class) as well as euclidean plane
  /// indices (main_index and mirror index in the MirrorVertexInfo
  /// structure), which are private concepts of this class.
  /// 
  /// In above picture our structures would hold the following info
  /// (the picture shows the euclidean-plane numbering)
  ///
  /// _mirror_info[cylinder_index = 0] = (0, 1)
  /// _mirror_info[cylinder_index = 1] = (2, INEXISTENT_VERTEX)
  ///
  /// We also need to be able to go from the euclidean plane indices
  /// back to the "true" cylinder index, and for this purpose we use
  /// the std::vector _cylinder_index_of_plane_vertex[...], which in the above example has
  /// the following contents
  ///
  /// _cylinder_index_of_plane_vertex[0] = 0
  /// _cylinder_index_of_plane_vertex[1] = 0
  /// _cylinder_index_of_plane_vertex[2] = 1
  ///

  /// 
  struct MirrorVertexInfo {
    /// index of the given point (appearing in the range 0--2pi) in the 
    /// 0--2pi euclidean plane structure (position will coincide with
    /// that on the 0--2pi cylinder, but index labelling it will be
    /// different)
    int main_index; 
    /// index of the mirror point (appearing in the range 2pi--3pi) in the
    /// 0--3pi euclidean plane structure
    int mirror_index; 
  };

  // for each "true" vertex we have reference to indices in the euclidean
  // plane structure
  std::vector<MirrorVertexInfo> _mirror_info;
  // for each index in the euclidean 0--2pi plane structure we want to
  // be able to get back to the "true" vertex index on the overall
  // 0--2pi cylinder structure
  std::vector<int> _cylinder_index_of_plane_vertex;

  // NB: we define POINTERS here because the initialisation gave
  //     us problems (things crashed!), perhaps because in practice
  //     we were making a copy without being careful and defining
  //     a proper copy constructor.
  DnnPlane * _DNN;

  /// given a phi value in the 0--pi range return one 
  /// in the 2pi--3pi range; whereas if it is in the pi-2pi range then
  /// remap it to be inthe range (-pi)--0.
  inline EtaPhi _remap_phi(const EtaPhi & point) {
    double phi = point.second;
    if (phi < pi) { phi += twopi ;} else {phi -= twopi;}
    return EtaPhi(point.first, phi);}


  //----------------------------------------------------------------------
  /// Actions here are similar to those in the
  /// Dnn3piCylinder::_RegisterCylinderPoint case, however here we do
  /// NOT create the mirror point -- instead we initialise the structure
  /// as if there were no need for the mirror point.
  ///
  /// ADDITIONALLY push the cylinder_point onto the vector plane_points.
  void _RegisterCylinderPoint (const EtaPhi & cylinder_point,
			       std::vector<EtaPhi> & plane_points);

  /// For each plane point specified in the vector plane_indices,
  /// establish whether there is a need to create a mirror point
  /// according to the following criteria:
  ///
  /// . phi < pi
  /// . mirror does not already exist
  /// . phi < NearestNeighbourDistance 
  ///   (if this is not true then there is no way that its mirror point
  ///   could have a nearer neighbour).
  ///
  /// If conditions all hold, then create the mirror point, insert it
  /// into the _DNN structure, adjusting any nearest neighbours, and
  /// return the list of plane points whose nearest neighbours have
  /// changed (this will include the new neighbours that have just been
  /// added)
  void _CreateNecessaryMirrorPoints(
			  const std::vector<int> & plane_indices,
			  std::vector<int> & updated_plane_points);

};


// here follow some inline implementations of the simpler of the
// functions defined above

//----------------------------------------------------------------------
/// Note: one of the difficulties of the 0--2pi mapping is that
/// a point may have its mirror copy as its own nearest neighbour
/// (if no other point is within a distance of 2pi). This does
/// not matter for the kt_algorithm with
/// reasonable values of radius, but might matter for a general use
/// of this algorithm -- depending on whether or not the user has
/// initialised the class with instructions to ignore this problem the
/// program will detect and ignore it, or crash.
inline int Dnn2piCylinder::NearestNeighbourIndex(const int & current) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
  int plane_index;
  if (mirror_index == INEXISTENT_VERTEX ) {
    plane_index = _DNN->NearestNeighbourIndex(main_index);
  } else {
    plane_index = (
	_DNN->NearestNeighbourDistance(main_index) < 
	_DNN->NearestNeighbourDistance(mirror_index)) ? 
      _DNN->NearestNeighbourIndex(main_index) : 
      _DNN->NearestNeighbourIndex(mirror_index) ; 
  }
  int this_cylinder_index = _cylinder_index_of_plane_vertex[plane_index];
  // either the user has acknowledged the fact that they may get the
  // mirror copy as the closest point, or crash if it should occur
  // that mirror copy is the closest point.
  assert(_ignore_nearest_is_mirror || this_cylinder_index != current);
  //if (this_cylinder_index == current) {
  //  cerr << "WARNING point "<<current<<
  //    " has its mirror copy as its own nearest neighbour"<<endl;
  //}
  return this_cylinder_index;
}

inline double Dnn2piCylinder::NearestNeighbourDistance(const int & current) const {
  int main_index = _mirror_info[current].main_index;
  int mirror_index = _mirror_info[current].mirror_index;
  if (mirror_index == INEXISTENT_VERTEX ) {
    return _DNN->NearestNeighbourDistance(main_index);
  } else {
    return (
	_DNN->NearestNeighbourDistance(main_index) < 
	_DNN->NearestNeighbourDistance(mirror_index)) ? 
      _DNN->NearestNeighbourDistance(main_index) : 
      _DNN->NearestNeighbourDistance(mirror_index) ; 
  }
 
}

inline bool Dnn2piCylinder::Valid(const int & index) const {
  return (_DNN->Valid(_mirror_info[index].main_index));
}


inline Dnn2piCylinder::~Dnn2piCylinder() {
  delete _DNN; 
}


FASTJET_END_NAMESPACE

#endif //  __FASTJET_DNN2PICYLINDER_HH__
#endif //DROP_CGAL 
