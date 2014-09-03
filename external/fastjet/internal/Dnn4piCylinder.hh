//FJSTARTHEADER
// $Id: Dnn4piCylinder.hh 3442 2014-07-24 07:20:49Z salam $
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
#ifndef __FASTJET_DNN4PICYLINDER_HH__
#define __FASTJET_DNN4PICYLINDER_HH__

#include "fastjet/internal/DynamicNearestNeighbours.hh"
#include "fastjet/internal/DnnPlane.hh"
#include "fastjet/internal/numconsts.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// \if internal_doc
/// @ingroup internal
/// \class Dnn4piCylinder
/// class derived from DynamicNearestNeighbours that provides an
/// implementation for the surface of cylinder (using two copies of
/// DnnPlane, one running from 0--2pi, the other from pi--3pi).
/// \endif
class Dnn4piCylinder : public DynamicNearestNeighbours {
 public:
  /// empty initaliser
  Dnn4piCylinder() {}

  /// Initialiser from a set of points on an Eta-Phi plane, where
  /// eta can have an arbitrary ranges and phi must be in range
  /// 0 <= phi < 2pi
  Dnn4piCylinder(const std::vector<EtaPhi> &, const bool & verbose = false );

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

  ~Dnn4piCylinder();

 private:

  bool _verbose;

  // NB: we define POINTERS here because the initialisation gave
  //     us problems (things crashed!), perhaps because in practice
  //     we were making a copy without being careful and defining
  //     a proper copy constructor.
  DnnPlane * _DNN1, * _DNN2;

  /// given a phi value in the 0--2pi range return one 
  /// in the pi--3pi range.
  inline EtaPhi _remap_phi(const EtaPhi & point) {
    double phi = point.second;
    if (phi < pi) { phi += twopi ;}
    return EtaPhi(point.first, phi);}

};


// here follow some inline implementations of the simpler of the
// functions defined above

inline int Dnn4piCylinder::NearestNeighbourIndex(const int current) const {
  return (_DNN1->NearestNeighbourDistance(current) < 
	  _DNN2->NearestNeighbourDistance(current)) ? 
    _DNN1->NearestNeighbourIndex(current) : 
    _DNN2->NearestNeighbourIndex(current) ; 
}

inline double Dnn4piCylinder::NearestNeighbourDistance(const int current) const {
  return (_DNN1->NearestNeighbourDistance(current) < 
	  _DNN2->NearestNeighbourDistance(current)) ? 
    _DNN1->NearestNeighbourDistance(current) : 
    _DNN2->NearestNeighbourDistance(current) ; 
}

inline bool Dnn4piCylinder::Valid(const int index) const {
  return (_DNN1->Valid(index) && _DNN2->Valid(index));
}


inline Dnn4piCylinder::~Dnn4piCylinder() {
  delete _DNN1; 
  delete _DNN2;
}


FASTJET_END_NAMESPACE

#endif //  __FASTJET_DNN4PICYLINDER_HH__
#endif //  DROP_CGAL 
