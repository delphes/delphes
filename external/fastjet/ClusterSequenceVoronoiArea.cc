//FJSTARTHEADER
// $Id: ClusterSequenceVoronoiArea.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2006-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/ClusterSequenceVoronoiArea.hh"
#include "fastjet/internal/Voronoi.hh"
#include <list>
#include <cassert>
#include <ostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <limits>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

typedef ClusterSequenceVoronoiArea::VoronoiAreaCalc VAC;

/// class for carrying out a voronoi area calculation on a set of
/// initial vectors
class ClusterSequenceVoronoiArea::VoronoiAreaCalc {
public:
  /// constructor that takes a range of a vector together with the
  /// effective radius for the intersection of discs with voronoi
  /// cells
  VoronoiAreaCalc(const vector<PseudoJet>::const_iterator &,
		  const vector<PseudoJet>::const_iterator &,
		  double effective_R);

  /// return the area of the particle associated with the given
  /// index
  inline double area (int index) const {return _areas[index];};

private:
  std::vector<double> _areas;     ///< areas, numbered as jets
  double _effective_R;            ///< effective radius
  double _effective_R_squared;    ///< effective radius squared

  /**
   * compute the intersection of one triangle with the circle
   * the area is returned
   */
  double edge_circle_intersection(const VPoint &p0,
				  const GraphEdge &edge);

  /// get the area of a circle of radius R centred on the point 0 with
  /// 1 and 2 on each "side" of the arc. dij is the distance between
  /// point i and point j and all distances are squared
  inline double circle_area(const double d12_2, double d01_2, double d02_2){
    return 0.5*_effective_R_squared
      *acos(min(1.0,(d01_2+d02_2-d12_2)/(2*sqrt(d01_2*d02_2))));
  }
};


/**
 * compute the intersection of one triangle with the circle
 * the area is returned
 */
double VAC::edge_circle_intersection(const VPoint &p0,
				     const GraphEdge &edge){
  VPoint p1(edge.x1-p0.x, edge.y1-p0.y);
  VPoint p2(edge.x2-p0.x, edge.y2-p0.y);
  VPoint pdiff = p2-p1;

  //fprintf(stdout, "\tpt(%f,%f)\n", p0.x, p0.y);

  double cross = vector_product(p1, p2);
  double d12_2 = norm(pdiff);
  double d01_2 = norm(p1);
  double d02_2 = norm(p2);

  // compute intersections between edge line and circle
  double delta = d12_2*_effective_R_squared - cross*cross;
  
  // if no intersection, area=area_circle
  if (delta<=0){
    return circle_area(d12_2, d01_2, d02_2);
  }

  // we'll only need delta's sqrt now
  delta = sqrt(delta);

  // b is the projection of 01 onto 12
  double b = scalar_product(pdiff, p1);

  // intersections with the circle:
  //   we compute the "coordinate along the line" of the intersection
  //   with t=0 (1) corresponding to p1 (p2)
  // points with 0<t<1 are within the circle others are outside

  // positive intersection
  double tp = (delta-b)/d12_2;

  // if tp is negative, tm also => inters = circle
  if (tp<0)
    return circle_area(d12_2, d01_2, d02_2);

  // we need the second intersection
  double tm = -(delta+b)/d12_2;

  // if tp<1, it lies in the circle
  if (tp<1){
    // if tm<0, the segment has one intersection
    // with the circle at p (t=tp)
    // the area is a triangle from 1 to p
    //        then a circle   from p to 2
    // several tricks can be used:
    //  - the area of the triangle is tp*area triangle
    //  - the lenght for the circle are easily obtained
    if (tm<0)
      return tp*0.5*fabs(cross)
        +circle_area((1-tp)*(1-tp)*d12_2, _effective_R_squared, d02_2);

    // now, 0 < tm < tp < 1
    // the segment intersects twice the circle
    //   area = 2 cirles at ends + a triangle in the middle
    // again, simplifications are staightforward
    return (tp-tm)*0.5*fabs(cross)
      + circle_area(tm*tm*d12_2, d01_2, _effective_R_squared)
      + circle_area((1-tp)*(1-tp)*d12_2, _effective_R_squared, d02_2);
  }

  // now, we have tp>1

  // if in addition tm>1, intersectino is a circle
  if (tm>1)
    return circle_area(d12_2, d01_2, d02_2);

  // if tm<0, the triangle is inside the circle
  if (tm<0)
    return 0.5*fabs(cross);

  // otherwise, only the "tm point" is on the segment
  //   area = circle from 1 to m and triangle from m to 2

  return (1-tm)*0.5*fabs(cross)
    +circle_area(tm*tm*d12_2, d01_2, _effective_R_squared);
}


// the constructor...
//----------------------------------------------------------------------
VAC::VoronoiAreaCalc(const vector<PseudoJet>::const_iterator &jet_begin,
		     const vector<PseudoJet>::const_iterator &jet_end,
		     double effective_R) {

  assert(effective_R < 0.5*pi);

  vector<VPoint> voronoi_particles;
  vector<int> voronoi_indices;

  _effective_R         = effective_R;
  _effective_R_squared = effective_R*effective_R;

  double minrap = numeric_limits<double>::max();
  double maxrap = -minrap;

  unsigned int n_tot = 0, n_added = 0;

  // loop over jets and create the triangulation, as well as cross-referencing
  // info
  for (vector<PseudoJet>::const_iterator jet_it = jet_begin; 
       jet_it != jet_end; jet_it++) {
    _areas.push_back(0.0);
    if ((jet_it->perp2()) != 0.0 || (jet_it->E() != jet_it->pz())){
      // generate the corresponding point
      double rap = jet_it->rap(), phi = jet_it->phi();
      voronoi_particles.push_back(VPoint(rap, phi));
      voronoi_indices.push_back(n_tot);
      n_added++;

      // insert a copy of the point if it falls within 2*_R_effective
      // of the 0,2pi borders (because we are interested in any
      // voronoi edge within _R_effective of the other border)
      if (phi < 2*_effective_R) {
	voronoi_particles.push_back(VPoint(rap,phi+twopi));
	voronoi_indices.push_back(-1);
	n_added++;
      } else if (twopi-phi < 2*_effective_R) {
	voronoi_particles.push_back(VPoint(rap,phi-twopi));
	voronoi_indices.push_back(-1);
	n_added++;
      }

      // track the rapidity range
      maxrap = max(maxrap,rap);
      minrap = min(minrap,rap);
    }
    n_tot++;
  }

  // allow for 0-particle case in graceful way
  if (n_added == 0) return;
  // assert(n_added > 0); // old (pre 2.4) non-graceful exit

  // add extreme cases (corner particles):
  double max_extend = 2*max(maxrap-minrap+4*_effective_R, twopi+8*_effective_R);
  voronoi_particles.push_back(VPoint(0.5*(minrap+maxrap)-max_extend, pi));
  voronoi_particles.push_back(VPoint(0.5*(minrap+maxrap)+max_extend, pi));
  voronoi_particles.push_back(VPoint(0.5*(minrap+maxrap), pi-max_extend));
  voronoi_particles.push_back(VPoint(0.5*(minrap+maxrap), pi+max_extend));

  // Build the VD
  VoronoiDiagramGenerator vdg;
  vdg.generateVoronoi(&voronoi_particles, 
		      0.5*(minrap+maxrap)-max_extend, 0.5*(minrap+maxrap)+max_extend,
		      pi-max_extend, pi+max_extend);

  vdg.resetIterator();
  GraphEdge *e=NULL;
  unsigned int v_index;
  int p_index;
  vector<PseudoJet>::const_iterator jet;

  while(vdg.getNext(&e)){
    v_index = e->point1;
    if (v_index<n_added){ // this removes the corner particles
      p_index = voronoi_indices[v_index];
      if (p_index!=-1){   // this removes the copies
	jet = jet_begin+voronoi_indices[v_index];
	_areas[p_index]+=
	  edge_circle_intersection(voronoi_particles[v_index], *e);
      }
    }
    v_index = e->point2;
    if (v_index<n_added){ // this removes the corner particles
      p_index = voronoi_indices[v_index];
      if (p_index!=-1){   // this removes the copies
	jet = jet_begin+voronoi_indices[v_index];
	_areas[p_index]+=
	  edge_circle_intersection(voronoi_particles[v_index], *e);
      }
    }
  }


}


//----------------------------------------------------------------------
///
void ClusterSequenceVoronoiArea::_initializeVA () {
  // run the VAC on our original particles
  _pa_calc = new VAC(_jets.begin(), 
		     _jets.begin()+n_particles(),
		     _effective_Rfact*_jet_def.R());

  // transfer the areas to our local structure
  //  -- first the initial ones
  _voronoi_area.reserve(2*n_particles());
  _voronoi_area_4vector.reserve(2*n_particles());
  for (unsigned int i=0; i<n_particles(); i++) {
    _voronoi_area.push_back(_pa_calc->area(i));
    // make a stab at a 4-vector area
    if (_jets[i].perp2() > 0) {
      _voronoi_area_4vector.push_back((_pa_calc->area(i)/_jets[i].perp())
                                      * _jets[i]);
    } else {
      // not sure what to do here -- just put zero (it won't be meaningful
      // anyway)
      _voronoi_area_4vector.push_back(PseudoJet(0.0,0.0,0.0,0.0));
    }
  }
	   
  //  -- then the combined areas that arise from the clustering
  for (unsigned int i = n_particles(); i < _history.size(); i++) {
    double area_local;
    PseudoJet area_4vect;
    if (_history[i].parent2 >= 0) {
      area_local = _voronoi_area[_history[i].parent1] + 
  	           _voronoi_area[_history[i].parent2];
      area_4vect = _voronoi_area_4vector[_history[i].parent1] + 
                   _voronoi_area_4vector[_history[i].parent2];
    } else {
      area_local = _voronoi_area[_history[i].parent1];
      area_4vect = _voronoi_area_4vector[_history[i].parent1];
    }
    _voronoi_area.push_back(area_local);
    _voronoi_area_4vector.push_back(area_4vect);
  }

}

//----------------------------------------------------------------------
ClusterSequenceVoronoiArea::~ClusterSequenceVoronoiArea() {
  delete _pa_calc;
}

FASTJET_END_NAMESPACE
