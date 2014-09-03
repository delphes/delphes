///////////////////////////////////////////////////////////////////////////////
// File: geom_2d.cpp                                                         //
// Description: source file for two-dimensional geometry tools               //
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
// $Revision:: 171                                                          $//
// $Date:: 2007-06-19 16:26:05 +0200 (Tue, 19 Jun 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "geom_2d.h"
#include <algorithm>

namespace siscone{

#define PHI_RANGE_MASK 0xFFFFFFFF

/*********************************************************
 * class Ceta_phi_range implementation                   *
 * class for holding a covering range in eta-phi         *
 *                                                       *
 * This class deals with ranges in the eta-phi plane. It *
 * implements methods to test if two ranges overlap and  *
 * to take the union of two overlapping intervals.       *
 *********************************************************/

using namespace std;

// static member default init
//----------------------------
double Ceta_phi_range::eta_min = -100.0;
double Ceta_phi_range::eta_max = 100.0;

// default ctor
//--------------
Ceta_phi_range::Ceta_phi_range(){
  eta_range = 0;
  phi_range = 0;
}

// ctor with initialisation
// we initialise with a centre (in eta,phi) and a radius
//  - c_eta   eta coordinate of the centre
//  - c_phi   phi coordinate of the centre
//  - R       radius
//-------------------------------------------------------
Ceta_phi_range::Ceta_phi_range(double c_eta, double c_phi, double R){
  // determination of the eta range
  //-------------------------------
  double xmin = max(c_eta-R,eta_min+0.0001);
  double xmax = min(c_eta+R,eta_max-0.0001);

  unsigned int cell_min = get_eta_cell(xmin);
  unsigned int cell_max = get_eta_cell(xmax);

  // warning: if cell_max==2^31, 2*cell_max==0 hence, 
  // even if the next formula is formally (2*cell_max-cell_min),
  // expressing it as (cell_max-cell_min)+cell_max is safe.
  eta_range = (cell_max-cell_min)+cell_max;

  // determination of the phi range
  // !! taking care of periodicity !!
  //---------------------------------
  xmin = phi_in_range(c_phi-R);
  xmax = phi_in_range(c_phi+R);

  cell_min = get_phi_cell(xmin);
  cell_max = get_phi_cell(xmax);

  // Also, if the interval goes through pi, inversion is needed
  if (xmax>xmin)
    phi_range = (cell_max-cell_min)+cell_max;
  else {
    phi_range = (cell_min==cell_max) 
      ? PHI_RANGE_MASK
      : ((PHI_RANGE_MASK^(cell_min-cell_max)) + cell_max);
  }
}

// assignment of range
//  - r   range to assign to current one
//---------------------------------------
Ceta_phi_range& Ceta_phi_range::operator = (const Ceta_phi_range &r){
  eta_range = r.eta_range;
  phi_range = r.phi_range;

  return *this;
}

// add a particle to the range
//  - eta  eta coordinate of the particle
//  - phi  phi coordinate of the particle
// \return 0 on success, 1 on error
//----------------------------------------
int Ceta_phi_range::add_particle(const double eta, const double phi){
  // deal with the eta coordinate
  eta_range |= get_eta_cell(eta);

  // deal with the phi coordinate
  phi_range |= get_phi_cell(phi);

  return 0;
}


// test overlap
//  - r1  first range
//  - r2  second range
// return true if overlap, false otherwise.
//------------------------------------------
bool is_range_overlap(const Ceta_phi_range &r1, const Ceta_phi_range &r2){
  // check overlap in eta AND phi
  return ((r1.eta_range & r2.eta_range) && (r1.phi_range & r2.phi_range));
}

// compute union
// Note: we assume that the two intervals overlap
//  - r1  first range
//  - r2  second range
// \return union of the two ranges
//------------------------------------------
const Ceta_phi_range range_union (const Ceta_phi_range &r1, const Ceta_phi_range &r2){
  Ceta_phi_range tmp;

  // compute union in eta
  tmp.eta_range = r1.eta_range | r2.eta_range;

  // compute union in phi
  tmp.phi_range = r1.phi_range | r2.phi_range;

  return tmp;
}

}
