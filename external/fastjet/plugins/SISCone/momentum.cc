///////////////////////////////////////////////////////////////////////////////
// File: momentum.cpp                                                        //
// Description: source file for 4-momentum class Cmomentum                   //
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

#include "momentum.h"
#include <math.h>
#include <stdlib.h>

namespace siscone{

/*************************************************************************
 * class Cmomentum                                                       *
 * This class contains the information for particle or group of          *
 * particles management.                                                 *
 * It includes all Lorentz properties as well as tools for summing them. *
 *************************************************************************/
 
// default ctor
//--------------
Cmomentum::Cmomentum(){
  eta = 0.0;
  phi = 0.0;
  px = py = pz = E = 0.0;
  ref = Creference();
  index = -1;
}

// ctor with initialisation
//--------------------------
Cmomentum::Cmomentum(double _px, double _py, double _pz, double _E){
  px = _px;
  py = _py;
  pz = _pz;
  E  = _E;

  // compute eta and phi
  build_etaphi();
  ref = Creference();
}

// ctor with detailed initialisation
//-----------------------------------
Cmomentum::Cmomentum(double _eta, double _phi, Creference _ref){
  eta = _eta;
  phi = _phi;

  ref = _ref;
}

// default dtor
//--------------
Cmomentum::~Cmomentum(){

}

// assignment of vectors
//-----------------------
Cmomentum& Cmomentum::operator = (const Cmomentum &v){
  px = v.px;
  py = v.py;
  pz = v.pz;
  E  = v.E;

  eta = v.eta;
  phi = v.phi;

  ref = v.ref;
  return *this;
}

// addition of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
const Cmomentum Cmomentum::operator + (const Cmomentum &v){
  Cmomentum tmp = *this;
  return tmp+=v;
}

// incrementation of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
Cmomentum& Cmomentum::operator += (const Cmomentum &v){
  px+=v.px;
  py+=v.py;
  pz+=v.pz;
  E +=v.E;

  ref+=v.ref;

  return *this;
}

// incrementation of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
Cmomentum& Cmomentum::operator -= (const Cmomentum &v){
  px-=v.px;
  py-=v.py;
  pz-=v.pz;
  E -=v.E;

  ref-=v.ref;
  return *this;
}

// build eta-phi from 4-momentum info
// !!!                WARNING                   !!!
// !!! computing eta and phi is time-consuming  !!!
// !!! use this whenever you need eta or phi    !!!
// !!! automatically called for single-particle !!!
//--------------------------------------------------
void Cmomentum::build_etaphi(){
  // note: the factor n (ref.nb) cancels in all expressions !!
  eta = 0.5*log((E+pz)/(E-pz));
  phi = atan2(py,px);
}


// ordering of two vectors
// the default ordering is w.r.t. their references
//-------------------------------------------------
bool operator < (const Cmomentum &v1, const Cmomentum &v2){
  return v1.ref < v2.ref;
}

// ordering of vectors in eta (e.g. used in collinear tests)
//-----------------------------------------------------------
bool momentum_eta_less(const Cmomentum &v1, const Cmomentum &v2){
  return v1.eta < v2.eta;
}

// ordering of vectors in pt
//---------------------------
bool momentum_pt_less(const Cmomentum &v1, const Cmomentum &v2){
  return v1.perp2() < v2.perp2();
}

}

