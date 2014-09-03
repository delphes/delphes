// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: area.h                                                              //
// Description: header file for the computation of jet area                  //
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
// $Revision:: 149                                                          $//
// $Date:: 2007-03-15 00:13:58 +0100 (Thu, 15 Mar 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "area.h"
#include "momentum.h"
#include <stdlib.h>
#include <iostream>

namespace siscone{
using namespace std;

/*******************************************************
 * Cjet_area implementation                            *
 * real Jet information, including its area(s)         *
 *                                                     *
 * This class contains information for one single jet. *
 * That is, first, its momentum carrying information   *
 * about its centre and pT, and second, its particle   *
 * contents (from CJeT).                               *
 * Compared to the Cjet class, it also includes the    *
 * passive and active areas of the jet computed using  *
 * the Carea class.                                    *
 *******************************************************/

// default ctor
//--------------
Cjet_area::Cjet_area(){
  active_area = passive_area = 0.0;
}

// jet-initiated ctor
//-------------------
Cjet_area::Cjet_area(Cjet &j){
  v = j.v;
  n = j.n;
  contents = j.contents;

  pass = j.pass;

  pt_tilde = j.pt_tilde;
  sm_var2 = j.sm_var2;

  active_area = passive_area = 0.0;
}

// default dtor
//--------------
Cjet_area::~Cjet_area(){

}


/******************************************************************
 * Csiscone_area implementation                                   *
 * class for the computation of jet areas.                        *
 *                                                                *
 * This is the class user should use whenever you want to compute *
 * the jet area (passive and active).                             *
 * It uses the SISCone algorithm to perform the jet analysis.     *
 ******************************************************************/

// default ctor
//-------------
Carea::Carea(){
  grid_size = 60;     // 3600 particles added
  grid_eta_max = 6.0; // maybe having twice more points in eta than in phi should be nice
  grid_shift = 0.5;   // 50% "shacking"

  pt_soft = 1e-100;
  pt_shift = 0.05;
  pt_soft_min = 1e-90;
}

// default dtor
//-------------
Carea::~Carea(){
  
}
  
/*
 * compute the jet areas from a given particle set.
 * The parameters of this method are the ones which control the jet clustering alghorithm.
 * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
 * particles/jets and thus is used internally.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _f           shared energy threshold for splitting&merging
 *  - _n_pass_max  maximum number of passes (0=full search, the default)
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *        NOTE: SM_pt leads to IR unsafety for some events with momentum conservation. 
 *              SM_Et is IR safe but not boost invariant and not implemented(!)
 *              SM_mt is IR safe for hadronic events, but not for decays of two 
 *                    back-to-back particles of identical mass
 *              SM_pttilde  
 *                    is always IR safe, and also boost invariant (default)
 *  - _hard_only   when this is set on, only hard jets are computed
 *                 and not the purely ghosted jets (default: false)
 * return the jets together with their areas
 * The return value is the number of jets (including pure-ghost ones if they are included)
 ********************************************************************************************/
int Carea::compute_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
			 int _n_pass_max, Esplit_merge_scale _split_merge_scale,
			 bool _hard_only){

  vector<Cmomentum> all_particles;

  // put "hardest cut-off" if needed
  // this avoids computation of ghosted jets when not required and 
  // significantly shortens the SM.
  if (_hard_only){
    SM_var2_hardest_cut_off = pt_soft_min*pt_soft_min;
  } 

  // clear potential previous runs
  jet_areas.clear();

  // put initial set of particles in the list
  int n_hard = _particles.size();
  all_particles = _particles;

  // build the set of ghost particles and add them to the particle list
  int i,j;
  double eta_g,phi_g,pt_g;

  for (i=0;i<grid_size;i++){
    for (j=0;j<grid_size;j++){
      eta_g = grid_eta_max*(-1.0+2.0*(i+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      phi_g = M_PI        *(-1.0+2.0*(j+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      pt_g  = pt_soft*(1.0+pt_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))));
      all_particles.push_back(Cmomentum(pt_g*cos(phi_g),pt_g*sin(phi_g),pt_g*sinh(eta_g),pt_g*cosh(eta_g)));
    }
  }
  
  // run clustering with all particles.
  // the split-merge here dynamically accounts for the purely soft jets.
  // we therefore end up with the active area for the jets
  int n_jets = compute_jets(all_particles, _radius, _f, _n_pass_max, 0.0, _split_merge_scale);

  // save jets in the Cjet_area structure
  // and determine their size
  // jet contents is ordered by increasing index of the initial
  // particles. Hence, we look for the first particle with index >= n_hard
  // and deduce the number of ghosts in the jet, hence the area.
  double area_factor = (2.0*grid_eta_max/grid_size)*(twopi/grid_size);
  
  for (i=0;i<(int) jets.size();i++){
    jet_areas.push_back(jets[i]);
    j=0;
    while ((j<jets[i].n) && (jets[i].contents[j]<n_hard)) j++;
    jet_areas[i].active_area = (jets[i].n-j)*area_factor;
  }

  // determine passive jet area
  // for that onem we use the pt_min cut in order to remove purely
  // soft jets from the SM procedure
  recompute_jets(_f, pt_soft_min);

  // for the area computation, we assume the jete order is the same!
  for (i=0;i<(int) jets.size();i++){
    j=0;
    while ((j<jets[i].n) && (jets[i].contents[j]<n_hard)) j++;
    jet_areas[i].passive_area = (jets[i].n-j)*area_factor;
  }

  // Note:
  // there surely remain purely soft je at the end
  // their active area is 0 by default (see ctor)

  jets.clear();

  return n_jets;
}

/*
 * compute the passive jet areas from a given particle set.
 * The parameters of this method are the ones which control the jet clustering alghorithm.
 * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
 * particles/jets and thus is used internally.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _f           shared energy threshold for splitting&merging
 *  - _n_pass_max  maximum number of passes (0=full search, the default)
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *        NOTE: SM_pt leads to IR unsafety for some events with momentum conservation. 
 *              SM_Et is IR safe but not boost invariant and not implemented(!)
 *              SM_mt is IR safe for hadronic events, but not for decays of two 
 *                    back-to-back particles of identical mass
 *              SM_pttilde  
 *                    is always IR safe, and also boost invariant (default)
 * return the jets together with their passive areas
 * The return value is the number of jets
 ********************************************************************************************/
int Carea::compute_passive_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
				 int _n_pass_max, Esplit_merge_scale _split_merge_scale){

  vector<Cmomentum> all_particles;

  // in the case of passive area, we do not need 
  // to put the ghosts in the stable-cone search
  // (they do no influence the list of stable cones)
  // Here's how it goes...
  stable_cone_soft_pt2_cutoff = pt_soft_min*pt_soft_min;

  // clear potential previous runs
  jet_areas.clear();

  // put initial set of particles in the list
  int n_hard = _particles.size();
  all_particles = _particles;

  // build the set of ghost particles and add them to the particle list
  int i,j;
  double eta_g,phi_g,pt_g;

  for (i=0;i<grid_size;i++){
    for (j=0;j<grid_size;j++){
      eta_g = grid_eta_max*(-1.0+2.0*(i+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      phi_g = M_PI        *(-1.0+2.0*(j+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      pt_g  = pt_soft*(1.0+pt_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))));
      all_particles.push_back(Cmomentum(pt_g*cos(phi_g),pt_g*sin(phi_g),pt_g*sinh(eta_g),pt_g*cosh(eta_g)));
    }
  }
  
  // determine passive jet area
  // for that onem we use the pt_min cut in order to remove purely
  // soft jets from the SM procedure
  int n_jets = compute_jets(all_particles, _radius, _f, _n_pass_max, pt_soft_min, _split_merge_scale);

  // save jets in the Cjet_area structure
  // and determine their size
  // jet contents is ordered by increasing index of the initial
  // particles. Hence, we look for the first particle with index >= n_hard
  // and deduce the number of ghosts in the jet, hence the area.
  double area_factor = (2.0*grid_eta_max/grid_size)*(twopi/grid_size);
  for (i=0;i<(int) jets.size();i++){
    j=0;
    while ((j<jets[i].n) && (jets[i].contents[j]<n_hard)) j++;
    jet_areas[i].passive_area = (jets[i].n-j)*area_factor;
  }

  jets.clear();

  return n_jets;
}

/*
 * compute the active jet areas from a given particle set.
 * The parameters of this method are the ones which control the jet clustering alghorithm.
 * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
 * particles/jets and thus is used internally.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _f           shared energy threshold for splitting&merging
 *  - _n_pass_max  maximum number of passes (0=full search, the default)
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *        NOTE: SM_pt leads to IR unsafety for some events with momentum conservation. 
 *              SM_Et is IR safe but not boost invariant and not implemented(!)
 *              SM_mt is IR safe for hadronic events, but not for decays of two 
 *                    back-to-back particles of identical mass
 *              SM_pttilde  
 *                    is always IR safe, and also boost invariant (default)
 *  - _hard_only   when this is set on, only hard jets are computed
 *                 and not the purely ghosted jets (default: false)
 * return the jets together with their active areas
 * The return value is the number of jets (including pure-ghost ones if they are included)
 ********************************************************************************************/
int Carea::compute_active_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
				int _n_pass_max, Esplit_merge_scale _split_merge_scale,
				bool _hard_only){

  vector<Cmomentum> all_particles;

  // put "hardest cut-off" if needed
  // this avoids computation of ghosted jets when not required and 
  // significantly shortens the SM.
  if (_hard_only){
    SM_var2_hardest_cut_off = pt_soft_min*pt_soft_min;
  }

  // clear potential previous runs
  jet_areas.clear();

  // put initial set of particles in the list
  int n_hard = _particles.size();
  all_particles = _particles;

  // build the set of ghost particles and add them to the particle list
  int i,j;
  double eta_g,phi_g,pt_g;

  for (i=0;i<grid_size;i++){
    for (j=0;j<grid_size;j++){
      eta_g = grid_eta_max*(-1.0+2.0*(i+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      phi_g = M_PI        *(-1.0+2.0*(j+0.5+grid_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))))/grid_size);
      pt_g  = pt_soft*(1.0+pt_shift*(-1.0+2.0*(rand()/(RAND_MAX+1.0))));
      all_particles.push_back(Cmomentum(pt_g*cos(phi_g),pt_g*sin(phi_g),pt_g*sinh(eta_g),pt_g*cosh(eta_g)));
    }
  }
  
  // run clustering with all particles.
  // the split-merge here dynamically accounts for the purely soft jets.
  // we therefore end up with the active area for the jets
  int n_jets = compute_jets(all_particles, _radius, _f, _n_pass_max, 0.0, _split_merge_scale);

  // save jets in the Cjet_area structure
  // and determine their size
  // jet contents is ordered by increasing index of the initial
  // particles. Hence, we look for the first particle with index >= n_hard
  // and deduce the number of ghosts in the jet, hence the area.
  double area_factor = (2.0*grid_eta_max/grid_size)*(twopi/grid_size);
  
  for (i=0;i<(int) jets.size();i++){
    jet_areas.push_back(jets[i]);
    j=0;
    while ((j<jets[i].n) && (jets[i].contents[j]<n_hard)) j++;
    jet_areas[i].active_area = (jets[i].n-j)*area_factor;
  }

  jets.clear();

  return n_jets;
}

}
