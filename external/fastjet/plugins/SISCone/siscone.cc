///////////////////////////////////////////////////////////////////////////////
// File: siscone.cpp                                                         //
// Description: source file for the main SISCone class                       //
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
// $Revision:: 371                                                          $//
// $Date:: 2014-09-09 10:05:32 +0200 (Tue, 09 Sep 2014)                     $//
///////////////////////////////////////////////////////////////////////////////

//#ifdef HAVE_CONFIG_H
#include "config.h"
//#else
//#define PACKAGE_NAME "SISCone"
//#define VERSION "3.0.0"
//#warning "No config.h file available, using preset values"
//#endif

#include "ranlux.h"
#include "momentum.h"
#include "defines.h"
#include "siscone.h"
#include "siscone_error.h"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace siscone{
using namespace std;

/***************************************************************
 * Csiscone implementation                                     *
 * final class: gather everything to compute the jet contents. *
 *                                                             *
 * This is the class user should use.                          *
 * It computes the jet contents of a list of particles         *
 * given a cone radius and a threshold for splitting/merging.  *
 ***************************************************************/

// default ctor
//--------------
Csiscone::Csiscone(){
  rerun_allowed = false;
}

// default dtor
//--------------
Csiscone::~Csiscone(){
  rerun_allowed = false;
}

bool Csiscone::init_done=false;
std::ostream* Csiscone::_banner_ostr = &cout;

/*
 * compute the jets from a given particle set doing multiple passes
 * such pass N looks for jets among all particles not put into jets
 * during previous passes.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _f           shared energy threshold for splitting&merging
 *  - _n_pass_max  maximum number of runs
 *  - _ptmin       minimum pT of the protojets
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *    NOTE: using pt leads to IR unsafety for some events with momentum
 *          conservation. So we strongly advise not to change the default
 *          value.
 * return the number of jets found.
 **********************************************************************/
int Csiscone::compute_jets(vector<Cmomentum> &_particles, double _radius, double _f, 
			   int _n_pass_max, double _ptmin,
			   Esplit_merge_scale _split_merge_scale){
  _initialise_if_needed();

  // run some general safety tests (NB: f will be checked in split-merge)
  if (_radius <= 0.0 || _radius >= 0.5*M_PI) {
    ostringstream message;
    message << "Illegal value for cone radius, R = " << _radius 
            << " (legal values are 0<R<pi/2)";
    throw Csiscone_error(message.str());
  }



  ptcomparison.split_merge_scale = _split_merge_scale;
  partial_clear(); // make sure some things are initialised properly

  // init the split_merge algorithm with the initial list of particles
  // this initialises particle list p_left of remaining particles to deal with
  init_particles(_particles);

  bool finished = false;

  rerun_allowed = false;
  protocones_list.clear();

#ifdef DEBUG_STABLE_CONES
  nb_hash_cones_total = 0;
  nb_hash_occupied_total = 0;
#endif

  do{
    // initialise stable_cone finder
    // here we use the list of remaining particles
    // AFTER COLLINEAR CLUSTERING !!!!!!
    Cstable_cones::init(p_uncol_hard);

    // get stable cones
    if (get_stable_cones(_radius)){
      // we have some new protocones; add them to candidates
      // Note that add_protocones has to be called first
      // if we want the 4-vect components to be available
      // on top of eta and phi.
      add_protocones(&protocones, R2, _ptmin);
      protocones_list.push_back(protocones);
#ifdef DEBUG_STABLE_CONES
      nb_hash_cones_total += nb_hash_cones;
      nb_hash_occupied_total += nb_hash_occupied;
#endif
    } else {
      // no new protocone: leave
      finished=true;
    }

    _n_pass_max--;
  } while ((!finished) && (n_left>0) && (_n_pass_max!=0));

  rerun_allowed = true;

  // split & merge
  return perform(_f, _ptmin);
}


/*
 * compute the jets from a given particle set doing multiple passes
 * such pass N looks for jets among all particles not put into jets
 * during previous passes.
 *  - _particles   list of particles
 *  - _radius      cone radius
 *  - _n_pass_max  maximum number of runs
 *  - _ptmin       minimum pT of the protojets
 *  - _ordering_scale    the ordering scale to decide which stable
 *                       cone is removed
 * return the number of jets found.
 **********************************************************************/
int Csiscone::compute_jets_progressive_removal(vector<Cmomentum> &_particles, double _radius, 
					       int _n_pass_max, double _ptmin,
					       Esplit_merge_scale _ordering_scale){
  _initialise_if_needed();

  // run some general safety tests (NB: f will be checked in split-merge)
  if (_radius <= 0.0 || _radius >= 0.5*M_PI) {
    ostringstream message;
    message << "Illegal value for cone radius, R = " << _radius 
            << " (legal values are 0<R<pi/2)";
    throw Csiscone_error(message.str());
  }

  ptcomparison.split_merge_scale = _ordering_scale;
  partial_clear(); // make sure some things are initialised properly

  // init the split_merge algorithm with the initial list of particles
  // this initialises particle list p_left of remaining particles to deal with
  //
  // this stores the "processed" particles in p_uncol_hard
  init_particles(_particles);
  jets.clear();

  bool unclustered_left;
  rerun_allowed = false;
  protocones_list.clear();

  do{
    //cout << n_left << " particle left" << endl; 

    // initialise stable_cone finder
    // here we use the list of remaining particles
    // AFTER COLLINEAR CLUSTERING !!!!!!
    Cstable_cones::init(p_uncol_hard);

    // get stable cones (stored in 'protocones')
    unclustered_left = get_stable_cones(_radius);

    // add the hardest stable cone to the list of jets
    if (add_hardest_protocone_to_jets(&protocones, R2, _ptmin)) break;
  
    _n_pass_max--;
  } while ((unclustered_left) && (n_left>0) && (_n_pass_max!=0));

  // split & merge
  return jets.size();
}


/*
 * recompute the jets with a different overlap parameter.
 * we use the same particles and R as in the preceeding call.
 *  - _f           shared energy threshold for splitting&merging
 *  - _ptmin       minimum pT of the protojets
 *  - _split_merge_scale    the scale choice for the split-merge procedure
 *    NOTE: using pt leads to IR unsafety for some events with momentum
 *          conservation. So we strongly advise not to change the default
 *          value.
 * return the number of jets found, -1 if recomputation not allowed.
 ********************************************************************/
int Csiscone::recompute_jets(double _f, double _ptmin,
			     Esplit_merge_scale _split_merge_scale){
  if (!rerun_allowed)
    return -1;

  ptcomparison.split_merge_scale = _split_merge_scale;

  // restore particle list
  partial_clear();
  init_pleft();

  // initialise split/merge algorithm
  unsigned int i;
  for (i=0;i<protocones_list.size();i++)
    add_protocones(&(protocones_list[i]), R2, _ptmin);

  // split & merge
  return perform(_f, _ptmin);
}  

// ensure things are initialised
void Csiscone::_initialise_if_needed(){
  // initialise random number generator
  if (init_done) return;

  // initialise random number generator
  ranlux_init();

  // do not do this again
  init_done=true;

  // print the banner
  if (_banner_ostr != 0){
    (*_banner_ostr) << "#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
    (*_banner_ostr) << "#                    SISCone   version " << setw(28) << left << siscone_version() << "o" << endl;
    (*_banner_ostr) << "#              http://projects.hepforge.org/siscone                o" << endl;
    (*_banner_ostr) << "#                                                                  o" << endl;
    (*_banner_ostr) << "# This is SISCone: the Seedless Infrared Safe Cone Jet Algorithm   o" << endl;
    (*_banner_ostr) << "# SISCone was written by Gavin Salam and Gregory Soyez             o" << endl;
    (*_banner_ostr) << "# It is released under the terms of the GNU General Public License o" << endl;
    (*_banner_ostr) << "#                                                                  o" << endl;
    (*_banner_ostr) << "# A description of the algorithm is available in the publication   o" << endl;
    (*_banner_ostr) << "# JHEP 05 (2007) 086 [arXiv:0704.0292 (hep-ph)].                   o" << endl;
    (*_banner_ostr) << "# Please cite it if you use SISCone.                               o" << endl;
    (*_banner_ostr) << "#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo" << endl;
    (*_banner_ostr) << endl;

    _banner_ostr->flush();
  }
}

// finally, a bunch of functions to access to 
// basic information (package name, version)
//---------------------------------------------

/* 
 * return SISCone package name.
 * This is nothing but "SISCone", it is a replacement to the
 * PACKAGE_NAME string defined in config.h and which is not
 * public by default.
 * return the SISCone name as a string
 */
string siscone_package_name(){
  return PACKAGE_NAME;
}

/* 
 * return SISCone version number.
 * return a string of the form "X.Y.Z" with possible additional tag
 *        (alpha, beta, devel) to mention stability status
 */
string siscone_version(){
  return VERSION;
}

}
