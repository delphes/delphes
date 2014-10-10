// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: siscone.h                                                           //
// Description: header file for the main SISCone class                       //
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
// $Revision:: 369                                                          $//
// $Date:: 2014-09-04 16:57:55 +0200 (Thu, 04 Sep 2014)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SISCONE_H__
#define __SISCONE_H__

#include "protocones.h"
#include "split_merge.h"

namespace siscone{

/**
 * \class Csiscone
 * final class: gather everything to compute the jet contents.
 * 
 * This is the class user should use.
 * It computes the jet contents of a list of particles
 * given a cone radius and a threshold for splitting/merging.
 *
 * After the call to 'perform', the vector jets is filled with
 * the jets found. the 'contents' field of each jets contains
 * the indices of the particles included in that jet. 
 */
class Csiscone : public Cstable_cones, public Csplit_merge{
 public:
  /// default ctor
  Csiscone();

  /// default dtor
  ~Csiscone();

  /**
   * compute the jets from a given particle set.
   * We are doing multiple passes such pass n_pass looks for jets among 
   * all particles not put into jets during previous passes.
   * By default the number of passes is infinite (0). 
   * \param _particles   list of particles
   * \param _radius      cone radius
   * \param _f           shared energy threshold for splitting&merging
   * \param _n_pass_max  maximum number of passes (0=full search)
   * \param _ptmin       minimum pT of the protojets
   * \param _split_merge_scale    the scale choice for the split-merge procedure
   *        NOTE: SM_pt leads to IR unsafety for some events with momentum conservation. 
   *              SM_Et is IR safe but not boost invariant and not implemented(!)
   *              SM_mt is IR safe for hadronic events, but not for decays of two 
   *                    back-to-back particles of identical mass
   *              SM_pttilde  
   *                    is always IR safe, and also boost invariant (default)
   *
   * \return the number of jets found.
   */
  int compute_jets(std::vector<Cmomentum> &_particles, double _radius, double _f, 
		   int _n_pass_max=0, double _ptmin=0.0,
		   Esplit_merge_scale _split_merge_scale=SM_pttilde);

  /**
   * compute the jets from a given particle set.
   * We are doing multiple passes such pass n_pass looks for jets among 
   * all particles not put into jets during previous passes.
   * By default the number of passes is infinite (0). 
   * \param _particles   list of particles
   * \param _radius      cone radius
   * \param _n_pass_max  maximum number of passes (0=full search)
   * \param _ptmin       minimum pT of the protojets
   * \param _ordering_scale    the ordering scale to decide which stable
   *                           cone is removed
   *
   * Note that the Csplit_merge::SM_var2_hardest_cut_off cut is not
   * used in the progressive removal variant.
   * 
   * \return the number of jets found.
   */
  int compute_jets_progressive_removal(std::vector<Cmomentum> &_particles, double _radius, 
				       int _n_pass_max=0, double _ptmin=0.0,
				       Esplit_merge_scale _ordering_scale=SM_pttilde);

  /**
   * recompute the jets with a different overlap parameter.
   * we use the same particles and R as in the preceeding call.
   * \param _f           shared energy threshold for splitting&merging
   * \param _ptmin       minimum pT of the protojets
   * \param _split_merge_scale    the scale choice for the split-merge procedure
   *                                           split--merge variable
   *        NOTE: using pt leads to IR unsafety for some events with momentum
   *              conservation. So we strongly advise not to change the default
   *              value.
   * \return the number of jets found, -1 if recomputation not allowed.
   */
  int recompute_jets(double _f, double _ptmin = 0.0,
		     Esplit_merge_scale _split_merge_scale=SM_pttilde);

  /// list of protocones found pass-by-pass (not filled by compute_jets_progressive_removal)
  std::vector<std::vector<Cmomentum> > protocones_list;

  // random number initialisation
  static bool init_done;      ///< check random generator initialisation

#ifdef DEBUG_STABLE_CONES
  int nb_hash_cones_total, nb_hash_occupied_total;
#endif

  /**
   * A call to this function modifies the stream
   * used to print banners (by default cout).
   *
   * Please note that if you distribute 3rd party code
   * that links with SISCone, that 3rd party code must not
   * use this call turn off the printing of thw banner
   * by default. This requirement reflects the spirit of
   * clause 2c of the GNU Public License (v2), under which
   * SISCone is distributed.
   */
  static void set_banner_stream(std::ostream * ostr) {_banner_ostr = ostr;}

  /**
   * returns a pointer to the stream to be used to print banners
   * (cout by default)
   */
  static std::ostream * banner_stream() {return _banner_ostr;}

 private:
  bool rerun_allowed;         ///< is recompute_jets allowed ?
  static std::ostream * _banner_ostr; ///< stream to use for banners

  /// ensure things are initialised
  void _initialise_if_needed();

};

  
// finally, a bunch of functions to access to 
// basic information (package name, version)
//---------------------------------------------

/** 
 * return SISCone package name.
 * This is nothing but "SISCone", it is a replacement to the
 * PACKAGE_NAME string defined in config.h and which is not
 * public by default.
 * \return the SISCone name as a string
 */
std::string siscone_package_name();

/** 
 * return SISCone version number.
 * \return a string of the form "X.Y.Z" with possible additional tag
 *         (alpha, beta, devel) to mention stability status
 */
std::string siscone_version();

}
#endif
