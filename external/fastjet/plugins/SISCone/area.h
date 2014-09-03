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

#ifndef __SISCONE_AREA_H__
#define __SISCONE_AREA_H__

#include "momentum.h"
#include "siscone.h"

namespace siscone{

/**
 * \class Cjet_area
 * real Jet information, including its area(s)
 *
 * This class contains information for one single jet. 
 * That is, first, its momentum carrying information
 * about its centre and pT, and second, its particle
 * contents.
 * Compared to the Cjet class, it also includes the 
 * passive and active areas of the jet computed using 
 * the Carea class.
 */
class Cjet_area : public Cjet{
 public:
  /// default ctor
  Cjet_area();

  /// jet-initialised ctor
  Cjet_area(Cjet &j);

  /// default dtor
  ~Cjet_area();

  // area information
  double passive_area;   ///< passive area
  double active_area;    ///< active area
};

/**
 * \class Carea
 * class for the computation of jet areas.
 * 
 * This is the class user should use whenever you want to compute
 * the jet area (passive and active). .
 * It uses the SISCone algorithm to perform the jet analysis.
 */
class Carea : public Csiscone{
 public:
  /// default ctor
  Carea();

  /// default dtor
  ~Carea();

  /**
   * compute the jet areas from a given particle set.
   * The parameters of this method are the ones which control the jet clustering alghorithn.
   * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
   * particles/jets and thus is used internally.
   * \param _particles   list of particles
   * \param _radius      cone radius
   * \param _f           shared energy threshold for splitting&merging
   * \param _n_pass_max  maximum number of passes (0=full search, the default)
   * \param _split_merge_scale    the scale choice for the split-merge procedure
   *        NOTE: SM_pt leads to IR unsafety for some events with momentum conservation. 
   *              SM_Et is IR safe but not boost invariant and not implemented(!)
   *              SM_mt is IR safe for hadronic events, but not for decays of two 
   *                    back-to-back particles of identical mass
   *              SM_pttilde  
   *                    is always IR safe, and also boost invariant (default)
   * \param _hard_only  when this is set on, only hard jets are computed
   *                    and not the purely ghosted jets (default: false)
   * \return the number of jets (including pure-ghost ones if they are included)
   */
  int compute_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
		    int _n_pass_max=0, Esplit_merge_scale _split_merge_scale=SM_pttilde,
		    bool _hard_only=false);

  /**
   * compute the jet active areas from a given particle set.
   * The parameters of this method are the ones which control the jet clustering alghorithn.
   * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
   * particles/jets and thus is used internally.
   * See compute_areas for paramters definition.
   * \return the number of jets (including pure-ghost ones if they are included)
   */
  int compute_active_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
			   int _n_pass_max=0, Esplit_merge_scale _split_merge_scale=SM_pttilde,
			   bool _hard_only=false);

  /**
   * compute the jet passive areas from a given particle set.
   * The parameters of this method are the ones which control the jet clustering alghorithn.
   * Note that the pt_min is not allowed here soince the jet-area determination involves soft 
   * particles/jets and thus is used internally.
   * See compute_areas for paramters definition.
   */
  int compute_passive_areas(std::vector<Cmomentum> &_particles, double _radius, double _f, 
			    int _n_pass_max=0, Esplit_merge_scale _split_merge_scale=SM_pttilde);

  int grid_size;        ///< size of the grid we add soft particles on (N_soft=(grid_size^2))
  double grid_eta_max;  ///< maximal value of eta we add soft particles on
  double grid_shift;    ///< fractional (random) displacement of the points om the grid

  double pt_soft;       ///< pt of the soft particles added
  double pt_shift;      ///< amplitude of the pt random shift
  double pt_soft_min;   ///< pt_min used in SM to compute passive areas

  /// jets with their areas
  std::vector<Cjet_area> jet_areas;
};

}
#endif
