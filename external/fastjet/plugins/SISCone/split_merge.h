// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: split_merge.h                                                       //
// Description: header file for splitting/merging (contains the CJet class)  //
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
// $Revision:: 367                                                          $//
// $Date:: 2014-09-04 15:57:37 +0200 (Thu, 04 Sep 2014)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __SPLIT_MERGE_H__
#define __SPLIT_MERGE_H__

#include "defines.h"
#include "geom_2d.h"
#include "momentum.h"
#include <stdio.h>
#include <vector>
#include <set>
#include <memory>
#include <string>

namespace siscone{

/**
 * \class Cjet
 * real Jet information.
 *
 * This class contains information for one single jet. 
 * That is, first, its momentum carrying information
 * about its centre and pT, and second, its particle
 * contents
 */
class Cjet{
 public:
  /// default ctor
  Cjet();

  /// default dtor
  ~Cjet();

  Cmomentum v;               ///< jet momentum
  double pt_tilde;           ///< p-scheme pt
  int n;                     ///< number of particles inside
  std::vector<int> contents; ///< particle contents (list of indices)

  /// ordering variable used for ordering and overlap in the
  /// split--merge. This variable is automatically set either to
  /// pt_tilde, or to mt or to pt, depending on the siscone
  /// parameter. Note that the default behaviour is pt_tilde and that
  /// other chices may lead to infrared unsafe situations.
  /// Note: we use the square of the varible rather than the variable itself
  double sm_var2;

  /// covered range in eta-phi
  Ceta_phi_range range;

  /// pass at which the jet has been found
  /// It starts at 0 (first pass), -1 means infinite rapidity
  int pass;
};

/// ordering of jets in pt (e.g. used in final jets ordering)
bool jets_pt_less(const Cjet &j1, const Cjet &j2);
  

/// the choices of scale variable that can be used in the split-merge
/// step, both for ordering the protojets and for measuing their
/// overlap; pt, Et and mt=sqrt(pt^2+m^2) are all defined in E-scheme
/// (4-momentum) recombination; pttilde = \sum_{i\in jet} |p_{t,i}|
///
/// NB: if one changes the order here, one _MUST_ also change the order
///     in the SISCone plugin
enum Esplit_merge_scale {
           SM_pt,     ///< transverse momentum (E-scheme), IR unsafe
           SM_Et,     ///< transverse energy (E-scheme), not long. boost inv.
                      ///< original run-II choice [may not be implemented]
           SM_mt,     ///< transverse mass (E-scheme), IR safe except
                      ///< in decays of two identical narrow heavy particles
           SM_pttilde ///< pt-scheme pt = \sum_{i in jet} |p_{ti}|, should
                      ///< be IR safe in all cases
};

/// return the name of the split-merge scale choice
std::string split_merge_scale_name(Esplit_merge_scale sms);

/**
 * \class Csplit_merge_ptcomparison
 * comparison of jets for split--merge ordering
 *
 * a class that allows us to carry out comparisons of pt of jets, using
 * information from exact particle contents where necessary.
 */
class Csplit_merge_ptcomparison{
public:
  /// default ctor
  Csplit_merge_ptcomparison() : 
    particles(0), split_merge_scale(SM_pttilde){};

  /// return the name corresponding to the SM scale variable
  std::string SM_scale_name() const {
    return split_merge_scale_name(split_merge_scale);}

  std::vector<Cmomentum> * particles; ///< pointer to the list of particles
  std::vector<double> * pt;           ///< pointer to the pt of the particles

  /// comparison between 2 jets
  bool operator()(const Cjet &jet1, const Cjet &jet2) const;

  /**
   * get the difference between 2 jets, calculated such that rounding
   * errors will not affect the result even if the two jets have
   * almost the same content (so that the difference is below the
   * rounding errors)
   *
   * \param j1        first jet
   * \param j2        second jet
   * \param v         jet1-jet2
   * \param pt_tilde  jet1-jet2 pt_tilde
   */
  void get_difference(const Cjet &j1, const Cjet &j2, Cmomentum *v, double *pt_tilde) const;

  /// the following parameter controls the variable we're using for 
  /// the split-merge process i.e. the variable we use for 
  ///  1. ordering jet candidates;
  ///  2. computing the overlap fraction of two candidates.
  /// The default value uses pttile (p-scheme pt). Other alternatives are
  /// pt, mt=sqrt(pt^2+m^2)=sqrt(E^2-pz^2) or Et. 
  /// NOTE: Modifying the default choice can have nasty effects:
  /// - using pt leads to some IR unsafety when we have two jets,
  ///   e.g. back-to-back, with the same pt. In that case, their ordering
  ///   in pt is random and can be affected by the addition of a
  ///   soft particle.  Hence, we highly recommand to keep this to
  ///   the default value i.e.  to use pt only for the purpose of
  ///   investigating the IR issue
  /// - using Et is safe but does not respect boost invariance
  /// - using mt solves the IR unsafety issues with the pt variable
  ///   for QCD jets but the IR unsafety remains for nack-to-back 
  ///   jets of unstable narrow-width particles (e.g. Higgs).
  /// Therefore, keeping the default value is strongly advised.
  Esplit_merge_scale split_merge_scale;
};


// iterator types
/// iterator definition for the jet candidates structure
typedef std::multiset<siscone::Cjet,Csplit_merge_ptcomparison>::iterator cjet_iterator;

/// iterator definition for the jet structure
typedef std::vector<siscone::Cjet>::iterator jet_iterator;



/**
 * \class Csplit_merge
 * Class used to split and merge jets.
 */
class Csplit_merge{
 public:
  /// default ctor
  Csplit_merge();

  /// default dtor
  ~Csplit_merge();


  //////////////////////////////
  // initialisation functions //
  //////////////////////////////

  /**
   * initialisation function
   * \param _particles  list of particles
   * \param protocones  list of protocones (initial jet candidates)
   * \param R2          cone radius (squared)
   * \param ptmin       minimal pT allowed for jets
   * \return 0 on success, 1 on error
   */
  int init(std::vector<Cmomentum> &_particles, std::vector<Cmomentum> *protocones, double R2, double ptmin=0.0);

  /**
   * initialisation function for particle list
   * \param _particles  list of particles
   * \return 0 on success, 1 on error
   */
  int init_particles(std::vector<Cmomentum> &_particles);

  /**
   * build initial list of left particles
   */
  int init_pleft();

  /**
   * use a pt-dependent boundary for splitting
   * When called with true, the criterium for splitting two protojets 
   * will be to compare D1^2/kt1^2 vs. D2^2/kt2^2, the (anti-)kt-weighted 
   * distance instead of the plain distance D1^2 vs. D2^2.
   * This can be set in order to produce more circular hard jets, 
   * with the same underlying philosophy as for the anti-kt algorithm.
   * We thus expect a behaviour closer to the IterativeCone one. 
   * By default, we use the standard D1^2 vs. D2^2 comparison and this 
   * function is not called.
   */
  inline int set_pt_weighted_splitting(bool _use_pt_weighted_splitting){
    use_pt_weighted_splitting = _use_pt_weighted_splitting;
    return 0;
  }

  ////////////////////////
  // cleaning functions //
  ////////////////////////

  /// partial clearance
  int partial_clear();

  /// full clearance
  int full_clear();

  ///////////////////////////////////////
  // user-defined stable-cone ordering //
  ///////////////////////////////////////

  /// \class Cuser_scale_base
  /// base class for user-defined ordering of stable cones
  ///
  /// derived classes have to implement the () operator that returns
  /// the scale associated with a given jet.
  class Cuser_scale_base{
  public:
    /// empty virtual dtor
    virtual ~Cuser_scale_base(){}

    /// the scale associated with a given jet
    ///
    /// "progressive removal" iteratively removes the stable cone with
    /// the largest scale
    virtual double operator()(const Cjet & jet) const = 0;

    /// returns true when the scale associated with jet a is larger than
    /// the scale associated with jet b
    ///
    /// By default this does a simple direct comparison but it can be
    /// overloaded for higher precision [recommended if possible]
    ///
    /// This function assumes that a.sm_var2 and b.sm_var2 have been
    /// correctly initialised with the signed squared output of
    /// operator(), as is by default the case when is_larger is called
    /// from within siscone.
    virtual bool is_larger(const Cjet & a, const Cjet & b) const{
      return (a.sm_var2 > b.sm_var2);
    }
  };

  /// associate a user-defined scale to order the stable cones
  ///
  /// Note that this is only used in "progressive-removal mode",
  /// e.g. in add_hardest_protocone_to_jets().
  void set_user_scale(const Cuser_scale_base * user_scale_in){
    _user_scale = user_scale_in;
  }

  /// return the user-defined scale (NULL if none)
  const Cuser_scale_base * user_scale() const { return _user_scale; }


  /////////////////////////////////
  // main parts of the algorithm //
  /////////////////////////////////
 
  /**
   * build the list 'p_uncol_hard' from p_remain by clustering
   * collinear particles and removing particles softer than
   * stable_cone_soft_pt2_cutoff
   * note that thins in only used for stable-cone detection 
   * so the parent_index field is unnecessary
   */
  int merge_collinear_and_remove_soft();

  /**
   * add a list of protocones
   * \param protocones  list of protocones (initial jet candidates)
   * \param R2          cone radius (squared)
   * \param ptmin       minimal pT allowed for jets
   * \return 0 on success, 1 on error
   */
  int add_protocones(std::vector<Cmomentum> *protocones, double R2, double ptmin=0.0);

  /**
   * remove the hardest protocone and declare it a jet 
   * \param protocones  list of protocones (initial jet candidates)
   * \param R2          cone radius (squared)
   * \param ptmin       minimal pT allowed for jets
   * \return 0 on success, 1 on error
   *
   * The list of remaining particles (and the uncollinear-hard ones)
   * is updated.
   */
  int add_hardest_protocone_to_jets(std::vector<Cmomentum> *protocones, double R2, double ptmin=0.0);

  /**
   * really do the splitting and merging
   * At the end, the vector jets is filled with the jets found.
   * the 'contents' field of each jets contains the indices
   * of the particles included in that jet. 
   * \param overlap_tshold  threshold for splitting/merging transition
   * \param ptmin           minimal pT allowed for jets
   * \return the number of jets is returned
   */
  int perform(double overlap_tshold, double ptmin=0.0);


  //////////////////////////////
  // save and debug functions //
  //////////////////////////////

  /// save final jets
  /// \param flux   stream to save the jet contentss
  int save_contents(FILE *flux);

  /// show jets/candidates status
  int show();

  // particle information
  int n;                               ///< number of particles
  std::vector<Cmomentum> particles;    ///< list of particles
  std::vector<double> pt;              ///< list of particles' pt
  int n_left;                          ///< numer of particles that does not belong to any jet
  std::vector<Cmomentum> p_remain;     ///< list of particles remaining to deal with
  std::vector<Cmomentum> p_uncol_hard; ///< list of particles remaining with collinear clustering
  int n_pass;                          ///< index of the run

  /// minimal difference in squared distance between a particle and
  /// two overlapping protojets when doing a split (useful when
  /// testing approx. collinear safety)
  double most_ambiguous_split; 

  // jets information
  std::vector<Cjet> jets;            ///< list of jets

  // working entries
  int *indices;                      ///< maximal size array for indices works
  int idx_size;                      ///< number of elements in indices1

  /// The following flag indicates that identical protocones
  /// are to be merged automatically each time around the split-merge
  /// loop and before anything else happens.
  ///
  /// This flag is only effective if ALLOW_MERGE_IDENTICAL_PROTOCONES
  /// is set in 'defines.h'
  /// Note that this lead to infrared-unsafety so it is disabled
  /// by default
  bool merge_identical_protocones;

  /// member used for detailed comparisons of pt's
  Csplit_merge_ptcomparison ptcomparison;

  /// stop split--merge or progressive-removal when the squared SM_var
  /// of the hardest protojet is below this cut-off. Note that this is
  /// a signed square (ie SM_var*|SM_var|) to be able to handle
  /// negative values.
  ///
  /// Note that the cut-off is set on the variable squared.
  double SM_var2_hardest_cut_off;

  /// pt cutoff for the particles to put in p_uncol_hard 
  /// this is meant to allow removing soft particles in the
  /// stable-cone search.
  ///
  /// This is not collinear-safe so you should not use this
  /// variable unless you really know what you are doing
  /// Note that the cut-off is set on the variable squared.
  double stable_cone_soft_pt2_cutoff;

 private:
  /**
   * get the overlap between 2 jets
   * \param j1   first jet
   * \param j2   second jet
   * \param v    returned overlap^2 (determined by the choice of SM variable)
   * \return true if overlapping, false if disjoint
   */
  bool get_overlap(const Cjet &j1, const Cjet &j2, double *v);


  /**
   * split the two given jets.
   * during this procedure, the jets j1 & j2 are replaced
   * by 2 new jets. Common particles are associted to the 
   * closest initial jet.
   * \param it_j1  iterator of the first jet in 'candidates'
   * \param it_j2  iterator of the second jet in 'candidates'
   * \param j1     first jet (Cjet instance)
   * \param j2     second jet (Cjet instance)
   * \return true on success, false on error
   */
  bool split(cjet_iterator &it_j1, cjet_iterator &it_j2);

  /**
   * merge the two given jet.
   * during this procedure, the jets j1 & j2 are replaced
   * by 1 single jets containing both of them.
   * \param it_j1  iterator of the first jet in 'candidates'
   * \param it_j2  iterator of the second jet in 'candidates'
   * \return true on success, false on error
   */
  bool merge(cjet_iterator &it_j1, cjet_iterator &it_j2);

  /**
   * Check whether or not a jet has to be inserted in the 
   * list of protojets. If it has, set its sm_variable and
   * insert it to the list of protojets.
   * \param jet    jet to insert
   */
  bool insert(Cjet &jet);

  /**
   * given a 4-momentum and its associated pT, return the 
   * variable tht has to be used for SM
   * \param v          4 momentum of the protojet
   * \param pt_tilde   pt_tilde of the protojet
   */
  double get_sm_var2(Cmomentum &v, double &pt_tilde);

  // jet information
  /// list of jet candidates
  std::auto_ptr<std::multiset<Cjet,Csplit_merge_ptcomparison> > candidates;

  /// minimal pt2
  double pt_min2;

  /**
   * do we have or not to use the pt-weighted splitting
   * (see description for set_pt_weighted_splitting)
   * This will be false by default
   */
  bool use_pt_weighted_splitting;

  /// use a user-defined scale to order the stable cones and jet
  /// candidates
  const Cuser_scale_base *_user_scale;

#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  /// checkxor for the candidates (to avoid having twice the same contents)
  std::set<Creference> cand_refs;
#endif
};

}


#endif
