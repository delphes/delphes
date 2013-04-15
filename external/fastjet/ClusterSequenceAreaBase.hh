//STARTHEADER
// $Id: ClusterSequenceAreaBase.hh 2687 2011-11-14 11:17:51Z soyez $
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

#ifndef __FASTJET_CLUSTERSEQUENCEAREABASE_HH__
#define __FASTJET_CLUSTERSEQUENCEAREABASE_HH__

#include "fastjet/ClusterSequence.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/Selector.hh"

FASTJET_BEGIN_NAMESPACE

/// @ingroup area_classes
/// \class ClusterSequenceAreaBase
/// base class that sets interface for extensions of ClusterSequence
/// that provide information about the area of each jet
///
/// the virtual functions here all return 0, since no area determination
/// is implemented.
class ClusterSequenceAreaBase : public ClusterSequence {
public:
  
  /// a constructor which just carries out the construction of the
  /// parent class
  template<class L> ClusterSequenceAreaBase
         (const std::vector<L> & pseudojets, 
	  const JetDefinition & jet_def_in,
	  const bool & writeout_combinations = false) :
           ClusterSequence(pseudojets, jet_def_in, writeout_combinations) {}


  /// default constructor
  ClusterSequenceAreaBase() {}


  /// destructor
  virtual ~ClusterSequenceAreaBase() {}


  /// return the area associated with the given jet; this base class
  /// returns 0.
  virtual double area       (const PseudoJet & ) const {return 0.0;}

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet; this base class returns 0.
  virtual double area_error (const PseudoJet & ) const {return 0.0;}

  /// return a PseudoJet whose 4-vector is defined by the following integral
  ///
  ///       \int drap d\phi PseudoJet("rap,phi,pt=one") *
  ///                           * Theta("rap,phi inside jet boundary")
  ///
  /// where PseudoJet("rap,phi,pt=one") is a 4-vector with the given
  /// rapidity (rap), azimuth (phi) and pt=1, while Theta("rap,phi
  /// inside jet boundary") is a function that is 1 when rap,phi
  /// define a direction inside the jet boundary and 0 otherwise.
  ///
  /// This base class returns a null 4-vector.
  virtual PseudoJet area_4vector(const PseudoJet & ) const {
    return PseudoJet(0.0,0.0,0.0,0.0);}

  /// true if a jet is made exclusively of ghosts
  ///
  /// NB: most area classes do not give any explicit ghost jets, but
  /// some do, and they should replace this function with their own
  /// version.
  virtual bool is_pure_ghost(const PseudoJet & ) const {
    return false;
  }

  /// returns true if ghosts are explicitly included within 
  /// jets for this ClusterSequence; 
  ///
  /// Derived classes that do include explicit ghosts should provide
  /// an alternative version of this routine and set it properly.
  virtual bool has_explicit_ghosts() const {
    return false;
  }

  /// return the total area, corresponding to the given Selector, that
  /// is free of jets, in general based on the inclusive jets.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  virtual double empty_area(const Selector & selector) const;

  /// return the total area, corresponding to the given Selector, that
  /// is free of jets, based on the supplied all_jets
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  double empty_area_from_jets(const std::vector<PseudoJet> & all_jets,
			      const Selector & selector) const;

  /// return something similar to the number of pure ghost jets
  /// in the given selector's range in an active area case.
  /// For the local implementation we return empty_area/(0.55 pi R^2),
  /// based on measured properties of ghost jets with kt and cam
  /// (cf arXiv:0802.1188).
  ///
  /// Note that the number returned is a double.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  virtual double n_empty_jets(const Selector & selector) const {
    double R = jet_def().R();
    return empty_area(selector)/(0.55*pi*R*R);
  }

  /// the median of (pt/area) for jets contained within the selector
  /// range, making use also of the info on n_empty_jets
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  double median_pt_per_unit_area(const Selector & selector) const;

  /// the median of (pt/area_4vector) for jets contained within the
  /// selector range, making use also of the info on n_empty_jets
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet
  double median_pt_per_unit_area_4vector(const Selector & selector) const;
  
  /// the function that does the work for median_pt_per_unit_area and 
  /// median_pt_per_unit_area_4vector: 
  /// - something_is_area_4vect = false -> use plain area
  /// - something_is_area_4vect = true  -> use 4-vector area
  double median_pt_per_unit_something(
                    const Selector & selector, bool use_area_4vector) const;

  /// using jets withing the selector range (and with 4-vector areas if
  /// use_area_4vector), calculate the median pt/area, as well as an
  /// "error" (uncertainty), which is defined as the 1-sigma
  /// half-width of the distribution of pt/A, obtained by looking for
  /// the point below which we have (1-0.6827)/2 of the jets
  /// (including empty jets).
  ///
  /// The subtraction for a jet with uncorrected pt pt^U and area A is
  ///
  ///   pt^S = pt^U - median*A +- sigma*sqrt(A)
  ///
  /// where the error is only that associated with the fluctuations
  /// in the noise and not that associated with the noise having 
  /// caused changes in the hard-particle content of the jet.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  ///
  /// NB: subtraction may also be done with 4-vector area of course,
  /// and this is recommended for jets with larger values of R, as
  /// long as rho has also been determined with a 4-vector area;
  /// using a scalar area causes one to neglect terms of relative
  /// order $R^2/8$ in the jet $p_t$.
  virtual void get_median_rho_and_sigma(const Selector & selector, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
                                        double & mean_area) const;

  /// a more advanced version of get_median_rho_and_sigma, which allows
  /// one to use any "view" of the event containing all jets (so that, 
  /// e.g. one might use Cam on a different resolution scale without
  /// have to rerun the algorithm).
  ///
  /// By default it will assume that "all" are not inclusive jets, 
  /// so that in dealing with empty area it has to calculate
  /// the number of empty jets based on the empty area and the
  /// the observed <area> of jets rather than a surmised area
  ///
  /// Note that for small effective radii, this can cause problems
  /// because the harder jets get an area >> <ghost-jet-area>
  /// and so the estimate comes out all wrong. In these situations
  /// it is highly advisable to use an area with explicit ghosts, since
  /// then the "empty" jets are actually visible.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  virtual void get_median_rho_and_sigma(const std::vector<PseudoJet> & all_jets,
					const Selector & selector, 
                                        bool use_area_4vector,
                                        double & median, double & sigma,
                                        double & mean_area,
					bool all_are_inclusive = false) const;

  /// same as the full version of get_median_rho_and_error, but without
  /// access to the mean_area
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  virtual void get_median_rho_and_sigma(const Selector & selector, 
                                bool use_area_4vector,
                                double & median, double & sigma) const {
    double mean_area;
    get_median_rho_and_sigma(selector,  use_area_4vector,
                             median,  sigma, mean_area);
  }
  

  /// fits a form pt_per_unit_area(y) = a + b*y^2 in the selector range. 
  /// exclude_above allows one to exclude large values of pt/area from fit. 
  ///               (if negative, the cut is discarded)
  /// use_area_4vector = true uses the 4vector areas.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  virtual void parabolic_pt_per_unit_area(double & a, double & b, 
                                          const Selector & selector, 
                                          double exclude_above=-1.0, 
                                          bool use_area_4vector=false) const;

  /// return a vector of all subtracted jets, using area_4vector, given rho.
  /// Only inclusive_jets above ptmin are subtracted and returned.
  /// the ordering is the same as that of sorted_by_pt(cs.inclusive_jets()),
  /// i.e. not necessarily ordered in pt once subtracted
  std::vector<PseudoJet> subtracted_jets(const double rho,
                                         const double ptmin=0.0) const;

  /// return a vector of subtracted jets, using area_4vector.
  /// Only inclusive_jets above ptmin are subtracted and returned.
  /// the ordering is the same as that of sorted_by_pt(cs.inclusive_jets()),
  /// i.e. not necessarily ordered in pt once subtracted
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  std::vector<PseudoJet> subtracted_jets(const Selector & selector, 
                                         const double ptmin=0.0) const;

  /// return a subtracted jet, using area_4vector, given rho
  PseudoJet subtracted_jet(const PseudoJet & jet,
                           const double rho) const;

  /// return a subtracted jet, using area_4vector; note
  /// that this is potentially inefficient if repeatedly used for many
  /// different jets, because rho will be recalculated each time
  /// around.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  PseudoJet subtracted_jet(const PseudoJet & jet,
                           const Selector & selector) const;

  /// return the subtracted pt, given rho
  double subtracted_pt(const PseudoJet & jet,
                       const double rho,
	               bool use_area_4vector=false) const;

  /// return the subtracted pt; note that this is
  /// potentially inefficient if repeatedly used for many different
  /// jets, because rho will be recalculated each time around.
  /// 
  /// The selector passed as an argument has to have a finite area and
  /// apply jet-by-jet (see the BackgroundEstimator and Subtractor
  /// tools for more generic usages)
  double subtracted_pt(const PseudoJet & jet,
                       const Selector & selector,
	               bool use_area_4vector=false) const;

protected:
  /// check the selector is suited for the computations i.e. applies jet by jet and has a finite area
  void _check_selector_good_for_median(const Selector &selector) const;


private:
  /// handle warning messages
  static LimitedWarning _warnings;
  static LimitedWarning _warnings_zero_area;
  static LimitedWarning _warnings_empty_area;

  /// check the jet algorithm is suitable (and if not issue a warning)
  void _check_jet_alg_good_for_median() const;
  
};



FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEAREABASE_HH__
