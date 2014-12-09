//FJSTARTHEADER
// $Id: ClusterSequenceActiveArea.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include<iostream>
#include<vector>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<valarray>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


using namespace std;


//int ClusterSequenceActiveArea::_n_seed_warnings = 0;
//const int _max_seed_warnings = 10;

//----------------------------------------------------------------------
/// global routine for running active area
void ClusterSequenceActiveArea::_initialise_and_run_AA (
		const JetDefinition & jet_def_in,
		const GhostedAreaSpec & ghost_spec,
		const bool & writeout_combinations) {

  bool continue_running;
  _initialise_AA(jet_def_in,  ghost_spec, writeout_combinations, continue_running);
  if (continue_running) {
    _run_AA(ghost_spec);
    _postprocess_AA(ghost_spec);
  }
}

//----------------------------------------------------------------------
void ClusterSequenceActiveArea::_resize_and_zero_AA () {
  // initialize our local area information
  _average_area.resize(2*_jets.size());  _average_area  = 0.0;
  _average_area2.resize(2*_jets.size()); _average_area2 = 0.0;
  _average_area_4vector.resize(2*_jets.size()); 
  _average_area_4vector = PseudoJet(0.0,0.0,0.0,0.0);
  _non_jet_area = 0.0; _non_jet_area2 = 0.0; _non_jet_number=0.0;
}

//---------------------------------a-------------------------------------
void ClusterSequenceActiveArea::_initialise_AA (
		const JetDefinition & jet_def_in,
		const GhostedAreaSpec & ghost_spec,
		const bool & writeout_combinations,
                bool & continue_running) 
{

  // store this for future use
  _ghost_spec_repeat = ghost_spec.repeat();

  // make sure placeholders are there & zeroed
  _resize_and_zero_AA();
     
  // for future reference...
  _maxrap_for_area = ghost_spec.ghost_maxrap();
  _safe_rap_for_area = _maxrap_for_area - jet_def_in.R();

  // Make sure we'll have at least one repetition -- then we can
  // deduce the unghosted clustering sequence from one of the ghosted
  // sequences. If we do not have any repetitions, then get the
  // unghosted sequence from the plain unghosted clustering.
  //
  // NB: all decanting and filling of initial history will then
  // be carried out by base-class routine
  if (ghost_spec.repeat() <= 0) {
    _initialise_and_run(jet_def_in, writeout_combinations);
    continue_running = false;
    return;
  }

  // transfer all relevant info into internal variables
  _decant_options(jet_def_in, writeout_combinations);

  // set up the history entries for the initial particles (those
  // currently in _jets)
  _fill_initial_history();

  // by default it does not...
  _has_dangerous_particles = false;
  
  continue_running = true;
}


//----------------------------------------------------------------------
void ClusterSequenceActiveArea::_run_AA (const GhostedAreaSpec & ghost_spec) {
  // record the input jets as they are currently
  vector<PseudoJet> input_jets(_jets);

  // code for testing the unique tree
  vector<int> unique_tree;

  // run the clustering multiple times so as to get areas of all the jets
  for (int irepeat = 0; irepeat < ghost_spec.repeat(); irepeat++) {

    ClusterSequenceActiveAreaExplicitGhosts clust_seq(input_jets, 
                                                      jet_def(), ghost_spec);

    _has_dangerous_particles |= clust_seq.has_dangerous_particles();
    if (irepeat == 0) {
      // take the non-ghost part of the history and put into our own
      // history.
      _transfer_ghost_free_history(clust_seq);
      // get the "unique" order that will be used for transferring all areas. 
      unique_tree = unique_history_order();
    }

    // transfer areas from clust_seq into our object
    _transfer_areas(unique_tree, clust_seq);
  }
}
  

//----------------------------------------------------------------------
/// run the postprocessing for the active area (and derived classes)
void ClusterSequenceActiveArea::_postprocess_AA (const GhostedAreaSpec & ghost_spec) {
  _average_area  /= ghost_spec.repeat();
  _average_area2 /= ghost_spec.repeat();
  if (ghost_spec.repeat() > 1) {
    // the VC compiler complains if one puts everything on a single line.
    // An alternative solution would be to use -1.0 (+single line)
    const double tmp = ghost_spec.repeat()-1;
    _average_area2 = sqrt(abs(_average_area2 - _average_area*_average_area)/tmp);
  } else {
    _average_area2 = 0.0;
  }

  _non_jet_area  /= ghost_spec.repeat();
  _non_jet_area2 /= ghost_spec.repeat();
  _non_jet_area2  = sqrt(abs(_non_jet_area2 - _non_jet_area*_non_jet_area)/
			 ghost_spec.repeat());
  _non_jet_number /= ghost_spec.repeat();

  // following bizarre way of writing things is related to 
  // poverty of operations on PseudoJet objects (as well as some confusion
  // in one or two places)
  for (unsigned i = 0; i < _average_area_4vector.size(); i++) {
    _average_area_4vector[i] = (1.0/ghost_spec.repeat()) * _average_area_4vector[i];
  }
  //cerr << "Non-jet area = " << _non_jet_area << " +- " << _non_jet_area2<<endl;
}


// //----------------------------------------------------------------------
// void ClusterSequenceActiveArea::_initialise_and_run_AA (
// 		const JetDefinition & jet_def,
// 		const GhostedAreaSpec & ghost_spec,
// 		const bool & writeout_combinations) 
// {
// 
//   // store this for future use
//   _ghost_spec_repeat = ghost_spec.repeat();
// 
//   // initialize our local area information
//   _average_area.resize(2*_jets.size());  _average_area  = 0.0;
//   _average_area2.resize(2*_jets.size()); _average_area2 = 0.0;
//   _average_area_4vector.resize(2*_jets.size()); 
//   _average_area_4vector = PseudoJet(0.0,0.0,0.0,0.0);
//   _non_jet_area = 0.0; _non_jet_area2 = 0.0; _non_jet_number=0.0;
//      
//   // for future reference...
//   _maxrap_for_area = ghost_spec.ghost_maxrap();
//   _safe_rap_for_area = _maxrap_for_area - jet_def.R();
// 
//   // Make sure we'll have at least one repetition -- then we can
//   // deduce the unghosted clustering sequence from one of the ghosted
//   // sequences. If we do not have any repetitions, then get the
//   // unghosted sequence from the plain unghosted clustering.
//   //
//   // NB: all decanting and filling of initial history will then
//   // be carried out by base-class routine
//   if (ghost_spec.repeat() <= 0) {
//     _initialise_and_run(jet_def, writeout_combinations);
//     return;
//   }
// 
//   // transfer all relevant info into internal variables
//   _decant_options(jet_def, writeout_combinations);
// 
//   // set up the history entries for the initial particles (those
//   // currently in _jets)
//   _fill_initial_history();
// 
//   // record the input jets as they are currently
//   vector<PseudoJet> input_jets(_jets);
// 
//   // code for testing the unique tree
//   vector<int> unique_tree;
// 
//   
//   
// 
//   // run the clustering multiple times so as to get areas of all the jets
//   for (int irepeat = 0; irepeat < ghost_spec.repeat(); irepeat++) {
// 
//     ClusterSequenceActiveAreaExplicitGhosts clust_seq(input_jets, 
//                                                       jet_def, ghost_spec);
// 
//     if (irepeat == 0) {
//       // take the non-ghost part of the history and put into our own
//       // history.
//       _transfer_ghost_free_history(clust_seq);
//       // get the "unique" order that will be used for transferring all areas. 
//       unique_tree = unique_history_order();
//     }
// 
//     // transfer areas from clust_seq into our object
//     _transfer_areas(unique_tree, clust_seq);
//   }
//   
//   _average_area  /= ghost_spec.repeat();
//   _average_area2 /= ghost_spec.repeat();
//   if (ghost_spec.repeat() > 1) {
//     _average_area2 = sqrt(abs(_average_area2 - _average_area*_average_area)/
//                           (ghost_spec.repeat()-1));
//   } else {
//     _average_area2 = 0.0;
//   }
// 
//   _non_jet_area  /= ghost_spec.repeat();
//   _non_jet_area2 /= ghost_spec.repeat();
//   _non_jet_area2  = sqrt(abs(_non_jet_area2 - _non_jet_area*_non_jet_area)/
// 			 ghost_spec.repeat());
//   _non_jet_number /= ghost_spec.repeat();
// 
//   // following bizarre way of writing things is related to 
//   // poverty of operations on PseudoJet objects (as well as some confusion
//   // in one or two places)
//   for (unsigned i = 0; i < _average_area_4vector.size(); i++) {
//     _average_area_4vector[i] = (1.0/ghost_spec.repeat()) * _average_area_4vector[i];
//   }
//   //cerr << "Non-jet area = " << _non_jet_area << " +- " << _non_jet_area2<<endl;
// 
//   
// }
// 


//----------------------------------------------------------------------
double ClusterSequenceActiveArea::pt_per_unit_area(
		       mean_pt_strategies strat, double range) const {
  
  vector<PseudoJet> incl_jets = inclusive_jets();
  vector<double> pt_over_areas;

  for (unsigned i = 0; i < incl_jets.size(); i++) {
    if (abs(incl_jets[i].rap()) < _safe_rap_for_area) {
      double this_area;
      if ( strat == median_4vector ) {
          this_area = area_4vector(incl_jets[i]).perp();
      } else {
          this_area = area(incl_jets[i]);
      }
      pt_over_areas.push_back(incl_jets[i].perp()/this_area);
    }
  }
  
  // there is nothing inside our region, so answer will always be zero
  if (pt_over_areas.size() == 0) {return 0.0;}
  
  // get median (pt/area) [this is the "old" median definition. It considers
  // only the "real" jets in calculating the median, i.e. excluding the
  // only-ghost ones]
  sort(pt_over_areas.begin(), pt_over_areas.end());
  double non_ghost_median_ratio = pt_over_areas[pt_over_areas.size()/2];

  // new median definition that takes into account non-jet area (i.e.
  // jets composed only of ghosts), and for fractional median position 
  // interpolates between the corresponding entries in the pt_over_areas array
  double nj_median_pos = (pt_over_areas.size()-1 - _non_jet_number)/2.0;
  double nj_median_ratio;
  if (nj_median_pos >= 0 && pt_over_areas.size() > 1) {
    int int_nj_median = int(nj_median_pos);
    nj_median_ratio = 
      pt_over_areas[int_nj_median] * (int_nj_median+1-nj_median_pos)
      + pt_over_areas[int_nj_median+1] * (nj_median_pos - int_nj_median);
  } else {
    nj_median_ratio = 0.0;
  }


  // get various forms of mean (pt/area)
  double pt_sum = 0.0, pt_sum_with_cut = 0.0;
  double area_sum = _non_jet_area, area_sum_with_cut = _non_jet_area;
  double ratio_sum = 0.0; 
  double ratio_n = _non_jet_number;
  for (unsigned i = 0; i < incl_jets.size(); i++) {
    if (abs(incl_jets[i].rap()) < _safe_rap_for_area) {
      double this_area;
      if ( strat == median_4vector ) {
          this_area = area_4vector(incl_jets[i]).perp();
      } else {
          this_area = area(incl_jets[i]);
      }
      pt_sum   += incl_jets[i].perp();
      area_sum += this_area;
      double ratio = incl_jets[i].perp()/this_area;
      if (ratio < range*nj_median_ratio) {
	pt_sum_with_cut   += incl_jets[i].perp();
	area_sum_with_cut += this_area;
	ratio_sum += ratio; ratio_n++;
      }
    }
  }
  
  if (strat == play) {
    double trunc_sum = 0, trunc_sumsqr = 0;
    vector<double> means(pt_over_areas.size()), sd(pt_over_areas.size());
    for (unsigned i = 0; i < pt_over_areas.size() ; i++ ) {
      double ratio = pt_over_areas[i];
      trunc_sum += ratio;
      trunc_sumsqr += ratio*ratio;
      means[i] = trunc_sum / (i+1);
      sd[i]    = sqrt(abs(means[i]*means[i]  - trunc_sumsqr/(i+1)));
      cerr << "i, means, sd: " <<i<<", "<< means[i] <<", "<<sd[i]<<", "<<
	sd[i]/sqrt(i+1.0)<<endl;
    }
    cout << "-----------------------------------"<<endl;
    for (unsigned i = 0; i <= pt_over_areas.size()/2 ; i++ ) {
      cout << "Median "<< i <<" = " << pt_over_areas[i]<<endl;
    }
    cout << "Number of non-jets: "<<_non_jet_number<<endl;
    cout << "Area of non-jets: "<<_non_jet_area<<endl;
    cout << "Default median position: " << (pt_over_areas.size()-1)/2.0<<endl;
    cout << "NJ median position: " << nj_median_pos <<endl;
    cout << "NJ median value: " << nj_median_ratio <<endl;
    return 0.0;
  }

  switch(strat) {
  case median:
  case median_4vector:
    return nj_median_ratio;
  case non_ghost_median:
    return non_ghost_median_ratio; 
  case pttot_over_areatot:
    return pt_sum / area_sum;
  case pttot_over_areatot_cut:
    return pt_sum_with_cut / area_sum_with_cut;
  case mean_ratio_cut:
    return ratio_sum/ratio_n;
  default:
    return nj_median_ratio;
  }

}


// The following functionality is now provided by the base class
// //----------------------------------------------------------------------
// // fit a parabola to pt/area as a function of rapidity, using the
// // formulae of CCN28-36 (which actually fits f = a+b*x^2)
// void ClusterSequenceActiveArea::parabolic_pt_per_unit_area(
//        double & a, double & b, double raprange, double exclude_above,
//        bool use_area_4vector) const {
//   
//   double this_raprange;
//   if (raprange <= 0) {this_raprange = _safe_rap_for_area;}
//   else {this_raprange = raprange;}
// 
//   int n=0;
//   int n_excluded = 0;
//   double mean_f=0, mean_x2=0, mean_x4=0, mean_fx2=0; 
// 
//   vector<PseudoJet> incl_jets = inclusive_jets();
// 
//   for (unsigned i = 0; i < incl_jets.size(); i++) {
//     if (abs(incl_jets[i].rap()) < this_raprange) {
//       double this_area;
//       if ( use_area_4vector ) {
//           this_area = area_4vector(incl_jets[i]).perp();     
//       } else {
//           this_area = area(incl_jets[i]);
//       }
//       double f = incl_jets[i].perp()/this_area;
//       if (exclude_above <= 0.0 || f < exclude_above) {
// 	double x = incl_jets[i].rap(); double x2 = x*x;
// 	mean_f   += f;
// 	mean_x2  += x2;
// 	mean_x4  += x2*x2;
// 	mean_fx2 += f*x2;
// 	n++;
//       } else {
// 	n_excluded++;
//       }
//     }
//   }
// 
//   if (n <= 1) {
//     // meaningful results require at least two jets inside the
//     // area -- mind you if there are empty jets we should be in 
//     // any case doing something special...
//     a = 0.0;
//     b = 0.0;
//   } else {
//     mean_f   /= n;
//     mean_x2  /= n;
//     mean_x4  /= n;
//     mean_fx2 /= n;
//     
//     b = (mean_f*mean_x2 - mean_fx2)/(mean_x2*mean_x2 - mean_x4);
//     a = mean_f - b*mean_x2;
//   }
//   //cerr << "n_excluded = "<< n_excluded << endl;
// }


//----------------------------------------------------------------------
double ClusterSequenceActiveArea::empty_area(const Selector & selector) const {
  // make sure that the selector applies jet by jet
  if (! selector.applies_jet_by_jet()){
    throw Error("ClusterSequenceActiveArea: empty area can only be computed from selectors applying jet by jet");
  }

  double empty = 0.0;
  // first deal with ghost jets
  for (unsigned  i = 0; i < _ghost_jets.size(); i++) {
    if (selector.pass(_ghost_jets[i])) {
      empty += _ghost_jets[i].area;
    }
  }
  // then deal with unclustered ghosts
  for (unsigned  i = 0; i < _unclustered_ghosts.size(); i++) {
    if (selector.pass(_unclustered_ghosts[i])) {
      empty += _unclustered_ghosts[i].area;
    }
  }
  empty /= _ghost_spec_repeat;
  return empty;
}

//----------------------------------------------------------------------
double ClusterSequenceActiveArea::n_empty_jets(const Selector & selector) const {
  _check_selector_good_for_median(selector);

  double inrange = 0;
  for (unsigned  i = 0; i < _ghost_jets.size(); i++) {
    if (selector.pass(_ghost_jets[i])) inrange++;
  }
  inrange /= _ghost_spec_repeat;
  return inrange;
}

//----------------------------------------------------------------------
/// transfer the history (and jet-momenta) from clust_seq to our
/// own internal structure while removing ghosts
void ClusterSequenceActiveArea::_transfer_ghost_free_history(
             const ClusterSequenceActiveAreaExplicitGhosts & ghosted_seq) {
  
  const vector<history_element> & gs_history  = ghosted_seq.history();
  vector<int> gs2self_hist_map(gs_history.size());

  // first transfer info about strategy used (which isn't necessarily
  // always the one that got asked for...)
  _strategy = ghosted_seq.strategy_used();

  // work our way through to first non-trivial combination
  unsigned igs = 0;
  unsigned iself = 0;
  while (igs < gs_history.size() && gs_history[igs].parent1 == InexistentParent) {
    // record correspondence 
    if (!ghosted_seq.is_pure_ghost(igs)) {
      gs2self_hist_map[igs] = iself++; 
    } else {
      gs2self_hist_map[igs] = Invalid; 
    }
    igs++;
  };

  // make sure the count of non-ghost initial jets is equal to
  // what we already have in terms of initial jets
  assert(iself == _history.size());

  // if there was no clustering in this event (e.g. SISCone passive
  // area with zero input particles, or with a pt cut on stable cones
  // that kills all jets), then don't bother with the rest (which
  // would crash!)
  if (igs == gs_history.size()) return;
  
  // now actually transfer things
  do  {
    // if we are a pure ghost, then go on to next round
    if (ghosted_seq.is_pure_ghost(igs)) {
      gs2self_hist_map[igs] = Invalid;
      continue;
    }

    const history_element & gs_hist_el = gs_history[igs];

    bool parent1_is_ghost = ghosted_seq.is_pure_ghost(gs_hist_el.parent1);
    bool parent2_is_ghost = ghosted_seq.is_pure_ghost(gs_hist_el.parent2);

    // if exactly one parent is a ghost then maintain info about the
    // non-ghost correspondence for this jet, and then go on to next
    // recombination in the ghosted sequence
    if (parent1_is_ghost && !parent2_is_ghost && gs_hist_el.parent2 >= 0) {
      gs2self_hist_map[igs] = gs2self_hist_map[gs_hist_el.parent2];
      continue;
    }
    if (!parent1_is_ghost && parent2_is_ghost) {
      gs2self_hist_map[igs] = gs2self_hist_map[gs_hist_el.parent1];
      continue;
    }

    // no parents are ghosts...
    if (gs_hist_el.parent2 >= 0) {
      // recombination of two non-ghosts
      gs2self_hist_map[igs] = _history.size();
      // record the recombination in our own sequence
      int newjet_k; // dummy var -- not used
      int jet_i = _history[gs2self_hist_map[gs_hist_el.parent1]].jetp_index;
      int jet_j = _history[gs2self_hist_map[gs_hist_el.parent2]].jetp_index;
      //cerr << "recombining "<< jet_i << " and "<< jet_j << endl;
      _do_ij_recombination_step(jet_i, jet_j, gs_hist_el.dij, newjet_k);
    } else {
      // we have a non-ghost that has become a beam-jet
      assert(gs_history[igs].parent2 == BeamJet);
      // record position
      gs2self_hist_map[igs] = _history.size();
      // record the recombination in our own sequence
      _do_iB_recombination_step(
             _history[gs2self_hist_map[gs_hist_el.parent1]].jetp_index,
             gs_hist_el.dij);
    }
  } while (++igs < gs_history.size());

}

//----------------------------------------------------------------------
void ClusterSequenceActiveArea::_transfer_areas(
	    const vector<int> & unique_hist_order,
    	    const ClusterSequenceActiveAreaExplicitGhosts & ghosted_seq  ) {

  const vector<history_element> & gs_history  = ghosted_seq.history();
  const vector<PseudoJet>       & gs_jets     = ghosted_seq.jets();
  vector<int>    gs_unique_hist_order = ghosted_seq.unique_history_order();

  const double tolerance = 1e-11; // to decide when two jets are the same

  int j = -1;
  int hist_index = -1;
  
  valarray<double> our_areas(_history.size());
  our_areas = 0.0;

  valarray<PseudoJet> our_area_4vectors(_history.size());
  our_area_4vectors = PseudoJet(0.0,0.0,0.0,0.0);

  for (unsigned i = 0; i < gs_history.size(); i++) {
    // only consider composite particles
    unsigned gs_hist_index = gs_unique_hist_order[i];
    if (gs_hist_index < ghosted_seq.n_particles()) continue;
    const history_element & gs_hist = gs_history[gs_unique_hist_order[i]];
    int parent1 = gs_hist.parent1;
    int parent2 = gs_hist.parent2;

    if (parent2 == BeamJet) {
      // need to look at parent to get the actual jet
      const PseudoJet & jet = 
  	  gs_jets[gs_history[parent1].jetp_index];
      double area_local = ghosted_seq.area(jet);
      PseudoJet ext_area = ghosted_seq.area_4vector(jet);

      if (ghosted_seq.is_pure_ghost(parent1)) {
        // record the existence of the pure ghost jet for future use
        _ghost_jets.push_back(GhostJet(jet,area_local));
	if (abs(jet.rap()) < _safe_rap_for_area) {
	  _non_jet_area  += area_local;
	  _non_jet_area2 += area_local*area_local;
	  _non_jet_number += 1;
	}
      } else {

	// get next "combined-particle" index in our own history
	// making sure we don't go beyond its bounds (if we do
	// then we're in big trouble anyway...)
	while (++j < int(_history.size())) {
	  hist_index = unique_hist_order[j];
	  if (hist_index >= _initial_n) break;}

        // sanity checking -- do not overrun
        if (j >= int(_history.size())) throw Error("ClusterSequenceActiveArea: overran reference array in diB matching");

        // sanity check -- make sure we are taking about the same 
        // jet in reference and new sequences
        int refjet_index = _history[_history[hist_index].parent1].jetp_index;
        assert(refjet_index >= 0 && refjet_index < int(_jets.size()));
	const PseudoJet & refjet = _jets[refjet_index];

      //cerr << "Inclusive" << endl;
      //cerr << gs_history[parent1].jetp_index << " " << gs_jets.size() << endl;
      //cerr << _history[_history[hist_index].parent1].jetp_index << " " << _jets.size() << endl;

        // If pt disagrees check E; if they both disagree there's a
        // problem here... NB: a massive particle with zero pt may
        // have its pt changed when a ghost is added -- this is why we
        // also require the energy to be wrong before complaining
        _throw_unless_jets_have_same_perp_or_E(jet, refjet, tolerance,
                                               ghosted_seq);

	// set the area at this clustering stage
	our_areas[hist_index]  = area_local; 
	our_area_4vectors[hist_index]  = ext_area; 

	// update the parent as well -- that way its area is the area
	// immediately before clustering (i.e. resolve an ambiguity in
	// the Cambridge case and ensure in the kt case that the original
	// particles get a correct area)
	our_areas[_history[hist_index].parent1] = area_local;
	our_area_4vectors[_history[hist_index].parent1] = ext_area;
	
      }
    }
    else if (!ghosted_seq.is_pure_ghost(parent1) && 
	     !ghosted_seq.is_pure_ghost(parent2)) {

      // get next "combined-particle" index in our own history
      while (++j < int(_history.size())) {
	hist_index = unique_hist_order[j];
	if (hist_index >= _initial_n) break;}
      
      // sanity checking -- do not overrun
      if (j >= int(_history.size())) throw Error("ClusterSequenceActiveArea: overran reference array in dij matching");

      // make sure that our reference history entry is also for
      // an exclusive (dij) clustering (otherwise the comparison jet
      // will not exist)
      if (_history[hist_index].parent2 == BeamJet) throw Error("ClusterSequenceActiveArea: could not match clustering sequences (encountered dij matched with diB)");

      //cerr << "Exclusive: hist_index,hist_size: " << hist_index << " " << _history.size()<< endl;
      //cerr << gs_hist.jetp_index << " " << gs_jets.size() << endl;
      //cerr << _history[hist_index].jetp_index << " " << _jets.size() << endl;

      const PseudoJet & jet = gs_jets[gs_hist.jetp_index];
      const PseudoJet & refjet = _jets[_history[hist_index].jetp_index];

      // run sanity check 
      _throw_unless_jets_have_same_perp_or_E(jet, refjet, tolerance,
                                             ghosted_seq);

      // update area and our local index (maybe redundant since later
      // the descendants will reupdate it?)
      double area_local  = ghosted_seq.area(jet);
      our_areas[hist_index]  += area_local; 

      PseudoJet ext_area = ghosted_seq.area_4vector(jet);

      // GPS TMP debugging (jetclu) -----------------------
      //ext_area = PseudoJet(1e-100,1e-100,1e-100,4e-100);
      //our_area_4vectors[hist_index] = ext_area;
      //cout << "aa " 
      //     << our_area_4vectors[hist_index].px() << " "
      //     << our_area_4vectors[hist_index].py() << " "
      //     << our_area_4vectors[hist_index].pz() << " "
      //     << our_area_4vectors[hist_index].E() << endl;
      //cout << "bb " 
      //     << ext_area.px() << " "
      //     << ext_area.py() << " "
      //     << ext_area.pz() << " "
      //     << ext_area.E() << endl;
      //---------------------------------------------------

      _jet_def.recombiner()->plus_equal(our_area_4vectors[hist_index], ext_area);

      // now update areas of parents (so that they becomes areas
      // immediately before clustering occurred). This is of use
      // because it allows us to set the areas of the original hard
      // particles in the kt algorithm; for the Cambridge case it
      // means a jet's area will be the area just before it clusters
      // with another hard jet.
      const PseudoJet & jet1 = gs_jets[gs_history[parent1].jetp_index];
      int our_parent1 = _history[hist_index].parent1;
      our_areas[our_parent1] = ghosted_seq.area(jet1);
      our_area_4vectors[our_parent1] = ghosted_seq.area_4vector(jet1);

      const PseudoJet & jet2 = gs_jets[gs_history[parent2].jetp_index];
      int our_parent2 = _history[hist_index].parent2;
      our_areas[our_parent2] = ghosted_seq.area(jet2);
      our_area_4vectors[our_parent2] = ghosted_seq.area_4vector(jet2);
    }

  }

  // now add unclustered ghosts to the relevant list so that we can
  // calculate empty area later.
  vector<PseudoJet> unclust = ghosted_seq.unclustered_particles();
  for (unsigned iu = 0; iu < unclust.size();  iu++) {
    if (ghosted_seq.is_pure_ghost(unclust[iu])) {
      double area_local = ghosted_seq.area(unclust[iu]);
      _unclustered_ghosts.push_back(GhostJet(unclust[iu],area_local));
    }
  }

  /*
   * WARNING:
   *   _average_area has explicitly been sized initially to 2*jets().size()
   *   which can be bigger than our_areas (of size _history.size()
   *   if there are some unclustered particles. 
   *   So we must take care about boundaries
   */

  for (unsigned int area_index = 0; area_index<our_areas.size(); area_index++){
    _average_area[area_index]  += our_areas[area_index]; 
    _average_area2[area_index] += our_areas[area_index]*our_areas[area_index]; 
  }

  //_average_area_4vector += our_area_4vectors;
  // Use the proper recombination scheme when averaging the area_4vectors
  // over multiple ghost runs (i.e. the repeat stage);
  for (unsigned i = 0; i < our_area_4vectors.size(); i++) {
    _jet_def.recombiner()->plus_equal(_average_area_4vector[i],
                                       our_area_4vectors[i]);
  }
}


/// check if two jets have the same momentum to within the
/// tolerance (and if pt's are not the same we're forgiving and
/// look to see if the energy is the same)
void ClusterSequenceActiveArea::_throw_unless_jets_have_same_perp_or_E(
                                const PseudoJet & jet, 
                                const PseudoJet & refjet, 
                                double tolerance,
          const ClusterSequenceActiveAreaExplicitGhosts & jets_ghosted_seq
) const {

  if (abs(jet.perp2()-refjet.perp2()) > 
      tolerance*max(jet.perp2(),refjet.perp2())
      && abs(jet.E()-refjet.E()) > tolerance*max(jet.E(),refjet.E())) {
    ostringstream ostr;
    ostr << "Could not match clustering sequence for an inclusive/exclusive jet when reconstructing areas. See FAQ for possible explanations." << endl;
    ostr << "  Ref-Jet: "
         << refjet.px() << " " 
         << refjet.py() << " " 
         << refjet.pz() << " " 
         << refjet.E() << endl;
    ostr << "  New-Jet: "
         << jet.px() << " " 
         << jet.py() << " " 
         << jet.pz() << " " 
         << jet.E() << endl;
    if (jets_ghosted_seq.has_dangerous_particles()) {
      ostr << "  NB: some particles have pt too low wrt ghosts -- this may be the cause" << endl;}
    //ostr << jet.perp() << " " << refjet.perp() << " "
    //     << jet.perp() - refjet.perp() << endl;
    throw Error(ostr.str());
  }
}

FASTJET_END_NAMESPACE

