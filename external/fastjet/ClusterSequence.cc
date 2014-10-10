//FJSTARTHEADER
// $Id: ClusterSequence.cc 3685 2014-09-11 20:15:00Z salam $
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

#include "fastjet/Error.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceStructure.hh"
#include "fastjet/version.hh" // stores the current version number
#include "fastjet/internal/LazyTiling9Alt.hh"
#include "fastjet/internal/LazyTiling9.hh"
#include "fastjet/internal/LazyTiling25.hh"
#ifndef __FJCORE__
#include "fastjet/internal/LazyTiling9SeparateGhosts.hh"
#endif  // __FJCORE__
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include<string>
#include<set>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
// here's where we put the main page for fastjet (as explained in the
// Doxygen FAQ)
// We put it inside the fastjet namespace to have the links without
// having to specify (fastjet::)
//......................................................................
/** \mainpage FastJet code documentation
 *
 * These pages provide automatically generated documentation for the 
 * FastJet package.
 * 
 * \section useful_classes The most useful classes
 *
 * Many of the facilities of FastJet can be accessed through the three
 * following classes:
 *
 * - PseudoJet: the basic class for holding the 4-momentum of a
 *   particle or a jet.
 *
 * - JetDefinition: the combination of a #JetAlgorithm and its
 *   associated parameters. Can also be initialised with a \ref plugins "plugin".  
 *
 * - ClusterSequence: constructed with a vector of input (PseudoJet)
 *   particles and a JetDefinition, it computes and stores the
 *   information on how the input particles are clustered into jets.
 *
 * \section advanced_classes Selected more advanced classes
 *
 * - ClusterSequenceArea: with the help of an AreaDefinition, provides
 *   jets that also contain information about their area.
 *
 * \section Tools Selected additional tools
 *
 * - JetMedianBackgroundEstimator: with the help of a Selector, a JetDefinition and
 *   an AreaDefinition, allows one to estimate the background noise density in an event; for a simpler, quicker, effective alternative, use GridMedianBackgroundEstimator
 *
 * - Transformer: class from which are derived various tools for
 *   manipulating jets and accessing their substructure. Examples are
 *   Subtractor, Filter, Pruner and various taggers (e.g. JHTopTagger
 *   and MassDropTagger).
 *
 * \section further_info Further information
 *
 * - Selected classes ordered by topics can be found under the <a
 * href="modules.html">modules</a> tab.
 *
 * - The complete list of classes is available under the  <a
 * href="annotated.html">classes</a> tab.
 * 
 * - For non-class material (<a href="namespacefastjet.html#enum-members">enums</a>,
 * <a href="namespacefastjet.html#typedef-members">typedefs</a>, 
 * <a href="namespacefastjet.html#func-members">functions</a>), see the 
 * #fastjet documentation
 * 
 * - For further information and normal documentation, see the main <a
 * href="http://fastjet.fr/">FastJet</a> page.
 *
 * \section examples Examples
 *   See our \subpage Examples page
 */

// define the doxygen groups
/// \defgroup basic_classes    Fundamental FastJet classes
/// \defgroup area_classes     Area-related classes
/// \defgroup sec_area_classes Secondary area-related classes
/// \defgroup plugins          Plugins for non-native jet definitions
/// \defgroup selectors        Selectors
/// \defgroup tools            FastJet tools
/// \{ \defgroup tools_generic     Generic tools
///    \defgroup tools_background  Background subtraction
///    \defgroup tools_taggers     Taggers
/// \}
/// \defgroup extra_info       Access to extra information
/// \defgroup error_handling   Error handling
/// \defgroup advanced_usage   Advanced usage
/// \if internal_doc
/// \defgroup internal
/// \endif

//----------------------------------------------------------------------


using namespace std;


// The following variable can be modified from within user code
// so as to redirect banners to an ostream other than cout.
//
// Please note that if you distribute 3rd party code
// that links with FastJet, that 3rd party code is NOT
// allowed to turn off the printing of FastJet banners
// by default. This requirement reflects the spirit of
// clause 2c of the GNU Public License (v2), under which
// FastJet and its plugins are distributed.
std::ostream * ClusterSequence::_fastjet_banner_ostr = &cout;


// destructor that guarantees proper bookkeeping for the CS Structure
ClusterSequence::~ClusterSequence () {
  // set the pointer in the wrapper to this object to NULL to say that
  // we're going out of scope
  if (_structure_shared_ptr()){
    ClusterSequenceStructure* csi = dynamic_cast<ClusterSequenceStructure*>(_structure_shared_ptr()); 
    // normally the csi is purely internal so it really should not be
    // NULL i.e assert should be OK
    // (we assert rather than throw an error, since failure here is a
    // sign of major internal problems)
    assert(csi != NULL);
    csi->set_associated_cs(NULL);

    // if the user had given the CS responsibility to delete itself,
    // but then deletes the CS themselves, the following lines of
    // code will ensure that the structure_shared_ptr will have
    // a proper object count (so that jets associated with the CS will
    // throw the correct error if the user tries to access their
    // constituents).
    if (_deletes_self_when_unused) {
      _structure_shared_ptr.set_count(_structure_shared_ptr.use_count() 
				        + _structure_use_count_after_construction);
    }
  }
}

//-----------
void ClusterSequence::signal_imminent_self_deletion() const {
  // normally if the destructor is called when
  // _deletes_self_when_unused is true, it assumes that it's been
  // called by the user (and it therefore resets the shared pointer
  // count to the true count).
  //
  // for self deletion (called from the destructor of the CSstructure,
  // the shared_ptr to which has just had its pointer -> 0) you do
  // _not_ want to reset the pointer count (otherwise you will end up
  // with a double delete on the shared pointer once you start
  // deleting the internal structure of the CS).
  //
  // the following modification ensures that the count reset will not
  // take place in the destructor
  assert(_deletes_self_when_unused);
  _deletes_self_when_unused = false;
}

//DEP //----------------------------------------------------------------------
//DEP void ClusterSequence::_initialise_and_run (
//DEP 				  const double R,
//DEP 				  const Strategy & strategy,
//DEP 				  const bool & writeout_combinations) {
//DEP 
//DEP   JetDefinition jet_def(_default_jet_algorithm, R, strategy);
//DEP   _initialise_and_run(jet_def, writeout_combinations);
//DEP }


//----------------------------------------------------------------------
void ClusterSequence::_initialise_and_run (
				  const JetDefinition & jet_def_in,
				  const bool & writeout_combinations) {

  // transfer all relevant info into internal variables
  _decant_options(jet_def_in, writeout_combinations);

  // now run
  _initialise_and_run_no_decant();
}

//----------------------------------------------------------------------
void ClusterSequence::_initialise_and_run_no_decant () {

  // set up the history entries for the initial particles (those
  // currently in _jets)
  _fill_initial_history();

  // don't run anything if the event is empty
  if (n_particles() == 0) return;

  // ----- deal with special cases: plugins & e+e- ------
  if (_jet_algorithm == plugin_algorithm) {
    // allows plugin_xyz() functions to modify cluster sequence
    _plugin_activated = true;
    // let the plugin do its work here
    _jet_def.plugin()->run_clustering( (*this) );
    _plugin_activated = false;
    _update_structure_use_count();
    return;
  } else if (_jet_algorithm == ee_kt_algorithm ||
	     _jet_algorithm == ee_genkt_algorithm) {
    // ignore requested strategy
    _strategy = N2Plain;
    if (_jet_algorithm == ee_kt_algorithm) {
      // make sure that R is large enough so that "beam" recomb only
      // occurs when a single particle is left
      // Normally, this should be automatically set to 4 from JetDefinition
      assert(_Rparam > 2.0); 
      // this is used to renormalise the dij to get a "standard" form
      // and our convention in e+e- will be different from that
      // in long.inv case; NB: _invR2 name should be changed -> _renorm_dij?
      _invR2 = 1.0;
    } else {
      // as of 2009-01-09, choose R to be an angular distance, in
      // radians.  Since the algorithm uses 2(1-cos(theta)) as its
      // squared angular measure, make sure that the _R2 is defined
      // in a similar way.
      if (_Rparam > pi) {
	// choose a value that ensures that back-to-back particles will
	// always recombine 
	//_R2 = 4.0000000000001;
	_R2 = 2 * ( 3.0 + cos(_Rparam) );
      } else {
	_R2    = 2 * ( 1.0 - cos(_Rparam) );
      }
      _invR2 = 1.0/_R2;
    }
    _simple_N2_cluster_EEBriefJet();
    return;
  } else if (_jet_algorithm == undefined_jet_algorithm) {
    throw Error("A ClusterSequence cannot be created with an uninitialised JetDefinition");
  }


  // automatically redefine the strategy according to N if that is
  // what the user requested -- transition points (and especially
  // their R-dependence) are based on empirical observations for a
  // R=0.4, 0.7 and 1.0, running on toth (3.4GHz, Pentium IV D [dual
  // core] with 2MB of cache).
  //-------------
  // 2011-11-15: lowered N2Plain -> N2Tiled switchover based on some
  //             new tests on an Intel Core 2 Duo T9400 @ 2.53 GHz
  //             with 6MB cache; tests performed with lines such as
  //             ./fastjet_timing_plugins -kt -nhardest 30 -repeat 50000 -strategy -3 -R 0.5 -nev 1  <  ../../data/Pythia-PtMin1000-LHC-1000ev.dat
  if (_strategy == Best) {
    _strategy = _best_strategy();
#ifdef DROP_CGAL
    // fall back strategy for large N when CGAL is missing
    if (_strategy == NlnN) _strategy = N2MHTLazy25;
#endif  // DROP_CGAL
  } else if (_strategy == BestFJ30) {
    int N = _jets.size();
    //if (N <= 55*max(0.5,min(1.0,_Rparam))) {// old empirical scaling with R
    //----------------------
    // 2011-11-15: new empirical scaling with R; NB: low-R N2Tiled
    // could be significantly improved at low N by limiting the
    // minimum size of tiles when R is small
    if (min(1.0,max(0.1,_Rparam)*3.3)*N <= 30) {
      _strategy = N2Plain;
    } else if (N > 6200/pow(_Rparam,2.0) && _jet_def.jet_algorithm() == cambridge_algorithm) {
      _strategy = NlnNCam;
#ifndef DROP_CGAL
    } else if ((N > 16000/pow(_Rparam,1.15) && _jet_def.jet_algorithm() != antikt_algorithm)
	       || N > 35000/pow(_Rparam,1.15)) {
      _strategy = NlnN;
#endif  // DROP_CGAL
    } else if (N <= 450) {
      _strategy = N2Tiled;
    } else {                   
      _strategy = N2MinHeapTiled;
    }
  }

  // R >= 2pi is not supported by all clustering strategies owing to
  // periodicity issues (a particle might cluster with itself). When
  // R>=2pi, we therefore automatically switch to a strategy that is
  // known to work.
  if (_Rparam >= twopi) {
    if (   _strategy == NlnN
	|| _strategy == NlnN3pi
	|| _strategy == NlnNCam
	|| _strategy == NlnNCam2pi2R
	|| _strategy == NlnNCam4pi) {
#ifdef DROP_CGAL
      _strategy = N2MinHeapTiled;
#else
      _strategy = NlnN4pi;
#endif    
    }
    if (_jet_def.strategy() != Best && _strategy != _jet_def.strategy()) {
      ostringstream oss;
      oss << "Cluster strategy " << strategy_string(_jet_def.strategy())
	  << " automatically changed to " << strategy_string()
	  << " because the former is not supported for R = " << _Rparam
	  << " >= 2pi";
      _changed_strategy_warning.warn(oss.str());
    }
  }


  // run the code containing the selected strategy
  // 
  // We order the strategies starting from the ones used by the Best
  // strategy in the order of increasing N, then the remaining ones
  // again in the order of increasing N.
  if (_strategy == N2Plain) {
    // BriefJet provides standard long.invariant kt alg.
    this->_simple_N2_cluster_BriefJet();
  } else if (_strategy == N2Tiled) {
    this->_faster_tiled_N2_cluster();
  } else if (_strategy == N2MinHeapTiled) {
    this->_minheap_faster_tiled_N2_cluster();
  } else if (_strategy == N2MHTLazy9Alt) {
    // attempt to use an external tiling routine -- it manipulates
    // the CS history via the plugin mechanism
    _plugin_activated = true;
    LazyTiling9Alt tiling(*this);
    tiling.run();
    _plugin_activated = false;

  } else if (_strategy == N2MHTLazy25) {
    // attempt to use an external tiling routine -- it manipulates
    // the CS history via the plugin mechanism
    _plugin_activated = true;
    LazyTiling25 tiling(*this);
    tiling.run();
    _plugin_activated = false;

  } else if (_strategy == N2MHTLazy9) {
    // attempt to use an external tiling routine -- it manipulates
    // the CS history via the plugin mechanism
    _plugin_activated = true;
    LazyTiling9 tiling(*this);
    tiling.run();
    _plugin_activated = false;

#ifndef __FJCORE__
  } else if (_strategy == N2MHTLazy9AntiKtSeparateGhosts) {
    // attempt to use an external tiling routine -- it manipulates
    // the CS history via the plugin mechanism
    _plugin_activated = true;
    LazyTiling9SeparateGhosts tiling(*this);
    tiling.run();
    _plugin_activated = false;
#else 
    throw Error("N2MHTLazy9AntiKtSeparateGhosts strategy not supported with FJCORE");
#endif  // __FJCORE__

  } else if (_strategy == NlnN) {
    this->_delaunay_cluster();
  } else if (_strategy == NlnNCam) {
    this->_CP2DChan_cluster_2piMultD();
  } else if (_strategy == NlnN3pi || _strategy == NlnN4pi ) {
    this->_delaunay_cluster();
  } else if (_strategy ==  N3Dumb ) {
    this->_really_dumb_cluster();
  } else if (_strategy == N2PoorTiled) {
    this->_tiled_N2_cluster();
  } else if (_strategy == NlnNCam4pi) {
    this->_CP2DChan_cluster();
  } else if (_strategy == NlnNCam2pi2R) {
    this->_CP2DChan_cluster_2pi2R();
  } else {
    ostringstream err;
    err << "Unrecognised value for strategy: "<<_strategy;
    throw Error(err.str());
  }

}


// these needs to be defined outside the class definition.
bool ClusterSequence::_first_time = true;
LimitedWarning ClusterSequence::_exclusive_warnings;


//----------------------------------------------------------------------
// the version string
string fastjet_version_string() {
  return "FastJet version "+string(fastjet_version);
}


//----------------------------------------------------------------------
// prints a banner on the first call
void ClusterSequence::print_banner() {

  if (!_first_time) {return;}
  _first_time = false;

  // make sure the user has not set the banner stream to NULL
  ostream * ostr = _fastjet_banner_ostr;
  if (!ostr) return;  

  (*ostr) << "#--------------------------------------------------------------------------\n";
  (*ostr) << "#                         FastJet release " << fastjet_version << endl;
  (*ostr) << "#                 M. Cacciari, G.P. Salam and G. Soyez                  \n"; 
  (*ostr) << "#     A software package for jet finding and analysis at colliders      \n";
  (*ostr) << "#                           http://fastjet.fr                           \n"; 
  (*ostr) << "#	                                                                      \n";
  (*ostr) << "# Please cite EPJC72(2012)1896 [arXiv:1111.6097] if you use this package\n";
  (*ostr) << "# for scientific work and optionally PLB641(2006)57 [hep-ph/0512210].   \n";
  (*ostr) << "#                                                                       \n";
  (*ostr) << "# FastJet is provided without warranty under the terms of the GNU GPLv2.\n";
  (*ostr) << "# It uses T. Chan's closest pair algorithm, S. Fortune's Voronoi code";
#ifndef DROP_CGAL
  (*ostr) << ",\n# CGAL ";
#else
  (*ostr) << "\n# ";
#endif  // DROP_CGAL
  (*ostr) << "and 3rd party plugin jet algorithms. See COPYING file for details.\n";
  (*ostr) << "#--------------------------------------------------------------------------\n";
  // make sure we really have the output done.
  ostr->flush();
}

//----------------------------------------------------------------------
// transfer all relevant info into internal variables
void ClusterSequence::_decant_options(const JetDefinition & jet_def_in,
                                      const bool & writeout_combinations) {
  // make a local copy of the jet definition (for future use)
  _jet_def = jet_def_in;
  _writeout_combinations = writeout_combinations;
  // initialised the wrapper to the current CS
  _structure_shared_ptr.reset(new ClusterSequenceStructure(this));

  _decant_options_partial();
}

//----------------------------------------------------------------------
// transfer all relevant info into internal variables
void ClusterSequence::_decant_options_partial() {
  // let the user know what's going on
  print_banner();
  
  _jet_algorithm = _jet_def.jet_algorithm();
  _Rparam = _jet_def.R();  _R2 = _Rparam*_Rparam; _invR2 = 1.0/_R2;
  _strategy = _jet_def.strategy();

  // disallow interference from the plugin
  _plugin_activated = false;

  // initialised the wrapper to the current CS
  //_structure_shared_ptr.reset(new ClusterSequenceStructure(this));
  _update_structure_use_count(); // make sure it's correct already here
}


//----------------------------------------------------------------------
// initialise the history in a standard way
void ClusterSequence::_fill_initial_history () {

  //if (_jets.size() == 0) {throw Error("Cannot run jet-finder on empty event");}

  // reserve sufficient space for everything
  _jets.reserve(_jets.size()*2);
  _history.reserve(_jets.size()*2);

  _Qtot = 0;

  for (int i = 0; i < static_cast<int>(_jets.size()) ; i++) {
    history_element element;
    element.parent1 = InexistentParent;
    element.parent2 = InexistentParent;
    element.child   = Invalid;
    element.jetp_index = i;
    element.dij     = 0.0;
    element.max_dij_so_far = 0.0;

    _history.push_back(element);
    
    // do any momentum preprocessing needed by the recombination scheme
    _jet_def.recombiner()->preprocess(_jets[i]);

    // get cross-referencing right from PseudoJets
    _jets[i].set_cluster_hist_index(i);
    _set_structure_shared_ptr(_jets[i]);

    // determine the total energy in the event
    _Qtot += _jets[i].E();
  }
  _initial_n = _jets.size();
  _deletes_self_when_unused = false;
}


//----------------------------------------------------------------------
string ClusterSequence::strategy_string (Strategy strategy_in)  const {
  string strategy;
  switch(strategy_in) {
  case NlnN:
    strategy = "NlnN"; break;
  case NlnN3pi:
    strategy = "NlnN3pi"; break;
  case NlnN4pi:
    strategy = "NlnN4pi"; break;
  case N2Plain:
    strategy = "N2Plain"; break;
  case N2Tiled:
    strategy = "N2Tiled"; break;
  case N2MinHeapTiled:
    strategy = "N2MinHeapTiled"; break;
  case N2PoorTiled:
    strategy = "N2PoorTiled"; break;
  case N2MHTLazy9:
    strategy = "N2MHTLazy9"; break;
  case N2MHTLazy9Alt:
    strategy = "N2MHTLazy9Alt"; break;
  case N2MHTLazy25:
    strategy = "N2MHTLazy25"; break;
  case N2MHTLazy9AntiKtSeparateGhosts:
    strategy = "N2MHTLazy9AntiKtSeparateGhosts"; break;
  case N3Dumb:
    strategy = "N3Dumb"; break;
  case NlnNCam4pi:
    strategy = "NlnNCam4pi"; break;
  case NlnNCam2pi2R:
    strategy = "NlnNCam2pi2R"; break;
  case NlnNCam:
    strategy = "NlnNCam"; break; // 2piMultD
  case plugin_strategy:
    strategy = "plugin strategy"; break;
  default:
    strategy = "Unrecognized";
  }
  return strategy;
}  


double ClusterSequence::jet_scale_for_algorithm(
				  const PseudoJet & jet) const {
  if (_jet_algorithm == kt_algorithm)             {return jet.kt2();}
  else if (_jet_algorithm == cambridge_algorithm) {return 1.0;}
  else if (_jet_algorithm == antikt_algorithm) {
    double kt2=jet.kt2();
    return kt2 > 1e-300 ? 1.0/kt2 : 1e300;
  } else if (_jet_algorithm == genkt_algorithm) {
    double kt2 = jet.kt2();
    double p   = jet_def().extra_param();
    if (p <= 0 && kt2 < 1e-300) kt2 = 1e-300; // dodgy safety check
    return pow(kt2, p);
  } else if (_jet_algorithm == cambridge_for_passive_algorithm) {
    double kt2 = jet.kt2();
    double lim = _jet_def.extra_param();
    if (kt2 < lim*lim && kt2 != 0.0) {
      return 1.0/kt2;
    } else {return 1.0;}
  } else {throw Error("Unrecognised jet algorithm");}
}

//----------------------------------------------------------------------
// returns a suggestion for the best strategy to use on event
// multiplicity, algorithm, R, etc.
//
// Some of the work to establish the best strategy is collected in
// issue-tracker/2014-07-auto-strategy-selection;
// transition_fit_v2.fit indicates the results of the fits that we're
// using here. (Automatically generated by transition_fit_v2.gp).
//
// The transition to NlnN is always present, and it is the the
// caller's responsibility to drop back down to N2MHTLazy25 if NlnN
// isn't available.
//
// This routine should be called only if the jet alg is one of kt,
// antikt, cam or genkt.
Strategy ClusterSequence::_best_strategy() const {
  int N = _jets.size();
  // define bounded R, always above 0.1, because we don't trust any
  // of our parametrizations below R = 0.1
  double bounded_R = max(_Rparam, 0.1);

  // the very first test thing is a quick hard-coded test to decide
  // if we immediately opt for N2Plain
  if (N <= 30 || N <= 39.0/(bounded_R + 0.6)) {
    return N2Plain;
  } 
  
  // Define objects that describe our various boundaries. A prefix N_
  // indicates that boundary is for N, while L_ means it's for log(N).
  //
  // Hopefully having them static will ensure minimal overhead
  // in creating them; collecting them in one place should
  // help with updates?
  //
  const static _Parabola N_Tiled_to_MHT_lowR             (-45.4947,54.3528,44.6283);
  const static _Parabola L_MHT_to_MHTLazy9_lowR          (0.677807,-1.05006,10.6994);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_akt_lowR(0.169967,-0.512589,12.1572);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_kt_lowR (0.16237,-0.484612,12.3373);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_cam_lowR = L_MHTLazy9_to_MHTLazy25_kt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_akt_lowR    (0.0472051,-0.22043,15.9196);
  const static _Parabola L_MHTLazy25_to_NlnN_kt_lowR     (0.118609,-0.326811,14.8287);
  const static _Parabola L_MHTLazy25_to_NlnN_cam_lowR    (0.10119,-0.295748,14.3924);

  const static _Line     L_Tiled_to_MHTLazy9_medR         (-1.31304,7.29621);
  const static _Parabola L_MHTLazy9_to_MHTLazy25_akt_medR = L_MHTLazy9_to_MHTLazy25_akt_lowR;
  const static _Parabola L_MHTLazy9_to_MHTLazy25_kt_medR  = L_MHTLazy9_to_MHTLazy25_kt_lowR;
  const static _Parabola L_MHTLazy9_to_MHTLazy25_cam_medR = L_MHTLazy9_to_MHTLazy25_cam_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_akt_medR     = L_MHTLazy25_to_NlnN_akt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_kt_medR      = L_MHTLazy25_to_NlnN_kt_lowR;
  const static _Parabola L_MHTLazy25_to_NlnN_cam_medR     = L_MHTLazy25_to_NlnN_cam_lowR;

  const static double    N_Plain_to_MHTLazy9_largeR         = 75;
  const static double    N_MHTLazy9_to_MHTLazy25_akt_largeR = 700;
  const static double    N_MHTLazy9_to_MHTLazy25_kt_largeR  = 1000;
  const static double    N_MHTLazy9_to_MHTLazy25_cam_largeR = 1000;
  const static double    N_MHTLazy25_to_NlnN_akt_largeR     = 100000;
  const static double    N_MHTLazy25_to_NlnN_kt_largeR      = 40000;
  const static double    N_MHTLazy25_to_NlnN_cam_largeR     = 15000;

  // We have timing studies only for kt, cam and antikt; for other
  // algorithms we set the local jet_algorithm variable to the one of
  // kt,cam,antikt that we think will be closest in behaviour to the
  // other alg.
  JetAlgorithm jet_algorithm;
  if (_jet_algorithm == genkt_algorithm) {
    // for genkt, then we set the local jet_algorithm variable (used
    // only for strategy choice) to be either kt or antikt, depending on
    // the p value.
    double p   = jet_def().extra_param();
    if (p < 0.0) jet_algorithm = antikt_algorithm;
    else         jet_algorithm =     kt_algorithm;
  } else if (_jet_algorithm == cambridge_for_passive_algorithm) {
    // we assume (but haven't tested) that using the kt-alg timing
    // transitions should be adequate for cambridge_for_passive_algorithm
    jet_algorithm = kt_algorithm;
  } else {
    jet_algorithm = _jet_algorithm;
  }

  if (bounded_R < 0.65) {
    // low R case
    if          (N    < N_Tiled_to_MHT_lowR(bounded_R))              return N2Tiled;
    double logN = log(double(N));
    if          (logN < L_MHT_to_MHTLazy9_lowR(bounded_R))           return N2MinHeapTiled;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_akt_lowR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_akt_lowR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_kt_lowR(bounded_R))  return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_kt_lowR(bounded_R))      return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (logN < L_MHTLazy9_to_MHTLazy25_cam_lowR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_cam_lowR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  } else if (bounded_R < 0.5*pi) {
    // medium R case
    double logN = log(double(N));
    if      (logN < L_Tiled_to_MHTLazy9_medR(bounded_R))             return N2Tiled;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_akt_medR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_akt_medR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (logN < L_MHTLazy9_to_MHTLazy25_kt_medR(bounded_R))  return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_kt_medR(bounded_R))      return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (logN < L_MHTLazy9_to_MHTLazy25_cam_medR(bounded_R)) return N2MHTLazy9;
        else if (logN < L_MHTLazy25_to_NlnN_cam_medR(bounded_R))     return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  } else {
    // large R case (R > pi/2)
    if      (N    < N_Plain_to_MHTLazy9_largeR)                      return N2Plain;
    else {
      if (jet_algorithm == antikt_algorithm){
        if      (N < N_MHTLazy9_to_MHTLazy25_akt_largeR)             return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_akt_largeR)                 return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == kt_algorithm){
        if      (N < N_MHTLazy9_to_MHTLazy25_kt_largeR)              return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_kt_largeR)                  return N2MHTLazy25;
        else                                                         return NlnN;
      } else if (jet_algorithm == cambridge_algorithm)  {
        if      (N < N_MHTLazy9_to_MHTLazy25_cam_largeR)             return N2MHTLazy9;
        else if (N < N_MHTLazy25_to_NlnN_cam_largeR)                 return N2MHTLazy25;
        else                                                         return NlnNCam;
      }
    }
  }
  
  bool code_should_never_reach_here = false;
  assert(code_should_never_reach_here); 
  return N2MHTLazy9;

}


// //----------------------------------------------------------------------
// /// transfer the sequence contained in other_seq into our own;
// /// any plugin "extras" contained in the from_seq will be lost
// /// from there.
// void ClusterSequence::transfer_from_sequence(ClusterSequence & from_seq) {
// 
//   if (will_delete_self_when_unused()) 
//     throw(Error("cannot use CS::transfer_from_sequence after a call to delete_self_when_unused()"));
// 
//   // the metadata
//   _jet_def                 = from_seq._jet_def                ;
//   _writeout_combinations   = from_seq._writeout_combinations  ;
//   _initial_n               = from_seq._initial_n              ;
//   _Rparam                  = from_seq._Rparam                 ;
//   _R2                      = from_seq._R2                     ;
//   _invR2                   = from_seq._invR2                  ;
//   _strategy                = from_seq._strategy               ;
//   _jet_algorithm           = from_seq._jet_algorithm          ;
//   _plugin_activated        = from_seq._plugin_activated       ;
// 
//   // the data
//   _jets     = from_seq._jets;
//   _history  = from_seq._history;
//   // the following transfers ownership of the extras from the from_seq
//   _extras   = from_seq._extras;
// 
//   // transfer of ownership
//   if (_structure_shared_ptr()) {
//     // anything that is currently associated with the cluster sequence
//     // should be told that its cluster sequence no longer exists
//     ClusterSequenceStructure* csi = dynamic_cast<ClusterSequenceStructure*>(_structure_shared_ptr()); 
//     assert(csi != NULL);
//     csi->set_associated_cs(NULL);
//   }
//   // create a new _structure_shared_ptr to reflect the fact that
//   // this CS is essentially a new one
//   _structure_shared_ptr.reset(new ClusterSequenceStructure(this));
//   _update_structure_use_count();
//   
//   for (vector<PseudoJet>::iterator jit = _jets.begin(); jit != _jets.end(); jit++)
//     _set_structure_shared_ptr(*jit);
// }


//----------------------------------------------------------------------
// transfer the sequence contained in other_seq into our own;
// any plugin "extras" contained in the from_seq will be lost
// from there.
//
// It also sets the ClusterSequence pointers of the PseudoJets in
// the history to point to this ClusterSequence
//
// The second argument is an action that will be applied on every
// jets in the resulting ClusterSequence
void ClusterSequence::transfer_from_sequence(const ClusterSequence & from_seq,
					     const FunctionOfPseudoJet<PseudoJet> * action_on_jets){

  if (will_delete_self_when_unused()) 
    throw(Error("cannot use CS::transfer_from_sequence after a call to delete_self_when_unused()"));

  // the metadata
  _jet_def                 = from_seq._jet_def                ;
  _writeout_combinations   = from_seq._writeout_combinations  ;
  _initial_n               = from_seq._initial_n              ;
  _Rparam                  = from_seq._Rparam                 ;
  _R2                      = from_seq._R2                     ;
  _invR2                   = from_seq._invR2                  ;
  _strategy                = from_seq._strategy               ;
  _jet_algorithm           = from_seq._jet_algorithm          ;
  _plugin_activated        = from_seq._plugin_activated       ;

  // the data

  // apply the transformation on the jets if needed
  if (action_on_jets)
    _jets     = (*action_on_jets)(from_seq._jets);
  else
    _jets     = from_seq._jets;
  _history  = from_seq._history;
  // the following shares ownership of the extras with the from_seq;
  // no transformations will be applied to the extras
  _extras   = from_seq._extras;

  // clean up existing structure
  if (_structure_shared_ptr()) {
    // If there are jets associated with an old version of the CS and
    // a new one, keeping track of when to delete the CS becomes more
    // complex; so we don't allow this situation to occur.
    if (_deletes_self_when_unused) throw Error("transfer_from_sequence cannot be used for a cluster sequence that deletes self when unused");
    
    // anything that is currently associated with the cluster sequence
    // should be told that its cluster sequence no longer exists
    ClusterSequenceStructure* csi = dynamic_cast<ClusterSequenceStructure*>(_structure_shared_ptr()); 
    assert(csi != NULL);
    csi->set_associated_cs(NULL);
  }
  // create a new _structure_shared_ptr to reflect the fact that
  // this CS is essentially a new one
  _structure_shared_ptr.reset(new ClusterSequenceStructure(this));
  _update_structure_use_count();
  
  for (unsigned int i=0; i<_jets.size(); i++){
    // we reset the cluster history index in case action_on_jets
    // messed up with it
    _jets[i].set_cluster_hist_index(from_seq._jets[i].cluster_hist_index());

    // reset the structure pointer
    _set_structure_shared_ptr(_jets[i]);
  }
}


//----------------------------------------------------------------------
// record an ij recombination and reset the _jets[newjet_k] momentum and
// user index to be those of newjet
void ClusterSequence::plugin_record_ij_recombination(
	   int jet_i, int jet_j, double dij, 
	   const PseudoJet & newjet, int & newjet_k) {

  plugin_record_ij_recombination(jet_i, jet_j, dij, newjet_k);

  // now transfer newjet into place
  int tmp_index = _jets[newjet_k].cluster_hist_index();
  _jets[newjet_k] = newjet;
  _jets[newjet_k].set_cluster_hist_index(tmp_index);
  _set_structure_shared_ptr(_jets[newjet_k]);
}


//----------------------------------------------------------------------
// return all inclusive jets with pt > ptmin
vector<PseudoJet> ClusterSequence::inclusive_jets (const double ptmin) const{
  double dcut = ptmin*ptmin;
  int i = _history.size() - 1; // last jet
  vector<PseudoJet> jets_local;
  if (_jet_algorithm == kt_algorithm) {
    while (i >= 0) {
      // with our specific definition of dij and diB (i.e. R appears only in 
      // dij), then dij==diB is the same as the jet.perp2() and we can exploit
      // this in selecting the jets...
      if (_history[i].max_dij_so_far < dcut) {break;}
      if (_history[i].parent2 == BeamJet && _history[i].dij >= dcut) {
	// for beam jets
	int parent1 = _history[i].parent1;
	jets_local.push_back(_jets[_history[parent1].jetp_index]);}
      i--;
    }
  } else if (_jet_algorithm == cambridge_algorithm) {
    while (i >= 0) {
      // inclusive jets are all at end of clustering sequence in the
      // Cambridge algorithm -- so if we find a non-exclusive jet, then
      // we can exit
      if (_history[i].parent2 != BeamJet) {break;}
      int parent1 = _history[i].parent1;
      const PseudoJet & jet = _jets[_history[parent1].jetp_index];
      if (jet.perp2() >= dcut) {jets_local.push_back(jet);}
      i--;
    }
  } else if (_jet_algorithm == plugin_algorithm 
             || _jet_algorithm == ee_kt_algorithm
             || _jet_algorithm == antikt_algorithm
             || _jet_algorithm == genkt_algorithm
             || _jet_algorithm == ee_genkt_algorithm
             || _jet_algorithm == cambridge_for_passive_algorithm) {
    // for inclusive jets with a plugin algorithm, we make no
    // assumptions about anything (relation of dij to momenta,
    // ordering of the dij, etc.)
    while (i >= 0) {
      if (_history[i].parent2 == BeamJet) {
	int parent1 = _history[i].parent1;
	const PseudoJet & jet = _jets[_history[parent1].jetp_index];
	if (jet.perp2() >= dcut) {jets_local.push_back(jet);}
      }
      i--;
    }
  } else {throw Error("cs::inclusive_jets(...): Unrecognized jet algorithm");}
  return jets_local;
}


//----------------------------------------------------------------------
// return the number of exclusive jets that would have been obtained
// running the algorithm in exclusive mode with the given dcut
int ClusterSequence::n_exclusive_jets (const double dcut) const {

  // first locate the point where clustering would have stopped (i.e. the
  // first time max_dij_so_far > dcut)
  int i = _history.size() - 1; // last jet
  while (i >= 0) {
    if (_history[i].max_dij_so_far <= dcut) {break;}
    i--;
  }
  int stop_point = i + 1;
  // relation between stop_point, njets assumes one extra jet disappears
  // at each clustering.
  int njets = 2*_initial_n - stop_point;
  return njets;
}

//----------------------------------------------------------------------
// return all exclusive jets that would have been obtained running
// the algorithm in exclusive mode with the given dcut
vector<PseudoJet> ClusterSequence::exclusive_jets (const double dcut) const {
  int njets = n_exclusive_jets(dcut);
  return exclusive_jets(njets);
}


//----------------------------------------------------------------------
// return the jets obtained by clustering the event to n jets.
// Throw an error if there are fewer than n particles.
vector<PseudoJet> ClusterSequence::exclusive_jets (const int njets) const {

  // make sure the user does not ask for more than jets than there
  // were particles in the first place.
  if (njets > _initial_n) {
    ostringstream err;
    err << "Requested " << njets << " exclusive jets, but there were only " 
	<< _initial_n << " particles in the event";
    throw Error(err.str());
  }

  return exclusive_jets_up_to(njets);
}

//----------------------------------------------------------------------
// return the jets obtained by clustering the event to n jets.
// If there are fewer than n particles, simply return all particles
vector<PseudoJet> ClusterSequence::exclusive_jets_up_to (const int njets) const {

  // provide a warning when extracting exclusive jets for algorithms 
  // that does not support it explicitly.
  // Native algorithm that support it are: kt, ee_kt, Cambridge/Aachen, 
  //   genkt and ee_genkt (both with p>=0)
  // For plugins, we check Plugin::exclusive_sequence_meaningful()
  if (( _jet_def.jet_algorithm() != kt_algorithm) &&
      ( _jet_def.jet_algorithm() != cambridge_algorithm) &&
      ( _jet_def.jet_algorithm() != ee_kt_algorithm) &&
      (((_jet_def.jet_algorithm() != genkt_algorithm) && 
	(_jet_def.jet_algorithm() != ee_genkt_algorithm)) || 
       (_jet_def.extra_param() <0)) &&
      ((_jet_def.jet_algorithm() != plugin_algorithm) ||
       (!_jet_def.plugin()->exclusive_sequence_meaningful()))) {
    _exclusive_warnings.warn("dcut and exclusive jets for jet-finders other than kt, C/A or genkt with p>=0 should be interpreted with care.");
  }


  // calculate the point where we have to stop the clustering.
  // relation between stop_point, njets assumes one extra jet disappears
  // at each clustering.
  int stop_point = 2*_initial_n - njets;
  // make sure it's safe when more jets are requested than there are particles
  if (stop_point < _initial_n) stop_point = _initial_n;

  // some sanity checking to make sure that e+e- does not give us
  // surprises (should we ever implement e+e-)...
  if (2*_initial_n != static_cast<int>(_history.size())) {
    ostringstream err;
    err << "2*_initial_n != _history.size() -- this endangers internal assumptions!\n";
    throw Error(err.str());
    //assert(false);
  }

  // now go forwards and reconstitute the jets that we have --
  // basically for any history element, see if the parent jets to
  // which it refers were created before the stopping point -- if they
  // were then add them to the list, otherwise they are subsequent
  // recombinations of the jets that we are looking for.
  vector<PseudoJet> jets_local;
  for (unsigned int i = stop_point; i < _history.size(); i++) {
    int parent1 = _history[i].parent1;
    if (parent1 < stop_point) {
      jets_local.push_back(_jets[_history[parent1].jetp_index]);
    }
    int parent2 = _history[i].parent2;
    if (parent2 < stop_point && parent2 > 0) {
      jets_local.push_back(_jets[_history[parent2].jetp_index]);
    }
    
  }

  // sanity check...
  if (int(jets_local.size()) != min(_initial_n, njets)) {
    ostringstream err;
    err << "ClusterSequence::exclusive_jets: size of returned vector ("
	 <<jets_local.size()<<") does not coincide with requested number of jets ("
	 <<njets<<")";
    throw Error(err.str());
  }

  return jets_local;
}

//----------------------------------------------------------------------
/// return the dmin corresponding to the recombination that went from
/// n+1 to n jets
double ClusterSequence::exclusive_dmerge (const int njets) const {
  assert(njets >= 0);
  if (njets >= _initial_n) {return 0.0;}
  return _history[2*_initial_n-njets-1].dij;
}


//----------------------------------------------------------------------
/// return the maximum of the dmin encountered during all recombinations 
/// up to the one that led to an n-jet final state; identical to
/// exclusive_dmerge, except in cases where the dmin do not increase
/// monotonically.
double ClusterSequence::exclusive_dmerge_max (const int njets) const {
  assert(njets >= 0);
  if (njets >= _initial_n) {return 0.0;}
  return _history[2*_initial_n-njets-1].max_dij_so_far;
}


//----------------------------------------------------------------------
/// return a vector of all subjets of the current jet (in the sense
/// of the exclusive algorithm) that would be obtained when running
/// the algorithm with the given dcut.
std::vector<PseudoJet> ClusterSequence::exclusive_subjets 
   (const PseudoJet & jet, const double dcut) const {

  set<const history_element*> subhist;

  // get the set of history elements that correspond to subjets at
  // scale dcut
  get_subhist_set(subhist, jet, dcut, 0);

  // now transfer this into a sequence of jets
  vector<PseudoJet> subjets;
  subjets.reserve(subhist.size());
  for (set<const history_element*>::iterator elem = subhist.begin(); 
       elem != subhist.end(); elem++) {
    subjets.push_back(_jets[(*elem)->jetp_index]);
  }
  return subjets;
}

//----------------------------------------------------------------------
/// return the size of exclusive_subjets(...); still n ln n with same
/// coefficient, but marginally more efficient than manually taking
/// exclusive_subjets.size()
int ClusterSequence::n_exclusive_subjets(const PseudoJet & jet, 
                        const double dcut) const {
  set<const history_element*> subhist;
  // get the set of history elements that correspond to subjets at
  // scale dcut
  get_subhist_set(subhist, jet, dcut, 0);
  return subhist.size();
}

//----------------------------------------------------------------------
/// return the list of subjets obtained by unclustering the supplied
/// jet down to nsub subjets. Throws an error if there are fewer than
/// nsub particles in the jet.
std::vector<PseudoJet> ClusterSequence::exclusive_subjets
   (const PseudoJet & jet, int nsub) const {
  vector<PseudoJet> subjets = exclusive_subjets_up_to(jet, nsub);
  if (int(subjets.size()) < nsub) {
    ostringstream err;
    err << "Requested " << nsub << " exclusive subjets, but there were only " 
	<< subjets.size() << " particles in the jet";
    throw Error(err.str());
  }
  return subjets;

}

//----------------------------------------------------------------------
/// return the list of subjets obtained by unclustering the supplied
/// jet down to nsub subjets (or all constituents if there are fewer
/// than nsub).
std::vector<PseudoJet> ClusterSequence::exclusive_subjets_up_to
   (const PseudoJet & jet, int nsub) const {

  set<const history_element*> subhist;

  // prepare the vector into which we'll put the result
  vector<PseudoJet> subjets;
  if (nsub <  0) throw Error("Requested a negative number of subjets. This is nonsensical.");
  if (nsub == 0) return subjets;

  // get the set of history elements that correspond to subjets at
  // scale dcut
  get_subhist_set(subhist, jet, -1.0, nsub);

  // now transfer this into a sequence of jets
  subjets.reserve(subhist.size());
  for (set<const history_element*>::iterator elem = subhist.begin(); 
       elem != subhist.end(); elem++) {
    subjets.push_back(_jets[(*elem)->jetp_index]);
  }
  return subjets;
}


//----------------------------------------------------------------------
/// return the dij that was present in the merging nsub+1 -> nsub 
/// subjets inside this jet.
/// 
/// If the jet has nsub or fewer constituents, it will return 0.
double ClusterSequence::exclusive_subdmerge(const PseudoJet & jet, int nsub) const {
  set<const history_element*> subhist;

  // get the set of history elements that correspond to subjets at
  // scale dcut
  get_subhist_set(subhist, jet, -1.0, nsub);
  
  set<const history_element*>::iterator highest = subhist.end();
  highest--;
  /// will be zero if nconst <= nsub, since highest will be an original 
  /// particle have zero dij
  return (*highest)->dij;
}


//----------------------------------------------------------------------
/// return the maximum dij that occurred in the whole event at the
/// stage that the nsub+1 -> nsub merge of subjets occurred inside 
/// this jet.
///
/// If the jet has nsub or fewer constituents, it will return 0.
double ClusterSequence::exclusive_subdmerge_max(const PseudoJet & jet, int nsub) const {

  set<const history_element*> subhist;

  // get the set of history elements that correspond to subjets at
  // scale dcut
  get_subhist_set(subhist, jet, -1.0, nsub);
  
  set<const history_element*>::iterator highest = subhist.end();
  highest--;
  /// will be zero if nconst <= nsub, since highest will be an original 
  /// particle have zero dij
  return (*highest)->max_dij_so_far;
}



//----------------------------------------------------------------------
/// return a set of pointers to history entries corresponding to the
/// subjets of this jet; one stops going working down through the
/// subjets either when 
///   - there is no further to go
///   - one has found maxjet entries
///   - max_dij_so_far <= dcut
void ClusterSequence::get_subhist_set(set<const history_element*> & subhist,
                                     const  PseudoJet & jet, 
                                     double dcut, int maxjet) const {
  assert(contains(jet));
  
  subhist.clear();
  subhist.insert(&(_history[jet.cluster_hist_index()]));

  // establish the set of jets that are relevant
  int njet = 1;
  while (true) {
    // first find out if we need to probe deeper into jet.
    // Get history element closest to end of sequence
    set<const history_element*>::iterator highest = subhist.end();
    assert (highest != subhist.begin()); 
    highest--;
    const history_element* elem = *highest;
    // make sure we haven't got too many jets
    if (njet == maxjet) break;
    // make sure it has parents
    if (elem->parent1 < 0)            break;
    // make sure that we still resolve it at scale dcut
    if (elem->max_dij_so_far <= dcut) break;

    // then do so: replace "highest" with its two parents
    subhist.erase(highest);
    subhist.insert(&(_history[elem->parent1]));
    subhist.insert(&(_history[elem->parent2]));
    njet++;
  }
}

//----------------------------------------------------------------------
// work through the object's history until
bool ClusterSequence::object_in_jet(const PseudoJet & object, 
                                    const PseudoJet & jet) const {

  // make sure the object conceivably belongs to this clustering
  // sequence
  assert(contains(object) && contains(jet));

  const PseudoJet * this_object = &object;
  const PseudoJet * childp;
  while(true) {
    if (this_object->cluster_hist_index() == jet.cluster_hist_index()) {
      return true;
    } else if (has_child(*this_object, childp)) {
      this_object = childp;
    } else {
      return false;
    }
  }
}

//----------------------------------------------------------------------
/// if the jet has parents in the clustering, it returns true
/// and sets parent1 and parent2 equal to them.
///
/// if it has no parents it returns false and sets parent1 and
/// parent2 to zero
bool ClusterSequence::has_parents(const PseudoJet & jet, PseudoJet & parent1, 
                              PseudoJet & parent2) const {

  const history_element & hist = _history[jet.cluster_hist_index()];

  // make sure we do not run into any unexpected situations --
  // i.e. both parents valid, or neither
  assert ((hist.parent1 >= 0 && hist.parent2 >= 0) || 
          (hist.parent1 < 0 && hist.parent2 < 0));

  if (hist.parent1 < 0) {
    parent1 = PseudoJet(0.0,0.0,0.0,0.0);
    parent2 = parent1;
    return false;
  } else {
    parent1 = _jets[_history[hist.parent1].jetp_index];
    parent2 = _jets[_history[hist.parent2].jetp_index];
    // order the parents in decreasing pt
    if (parent1.perp2() < parent2.perp2()) std::swap(parent1,parent2);
    return true;
  }
}

//----------------------------------------------------------------------
/// if the jet has a child then return true and give the child jet
/// otherwise return false and set the child to zero
bool ClusterSequence::has_child(const PseudoJet & jet, PseudoJet & child) const {

  //const history_element & hist = _history[jet.cluster_hist_index()];
  //
  //if (hist.child >= 0) {
  //  child = _jets[_history[hist.child].jetp_index];
  //  return true;
  //} else {
  //  child = PseudoJet(0.0,0.0,0.0,0.0);
  //  return false;
  //}
  const PseudoJet * childp;
  bool res = has_child(jet, childp);
  if (res) {
    child = *childp;
    return true;
  } else {
    child = PseudoJet(0.0,0.0,0.0,0.0);
    return false;
  }
}

bool ClusterSequence::has_child(const PseudoJet & jet, const PseudoJet * & childp) const {

  const history_element & hist = _history[jet.cluster_hist_index()];

  // check that this jet has a child and that the child corresponds to
  // a true jet [RETHINK-IF-CHANGE-NUMBERING: what is the right
  // behaviour if the child is the same jet but made inclusive...?]
  if (hist.child >= 0 && _history[hist.child].jetp_index >= 0) {
    childp = &(_jets[_history[hist.child].jetp_index]);
    return true;
  } else {
    childp = NULL;
    return false;
  }
}


//----------------------------------------------------------------------
/// if this jet has a child (and so a partner) return true
/// and give the partner, otherwise return false and set the
/// partner to zero
bool ClusterSequence::has_partner(const PseudoJet & jet, 
                              PseudoJet & partner) const {

  const history_element & hist = _history[jet.cluster_hist_index()];

  // make sure we have a child and that the child does not correspond
  // to a clustering with the beam (or some other invalid quantity)
  if (hist.child >= 0 && _history[hist.child].parent2 >= 0) {
    const history_element & child_hist = _history[hist.child];
    if (child_hist.parent1 == jet.cluster_hist_index()) {
      // partner will be child's parent2 -- for iB clustering
      // parent2 will not be valid
      partner = _jets[_history[child_hist.parent2].jetp_index];
    } else {
      // partner will be child's parent1
      partner = _jets[_history[child_hist.parent1].jetp_index];
    }
    return true;
  } else {
    partner = PseudoJet(0.0,0.0,0.0,0.0);
    return false;
  }
}


//----------------------------------------------------------------------
// return a vector of the particles that make up a jet
vector<PseudoJet> ClusterSequence::constituents (const PseudoJet & jet) const {
  vector<PseudoJet> subjets;
  add_constituents(jet, subjets);
  return subjets;
}

//----------------------------------------------------------------------
/// output the supplied vector of jets in a format that can be read
/// by an appropriate root script; the format is:
/// jet-n jet-px jet-py jet-pz jet-E 
///   particle-n particle-rap particle-phi particle-pt
///   particle-n particle-rap particle-phi particle-pt
///   ...
/// #END
/// ... [i.e. above repeated]
void ClusterSequence::print_jets_for_root(const std::vector<PseudoJet> & jets_in, 
                                          ostream & ostr) const {
  for (unsigned i = 0; i < jets_in.size(); i++) {
    ostr << i  << " "
         << jets_in[i].px() << " "
         << jets_in[i].py() << " "
         << jets_in[i].pz() << " "
         << jets_in[i].E() << endl;
    vector<PseudoJet> cst = constituents(jets_in[i]);
    for (unsigned j = 0; j < cst.size() ; j++) {
      ostr << " " << j << " "
           << cst[j].rap() << " "
           << cst[j].phi() << " "
           << cst[j].perp() << endl;
    }
    ostr << "#END" << endl;
  }
}

void ClusterSequence::print_jets_for_root(const std::vector<PseudoJet> & jets_in, 
					  const std::string & filename,
					  const std::string & comment ) const {
  std::ofstream ostr(filename.c_str());
  if (comment != "") ostr << "# " << comment << endl;
  print_jets_for_root(jets_in, ostr);
}


// Not yet. Perhaps in a future release
// //----------------------------------------------------------------------
// // print out all inclusive jets with pt > ptmin
// void ClusterSequence::print_jets (const double ptmin) const{
//     vector<PseudoJet> jets = sorted_by_pt(inclusive_jets(ptmin));
// 
//     for (size_t j = 0; j < jets.size(); j++) {
//        printf("%5u %7.3f %7.3f %9.3f\n",
//        j,jets[j].rap(),jets[j].phi(),jets[j].perp());
//     }
// }

//----------------------------------------------------------------------
/// returns a vector of size n_particles() which indicates, for 
/// each of the initial particles (in the order in which they were
/// supplied), which of the supplied jets it belongs to; if it does
/// not belong to any of the supplied jets, the index is set to -1;
vector<int> ClusterSequence::particle_jet_indices(
                        const vector<PseudoJet> & jets_in) const {

  vector<int> indices(n_particles());

  // first label all particles as not belonging to any jets
  for (unsigned ipart = 0; ipart < n_particles(); ipart++) 
    indices[ipart] = -1;

  // then for each of the jets relabel its consituents as belonging to
  // that jet
  for (unsigned ijet = 0; ijet < jets_in.size(); ijet++) {

    vector<PseudoJet> jet_constituents(constituents(jets_in[ijet]));

    for (unsigned ip = 0; ip < jet_constituents.size(); ip++) {
      // a safe (if slightly redundant) way of getting the particle
      // index (for initial particles it is actually safe to assume
      // ipart=iclust).
      unsigned iclust = jet_constituents[ip].cluster_hist_index();
      unsigned ipart = history()[iclust].jetp_index;
      indices[ipart] = ijet;
    }
  }

  return indices;
}


//----------------------------------------------------------------------
// recursive routine that adds on constituents of jet to the subjet_vector
void ClusterSequence::add_constituents (
           const PseudoJet & jet, vector<PseudoJet> & subjet_vector) const {
  // find out position in cluster history
  int i = jet.cluster_hist_index();
  int parent1 = _history[i].parent1;
  int parent2 = _history[i].parent2;

  if (parent1 == InexistentParent) {
    // It is an original particle (labelled by its parent having value
    // InexistentParent), therefore add it on to the subjet vector
    // Note: we add the initial particle and not simply 'jet' so that
    //       calling add_constituents with a subtracted jet containing
    //       only one particle will work.
    subjet_vector.push_back(_jets[i]);
    return;
  } 

  // add parent 1
  add_constituents(_jets[_history[parent1].jetp_index], subjet_vector);

  // see if parent2 is a real jet; if it is then add its constituents
  if (parent2 != BeamJet) {
    add_constituents(_jets[_history[parent2].jetp_index], subjet_vector);
  }
}



//----------------------------------------------------------------------
// initialise the history in a standard way
void ClusterSequence::_add_step_to_history (
	       const int step_number, const int parent1, 
	       const int parent2, const int jetp_index,
	       const double dij) {

  history_element element;
  element.parent1 = parent1;
  element.parent2 = parent2;
  element.jetp_index = jetp_index;
  element.child = Invalid;
  element.dij   = dij;
  element.max_dij_so_far = max(dij,_history[_history.size()-1].max_dij_so_far);
  _history.push_back(element);

  int local_step = _history.size()-1;
  assert(local_step == step_number);

  assert(parent1 >= 0);
  _history[parent1].child = local_step;
  if (parent2 >= 0) {_history[parent2].child = local_step;}

  // get cross-referencing right from PseudoJets
  if (jetp_index != Invalid) {
    assert(jetp_index >= 0);
    //cout << _jets.size() <<" "<<jetp_index<<"\n";
    _jets[jetp_index].set_cluster_hist_index(local_step);
    _set_structure_shared_ptr(_jets[jetp_index]);
  }

  if (_writeout_combinations) {
    cout << local_step << ": " 
	 << parent1 << " with " << parent2
	 << "; y = "<< dij<<endl;
  }

}




//======================================================================
// Return an order in which to read the history such that _history[order[i]] 
// will always correspond to the same set of consituent particles if 
// two branching histories are equivalent in terms of the particles
// contained in any given pseudojet.
vector<int> ClusterSequence::unique_history_order() const {

  // first construct an array that will tell us the lowest constituent
  // of a given jet -- this will always be one of the original
  // particles, whose order is well defined and so will help us to
  // follow the tree in a unique manner.
  valarray<int> lowest_constituent(_history.size());
  int hist_n = _history.size();
  lowest_constituent = hist_n; // give it a large number
  for (int i = 0; i < hist_n; i++) {
    // sets things up for the initial partons
    lowest_constituent[i] = min(lowest_constituent[i],i); 
    // propagates them through to the children of this parton
    if (_history[i].child > 0) lowest_constituent[_history[i].child] 
      = min(lowest_constituent[_history[i].child],lowest_constituent[i]);
  }

  // establish an array for what we have and have not extracted so far
  valarray<bool> extracted(_history.size()); extracted = false;
  vector<int> unique_tree;
  unique_tree.reserve(_history.size());

  // now work our way through the tree
  for (unsigned i = 0; i < n_particles(); i++) {
    if (!extracted[i]) {
      unique_tree.push_back(i);
      extracted[i] = true;
      _extract_tree_children(i, extracted, lowest_constituent, unique_tree);
    }
  }

  return unique_tree;
}

//======================================================================
// helper for unique_history_order
void ClusterSequence::_extract_tree_children(
       int position, 
       valarray<bool> & extracted, 
       const valarray<int> & lowest_constituent,
       vector<int> & unique_tree) const {
  if (!extracted[position]) {
    // that means we may have unidentified parents around, so go and
    // collect them (extracted[position]) will then be made true)
    _extract_tree_parents(position,extracted,lowest_constituent,unique_tree);
  } 
  
  // now look after the children...
  int child = _history[position].child;
  if (child  >= 0) _extract_tree_children(child,extracted,lowest_constituent,unique_tree);
}


//======================================================================
// return the list of unclustered particles
vector<PseudoJet> ClusterSequence::unclustered_particles() const {
  vector<PseudoJet> unclustered;
  for (unsigned i = 0; i < n_particles() ; i++) {
    if (_history[i].child == Invalid) 
      unclustered.push_back(_jets[_history[i].jetp_index]);
  }
  return unclustered;
}

//======================================================================
/// Return the list of pseudojets in the ClusterSequence that do not
/// have children (and are not among the inclusive jets). They may
/// result from a clustering step or may be one of the pseudojets
/// returned by unclustered_particles().
vector<PseudoJet> ClusterSequence::childless_pseudojets() const {
  vector<PseudoJet> unclustered;
  for (unsigned i = 0; i < _history.size() ; i++) {
    if ((_history[i].child == Invalid) && (_history[i].parent2 != BeamJet))
      unclustered.push_back(_jets[_history[i].jetp_index]);
  }
  return unclustered;
}



//----------------------------------------------------------------------
// returns true if the cluster sequence contains this jet (i.e. jet's
// structure is this cluster sequence's and the cluster history index
// is in a consistent range)
bool ClusterSequence::contains(const PseudoJet & jet) const {
  return jet.cluster_hist_index() >= 0 
    &&   jet.cluster_hist_index() < int(_history.size())
    &&   jet.has_valid_cluster_sequence()
    &&   jet.associated_cluster_sequence() == this;
}



//======================================================================
// helper for unique_history_order
void ClusterSequence::_extract_tree_parents(
       int position, 
       valarray<bool> & extracted, 
       const valarray<int> & lowest_constituent,
       vector<int> & unique_tree) const {

  if (!extracted[position]) {
    int parent1 = _history[position].parent1;
    int parent2 = _history[position].parent2;
    // where relevant order parents so that we will first treat the
    // one containing the smaller "lowest_constituent"
    if (parent1 >= 0 && parent2 >= 0) {
      if (lowest_constituent[parent1] > lowest_constituent[parent2]) 
	std::swap(parent1, parent2);
    }
    // then actually run through the parents to extract the constituents...
    if (parent1 >= 0 && !extracted[parent1]) 
      _extract_tree_parents(parent1,extracted,lowest_constituent,unique_tree);
    if (parent2 >= 0 && !extracted[parent2]) 
      _extract_tree_parents(parent2,extracted,lowest_constituent,unique_tree);
    // finally declare this position to be accounted for and push it
    // onto our list.
    unique_tree.push_back(position);
    extracted[position] = true;
  }
}


//======================================================================
/// carries out the bookkeeping associated with the step of recombining
/// jet_i and jet_j (assuming a distance dij) and returns the index
/// of the recombined jet, newjet_k.
void ClusterSequence::_do_ij_recombination_step(
                               const int jet_i, const int jet_j, 
			       const double dij, 
			       int & newjet_k) {

  // Create the new jet by recombining the first two.
  //
  // For efficiency reasons, use a ctr that initialises only the
  // shared pointers, since the rest of the info will anyway be dealt
  // with by the recombiner.
  PseudoJet newjet(false); 
  _jet_def.recombiner()->recombine(_jets[jet_i], _jets[jet_j], newjet);
  _jets.push_back(newjet);
  // original version...
  //_jets.push_back(_jets[jet_i] + _jets[jet_j]);

  // get its index
  newjet_k = _jets.size()-1;

  // get history index
  int newstep_k = _history.size();
  // and provide jet with the info
  _jets[newjet_k].set_cluster_hist_index(newstep_k);

  // finally sort out the history 
  int hist_i = _jets[jet_i].cluster_hist_index();
  int hist_j = _jets[jet_j].cluster_hist_index();

  _add_step_to_history(newstep_k, min(hist_i, hist_j), max(hist_i,hist_j),
		       newjet_k, dij);

}


//======================================================================
/// carries out the bookkeeping associated with the step of recombining
/// jet_i with the beam
void ClusterSequence::_do_iB_recombination_step(
				  const int jet_i, const double diB) {
  // get history index
  int newstep_k = _history.size();

  // recombine the jet with the beam
  _add_step_to_history(newstep_k,_jets[jet_i].cluster_hist_index(),BeamJet,
		       Invalid, diB);

}


// make sure the static member _changed_strategy_warning is defined. 
LimitedWarning ClusterSequence::_changed_strategy_warning;


//----------------------------------------------------------------------
void ClusterSequence::_set_structure_shared_ptr(PseudoJet & j) {
  j.set_structure_shared_ptr(_structure_shared_ptr);
  // record the use count of the structure shared point to help
  // in case we want to ask the CS to handle its own memory
  _update_structure_use_count();
}


//----------------------------------------------------------------------
void ClusterSequence::_update_structure_use_count() {
  // record the use count of the structure shared point to help
  // in case we want to ask the CS to handle its own memory
  _structure_use_count_after_construction = _structure_shared_ptr.use_count();
}

//----------------------------------------------------------------------
/// by calling this routine you tell the ClusterSequence to delete
/// itself when all the Pseudojets associated with it have gone out
/// of scope. 
void ClusterSequence::delete_self_when_unused() {
  // the trick we use to handle this is to modify the use count; 
  // that way the structure will be deleted when there are no external
  // objects left associated the CS and the structure's destructor will then
  // look after deleting the cluster sequence
  
  // first make sure that there is at least one other object
  // associated with the CS
  int new_count = _structure_shared_ptr.use_count() - _structure_use_count_after_construction;
  if (new_count <= 0) {
    throw Error("delete_self_when_unused may only be called if at least one object outside the CS (e.g. a jet) is already associated with the CS");
  }

  _structure_shared_ptr.set_count(new_count);
  _deletes_self_when_unused = true;
}


FASTJET_END_NAMESPACE

