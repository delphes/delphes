//FJSTARTHEADER
// $Id: JetDefinition.cc 3677 2014-09-09 22:45:25Z soyez $
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

#include "fastjet/JetDefinition.hh"
#include "fastjet/Error.hh"
#include "fastjet/CompositeJetStructure.hh"
#include<sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

const double JetDefinition::max_allowable_R = 1000.0;

//----------------------------------------------------------------------
// [NB: implementation was getting complex, so in 2.4-devel moved it
//  from .hh to .cc]
JetDefinition::JetDefinition(JetAlgorithm jet_algorithm_in, 
			     double R_in, 
			     Strategy strategy_in,
			     RecombinationScheme recomb_scheme_in,
                             int nparameters) :
  _jet_algorithm(jet_algorithm_in), _Rparam(R_in), _strategy(strategy_in) {

  // set R parameter or ensure its sensibleness, as appropriate
  if (_jet_algorithm == ee_kt_algorithm) {
    _Rparam = 4.0; // introduce a fictional R that ensures that
                   // our clustering sequence will not produce
                   // "beam" jets except when only a single particle remains.
                   // Any value > 2 would have done here
  } else {
    // We maintain some limit on R because particles with pt=0, m=0
    // can have rapidities O(100000) and one doesn't want the
    // clustering to start including them as if their rapidities were
    // physical.
    if (R_in > max_allowable_R) {
      ostringstream oss;
      oss << "Requested R = " << R_in << " for jet definition is larger than max_allowable_R = " << max_allowable_R;
      throw Error(oss.str());
    }
  }

  // cross-check the number of parameters that were declared in setting up the
  // algorithm (passed internally from the public constructors)
  unsigned int nparameters_expected = n_parameters_for_algorithm(jet_algorithm_in);
  if (nparameters != (int) nparameters_expected){
    ostringstream oss;
    oss << "The jet algorithm you requested ("
        << jet_algorithm_in << ") should be constructed with " << nparameters_expected 
        << " parameter(s) but was called with " << nparameters << " parameter(s)\n";
    throw Error(oss.str()); 
  }

  // make sure the strategy requested is sensible
  assert (_strategy  != plugin_strategy);

  _plugin = NULL;
  set_recombination_scheme(recomb_scheme_in);
  set_extra_param(0.0); // make sure it's defined
}


//----------------------------------------------------------------------
// returns true if the jet definition involves an algorithm
// intended for use on a spherical geometry (e.g. e+e- algorithms,
// as opposed to most pp algorithms, which use a cylindrical,
// rapidity-phi geometry).
bool JetDefinition::is_spherical() const {
  if (jet_algorithm() == plugin_algorithm) {
    return plugin()->is_spherical();
  } else {
    return (jet_algorithm() == ee_kt_algorithm ||  // as of 2013-02-14, the two
            jet_algorithm() == ee_genkt_algorithm  // native spherical algorithms
            );
  }
}

//----------------------------------------------------------------------
string JetDefinition::description() const {
  ostringstream name;
  
  name << description_no_recombiner();

  if ((jet_algorithm() == plugin_algorithm) || (jet_algorithm() == undefined_jet_algorithm)){
    return name.str();
  }

  if (n_parameters_for_algorithm(jet_algorithm()) == 0)
    name << " with ";
  else 
    name << " and ";
  name << recombiner()->description();

  return name.str();
}

//----------------------------------------------------------------------
string JetDefinition::description_no_recombiner() const {
  
  ostringstream name;
  if (jet_algorithm() == plugin_algorithm) {
    return plugin()->description();
  } else if (jet_algorithm() == undefined_jet_algorithm) {
    return "uninitialised JetDefinition (jet_algorithm=undefined_jet_algorithm)" ;
  }

  name << algorithm_description(jet_algorithm());
  switch (n_parameters_for_algorithm(jet_algorithm())){
  case 0: name << " (NB: no R)"; break;
  case 1: name << " with R = " << R(); break; // the parameter is always R
  case 2: 
    // the 1st parameter is always R
    name << " with R = " << R();
    // the 2nd depends on the algorithm
    if (jet_algorithm() == cambridge_for_passive_algorithm){
      name << "and a special hack whereby particles with kt < " 
           << extra_param() << "are treated as passive ghosts";
    } else {
      name << ", p = " << extra_param();
    }
  };

  return name.str();
}

//----------------------------------------------------------------------
string JetDefinition::algorithm_description(const JetAlgorithm jet_alg){
  ostringstream name;
  switch (jet_alg){
  case plugin_algorithm:                return "plugin algorithm";
  case kt_algorithm:                    return "Longitudinally invariant kt algorithm";
  case cambridge_algorithm:             return "Longitudinally invariant Cambridge/Aachen algorithm";
  case antikt_algorithm:                return "Longitudinally invariant anti-kt algorithm";
  case genkt_algorithm:                 return "Longitudinally invariant generalised kt algorithm";
  case cambridge_for_passive_algorithm: return "Longitudinally invariant Cambridge/Aachen algorithm";
  case ee_kt_algorithm:                 return "e+e- kt (Durham) algorithm (NB: no R)";
  case ee_genkt_algorithm:              return "e+e- generalised kt algorithm";
  case undefined_jet_algorithm:         return "undefined jet algorithm";
  default:
    throw Error("JetDefinition::algorithm_description(): unrecognized jet_algorithm");
  };
}

//----------------------------------------------------------------------
unsigned int JetDefinition::n_parameters_for_algorithm(const JetAlgorithm jet_alg){
  switch (jet_alg) {
  case ee_kt_algorithm:    return 0;
  case genkt_algorithm:
  case ee_genkt_algorithm: return 2;
  default:                 return 1;
  };
}

//----------------------------------------------------------------------
void JetDefinition::set_recombination_scheme(
                               RecombinationScheme recomb_scheme) {
  _default_recombiner = JetDefinition::DefaultRecombiner(recomb_scheme);

  // do not forget to delete the existing recombiner if needed
  if (_shared_recombiner()) _shared_recombiner.reset();

  _recombiner = 0;
}

void JetDefinition::set_recombiner(const JetDefinition &other_jet_def){
  // make sure the "invariants" of the other jet def are sensible
  assert(other_jet_def._recombiner || 
         other_jet_def.recombination_scheme() != external_scheme);

  // first treat the situation where we're using the default recombiner
  if (other_jet_def._recombiner == 0){
    set_recombination_scheme(other_jet_def.recombination_scheme());
    return;
  }

  // in other cases, copy the pointer to the recombiner
  _recombiner = other_jet_def._recombiner;
  // set the default recombiner appropriately
  _default_recombiner = DefaultRecombiner(external_scheme);
  // and set the _shared_recombiner to the same state
  // as in the other_jet_def, whatever that was
  _shared_recombiner.reset(other_jet_def._shared_recombiner);

  // NB: it is tempting to go via set_recombiner and then to sort
  // out the shared part, but this would be dangerous in the
  // specific (rare?) case where other_jet_def is the same as this
  // it deletes_recombiner_when_unused. In that case the shared
  // pointer reset would delete the recombiner.
}


// returns true if the current jet definitions shares the same
// recombiner as teh one passed as an argument
bool JetDefinition::has_same_recombiner(const JetDefinition &other_jd) const{
  // first make sure that they have the same recombination scheme
  const RecombinationScheme & scheme = recombination_scheme();
  if (other_jd.recombination_scheme() != scheme) return false;

  // if the scheme is "external", also check that they have the same
  // recombiner
  return (scheme != external_scheme) 
    || (recombiner() == other_jd.recombiner());
}

/// causes the JetDefinition to handle the deletion of the
/// recombiner when it is no longer used
void JetDefinition::delete_recombiner_when_unused(){
  if (_recombiner == 0){
    throw Error("tried to call JetDefinition::delete_recombiner_when_unused() for a JetDefinition without a user-defined recombination scheme");
  } else if (_shared_recombiner.get()) {
    throw Error("Error in JetDefinition::delete_recombiner_when_unused: the recombiner is already scheduled for deletion when unused (or was already set as shared)");
  }

  _shared_recombiner.reset(_recombiner);
}

/// allows to let the JetDefinition handle the deletion of the
/// plugin when it is no longer used
void JetDefinition::delete_plugin_when_unused(){
  if (_plugin == 0){
    throw Error("tried to call JetDefinition::delete_plugin_when_unused() for a JetDefinition without a plugin");
  }

  _plugin_shared.reset(_plugin);
}



string JetDefinition::DefaultRecombiner::description() const {
  switch(_recomb_scheme) {
  case E_scheme:
    return "E scheme recombination";
  case pt_scheme:
    return "pt scheme recombination";
  case pt2_scheme:
    return "pt2 scheme recombination";
  case Et_scheme:
    return "Et scheme recombination";
  case Et2_scheme:
    return "Et2 scheme recombination";
  case BIpt_scheme:
    return "boost-invariant pt scheme recombination";
  case BIpt2_scheme:
    return "boost-invariant pt2 scheme recombination";
  case WTA_pt_scheme:
    return "pt-ordered Winner-Takes-All recombination";
  // Energy-ordering can lead to dangerous situations with particles at
  // rest. We instead implement the WTA_modp_scheme
  //
  //   case WTA_E_scheme:
  //     return "energy-ordered Winner-Takes-All recombination";
  case WTA_modp_scheme:
    return "|3-momentum|-ordered Winner-Takes-All recombination";
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }
}


void JetDefinition::DefaultRecombiner::recombine(
           const PseudoJet & pa, const PseudoJet & pb,
           PseudoJet & pab) const {
  
  double weighta, weightb;

  switch(_recomb_scheme) {
  case E_scheme:
    // a call to reset turns out to be somewhat more efficient
    // than a sum and assignment
    //pab = pa + pb; 
    pab.reset(pa.px()+pb.px(),
    	      pa.py()+pb.py(),
    	      pa.pz()+pb.pz(),
    	      pa.E ()+pb.E ());
    return;
  // all remaining schemes are massless recombinations and locally
  // we just set weights, while the hard work is done below...
  case pt_scheme:
  case Et_scheme:
  case BIpt_scheme:
    weighta = pa.perp(); 
    weightb = pb.perp();
    break;
  case pt2_scheme:
  case Et2_scheme:
  case BIpt2_scheme:
    weighta = pa.perp2(); 
    weightb = pb.perp2();
    break;
  case WTA_pt_scheme:{
    const PseudoJet & phard = (pa.pt2() >= pb.pt2()) ? pa : pb;
    /// keep y,phi and m from the hardest, sum pt
    pab.reset_PtYPhiM(pa.pt()+pb.pt(), 
                      phard.rap(), phard.phi(), phard.m());
    return;}
  // Energy-ordering can lead to dangerous situations with particles at
  // rest. We instead implement the WTA_modp_scheme
  //
  //   case WTA_E_scheme:{
  //     const PseudoJet & phard = (pa.E() >= pb.E()) ? pa : pb;
  //     /// keep 3-momentum direction and mass from the hardest, sum energies
  //     ///
  //     /// If the particle with the largest energy is at rest, the sum
  //     /// remains at rest, implying that the mass of the sum is larger
  //     /// than the mass of pa.
  //     double Eab = pa.E() + pb.E();
  //     double scale = (phard.modp2()==0.0)
  //       ? 0.0
  //       : sqrt((Eab*Eab - phard.m2())/phard.modp2());
  //     pab.reset(phard.px()*scale, phard.py()*scale, phard.pz()*scale, Eab);
  //     return;}
  case WTA_modp_scheme:{
    // Note: we need to compute both a and b modp. And we need pthard
    // and its modp. If we want to avoid repeating the test and do
    // only 2 modp calculations, we'd have to duplicate the code (or
    // use a pair<const PJ&>). An alternative is to write modp_soft as
    // modp_ab-modp_hard but this could suffer from larger rounding
    // errors
    bool a_hardest = (pa.modp2() >= pb.modp2());
    const PseudoJet & phard = a_hardest ? pa : pb;
    const PseudoJet & psoft = a_hardest ? pb : pa;
    /// keep 3-momentum direction and mass from the hardest, sum modp
    ///
    /// If the hardest particle is at rest, the sum remains at rest
    /// (the energy of the sum is therefore the mass of pa)
    double modp_hard = phard.modp();
    double modp_ab = modp_hard + psoft.modp();
    if (phard.modp2()==0.0){
      pab.reset(0.0, 0.0, 0.0, phard.m());
    } else {
      double scale = modp_ab/modp_hard;
      pab.reset(phard.px()*scale, phard.py()*scale, phard.pz()*scale,
                sqrt(modp_ab*modp_ab + phard.m2()));
    }
    return;}
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }

  double perp_ab = pa.perp() + pb.perp();
  if (perp_ab != 0.0) { // weights also non-zero...
    double y_ab    = (weighta * pa.rap() + weightb * pb.rap())/(weighta+weightb);
    
    // take care with periodicity in phi...
    double phi_a = pa.phi(), phi_b = pb.phi();
    if (phi_a - phi_b > pi)  phi_b += twopi;
    if (phi_a - phi_b < -pi) phi_b -= twopi;
    double phi_ab = (weighta * phi_a + weightb * phi_b)/(weighta+weightb);

    // this is much more efficient...
    pab.reset_PtYPhiM(perp_ab,y_ab,phi_ab);
    // pab = PseudoJet(perp_ab*cos(phi_ab),
    // 		    perp_ab*sin(phi_ab),
    // 		    perp_ab*sinh(y_ab),
    // 		    perp_ab*cosh(y_ab));
  } else { // weights are zero
    //pab = PseudoJet(0.0,0.0,0.0,0.0);
    pab.reset(0.0, 0.0, 0.0, 0.0);
  }
}


void JetDefinition::DefaultRecombiner::preprocess(PseudoJet & p) const {
  switch(_recomb_scheme) {
  case E_scheme:
  case BIpt_scheme:
  case BIpt2_scheme:
  case WTA_pt_scheme:
  //case WTA_E_scheme:
  case WTA_modp_scheme:
    break;
  case pt_scheme:
  case pt2_scheme:
    {
      // these schemes (as in the ktjet implementation) need massless
      // initial 4-vectors with essentially E=|p|.
      double newE = sqrt(p.perp2()+p.pz()*p.pz());
      p.reset_momentum(p.px(), p.py(), p.pz(), newE);
      // FJ2.x version
      // int    user_index = p.user_index();
      // p = PseudoJet(p.px(), p.py(), p.pz(), newE);
      // p.set_user_index(user_index);
    }
    break;
  case Et_scheme:
  case Et2_scheme:
    {
      // these schemes (as in the ktjet implementation) need massless
      // initial 4-vectors with essentially E=|p|.
      double rescale = p.E()/sqrt(p.perp2()+p.pz()*p.pz());
      p.reset_momentum(rescale*p.px(), rescale*p.py(), rescale*p.pz(), p.E());
      // FJ2.x version
      // int    user_index = p.user_index();
      // p = PseudoJet(rescale*p.px(), rescale*p.py(), rescale*p.pz(), p.E());
      // p.set_user_index(user_index);
    }
    break;
  default:
    ostringstream err;
    err << "DefaultRecombiner: unrecognized recombination scheme " 
        << _recomb_scheme;
    throw Error(err.str());
  }
}

void JetDefinition::Plugin::set_ghost_separation_scale(double /*scale*/) const {
  throw Error("set_ghost_separation_scale not supported");
}



//-------------------------------------------------------------------------------
// helper functions to build a jet made of pieces
//
// This is the extended version with support for a user-defined
// recombination-scheme
// -------------------------------------------------------------------------------

// build a "CompositeJet" from the vector of its pieces
//
// the user passes the reciombination scheme used to "sum" the pieces.
PseudoJet join(const vector<PseudoJet> & pieces, const JetDefinition::Recombiner & recombiner){
  // compute the total momentum
  //--------------------------------------------------
  PseudoJet result;  // automatically initialised to 0
  if (pieces.size()>0){
    result = pieces[0];
    for (unsigned int i=1; i<pieces.size(); i++)
      recombiner.plus_equal(result, pieces[i]);
  }

  // attach a CompositeJetStructure to the result
  //--------------------------------------------------
  CompositeJetStructure *cj_struct = new CompositeJetStructure(pieces, &recombiner);

  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));

  return result;
}

// build a "CompositeJet" from a single PseudoJet
PseudoJet join(const PseudoJet & j1, 
	       const JetDefinition::Recombiner & recombiner){
  return join(vector<PseudoJet>(1,j1), recombiner);
}

// build a "CompositeJet" from two PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join(pieces, recombiner);
}

// build a "CompositeJet" from 3 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, 
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join(pieces, recombiner);
}

// build a "CompositeJet" from 4 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4,
	       const JetDefinition::Recombiner & recombiner){
  vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join(pieces, recombiner);
}


 

FASTJET_END_NAMESPACE
