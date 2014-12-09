///////////////////////////////////////////////////////////////////////////////
// File: split_merge.cpp                                                     //
// Description: source file for splitting/merging (contains the CJet class)  //
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
// $Revision:: 370                                                          $//
// $Date:: 2014-09-04 17:03:15 +0200 (Thu, 04 Sep 2014)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "split_merge.h"
#include "siscone_error.h"
#include "momentum.h"
#include <limits>   // for max
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>

namespace siscone{

using namespace std;

/********************************************************
 * class Cjet implementation                            *
 * real Jet information.                                *
 * This class contains information for one single jet.  *
 * That is, first, its momentum carrying information    *
 * about its centre and pT, and second, its particle    *
 * contents                                             *
 ********************************************************/
// default ctor
//--------------
Cjet::Cjet(){
  n = 0;
  v = Cmomentum();
  pt_tilde = 0.0;
  sm_var2 = 0.0;
}

// default dtor
//--------------
Cjet::~Cjet(){

}

// ordering of jets in pt (e.g. used in final jets ordering)
//-----------------------------------------------------------
bool jets_pt_less(const Cjet &j1, const Cjet &j2){
  return j1.v.perp2() > j2.v.perp2();
}


/********************************************************
 * Csplit_merge_ptcomparison implementation             *
 * This deals with the ordering of the jets candidates  *
 ********************************************************/

// odering of two jets
// The variable of the ordering is pt or mt 
// depending on 'split_merge_scale' choice
//
// with EPSILON_SPLITMERGE defined, this attempts to identify
// delicate cases where two jets have identical momenta except for
// softish particles -- the difference of pt's may not be correctly
// identified normally and so lead to problems for the fate of the
// softish particle.
//
// NB: there is a potential issue in momentum-conserving events,
// whereby the harder of two jets becomes ill-defined when a soft
// particle is emitted --- this may have a knock-on effect on
// subsequent split-merge steps in cases with sufficiently large R
// (but we don't know what the limit is...)
//------------------------------------------------------------------
bool Csplit_merge_ptcomparison::operator ()(const Cjet &jet1, const Cjet &jet2) const{
  double q1, q2;

  // compute the value for comparison for both jets
  // This depends on the choice of variable (mt is the default)
  q1 = jet1.sm_var2;
  q2 = jet2.sm_var2;

  bool res = q1 > q2;

  // if we enable the refined version of the comparison (see defines.h),
  // we compute the difference more precisely when the two jets are very
  // close in the ordering variable.
#ifdef EPSILON_SPLITMERGE
  if ( (fabs(q1-q2) < EPSILON_SPLITMERGE*max(q1,q2)) &&
       (jet1.v.ref != jet2.v.ref) ) {
    // get the momentum of the difference
    Cmomentum difference;
    double pt_tilde_difference;
    get_difference(jet1,jet2,&difference,&pt_tilde_difference);
    
    // use the following relation: pt1^2 - pt2^2 = (pt1+pt2)*(pt1-pt2)
    double qdiff;
    Cmomentum sum = jet1.v ;
    sum +=  jet2.v;
    double pt_tilde_sum = jet1.pt_tilde + jet2.pt_tilde;
    
    // depending on the choice of ordering variable, set the result
    switch (split_merge_scale){
    case SM_mt:
      qdiff = sum.E*difference.E - sum.pz*difference.pz;
      break;
    case SM_pt:
      qdiff = sum.px*difference.px + sum.py*difference.py;
      break;
    case SM_pttilde:  
      qdiff = pt_tilde_sum*pt_tilde_difference;
      break;
    case SM_Et:
      // diff = E^2 (dpt^2 pz^2- pt^2 dpz^2)
      //      + dE^2 (pt^2+pz^2) pt2^2
      // where, unless explicitely specified the variables
      // refer to the first jet or differences jet1-jet2.
      qdiff = jet1.v.E*jet1.v.E*
	((sum.px*difference.px + sum.py*difference.py)*jet1.v.pz*jet1.v.pz
	 -jet1.v.perp2()*sum.pz*difference.pz)
	+sum.E*difference.E*(jet1.v.perp2()+jet1.v.pz*jet1.v.pz)*jet2.v.perp2();
      break;
    default:
      throw Csiscone_error("Unsupported split-merge scale choice: "
			   + SM_scale_name());
    }
    res = qdiff > 0;
  }
#endif  // EPSILON_SPLITMERGE

  return res;
}


/// return a name for the sm scale choice 
/// NB: used internally and also by fastjet
std::string split_merge_scale_name(Esplit_merge_scale sms) {
  switch(sms) {
  case SM_pt:
    return "pt (IR unsafe)";
  case SM_Et:
    return "Et (boost dep.)";
  case SM_mt:
    return "mt (IR safe except for pairs of identical decayed heavy particles)";
  case SM_pttilde:
    return "pttilde (scalar sum of pt's)";
  default:
    return "[SM scale without a name]";
  }
}


// get the difference between 2 jets
//  - j1         first jet
//  - j2         second jet
//  - v          jet1-jet2
//  - pt_tilde   jet1-jet2 pt_tilde
// return true if overlapping, false if disjoint
//-----------------------------------------------
void Csplit_merge_ptcomparison::get_difference(const Cjet &j1, const Cjet &j2, Cmomentum *v, double *pt_tilde) const {
  int i1,i2;

  // initialise
  i1=i2=0;
  *v = Cmomentum();
  *pt_tilde = 0.0;

  // compute overlap
  // at the same time, we store union in indices
  do{
    if (j1.contents[i1]==j2.contents[i2]) {
      i1++;
      i2++;
    } else if (j1.contents[i1]<j2.contents[i2]){
      (*v) += (*particles)[j1.contents[i1]];
      (*pt_tilde) += (*pt)[j1.contents[i1]];
      i1++;
    } else if (j1.contents[i1]>j2.contents[i2]){
      (*v) -= (*particles)[j2.contents[i2]];
      (*pt_tilde) -= (*pt)[j2.contents[i2]];
      i2++;
    } else {
      throw Csiscone_error("get_non_overlap reached part it should never have seen...");
    }
  } while ((i1<j1.n) && (i2<j2.n));

  // deal with particles at the end of the list...
  while (i1 < j1.n) {
    (*v) += (*particles)[j1.contents[i1]];
    (*pt_tilde) += (*pt)[j1.contents[i1]];
    i1++;
  }
  while (i2 < j2.n) {
    (*v) -= (*particles)[j2.contents[i2]];
    (*pt_tilde) -= (*pt)[j2.contents[i2]];
    i2++;
  }
}


/********************************************************
 * class Csplit_merge implementation                    *
 * Class used to split and merge jets.                  *
 ********************************************************/
// default ctor
//--------------
Csplit_merge::Csplit_merge(){
  merge_identical_protocones = false;
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
#ifdef MERGE_IDENTICAL_PROTOCONES_DEFAULT_TRUE
  merge_identical_protocones = true;
#endif
#endif
  _user_scale = NULL;
  indices = NULL;

  // ensure that ptcomparison points to our set of particles (though params not correct)
  ptcomparison.particles = &particles;
  ptcomparison.pt = &pt;
  candidates.reset(new multiset<Cjet,Csplit_merge_ptcomparison>(ptcomparison));

  // no hardest cut (col-unsafe)
  SM_var2_hardest_cut_off = -numeric_limits<double>::max();

  // no pt cutoff for the particles to put in p_uncol_hard
  stable_cone_soft_pt2_cutoff = -1.0;

  // no pt-weighted splitting
  use_pt_weighted_splitting = false;
}


// default dtor
//--------------
Csplit_merge::~Csplit_merge(){
  full_clear();
}


// initialisation function
//  - _particles  list of particles
//  - protocones  list of protocones (initial jet candidates)
//  - R2          cone radius (squared)
//  - ptmin       minimal pT allowed for jets
//-------------------------------------------------------------
int Csplit_merge::init(vector<Cmomentum> & /*_particles*/, vector<Cmomentum> *protocones, double R2, double ptmin){
  // browse protocones
  return add_protocones(protocones, R2, ptmin);
}


// initialisation function for particle list
//  - _particles  list of particles
//-------------------------------------------------------------
int Csplit_merge::init_particles(vector<Cmomentum> &_particles){
  full_clear();

  // compute the list of particles
  // here, we do not need to throw away particles 
  // with infinite rapidity (colinear with the beam)
  particles = _particles;
  n = particles.size();

  // build the vector of particles' pt
  pt.resize(n);
  for (int i=0;i<n;i++)
    pt[i] = particles[i].perp();

  // ensure that ptcomparison points to our set of particles (though params not correct)
  ptcomparison.particles = &particles;
  ptcomparison.pt = &pt;

  // set up the list of particles left.
  init_pleft();

  indices = new int[n];

  return 0;
}


// build initial list of remaining particles
//------------------------------------------
int Csplit_merge::init_pleft(){
  // at this level, we only rule out particles with 
  // infinite rapidity
  // in the parent particle list, index contain the run 
  // at which particles are puts in jets:
  //  - -1 means infinity rapidity
  //  -  0 means not included
  //  -  i mean included at run i
  int i,j;

  // copy particles removing the ones with infinite rapidity
  // determine min,max eta
  j=0;
  double eta_min=0.0;  /// for the Ceta_phi_range static member!
  double eta_max=0.0;  /// for the Ceta_phi_range static member!
  p_remain.clear();
  for (i=0;i<n;i++){
    // set ref for checkxor
    particles[i].ref.randomize();

    // check if rapidity is not infinite or ill-defined
    if (fabs(particles[i].pz) < (particles[i].E)){
      p_remain.push_back(particles[i]);
      // set up parent index for tracability
      p_remain[j].parent_index = i;
      // left particles are marked with a 1
      // IMPORTANT NOTE: the meaning of index in p_remain is
      //   somehow different as in the initial particle list.
      //   here, within one pass, we use the index to track whether
      //   a particle is included in the current pass (index set to 0
      //   in add_protocones) or still remain (index still 1)
      p_remain[j].index = 1;

      j++;
      // set up parent-particle index
      particles[i].index = 0;

      eta_min = min(eta_min, particles[i].eta);
      eta_max = max(eta_max, particles[i].eta);
    } else {
      particles[i].index = -1;
    }
  }
  n_left = p_remain.size();
  n_pass = 0;

  Ceta_phi_range epr;
  epr.eta_min = eta_min-0.01;
  epr.eta_max = eta_max+0.01;

  merge_collinear_and_remove_soft();

  return 0;
}


// partial clearance
// we want to keep   particle list and indices
// for future usage, so do not clear it !
// this is done in full_clear
//----------------------------------------
int Csplit_merge::partial_clear(){
  // release jets

  // set up the auto_ptr for the multiset with the _current_ state of
  // ptcomparison (which may have changed since we constructed the
  // class)
  candidates.reset(new multiset<Cjet,Csplit_merge_ptcomparison>(ptcomparison));

  // start off with huge number
  most_ambiguous_split = numeric_limits<double>::max();

  jets.clear();
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if (merge_identical_protocones)
    cand_refs.clear();
#endif

  p_remain.clear();

  return 0;
}


// full clearance
//----------------
int Csplit_merge::full_clear(){
  partial_clear();

  // clear previously allocated memory
  if (indices != NULL){
    delete[] indices;
  }
  particles.clear();

  return 0;
}


// build the list 'p_uncol_hard' from p_remain by clustering collinear particles
// note that thins in only used for stable-cone detection 
// so the parent_index field is unnecessary
//-------------------------------------------------------------------------
int Csplit_merge::merge_collinear_and_remove_soft(){
  int i,j;
  vector<Cmomentum> p_sorted;
  bool collinear;
  double dphi;

  p_uncol_hard.clear();

  // we first sort the particles according to their rapidity
  for (i=0;i<n_left;i++)
    p_sorted.push_back(p_remain[i]);
  sort(p_sorted.begin(), p_sorted.end(), momentum_eta_less);

  // then we cluster particles looping over the particles in the following way
  // if (a particle i has same eta-phi a one after (j))
  // then add momentum i to j
  // else add i to the p_uncol_hard list
  i = 0;
  while (i<n_left){
    // check if the particle passes the stable_cone_soft_pt2_cutoff
    if (p_sorted[i].perp2()<stable_cone_soft_pt2_cutoff) {
      i++;
      continue;
    }

    // check if there is(are) particle(s) with the 'same' eta
    collinear = false;
    j=i+1;
    while ((j<n_left) && (fabs(p_sorted[j].eta-p_sorted[i].eta)<EPSILON_COLLINEAR) && (!collinear)){
      dphi = fabs(p_sorted[j].phi-p_sorted[i].phi);
      if (dphi>M_PI) dphi = twopi-dphi;
      if (dphi<EPSILON_COLLINEAR){
	// i is collinear with j; add the momentum (i) to the momentum (j) 
	p_sorted[j] += p_sorted[i];
	// set collinearity test to true
	collinear = true;
      }
      j++;
    }
    // if no collinearity detected, add the particle to our list
    if (!collinear)
      p_uncol_hard.push_back(p_sorted[i]);
    i++;
  }

  return 0;
}


// add a list of protocones
//  - protocones  list of protocones (initial jet candidates)
//  - R2          cone radius (squared)
//  - ptmin       minimal pT allowed for jets
//-------------------------------------------------------------
int Csplit_merge::add_protocones(vector<Cmomentum> *protocones, double R2, double ptmin){
  int i;
  Cmomentum *c;
  Cmomentum *v;
  double eta, phi;
  double dx, dy;
  double R;
  Cjet jet;

  if (protocones->size()==0)
    return 1;

  pt_min2 = ptmin*ptmin;
  R = sqrt(R2);

  // browse protocones
  // for each of them, build the list of particles in them
  for (vector<Cmomentum>::iterator p_it = protocones->begin();p_it != protocones->end();p_it++){
    // initialise variables
    c = &(*p_it);

    // note: cones have been tested => their (eta,phi) coordinates are computed
    eta = c->eta;
    phi = c->phi;

    // browse particles to create cone contents
    // note that jet is always initialised with default values at this level
    jet.v = Cmomentum();
    jet.pt_tilde=0;
    jet.contents.clear();
    for (i=0;i<n_left;i++){
      v = &(p_remain[i]);
      // for security, avoid including particles with infinite rapidity)
      // NO NEEDED ANYMORE SINCE REMOVED FROM p_remain FROM THE BEGINNING
      //if (fabs(v->pz)!=v->E){
      dx = eta - v->eta;
      dy = fabs(phi - v->phi);
      if (dy>M_PI) 
	dy -= twopi;
      if (dx*dx+dy*dy<R2){
	jet.contents.push_back(v->parent_index);
	jet.v+= *v;
	jet.pt_tilde+= pt[v->parent_index];
	v->index=0;
      }
    }
    jet.n=jet.contents.size();

    // set the momentum in protocones 
    // (it was only known through eta and phi up to now)
    *c = jet.v;
    c->eta = eta; // restore exact original coords
    c->phi = phi; // to avoid rounding error inconsistencies

    // set the jet range
    jet.range=Ceta_phi_range(eta,phi,R);

#ifdef DEBUG_SPLIT_MERGE
    cout << "adding jet: ";
    for (int i2=0;i2<jet.n;i2++)
      cout << jet.contents[i2] << " ";
    cout << endl;
#endif

    // add it to the list of jets
    insert(jet);
  }
  
  // update list of included particles
  n_pass++;

#ifdef DEBUG_SPLIT_MERGE
  cout << "remaining particles: "; 
#endif
  int j=0;
  for (i=0;i<n_left;i++){
    if (p_remain[i].index){
      // copy particle
      p_remain[j]=p_remain[i];
      p_remain[j].parent_index = p_remain[i].parent_index;
      p_remain[j].index=1;
      // set run in initial list
      particles[p_remain[j].parent_index].index = n_pass;
#ifdef DEBUG_SPLIT_MERGE
      cout << p_remain[j].parent_index << " ";
#endif
      j++;
    }
  }
#ifdef DEBUG_SPLIT_MERGE
  cout << endl;
#endif
  n_left = j;
  p_remain.resize(j);

  merge_collinear_and_remove_soft();

  return 0;
}


/*
 * remove the hardest protocone and declare it a a jet 
 *  - protocones  list of protocones (initial jet candidates)
 *  - R2          cone radius (squared)
 *  - ptmin       minimal pT allowed for jets
 * return 0 on success, 1 on error
 *
 * The list of remaining particles (and the uncollinear-hard ones)
 * is updated.
 */
int Csplit_merge::add_hardest_protocone_to_jets(std::vector<Cmomentum> *protocones, double R2, double ptmin){

  int i;
  Cmomentum *c;
  Cmomentum *v;
  double eta, phi;
  double dx, dy;
  double R;
  Cjet jet, jet_candidate;
  bool found_jet = false;

  if (protocones->size()==0)
    return 1;

  pt_min2 = ptmin*ptmin;
  R = sqrt(R2);

  // browse protocones
  // for each of them, build the list of particles in them
  for (vector<Cmomentum>::iterator p_it = protocones->begin();p_it != protocones->end();p_it++){
    // initialise variables
    c = &(*p_it);

    // note: cones have been tested => their (eta,phi) coordinates are computed
    eta = c->eta;
    phi = c->phi;

    // NOTE: this is common to this method and add_protocones, so it
    // could be moved into a 'build_jet_from_protocone' method
    //
    // browse particles to create cone contents
    jet_candidate.v = Cmomentum();
    jet_candidate.pt_tilde=0;
    jet_candidate.contents.clear();
    for (i=0;i<n_left;i++){
      v = &(p_remain[i]);
      // for security, avoid including particles with infinite rapidity)
      // NO NEEDED ANYMORE SINCE REMOVED FROM p_remain FROM THE BEGINNING
      //if (fabs(v->pz)!=v->E){
      dx = eta - v->eta;
      dy = fabs(phi - v->phi);
      if (dy>M_PI) 
	dy -= twopi;
      if (dx*dx+dy*dy<R2){
	jet_candidate.contents.push_back(v->parent_index);
	jet_candidate.v+= *v;
	jet_candidate.pt_tilde+= pt[v->parent_index];
	v->index=0;
      }
    }
    jet_candidate.n=jet_candidate.contents.size();

    // set the momentum in protocones 
    // (it was only known through eta and phi up to now)
    *c = jet_candidate.v;
    c->eta = eta; // restore exact original coords
    c->phi = phi; // to avoid rounding error inconsistencies

    // set the jet range
    jet_candidate.range=Ceta_phi_range(eta,phi,R);

    // check that the protojet has large enough pt
    if (jet_candidate.v.perp2()<pt_min2)
      continue;

    // assign the split-merge (or progressive-removal) squared scale variable
    if (_user_scale) {
      // sm_var2 is the signed square of the user scale returned
      // for the jet candidate
      jet_candidate.sm_var2 = (*_user_scale)(jet_candidate);
      jet_candidate.sm_var2 *= abs(jet_candidate.sm_var2);
    } else {
      jet_candidate.sm_var2 = get_sm_var2(jet_candidate.v, jet_candidate.pt_tilde);
    }

    // now check if it is possibly the hardest
    if ((! found_jet) ||
	(_user_scale ? _user_scale->is_larger(jet_candidate, jet)
	             : ptcomparison(jet_candidate, jet))){
      jet = jet_candidate;
      found_jet = true;
    }
  }

  // make sure at least one of the jets has passed the selection
  if (!found_jet) return 1;  
  
  // add the jet to the list of jets
  jets.push_back(jet);
  jets[jets.size()-1].v.build_etaphi();

#ifdef DEBUG_SPLIT_MERGE
  cout << "PR-Jet " << jets.size() << " [size " << next_jet.contents.size() << "]:";
#endif
    
  // update the list of what particles are left
  int p_remain_index = 0;
  int contents_index = 0;
  //sort(next_jet.contents.begin(),next_jet.contents.end());
  for (int index=0;index<n_left;index++){
    if ((contents_index<(int) jet.contents.size()) &&
	(p_remain[index].parent_index == jet.contents[contents_index])){
      // this particle belongs to the newly found jet
      // set pass in initial list
      particles[p_remain[index].parent_index].index = n_pass;
#ifdef DEBUG_SPLIT_MERGE
      cout << " " << jet.contents[contents_index];
#endif
      contents_index++;
    } else {
      // this particle still has to be clustered
      p_remain[p_remain_index] = p_remain[index];
      p_remain[p_remain_index].parent_index = p_remain[index].parent_index;
      p_remain[p_remain_index].index=1;
      p_remain_index++;
    }
  }
  p_remain.resize(n_left-jet.contents.size());
  n_left = p_remain.size();
  jets[jets.size()-1].pass = particles[jet.contents[0]].index;

  // increase the pass index
  n_pass++;

#ifdef DEBUG_SPLIT_MERGE
  cout << endl;
#endif

  // male sure the list of uncol_hard particles (used for the next
  // stable cone finding) is updated [could probably be more
  // efficient]
  merge_collinear_and_remove_soft();
  
  return 0;
}

/*
 * really do the splitting and merging
 * At the end, the vector jets is filled with the jets found.
 * the 'contents' field of each jets contains the indices
 * of the particles included in that jet. 
 *  - overlap_tshold    threshold for splitting/merging transition
 *  - ptmin             minimal pT allowed for jets
 * return the number of jets is returned
 ******************************************************************/
int Csplit_merge::perform(double overlap_tshold, double ptmin){
  // iterators for the 2 jets
  cjet_iterator j1;
  cjet_iterator j2;

  pt_min2 = ptmin*ptmin;

  if (candidates->size()==0)
    return 0;

  if (overlap_tshold>=1.0 || overlap_tshold <= 0) {
    ostringstream message;
    message << "Illegal value for overlap_tshold, f = " << overlap_tshold;
    message << "  (legal values are 0<f<1)";
    throw Csiscone_error(message.str());
  }

  // overlap (the contents of this variable depends on the choice for
  // the split--merge variable.)
  // Note that the square of the ovelap is used
  double overlap2;

  // avoid to compute tshold*tshold at each overlap
  double overlap_tshold2 = overlap_tshold*overlap_tshold;

  do{
    if (candidates->size()>0){
      // browse for the first jet
      j1 = candidates->begin();
      
      // if hardest jet does not pass threshold then nothing else will
      // either so one stops the split merge.
      if (j1->sm_var2<SM_var2_hardest_cut_off) {break;}

      // browse for the second jet
      j2 = j1;
      j2++;
      int j2_relindex = 1; // used only in ifdef, but costs little so keep it outside

      while (j2 != candidates->end()){
#ifdef DEBUG_SPLIT_MERGE
	show();
#endif
	// check overlapping
	if (get_overlap(*j1, *j2, &overlap2)){
	  // check if overlapping energy passes threshold
	  // Note that this depends on the ordering variable
#ifdef DEBUG_SPLIT_MERGE
          cout << "overlap between cdt 1 and cdt " << j2_relindex+1 << " with overlap " 
               << sqrt(overlap2/j2->sm_var2) << endl<<endl;
#endif
	  if (overlap2<overlap_tshold2*j2->sm_var2){
	    // split jets
	    split(j1, j2);
	    
	    // update iterators
	    j2 = j1 = candidates->begin();
            j2_relindex = 0;
	  } else {
	    // merge jets
	    merge(j1, j2);
	    
	    // update iterators
	    j2 = j1 = candidates->begin();
            j2_relindex = 0;
	  }
	}
        // watch out: split/merge might have caused new jets with pt <
        // ptmin to disappear, so the total number of jets may
        // have changed by more than expected and j2 might already by
        // the end of the candidates list...
        j2_relindex++;
	if (j2 != candidates->end()) j2++;
      } // end of loop on the second jet
      
      if (j1 != candidates->end()) {
        // all "second jet" passed without overlapping
        // (otherwise we won't leave the j2 loop)
        // => store jet 1 as real jet
        jets.push_back(*j1);
        jets[jets.size()-1].v.build_etaphi();
        // a bug where the contents has non-zero size has been cropping
        // up in many contexts -- so catch it!
        assert(j1->contents.size() > 0);
        jets[jets.size()-1].pass = particles[j1->contents[0]].index;
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
        cand_refs.erase(j1->v.ref);
#endif
        candidates->erase(j1);

	//// test that the hardest jet pass the potential cut-off
	//if ((candidates->size()!=0) && 
	//    (candidates->begin()->sm_var2<SM_var2_hardest_cut_off)){
	//  candidates->clear();
	//}
      }
    }
  } while (candidates->size()>0);

  // sort jets by pT
  sort(jets.begin(), jets.end(), jets_pt_less);
#ifdef DEBUG_SPLIT_MERGE
      show();
#endif

  return jets.size();
}



// save the event on disk
//  - flux   stream used to save jet contents
//--------------------------------------------
int Csplit_merge::save_contents(FILE *flux){
  jet_iterator it_j;
  Cjet *j1;
  int i1, i2;

  fprintf(flux, "# %d jets found\n", (int) jets.size());
  fprintf(flux, "# columns are: eta, phi, pt and number of particles for each jet\n");
  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j1 = &(*it_j);
    j1->v.build_etaphi();
    fprintf(flux, "%f\t%f\t%e\t%d\n", 
	    j1->v.eta, j1->v.phi, j1->v.perp(), j1->n);
  }
  
  fprintf(flux, "# jet contents\n");
  fprintf(flux, "# columns are: eta, phi, pt, particle index and jet number\n");
  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j1 = &(*it_j);
    for (i2=0;i2<j1->n;i2++)
      fprintf(flux, "%f\t%f\t%e\t%d\t%d\n", 
      	      particles[j1->contents[i2]].eta, particles[j1->contents[i2]].phi,
      	      particles[j1->contents[i2]].perp(), j1->contents[i2], i1);
  }
  
  return 0;
}


// show current jets/candidate status
//------------------------------------
int Csplit_merge::show(){
  jet_iterator it_j;
  cjet_iterator it_c;
  Cjet *j;
  const Cjet *c;
  int i1, i2;

  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j = &(*it_j);
    fprintf(stdout, "jet %2d: %e\t%e\t%e\t%e\t", i1+1,
	    j->v.px, j->v.py, j->v.pz, j->v.E);
    for (i2=0;i2<j->n;i2++)
      fprintf(stdout, "%d ", j->contents[i2]);
    fprintf(stdout, "\n");
  }
  
  for (it_c = candidates->begin(), i1=0 ; it_c != candidates->end() ; it_c++, i1++){
    c = &(*it_c);
    fprintf(stdout, "cdt %2d: %e\t%e\t%e\t%e\t%e\t", i1+1,
	    c->v.px, c->v.py, c->v.pz, c->v.E, sqrt(c->sm_var2));
    for (i2=0;i2<c->n;i2++)
      fprintf(stdout, "%d ", c->contents[i2]);
    fprintf(stdout, "\n");
  }
  
  fprintf(stdout, "\n");
  return 0;
}


// get the overlap between 2 jets
//  - j1        first jet
//  - j2        second jet
//  - overlap2  returned overlap^2 (determined by the choice of SM variable)
// return true if overlapping, false if disjoint
//---------------------------------------------------------------------
bool Csplit_merge::get_overlap(const Cjet &j1, const Cjet &j2, double *overlap2){
  // check if ranges overlap
  if (!is_range_overlap(j1.range,j2.range))
    return false;

  int i1,i2;
  bool is_overlap;

  // initialise
  i1=i2=idx_size=0;
  is_overlap = false;
  Cmomentum v;
  double pt_tilde=0.0;

  // compute overlap
  // at the same time, we store union in indices
  do{
    if (j1.contents[i1]<j2.contents[i2]){
      indices[idx_size] = j1.contents[i1];
      i1++;
    } else if (j1.contents[i1]>j2.contents[i2]){
      indices[idx_size] = j2.contents[i2];
      i2++;
    } else { // (j1.contents[i1]==j2.contents[i2])
      v += particles[j1.contents[i1]];
      pt_tilde += pt[j1.contents[i1]];
      indices[idx_size] = j1.contents[i1];
      i1++;
      i2++;
      is_overlap = true;
    }
    idx_size++;
  } while ((i1<j1.n) && (i2<j2.n));

  // finish computing union
  // (only needed if overlap !)
  if (is_overlap){
    while (i1<j1.n){
      indices[idx_size] = j1.contents[i1];
      i1++;
      idx_size++;
    }
    while (i2<j2.n){
      indices[idx_size] = j2.contents[i2];
      i2++;
      idx_size++;
    }
  }

  // assign the overlapping var as return variable
  (*overlap2) = get_sm_var2(v, pt_tilde);

  return is_overlap;
}



// split the two given jet.
// during this procedure, the jets j1 & j2 are replaced
// by 2 new jets. Common particles are associted to the 
// closest initial jet.
//  - it_j1  iterator of the first jet in 'candidates'
//  - it_j2  iterator of the second jet in 'candidates'
//  - j1     first jet (Cjet instance)
//  - j2     second jet (Cjet instance)
// return true on success, false on error
////////////////////////////////////////////////////////////////
bool Csplit_merge::split(cjet_iterator &it_j1, cjet_iterator &it_j2){
  int i1, i2;
  Cjet jet1, jet2;
  double eta1, phi1, pt1_weight, eta2, phi2, pt2_weight;
  double dx1, dy1, dx2, dy2;
  Cmomentum tmp;
  Cmomentum *v;

  // shorthand to avoid having "->" all over the place
  const Cjet & j1 = * it_j1;
  const Cjet & j2 = * it_j2;

  i1=i2=0;
  jet2.v = jet1.v = Cmomentum();
  jet2.pt_tilde = jet1.pt_tilde = 0.0;

  // compute centroids
  // When use_pt_weighted_splitting is activated, the
  // "geometrical" distance is weighted by the inverse
  // of the pt of the protojet
  // This is stored in pt{1,2}_weight
  tmp = j1.v;
  tmp.build_etaphi();
  eta1 = tmp.eta;
  phi1 = tmp.phi;
  pt1_weight = (use_pt_weighted_splitting) ? 1.0/tmp.perp2() : 1.0;

  tmp = j2.v;
  tmp.build_etaphi();
  eta2 = tmp.eta;
  phi2 = tmp.phi;
  pt2_weight = (use_pt_weighted_splitting) ? 1.0/tmp.perp2() : 1.0;

  jet1.v = jet2.v = Cmomentum();  

  // compute jet splitting
  do{
    if (j1.contents[i1]<j2.contents[i2]){
      // particle i1 belong only to jet 1
      v = &(particles[j1.contents[i1]]);
      jet1.contents.push_back(j1.contents[i1]);
      jet1.v += *v;
      jet1.pt_tilde += pt[j1.contents[i1]];
      i1++;
      jet1.range.add_particle(v->eta,v->phi);
    } else if (j1.contents[i1]>j2.contents[i2]){
      // particle i2 belong only to jet 2
      v = &(particles[j2.contents[i2]]);
      jet2.contents.push_back(j2.contents[i2]);
      jet2.v += *v;
      jet2.pt_tilde += pt[j2.contents[i2]];
      i2++;
      jet2.range.add_particle(v->eta,v->phi);
    } else { // (j1.contents[i1]==j2.contents[i2])
      // common particle, decide which is the closest centre
      v = &(particles[j1.contents[i1]]);

      // distance w.r.t. centroid 1
      dx1 = eta1 - v->eta;
      dy1 = fabs(phi1 - v->phi);
      if (dy1>M_PI) 
	dy1 -= twopi;

      // distance w.r.t. centroid 2
      dx2 = eta2 - v->eta;
      dy2 = fabs(phi2 - v->phi);
      if (dy2>M_PI) 
	dy2 -= twopi;

      //? what when == ? 
      // When use_pt_weighted_splitting is activated, the
      // "geometrical" distance is weighted by the inverse
      // of the pt of the protojet
      double d1sq = (dx1*dx1+dy1*dy1)*pt1_weight;
      double d2sq = (dx2*dx2+dy2*dy2)*pt2_weight;
      // do bookkeeping on most ambiguous split
      if (fabs(d1sq-d2sq) < most_ambiguous_split) 
        most_ambiguous_split = fabs(d1sq-d2sq);

      if (d1sq<d2sq){
	// particle i1 belong only to jet 1
	jet1.contents.push_back(j1.contents[i1]);
	jet1.v += *v;
	jet1.pt_tilde += pt[j1.contents[i1]];
	jet1.range.add_particle(v->eta,v->phi);
      } else {
	// particle i2 belong only to jet 2
	jet2.contents.push_back(j2.contents[i2]);
	jet2.v += *v;
	jet2.pt_tilde += pt[j2.contents[i2]];
	jet2.range.add_particle(v->eta,v->phi);
      }      

      i1++;
      i2++;
    }
  } while ((i1<j1.n) && (i2<j2.n));
  
  while (i1<j1.n){
    v = &(particles[j1.contents[i1]]);
    jet1.contents.push_back(j1.contents[i1]);
    jet1.v += *v;
    jet1.pt_tilde += pt[j1.contents[i1]];
    i1++;
    jet1.range.add_particle(v->eta,v->phi);
  }
  while (i2<j2.n){
    v = &(particles[j2.contents[i2]]);
    jet2.contents.push_back(j2.contents[i2]);
    jet2.v += *v;
    jet2.pt_tilde += pt[j2.contents[i2]];
    i2++;
    jet2.range.add_particle(v->eta,v->phi);
  }

  // finalise jets
  jet1.n = jet1.contents.size();
  jet2.n = jet2.contents.size();

  //jet1.range = j1.range;
  //jet2.range = j2.range;

  // remove previous jets
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  cand_refs.erase(j1.v.ref);
  cand_refs.erase(j2.v.ref);
#endif
  candidates->erase(it_j1);
  candidates->erase(it_j2);

  // reinsert new ones
  insert(jet1);
  insert(jet2);

  return true;
}

// merge the two given jet.
// during this procedure, the jets j1 & j2 are replaced
// by 1 single jets containing both of them.
//  - it_j1  iterator of the first jet in 'candidates'
//  - it_j2  iterator of the second jet in 'candidates'
// return true on success, false on error
////////////////////////////////////////////////////////////////
bool Csplit_merge::merge(cjet_iterator &it_j1, cjet_iterator &it_j2){
  Cjet jet;
  int i;

  // build new jet
  // note: particles within j1 & j2 have already been stored in indices
  for (i=0;i<idx_size;i++){
    jet.contents.push_back(indices[i]);
    jet.v += particles[indices[i]];
    jet.pt_tilde += pt[indices[i]];
  }
  jet.n = jet.contents.size();

  // deal with ranges
  jet.range = range_union(it_j1->range, it_j2->range);

  // remove old candidates
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if (merge_identical_protocones){
    cand_refs.erase(it_j1->v.ref);
    cand_refs.erase(it_j2->v.ref);
  }
#endif
  candidates->erase(it_j1);
  candidates->erase(it_j2);

  // reinsert new candidate
  insert(jet);

  return true;
}

/**
 * Check whether or not a jet has to be inserted in the 
 * list of protojets. If it has, set its sm_variable and
 * insert it to the list of protojets.
 */
bool Csplit_merge::insert(Cjet &jet){

  // eventually check that no other candidate are present with the
  // same cone contents. We recall that this automatic merging of
  // identical protocones can lead to infrared-unsafe situations.
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if ((merge_identical_protocones) && (!cand_refs.insert(jet.v.ref).second))
    return false;
#endif

  // check that the protojet has large enough pt
  if (jet.v.perp2()<pt_min2)
    return false;

  // assign SM variable
  jet.sm_var2 = get_sm_var2(jet.v, jet.pt_tilde);

  // insert the jet.
  candidates->insert(jet);

  return true;
}

/**
 * given a 4-momentum and its associated pT, return the 
 * variable that has to be used for SM
 * \param v          4 momentum of the protojet
 * \param pt_tilde   pt_tilde of the protojet
 */
double Csplit_merge::get_sm_var2(Cmomentum &v, double &pt_tilde){
  switch(ptcomparison.split_merge_scale) {
  case SM_pt:      return v.perp2();            
  case SM_mt:      return v.perpmass2();        
  case SM_pttilde: return pt_tilde*pt_tilde;
  case SM_Et:      return v.Et2();
  default:
    throw Csiscone_error("Unsupported split-merge scale choice: "
                                 + ptcomparison.SM_scale_name());
  }

  return 0.0;
}

}
