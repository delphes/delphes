//FJSTARTHEADER
// $Id: ClusterSequence_DumbN3.cc 3433 2014-07-23 08:17:03Z salam $
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
#include<iostream>
#include<cmath>
#include <cstdlib>
#include<cassert>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


using namespace std;


//----------------------------------------------------------------------
/// Run the clustering in a very slow variant of the N^3 algorithm. 
///
/// The only thing this routine has going for it is that memory usage
/// is O(N)!
void ClusterSequence::_really_dumb_cluster () {

  // the array that will be overwritten here will be one
  // of pointers to jets.
  vector<PseudoJet *> jetsp(_jets.size());
  vector<int>         indices(_jets.size());

  for (size_t i = 0; i<_jets.size(); i++) {
    jetsp[i] = & _jets[i];
    indices[i] = i;
  }

  for (int n = jetsp.size(); n > 0; n--) {
    int ii, jj;
    // find smallest beam distance [remember jet_scale_for_algorithm 
    // will return kt^2 for the kt-algorithm and 1 for the Cambridge/Aachen]
    double ymin = jet_scale_for_algorithm(*(jetsp[0]));
    ii = 0; jj = -2;
    for (int i = 0; i < n; i++) {
      double yiB = jet_scale_for_algorithm(*(jetsp[i]));
      if (yiB < ymin) {
	ymin = yiB; ii = i; jj = -2;}
    }

    // find smallest distance between pair of jetsp
    for (int i = 0; i < n-1; i++) {
      for (int j = i+1; j < n; j++) {
	//double y = jetsp[i]->kt_distance(*jetsp[j])*_invR2;
	double y = min(jet_scale_for_algorithm(*(jetsp[i])), 
		       jet_scale_for_algorithm(*(jetsp[j])))
	            * jetsp[i]->plain_distance(*jetsp[j])*_invR2;
	if (y < ymin) {ymin = y; ii = i; jj = j;}
      }
    }

    // output recombination sequence
    // old "ktclus" way of labelling
    //cout <<n<< " "<< ii+1 << " with " << jj+1 << "; y = "<< ymin<<endl;
    //OBS // new delaunay way of labelling
    //OBS int jjindex_or_beam, iiindex;
    //OBS if (jj < 0) {jjindex_or_beam = BeamJet; iiindex = indices[ii];} 
    //OBS else {
    //OBS   jjindex_or_beam = max(indices[ii],indices[jj]);
    //OBS   iiindex =         min(indices[ii],indices[jj]);
    //OBS }

    // now recombine
    int newn = 2*jetsp.size() - n;
    if (jj >= 0) {
      // combine pair
      int nn; // new jet index
      _do_ij_recombination_step(jetsp[ii]-&_jets[0], 
				jetsp[jj]-&_jets[0], ymin, nn);
      
      // sort out internal bookkeeping
      jetsp[ii] = &_jets[nn];
      // have jj point to jet that was pointed at by n-1 
      // (since original jj is no longer current, so put n-1 into jj)
      jetsp[jj] = jetsp[n-1];
      indices[ii] = newn;
      indices[jj] = indices[n-1];

      //OBS_jets.push_back(*jetsp[ii] + *jetsp[jj]);
      //OBSjetsp[ii] = &_jets[_jets.size()-1];
      //OBS// have jj point to jet that was pointed at by n-1 
      //OBS// (since original jj is no longer current, so put n-1 into jj)
      //OBSjetsp[jj] = jetsp[n-1];
      //OBS
      //OBSindices[ii] = newn;
      //OBSindices[jj] = indices[n-1];
      //OBS_add_step_to_history(newn,iiindex,
      //OBS			      jjindex_or_beam,_jets.size()-1,ymin);
    } else {
      // combine ii with beam
      _do_iB_recombination_step(jetsp[ii]-&_jets[0], ymin);
      // put last jet (pointer) in place of ii (which has disappeared)
      jetsp[ii] = jetsp[n-1];
      indices[ii] = indices[n-1];
      //OBS_add_step_to_history(newn,iiindex,jjindex_or_beam,Invalid, ymin);
    }
  }

}

FASTJET_END_NAMESPACE

