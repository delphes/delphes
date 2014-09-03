//FJSTARTHEADER
// $Id: CDFMidPointPlugin.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Error.hh"
#include <sstream>

// CDF stuff
#include "MidPointAlgorithm.hh"
#include "PhysicsTower.hh"
#include "Cluster.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace cdf;

bool CDFMidPointPlugin::_first_time = true;

string CDFMidPointPlugin::description () const {
  ostringstream desc;
  
  string sm_scale_string = "split-merge uses ";
  switch(_sm_scale) {
  case SM_pt:
    sm_scale_string += "pt";
    break;
  case SM_Et:
    sm_scale_string += "Et";
    break;
  case SM_mt:
    sm_scale_string += "mt";
    break;
  case SM_pttilde:
    sm_scale_string += "pttilde (scalar sum of pts)";
    break;
  default:
    ostringstream err;
    err << "Unrecognized split-merge scale choice = " << _sm_scale;
    throw Error(err.str());
  }

  
  if (cone_area_fraction() == 1) {
    desc << "CDF MidPoint jet algorithm, with " ;
  } else {
    desc << "CDF MidPoint+Searchcone jet algorithm, with ";
  }
  desc << "seed_threshold = "     << seed_threshold     () << ", "
       << "cone_radius = "        << cone_radius        () << ", "
       << "cone_area_fraction = " << cone_area_fraction () << ", " 
       << "max_pair_size = "      << max_pair_size      () << ", "
       << "max_iterations = "     << max_iterations     () << ", "
       << "overlap_threshold  = " << overlap_threshold  () << ", "
       << sm_scale_string ;

  return desc.str();
}


void CDFMidPointPlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());
 
  // create the physics towers needed by the CDF code
  vector<PhysicsTower> towers;
  towers.reserve(clust_seq.jets().size());
  for (unsigned i = 0; i < clust_seq.jets().size(); i++) {
    LorentzVector fourvect(clust_seq.jets()[i].px(),
			   clust_seq.jets()[i].py(),
			   clust_seq.jets()[i].pz(),
			   clust_seq.jets()[i].E());
    PhysicsTower tower(fourvect);
    // misuse one of the indices for tracking, since the MidPoint
    // implementation doesn't seem to make use of these indices
    tower.calTower.iEta = i;
    towers.push_back(tower);
  }

  // prepare the CDF algorithm
  MidPointAlgorithm m(_seed_threshold,_cone_radius,_cone_area_fraction,
		      _max_pair_size,_max_iterations,_overlap_threshold,
                      MidPointAlgorithm::SplitMergeScale(_sm_scale));
    
  // run the CDF algorithm
  std::vector<Cluster> jets;
  m.run(towers,jets);


  // now transfer the jets back into our own structure -- we will
  // mimic the cone code with a sequential recombination sequence in
  // which the jets are built up by adding one particle at a time
  for(vector<Cluster>::const_iterator jetIter = jets.begin(); 
                                      jetIter != jets.end(); jetIter++) {
    const vector<PhysicsTower> & tower_list = jetIter->towerList;
    int jet_k = tower_list[0].calTower.iEta;
  
    int ntow = int(jetIter->towerList.size());
    for (int itow = 1; itow < ntow; itow++) {
      int jet_i = jet_k;
      // retrieve our misappropriated index for the jet
      int jet_j = tower_list[itow].calTower.iEta;
      // do a fake recombination step with dij=0
      double dij = 0.0;
      clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
    }
  
    // NB: put a sensible looking d_iB just to be nice...
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);
  }


  // following code is for testing only
  //cout << endl;
  //for(vector<Cluster>::const_iterator jetIter = jets.begin(); 
  //                                    jetIter != jets.end(); jetIter++) {
  //  cout << jetIter->fourVector.pt() << " " << jetIter->fourVector.y() << endl;
  //}
  //cout << "-----------------------------------------------------\n";
  //vector<PseudoJet> ourjets(clust_seq.inclusive_jets());
  //for (vector<PseudoJet>::const_reverse_iterator ourjet = ourjets.rbegin();
  //     ourjet != ourjets.rend(); ourjet++) {
  //  cout << ourjet->perp() << " " << ourjet->rap() << endl;
  //}
  //cout << endl;
}

// print a banner for reference to the 3rd-party code
void CDFMidPointPlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the CDF MidPoint plugin for FastJet                     " << endl;
  (*ostr) << "# This is based on an implementation provided by Joey Huston.             " << endl;
  (*ostr) << "# If you use this plugin, please cite                                     " << endl;
  (*ostr) << "#   G. C. Blazey et al., hep-ex/0005012.                                  " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                             " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
