//FJSTARTHEADER
// $Id: CDFJetCluPlugin.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include <sstream>
#include <cassert>

// CDF stuff
#include "JetCluAlgorithm.hh"
#include "PhysicsTower.hh"
#include "Cluster.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace cdf;

bool CDFJetCluPlugin::_first_time = true;

string CDFJetCluPlugin::description () const {
  ostringstream desc;
  
  desc << "CDF JetClu jet algorithm with " 
       << "seed_threshold = "     << seed_threshold    () << ", "
       << "cone_radius = "        << cone_radius       () << ", "
       << "adjacency_cut = "      << adjacency_cut     () << ", " 
       << "max_iterations = "     << max_iterations    () << ", "
       << "iratch = "             << iratch            () << ", "
       << "overlap_threshold = "  << overlap_threshold () ;

  return desc.str();
}


void CDFJetCluPlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());
 
  // create the physics towers needed by the CDF code
  vector<PhysicsTower> towers;
  towers.reserve(clust_seq.jets().size());

  // create a map to identify jets (actually just the input particles)...
  //map<double,int> jetmap;

  for (unsigned i = 0; i < clust_seq.jets().size(); i++) {
    PseudoJet particle(clust_seq.jets()[i]);
    //_insert_unique(particle, jetmap);
    LorentzVector fourvect(particle.px(), particle.py(),
			   particle.pz(), particle.E());
    PhysicsTower tower(fourvect);
    // add tracking information for later
    tower.fjindex = i;
    towers.push_back(tower);
  }

  // prepare the CDF algorithm
  JetCluAlgorithm j(seed_threshold(), cone_radius(), adjacency_cut(),
                    max_iterations(), iratch(), overlap_threshold());
    
  // run the CDF algorithm
  std::vector<Cluster> jets;
  j.run(towers,jets);


  // now transfer the jets back into our own structure -- we will
  // mimic the cone code with a sequential recombination sequence in
  // which the jets are built up by adding one particle at a time

  // NB: with g++-4.0, the reverse iterator code gave problems, so switch
  //     to indices instead
  //for(vector<Cluster>::const_reverse_iterator jetIter = jets.rbegin(); 
  //                                    jetIter != jets.rend(); jetIter++) {
  //  const vector<PhysicsTower> & tower_list = jetIter->towerList;
  //  int jet_k = jetmap[tower_list[0].fourVector.E];
  //
  //  int ntow = int(jetIter->towerList.size());

  for(int iCDFjets = jets.size()-1; iCDFjets >= 0; iCDFjets--) {

    const vector<PhysicsTower> & tower_list = jets[iCDFjets].towerList;
    int ntow = int(tower_list.size());
    
    // 2008-09-04: sort the towerList (according to fjindex) so as
    //             to have a consistent order for particles in jet
    //             (necessary because addition of ultra-soft particles
    //             sometimes often modifies the order, while maintaining
    //             the same overall set)
    vector<int>    jc_indices(ntow);
    vector<double> fj_indices(ntow); // use double: benefit from existing routine
    for (int itow = 0; itow < ntow; itow++) {
      jc_indices[itow] = itow;
      fj_indices[itow] = tower_list[itow].fjindex;
    }
    sort_indices(jc_indices, fj_indices);

    int jet_k = tower_list[jc_indices[0]].fjindex;
  
    for (int itow = 1; itow < ntow; itow++) {
      if (tower_list[jc_indices[itow]].Et() > 1e-50) {
      }
      int jet_i = jet_k;
      // retrieve our index for the jet
      int jet_j;
      jet_j = tower_list[jc_indices[itow]].fjindex;

      // safety check
      assert (jet_j >= 0 && jet_j < int(towers.size()));

      // do a fake recombination step with dij=0
      double dij = 0.0;

      // JetClu does E-scheme recombination so we can stick with the
      // simple option
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
  //for (vector<PseudoJet>::const_iterator ourjet = ourjets.begin();
  //     ourjet != ourjets.end(); ourjet++) {
  //  cout << ourjet->perp() << " " << ourjet->rap() << endl;
  //}
  //cout << endl;
}


// print a banner for reference to the 3rd-party code
void CDFJetCluPlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the CDF JetClu plugin for FastJet                       " << endl;
  (*ostr) << "# This is based on an implementation provided by Joey Huston.             " << endl;
  (*ostr) << "# If you use this plugin, please cite                                     " << endl;
  (*ostr) << "#   F. Abe et al. [CDF Collaboration], Phys. Rev. D 45 (1992) 1448.       " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                             " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
