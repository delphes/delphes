//FJSTARTHEADER
// $Id: ATLASConePlugin.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2007-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ATLASConePlugin.hh"

// SpartyJet stuff
#include "CommonUtils.hh"
#include "JetConeFinderTool.hh"
#include "JetSplitMergeTool.hh"

// other stuff
#include <vector>
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

bool ATLASConePlugin::_first_time = true;

string ATLASConePlugin::description () const {
  ostringstream desc;
  desc << "ATLASCone plugin with R = "<< _radius 
       << ", seed threshold = " << _seedPt
       << ", overlap threshold f = " << _f;
  return desc.str();
}

void ATLASConePlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());


  // transfer the list of PseudoJet into a atlas::Jet::jet_list_t
  //  jet_list_t is a vector<Jet*>
  // We set the index of the 4-vect to trace the constituents at the end
  //------------------------------------------------------------------
  // cout << "ATLASConePlugin: transferring vectors from ClusterSequence" << endl;
  atlas::JetConeFinderTool::jetcollection_t jets_ptr;
  vector<atlas::Jet*> particles_ptr;

  for (unsigned int i=0 ; i<clust_seq.jets().size() ; i++) {
    const PseudoJet & mom = clust_seq.jets()[i];
    
    // first create the particle
    atlas::Jet *particle = new atlas::Jet(mom.px(), mom.py(), mom.pz(), mom.E(), i);
    particles_ptr.push_back(particle);

    // then add it to the list of particles we'll use for the clustering
    atlas::Jet *jet = new atlas::Jet;
    jet->set_index(particle->index());
    jet->addConstituent(particle);

    // and finally add that one to the list of jets
    jets_ptr.push_back(jet);
  }
  // cout << "ATLASCone: " << jets_ptr.size() << " particles to cluster" << endl;

  // search the stable cones
  //------------------------------------------------------------------
  // cout << "ATLASConePlugin: searching for stable cones" << endl;
  atlas::JetConeFinderTool stable_cone_finder;

  // set the parameters
  stable_cone_finder.m_coneR  = _radius;
  stable_cone_finder.m_seedPt = _seedPt;

  // really do the search.
  // Note that the list of protocones is returned 
  // through the argument
  stable_cone_finder.execute(jets_ptr);
  // cout << "ATLASCone: " << jets_ptr.size() << " stable cones found" << endl;

  // perform the split-merge
  //------------------------------------------------------------------
  // cout << "ATLASConePlugin: running the split-merge" << endl;
  atlas::JetSplitMergeTool split_merge;
  
  // set the parameters
  split_merge.m_f = _f;

  // do the work
  // again, the list of jets is returned through the argument
  split_merge.execute(&jets_ptr);
  // cout << "ATLASCone: " << jets_ptr.size() << " jets after split--merge" << endl;

  // build the FastJet jets (a la SISConePlugin)
  //------------------------------------------------------------------
  // cout << "ATLASConePlugin: backporting jets to the ClusterSequence" << endl;
  for (atlas::Jet::jet_list_t::iterator jet_it = jets_ptr.begin() ;
	 jet_it != jets_ptr.end(); jet_it++){
    // iterate over the constituents, starting from the first one
    // that we just take as a reference
    atlas::Jet::constit_vect_t::iterator constit_it = (*jet_it)->firstConstituent();
    // cout << " atlas: jet has " << (*jet_it)->getConstituentNum() << " constituents" << endl;
    int jet_k = (*constit_it)->index();
    constit_it++;
    
    // loop over the remaining particles
    while (constit_it != (*jet_it)->lastConstituent()){
      // take the last result of the merge
      int jet_i = jet_k;
      // and the next element of the jet
      int jet_j = (*constit_it)->index();
      // and merge them (with a fake dij)
      double dij = 0.0;

      // create the new jet by hand so that we can adjust its user index
      // Note again the use of the E-scheme recombination here!
      PseudoJet newjet = clust_seq.jets()[jet_i] + clust_seq.jets()[jet_j];
      clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, newjet, jet_k);

      // jump to the next constituent
      constit_it++;
    }

    // we have merged all the jet's particles into a single object, so now
    // "declare" it to be a beam (inclusive) jet.
    // [NB: put a sensible looking d_iB just to be nice...]
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);
  }
    
  // cout << "ATLASConePlugin: Bye" << endl;
  clear_list(particles_ptr);
  clear_list(jets_ptr);
}

// print a banner for reference to the 3rd-party code
void ATLASConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the ATLAS Cone plugin for FastJet                       " << endl;
  (*ostr) << "# Original code from SpartyJet; interface by the FastJet authors          " << endl;
  (*ostr) << "# If you use this plugin, please cite                                     " << endl;
  (*ostr) << "#   P.A. Delsart, K. Geerlings, J. Huston, B. Martin and C. Vermilion,    " << endl;
  (*ostr) << "#   SpartyJet, http://projects.hepforge.org/spartyjet                     " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                             " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
