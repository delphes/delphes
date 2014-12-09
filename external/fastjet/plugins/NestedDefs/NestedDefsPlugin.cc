//FJSTARTHEADER
// $Id: NestedDefsPlugin.cc 3433 2014-07-23 08:17:03Z salam $
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

// TODO
// ? Maybe one could provide additional recomb. dists as an "extra".;

// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/NestedDefsPlugin.hh"

// other stuff
#include <vector>
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

string NestedDefsPlugin::description () const {
  ostringstream desc;
  
  desc << "NestedDefs: successive application of " ;
  unsigned int i=1;
  for (list<JetDefinition>::const_iterator it=_defs.begin();it!=_defs.end();it++){
    desc << "Definition " << i++ << " [" << it->description() << "] - ";
  }

  return desc.str();
}

void NestedDefsPlugin::run_clustering(ClusterSequence & clust_seq) const {
  vector<PseudoJet> momenta;

  // build the initial list of particles
  momenta = clust_seq.jets();
  unsigned int step_n = momenta.size();

  // initialise the conversion table, which works as follows
  // conversion_table[step_cs_jet_index] = main_cs_jet_index
  vector<unsigned int> conversion_table(2*step_n);
  vector<unsigned int> new_conversion_table;
  for (unsigned int i=0;i<step_n;i++)
    conversion_table[i]=i;

  // Now the steps go as follows:
  // for each definition in the list, 
  //  - do the clustering,
  //  - copy the history into the main one
  //  - update the list of momenta and the index conversion table
  list<JetDefinition>::const_iterator def_iterator = _defs.begin();
  unsigned int def_index=0;
  bool last_def=false;

  while (def_iterator!=_defs.end()){
    last_def = (def_index == (_defs.size()-1));

    // do the clustering
    ClusterSequence step_cs(momenta, *def_iterator);

    // clear the momenta as we shall fill them again
    momenta.clear();
    new_conversion_table.clear();

    // retrieve the history
    const vector<ClusterSequence::history_element> & step_history = step_cs.history();

    // copy the history
    // note that we skip the initial steps which are just the 
    // declaration of the particles.
    vector<ClusterSequence::history_element>::const_iterator 
      hist_iterator = step_history.begin();

    for (unsigned int i=step_n;i!=0;i--)
      hist_iterator++;

    while (hist_iterator != step_history.end()){
      // check if it is a recombination with the beam or a simple recombination
      if (hist_iterator->parent2 == ClusterSequence::BeamJet){
	// save this jet for future clustering
	// unless we've reached the last def in which case, record the clustering
	unsigned int step_jet_index = step_cs.history()[hist_iterator->parent1].jetp_index;
	if (last_def){
	  clust_seq.plugin_record_iB_recombination(conversion_table[step_jet_index], 
						   hist_iterator->dij);
	} else {
	  momenta.push_back(step_cs.jets()[step_jet_index]);
	  new_conversion_table.push_back(conversion_table[step_jet_index]);
	}
      } else {
	// record combination
	// note that we set the recombination distance to 0 except for the last alg
	unsigned int step_jet1_index = step_cs.history()[hist_iterator->parent1].jetp_index;
	unsigned int step_jet2_index = step_cs.history()[hist_iterator->parent2].jetp_index;
	PseudoJet newjet = step_cs.jets()[hist_iterator->jetp_index];
	int jet_k;
	clust_seq.plugin_record_ij_recombination(conversion_table[step_jet1_index], 
						 conversion_table[step_jet2_index],
						 last_def ? hist_iterator->dij : 0.0,
						 newjet, jet_k);

	// save info in the conversion table for tracking purposes
	conversion_table[hist_iterator->jetp_index]=jet_k;
      }

      // go to the next history element
      hist_iterator++;
    }

    // finalise this step:
    //  - update nr of particles
    //  - update conversion table
    step_n = momenta.size();
    for (unsigned int i=0;i<step_n;i++)
      conversion_table[i] = new_conversion_table[i];

    // go to the next alg
    def_index++;
    def_iterator++;
  }

}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
