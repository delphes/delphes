//FJSTARTHEADER
// $Id: D0RunIIConePlugin.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/D0RunIIConePlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Error.hh"
#include <sstream>

// D0 stuff
#include <list>
#include "ILConeAlgorithm.hpp"
#include "HepEntity.h"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
using namespace d0;

const double D0RunIIConePlugin::_DEFAULT_split_ratio              = 0.5  ; // overlap threshold
const double D0RunIIConePlugin::_DEFAULT_far_ratio                = 0.5  ;
const double D0RunIIConePlugin::_DEFAULT_Et_min_ratio             = 0.5  ;
const bool   D0RunIIConePlugin::_DEFAULT_kill_duplicate           = true ;
const double D0RunIIConePlugin::_DEFAULT_duplicate_dR             = 0.005; 
const double D0RunIIConePlugin::_DEFAULT_duplicate_dPT            = 0.01 ; 
const double D0RunIIConePlugin::_DEFAULT_search_factor            = 1.0  ; 
const double D0RunIIConePlugin::_DEFAULT_pT_min_leading_protojet  = 0.   ; 
const double D0RunIIConePlugin::_DEFAULT_pT_min_second_protojet   = 0.   ;
const int    D0RunIIConePlugin::_DEFAULT_merge_max                = 10000; 
const double D0RunIIConePlugin::_DEFAULT_pT_min_nomerge           = 0.   ;

bool D0RunIIConePlugin::_first_time = true;

string D0RunIIConePlugin::description () const {
  ostringstream desc;
  
  desc << "D0 Run II Improved Legacy (midpoint) cone jet algorithm, with ";
  desc << "cone_radius = "        << cone_radius        () << ", "
       << "min_jet_Et = "         << min_jet_Et         () << ", " 
       << "split_ratio = "        << split_ratio        ();

  return desc.str();
}


void D0RunIIConePlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());
 
  // create the entities needed by the D0 code
  vector<HepEntity> entities(clust_seq.jets().size());
  list<const HepEntity * > ensemble;
  for (unsigned i = 0; i < clust_seq.jets().size(); i++) {
    entities[i].Fill(clust_seq.jets()[i].E(),
                     clust_seq.jets()[i].px(),
                     clust_seq.jets()[i].py(),
                     clust_seq.jets()[i].pz(),
                     i);
    // use only the particles that do not have infinite rapidity
    if (abs(entities[i].pz) < entities[i].E) {
      ensemble.push_back(& (entities[i]));
    }
  }

  // prepare the D0 algorithm
  ILConeAlgorithm<HepEntity> 
    ilegac(cone_radius(), 
           min_jet_Et(), 
           split_ratio(),
	   far_ratio(), 
           Et_min_ratio(), 
           kill_duplicate(), 
           duplicate_dR(), 
	   duplicate_dPT(), 
           search_factor(), 
           pT_min_leading_protojet(), 
	   pT_min_second_protojet(), 
           merge_max(), 
           pT_min_nomerge());

  // run the algorithm
  float Item_ET_Threshold = 0.;
  list<HepEntity> jets;
  ilegac.makeClusters(jets, ensemble, Item_ET_Threshold);

  // now transfer the information about the jets into the
  // FastJet structure
  for(int i = ilegac.ilcv.size()-1; i >= 0; i--) {
    
    std::list<const HepEntity*> tlist = ilegac.ilcv[i].LItems();
    std::list<const HepEntity*>::iterator tk;
    
    // get first particle in list
    tk = tlist.begin();

    // if there is no particle, just discard it
    // Note: this unexpected behaviour has been observed when the
    //       min_jet_Et parameter was set to 0
    if (tk==tlist.end())
      continue;

    int jet_k = (*tk)->index;
    // now merge with remaining particles in list
    tk++;
    for (; tk != tlist.end(); tk++) {
      int jet_i = jet_k;
      int jet_j = (*tk)->index;
      // do a fake recombination step with dij=0
      double dij = 0.0;
      clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
    }

    // NB: put a sensible looking d_iB just to be nice...
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);

  }
}

// print a banner for reference to the 3rd-party code
void D0RunIIConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#--------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the D0 Run II Cone plugin for FastJet                    " << endl;
  (*ostr) << "# Original code by the D0 collaboration, provided by Lars Sonnenschein;    " << endl;
  (*ostr) << "# interface by FastJet authors                                             " << endl;
  (*ostr) << "# If you use this plugin, please cite                                      " << endl;
  (*ostr) << "#   G. C. Blazey et al., hep-ex/0005012                                    " << endl;
  (*ostr) << "#   V. M. Abazov et al. [D0 Collaboration], arXiv:1110.3771 [hep-ex]       " << endl; 
  (*ostr) << "# in addition to the usual FastJet reference.                              " << endl;
  (*ostr) << "#--------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
