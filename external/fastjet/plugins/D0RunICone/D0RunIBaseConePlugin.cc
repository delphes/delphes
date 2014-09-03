//FJSTARTHEADER
// $Id: D0RunIBaseConePlugin.cc 1779 2010-10-25 10:32:59Z soyez $
//
// Copyright (c) 2009-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

// D0 stuff
// apparently this has to go first to avoid a problem with gcc-4.0.1 builds on Macs
#include <list>
#include "ConeClusterAlgo.hpp"
#include "HepEntityIpre96.h"
#include "HepEntityI.h"

#include "fastjet/D0RunIBaseConePlugin.hh"
#include "fastjet/D0RunIpre96ConePlugin.hh"
#include "fastjet/D0RunIConePlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Error.hh"
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std; 
using namespace d0runi;

/////////////////////////////////////////////
//                                         //
// D0RunIBaseConePlugin implementation     //
//                                         //
/////////////////////////////////////////////

const double D0RunIBaseConePlugin::_DEFAULT_SPLifr             = 0.5; //shared Et fraction threshold
const double D0RunIBaseConePlugin::_DEFAULT_TWOrad             = 0.; 
const bool   D0RunIBaseConePlugin::_DEFAULT_D0_Angle           = false;
const bool   D0RunIBaseConePlugin::_DEFAULT_Increase_Delta_R   = true;
const bool   D0RunIBaseConePlugin::_DEFAULT_Kill_Far_Clusters  = true;
const bool   D0RunIBaseConePlugin::_DEFAULT_Jet_Et_Min_On_Iter = true;
const double D0RunIBaseConePlugin::_DEFAULT_Far_Ratio          = 0.5;
const double D0RunIBaseConePlugin::_DEFAULT_Eitem_Negdrop      = -1.0;
const double D0RunIBaseConePlugin::_DEFAULT_Et_Min_Ratio       = 0.5;
const double D0RunIBaseConePlugin::_DEFAULT_Thresh_Diff_Et     = 0.01;


// for the real work, we write a template class that decides which
// HepEntity type to use
template<typename HepEntityType>
void D0RunIBaseConePlugin::run_clustering_worker(ClusterSequence & clust_seq) const{
  // create the entities needed by the D0 code
  vector<HepEntityType> entities(clust_seq.jets().size());
  list<const HepEntityType * > ensemble;
  for (unsigned i = 0; i < clust_seq.jets().size(); i++) {
    entities[i].Fill(clust_seq.jets()[i].E(),
		     clust_seq.jets()[i].px(),
		     clust_seq.jets()[i].py(),
		     clust_seq.jets()[i].pz(),
		     i);
    // use only the particles that do not have infinite rapidity
    if (abs(entities[i].pz() ) < entities[i].E() ) {
      ensemble.push_back(& (entities[i]));
    }
  }

  // prepare the D0 algorithm
  ConeClusterAlgo<HepEntityType> 
    RunIconeAlgo(CONErad(), 
		 JETmne(),
		 TWOrad(),
		 SPLifr(),
		 D0_Angle(),
		 Increase_Delta_R(),
		 Kill_Far_Clusters(),
		 Jet_Et_Min_On_Iter(),
		 Far_Ratio(),
		 Eitem_Negdrop(),
		 Et_Min_Ratio(),
		 Thresh_Diff_Et());


  // run the algorithm
  float Zvertex = 0.;
  list<HepEntityType> jets;
  RunIconeAlgo.makeClusters(jets, ensemble, Zvertex);

  // now transfer the information about the jets into the
  // FastJet structure
  for(int i = RunIconeAlgo.TempColl.size()-1; i >= 0; i--) {
    
    std::list<const HepEntityType*> tlist = RunIconeAlgo.TempColl[i].LItems();
    typename std::list<const HepEntityType*>::iterator tk;
        
    // get first particle in list
    tk = tlist.begin();
    int jet_k = (*tk)->index;

    // GS addition: in order to use the proper recombination scheme
    // used by D0 we need to keep track of the sum as a
    // "HepEntityType"
    HepEntityType jet_current_momentum = *(*tk);

    // now merge with remaining particles in list
    tk++;
    for (; tk != tlist.end(); tk++) {
      int jet_i = jet_k;
      int jet_j = (*tk)->index;
      // do a fake recombination step with dij=0
      double dij = 0.0;

      // GS addition: find the new momentum and convert that into a
      // pseudo-jet
      jet_current_momentum.Add(**tk);
      PseudoJet new_mom(jet_current_momentum.px(), jet_current_momentum.py(), 
			jet_current_momentum.pz(), jet_current_momentum.E());

      clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, new_mom, jet_k);
    }
    
    // NB: put a sensible looking d_iB just to be nice...
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);    
  }
}


/////////////////////////////////////////////
//                                         //
// D0RunIpre96ConePlugin implementation    //
//                                         //
/////////////////////////////////////////////

bool D0RunIpre96ConePlugin::_first_time=true;

string D0RunIpre96ConePlugin::description () const {
  ostringstream desc;
  
  desc << "D0 Run I (pre 96) cone jet algorithm, with ";
  desc << "cone_radius = "        << CONErad        () << ", "
       << "min_jet_Et = "         << JETmne         () << ", " 
       << "split_fraction = "     << SPLifr         ();

  return desc.str();
}

void D0RunIpre96ConePlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());

  run_clustering_worker<HepEntityIpre96>(clust_seq);
}

// print a banner for reference to the 3rd-party code
void D0RunIpre96ConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#--------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the D0 Run I (pre96) Cone plugin for FastJet             " << endl;
  (*ostr) << "# Original code by the D0 collaboration, provided by Lars Sonnenschein;    " << endl;
  (*ostr) << "# interface by FastJet authors                                             " << endl;
  (*ostr) << "# If you use this plugin, please cite                                      " << endl;
  (*ostr) << "#   B. Abbott et al. [D0 Collaboration], FERMILAB-PUB-97-242-E.            " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                              " << endl;
  (*ostr) << "#--------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}


/////////////////////////////////////////////
//                                         //
// D0RunIConePlugin implementation         //
//                                         //
/////////////////////////////////////////////

bool D0RunIConePlugin::_first_time=true;

string D0RunIConePlugin::description () const {
  ostringstream desc;
  
  desc << "D0 Run I cone jet algorithm, with ";
  desc << "cone_radius = "        << CONErad        () << ", "
       << "min_jet_Et = "         << JETmne         () << ", " 
       << "split_fraction = "     << SPLifr         ();

  return desc.str();
}

void D0RunIConePlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  _print_banner(clust_seq.fastjet_banner_stream());

  run_clustering_worker<HepEntityI>(clust_seq);
}

// print a banner for reference to the 3rd-party code
void D0RunIConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time) return;
  _first_time=false;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#--------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the D0 Run I Cone plugin for FastJet                     " << endl;
  (*ostr) << "# Original code provided by Lars Sonnenschein; interface by FastJet authors" << endl;
  (*ostr) << "# If you use this plugin, please cite                                      " << endl;
  (*ostr) << "#   B. Abbott et al. [D0 Collaboration], FERMILAB-PUB-97-242-E.            " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                              " << endl;
  (*ostr) << "#--------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}



FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
