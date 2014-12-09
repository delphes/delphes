//FJSTARTHEADER
// $Id: EECambridgePlugin.cc 3433 2014-07-23 08:17:03Z salam $
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
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/NNH.hh"

// other stuff
#include <sstream>
#include <limits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

//----------------------------------------------------------------------
/// class to help run an e+e- Cambridge algorithm
class EECamBriefJet {
public:
  void init(const PseudoJet & jet) {
    double norm = 1.0/sqrt(jet.modp2());
    nx = jet.px() * norm;
    ny = jet.py() * norm;
    nz = jet.pz() * norm;
  }

  double distance(const EECamBriefJet * jet) const {
    double dij = 1 - nx*jet->nx
                   - ny*jet->ny
                   - nz*jet->nz;
    return dij;
  }

  double beam_distance() const {
    return numeric_limits<double>::max();
  }

private:
  double nx, ny, nz;
};


string EECambridgePlugin::description () const {
  ostringstream desc;
  desc << "EECambridge plugin with ycut = " << ycut() ;
  return desc.str();
}

void EECambridgePlugin::run_clustering(ClusterSequence & cs) const {
  int njets = cs.jets().size();
  NNH<EECamBriefJet> nnh(cs.jets());

  double Q2 = cs.Q2(); 

  while (njets > 0) {
    int i, j, k;
    // here we get a minimum based on the purely angular variable from the NNH class
    // (called dij there, but vij in the Cambridge article (which uses dij for 
    // a kt distance...)
    double vij = nnh.dij_min(i, j); // i,j are return values...

    // next we work out the dij (ee kt distance), and based on its
    // value decide whether we have soft-freezing (represented here by
    // a "Beam" clustering) or not
    double dij;
    if (j >= 0) {
      double scale = min(cs.jets()[i].E(), cs.jets()[j].E());
      dij = 2 * vij * scale * scale;
      if (dij > Q2 * ycut()) {
	// we'll call the softer partner a "beam" jet
	if (cs.jets()[i].E() > cs.jets()[j].E()) std::swap(i,j);
	j = -1;
      }
    } else {
      // for the last particle left, just use yij = 1
      dij = Q2;
    }
    
    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
    } else {
      cs.plugin_record_iB_recombination(i, dij);
      nnh.remove_jet(i);
    }
    njets--;
  }
    
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
