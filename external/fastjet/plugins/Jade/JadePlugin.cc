//FJSTARTHEADER
// $Id: JadePlugin.cc 3433 2014-07-23 08:17:03Z salam $
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
#include "fastjet/JadePlugin.hh"
#include <iostream>
//#include "fastjet/internal/ClusterSequence_N2.icc"
#include "fastjet/NNH.hh"

// other stuff
#include <vector>
#include <sstream>
#include <limits>




using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
/// class to help run a JADE algorithm
class JadeBriefJet {
public:
  void init(const PseudoJet & jet) {
    double norm = 1.0/sqrt(jet.modp2());
    nx = jet.px() * norm;
    ny = jet.py() * norm;
    nz = jet.pz() * norm;
    rt2E = sqrt(2.0)*jet.E();
  }

  double distance(const JadeBriefJet * jet) const {
    double dij = 1 - nx*jet->nx
                   - ny*jet->ny
                   - nz*jet->nz;
    dij *= rt2E*jet->rt2E;
    return dij;
  }

  double beam_distance() const {
    return numeric_limits<double>::max();
  }

private:
  double rt2E, nx, ny, nz;
};


//----------------------------------------------------------------------
string JadePlugin::description () const {
  ostringstream desc;
  desc << "e+e- JADE algorithm plugin";
  return desc.str();
}

//----------------------------------------------------------------------
void JadePlugin::run_clustering(ClusterSequence & cs) const {
  int njets = cs.jets().size();
  NNH<JadeBriefJet> nnh(cs.jets());

  // if testing against Hoeth's implementation, need to rescale the
  // dij by Q^2.
  //double Q2 = cs.Q2(); 

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j);

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
    } else {
      double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
      cs.plugin_record_iB_recombination(i, diB);
      nnh.remove_jet(i);
    }
    njets--;
  }
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
