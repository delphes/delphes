//FJSTARTHEADER
// $Id: JadePlugin.cc 4354 2018-04-22 07:12:37Z salam $
//
// Copyright (c) 2007-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
#include "fastjet/NNFJN2Plain.hh"

// other stuff
#include <vector>
#include <sstream>
#include <limits>




using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
/// class to help run a JADE algorithm
///
/// This class works both with NNH and NNFJN2Plain clustering
/// helpers. They both use the same init(...) call, but for the
/// clustering:
///
/// - NNH uses distance(...) and beam_distance()
/// - NNFJPlainN2 uses geometrical_distance(...), momentum_factor()
///   and geometrical_beam_distance()
///
/// For NNFJPlainN2 the 2 E_i E_j (1-cos theta_{ij}) factor
/// gets broken up into
///
///     sqrt(2)*min(E_i,E_j) * [sqrt(2)*max(E_i,E_j) (1 - cos \theta_{ij})]
///
/// The second factor is what we call the "geometrical_distance" even
/// though it isn't actually purely geometrical. But the fact that it
/// gets multiplied by min(E_i,E_j) to get the full distance is
/// sufficient for the validity of the FJ lemma, allowing for the use
/// of NNFJN2Plain.
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

  double geometrical_distance(const JadeBriefJet * jet) const {
    double dij = 1 - nx*jet->nx
                   - ny*jet->ny
                   - nz*jet->nz;
    dij *= max(rt2E,jet->rt2E);
    return dij;
  }

  double momentum_factor() const {
    return rt2E;
  }
  
  double beam_distance() const {
    return numeric_limits<double>::max();
  }

  double geometrical_beam_distance() const {
    // get a number that is almost the same as max(), just a little
    // smaller so as to ensure that when we divide it by rt2E and then
    // multiply it again, we won't get an overflow
    const double almost_max = numeric_limits<double>::max() * (1 - 1e-13);
    return almost_max / rt2E;
  }
  
private:
  double rt2E, nx, ny, nz;
};


//----------------------------------------------------------------------
string JadePlugin::description () const {
  ostringstream desc;
  desc << "e+e- JADE algorithm plugin";
  switch(_strategy) {
  case strategy_NNH:
    desc << ", using NNH strategy"; break;
  case strategy_NNFJN2Plain:
    desc << ", using NNFJN2Plain strategy"; break;
  default:
    throw Error("Unrecognized strategy in JadePlugin");
  }

  return desc.str();
}

// //----------------------------------------------------------------------
// void JadePlugin::run_clustering(ClusterSequence & cs) const {
//   int njets = cs.jets().size();
// 
//   //SharedPtr<NNBase<> > nn;
//   NNBase<> * nn;
//   switch(_strategy) {
//   case strategy_NNH:
//     //nn.reset(new NNH<JadeBriefJet>(cs.jets()));
//     nn = new NNH<JadeBriefJet>(cs.jets());
//     break;
//   case strategy_NNFJN2Plain:
//     //nn.reset(new NNFJN2Plain<JadeBriefJet>(cs.jets()));
//     nn = new NNFJN2Plain<JadeBriefJet>(cs.jets());
//     break;
//   default:
//     throw Error("Unrecognized strategy in JadePlugin");
//   }
//   //NNH<JadeBriefJet> nnh(cs.jets());
//   //NNFJN2Plain<JadeBriefJet> nnh(cs.jets());
// 
//   // if testing against Hoeth's implementation, need to rescale the
//   // dij by Q^2.
//   //double Q2 = cs.Q2(); 
// 
//   while (njets > 0) {
//     int i, j, k;
//     double dij = nn->dij_min(i, j);
// 
//     if (j >= 0) {
//       cs.plugin_record_ij_recombination(i, j, dij, k);
//       nn->merge_jets(i, j, cs.jets()[k], k);
//     } else {
//       double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
//       cs.plugin_record_iB_recombination(i, diB);
//       nn->remove_jet(i);
//     }
//     njets--;
//   }
//   delete nn;
// }


template<class N> void JadePlugin::_actual_run_clustering(ClusterSequence & cs) const {

  int njets = cs.jets().size();

  N nn(cs.jets());

  // if testing against Hoeth's implementation, need to rescale the
  // dij by Q^2.
  //double Q2 = cs.Q2(); 

  while (njets > 0) {
    int i, j, k;
    double dij = nn.dij_min(i, j);

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nn.merge_jets(i, j, cs.jets()[k], k);
    } else {
      double diB = cs.jets()[i].E()*cs.jets()[i].E(); // get new diB
      cs.plugin_record_iB_recombination(i, diB);
      nn.remove_jet(i);
    }
    njets--;
  }

}

//----------------------------------------------------------------------
void JadePlugin::run_clustering(ClusterSequence & cs) const {

  switch(_strategy) {
  case strategy_NNH:
    _actual_run_clustering<NNH<JadeBriefJet> >(cs);
    break;
  case strategy_NNFJN2Plain:
    _actual_run_clustering<NNFJN2Plain<JadeBriefJet> >(cs);
    break;
  default:
    throw Error("Unrecognized strategy in JadePlugin");
  }
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
