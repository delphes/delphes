// $Id: ValenciaPlugin.cc 1209 2018-12-05 16:18:01Z vos $
//
// Copyright (c) 2014, Marcel Vos and Ignacio Garcia 
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "ValenciaPlugin.hh"
#include "fastjet/NNH.hh"

// strings and streams
#include <sstream>
#include <limits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
/// class that contains the algorithm parameters R and beta
//  needed to run the Valencia algorithm
class ValenciaInfo {

public: 
  ValenciaInfo(double Ri, double betai, double gammai) 
  { R_ = Ri; beta_ = betai; gamma_ = gammai;}

  double beta() { return beta_; }
  double gamma() { return gamma_; }
  double R() { return R_; }

private:

  double R_, beta_, gamma_;
};

class ValenciaBriefJet {
public:
  void init(const PseudoJet & jet, ValenciaInfo * info) {
   double norm = 1.0/sqrt(jet.modp2());
   // pseudo-jet information needed to calculate distance
    nx = jet.px() * norm;
    ny = jet.py() * norm;
    nz = jet.pz() * norm;
    E = jet.E();
   
    // algorithm parameters needed in distance calculation
    R = info->R();
    beta = info->beta();

    // beam distance
    // E-to-the-2*beta times sin(polar angle)-to-the-2*gamma
    if (E==0. || jet.perp()==0.) diB=0.;
    // modified diB in release 2.0.1    
    diB = pow(E,2*beta) * pow(jet.perp()/sqrt(jet.perp2()+jet.pz()*jet.pz()),2*info->gamma());
  }

  double distance(const ValenciaBriefJet * jet) const {
    // inter-particle distance is similar to Durham and e+e- kt
    // two differences: 
    //              R is redefined for greater flexibility
    //              the energy is elevated to 2*beta
    double dij = 1 - nx*jet->nx
                   - ny*jet->ny
                   - nz*jet->nz;

    if (pow(jet->E,2*beta) < pow(E,2*beta))
      dij *= 2 * pow(jet->E,2*beta);
    else 
      dij *= 2 * pow(E,2*beta);

    dij/=pow(R,2);

    return dij;
  }

  double beam_distance() const {
    return diB;
  }


  double E, nx, ny, nz;
  double diB;
  double R, beta;
};


std::string ValenciaPlugin::description () const {
  std::ostringstream desc;
  desc << "Valencia plugin with R = " << R() << ", beta = " << beta() << " and gamma = " << gamma();
  return desc.str();
}

void ValenciaPlugin::run_clustering(fastjet::ClusterSequence & cs) const {
  int njets = cs.jets().size();
  ValenciaInfo vinfo(R(),beta(),gamma());
  
  NNH<ValenciaBriefJet,ValenciaInfo> nnh(cs.jets(),&vinfo);

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j); 

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


} // namespace contrib

FASTJET_END_NAMESPACE
