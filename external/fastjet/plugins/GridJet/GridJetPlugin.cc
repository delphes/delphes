//FJSTARTHEADER
// $Id: GridJetPlugin.cc 2268 2011-06-20 15:12:26Z salam $
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

// fastjet stuff
#include "fastjet/ClusterSequence.hh"
#include "fastjet/GridJetPlugin.hh"

// other stuff
#include <vector>
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

//----------------------------------------------------------------------
GridJetPlugin::GridJetPlugin (double ymax,
                              double requested_grid_spacing, 
                              const JetDefinition & post_jet_def) :
#ifdef FASTJET_GRIDJET_USEFJGRID
  RectangularGrid(ymax, requested_grid_spacing), _post_jet_def(post_jet_def) {
}
#else
  _ymin(-ymax), _ymax(ymax), 
  _requested_grid_spacing(requested_grid_spacing) ,
  _post_jet_def(post_jet_def)
{
  setup_grid();
}
#endif

#ifdef FASTJET_GRIDJET_USEFJGRID
GridJetPlugin::GridJetPlugin (const RectangularGrid & grid,
                              const JetDefinition & post_jet_def) : 
  RectangularGrid(grid), _post_jet_def(post_jet_def) {
  if (!RectangularGrid::is_initialised()) 
    throw Error("attempt to construct GridJetPlugin with uninitialised RectangularGrid");
}
#endif // FASTJET_GRIDJET_USEFJGRID

#ifndef FASTJET_GRIDJET_USEFJGRID
void GridJetPlugin::setup_grid() {
  // since we've exchanged the arguments of the constructor,
  // there's a danger of calls with exchanged ymax,spacing arguments -- 
  // the following check should catch most such situations.
  assert(_ymax>0 && _ymax - _ymin >= _requested_grid_spacing);

  double ny_double = (_ymax-_ymin) / _requested_grid_spacing;
  _ny = int(ny_double+0.49999);
  _dy = (_ymax-_ymin) / _ny;
  
  _nphi = int (twopi / _requested_grid_spacing + 0.5);
  _dphi = twopi / _nphi;

  // some sanity checking (could throw a fastjet::Error)
  assert(_ny >= 1 && _nphi >= 1);

  _ntotal = _nphi * _ny;
}

//----------------------------------------------------------------------
int GridJetPlugin::tile_index(const PseudoJet & p) const {
  // directly taking int does not work for values between -1 and 0
  // so use floor instead
  // double iy_double = (p.rap() - _ymin) / _dy;
  // if (iy_double < 0.0) return -1;
  // int iy = int(iy_double);
  // if (iy >= _ny) return -1;

  // writing it as below gives a huge speed gain (factor two!). Even
  // though answers are identical and the routine here is not the
  // speed-critical step. It's not at all clear why.
  int iy = int(floor( (p.rap() - _ymin) / _dy ));
  if (iy < 0 || iy >= _ny) return -1;

  int iphi = int( p.phi()/_dphi );
  assert(iphi >= 0 && iphi <= _nphi);
  if (iphi == _nphi) iphi = 0; // just in case of rounding errors

  int igrid_res = iy*_nphi + iphi;
  assert (igrid_res >= 0 && igrid_res < _ny*_nphi);
  return igrid_res;
}
#endif // not FASTJET_GRIDJET_USEFJGRID


//----------------------------------------------------------------------
string GridJetPlugin::description () const {
  ostringstream desc;
  desc << "GridJetPlugin plugin with ";
#ifndef FASTJET_GRIDJET_USEFJGRID
  desc << "ymax = " << _ymax << ", dy = " << _dy << ", dphi = " << _dphi << " (requested grid spacing was " << _requested_grid_spacing << ")";
#else
  desc << RectangularGrid::description();
#endif
  if (_post_jet_def.jet_algorithm() != undefined_jet_algorithm) {
    desc << ", followed by " << _post_jet_def.description();
  }
  return desc.str();
}


//----------------------------------------------------------------------
double GridJetPlugin::R() const {return sqrt(drap()*dphi()/pi);}


//----------------------------------------------------------------------
void GridJetPlugin::run_clustering(ClusterSequence & cs) const {
  
  // we will create a grid; 
  //  * -1 will indicate there is no jet here currently
  //  * a number >= 0 will mean that particle indicated by the index
  //    is currently the jet on the grid
  vector<int> grid(n_tiles(), -1);
 
  int nparticles = cs.jets().size();
  double dij_or_diB = 1.0;

  int ngrid_active = 0;

  // combine particles with whatever is in the grid
  for (int i = 0; i < nparticles; i++) {
    int igrd = tile_index(cs.jets()[i]);
    //cout << i << " " << cs.jets()[i].rap() << " " << cs.jets()[i].phi() 
    // 	 << " " << igrd << " " << grid.size() << " " << _ntotal << endl;
    if (igrd < 0) continue;
    assert(igrd <= n_tiles());
    if (grid[igrd] == -1) {
      grid[igrd] = i; // jet index of initial particle i is i
      ngrid_active++;
    } else {
      int k;
      cs.plugin_record_ij_recombination(grid[igrd], i, dij_or_diB, k);
      grid[igrd] = k; // grid takes jet index of new particle
      //cout << "  res: " << cs.jets()[k].rap() << " " << cs.jets()[k].phi() << endl;
    }
  }

  if (_post_jet_def.jet_algorithm() == undefined_jet_algorithm) {
    // make the final jets via iB recombinations
    for (unsigned igrd = 0; igrd < grid.size(); igrd++) {
      if (grid[igrd] != -1 && tile_is_good(igrd)) 
                   cs.plugin_record_iB_recombination(grid[igrd], dij_or_diB);
    }
  } else {
    // otherwise post-cluster the grid elements with a normal jet algorithm
    vector<PseudoJet> inputs;
    vector<int>       cs_indices;
    inputs.reserve(ngrid_active);
    cs_indices.reserve(2*ngrid_active);
    for (unsigned igrd = 0; igrd < grid.size(); igrd++) {
      if (grid[igrd] != -1) {
	inputs.push_back(cs.jets()[grid[igrd]]);
	cs_indices.push_back(grid[igrd]);
      }
    }
    ClusterSequence post_cs(inputs, _post_jet_def);
    const vector<ClusterSequence::history_element> & post_history = post_cs.history();
    const vector<PseudoJet>                        & post_jets = post_cs.jets();
    for (unsigned ihist = ngrid_active; ihist < post_history.size(); ihist++) {
      const ClusterSequence::history_element & hist = post_history[ihist];
      int post_ij1 = post_history[hist.parent1].jetp_index;
      int ij1 = cs_indices[post_ij1];
      if (hist.parent2 >= 0) {
	int post_ij2 = post_history[hist.parent2].jetp_index;
	int ij2 = cs_indices[post_ij2];
	int k;
	cs.plugin_record_ij_recombination(ij1, ij2, hist.dij, post_jets[hist.jetp_index], k); 
	assert(int(cs_indices.size()) == hist.jetp_index);
	cs_indices.push_back(k);
      } else {
	cs.plugin_record_iB_recombination(ij1, hist.dij);
      }
    }
    
  }
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
