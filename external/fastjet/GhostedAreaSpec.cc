//FJSTARTHEADER
// $Id: GhostedAreaSpec.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Error.hh"
#include<iostream>
#include<sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

BasicRandom<double> GhostedAreaSpec::_random_generator;
LimitedWarning GhostedAreaSpec::_warn_fj2_placement_deprecated;

/// explicit constructor
GhostedAreaSpec::GhostedAreaSpec(
                           const Selector & selector,
                           int    repeat_in        ,
                           double ghost_area_in    ,   
                           double grid_scatter_in  , 
                           double pt_scatter_in    ,   
                           double mean_ghost_pt_in 
                          ): 
    _repeat(repeat_in), 
    _ghost_area(ghost_area_in), 
    _grid_scatter(grid_scatter_in),  
    _pt_scatter(pt_scatter_in), 
    _mean_ghost_pt(mean_ghost_pt_in),
    _fj2_placement(false),
    _selector(selector),
    _actual_ghost_area(-1.0)
  {
    // check the selector has the properties needed -- an area and
    // applicability jet-by-jet (the latter follows automatically from
    // the former?)
    if (!_selector.has_finite_area()) throw Error("To construct a GhostedAreaSpec with a Selector, the selector must have a finite area");
    if (!_selector.applies_jet_by_jet()) throw Error("To construct a GhostedAreaSpec with a Selector, the selector must apply jet-by-jet");
    // get the internal rapidity extent from the selector
    double ghost_maxrap_local, ghost_minrap_local;
    _selector.get_rapidity_extent(ghost_minrap_local, ghost_maxrap_local);
    _ghost_maxrap     = 0.5*(ghost_maxrap_local - ghost_minrap_local); 
    _ghost_rap_offset = 0.5*(ghost_maxrap_local + ghost_minrap_local);
    
    _initialize();
  
}

//======================================================================
// sets fj2 ghost placement
void GhostedAreaSpec::set_fj2_placement(bool val) {
  _fj2_placement  = val; _initialize();
  if (val) _warn_fj2_placement_deprecated.warn("FJ2 placement of ghosts can lead to systematic edge effects in area evaluation and is deprecated. Prefer new (default) FJ3 placement.");
}

//======================================================================
/// sets the detailed parameters for the ghosts (which may not be quite
/// the same as those requested -- this is in order for things to fit
/// in nicely into 2pi etc...
void GhostedAreaSpec::_initialize() {
  // add on area-measuring dummy particles
  _drap = sqrt(_ghost_area);
  _dphi = _drap;
  if (_fj2_placement) {
    _nphi = int(ceil(twopi/_dphi)); _dphi = twopi/_nphi;
    _nrap = int(ceil(_ghost_maxrap/_drap)); _drap = _ghost_maxrap / _nrap;
    _actual_ghost_area = _dphi * _drap;
    _n_ghosts   = (2*_nrap+1)*_nphi;
  } else {
    // for FJ3, update the ghost placement as follows
    // - use nearest int rather than ceiling in determining number of
    //   phi and rapidity locations, because this is more stable when
    //   the user is trying to get an exact number based on the area
    // - rather than placing ghosts up to maximum rapidity
    _nphi = int(twopi/_dphi + 0.5); _dphi = twopi/_nphi;
    _nrap = int(_ghost_maxrap/_drap + 0.5); _drap = _ghost_maxrap / _nrap;
    _actual_ghost_area = _dphi * _drap;
    _n_ghosts   = (2*_nrap)*_nphi;
  }
  // checkpoint the status of the random number generator.
  checkpoint_random();
  //_random_generator.info(cerr);
}

//----------------------------------------------------------------------
/// adds the ghost 4-momenta to the vector of PseudoJet's
void GhostedAreaSpec::add_ghosts(vector<PseudoJet> & event) const {

  double rap_offset;
  int nrap_upper;
  if (_fj2_placement) {
    rap_offset  = 0.0;
    nrap_upper  = _nrap;
  } else {
    rap_offset  = 0.5;
    nrap_upper  = _nrap-1;
  }

  // add momenta for ghosts
  for (int irap = -_nrap; irap <= nrap_upper; irap++) {
    for (int iphi = 0; iphi < _nphi; iphi++) {
     
      // include random offsets for all quantities
      //----------------------------------------------
      // NB: in FJ2 we'd exchanged the px and py components relative to a
      // standard definition of phi; to preserve the same areas as fj2
      // we now generate a "phi_fj2", and then convert to a standard phi
      double phi_fj2 = (iphi+0.5) * _dphi + _dphi*(_our_rand()-0.5)*_grid_scatter;
      double phi;
      if (_fj2_placement) phi = 0.5*pi - phi_fj2;
      else                phi = phi_fj2;
      double rap = (irap+rap_offset) * _drap + _drap*(_our_rand()-0.5)*_grid_scatter
	                                                 + _ghost_rap_offset ;
      double pt = _mean_ghost_pt*(1+(_our_rand()-0.5)*_pt_scatter);

      double exprap = exp(+rap);
      double pminus = pt/exprap;
      double pplus  = pt*exprap;
      double px = pt*cos(phi);
      double py = pt*sin(phi);
      PseudoJet mom(px,py,0.5*(pplus-pminus),0.5*(pplus+pminus));
      // this call fills in the PseudoJet's cached rap,phi information,
      // based on pre-existing knowledge. Watch out: if you get the hint
      // wrong nobody will tell you, but you will certainly mess up
      // your results.
      mom.set_cached_rap_phi(rap,phi);

      // if we have an active selector and the particle does not pass the 
      // selection condition, move on to the next momentum
      if (_selector.worker().get() && !_selector.pass(mom)) continue;
      event.push_back(mom);
    }
  }
}

string GhostedAreaSpec::description() const {

  ostringstream ostr;
  ostr << "ghosts of area " << actual_ghost_area() 
       << " (had requested " << ghost_area() << ")";
  if (_selector.worker().get()) 
    ostr << ", placed according to selector (" << _selector.description() << ")";
  else
    ostr << ", placed up to y = " << ghost_maxrap() ;
  ostr << ", scattered wrt to perfect grid by (rel) " << grid_scatter() 
       << ", mean_ghost_pt = " << mean_ghost_pt()
       << ", rel pt_scatter =  " << pt_scatter()
       << ", n repetitions of ghost distributions =  " << repeat();
  return ostr.str();
}

FASTJET_END_NAMESPACE

