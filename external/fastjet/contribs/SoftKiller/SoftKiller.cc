// $Id$
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin. P. Salam and Gregory Soyez
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

#include "SoftKiller.hh"
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// ctor with simple initialisation
///  \param rapmax     the maximal absolute rapidity extent of the grid
///  \param cell_size  the grid spacing (equivalently, cell size)
///  \param sifter     when provided, the soft killer is applied
///                    only to particles that pass the sifter (the
///                    others are kept untouched)
SoftKiller::SoftKiller(double rapmax, double cell_size,
                       Selector sifter) :
#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
  RectangularGrid(rapmax, cell_size), _sifter(sifter) {}
#else // not FJCONTRIB_SOFTKILLER_USEFJGRID
  _ymax(rapmax), _ymin(-rapmax), 
  _requested_drap(cell_size), _requested_dphi(cell_size),
  _sifter(sifter) {
    _setup_grid();
}
#endif

/// ctor with more control over initialisation
///  \param rapmin     the minimum rapidity extent of the grid
///  \param rapmax     the maximum rapidity extent of the grid
///  \param drap       the grid spacing in rapidity
///  \param dphi       the grid spacing in azimuth
///  \param sifter     when provided, the soft killer is applied
///                    onyl to particles that pass the sifter (the
///                    others are kept untouched)
SoftKiller::SoftKiller(double rapmin, double rapmax, double drap, double dphi,
                       Selector sifter) :
#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
  RectangularGrid(rapmin, rapmax, drap, dphi), _sifter(sifter) {}
#else
   _ymax(rapmax), _ymin(rapmin), 
    _requested_drap(drap), _requested_dphi(dphi),
    _sifter(sifter) {
  _setup_grid();
}
#endif

#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
SoftKiller::SoftKiller(const RectangularGrid & grid, Selector sifter) : 
            RectangularGrid(grid), _sifter(sifter) {}
#endif 

/// dummy ctor (will give an unusable SoftKiller)
SoftKiller::SoftKiller() 
#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
  {}  
#else
  : _ymax(-1.0), _ymin(1.0), _requested_drap(-1.0), _requested_dphi(-1.0) {
  _ntotal = 0;
}
#endif // FJCONTRIB_SOFTKILLER_USEFJGRID


//------------------------------------------------------------------------
// description of the soft killer
std::string SoftKiller::description() const{
  ostringstream oss;


#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
  oss << "SoftKiller with " << RectangularGrid::description();
#else
  if (_requested_drap < 0 || _requested_dphi < 0)
    return "Uninitialised SoftKiller";

  oss << "SoftKiller with rapidity extent " << _ymin << " < rap < " << _ymax
      << ", cell size drap x dphi = " << _dy << " x " << _dphi;
#endif
  if (_sifter.worker()) {
    oss << " and applied to particles passing the selection (" 
        << _sifter.description() << ")";
  }
  return oss.str();
}

//------------------------------------------------------------------------
// similarly to Transformers in FastJet, introduce a 'result'
// method equivalent to the () operator, i.e. returns the event
// after the soft killer has been applied
//vector<PseudoJet> SoftKiller::result(const vector<PseudoJet> & event) const {
void SoftKiller::apply(const vector<PseudoJet> & event, 
                       vector<PseudoJet> & reduced_event, 
                         double & pt_threshold) const {
  // a safety check: we impose at least 2 cells (otherwise, this is
  // equivalent to asking an empty event and there are more
  // efficient ways to do that)
  if (n_tiles()<2){
    throw Error("SoftKiller not properly initialised.");
  }

  // we're not set up to handle the case where the event and reduced
  // event are the same vector; so crash in that case.
  assert(&event != &reduced_event);
  // currently we can only handle the case where all tiles have equal
  // area; that is the case for Rectangular tilings, but in the future
  // one might imagine having non-rectangular tilings.
  assert(all_tiles_equal_area());

  //profiling: CPUTimer t;
  //profiling: t.start();

  // init an array to hold the max pt in each grid cell
  //
  // This is better (stack) but only C99 (fails with pedantic):
  //   double max_pt2[n_tiles()];
  // See 
  //   http://binglongx.wordpress.com/2011/05/08/create-variable-length-array-on-the-stack-in-c/
  // for a possible workaround
  vector<double> max_pt2(n_tiles(), 0.0);
  //double max_pt2[n_tiles()];
  //memset(max_pt2, 0, n_tiles()*sizeof(double));
  //double *max_pt2 = new double[n_tiles()]; // if this is reinstated, reinstate also the delete, below
  //for (int i = 0; i < n_tiles(); i++) {max_pt2[i] = 0;}

  //profiling: t.stop();
  //profiling: cout << "   copy ptrs: " << t.total() << endl;
  //profiling: t.start();

  // leave away particles that are sifted
  // vector<PseudoJet> reduced_event, saved_particles;
  // _sifter.sift(event, reduced_event, saved_particles);

  vector<const PseudoJet *> event_ptrs(event.size());
  for (unsigned i = 0; i < event.size(); i++) {
    event_ptrs[i] = & event[i];
  }
  // only run the sifter if it serves a purpose
  if (_sifter.worker()) _sifter.nullify_non_selected(event_ptrs);
    
  //profiling: t.stop();
  //profiling: cout << "   sifter   : " << t.total() << endl;
  //profiling: t.start();
    
  //vector<PseudoJet> reduced_event = event;

  // browse the particles and figure which is the min pt in each cell
  for(unsigned int i=0;i<event.size(); i++){
    if (event_ptrs[i] == 0) continue;
    int idx = tile_index(event[i]);
    if (idx<0) continue;
    max_pt2[idx] = max(max_pt2[idx], event[i].pt2());
  }
    
  //profiling: t.stop();
  //profiling: cout << "   browsing : " << t.total() << endl;
  //profiling: t.start();
  
  
  // if there are some "bad" tiles, then we need to exclude them from
  // the calculation of the median. We'll do this by condensing the
  // max_pt2 vector down to just the values for the tiles that are
  // good.
  //
  // tested answers look right in "issue" 2014-08-08-testing-rect-grid
  if (n_good_tiles() != n_tiles()) {
    int newn = 0;
    for (unsigned i = 0; i < max_pt2.size(); i++) {
      if (tile_is_good(i)) {
        // clang gets confused with the SharedPtr swap if we don't
        // have std:: here
        std::swap(max_pt2[i],max_pt2[newn]);
        newn++;
      }
    }
    max_pt2.resize(newn);
  }

  // sort the list of max values (sort works with arrays)
  //sort(max_pt2, max_pt2+n_tiles()); // use this one if we have the allocated C-style array above
  sort(max_pt2.begin(), max_pt2.end());

  // Note that if we ask that half of the event is empty that means
  // that we need at least half of the cells empty... so we need to
  // round up
  //
  // For a potentially faster median search, see
  //
  //   http://ndevilla.free.fr/median/median/index.html?utm_source=feedburner&utm_medium=twitter&utm_campaign=Feed%3A+hnycombinator+%28HN+-+hnycombinator%29
  int int_median_pos = max_pt2.size()/2;

  //--------------------------------------------------------------
  // alternative to what is below:
  //
  // apply the pt cut manually (always on pt^2): first get a vector of
  // indices that we'll want to keep and then when we know the size of
  // the resulting vector we actually do the transfer of the
  // PseudoJets.  
  //
  // This is almost a factor of two faster than doing a
  // push_back on the PseudoJets themselves (which is logical, since a
  // push_back probably averages out as doing the copy twice).
  // (2014-08-11, also tried it with pointers, which seemed 
  // marginally slower).
  double pt2cut = (1+1e-12)*max_pt2[int_median_pos];

  vector<int> indices;
  for(unsigned int i=0;i<event.size(); i++){
    if ((event_ptrs[i] == 0) || (event[i].pt2() >= pt2cut))
      indices.push_back(i);
  }
  reduced_event.reserve(indices.size());
  for(unsigned int i=0;i<indices.size(); i++){
    reduced_event.push_back(event[indices[i]]);
  }
  
  // //vector<PseudoJet> reduced_event;
  // reduced_event.clear();
  // for(unsigned int i=0;i<event.size(); i++){
  //   if ((event_ptrs[i] == 0) ||
  //       (event[i].pt2() >= pt2cut))
  //     reduced_event.push_back(event[i]);
  // }

  // free memory: reinstate this is max_pt2 becomes an allocated variable again
  //delete max_pt2;

  //return reduced_event;
  pt_threshold = sqrt(pt2cut);
    
  // end of alternative
  //--------------------------------------------------------------
    
  // double median_maxpt = sqrt(max_pt2[int_median_pos]);
  // 
  // //profiling: t.stop();
  // //profiling: cout << "   median   : " << t.total() << endl;
  // //profiling: t.start();
  // 
  // // apply a cut on pt using a selector
  // //
  // // Watch out that the Selector checks pt >= ptcut, and, since it
  // // uses pt2 to perform the comparison, may lead to rounding
  // // errors. By multiplying the ptcut by a small amount, we make
  // // sure the particle with pt=ptcut is killed.
  // //
  // // We're actually going to set to null the ones that will be kept
  // SelectorPtMax((1+1e-12)*median_maxpt).nullify_non_selected(event_ptrs);
  // 
  // // basic information
  // _ptcut = median_maxpt; // a good first start
  // 
  // //profiling: t.stop();
  // //profiling: cout << "   cutting  : " << t.total() << endl;
  // //profiling: t.start();
  // 
  // // then put back in the saved particles
  // //
  // // Note that here we may want to use the killer independently on
  // // the 2 sets of particles but then the best woulb be to move most
  // // of the above in a separate method and call it twice (left for
  // // future experimentation)
  // //copy(saved_particles.begin(), saved_particles.end(), back_inserter(reduced_event));
  // 
  // vector<PseudoJet> reduced_event;
  // // no sizeable speed gain: reduced_event.reserve(event.size());
  // for (unsigned int i=0;i<event.size();i++)
  //   if (event_ptrs[i] == 0) reduced_event.push_back(event[i]);
  // 
  // // free memory
  // delete max_pt2;
  // 
  // //profiling: t.stop();
  // //profiling: cout << "   finishing: " << t.total() << endl;
  // 
  // return reduced_event;
}

#ifndef FJCONTRIB_SOFTKILLER_USEFJGRID

// configure the grid
void SoftKiller::_setup_grid() {
  // this grid-definition code is becoming repetitive -- it should
  // probably be moved somewhere central...
  double ny_double = (_ymax-_ymin) / _requested_drap;
  _ny = max(int(ny_double+0.5),1);
  _dy = (_ymax-_ymin) / _ny;
  _inverse_dy = _ny/(_ymax-_ymin);
  
  _nphi = int (twopi / _requested_dphi + 0.5);
  _dphi = twopi / _nphi;
  _inverse_dphi = _nphi/twopi;

  // some sanity checking (could throw a fastjet::Error)
  assert(_ny >= 1 && _nphi >= 1);

  _ntotal = _nphi * _ny;
  //_max_pt.resize(_ntotal);
  _cell_area = _dy * _dphi;
}

// retrieve the grid cell index for a given PseudoJet
inline int SoftKiller::tile_index(const PseudoJet & p) const {
  // writing it as below gives a huge speed gain (factor two!). Even
  // though answers are identical and the routine here is not the
  // speed-critical step. It's not at all clear why.
  int iy = int(floor( (p.rap() - _ymin) * _inverse_dy ));
  if (iy < 0 || iy >= _ny) return -1;

  int iphi = int( p.phi() * _inverse_dphi );
  //assert(iphi >= 0 && iphi <= _nphi);
  if (iphi == _nphi) iphi = 0; // just in case of rounding errors

  //int igrid_res = iy*_nphi + iphi;
  //assert (igrid_res >= 0 && igrid_res < _ny*_nphi);
  return iy*_nphi + iphi; //igrid_res;
}

#endif // not FJCONTRIB_SOFTKILLER_USEFJGRID



} // namespace contrib

FASTJET_END_NAMESPACE
