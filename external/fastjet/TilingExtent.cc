//FJSTARTHEADER
// $Id: TilingExtent.cc 3433 2014-07-23 08:17:03Z salam $
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

#include <iomanip>
#include <limits>
#include <cmath>
#include "fastjet/internal/TilingExtent.hh"
using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

TilingExtent::TilingExtent(ClusterSequence & cs) {
  _determine_rapidity_extent(cs.jets());
}
  
void TilingExtent::_determine_rapidity_extent(const vector<PseudoJet> & particles) {
  // have a binning of rapidity that goes from -nrap to nrap
  // in bins of size 1; the left and right-most bins include
  // include overflows from smaller/larger rapidities
  int nrap = 20; 
  int nbins = 2*nrap;
  vector<double> counts(nbins, 0);
  
  // get the minimum and maximum rapidities and at the same time bin
  // the multiplicities as a function of rapidity to help decide how
  // far out it's worth going
  _minrap =  numeric_limits<double>::max();
  _maxrap = -numeric_limits<double>::max();
  int ibin;
  for (unsigned i = 0; i < particles.size(); i++) {
    // ignore particles with infinite rapidity
    if (particles[i].E() == abs(particles[i].pz())) continue;
    double rap = particles[i].rap();
    if (rap < _minrap) _minrap = rap;
    if (rap > _maxrap) _maxrap = rap;
    // now bin the rapidity to decide how far to go with the tiling.
    // Remember the bins go from ibin=0 (rap=-infinity..-19)
    // to ibin = nbins-1 (rap=19..infinity for nrap=20)
    ibin = int(rap+nrap); 
    if (ibin < 0) ibin = 0;
    if (ibin >= nbins) ibin = nbins - 1;
    counts[ibin]++;
  }

  // now figure out the particle count in the busiest bin
  double max_in_bin = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    if (max_in_bin < counts[ibin]) max_in_bin = counts[ibin];
  }
  
  // and find _minrap, _maxrap such that edge bin never contains more
  // than some fraction of busiest, and at least a few particles; first do
  // it from left. NB: the thresholds chosen here are largely
  // guesstimates as to what might work.
  //
  // 2014-07-17: in some tests at high multiplicity (100k) and particles going up to
  //             about 7.3, anti-kt R=0.4, we found that 0.25 gave 20% better run times
  //             than the original value of 0.5.
  const double allowed_max_fraction = 0.25;
  // the edge bins should also contain at least min_multiplicity particles
  const double min_multiplicity = 4;
  // now calculate how much we can accumulate into an edge bin
  double allowed_max_cumul = floor(max(max_in_bin * allowed_max_fraction, min_multiplicity));
  // make sure we don't require more particles in a bin than max_in_bin
  if (allowed_max_cumul > max_in_bin) allowed_max_cumul = max_in_bin;

  // start scan over rapidity bins from the left, to find out minimum rapidity of tiling
  double cumul_lo = 0;
  _cumul2 = 0;
  for (ibin = 0; ibin < nbins; ibin++) {
    cumul_lo += counts[ibin];
    if (cumul_lo >= allowed_max_cumul) {
      double y = ibin-nrap;
      if (y > _minrap) _minrap = y;
      break;
    }
  }
  assert(ibin != nbins); // internal consistency check that you found a bin
  _cumul2 += cumul_lo*cumul_lo;

  // ibin_lo is the index of the leftmost bin that should be considered
  int ibin_lo = ibin;

  // then do it from right, to find out maximum rapidity of tiling
  double cumul_hi = 0;
  for (ibin = nbins-1; ibin >= 0; ibin--) {
    cumul_hi += counts[ibin];
    if (cumul_hi >= allowed_max_cumul) {
      double y = ibin-nrap+1; // +1 here is the rapidity bin width
      if (y < _maxrap) _maxrap = y;
      break;
    }
  }
  assert(ibin >= 0); // internal consistency check that you found a bin

  // ibin_hi is the index of the rightmost bin that should be considered
  int ibin_hi = ibin;

  // consistency check 
  assert(ibin_hi >= ibin_lo); 

  // now work out cumul2
  if (ibin_hi == ibin_lo) {
    // if there is a single bin (potentially including overflows
    // from both sides), cumul2 is the square of the total contents
    // of that bin, which we obtain from cumul_lo and cumul_hi minus
    // the double counting of part that is contained in both
    // (putting double 
    _cumul2 = pow(double(cumul_lo + cumul_hi - counts[ibin_hi]), 2);
  } else {
    // otherwise we have a straightforward sum of squares of bin
    // contents
    _cumul2 += cumul_hi*cumul_hi;

    // now get the rest of the squared bin contents
    for (ibin = ibin_lo+1; ibin < ibin_hi; ibin++) {
      _cumul2 += counts[ibin]*counts[ibin];
    }
  }

}


FASTJET_END_NAMESPACE
