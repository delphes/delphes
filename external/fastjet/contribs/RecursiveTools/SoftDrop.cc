// $Id: SoftDrop.cc 686 2014-06-14 03:25:09Z jthaler $
//
// Copyright (c) 2014-, Gregory Soyez, Jesse. Thaler
// based on arXiv:1402.2657 by Andrew J. Larkoski, Simone Marzani,
// Gregory Soyez, Jesse Thaler
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

#include "SoftDrop.hh"
#include <cmath>
#include <sstream>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
// TODO:
//
// - implement reclustering (at the moment we assume it's C/A as for mMDT
//
// - what is returned if no substructure is found? [and for negativeve
//   pt, m2 or other situations where mMDT currentlty returns an empty
//   PseudoJet]
//
//   GS.: mMDT seeks substructure so it makes sense for it to return
//   an empry PseudoJet when no substructure is found (+ it is the
//   original behaviour). For SoftDrop, in grooming mode (beta>0), if
//   would make sense to return a jet with a single particle. In
//   tagging mode (beta<0), the situation is less clear. At the level
//   of the implementation, having a virtual function could work
//   (with a bit of care to cover the -ve pt or m2 cases)
//
// - Do we include Andrew and Simone in the "contrib-author" list
//   since they are on the paper? Do we include Gavin in the author's
//   list since he started this contrib?
//
//----------------------------------------------------------------------

//----------------------------------------------------------------------
double SoftDrop::symmetry_cut_fn(const PseudoJet & p1, 
                                 const PseudoJet & p2) const{
  return _symmetry_cut * pow(p1.squared_distance(p2)/_R0sqr, 0.5*_beta);
}

//----------------------------------------------------------------------
string SoftDrop::symmetry_cut_description() const {
  ostringstream ostr;
  ostr << _symmetry_cut << " (theta/" << sqrt(_R0sqr) << ")^" << _beta << " [SoftDrop]";
  return ostr.str();
}

} // namespace contrib

FASTJET_END_NAMESPACE
