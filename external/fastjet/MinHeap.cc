//FJSTARTHEADER
// $Id: MinHeap.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/internal/MinHeap.hh"
#include<iostream>
#include<cmath>
#include<limits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

//----------------------------------------------------------------------
/// construct the MinHeap; structure will be as follows:
///   . _heap[0].minloc points to globally smallest entry
///     _heap[1].minloc points to smallest entry in one half of heap
///     _heap[2].minloc points to smallest entry in other half of heap
///
///   . for _heap[i], the "parent" is to be found at (i-1)/2
void MinHeap::initialise(const std::vector<double> & values){
  
  // fill the high-range of the heap with the largest possible value
  // (minloc of each entry is itself)
  for (unsigned i = values.size(); i < _heap.size(); i++) {
    _heap[i].value = std::numeric_limits<double>::max();
    _heap[i].minloc = &(_heap[i]);
  }

  // fill the rest of the heap with the actual values
  // (minloc of each entry is itself)
  for (unsigned i = 0; i < values.size(); i++) {
    _heap[i].value = values[i];
    _heap[i].minloc = &(_heap[i]);
  }
  
  // now adjust the minlocs so that everything is OK...
  for (unsigned i = _heap.size()-1; i > 0; i--) {
    ValueLoc * parent = &(_heap[(i-1)/2]);
    ValueLoc * here   = &(_heap[i]);
    if (here->minloc->value < parent->minloc->value) {
      parent->minloc = here->minloc;
    }
  }
  //cout << minloc() << " "<<sqrt(minval())<<endl;
  //cout << sqrt(_heap[47].value)<<endl;
  //cout << sqrt(_heap[48].value)<<endl;
  //cout << sqrt(_heap[25].value)<<endl;
}


//----------------------------------------------------------------------
void MinHeap::update(unsigned int loc, double new_value) {
  

  assert(loc < _heap.size());
  ValueLoc * start = &(_heap[loc]);

  // if the minloc is somewhere below us and our value is no smaller
  // than the previous value, we can't possibly change the minloc
  if (start->minloc != start && !(new_value < start->minloc->value)) {
    start->value = new_value;
    //std::cout << "                     had easy exit\n";
    return;
  }

  // update the value and put a temporary location
  start->value = new_value;
  start->minloc = start;
  // warn that a change has been made at this place
  bool change_made = true;
  // to make sure we don't go off edge...
  ValueLoc * heap_end = (&(_heap[0])) + _heap.size();

  // now work our way up the heap
  while(change_made) {
    ValueLoc * here = &(_heap[loc]);
    change_made     = false;

    // if we were pointing to start, then we must re-initialise things
    if (here->minloc == start) {
      here->minloc = here; change_made = true;
    }

    // now compare current location to children (at 2*loc+1, 2*loc+2)
    ValueLoc * child = &(_heap[2*loc+1]);
    if (child < heap_end && child->minloc->value < here->minloc->value ) {
      here->minloc = child->minloc;
      change_made = true;}
    child++;
    if (child < heap_end && child->minloc->value < here->minloc->value ) {
      here->minloc = child->minloc;
      change_made = true;}
    
    // then we move up (loc ->(loc-1)/2) if there's anywhere to go 
    if (loc == 0) {break;}
    loc = (loc-1)/2;
  }

}

FASTJET_END_NAMESPACE

