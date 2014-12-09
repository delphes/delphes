//FJSTARTHEADER
// $Id: MinHeap.hh 3433 2014-07-23 08:17:03Z salam $
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

#ifndef __FASTJET_MINHEAP__HH__
#define __FASTJET_MINHEAP__HH__

#include<vector>
#include<cassert>
#include<memory>
#include<limits>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//======================================================================
/// \if internal_doc
/// @ingroup internal
/// \class MinHeap
/// A class which provides a "heap"-like structure that allows
/// access to a the minimal value of a dynamically changing set of numbers
/// \endif
class MinHeap {
public:
  /// construct a MinHeap from the vector of values, allowing future
  /// expansion to a maximum size max_size;
  MinHeap (const std::vector<double> & values, unsigned int max_size) :
    _heap(max_size) {initialise(values);}

  /// do the minimal setup for a MinHeap that can reach max_size;
  /// initialisation must be performed later with the actual values.
  MinHeap (unsigned int max_size) : _heap(max_size) {}

  /// constructor in which the the maximum size is the size of the values array
  MinHeap (const std::vector<double> & values) :
    _heap(values.size()) {initialise(values);}

  /// initialise the heap with the supplied values. Should only be called if
  /// the constructor did not supply values.
  void initialise(const std::vector<double> & values);

  /// return the location of the minimal value on the heap
  inline unsigned int minloc() const {
    return (_heap[0].minloc) - &(_heap[0]);}
  
  /// return the minimal value on the heap
  inline double       minval() const {return _heap[0].minloc->value;}

  inline double operator[](int i) const {return _heap[i].value;}

  /// remove the value at the specified location (i.e. replace it with
  /// the largest possible value).
  void remove(unsigned int loc) {
    update(loc,std::numeric_limits<double>::max());};

  /// update the value at the specified location
  void update(unsigned int, double);

private:

  struct ValueLoc{
    double value;
    ValueLoc * minloc;
  };
      
  std::vector<ValueLoc> _heap;



};


FASTJET_END_NAMESPACE

#endif // __FASTJET_MINHEAP__HH__
