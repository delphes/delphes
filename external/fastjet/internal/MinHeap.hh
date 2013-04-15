//STARTHEADER
// $Id: MinHeap.hh 2577 2011-09-13 15:11:38Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

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
    _heap(max_size) {_initialise(values);};

  /// constructor in which the the maximum size is the size of the values array
  MinHeap (const std::vector<double> & values) :
    _heap(values.size()) {_initialise(values);};
  
  /// return the location of the minimal value on the heap
  inline unsigned int minloc() const {
    return (_heap[0].minloc) - &(_heap[0]);};
  
  /// return the minimal value on the heap
  inline double       minval() const {return _heap[0].minloc->value;};

  inline double operator[](int i) const {return _heap[i].value;};

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

  void _initialise(const std::vector<double> & values);


};


FASTJET_END_NAMESPACE

#endif // __FASTJET_MINHEAP__HH__
