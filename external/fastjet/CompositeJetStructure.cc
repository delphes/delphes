//FJSTARTHEADER
// $Id: CompositeJetStructure.cc 3433 2014-07-23 08:17:03Z salam $
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

#include <fastjet/CompositeJetStructure.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;


//-------------------------------------------------------------------------------
// \class CompositeJetStructure
// The structure for a jet made of pieces
//
// This stores the vector of the pieces that make the jet and provide
// the methods to access them
// -------------------------------------------------------------------------------

CompositeJetStructure::CompositeJetStructure(const std::vector<PseudoJet> & initial_pieces, 
					     const JetDefinition::Recombiner * recombiner)
  : _pieces(initial_pieces){

#ifndef __FJCORE__
  // deal with area support (cache the area if needed)
  //--------------------------------------------------
  // check if all the pieces have area, in which case store it
  bool has_area_local = true;
  for (vector<PseudoJet>::const_iterator pit=_pieces.begin(); pit!=_pieces.end(); pit++){
    if (!pit->has_area()){
      has_area_local = false;
      continue;
    }
  }

  if (has_area_local){
    _area_4vector_ptr = new PseudoJet();
    for (unsigned int i=0; i<_pieces.size(); i++){
      const PseudoJet & p = _pieces[i];
      if (recombiner)
	recombiner->plus_equal(*_area_4vector_ptr, p.area_4vector());
      else
	*_area_4vector_ptr += p.area_4vector();
    } 
  } else {
    _area_4vector_ptr = 0;
  }
#else
  if (recombiner){};  // ugly trick to prevent a gcc warning
  _area_4vector_ptr = 0;
#endif

}


// description
std::string CompositeJetStructure::description() const{ 
  string str = "Composite PseudoJet";
  return str; 
}



// things reimplemented from the base structure
//------------------------------------------------------------------------------
bool CompositeJetStructure::has_constituents() const{
  //for (vector<PseudoJet>::const_iterator pit=_pieces.begin(); pit!=_pieces.end(); pit++)
  //  if (!pit->has_constituents()) return false;
  //
  //return true;

  // the only case where we do not have constituents is the case where
  // there is no pieces!
  return _pieces.size()!=0;
}

std::vector<PseudoJet> CompositeJetStructure::constituents(const PseudoJet & /*jet*/) const{
  // recurse into the pieces that ahve constituents, just append the others
  // the following code automatically throws an Error if any of the
  // pieces has no constituents
  vector<PseudoJet> all_constituents;
  for (unsigned i = 0; i < _pieces.size(); i++) {
    if (_pieces[i].has_constituents()){
      vector<PseudoJet> constits = _pieces[i].constituents();
      copy(constits.begin(), constits.end(), back_inserter(all_constituents));
    } else {
      all_constituents.push_back(_pieces[i]);
    }
  }
 
  return all_constituents;
}

std::vector<PseudoJet> CompositeJetStructure::pieces(const PseudoJet & /*jet*/) const{
  return _pieces;
}


#ifndef __FJCORE__
// area-related material

// check if it has a well-defined area
bool CompositeJetStructure::has_area() const{
  return (_area_4vector_ptr != 0);
}

// return the jet (scalar) area.
double CompositeJetStructure::area(const PseudoJet & /*reference*/) const{
  if (! has_area())
    throw Error("One or more of this composite jet's pieces does not support area");

  double a=0;
  for (unsigned i = 0; i < _pieces.size(); i++)
    a += _pieces[i].area();

  return a;
}

// return the error (uncertainty) associated with the determination
// of the area of this jet.
// 
// Be conservative: return the sum of the errors
double CompositeJetStructure::area_error(const PseudoJet & /*reference*/) const{
  if (! has_area())
    throw Error("One or more of this composite jet's pieces does not support area");

  double a_err=0;
  for (unsigned i = 0; i < _pieces.size(); i++)
    a_err += _pieces[i].area_error();

  return a_err;
}

// return the jet 4-vector area.
PseudoJet CompositeJetStructure::area_4vector(const PseudoJet & /*reference*/) const{
  if (! has_area())
    throw Error("One or more of this composite jet's pieces does not support area");

  return *_area_4vector_ptr; // one is supposed to call has_area before!
}

// true if this jet is made exclusively of ghosts.
//
// In this case, it will be true if all pieces are pure ghost
bool CompositeJetStructure::is_pure_ghost(const PseudoJet & /*reference*/) const{
  for (unsigned i = 0; i < _pieces.size(); i++)
    if (! _pieces[i].is_pure_ghost()) return false;

  return true;
}

#endif  // __FJCORE__


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
