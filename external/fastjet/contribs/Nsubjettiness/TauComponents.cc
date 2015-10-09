//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessDefinition.cc 704 2014-07-07 14:30:43Z jthaler $
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

#include "TauComponents.hh"
#include "MeasureDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

// This constructor takes input vector and double and calculates all necessary tau components
TauComponents::TauComponents(TauMode tau_mode,
              const std::vector<double> & jet_pieces_numerator,
              double beam_piece_numerator,
              double denominator,
              const std::vector<PseudoJet> & jets,
              const std::vector<PseudoJet> & axes
              )
:  _tau_mode(tau_mode),
_jet_pieces_numerator(jet_pieces_numerator),
_beam_piece_numerator(beam_piece_numerator),
_denominator(denominator),
_jets(jets),
_axes(axes)
{
   
   if (!has_denominator()) assert(_denominator == 1.0); //make sure no effect from _denominator if _has_denominator is false
   if (!has_beam()) assert (_beam_piece_numerator == 0.0); //make sure no effect from _beam_piece_numerator if _has_beam is false
   
   // Put the pieces together
   _numerator = _beam_piece_numerator;
   _jet_pieces.resize(_jet_pieces_numerator.size(),0.0);
   for (unsigned j = 0; j < _jet_pieces_numerator.size(); j++) {
      _jet_pieces[j] = _jet_pieces_numerator[j]/_denominator;
      _numerator += _jet_pieces_numerator[j];
      
      // Add structural information to jets
      StructureType * structure = new StructureType(_jets[j]);
      structure->_tau_piece = _jet_pieces[j];
      _jets[j].set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
   }
   
   _beam_piece = _beam_piece_numerator/_denominator;
   _tau = _numerator/_denominator;
   
   // Add total_jet with structural information
   _total_jet = join(_jets);
   StructureType * total_structure = new StructureType(_total_jet);
   total_structure->_tau_piece = _tau;
   _total_jet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(total_structure));
}
   
   
   
// test for denominator/beams
bool TauComponents::has_denominator() const {
   return (_tau_mode == NORMALIZED_JET_SHAPE
           || _tau_mode == NORMALIZED_EVENT_SHAPE);
}

bool TauComponents::has_beam() const {
   return (_tau_mode == UNNORMALIZED_EVENT_SHAPE
           || _tau_mode == NORMALIZED_EVENT_SHAPE);
}

} // namespace contrib

FASTJET_END_NAMESPACE
