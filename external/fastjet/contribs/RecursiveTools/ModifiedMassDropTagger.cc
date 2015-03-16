// $Id: ModifiedMassDropTagger.cc 683 2014-06-13 14:38:38Z gsoyez $
//
// Copyright (c) 2014-, Gavin P. Salam
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

#include "ModifiedMassDropTagger.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include <algorithm> 
#include <cstdlib> 

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
string ModifiedMassDropTagger::symmetry_cut_description() const {
  ostringstream ostr;
  ostr << _symmetry_cut << " [ModifiedMassDropTagger]";
  return ostr.str();
}

} // namespace contrib

FASTJET_END_NAMESPACE
