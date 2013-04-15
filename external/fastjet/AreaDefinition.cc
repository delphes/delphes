//STARTHEADER
// $Id$
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

#include "fastjet/AreaDefinition.hh"
#include<sstream>
#include<string>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

string VoronoiAreaSpec::description() const {
  ostringstream ostr;
  ostr << "Voronoi area with effective_Rfact = " << effective_Rfact() ;
  return ostr.str();
}


//----------------------------------------------------------------------
///  return info about the type of area being used by this defn
string AreaDefinition::description() const {
  ostringstream ostr;

  switch(area_type()) {
  case active_area:
    ostr << "Active area (hidden ghosts) with " ;
    ostr << ghost_spec().description();
    break;
  case active_area_explicit_ghosts:
    ostr << "Active area (explicit ghosts) with " ;
    ostr << ghost_spec().description();
    break;
  case one_ghost_passive_area:
    ostr << "Passive area (one ghost at a time) with " ;
    ostr << ghost_spec().description();
    break;
  case passive_area:
    ostr << "Passive area (optimal alg. based on jet.def.), where relevant with " ;
    ostr << ghost_spec().description()  ;
    break;
  case voronoi_area:
    ostr << voronoi_spec().description();
    break;
  default:
    cerr << "Error: unrecognized area_type in AreaDefinition::description():" 
         << area_type() << endl;
    exit(-1);
  }
  return ostr.str();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
