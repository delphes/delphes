/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2026  Universite catholique de Louvain (UCL), Belgium
 *                           AGH University of Krakow, Poland
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class PyDelpheEvent
 *
 *  Event I/O helper for Python bindings
 *
 *  \author L. Forthomme - AGH, Kraków
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

#include "PyDelphesEvent.h"

std::vector<Candidate *> *PyDelphesEvent::Get(std::string_view collName) const
{
  if(CandidatesCollection collObj = fFactory.Attach<std::vector<Candidate *> >(collName); collObj)
    return collObj.get();
  return nullptr;
}

//------------------------------------------------------------------------------
