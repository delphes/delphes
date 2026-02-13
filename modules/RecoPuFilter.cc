/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
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

/*  \class RecoPuFilter
 *
 *  Removes particles with RecoPU flag = true.
 *  Input collection needs to pass by TrackPileUpSubtractor first)
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/RecoPuFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"

//------------------------------------------------------------------------------

void RecoPuFilter::Init()
{
  // import input array
  GetFactory()->EventModel()->Attach(GetString("InputArray", "Delphes/allParticles"), fInputArray);
  // create output array
  GetFactory()->EventModel()->Book(fOutputArray, GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void RecoPuFilter::Finish()
{
}

//------------------------------------------------------------------------------

void RecoPuFilter::Process()
{
  for(const auto &candidate : *fInputArray)
  {
    if(candidate.IsRecoPU) continue;
    fOutputArray->emplace_back(candidate);
  }
}
