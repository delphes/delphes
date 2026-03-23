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

/** \class RecoPuFilter
 *
 *  Removes particles with RecoPU flag = true.
 *  Input collection needs to pass by TrackPileUpSubtractor first)
 *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

class RecoPuFilter: public DelphesModule
{
public:
  using DelphesModule::DelphesModule;

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "filteredParticles"));
  }
  void Process() override
  {
    fOutputArray->clear();
    for(Candidate *const &candidate : *fInputArray)
    {
      if(candidate->IsRecoPU) continue;
      fOutputArray->emplace_back(candidate);
    }
  }

private:
  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

REGISTER_MODULE("RecoPuFilter", RecoPuFilter);
