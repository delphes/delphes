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

/** \class Clone
 *
 *  Clone candidate array
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

#include <TObjArray.h>

class Cloner: public DelphesModule
{
public:
  Cloner() = default;

  void Init() override
  {
    // import input array(s)
    fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
    fItInputArray.reset(fInputArray->MakeIterator());

    // create output array(s)
    fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  }
  void Process() override
  {
    Candidate *candidate;

    // loop over all input candidates
    fItInputArray->Reset();
    while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
    {
      candidate = static_cast<Candidate *>(candidate->Clone());
      fOutputArray->Add(candidate);
    }
  }

private:
  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

REGISTER_MODULE("Cloner", Cloner);
