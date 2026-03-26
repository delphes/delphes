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

/** \class DelphesReader
 *
 *  Base object to read event content file
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve, L. Forthomme - AGH, Krakow
 *
 */

#include <ExRootAnalysis/ExRootProgressBar.h>

#include "classes/DelphesFactory.h"
#include "classes/DelphesReader.h"

DelphesReader::~DelphesReader()
{
  if(fInputFile)
  {
    fseek(fInputFile, 0L, SEEK_END);
    fProgressBar->Update(ftello(fInputFile), fEventCounter, true);
    fProgressBar->Finish();

    if(fInputFile != stdin) fclose(fInputFile);
  }
}

//---------------------------------------------------------------------------

void DelphesReader::SetFactory(DelphesFactory *factory)
{
  fFactory = factory;
  fAllParticleOutputArray = GetFactory()->Book<std::vector<Candidate *> >("Delphes/allParticles", false);
  fStableParticleOutputArray = GetFactory()->Book<std::vector<Candidate *> >("Delphes/stableParticles", false);
  fPartonOutputArray = GetFactory()->Book<std::vector<Candidate *> >("Delphes/partons", false);
}

//---------------------------------------------------------------------------

DelphesFactory *DelphesReader::GetFactory() const
{
  if(!fFactory)
    throw std::runtime_error("Failed to retrieve the objects factory for the reader module.");
  return fFactory;
}

//---------------------------------------------------------------------------

void DelphesReader::SetSkipEvents(long long skipEvents)
{
  if(skipEvents < 0)
    throw std::runtime_error("SkipEvents must be zero or positive");

  fSkipEvents = skipEvents;
}

//---------------------------------------------------------------------------

void DelphesReader::SetMaxEvents(long long maxEvents)
{
  if(maxEvents < 0)
    throw std::runtime_error("MaxEvents must be zero or positive");

  fMaxEvents = maxEvents;
}

//---------------------------------------------------------------------------

bool DelphesReader::ReadEvent()
{
  if(fMaxEvents > 0 && fEventCounter - fSkipEvents >= fMaxEvents)
    return false;
  // initialise the collections if not already done

  fReadStopWatch.Start();
  Clear();
  fAllParticleOutputArray->clear();
  fStableParticleOutputArray->clear();
  fPartonOutputArray->clear();
  while(ReadBlock())
  {
    if(EventReady())
    {
      ++fEventCounter;
      if(fEventCounter > fSkipEvents)
      {
        if(fInputFile)
          fProgressBar->Update(ftello(fInputFile), fEventCounter);
        fReadStopWatch.Stop();
        return true;
      }
      else if(fSkipEvents > 0)
      {
        fAllParticleOutputArray->clear();
        fStableParticleOutputArray->clear();
        fPartonOutputArray->clear();
      }
      continue;
    }
  }
  return false;
}

//---------------------------------------------------------------------------
