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

#ifndef HighMassVertexRecover_h
#define HighMassVertexRecover_h

/** \class HighMassVertexRecover
 *
 *  Performs transverse time smearing.
 *
 *  \author Olmo Cerri, Caltech
 *
 */

 #include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"

class TIterator;
class TObjArray;

class HighMassVertexRecover: public DelphesModule
{
public:

  HighMassVertexRecover();
  ~HighMassVertexRecover();

  void Init();
  void Process();
  void Finish();

  std::pair<Double_t, Double_t> ComputeCATime(Candidate * tk, Double_t m);

private:

  UInt_t fVerbose;
  Double_t fSigmaCompatibility;
  std::vector<Double_t> fMassList;

  TIterator *fItTrackInputArray;
  const TObjArray *fTrackInputArray;

  TIterator *fItVertexInputArray;
  const TObjArray *fVertexInputArray;

  TObjArray *fTrackOutputArray;
  TObjArray *fVertexOutputArray;
  TIterator *fItVertexOutputArray;

  ClassDef(HighMassVertexRecover, 1)
};

#endif
