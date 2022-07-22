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

#ifndef UnstablePropagator_h
#define UnstablePropagator_h

/** \class UnstablePropagator
 *
 *  Propagates charged unstable particles in magnetic field
 *  and updates coordinates of its daughters iteratively
 *
 *  \author M. Selvaggi - CERN
 *
 */

#include "classes/DelphesModule.h"

class TClonesArray;
class TIterator;
class TLorentzVector;
class Candidate;

class UnstablePropagator: public DelphesModule
{
public:
  UnstablePropagator();
  ~UnstablePropagator();

  void Init();
  void Process();
  void Finish();

private:
  Double_t fRadius, fRadius2, fRadiusMax, fHalfLength, fHalfLengthMax;
  Double_t fBz;
  Double_t fLmin; // minimum

  Bool_t fDebug;
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  std::vector < Int_t > DaughterIndices(Candidate *candidate);
  void PrintPart(TString prefix, Candidate *candidate);
  Double_t FlightDistance(Candidate *mother, Candidate *daughter);
  Int_t Index(Candidate *candidate);
  void ComputeChainFlightDistances(TString prefix, Candidate *candidate);
  void PropagateAndUpdateChain(TString prefix, Candidate *candidate);
  TLorentzVector PropagatedPosition(Candidate *candidate);

  ClassDef(UnstablePropagator, 1)
};

#endif
