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

#ifndef DelphesTreeWriter_h
#define DelphesTreeWriter_h

/** \class DelphesTreeWriter
 *
 *  ExRootTreeWriter with event-level filtering and normalisation.
 *
 *  \author Sihyun Jeon (Boston University)
 *
 */

#include "ExRootAnalysis/ExRootTreeWriter.h"

class TClonesArray;
class TFile;

class DelphesTreeWriter: public ExRootTreeWriter
{
public:
  DelphesTreeWriter(TFile *file = 0, const char *treeName = "Analysis");

  void SetEventVeto() { fEventVeto = kTRUE; }

  void Fill();
  void Write();

private:
  TFile *fOutputFile;

  Bool_t fEventVeto;
  Bool_t fEventArrayResolved;
  TClonesArray *fEventArray;
  Long64_t fInputEvents;
  Long64_t fWrittenEvents;
  Double_t fSumWeightsInput;
  Double_t fSumWeightsWritten;
};

#endif /* DelphesTreeWriter */
