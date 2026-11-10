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

/** \class DelphesTreeWriter
 *
 *  ExRootTreeWriter with event-level filtering and normalisation.
 *
 *  \author Sihyun Jeon (Boston University)
 *
 */

#include "classes/DelphesTreeWriter.h"

#include "classes/DelphesClasses.h"

#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

//------------------------------------------------------------------------------

DelphesTreeWriter::DelphesTreeWriter(TFile *file, const char *treeName) :
  ExRootTreeWriter(file, treeName),
  fOutputFile(file),
  fEventVeto(kFALSE), fEventArrayResolved(kFALSE), fEventArray(0),
  fInputEvents(0), fWrittenEvents(0),
  fSumWeightsInput(0.0), fSumWeightsWritten(0.0)
{
}

//------------------------------------------------------------------------------

void DelphesTreeWriter::Fill()
{
  Double_t weight = 1.0;
  const HepMCEvent *hepMCEvent;
  const LHEFEvent *lhefEvent;

  if(!fEventArrayResolved)
  {
    fEventArrayResolved = kTRUE;
    TBranch *branch = GetTree() ? GetTree()->GetBranch("Event") : 0;
    if(branch && branch->GetAddress()) fEventArray = *(TClonesArray **)branch->GetAddress();
  }

  if(fEventArray && fEventArray->GetEntriesFast() > 0)
  {
    hepMCEvent = dynamic_cast<HepMCEvent *>(fEventArray->At(0));
    lhefEvent = dynamic_cast<LHEFEvent *>(fEventArray->At(0));
    if(hepMCEvent)
      weight = hepMCEvent->Weight;
    else if(lhefEvent)
      weight = lhefEvent->Weight;
  }

  ++fInputEvents;
  fSumWeightsInput += weight;

  if(fEventVeto)
  {
    fEventVeto = kFALSE;
    return;
  }

  ++fWrittenEvents;
  fSumWeightsWritten += weight;

  ExRootTreeWriter::Fill();
}

//------------------------------------------------------------------------------

void DelphesTreeWriter::Write()
{
  TFile *file = fOutputFile;
  if(!file && GetTree()) file = GetTree()->GetCurrentFile();

  if(file)
  {
    TH1D eventCount("EventCount", "event filter", 4, 0.0, 4.0);
    eventCount.SetDirectory(0);
    eventCount.GetXaxis()->SetBinLabel(1, "input events");
    eventCount.GetXaxis()->SetBinLabel(2, "written events");
    eventCount.GetXaxis()->SetBinLabel(3, "sum weights input");
    eventCount.GetXaxis()->SetBinLabel(4, "sum weights written");
    eventCount.SetBinContent(1, fInputEvents);
    eventCount.SetBinContent(2, fWrittenEvents);
    eventCount.SetBinContent(3, fSumWeightsInput);
    eventCount.SetBinContent(4, fSumWeightsWritten);
    file->cd();
    eventCount.Write();
  }

  ExRootTreeWriter::Write();
}

//------------------------------------------------------------------------------
