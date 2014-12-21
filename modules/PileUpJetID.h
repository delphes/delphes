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

#ifndef PileUpJetID_h
#define PileUpJetID_h

/** \class PileUpJetID
 *
 *  CMS PileUp Jet ID Variables, based on http://cds.cern.ch/record/1581583
 *
 *  \author S. Zenz, December 2013
 *
 */


#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;

class PileUpJetID: public DelphesModule
{
public:

  PileUpJetID();
  ~PileUpJetID();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;
  Double_t fParameterR;

  // If set to true, may have weird results for PFCHS
  // If set to false, uses everything within dR < fParameterR even if in other jets &c.
  // Results should be very similar for PF
  Int_t fUseConstituents;

  Bool_t fAverageEachTower;

  TIterator *fItJetInputArray; //!

  const TObjArray *fJetInputArray; //!

  const TObjArray *fTrackInputArray; //!
  const TObjArray *fNeutralInputArray; //!

  TIterator *fItTrackInputArray; //!
  TIterator *fItNeutralInputArray; //!

  TObjArray *fOutputArray; //!

  TIterator *fItVertexInputArray; //!
  const TObjArray *fVertexInputArray; //!

  Double_t fZVertexResolution;

  ClassDef(PileUpJetID, 1)
};

#endif

