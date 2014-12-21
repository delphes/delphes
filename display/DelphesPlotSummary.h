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

#ifndef DelphesPlotSummary_h
#define DelphesPlotSummary_h

#include <vector>
#include <map>
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEveWindow.h"
#include "DelphesBranchElement.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <RQ_OBJECT.h>


class DelphesPlotSummary
{
    RQ_OBJECT("DelphesPlotSummary")
  public:
    DelphesPlotSummary(TEveWindowTab* tab);
    virtual ~DelphesPlotSummary();
    void Init(std::vector<DelphesBranchBase*>& elements);
    void FillSample(ExRootTreeReader* treeReader, Int_t event_id);
    void FillEvent();
    void Draw();
    void Progress(Int_t); // *SIGNAL*

  private:
    TEveWindowTab* tab_;
    std::map< TString, TCanvas* >           canvases_;
    std::map< TString, std::vector<TH1F*> > histograms_;
    std::vector<DelphesBranchBase*>* elements_;
    std::map< TString, std::vector<TMarker*> > eventMarkers_;
    std::map< TString, std::vector<TH1F*> > eventProfiles_;

};

#endif // DelphesPlotSummary_h

