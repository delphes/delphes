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

#ifndef DelphesEventDisplay_h
#define DelphesEventDisplay_h

#include <vector>
#include <iostream>
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "display/DelphesDisplay.h"
#include "display/Delphes3DGeometry.h"
#include "display/DelphesHtmlSummary.h"
#include "display/DelphesPlotSummary.h"
#include "TSystem.h"
#include "TChain.h"
#include "TAxis.h"
#include "TGHtml.h"
#include "TClonesArray.h"
#include "TGStatusBar.h"
#include "TGNumberEntry.h"
#include "TGProgressBar.h"
#include <RQ_OBJECT.h>



/*
 * assembly.C: sauvegarde as shape-extract -> implement in the geometry class (read/write)
 * histobrowser.C: intégration d'histogrammes dans le display (on pourrait avoir Pt, eta, phi pour les principales collections)
 * also from alice_esd: summary html table
 * 
 */

class DelphesEventDisplay
{
    RQ_OBJECT("DelphesEventDisplay")
  public:
    DelphesEventDisplay();
    DelphesEventDisplay(const char *configFile, const char *inputFile, Delphes3DGeometry& det3D);
    ~DelphesEventDisplay();
    void EventChanged(Int_t); // *SIGNAL*

  private:
    void update_html_summary();
    void make_gui();
    void load_event();
    void readConfig(const char *configFile, std::vector<DelphesBranchBase*>& elements);

    // Configuration and global variables.
    Int_t event_id_;
    Int_t event_id_tmp_;
    ExRootTreeReader *treeReader_;
    Double_t tkRadius_, totRadius_, tkHalfLength_, muHalfLength_, bz_;
    TAxis *etaAxis_, *phiAxis_;
    TChain* chain_;
    std::vector<DelphesBranchBase*> elements_;
    DelphesDisplay *delphesDisplay_;
    DelphesHtmlSummary *htmlSummary_;
    TGHtml *gHtml_;
    DelphesPlotSummary *plotSummary_;
    TGStatusBar* fStatusBar_;
    
    // gui controls
  public:
     void Fwd() {  
        if (event_id_ < treeReader_->GetEntries() - 2) {
           EventChanged(event_id_+1);
        } else {
           printf("Already at last event.\n");
        }
     }

     void Bck() {
        if (event_id_ > 0) {
           EventChanged(event_id_-1);
        } else {
           printf("Already at first event.\n");
        }
     }

    void PreSetEv(char* ev) {
      event_id_tmp_ = Int_t(atoi(ev));
    }

    void GoTo() {
      if (event_id_tmp_>=0 && event_id_tmp_ < treeReader_->GetEntries()-1) {
        EventChanged(event_id_tmp_);
      } else {
        printf("Error: no such event.\n");
      }
    }

    void InitSummaryPlots() {
      plotSummary_->FillSample(treeReader_, event_id_);
      plotSummary_->FillEvent();
      plotSummary_->Draw();
    }

    void DisplayProgress(Int_t p) { 
         fStatusBar_->SetText(Form("Processing... %d %%",p), 1);
         gSystem->ProcessEvents();
    }
};

#endif //DelphesEventDisplay_h
