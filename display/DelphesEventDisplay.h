/*
  * Delphes: a framework for fast simulation of a generic collider experiment
  * Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
 *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
 *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DelphesEventDisplay_h
#define DelphesEventDisplay_h

#include <vector>

#include "Rtypes.h"
#include "RQ_OBJECT.h"

class TAxis;
class TChain;
class TGHtml;
class TGStatusBar;
class DelphesDisplay;
class Delphes3DGeometry;
class DelphesBranchBase;
class DelphesHtmlSummary;
class DelphesPlotSummary;
class ExRootTreeReader;

/*
  *assembly.C: sauvegarde as shape-extract -> implement in the geometry class (read/write)
  *histobrowser.C: intégration d'histogrammes dans le display (on pourrait avoir Pt, eta, phi pour les principales collections)
  *also from alice_esd: summary html table
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
    void readConfig(const char *configFile, std::vector<DelphesBranchBase *>& elements);

    // Configuration and global variables.
    Int_t event_id_;
    Int_t event_id_tmp_;
    ExRootTreeReader *treeReader_;
    Double_t tkRadius_, totRadius_, tkHalfLength_, muHalfLength_, bz_;
    TAxis *etaAxis_, *phiAxis_;
    TChain *chain_;
    std::vector<DelphesBranchBase *> elements_;
    DelphesDisplay *delphesDisplay_;
    DelphesHtmlSummary *htmlSummary_;
    TGHtml *gHtml_;
    DelphesPlotSummary *plotSummary_;
    TGStatusBar *fStatusBar_;
    
    // gui controls
  public:
     void Fwd();

     void Bck();

    void PreSetEv(char *ev);

    void GoTo();

    void InitSummaryPlots();

    void DisplayProgress(Int_t p);
};

#endif //DelphesEventDisplay_h
