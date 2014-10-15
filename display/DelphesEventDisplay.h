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
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "display/DelphesDisplay.h"
#include "display/Delphes3DGeometry.h"
#include "TChain.h"
#include "TClonesArray.h"

/*
 * assembly.C: sauvegarde as shape-extract -> implement in the geometry class (read/write)
 * histobrowser.C: intégration d'histogrammes dans le display (on pourrait avoir Pt, eta, phi pour les principales collections)
 * also from alice_esd: summary html table
 * 
 */

class DelphesEventDisplay
{
  public:
    DelphesEventDisplay();
    DelphesEventDisplay(const char *configFile, const char *inputFile, Delphes3DGeometry& det3D);
    ~DelphesEventDisplay();
    Int_t event_id;
    ExRootTreeReader *gTreeReader;
    void load_event();

  private:
    void make_gui();
    void delphes_read();
    void delphes_read_towers(TClonesArray* data, DelphesBranchBase* element);
    void delphes_read_tracks(TClonesArray* data, DelphesBranchBase* element);
    void delphes_read_jets(TClonesArray* data, DelphesBranchBase* element);
    void delphes_read_vectors(TClonesArray* data, DelphesBranchBase* element);
    void readConfig(const char *configFile, Delphes3DGeometry& det3D, std::vector<DelphesBranchBase*>& elements, std::vector<TClonesArray*>& arrays);

    // Configuration and global variables.
    Double_t tkRadius_, totRadius_, tkHalfLength_, bz_;
    TChain* chain_;
    std::vector<DelphesBranchBase*> gElements;
    std::vector<TClonesArray*> gArrays;
    DelphesDisplay *gDelphesDisplay;

    // EvNavHandler class is needed to connect GUI signals.
    class EvNavHandler
    {     
      public:
      
         EvNavHandler(DelphesEventDisplay* display):display_(display) {}

         ~EvNavHandler() {}

         void Fwd() {  
            if (display_->event_id < display_->gTreeReader->GetEntries() - 1) {
               ++(display_->event_id);
               display_->load_event();
            } else {
               printf("Already at last event.\n");
            }
         }

         void Bck() {
            if (display_->event_id > 0) {
               --(display_->event_id);
               display_->load_event();
            } else {
               printf("Already at first event.\n");
            }
         }
       private:
         DelphesEventDisplay* display_;
    };
};

#endif //DelphesEventDisplay_h
