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

class DelphesPlotSummary
{
  public:
    DelphesPlotSummary(TEveWindowTab* tab);
    virtual ~DelphesPlotSummary();
    void Init(std::vector<DelphesBranchBase*>& elements);
    void FillSample(ExRootTreeReader* treeReader, Int_t event_id);
    void FillEvent();
    void Draw();

  private:
    TEveWindowTab* tab_;
    std::map< TString, TCanvas* >           canvases_;
    std::map< TString, std::vector<TH1F*> > histograms_;
    std::vector<DelphesBranchBase*>* elements_;
    std::map< TString, std::vector<TMarker*> > eventMarkers_;
    std::map< TString, std::vector<TH1F*> > eventProfiles_;

};

#endif // DelphesPlotSummary_h

