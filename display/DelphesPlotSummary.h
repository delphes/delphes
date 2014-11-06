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

