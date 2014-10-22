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
    void FillSample(ExRootTreeReader* treeReader);
    void FillEvent();

  private:
    TEveWindowTab* tab_;
    std::vector<TCanvas*> canvases;
    std::map< TString, std::vector<TH1F*> > histograms;

};

#endif // DelphesPlotSummary_h

