#ifndef ALGOOBJ_HH
#define ALGOOBJ_HH
#include <vector>

class AlgoSubObj
{
public:
  int    metricId;
  bool   useCharged;
  bool   applyLowPUCorr;
  int    combId;
  double coneSize;
  double rmsPtMin;
  double rmsScaleFactor;
};
class AlgoObj 
{
public:
      AlgoObj():
	etaMin(0), etaMax(0), 
	ptMin(0),minNeutralPt(0),minNeutralPtSlope(0),
	rmsEtaSF(0),medEtaSF(0),etaMaxExtrap(0)
    {}
    ~AlgoObj(){}
  
    float etaMin;
    float etaMax;
    float ptMin;
    double minNeutralPt;
    double minNeutralPtSlope;
    double rmsEtaSF;
    double medEtaSF;
    double etaMaxExtrap;  
    std::vector<AlgoSubObj> subAlgos; 
};
#endif
