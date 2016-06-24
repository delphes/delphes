#include "PuppiAlgo.hh"
#include "RecoObj2.hh"
#include "fastjet/PseudoJet.hh"
#include <vector>

using namespace std;

class PuppiContainer{
public:
    //PuppiContainer(const edm::ParameterSet &iConfig);
    PuppiContainer(const std::string &iConfig);
    PuppiContainer(bool iApplyCHS, bool iUseExp,double iPuppiWeightCut,std::vector<AlgoObj> &iAlgos);
    ~PuppiContainer(); 
    void initialize(const std::vector<RecoObj> &iRecoObjects);
    std::vector<fastjet::PseudoJet> pfParticles(){ return fPFParticles; }    
    std::vector<fastjet::PseudoJet> pvParticles(){ return fChargedPV; }        
    const std::vector<double> puppiWeights();
    std::vector<fastjet::PseudoJet> puppiParticles() { return fPupParticles;}

protected:
    double  goodVar      (fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt,double iRCone);
    void    getRMSAvg    (int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,std::vector<fastjet::PseudoJet> &iChargeParticles);
    double  getChi2FromdZ(double iDZ);
    int     getPuppiId   (const float &iPt,const float &iEta);
    double  var_within_R (int iId, const std::vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, double R);  
    
    std::vector<RecoObj>  fRecoParticles;
    std::vector<fastjet::PseudoJet> fPFParticles;
    std::vector<fastjet::PseudoJet> fChargedPV;
    std::vector<fastjet::PseudoJet> fPupParticles;
    std::vector<double>    fWeights;
    std::vector<double>    fVals;
    bool   fApplyCHS;
    bool   fUseExp;
    double fNeutralMinPt;
    double fNeutralSlope;
    double fPuppiWeightCut;
    int    fNAlgos;
    int    fNPV;
    double fPVFrac;
    std::vector<PuppiAlgo> fPuppiAlgo;
};


