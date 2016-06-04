#ifndef PUPPICLEANCONTAINER_HH
#define PUPPICLEANCONTAINER_HH


#include "PUPPI/RecoObj.hh"
#include "PUPPI/puppiParticle.hh"
#include "PUPPI/puppiAlgoBin.hh"

#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include <algorithm>

using namespace std;

//......................

class puppiCleanContainer{

 public:

  // basic constructor which takes input particles as RecoObj, the tracker eta extension and other two boolean info
  puppiCleanContainer(std::vector<RecoObj> inParticles,    // incoming particles of the event
                      std::vector<puppiAlgoBin> puppiAlgo, // vector with the definition of the puppi algorithm in different eta region (one for each eta) 
                      float minPuppiWeight  = 0.01,        // min puppi weight cut
                      bool  useExp = false                 // useDz vertex probability
  ); 

  ~puppiCleanContainer(); 

  // ----- get methods

  // get all the PF particles
  std::vector<fastjet::PseudoJet> pfParticles()  { return  fPFParticles_; }    
  // get all the PF charged from PV
  std::vector<fastjet::PseudoJet> pvParticles()  { return  fChargedPV_; }        
  // get all the PF charged from PU
  std::vector<fastjet::PseudoJet> puParticles()  { return  fChargedNoPV_; }    
  // get CHS particle collection
  std::vector<fastjet::PseudoJet> pfchsParticles(){ return fPFchsParticles_; }    
  // get puppi weight for all particles 
  std::vector<float> getPuppiWeights() { return fPuppiWeights_; };

  // process puppi
  std::vector<fastjet::PseudoJet> puppiEvent();
 
 protected:

   void    getRMSAvg(const int &, std::vector<fastjet::PseudoJet> &, std::vector<fastjet::PseudoJet> &);        
   float   goodVar  (const fastjet::PseudoJet &, const std::vector<fastjet::PseudoJet> &, const int &, const float &);    
   void    computeMedRMS(const int &);  
   float   compute(const float &, const std::vector<puppiParticle> &, const std::vector<puppiAlgoBin> &, const std::vector<int> &);

   // some get functions
   float getNeutralPtCut(const float&, const float&, const int&);
   std::vector<int> getPuppiId(const float &, const float &, const std::vector<puppiAlgoBin> &);
   bool  isGoodPuppiId(const float &, const float &, const puppiAlgoBin &);
   float getChi2FromdZ(float);
   // other functions
   float  var_within_R(const int &, const vector<fastjet::PseudoJet> &, const fastjet::PseudoJet &, const float &);
   float  pt_within_R(const std::vector<fastjet::PseudoJet> &, const fastjet::PseudoJet &, const float &);
   fastjet::PseudoJet flow_within_R(const vector<fastjet::PseudoJet> &, const fastjet::PseudoJet &, const float &);
   
  
 private:    
 
  std::vector<RecoObj>            fRecoParticles_;
  std::vector<fastjet::PseudoJet> fPFParticles_;
  std::vector<fastjet::PseudoJet> fPFchsParticles_;    
  std::vector<fastjet::PseudoJet> fChargedPV_;
  std::vector<fastjet::PseudoJet> fChargedNoPV_;

  std::vector<puppiAlgoBin> puppiAlgo_;
  std::vector<float> fPuppiWeights_;

  float  fMinPuppiWeight_;
  float  fPVFrac_;

  int    fNPV_;  
  bool   fUseExp_ ;
    
};

#endif
