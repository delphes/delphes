#ifndef PUPPIPARTICLE_HH
#define PUPPIPARTICLE_HH

#include <functional>
#include <algorithm>


//............ class to keep track of particles used in the puppi code storing pt,eta,metric value, id (LV,PU) and position in the orignal vector of particles 
class puppiParticle {

 public:

  puppiParticle(){
    fPt_         = 0.;
    fEta_        = 0.;
    fPval_       = 1.;
    fParticleId_ = -1;
    fPosition_   = -1; 
  }; 

  puppiParticle(const float & pt, const float & eta, const float & pval, const int & particleId, const int & position):
    fPt_(pt),
    fEta_(eta),
    fPval_(pval),
    fParticleId_(particleId),
    fPosition_(position)
  {};

  puppiParticle(const puppiParticle & particle) {
    fPt_         = particle.fPt_;
    fEta_        = particle.fEta_;
    fPval_       = particle.fPval_;
    fParticleId_ = particle.fParticleId_;
    fPosition_   = particle.fPosition_;

  }

  virtual ~puppiParticle(){};

  // sort asaf of pt
  bool operator < (const puppiParticle & particle2) const {
    if(fPt_ < particle2.fPt_) return true;
    else return false;
  }; 


  // oper operator =
  bool operator == (const puppiParticle & particle2) const {
    if( fPt_ == particle2.fPt_     and fEta_ == particle2.fEta_ and 
        fPval_ == particle2.fPval_ and fParticleId_ == particle2.fParticleId_ and 
        fPosition_ == particle2.fPosition_) return true;
    else return false ;
  };

  float fPt_;         // pt of the candidate                                                                                                                                              
  float fEta_;        // eta of the candidate                                                                                                                                             
  float fPval_;       // metric value                                                                                                                                                     
  int   fParticleId_; // particle id means user_index                                                                                                                                     
  int   fPosition_;   // position in the original particle vector                                                                                                                        

};

class puppiValSort : public std::binary_function<int,int,bool> {
 public:

  puppiValSort(){};

  ~puppiValSort(){};

  bool operator() (const puppiParticle & x, const puppiParticle & y){
    if(x.fPval_ < y.fPval_ ) return true;
    return false;
  }
};

#endif
