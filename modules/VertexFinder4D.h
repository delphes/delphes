#ifndef VertexFinder4D_h
#define VertexFinder4D_h

/** \class VertexFinder4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors M. Selvaggi
 *
 */



#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

class TObjArray;
class Candidate;
class TVector3;

class VertexFinder4D: public DelphesModule
{
public:

  VertexFinder4D();
  ~VertexFinder4D();

  void Init();
  void Process();
  void Finish();

  struct track_t{
  double z;              // z-coordinate at point of closest approach to the beamline
  double t;              // t-coordinate at point of closest approach to the beamline
  double dz2;            // square of the error of z(pca)
  double dtz;            // covariance of z-t
  double dt2;            // square of the error of t(pca)
  //const reco::TransientTrack* tt;  // a pointer to the Transient Track
  Candidate* tt;  // a pointer to the Candidate Track
  double Z;              // Z[i]   for DA clustering
  double pi;             // track weight
};


struct vertex_t{
  double z;    //           z coordinate
  double t;    //           t coordinate
  double pk;   //           vertex weight for "constrained" clustering
  // --- temporary numbers, used during update
  double ei;
  double sw;
  double swz;
  double swt;
  double se;
  // ---for Tc
  double swE;
  double Tc;
};

  

  void clusterize(const TObjArray & tracks, TObjArray & clusters);

  std::vector< Candidate >
    vertices(const TObjArray & tracks, const int verbosity=0);

  std::vector<track_t> fill(const TObjArray & tracks)const;

  bool split( double beta,
             std::vector<track_t> & tks,
             std::vector<vertex_t> & y,
             double threshold ) const;

  double update( double beta,
                std::vector<track_t> & tks,
                std::vector<vertex_t> & y ) const;

  double update(double beta,
               std::vector<track_t> & tks,
               std::vector<vertex_t> & y,
               double & )const;

  void dump(const double beta, const std::vector<vertex_t> & y, const std::vector<track_t> & tks, const int verbosity=0) const;
  bool merge(std::vector<vertex_t> &,int ) const;
  bool merge(std::vector<vertex_t> &,double & ) const;
  bool purge(std::vector<vertex_t> &, std::vector<track_t> & , double &, const double ) const;

  void splitAll( std::vector<vertex_t> & y ) const;

  double beta0(const double betamax,
	       std::vector<track_t> & tks,
	       std::vector<vertex_t> & y )const;

  double Eik(const track_t & t, const vertex_t & k)const;


private:

  Double_t fSigma;
  Double_t fMinPT;
  Double_t fMaxEta;
  Double_t fSeedMinPT;
  Int_t fMinNDF;
  Int_t fGrowSeeds;
  Bool_t fVerbose;

  bool useTc_;
  float vertexSize_;
  int maxIterations_;
  double coolingFactor_;
  float betamax_;
  float betastop_;
  double dzCutOff_;
  double d0CutOff_;
  double dtCutOff_; // for when the beamspot has time

  TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fOutputArray;
  TObjArray *fVertexOutputArray;

  ClassDef(VertexFinder4D, 1)
};

#endif
