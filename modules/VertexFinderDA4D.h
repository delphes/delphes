#ifndef VertexFinderDA4D_h
#define VertexFinderDA4D_h

/** \class VertexFinderDA4D
 *
 *  Cluster vertices from tracks using deterministic annealing and timing information
 *
 *  \authors M. Selvaggi, L. Gray
 *
 */


#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;
class TIterator;
class Candidate;

class VertexFinderDA4D: public DelphesModule
{
public:

  VertexFinderDA4D();
  ~VertexFinderDA4D();

  void Init();
  void Process();
  void Finish();

  void clusterize(const TObjArray &tracks, TObjArray &clusters);
  std::vector< Candidate* > vertices();

private:

  Bool_t fVerbose;
  Double_t fMinPT;

  Float_t fVertexSpaceSize;
  Float_t fVertexTimeSize;
  Bool_t fUseTc;
  Float_t fBetaMax;
  Float_t fBetaStop;
  Double_t fCoolingFactor;
  Int_t fMaxIterations;
  Double_t fDzCutOff;
  Double_t fD0CutOff;
  Double_t fDtCutOff; // for when the beamspot has time

  TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fOutputArray;
  TObjArray *fVertexOutputArray;

  ClassDef(VertexFinderDA4D, 1)
};

#endif
