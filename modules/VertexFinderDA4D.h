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
  std::vector<Candidate *> vertices();

private:
  Bool_t fVerbose{0};
  Double_t fMinPT{0.};

  Float_t fVertexSpaceSize{0.};
  Float_t fVertexTimeSize{0.};
  Bool_t fUseTc{false};
  Float_t fBetaMax{0.};
  Float_t fBetaStop{0.};
  Double_t fCoolingFactor{0.};
  Int_t fMaxIterations{0};
  Double_t fDzCutOff{0.};
  Double_t fD0CutOff{0.};
  Double_t fDtCutOff{0.}; // for when the beamspot has time

  TObjArray *fInputArray{nullptr};
  TIterator *fItInputArray{nullptr};

  TObjArray *fOutputArray{nullptr};
  TObjArray *fVertexOutputArray{nullptr};

  ClassDef(VertexFinderDA4D, 1)
};

#endif
