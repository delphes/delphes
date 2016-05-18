#ifndef VertexFinder_h
#define VertexFinder_h

/** \class VertexFinder
 *
 *  Cluster vertices from tracks
 *
 *  \authors A. Hart, M. Selvaggi
 *
 */


#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

class TObjArray;
class Candidate;
class TVector3;

class VertexFinder: public DelphesModule
{
public:

  VertexFinder();
  ~VertexFinder();

  void Init();
  void Process();
  void Finish();

  static Bool_t secondDescending (pair<UInt_t, Double_t>, pair<UInt_t, Double_t>);
  static Bool_t secondAscending (pair<UInt_t, Double_t>, pair<UInt_t, Double_t>);

private:

  void createSeeds ();
  void growCluster (const UInt_t);
  Double_t weight (const UInt_t);
  void addTrackToCluster (const UInt_t, const UInt_t);
  void removeTrackFromCluster (const UInt_t, const UInt_t);

  Double_t fSigma;
  Double_t fMinPT;
  Double_t fMaxEta;
  Double_t fSeedMinPT;
  Int_t fMinNDF;
  Int_t fGrowSeeds;

  TObjArray *fInputArray;
  TIterator *fItInputArray;

  TObjArray *fOutputArray;
  TObjArray *fVertexOutputArray;

  map<UInt_t, map<string, Double_t> > trackIDToDouble;
  map<UInt_t, map<string, Int_t> > trackIDToInt;
  map<UInt_t, map<string, Bool_t> > trackIDToBool;

  map<UInt_t, map<string, Double_t> > clusterIDToDouble;
  map<UInt_t, map<string, Int_t> > clusterIDToInt;
  map<UInt_t, map<string, Bool_t> > clusterIDToBool;
  vector<pair<UInt_t, Double_t> > trackPT;
  vector<pair<UInt_t, Double_t> > clusterSumPT2;

  ClassDef(VertexFinder, 1)
};

#endif
