#ifndef VertexFinder_h
#define VertexFinder_h

/** \class VertexFinder
 *
 *  Cluster vertices from tracks
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
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

  static bool secondDescending (pair<unsigned, double>, pair<unsigned, double>);
  static bool secondAscending (pair<unsigned, double>, pair<unsigned, double>);

private:

  void createSeeds ();
  void growCluster (const unsigned);
  double weight (const unsigned);
  void addTrackToCluster (const unsigned, const unsigned);
  void removeTrackFromCluster (const unsigned, const unsigned);

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

  map<unsigned, map<string, double> > trackIDToDouble;
  map<unsigned, map<string, int> > trackIDToInt;
  map<unsigned, map<string, bool> > trackIDToBool;

  map<unsigned, map<string, double> > clusterIDToDouble;
  map<unsigned, map<string, int> > clusterIDToInt;
  map<unsigned, map<string, bool> > clusterIDToBool;
  vector<pair<unsigned, double> > trackPT;
  vector<pair<unsigned, double> > clusterSumPT2;

  ClassDef(VertexFinder, 1)
};

#endif
