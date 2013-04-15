
/** \class BTagging
 *
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags 
 *
 *  $Date$
 *  $Revision$
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/BTagging.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

class BTaggingPartonClassifier : public ExRootClassifier
{
public:

  BTaggingPartonClassifier() {}

  Int_t GetCategory(TObject *object);
  
  Double_t fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------

Int_t BTaggingPartonClassifier::GetCategory(TObject *object)
{
  Candidate *parton = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1;
  
  pdgCode = TMath::Abs(parton->PID);
  if(pdgCode != 21 && pdgCode > 5) return -1;

  return 0;
}

//------------------------------------------------------------------------------

BTagging::BTagging() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0), fItJetInputArray(0)
{
  fClassifier = new BTaggingPartonClassifier;
}

//------------------------------------------------------------------------------

BTagging::~BTagging()
{
  if(fClassifier) delete fClassifier;
}

//------------------------------------------------------------------------------

void BTagging::Init()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fDeltaR = GetDouble("DeltaR", 0.5);

  fClassifier->fPTMin = GetDouble("PartonPTMin", 1.0);
  fClassifier->fEtaMax = GetDouble("PartonEtaMax", 2.5);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();
  
  fEfficiencyMap.clear();
  for(i = 0; i < size/2; ++i)
  {
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());

    fEfficiencyMap[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end())
  {
    formula = new DelphesFormula;
    formula->Compile("0.0");

    fEfficiencyMap[0] = formula;
  }

  // import input array(s)

  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fFilter = new ExRootFilter(fPartonInputArray);
  
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();
}

//------------------------------------------------------------------------------

void BTagging::Finish()
{
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  if(fFilter) delete fFilter;
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;

  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap)
  {
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void BTagging::Process()
{
  Candidate *jet, *parton;
  Double_t pt, eta, phi;
  TObjArray *partonArray;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  Int_t pdgCode, pdgCodeMax;

  // select quark and gluons
  fFilter->Reset();
  partonArray = fFilter->GetSubArray(fClassifier, 0);
  
  if(partonArray == 0) return;

  TIter itPartonArray(partonArray);
  
  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
  {
    const TLorentzVector &jetMomentum = jet->Momentum;
    pdgCodeMax = -1;
    eta = jetMomentum.Eta();
    phi = jetMomentum.Phi();
    pt = jetMomentum.Pt();

    // loop over all input partons
    itPartonArray.Reset();
    while((parton = static_cast<Candidate*>(itPartonArray.Next())))
    {
      pdgCode = TMath::Abs(parton->PID);
      if(pdgCode == 21) pdgCode = 0;
      if(jetMomentum.DeltaR(parton->Momentum) <= fDeltaR)
      {
        if(pdgCodeMax < pdgCode) pdgCodeMax = pdgCode;
      }
    }
    if(pdgCodeMax == 0) pdgCodeMax = 21;
    if(pdgCodeMax == -1) pdgCodeMax = 0;

    // find an efficency formula
    itEfficiencyMap = fEfficiencyMap.find(pdgCodeMax);
    if(itEfficiencyMap == fEfficiencyMap.end())
    {
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;

    // apply an efficency formula
    jet->BTag = gRandom->Uniform() <= formula->Eval(pt, eta);
  }
}

//------------------------------------------------------------------------------
