
/** \class Clone
 *
 *  Clone candidate array
 *
 *  $Date$
 *  $Revision$
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Cloner.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

Cloner::Cloner() :
  fItInputArray(0)
{

}

//------------------------------------------------------------------------------

Cloner::~Cloner()
{

}

//------------------------------------------------------------------------------

void Cloner::Init()
{
  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void Cloner::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void Cloner::Process()
{
  Candidate *candidate;
 
 // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidate = static_cast<Candidate*>(candidate->Clone());
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
