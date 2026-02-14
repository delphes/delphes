/** \class BeamSpotFilter
 *
 *  Extracts beam spot
 *
 *  \author Michele Selvaggi
 *
 */

#include "modules/BeamSpotFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{
  // import input array
  ImportArray(GetString("InputArray", "Delphes/allParticles"), fInputArray);
  // create output array
  ExportArray(fOutputArray, GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Finish()
{
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  for(const auto &candidate : *fInputArray)
  {
    fOutputArray->emplace_back(candidate);
    if(candidate.IsPU == 0) break;
  }
}
