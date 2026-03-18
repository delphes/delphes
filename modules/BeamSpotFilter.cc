/** \class BeamSpotFilter
 *
 *  Extracts beam spot
 *
 *  \author Michele Selvaggi
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesModuleFactory.h"

class BeamSpotFilter: public DelphesModule
{
public:
  BeamSpotFilter() = default;

  void Init() override;
  void Process() override;

private:
  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  fOutputArray->clear();
  Bool_t passed = false;
  for(const auto &candidate : *fInputArray)
  {
    if(passed) break;
    if(candidate->IsPU == 0) passed = true;
    fOutputArray->emplace_back(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("BeamSpotFilter", BeamSpotFilter);
