/** \class BeamSpotFilter
 *
 *  Extracts beam spot
 *
 *  \author Michele Selvaggi
 *
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesModule.h"

class BeamSpotFilter: public DelphesModule
{
public:
  BeamSpotFilter() = default;

  void Init() override
  {
    // import input array
    fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));

    // create output array
    fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
  }
  void Process() override
  {
    fOutputArray->clear();
    Bool_t passed = false;
    for(Candidate *const &candidate : *fInputArray)
    {
      if(passed) break;
      if(candidate->IsPU == 0) passed = true;
      fOutputArray->emplace_back(candidate);
    }
  }

private:
  CandidatesCollection fInputArray; //!
  CandidatesCollection fOutputArray; //!
};

//------------------------------------------------------------------------------

REGISTER_MODULE("BeamSpotFilter", BeamSpotFilter);
