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
  using DelphesModule::DelphesModule;

  void Init() override
  {
    fInputArray = ImportArray(Steer<std::string>("InputArray", "Delphes/allParticles"));
    fOutputArray = ExportArray(Steer<std::string>("OutputArray", "filteredParticles"));
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
