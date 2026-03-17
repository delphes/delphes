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

#include "TObjArray.h"

class BeamSpotFilter: public DelphesModule
{
public:
  BeamSpotFilter() = default;

  void Init() override;
  void Process() override;

private:
  const TObjArray *fInputArray{nullptr}; //!
  std::unique_ptr<TIterator> fItInputArray; //!

  TObjArray *fOutputArray{nullptr}; //!
};

//------------------------------------------------------------------------------

void BeamSpotFilter::Init()
{
  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray.reset(fInputArray->MakeIterator());

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void BeamSpotFilter::Process()
{
  Candidate *candidate = nullptr;
  Bool_t passed = false;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())) && !passed)
  {
    if(candidate->IsPU == 0) passed = true;
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

REGISTER_MODULE("BeamSpotFilter", BeamSpotFilter);
