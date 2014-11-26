#ifndef FastJetGridMedianEstimator_h
#define FastJetGridMedianEstimator_h


/** \class FastJetGridMedianEstimator
 *
 *  Computes median energy density per event using a fixed grid.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include <map>
#include <utility>


class TObjArray;
class TIterator;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class Selector;
  namespace contrib {
    class NjettinessPlugin;
  }
}

class FastJetGridMedianEstimator: public DelphesModule
{
public:

  FastJetGridMedianEstimator();
  ~FastJetGridMedianEstimator();

  void Init();
  void Process();
  void Finish();

private:
  
  typedef std::map< std::pair< Double_t , Double_t > , std::pair< Double_t , Double_t > > TGrid; //!
    
  TGrid fGrid; //!
 
  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fRhoOutputArray; //!

  ClassDef(FastJetGridMedianEstimator, 1)
};

#endif
