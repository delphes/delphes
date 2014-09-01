
/** \class
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2014-04-16 17:17:35 +0200 (Wed, 16 Apr 2014) $
 *  $Revision: 1369 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Delphes.h"

#include "modules/FastJetFinder.h"
#include "modules/ParticlePropagator.h"
#include "modules/Efficiency.h"
#include "modules/IdentificationMap.h"
#include "modules/EnergySmearing.h"
#include "modules/MomentumSmearing.h"
#include "modules/ImpactParameterSmearing.h"
#include "modules/TimeSmearing.h"
#include "modules/SimpleCalorimeter.h"
#include "modules/Calorimeter.h"
#include "modules/Isolation.h"
#include "modules/EnergyScale.h"
#include "modules/UniqueObjectFinder.h"
#include "modules/TrackCountingBTagging.h"
#include "modules/BTagging.h"
#include "modules/TauTagging.h"
#include "modules/TreeWriter.h"
#include "modules/Merger.h"
#include "modules/LeptonDressing.h"
#include "modules/PileUpMerger.h"
#include "modules/JetPileUpSubtractor.h"
#include "modules/TrackPileUpSubtractor.h"
#include "modules/PileUpJetID.h"
#include "modules/ConstituentFilter.h"
#include "modules/StatusPidFilter.h"
#include "modules/Cloner.h"
#include "modules/Weighter.h"
#include "modules/Hector.h"
#include "modules/ExampleModule.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class Delphes+;

#pragma link C++ class FastJetFinder+;
#pragma link C++ class ParticlePropagator+;
#pragma link C++ class Efficiency+;
#pragma link C++ class IdentificationMap+;
#pragma link C++ class EnergySmearing+;
#pragma link C++ class MomentumSmearing+;
#pragma link C++ class ImpactParameterSmearing+;
#pragma link C++ class TimeSmearing+;
#pragma link C++ class SimpleCalorimeter+;
#pragma link C++ class Calorimeter+;
#pragma link C++ class Isolation+;
#pragma link C++ class EnergyScale+;
#pragma link C++ class UniqueObjectFinder+;
#pragma link C++ class TrackCountingBTagging+;
#pragma link C++ class BTagging+;
#pragma link C++ class TauTagging+;
#pragma link C++ class TreeWriter+;
#pragma link C++ class Merger+;
#pragma link C++ class LeptonDressing+;
#pragma link C++ class PileUpMerger+;
#pragma link C++ class JetPileUpSubtractor+;
#pragma link C++ class TrackPileUpSubtractor+;
#pragma link C++ class PileUpJetID+;
#pragma link C++ class ConstituentFilter+;
#pragma link C++ class StatusPidFilter+;
#pragma link C++ class Cloner+;
#pragma link C++ class Weighter+;
#pragma link C++ class Hector+;
#pragma link C++ class ExampleModule+;

#endif
