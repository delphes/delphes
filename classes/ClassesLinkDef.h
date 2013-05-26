
/** ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date: 2008-07-08 12:01:52 $
 *  $Revision: 1.2 $
 *
 *  
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesFactory.h"

#include "classes/SortableObject.h"
#include "classes/DelphesClasses.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DelphesModule+;
#pragma link C++ class DelphesFactory+;

#pragma link C++ class SortableObject+;

#pragma link C++ class Event+;
#pragma link C++ class LHCOEvent+;
#pragma link C++ class LHEFEvent+;
#pragma link C++ class HepMCEvent+;
#pragma link C++ class GenParticle+;
#pragma link C++ class MissingET+;
#pragma link C++ class ScalarHT+;
#pragma link C++ class Rho+;
#pragma link C++ class Weight+;
#pragma link C++ class Photon+;
#pragma link C++ class Electron+;
#pragma link C++ class Muon+;
#pragma link C++ class Jet+;
#pragma link C++ class Track+;
#pragma link C++ class Tower+;

#pragma link C++ class Candidate+;

#endif

