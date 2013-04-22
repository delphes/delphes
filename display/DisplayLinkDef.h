
/** ExRootAnalysisLinkDef
 *
 *  Lists classes to be included in cint dicitonary
 *
 *  $Date$
 *  $Revision$
 *
 *  
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "display/DelphesDisplay.h"
#include "display/DelphesCaloData.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class DelphesDisplay+;
#pragma link C++ class DelphesCaloData+;

#endif

