
/** \class ExRootAnalysisLinkDef
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

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTask.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class ExRootTreeReader+;
#pragma link C++ class ExRootTreeBranch+;
#pragma link C++ class ExRootTreeWriter+;
#pragma link C++ class ExRootResult+;
#pragma link C++ class ExRootClassifier+;
#pragma link C++ class ExRootFilter+;

#pragma link C++ class ExRootProgressBar+;
#pragma link C++ class ExRootConfReader+;
#pragma link C++ class ExRootConfParam+;
#pragma link C++ class ExRootTask+;

#pragma link C++ function HistStyle;
#pragma link C++ function FillChain;

#endif

