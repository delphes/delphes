#ifndef ExRootUtilities_h
#define ExRootUtilities_h

/** \class ExRootUtilities
 *
 *  Functions simplifying ROOT tree analysis
 *
 *  $Date: 2008-06-04 13:57:28 $
 *  $Revision: 1.1 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "Rtypes.h"

class TH1;
class TChain;

void HistStyle(TH1 *hist, Bool_t stats = kTRUE);

Bool_t FillChain(TChain *chain, const char *inputFileList);

#endif // ExRootUtilities_h
