#include <TROOT.h>

void goBs()
{
gROOT->Reset();
gROOT->ProcessLine(".L examples/classes/VState.cc+");
gROOT->ProcessLine(".L examples/classes/VList.cc+");
gROOT->ProcessLine(".L examples/classes/BsPulls.cc+");
}
