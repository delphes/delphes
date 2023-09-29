#include <TROOT.h>

void goKs()
{
gROOT->Reset();
gROOT->ProcessLine(".L examples/classes/VState.cc+");
gROOT->ProcessLine(".L examples/classes/VList.cc+");
gROOT->ProcessLine(".L examples/classes/KsPulls.cc+");
}
