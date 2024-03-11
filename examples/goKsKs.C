#include <TROOT.h>

void goKsKs()
{
gROOT->Reset();
gROOT->ProcessLine(".L examples/classes/VState.cc+");
gROOT->ProcessLine(".L examples/classes/VList.cc+");
gROOT->ProcessLine(".L examples/classes/B0KsKsPulls.cc+");
gROOT->ProcessLine(".L examples/classes/KsPulls.cc+");
}
