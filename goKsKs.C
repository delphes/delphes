#include <TROOT.h>
#include <TString.h>

void goKsKs()
{
gROOT->Reset();
gROOT->ProcessLine(".I /home/users/bedeschi/git/fbedesch/delphes2/delphes");
gROOT->ProcessLine(".L examples/VState.cc+");
gROOT->ProcessLine(".L examples/VList.cc+");
gROOT->ProcessLine(".L examples/B0KsKsPulls.cc+");
}
