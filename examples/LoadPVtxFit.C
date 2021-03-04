#include <TROOT.h>

void LoadPVtxFit()
{
gROOT->Reset();
gROOT->ProcessLine(".L libDelphes.so");
gROOT->ProcessLine(".L external/TrackCovariance/VertexFit.cc+");
gROOT->ProcessLine(".L examples/ExamplePVtxFit.C+");
}
