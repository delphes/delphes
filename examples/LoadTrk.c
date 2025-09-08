#include <TROOT.h>

void LoadTrk()
{
gROOT->Reset();
gROOT->ProcessLine(".L libDelphes.so");
gROOT->ProcessLine(".I ~/git/fbedesch/Delphes/");
gROOT->ProcessLine(".L examples/classes/SolGeom.cc+");
gROOT->ProcessLine(".L examples/classes/KalmanCk.cc+");
gROOT->ProcessLine(".L examples/KalmanCheck.c+");
}
