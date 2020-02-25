#include <iostream>

#include <TString.h>

#include "SolGeom.h"
#include "SolTrack.h"

using namespace std;

SolGeom::SolGeom()
{
  //
  // Magnetic field
  //
  fB = 2.0;
  //
  // Create arrays
  //
  ftyLay = new Int_t[fNlMax];     // Layer type 1 = R (barrel) or 2 = z (forward/backward)
  fLyLabl = new TString[fNlMax];  // Layer label
  fxMin = new Double_t[fNlMax];   // Minimum dimension z for barrel  or R for forward
  fxMax = new Double_t[fNlMax];   // Maximum dimension z for barrel  or R for forward
  frPos = new Double_t[fNlMax];   // R/z location of layer
  fthLay = new Double_t[fNlMax];  // Thickness (meters)
  frlLay = new Double_t[fNlMax];  // Radiation length (meters)
  fnmLay = new Int_t[fNlMax];     // Number of measurements in layers (1D or 2D)
  fstLayU = new Double_t[fNlMax]; // Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
  fstLayL = new Double_t[fNlMax]; // Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
  fsgLayU = new Double_t[fNlMax]; // Resolution Upper side (meters) - 0 = no measurement
  fsgLayL = new Double_t[fNlMax]; // Resolution Lower side (meters) - 0 = no measurement
  fflLay = new Bool_t[fNlMax];    // measurement flag = T, scattering only = F
  //
  // Load geometry info in SolGeom.h
  //
  fNlay = 0; // Actual number of layers
  fBlay = 0; // Nr. of barrel layers
  fFlay = 0; // Nr. of forward/backward layers
  fNm = 0;   // Nr. of measuring layers
}

void SolGeom::Read(const char *data)
{
  char strng[200];
  int nbytes = 200;
  FILE *fdata = fopen(data, "r");
  if (!fdata)
  {
    cout << "SolGeom::GeoRead - can't open input file" << endl;
    return;
  }
  Int_t tyLay;
  char LyLabl[20];
  float xMin;
  float xMax;
  float rPos;
  float thLay;
  float rlLay;
  Int_t nmLay;
  float stLayU;
  float stLayL;
  float sgLayU;
  float sgLayL;
  Int_t flLay;

  fNlay = 0;
  while (fgets(strng, nbytes, fdata) != NULL)
  {
    cout << strng;
    int status = sscanf(strng, "%d %s %g %g %g %g %g %d %g %g %g %g %d",
      &tyLay, LyLabl, &xMin, &xMax, &rPos, &thLay,
      &rlLay, &nmLay, &stLayU, &stLayL, &sgLayU, &sgLayL, &flLay);
    ftyLay[fNlay] = tyLay;
    fLyLabl[fNlay] = LyLabl;
    fxMin[fNlay] = xMin;
    fxMax[fNlay] = xMax;
    frPos[fNlay] = rPos;
    fthLay[fNlay] = thLay;
    frlLay[fNlay] = rlLay;
    fnmLay[fNlay] = nmLay;
    fstLayU[fNlay] = stLayU;
    fstLayL[fNlay] = stLayL;
    fsgLayU[fNlay] = sgLayU;
    fsgLayL[fNlay] = sgLayL;
    fflLay[fNlay] = flLay;

    fNlay++;
    if (tyLay == 1) fBlay++;
    if (tyLay == 2) fFlay++;
    if (flLay == 1) fNm++;

  }

  cout << "SolGeom::GeoRead completed with " << fNlay << " layers input" << endl;
}

SolGeom::~SolGeom()
{
  delete[] ftyLay;
  delete[] fxMin;
  delete[] fxMax;
  delete[] frPos;
  delete[] fthLay;
  delete[] frlLay;
  delete[] fnmLay;
  delete[] fstLayU;
  delete[] fstLayL;
  delete[] fsgLayU;
  delete[] fsgLayL;
  delete[] fflLay;
}
