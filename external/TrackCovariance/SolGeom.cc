#include <iostream>
#include <sstream>

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
  Int_t tyLay;
  string LyLabl;
  Double_t xMin;
  Double_t xMax;
  Double_t rPos;
  Double_t thLay;
  Double_t rlLay;
  Int_t nmLay;
  Double_t stLayU;
  Double_t stLayL;
  Double_t sgLayU;
  Double_t sgLayL;
  Int_t flLay;

  stringstream data_stream(data);
  string line;

  fNlay = 0;
  while(getline(data_stream, line))
  {
    stringstream line_stream(line);

    line_stream >> tyLay >> LyLabl >> xMin >> xMax >> rPos >> thLay >> rlLay >> nmLay >> stLayU >> stLayL >> sgLayU >> sgLayL >> flLay;

    if(line_stream.fail()) continue;

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
//
// Define inner box for fast tracking
//
    SetMinBoundaries();
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

//
// Get inner boundaries of cylindrical box for fast simulation
//
void SolGeom::SetMinBoundaries()
{
	// Get radius of first barrel layer
	fRmin = 1000000.0;
	fZminPos = 1000000.0;
	fZminNeg = -1000000.0;
	for (Int_t i = 0; i < fNlay; i++){
		if (ftyLay[i] == 1) {				// Cylinders
			if (frPos[i] < fRmin) fRmin = frPos[i];
		}
		if (ftyLay[i] == 2) {				// Disks
			if (frPos[i] > 0.0 && frPos[i] < fZminPos) fZminPos = frPos[i];	// Positive direction
			if (frPos[i] < 0.0 && frPos[i] > fZminNeg) fZminNeg = frPos[i];	// Negative direction
		}
	}
}

//
// Set magnetic field strength manually 
//
void SolGeom::SetBz(const Double_t Bz)
{
	fB = Bz;
}
