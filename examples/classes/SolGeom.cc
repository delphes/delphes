#include <TMath.h>
#include <TGraph.h>
#include <iostream>
#include <TCanvas.h>
#include <TPave.h>
#include <TLine.h>
#include <TPolyLine.h>
#include <TF1.h>
#include <TString.h>
#include "examples/classes/SolGeom.h"

SolGeom::SolGeom()
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = kTRUE;	// default is everything enabled
	SolGeoFill();
	SetMinBoundaries();
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer "<<i<<" reached" << endl;
				return;
			}
		}
	}
	for (Int_t i = 0; i < k; i++)cout << "i = " << i << ", Detector = " << fDtype[i] << endl;
}
//
SolGeom::SolGeom(Bool_t *OK)
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = OK[i];	// User defined list
	SolGeoFill(); 
	SetMinBoundaries();
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer " << i << " reached" << endl;
				return;
			}
		}
	}
}
SolGeom::SolGeom(char* fname)
{
	SolGeoInit();
	for (Int_t i = 0; i < fNdet; i++)fEnable[i] = kTRUE;	// default is everything enabled
	GeoRead(fname);
	SetMinBoundaries();
	TString OldLab = " "; Int_t k = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] != OldLab)
		{
			fDtype[k] = fLyLabl[i];
			fDfstLay[k] = i;
			OldLab = fLyLabl[i];
			k++;
			if (k > fNdty)
			{
				cout << "SolGeom::SolGeom : Too many detector types! Layer " << i << " reached" << endl;
				return;
			}
		}
	}
	for (Int_t i = 0; i < k; i++)cout << "i = " << i << ", Detector = " << fDtype[i] << endl;
}
void SolGeom::SolGeoInit()
{
	//
	// Magnetic field
	//
	fB = 2.0;
	//
	// Create arrays
	//
	ftyLay = new Int_t[fNlMax];		// Layer type 1 = R (barrel) or 2 = z (forward/backward)
	fLyLabl = new TString[fNlMax];	// Layer label
	fxMin = new Double_t[fNlMax];	// Minimum dimension z for barrel  or R for forward
	fxMax = new Double_t[fNlMax];	// Maximum dimension z for barrel  or R for forward
	frPos = new Double_t[fNlMax];	// R/z location of layer
	fthLay = new Double_t[fNlMax];	// Thickness (meters)
	frlLay = new Double_t[fNlMax];	// Radiation length (meters)
	fnmLay = new Int_t[fNlMax];	// Number of measurements in layers (1D or 2D)
	fstLayU = new Double_t[fNlMax];	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
	fstLayL = new Double_t[fNlMax];	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
	fsgLayU = new Double_t[fNlMax];	// Resolution Upper side (meters) - 0 = no measurement
	fsgLayL = new Double_t[fNlMax];	// Resolution Lower side (meters) - 0 = no measurement
	fflLay = new Bool_t[fNlMax];		// measurement flag = T, scattering only = F
	fEnable = new Bool_t[fNdet];		// list of enabled detectors
	fDtype = new TString[fNdty];		// Array with layer labels 
	fDfstLay = new Int_t[fNdty];		// Array with start layer
	//
	// Load geometry info in SolGeom.h
	//
	fNlay = 0;	// Actual number of layers
	fBlay = 0;	// Nr. of barrel layers
	fFlay = 0;	// Nr. of forward/backward layers
	fNm = 0;  // Nr. of measuring layers
}
	//
void SolGeom::SolGeoFill()
{
	//===================================================================================
	//		BARREL REGION
	//===================================================================================
	//
	Double_t R12 = TMath::Sqrt(12);
	//
	// Beam pipe
	//
	if (fEnable[0])
	{
		ftyLay[fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
		fLyLabl[fNlay] = "PIPE";
		fxMin[fNlay] = -100.;		// Minimum dimension z for barrel  or R for forward
		fxMax[fNlay] = 100.;		// Maximum dimension z for barrel  or R for forward
		frPos[fNlay] = 0.015;		// R/z location of layer
		fthLay[fNlay] = 0.0012;		// Thickness (meters)
		frlLay[fNlay] = 35.276e-2;	// Radiation length (meters)
		fnmLay[fNlay] = 0;			// Number of measurements in layers (1D or 2D)
		fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
		fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
		fsgLayU[fNlay] = 0.;			// Resolution Upper side (meters) - 0 = no measurement
		fsgLayL[fNlay] = 0.;			// Resolution Lower side (meters) - 0 = no measurement
		fflLay[fNlay] = kFALSE;		// measurement flag = T, scattering only = F
		fNlay++; fBlay++;
	}
	//
	// Vertex  detector (inner)
	if (fEnable[1])
	{
		const Int_t NlVtx = 3;	// Assume 3 vertex pixel layers
		Double_t rVtx[NlVtx] = { 1.7, 2.3, 3.1 };		// Vertex layer radii in cm
		Double_t lVtx[NlVtx] = { 12.0, 16.0, 16.0 };		// Vertex layer half length in cm
		for (Int_t i = 0; i < NlVtx; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXLOW";			// Layer label
			fxMin[fNlay] = -lVtx[i] * 1.e-2;		// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lVtx[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rVtx[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 280.E-6;			// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;		// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 3.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 3.E-6;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
			fNm++;
		}
	}
	//
	// Vertex  detector (outer)
	if (fEnable[2])
	{
		const Int_t NlVtxo = 2;	// Assume 2 vertex strip layers
		Double_t rVtxo[NlVtxo] = { 32., 34. };		// Vertex layer radii in cm
		Double_t lVtxo[NlVtxo] = { 100., 105 };		// Vertex layer half length in cm
		for (Int_t i = 0; i < NlVtxo; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "VTXHIGH";			// Layer label
			fxMin[fNlay] = -lVtxo[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lVtxo[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rVtxo[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 470.E-6;			// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			if (i == 0) fnmLay[fNlay] = 2;
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0 = axial layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - pi/2 = z layer - Lower side
			fsgLayU[fNlay] = 7.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 7.E-6;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			//if (i == 0) fflLay[fNlay] = kFALSE;
			fNlay++; fBlay++;
			fNm++;
		}
	}
	//
	// Tracker inner can wall
	//
	if (fEnable[3])
	{
		ftyLay[fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
		fLyLabl[fNlay] = "DCHCANI";	// Layer label
		fxMin[fNlay] = -2.125;		// Minimum dimension z for barrel  or R for forward
		fxMax[fNlay] = 2.125;			// Maximum dimension z for barrel  or R for forward
		frPos[fNlay] = 0.345;		// R/z location of layer
		fthLay[fNlay] = 0.0002;		// Thickness (meters)
		frlLay[fNlay] = 23.72226e-2;	// Radiation length (meters)
		fnmLay[fNlay] = 0;			// Number of measurements in layers (1D or 2D)
		fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
		fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
		fsgLayU[fNlay] = 0.;			// Resolution Upper side (meters) - 0 = no measurement
		fsgLayL[fNlay] = 0.;			// Resolution Lower side (meters) - 0 = no measurement
		fflLay[fNlay] = kFALSE;		// measurement flag = T, scattering only = F
		fNlay++; fBlay++;
		//
		// Drift chamber  detector
		//const Int_t NlDch = 112;			// Assume 112 DCH layers
		const Int_t NlDch = 112;
		Double_t sgDch = 100.E-6;		// Drift chamber resolution in m
		Double_t rMin = frPos[fNlay - 1] + 1.5e-2;	// Min radius @ z = 2 meters
		Double_t rMax = 2.0;						    // Max radius @ z = 2 meters
		Double_t rStp = (rMax - rMin) / (Double_t)(NlDch - 1); // layer spacing
		Double_t stAng = 13.*TMath::Pi() / 180.; // Stereo angle expressed as plate rotation angle
		//Double_t stAng = 0;
		//
		for (Int_t i = 0; i < NlDch; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "DCH";				// Layer label
			fxMin[fNlay] = -2.;					// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = 2.;					// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rMin + i*rStp;		// R/z location of layer @ z = 0
			fthLay[fNlay] = rStp;				// Thickness (meters)
			frlLay[fNlay] = 1400.;				// Radiation length (meters) - combination of gas and wires
			fnmLay[fNlay] = 1;					// Number of measurements in layers (1D or 2D)
			Double_t tgt = pow(-1, i) * 2 * frPos[fNlay] *
				TMath::Sin(stAng / 2.) / (2 * fxMax[fNlay]);	// Stereo angle (rad)  - Upper side
			fstLayU[fNlay] = TMath::ATan(tgt);
			fstLayL[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 100.E-6;			// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;					// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fBlay++;
			fNm++;
		}
		//
		// Tracker outerer can wall
		//
		ftyLay[fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
		fLyLabl[fNlay] = "DCHCANO";	// Layer label
		fxMin[fNlay] = -2.125;		// Minimum dimension z for barrel  or R for forward
		fxMax[fNlay] = 2.125;			// Maximum dimension z for barrel  or R for forward
		frPos[fNlay] = 2.02;		// R/z location of layer
		fthLay[fNlay] = 0.02;		// Thickness (meters)
		frlLay[fNlay] = 166.7e-2;	// Radiation length (meters)
		fnmLay[fNlay] = 0;			// Number of measurements in layers (1D or 2D)
		fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
		fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
		fsgLayU[fNlay] = 0.;			// Resolution Upper side (meters) - 0 = no measurement
		fsgLayL[fNlay] = 0.;			// Resolution Lower side (meters) - 0 = no measurement
		fflLay[fNlay] = kFALSE;		// measurement flag = T, scattering only = F
		fNlay++; fBlay++;
	}
	//
	// Silicon wrapper barrel region
	//
	if (fEnable[4])
	{
		const Int_t NlVtxw = 2;							// Assume 2 strip layers
		Double_t rVtxw[NlVtxw] = { 204., 206. };		// Si wrapper layer radii in cm
		Double_t lVtxw[NlVtxw] = { 235., 235. };		// Si wrapper layer half length in cm
		for (Int_t i = 0; i < NlVtxw; i++)
		{
			ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "BSILWRP";			// Layer label
			fxMin[fNlay] = -lVtxw[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = lVtxw[i] * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = rVtxw[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 470.E-6;			// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 7.E-6;				// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 90.E-6;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			//if (i == 1)fflLay[fNlay] = kFALSE;
			fNlay++; fBlay++;
			fNm++;
		}
	}
	//
	// Magnet
	//
	ftyLay [fNlay] = 1;			// Layer type 1 = R (barrel) or 2 = z (forward/backward)
	fLyLabl[fNlay] = "MAG";		// Layer label
	fxMin  [fNlay] = -2.5;		// Minimum dimension z for barrel  or R for forward
	fxMax  [fNlay] = 2.5;		// Maximum dimension z for barrel  or R for forward
	frPos  [fNlay] = 2.25;		// R/z location of layer
	fthLay [fNlay] = 0.05;		// Thickness (meters)
	frlLay [fNlay] = 6.58e-2;	// Radiation length (meters)
	fnmLay [fNlay] = 0;			// Number of measurements in layers (1D or 2D)
	fstLayU[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
	fstLayL[fNlay] = 0;			// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
	fsgLayU[fNlay] = 0.;			// Resolution Upper side (meters) - 0 = no measurement
	fsgLayL[fNlay] = 0.;			// Resolution Lower side (meters) - 0 = no measurement
	fflLay [fNlay] = kFALSE;		// measurement flag = T, scattering only = F
	fNlay++; fBlay++;
	//
	// Preshower
	if (fEnable[5])
	{
		ftyLay[fNlay] = 1;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
		fLyLabl[fNlay] = "BPRESH";			// Layer label
		fxMin[fNlay] = -2.55;				// Minimum dimension z for barrel  or R for forward
		fxMax[fNlay] = 2.55;				// Maximum dimension z for barrel  or R for forward
		frPos[fNlay] = 2.45;				// R/z location of layer
		fthLay[fNlay] = .02;				// Thickness (meters)
		frlLay[fNlay] = 100e-2;			// Radiation length (meters)
		fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
		fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
		fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
		fsgLayU[fNlay] = 70.E-6;				// Resolution Upper side (meters) - 0 = no measurement
		fsgLayL[fNlay] = 1.E-2;				// Resolution Lower side (meters) - 0 = no measurement
		fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
		fNlay++; fBlay++;
		fNm++;
	}

	//================================================================================================
	//		FORWARD/BACKWARD
	//================================================================================================
	//
	// Vertex disks
	if (fEnable[6])
	{
		const Int_t NlVtxd = 8;							// Assume 8 pixel disk layers
		Double_t zVtxd[NlVtxd] = { -92., -90., -42., -40., 40., 42., 90., 92. };		// z location in cm
		Double_t rinVtxd[NlVtxd] = { 14.1, 13.8, 6.5, 6.2, 6.2, 6.5, 13.8, 14.1 };      // Lower radius in cm
		Double_t rotVtxd = 30.0;			// Outer radius in cm
		for (Int_t i = 0; i < NlVtxd; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward
			fLyLabl[fNlay] = "VTXDSK";			// Layer label
			fxMin[fNlay] = rinVtxd[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = rotVtxd * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = zVtxd[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 280.E-6;			// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 7.E-6;			// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 7.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
	}
	//
	// DCH side wall	
	if (fEnable[3])
	{
		const Int_t nW = 2;
		for (Int_t i = 0; i < nW; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "DCHWALL";			// Layer label
			fxMin[fNlay] = 34.5*1e-2;			// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = 202.0 * 1.e-2;		// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = pow(-1, i)*212.5 * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 25.E-2;				// Thickness (meters)
			frlLay[fNlay] = 555.0e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 0;					// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;					// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
	}
	//
	// Forw/Backw. Si wrapper
	if (fEnable[7])
	{
		const Int_t NlSid = 4;										// Assume 4 disk layers
		Double_t zSid[NlSid] = { -232., -230., 230., 232. };		// z location in cm
		Double_t rinSid[NlSid] = { 35.4, 35.0, 35.0, 35.4 };		// Lower radius in cm
		Double_t rotSid = 202.0;										// Outer radius in cm
		for (Int_t i = 0; i < NlSid; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "FSILWRP";			// Layer label
			fxMin[fNlay] = rinSid[i] * 1.e-2;	// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = rotSid * 1.e-2;	// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = zSid[i] * 1.e-2;	// R/z location of layer
			fthLay[fNlay] = 470.E-6;			// Thickness (meters)
			frlLay[fNlay] = 9.370e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 7.E-6;			// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 90.E-6;			// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
	}
	//
	// pre-shower radiators
	if (fEnable[8])
	{
		const Int_t nPS = 2;
		for (Int_t i = 0; i < nPS; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "FRAD";				// Layer label
			fxMin[fNlay] = 38 * 1e-2;			// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = 209.0 * 1.e-2;		// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = pow(-1, i)*249.* 1.e-2;	// R/z location of layer
			fthLay[fNlay] = .43E-2;				// Thickness (meters)
			frlLay[fNlay] = 0.5612e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 0;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = 0;					// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 0;					// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 0;					// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kFALSE;			// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
		//
		// Forw/backw. preshower
		const Int_t NlPSd = 2;										// Assume 2 disk layers
		Double_t zPSd[NlPSd] = { -2.55, 2.55 };					// z location in m
		Double_t rinPSd[NlPSd] = { 0.39, 0.39 };					// Lower radius in m
		Double_t rotPSd = 2.43;										// Outer radius in m
		for (Int_t i = 0; i < NlPSd; i++)
		{
			ftyLay[fNlay] = 2;					// Layer type 1 = R (barrel) or 2 = z (forward/backward)
			fLyLabl[fNlay] = "FPRESH";			// Layer label
			fxMin[fNlay] = rinPSd[i];			// Minimum dimension z for barrel  or R for forward
			fxMax[fNlay] = rotPSd;			// Maximum dimension z for barrel  or R for forward
			frPos[fNlay] = zPSd[i];			// R/z location of layer
			fthLay[fNlay] = 0.02;				// Thickness (meters)
			frlLay[fNlay] = 100.0e-2;			// Radiation length (meters)
			fnmLay[fNlay] = 2;					// Number of measurements in layers (1D or 2D)
			fstLayU[fNlay] = 0.0;				// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
			fstLayL[fNlay] = TMath::Pi() / 2.;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
			fsgLayU[fNlay] = 70.E-6;			// Resolution Upper side (meters) - 0 = no measurement
			fsgLayL[fNlay] = 1.E-2;				// Resolution Lower side (meters) - 0 = no measurement
			fflLay[fNlay] = kTRUE;				// measurement flag = T, scattering only = F
			fNlay++; fFlay++;
			fNm++;
		}
	}

	cout << "Geometry created with " << fNlay << "/" << fNm << " layers" << endl;
}
//
// Get inner boundaries
//
void SolGeom::SetMinBoundaries()
{
	// Get radius of first barrel measurement layer
	fRmin = 1000000.0;
	fZminPos = 1000000.0;
	fZminNeg = -1000000.0;
	for (Int_t i = 0; i < fNlay; i++){
		if (ftyLay[i] == 1) {						// Cylinders
			if (frPos[i] < fRmin) fRmin = frPos[i];
		}
		if (ftyLay[i] == 2) {						// Disks
			if (frPos[i] > 0.0 && frPos[i] < fZminPos) fZminPos = frPos[i];	// Positive direction
			if (frPos[i] < 0.0 && frPos[i] > fZminNeg) fZminNeg = frPos[i];	// Negative direction
		}
	}
}
//
// Print geometry
void SolGeom::GeoPrint(char* fname)
{
	FILE *fdata = fopen(fname, "w");
	if (!fdata)
	{
		cout << "SolGeom::GeoPrint - can't open output file" << endl;
		return;
	}
	for (Int_t l = 0; l < fNlay; l++)
	{
		fprintf(fdata, "%d %s %g %g %g %g %g %d %g %g %g %g %d\n",
		ftyLay[l], fLyLabl[l].Data(), fxMin[l], fxMax[l], frPos[l], fthLay[l],
		frlLay[l], fnmLay[l], fstLayU[l], fstLayL[l], fsgLayU[l], fsgLayL[l], fflLay[l]);
		//cout << strng << endl<< endl;
	}	
	fclose(fdata);
}
//
// Material counter (Units are fraction of X0)
Double_t *SolGeom::FracX0(Double_t theta)
{
	//
	// Calculates amount of material crossed by a straight track at a polar angle theta
	// for each subdetector:
	// 0: Pipe, 1: VTXLOW, 2: VTXHIGH, 3: DCHCANI, 4: DCH, 5: DCHCANO, 6: BSILWRP, 7: MAG,
	// 8: BPRESH, 9: VTXDSK, 10: DCHWALL, 11: FSILWRP, 12: FRAD, 13: FPRESH
	//
	Double_t *Mat;
	Mat = new Double_t[fNdty];
	for (Int_t i = 0; i < fNdty; i++)Mat[i] = 0;
	if (fNlay <= 0)
	{
		cout << "SolGeom::FracX0 : No geometry available. # layers = " << fNlay << endl;
		return Mat;
	}
	//
	// Loop over all layers
	Double_t lmb = 0.0;
	if (TMath::Abs(theta - TMath::PiOver2()) > 1.0e-10 && 
		TMath::Abs(theta)  > 1.0e-10) lmb = 1.0 / TMath::Tan(theta);	// Cot(theta)
	if (theta == 0.0) lmb = 1.e10;
	for (Int_t il = 0; il < fNlay; il++)
	{
		Int_t dNum;
		for (Int_t i = 0; i<fNdty; i++) if (fLyLabl[il] == fDtype[i])dNum = i;
		//cout << "dnum = " << dNum << ", detector: "<<fDtype[dNum]<<endl;
		if (ftyLay[il] == 1)		// Cylinder at constant R
		{
			Double_t R = frPos[il];
			Double_t z = lmb*R;
			//cout << "l num: " << il << ", R = " << R << ", z min: "<<fxMin[il]<<", z = " << z<<" , z max: "<<fxMax[il] << endl;
			if (z>fxMin[il] && z < fxMax[il])	// the layer is hit
			{
				Mat[dNum] += fthLay[il] / (TMath::Sin(theta)*frlLay[il]);
			}
		}
		else if (ftyLay[il] == 2) // disk at constant z
		{
			Double_t z = frPos[il];
			Double_t R = z / lmb;
			//cout << "l num: " << il << ", z = " << z << ", R min: " << fxMin[il] << ", R = " << R << " , R max: " << fxMax[il] << endl;
			if (R>fxMin[il] && R < fxMax[il])	// the layer is hit
			{
				Mat[dNum] += fthLay[il] / (TMath::Cos(theta)*frlLay[il]);
			}
		}
	}
	//
	return Mat;
}
//
// Read geometry
void SolGeom::GeoRead(char* fname)
{
	cout << "SolGeom::GeoRead: " << fname << endl;
	char strng[200];
	int nbytes = 200;
	FILE *fdata = fopen(fname, "r");
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
	//
	cout << "SolGeom::GeoRead: before fgets" << endl;
	while (fgets(strng, nbytes, fdata) != NULL)
	{
		cout << strng;
		int status = sscanf(strng, "%d %s %g %g %g %g %g %d %g %g %g %g %d",
			&tyLay, LyLabl, &xMin, &xMax, &rPos, &thLay,
			&rlLay, &nmLay, &stLayU, &stLayL, &sgLayU, &sgLayL, &flLay);
		ftyLay[fNlay] = tyLay; 
		fLyLabl[fNlay] = LyLabl; 
		fxMin[fNlay] = (Double_t) xMin; 
		fxMax[fNlay] = (Double_t) xMax; 
		frPos[fNlay] = (Double_t) rPos; 
		fthLay[fNlay] = (Double_t) thLay; 
		frlLay[fNlay] = (Double_t) rlLay; 
		fnmLay[fNlay] = nmLay; 
		fstLayU[fNlay] = (Double_t) stLayU; 
		fstLayL[fNlay] = (Double_t) stLayL; 
		fsgLayU[fNlay] = (Double_t) sgLayU; 
		fsgLayL[fNlay] = (Double_t) sgLayL; 
		fflLay[fNlay] = (Bool_t) flLay; 
		cout << "Layer # " << fNlay << ": " << fLyLabl[fNlay] << ", Position: " << frPos[fNlay]
			<< ", Measurement: " << fflLay[fNlay] << endl;
		
		fNlay++;
		if (tyLay == 1)fBlay++;
		if (flLay == 1)fNm++;

	}
	fclose(fdata);
	cout << "SolGeom::GeoRead completed with " << fNlay << " layers input" << endl;
}
//
// Destructor
SolGeom::~SolGeom()
{
	fNlay = 0;
	fBlay = 0;
	fNm = 0;

	delete[] & ftyLay;
	delete[] & fxMin;
	delete[] & fxMax;
	delete[] & frPos;
	delete[] & fthLay;
	delete[] & frlLay;
	delete[] & fnmLay;
	delete[] & fstLayU;
	delete[] & fstLayL;
	delete[] & fsgLayU;
	delete[] & fsgLayL;
	delete[] & fflLay;
	delete[] & fEnable;
}
//
// Draw the geometry (just a sketch)
//
void SolGeom::Draw()
{
	Double_t zMin = -2.75; Double_t zMax = 2.75;
	Double_t rMax = 2.6;
	fcnv = new TCanvas("cnv", "Geometry sketch", 10, 10, 950, 550);
	fcnv->Range(zMin, -0.1, zMax, rMax);
	// 
	// beam pipe
	if (fEnable[0])
	{
		TPave *pipe = new TPave(zMin, -frPos[0], zMax, frPos[0], 0, "");
		pipe->SetFillColor(kYellow);
		pipe->Draw();
	}
	// Beamline
	TLine *beam = new TLine(zMin, 0.0, zMax, 0.0);
	beam->SetLineColor(kBlack);
	beam->SetLineWidth(1);
	beam->SetLineStyle(9);
	beam->Draw("SAME");
	// Magnet
	TPave *sol = new TPave(-2.5, 2.2, 2.5, 2.3, 0, "");
	sol->SetFillColor(30);
	sol->Draw("SAME");
	//
	// Draw Calorimeter
	// Barrel
	const Int_t nP = 5;
	Double_t brCalX[nP] = { -2.6, 2.6, 4.6, -4.6, -2.6 };
	Double_t brCalY[nP] = { 2.5, 2.5, 4.5, 4.5, 2.5 };
	TPolyLine *brCalor = new TPolyLine(nP, brCalX, brCalY, "F");
	brCalor->SetFillColor(38);
	brCalor->SetLineColor(kBlack);
	brCalor->Draw("FSAME");
	// Backward
	Double_t bkCalX[nP] = { -4.6, -2.6, -2.6, -4.6, -4.6 };
	Double_t bkCalY[nP] = { 0.68, 0.39, 2.5, 4.5, 0.68 };
	TPolyLine *bkCalor = new TPolyLine(nP, bkCalX, bkCalY, "F");
	bkCalor->SetFillColor(38);
	bkCalor->SetLineColor(kBlack);
	bkCalor->Draw("FSAME");
	// Forward
	Double_t bfCalX[nP] = { 2.6, 4.6, 4.6, 2.6, 2.6 };
	Double_t bfCalY[nP] = { 0.39, 0.68, 4.5, 2.5, 0.39 };
	TPolyLine *bfCalor = new TPolyLine(nP, bfCalX, bfCalY, "F");
	bfCalor->SetFillColor(38);
	bfCalor->SetLineColor(kBlack);
	bfCalor->Draw("FSAME");
	// All other layers
	// Measurement silicon (red), blue (DCH), scattering black
	//
	const Int_t lMax = 200;
	TLine *ln[lMax];
	TF1   *fn[lMax];
	Int_t il = 0;
	Int_t ig = 0;
	for (Int_t i = 0; i < fNlay; i++)
	{
		if (fLyLabl[i] == "DCH")		// Drift chamber layers (hypeboloids)
		{
			char lab[10];
			Int_t stat;
			stat = sprintf(lab, "fun%d", ig);
			fn[ig] = new TF1(lab, this, &SolGeom::StereoHyp, lxMin(i), lxMax(i), 3, "SolGeom", "StereoHyp");
			fn[ig]->SetParameter(0, lPos(i));
			fn[ig]->SetParameter(1, lStU(i));
			fn[ig]->SetParameter(2, (Double_t)i);
			fn[ig]->SetLineColor(kBlue);
			fn[ig]->Draw("SAME");
			ig++;
		}
		else
		{
			if (ftyLay[i] == 1)ln[il] = new TLine(lxMin(i), lPos(i), lxMax(i), lPos(i));
			else ln[il] = new TLine(lPos(i), lxMin(i), lPos(i), lxMax(i));
			ln[il]->SetLineColor(kBlack);
			if (isMeasure(i))ln[il]->SetLineColor(kRed);
			ln[il]->Draw("SAME");
			il++;
		}
	}
}
//
Double_t SolGeom::StereoHyp(Double_t *x, Double_t *p)
{
	Double_t R   = p[0];
	Double_t tg  = TMath::Tan(p[1]);
	Int_t i = (Int_t)p[2];
	Double_t r = TMath::Sqrt(R*R - lxMax(i)*lxMax(i)*tg*tg);
	Double_t z   = x[0];
	//
	return TMath::Sqrt(r*r + z*z*tg*tg);
}
