//
#ifndef G__SOLGEOM_H
#define G__SOLGEOM_H
#include <TCanvas.h>
#include <TString.h>
//
// Class to create geometry for solenoid geometry

class SolGeom{
	//
	// Units are m
	//
private:	
	const Int_t fNlMax = 200;	// Maximum number of layers
	// B field
	Double_t fB;			// B field in Tesla
	//
	// Barrel layer properties
	Int_t fNlay;			// Total number of layers
	Int_t fBlay;			// Number of barrel layers
	Int_t fFlay;			// Number of forward/backward layers
	Int_t fNm;			// Nr. measurement layers
	Int_t    *ftyLay;	// Layer type 1 = R (barrel) or 2 = z (forward/backward)
	TString  *fLyLabl;	// Layer label 
				// Barrel: PIPE, VTXLOW, VTXHIGH, DCHCANI, DCH, DCHCANO, BSILWRP, MAG, BPRESH
				// Fw/Bw: VTXDSK, DCHWALL, FSILWRP, FRAD, FPRESH
	Double_t *fxMin;		// Minimum dimension z for barrel  or R for forward
	Double_t *fxMax;		// Maximum dimension z for barrel  or R for forward
	Double_t *frPos;		// R/z location of layer
	Double_t *fthLay;	// Thickness (meters)
	Double_t *frlLay;	// Radiation length (meters)
	Int_t    *fnmLay;	// Number of measurements in layers (1D or 2D)
	Double_t *fstLayU;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
	Double_t *fstLayL;	// Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
	Double_t *fsgLayU;	// Resolution Upper side (meters) - 0 = no measurement
	Double_t *fsgLayL;	// Resolution Lower side (meters) - 0 = no measurement
	Bool_t   *fflLay;	// measurement flag = T, scattering only = F
	TCanvas  *fcnv;		// pointer to canvas with geo drawing
	//
	Double_t fRmin;		// Radius of first measurement layer
	Double_t fZminPos;	// Z of first disk in positive direction
	Double_t fZminNeg;	// Z of first disk in positive direction
	//
	const Int_t fNdet = 10;  // Max number of tracking/passive detectors 
	Bool_t *fEnable; 	// Array of enabled detector
	// 0: beam pipe, 1: Inner VTX pixels, 2: Outer VTX silicon, 3: Drift chamber
	// 4: Outer silicon wrapper, 5: barrel pre-shower, 6: forward VTX pixels
	// 7: forward Silicon wrapper, 8: forward pre-shower
	const Int_t fNdty = 14;	// Number of detector types
	TString *fDtype;		// Array with layer labels 
	Int_t   *fDfstLay;	// Array with start layer number for given detector types
	Double_t StereoHyp(Double_t *x, Double_t *p);
	void SolGeoInit();
	void SolGeoFill();
	void GeoRead(char* fname);
	void SetMinBoundaries();
public:
	//
	// Constructors
	SolGeom();				// Create standard geometry from scratch
	SolGeom(Bool_t *OK);		// Create standard geometry with detector selection
	SolGeom(char * fname);	// Read in geometry from file
	// Destructor
	~SolGeom();
	// Accessors
	//const Int_t Ndet()				{ return fNdet; }
	Double_t B()						{ return fB; }
	Int_t    Nl()					{ return fNlay; }
	Int_t    Nm()					{ return fNm; }
	Int_t    NBl()					{ return fBlay; }
	TString  dType(Int_t ind)		{ return fDtype[ind]; }
	TString  lLabl(Int_t nlayer)	{ return fLyLabl[nlayer]; }
	Int_t    lTyp(Int_t nlayer)		{ return ftyLay[nlayer]; }
	Double_t lxMin(Int_t nlayer)		{ return fxMin[nlayer]; }
	Double_t lxMax(Int_t nlayer)		{ return fxMax[nlayer]; }
	Double_t lPos(Int_t nlayer)		{ return frPos[nlayer]; }
	Double_t lTh(Int_t nlayer)		{ return fthLay[nlayer]; }
	Double_t lX0(Int_t nlayer)		{ return frlLay[nlayer]; }
	Int_t    lND(Int_t nlayer)		{ return fnmLay[nlayer]; }
	Double_t lStU(Int_t nlayer)		{ return fstLayU[nlayer]; }
	Double_t lStL(Int_t nlayer)		{ return fstLayL[nlayer]; }
	Double_t lSgU(Int_t nlayer)		{ return fsgLayU[nlayer]; }
	Double_t lSgL(Int_t nlayer)		{ return fsgLayL[nlayer]; }
	Bool_t   isMeasure(Int_t nlayer)	{ return fflLay[nlayer]; }
	TCanvas *cnv()					{ return fcnv; }
	//
	Double_t GetRmin() { return fRmin; }
	Double_t GetZminPos() { return fZminPos; }
	Double_t GetZminNeg() { return fZminNeg; }
	//
	// B field control
	//
	void SetBfield(Double_t B) {fB = B; };
	//
	// I/O functions
	void GeoPrint(char * fname);
	//
	// Material counter
	Double_t *FracX0(Double_t theta);
	//
	// Draw geometry R - z
	void Draw();
	//
};

#endif
