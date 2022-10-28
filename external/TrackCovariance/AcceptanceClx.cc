#include <TVector3.h>
#include "AcceptanceClx.h"
//
// Pt splitting routine
//
void AcceptanceClx::VecInsert(Int_t i, Float_t x, TVectorF& Vec)
{
	// Insert a new element x in location i+1 of vector Vec
	//
	Int_t N = Vec.GetNrows();	// Get vector nitial size
	Vec.ResizeTo(N + 1);		// Increase size
	for (Int_t j = N - 1; j > i; j--)Vec(j + 1) = Vec(j);	// Shift all elements above i
	Vec(i+1) = x;
}
//
void AcceptanceClx::SplitPt(Int_t i, TVectorF &AccPt)
{
	Int_t Nrows = fAcc.GetNrows();
	Int_t Ncols = fAcc.GetNcols();
	TMatrixF AccMod(Nrows + 1, Ncols); // Size output matrix
	for (Int_t ith = 0; ith < Ncols; ith++)
	{
		AccMod(i + 1, ith) = AccPt(ith);
		for (Int_t ipt = 0; ipt <= i; ipt++) AccMod(ipt, ith) = fAcc(ipt, ith);
		for (Int_t ipt = i + 1; ipt < Nrows; ipt++) AccMod(ipt + 1, ith) = fAcc(ipt, ith);
	}
	//
	fAcc.ResizeTo(Nrows + 1, Ncols);
	fAcc = AccMod;
}
//
// Theta splitting routine
void AcceptanceClx::SplitTh(Int_t i, TVectorF &AccTh)
{
	Int_t Nrows = fAcc.GetNrows();
	Int_t Ncols = fAcc.GetNcols();
	TMatrixF AccMod(Nrows, Ncols + 1); // Size output matrix
	for (Int_t ipt = 0; ipt < Nrows; ipt++)
	{
		AccMod(ipt, i + 1) = AccTh(ipt);
		for (Int_t ith = 0; ith <= i; ith++) AccMod(ipt, ith) = fAcc(ipt, ith);
		for (Int_t ith = i + 1; ith < Ncols; ith++) AccMod(ipt, ith + 1) = fAcc(ipt, ith);
	}
	//
	fAcc.ResizeTo(Nrows, Ncols + 1);
	fAcc = AccMod;
}
//
// Constructors
//
AcceptanceClx::AcceptanceClx(TString InFile)
{
	ReadAcceptance(InFile);
}
//
AcceptanceClx::AcceptanceClx(SolGeom* InGeo)
{
	// Initializations
	//
	//cout << "Entered constructor of AccpeptanceClx" << endl;
	// Setup grid
	// Start grid parameters
	//
	// Pt nodes
	const Int_t NpPtInp = 10;
	Float_t PtInit[NpPtInp] = { 0., 1., 10., 100., 250.,
		 		   500., 1000., 2000., 10000, 50000. };
	TVectorF Pta(NpPtInp, PtInit);
	Int_t NpPt = Pta.GetNrows();	// Nr. of starting pt points
	// Theta nodes
	const Int_t NpThInp = 15;
	Float_t ThInit[NpThInp] = { 0.,5.,10.,20.,30.,40.,50.,90.,
	130.,140.,150.,160.,170.,175., 180. };
	TVectorF Tha(NpThInp, ThInit);
	Int_t NpTh = Tha.GetNrows();	// Nr. of starting theta points
	//cout << "AcceptanceClv:: Pta and Tha arrays defined" << endl;
	//
	// Grid splitting parameters
	Float_t dPtMin = 0.2;		// Grid Pt resolution (GeV)
	Float_t dThMin = 2.0;		// Grid Theta resolution (degrees)
	Float_t dAmin  = 1.0;		// Minimum # hits step
	//
	fAcc.ResizeTo(NpPt, NpTh);
	//
	//
	// Event loop: fill matrix starting nodes
	//
	TVector3 xv(0., 0., 0.);
	for (Int_t ipt = 0; ipt < NpPt; ipt++)	// Scan pt bins
	{
		// Momentum in GeV
		Double_t pt = Pta(ipt);
		//
		for (Int_t ith = 0; ith < NpTh; ith++)	// Scan theta bins
		{
			//
			// Theta in from degrees to radians
			Double_t th = TMath::Pi() * Tha(ith) / 180.;
			Double_t pz = pt / TMath::Tan(th);
			TVector3 tp(pt, 0., pz);
			//
			// Get number of measurement hits
			//
			SolTrack* gTrk = new SolTrack(xv, tp, InGeo);	// Generated track
			Int_t Mhits = gTrk->nmHit();			// Nr. Measurement hits
			fAcc(ipt, ith) = (Float_t)Mhits;
		}
	}
	//
	// Scan nodes and split if needed
	//
	Int_t Nsplits = 1;		// Number of split per iteration
	Int_t Ncycles = 0;		// Number of iterations
	Int_t MaxSplits = 200;		// Maximum number of splits
	Int_t NsplitCnt = 0;
	//
	while (Nsplits > 0)
	{
		Nsplits = 0;
		// Scan nodes
		for (Int_t ipt = 0; ipt < NpPt - 1; ipt++)	// Scan pt bins
		{
			for (Int_t ith = 0; ith < NpTh - 1; ith++)	// Scan theta bins
			{
				Float_t dAp = TMath::Abs(fAcc(ipt + 1, ith) - fAcc(ipt, ith));
				Float_t dPt = TMath::Abs(Pta(ipt + 1) - Pta(ipt));
				//
				// Pt split
				if (dPt > dPtMin && dAp > dAmin && Nsplits < MaxSplits) {
					NsplitCnt++;	// Total splits counter
					Nsplits++;	// Increase splits/cycle
					NpPt++;		// Increase #pt points 
					Float_t newPt = 0.5 * (Pta(ipt + 1) + Pta(ipt));
					VecInsert(ipt, newPt, Pta);
					TVectorF AccPt(NpTh);
					for (Int_t i = 0; i < NpTh; i++)
					{
						Double_t pt = newPt;
						Double_t th = TMath::Pi() * Tha[i] / 180.;
						Double_t pz = pt / TMath::Tan(th);
						TVector3 tp(pt, 0., pz);
						SolTrack* gTrk = new SolTrack(xv, tp, InGeo);	// Generated track
						Int_t Mhits = gTrk->nmHit();			// Nr. Measurement hits
						AccPt(i) = (Float_t)Mhits;
					}
					SplitPt(ipt, AccPt);
					// Completed Pt split
				}
				//
				Float_t dAt = TMath::Abs(fAcc(ipt, ith + 1) - fAcc(ipt, ith));
				Float_t dTh = TMath::Abs(Tha(ith + 1) - Tha(ith));
				//
				// Theta split
				if (dTh > dThMin && dAt > dAmin && Nsplits < MaxSplits) {
					//cout << "Th(" << ith << ") = " << Tha(ith) << ", dAt = " << dAt << endl;
					NsplitCnt++;	// Total splits counter
					Nsplits++;	// Increase splits
					NpTh++;		// Increase #pt points 
					Float_t newTh = 0.5 * (Tha(ith + 1) + Tha(ith));
					VecInsert(ith, newTh, Tha);
					TVectorF AccTh(NpPt);
					for (Int_t i = 0; i < NpPt; i++)
					{
						Double_t pt = Pta(i);
						Double_t th = TMath::Pi() * newTh / 180.;
						Double_t pz = pt / TMath::Tan(th);
						TVector3 tp(pt, 0., pz);
						SolTrack* gTrk = new SolTrack(xv, tp, InGeo);	// Generated track
						Int_t Mhits = gTrk->nmHit();			// Nr. Measurement hits
						AccTh(i) = (Float_t)Mhits;
					}
					SplitTh(ith, AccTh);
					// Theta splits completed
				}
			}  // End loop on theta nodes
		}  // End loop on pt nodes
		Ncycles++;
	}			// End loop on iteration
	//
	std::cout<<"AcceptanceClx:: Acceptance generation completed after "<<Ncycles<<" cycles"
		<<" and a total number of "<<NsplitCnt<< " splits"<<std::endl;
	//
	// Store final variables and parameters
	//
	fNPtNodes = NpPt;
	Int_t PtCk = Pta.GetNrows();
	if (PtCk != NpPt)
		std::cout << "AcceptanceClx:: Error in grid generation NpPt=" << NpPt << ", Pta size= " << PtCk << std::endl;
	fPtArray.ResizeTo(NpPt);
	fPtArray = Pta;	// Array of Pt nodes
	//
	fNThNodes = NpTh;
	Int_t ThCk = Tha.GetNrows();
	if (ThCk != NpTh)
		std::cout << "AcceptanceClx:: Error in grid generation NpTh=" << NpTh << ", Tha size= " << ThCk << std::endl;
	fThArray.ResizeTo(NpTh);
	fThArray = Tha;	// Array of Theta nodes
	//
		std::cout << "AcceptanceClx:: Acceptance encoding with " << fNPtNodes
		<<" pt nodes and "<< fNThNodes <<" theta nodes"<< std::endl;
	Int_t Nrows = fAcc.GetNrows();
	Int_t Ncols = fAcc.GetNcols();
}

// Destructor
AcceptanceClx::~AcceptanceClx()
{
	fNPtNodes = 0;
	fNThNodes = 0;
	fAcc.Clear();
	fPtArray.Clear();
	fThArray.Clear();
}
//
void AcceptanceClx::WriteAcceptance(TFile *fout)
{
	//
	// Write out data
	TTree* tree = new TTree("treeAcc", "Acceptance tree");
	TMatrixF* pntMatF = &fAcc;
	TVectorF* pntVecP = &fPtArray;
	TVectorF* pntVecT = &fThArray;
	tree->Branch("AcceptanceMatrix", "TMatrixF", &pntMatF, 64000, 0);
	tree->Branch("AcceptancePtVec",  "TVectorF", &pntVecP, 64000, 0);
	tree->Branch("AcceptanceThVec",  "TVectorF", &pntVecT, 64000, 0);
	tree->Fill();
	fout->Write();
}
//
void AcceptanceClx::WriteAcceptance(TString OutFile)
{
	//
	// Write out data
	TFile* fout = new TFile(OutFile,"RECREATE");
	WriteAcceptance(fout);
	fout->Close();
	delete fout;
}
//
void AcceptanceClx::ReadAcceptance(TString InFile)
{
	//
	// Read in data
	TFile* f = new TFile(InFile, "READ");
	//
	// Import TTree
	TTree* T = (TTree*)f->Get("treeAcc");
	//
	// Get matrix
	TMatrixF* pAcc = new TMatrixF();
	T->SetBranchAddress("AcceptanceMatrix", &pAcc);
	T->GetEntry(0);
	//
	fNPtNodes = pAcc->GetNrows();
	fNThNodes = pAcc->GetNcols();
	fAcc.ResizeTo(fNPtNodes, fNThNodes);
	fAcc = *pAcc;
	//
	// Get Pt grid
	TVectorF* pPtArray = new TVectorF();
	T->SetBranchAddress("AcceptancePtVec", &pPtArray);
	T->GetEntry(0);
	fPtArray.ResizeTo(fNPtNodes);
	fPtArray = *pPtArray;
	//
	// Get Theta array
	TVectorF* pThArray = new TVectorF();
	T->SetBranchAddress("AcceptanceThVec", &pThArray);
	T->GetEntry(0);
	fThArray.ResizeTo(fNThNodes);
	fThArray = *pThArray;
	//
	std::cout << "AcceptanceClx::Read complete: Npt= " << fNPtNodes << ", Nth= " << fNThNodes << std::endl;
	//
	f->Close();
	delete f;
}
//
Double_t AcceptanceClx::HitNumber(Double_t pt, Double_t theta)
{
	//
	// Protect against values out of range
	Float_t eps = 1.0e-4;
	Float_t pt0 = (Float_t)pt;
	if (pt0 <= fPtArray(0)) pt0 = fPtArray(0) + eps;
	else if (pt0 >= fPtArray(fNPtNodes-1)) pt0 = fPtArray(fNPtNodes-1) - eps;
	Float_t th0 = (Float_t)theta;
	if (th0 <= fThArray(0)) th0 = fThArray(0) + eps;
	else if (th0 >= fThArray(fNThNodes - 1)) th0 = fThArray(fNThNodes - 1) - eps;
	//
	// Find cell
	Float_t* parray = fPtArray.GetMatrixArray();
	Int_t ip = TMath::BinarySearch(fNPtNodes, parray, pt0);
	Float_t* tarray = fThArray.GetMatrixArray();
	Int_t it = TMath::BinarySearch(fNThNodes, tarray, th0);
	//
	if (ip<0 || ip >fNPtNodes - 2)
	{
		std::cout << "Search error: (ip, pt) = (" << ip << ", " << pt << "), pt0 = " << pt0 << std::endl;
		std::cout << "Search error: pt nodes = " << fNPtNodes 
			<< " , last value = " << fPtArray(fNPtNodes - 1) << std::endl;
	}
	if (it<0 || ip >fNThNodes - 2)
	{
		std::cout << "Search error: (it, th) = (" << it << ", " << theta << "), th0 = " << th0 << std::endl;
		std::cout << "Search error: th nodes = " << fNThNodes
			<< " , last value = " << fThArray(fNThNodes - 1) << std::endl;
	}
	//
	// Bilinear interpolation
	//
	Double_t spt = (pt0 - fPtArray(ip)) / (fPtArray(ip + 1) - fPtArray(ip));
	Double_t sth = (th0 - fThArray(it)) / (fThArray(it + 1) - fThArray(it));
	Double_t A11 = fAcc(ip, it);
	Double_t A12 = fAcc(ip, it+1);
	Double_t A21 = fAcc(ip+1, it);
	Double_t A22 = fAcc(ip+1, it+1);
	Double_t A = A11 * (1 - spt) * (1 - sth) + A12 * (1 - spt) * sth +
		A21 * spt * (1 - sth) + A22 * spt * sth;
	//
	return A;
}
//
Double_t AcceptanceClx::HitNum(Double_t *x, Double_t *p) // Theta in degrees
{
	Double_t pt = x[0];
	Double_t th = x[1];
	//
	return HitNumber(pt, th);
}
//

