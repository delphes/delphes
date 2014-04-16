/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_RecRPObject.cc
/// \brief Classes aiming at reconstruction particle properties

// C++ #includes
#include <iostream>
#include <iomanip>

// ROOT #includes
//#include "TGraph.h"
//#include "TF1.h"
//#include "TCanvas.h"

// local #includes
#include "H_RecRPObject.h"
#include "H_RomanPot.h"
#include "H_BeamParticle.h"
using namespace std;

H_RecRPObject::H_RecRPObject() {
	x1 = 0.;
	x2 = 0.;
	y1 = 0.;
	y2 = 0.;
	s1 = 0.;
	s2 = 0.;
	corr1_TM = 0;
	corr2_TM = 0;
	corr1_AM = 0;
	corr2_AM = 0;
	thx = NOT_YET_COMPUTED;
	thy = NOT_YET_COMPUTED;
	x0  = NOT_YET_COMPUTED;
	y0  = NOT_YET_COMPUTED;
	energy = NOT_YET_COMPUTED;
	virtuality = NOT_YET_COMPUTED;
	matrp1 = new TMatrix(MDIM,MDIM);
	matrp2 = new TMatrix(MDIM,MDIM);
	thebeam = new H_AbstractBeamLine();
}

H_RecRPObject::H_RecRPObject(const float S1, const float S2, const H_AbstractBeamLine& beamline) {
	x1 = 0;
	x2 = 0;
	y1 = 0;
	y2 = 0;
	s1 = S1;
	s2 = S2;
	thx = NOT_YET_COMPUTED;
	thy = NOT_YET_COMPUTED;
	x0  = NOT_YET_COMPUTED;
	y0  = NOT_YET_COMPUTED;
	energy = NOT_YET_COMPUTED;
	virtuality = NOT_YET_COMPUTED;
//	matrp1 = new TMatrix(MDIM,MDIM);
//	matrp2 = new TMatrix(MDIM,MDIM);
	thebeam = new H_AbstractBeamLine(beamline);
	H_RomanPot * rp1 = new H_RomanPot("rp1",s1,0);
	thebeam->add(rp1);
	H_RomanPot * rp2 = new H_RomanPot("rp2",s2,0);
	thebeam->add(rp2);
	matrp1 = new TMatrix(*(thebeam->getPartialMatrix("rp1",0.,MP,QP)));
	matrp2 = new TMatrix(*(thebeam->getPartialMatrix("rp2",0.,MP,QP)));

	corr1_TM = getECorrectionFactor(0,TM);
	corr2_TM = getECorrectionFactor(1,TM);
	corr1_AM = getECorrectionFactor(0,AM);
	corr2_AM = getECorrectionFactor(1,AM);
//	cout << corr1_TM << " " << corr2_TM << endl;
//	cout << corr1_AM << " " << corr2_AM << endl;
}

H_RecRPObject::H_RecRPObject(const H_RecRPObject& r) {
	x1 = r.x1;
	x2 = r.x2;
	y1 = r.y1;
	y2 = r.y2;
	s1 = r.s1;
	s2 = r.s2;
	x0 = r.x0;
	y0 = r.y0;
	thx = r.thx;
	thy = r.thy;
	energy = r.energy;
	virtuality = r.virtuality;
	matrp1 = new TMatrix(*(r.matrp1));
	matrp2 = new TMatrix(*(r.matrp2));
	corr1_TM = r.corr1_TM;
	corr2_TM = r.corr2_TM;
	corr1_AM = r.corr1_AM;
	corr2_AM = r.corr2_AM;
	thebeam = new H_AbstractBeamLine(*(r.thebeam));
}

H_RecRPObject& H_RecRPObject::operator=(const H_RecRPObject& r) {
	if(this==&r) return *this;
	x1 = r.x1;
	x2 = r.x2;
	y1 = r.y1;
	y2 = r.y2;
	s1 = r.s1;
	s2 = r.s2;
	x0 = r.x0;
	y0 = r.y0;
	thx = r.thx;
	thy = r.thy;
	energy = r.energy;
	virtuality = r.virtuality;
	matrp1 = new TMatrix(*(r.matrp1));
	matrp2 = new TMatrix(*(r.matrp2));
	corr1_TM = r.corr1_TM;
	corr2_TM = r.corr2_TM;
	corr1_AM = r.corr1_AM;
	corr2_AM = r.corr2_AM;
	thebeam = new H_AbstractBeamLine(*(r.thebeam));
	return *this;
}

float H_RecRPObject::getECorrectionFactor(const unsigned int facn, const unsigned int method) {
/*
 * commented out because CMSSW does not want any TGraph/TCanvas/TFit
 *
 * to be fixed !
 * X.R. 07/05/2009
 *
	float beta1 = ((thebeam->getPartialMatrix("rp1",0,MP,QP))->GetMatrixArray())[1*MDIM];
	float beta2 = ((thebeam->getPartialMatrix("rp2",0,MP,QP))->GetMatrixArray())[1*MDIM];
	float disp1 = ((thebeam->getPartialMatrix("rp1",0,MP,QP))->GetMatrixArray())[4*MDIM]*URAD;
	float disp2 = ((thebeam->getPartialMatrix("rp2",0,MP,QP))->GetMatrixArray())[4*MDIM]*URAD;
	const int n = 20; //using 20 points to get a good quadratic fit
	float ee[n], rece[n];

	for(int i = 0; i < n; i++) {
		ee[i] = 10 + i*200./((float)n-1);
		H_BeamParticle p1;
		p1.emitGamma(ee[i],0.);
		p1.computePath(thebeam,1);
		p1.propagate(s1);
		float x1 = p1.getX();
		p1.propagate(s2);
		float x2 = p1.getX();
		switch (method) {
			case TM: { rece[i] = -x1/disp1; }; break;
			case AM: { rece[i] = -(beta2*x1-beta1*x2)/(beta2*disp1-beta1*disp2); }; break;
			case PM: { rece[i] = -x1/disp1; cout<<"this method has not been implemented, using trivial reconstruction"<<endl; } break;
			default: { rece[i] = -(beta2*x1-beta1*x2)/(beta2*disp1-beta1*disp2); }; break;
		}
		ee[i] = ee[i] - rece[i];
	}
	char mytitle[50];
	sprintf(mytitle,"c_%d",method);
//	TCanvas*c = new TCanvas();
//	c->SetTitle(mytitle);
	TGraph* g1 = new TGraph(n,rece,ee);
	TF1* fit1 = new TF1("fit1","[0]*x + [1]*x*x",0,100);
	g1->Fit("fit1","Q");
	float xfact = fit1->GetParameter(facn);
//	g1->Draw("APL");
	delete g1;
	delete fit1;

	return xfact;
*/
        return 1.;
}


void H_RecRPObject::setPositions(const float X1, const float Y1, const float X2, const float Y2) {
	thx = NOT_YET_COMPUTED;
	thy = NOT_YET_COMPUTED;
	x0  = NOT_YET_COMPUTED;
	y0  = NOT_YET_COMPUTED;
	energy = NOT_YET_COMPUTED;
	virtuality = NOT_YET_COMPUTED;
	x1 = X1;
	x2 = X2;
	y1 = Y1;
	y2 = Y2;
	
	return;
}

void H_RecRPObject::printProperties() const {
	cout << "Roman pot variables :" << endl;
	cout << "\t pot 1 : (x,y,s) = (" << x1 << " , " << y1 << " , " << s1 << " )" << endl;
	cout << "\t pot 2 : (x,y,s) = (" << x2 << " , " << y2 << " , " << s2 << " )" << endl;
	cout << endl << "Reconstructed variables :" << endl;
	cout << "\t IP : (x,y) = (" << x0 << " , " << y0 << ") and (theta_x, theta_y) = (" << thx << " , " << thy << " )" << endl;
	if (energy==NOT_YET_COMPUTED) cout << "\t Energy not yet computed" << endl;
	else cout << "\t Energy = " << energy << " GeV" << endl;	
	if (virtuality==NOT_YET_COMPUTED) cout << "\t Virtuality not yet computed" << endl;
	else cout << "\t Virtuality = " << virtuality << " GeV^2" << endl;
	cout << endl;
}

float H_RecRPObject::getE() {
	if(energy==NOT_YET_COMPUTED) {
		cout<<"Please first compute energy using your favourite method"<<endl;
		return NOT_YET_COMPUTED;
	}
	return energy;
}

float H_RecRPObject::getE(const unsigned int method) {
	switch (method) {
		case TM: {energy = computeE_TM();} break;
		case AM: {energy = computeE_AM();} break;
		case PM: {energy = computeE_PM();} break;
		default: {energy = computeE_AM();} break;
	};
	return energy;
}

float H_RecRPObject::computeX0() {
	if(energy==NOT_YET_COMPUTED) {
		cout<<"Please first compute energy using your favourite method"<<endl;
		return NOT_YET_COMPUTED;
	}
	float alpha1 = (matrp1->GetMatrixArray())[0];
	float alpha2 = (matrp2->GetMatrixArray())[0];
	float disp1 = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
	float disp2 = (matrp2->GetMatrixArray())[4*MDIM]*URAD;
	x0 = (disp2*x1-disp1*x2)/(disp2*alpha1-disp1*alpha2);
	return x0;
}

float H_RecRPObject::computeY0() {
	if(energy==NOT_YET_COMPUTED) {
		cout<<"Please first compute energy using your favourite method"<<endl;
		return NOT_YET_COMPUTED;
	}
	float gamma1 = (matrp1->GetMatrixArray())[2*MDIM+2];
	float gamma2 = (matrp2->GetMatrixArray())[2*MDIM+2];
	float delta1 = (matrp1->GetMatrixArray())[3*MDIM+2];
	float delta2 = (matrp2->GetMatrixArray())[3*MDIM+2];
	y0 = (delta2*y1-delta1*y2)/(delta2*gamma1-delta1*gamma2);
	return y0;
}

float H_RecRPObject::computeTX() {
	if(energy==NOT_YET_COMPUTED) {
	        cout<<"Please first compute energy using your favourite method"<<endl;
	        return NOT_YET_COMPUTED;
	} 
	float beta1 = (matrp1->GetMatrixArray())[1*MDIM];
	float beta2 = (matrp2->GetMatrixArray())[1*MDIM];
	float disp1 = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
	float disp2 = (matrp2->GetMatrixArray())[4*MDIM]*URAD;
	// computes thx (murad)
	thx = (x1*disp2-x2*disp1)/(beta1*disp2-beta2*disp1);
	return thx;
}

float H_RecRPObject::computeTY() {
	if(energy==NOT_YET_COMPUTED) {
		cout<<"Please first compute energy using your favourite method"<<endl;
		return NOT_YET_COMPUTED;
	}
	float gamma1 = (matrp1->GetMatrixArray())[2*MDIM+2];
	float gamma2 = (matrp2->GetMatrixArray())[2*MDIM+2];
	float delta1 = (matrp1->GetMatrixArray())[3*MDIM+2];
	float delta2 = (matrp2->GetMatrixArray())[3*MDIM+2];
	// computes thy (murad)
	thy = (y1*gamma2-y2*gamma1)/(delta1*gamma2-delta2*gamma1);
	return thy;
}

float H_RecRPObject::computeE_TM() {
	// computes the emitted particle energy, from the trivial method
	float disp = ((thebeam->getPartialMatrix("rp1",0,MP,QP))->GetMatrixArray())[4*MDIM]*URAD;
	energy = -x1/disp;
	// corrects for nonlinear effects
	energy = (1+corr1_TM)*energy + corr2_TM*energy*energy;
	// sets the rp matrices at obtained energy
	delete matrp1;
	delete matrp2;
	matrp1 = new TMatrix(*(thebeam->getPartialMatrix("rp1",energy,MP,QP)));
	matrp2 = new TMatrix(*(thebeam->getPartialMatrix("rp2",energy,MP,QP)));
	// returns ...
	return energy;
}

float H_RecRPObject::computeE_AM() {
	// computes the emitted particle energy, from the angle compensation method, iterative way
	const int N = 10;
	delete matrp1;
	delete matrp2;
	matrp1 = new TMatrix(*(thebeam->getPartialMatrix("rp1",0,MP,QP)));
	float disp = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
	delete matrp1;
	energy = -x1/disp;
	matrp1 = new TMatrix(*(thebeam->getPartialMatrix("rp1",energy,MP,QP)));
	matrp2 = new TMatrix(*(thebeam->getPartialMatrix("rp2",energy,MP,QP)));
	float beta1 = (matrp1->GetMatrixArray())[1*MDIM];
	float beta2 = (matrp2->GetMatrixArray())[1*MDIM];
	float disp1 = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
	float disp2 = (matrp2->GetMatrixArray())[4*MDIM]*URAD;
	for(int i = 0; i < N; i++) {
		energy = -(beta2*x1-beta1*x2)/(beta2*disp1-beta1*disp2);
		delete matrp1;
		delete matrp2;
		matrp1 = new TMatrix(*(thebeam->getPartialMatrix("rp1",energy,MP,QP)));
		matrp2 = new TMatrix(*(thebeam->getPartialMatrix("rp2",energy,MP,QP)));
		beta1 = (matrp1->GetMatrixArray())[1*MDIM];
		beta2 = (matrp2->GetMatrixArray())[1*MDIM];
		disp1 = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
		disp2 = (matrp2->GetMatrixArray())[4*MDIM]*URAD;
	}
	// returns ...
	return energy;
}

float H_RecRPObject::computeE_PM() {
	cout<<"Not yet implemented, nothing done"<<endl;
	energy = NOT_YET_COMPUTED;
	return energy;
}

float H_RecRPObject::computeQ2() {
	// computes the emitted particle virtuality
	// energy should be teconstructed first
	if(energy==NOT_YET_COMPUTED) {
		cout<<"Please first compute energy using your favourite method"<<endl;
		return NOT_YET_COMPUTED;
	} 
	// getting parameters for reconstructed particle energy
    float beta1 = (matrp1->GetMatrixArray())[1*MDIM];
    float beta2 = (matrp2->GetMatrixArray())[1*MDIM];
    float gamma1 = (matrp1->GetMatrixArray())[2*MDIM+2];
    float gamma2 = (matrp2->GetMatrixArray())[2*MDIM+2];
    float delta1 = (matrp1->GetMatrixArray())[3*MDIM+2];
    float delta2 = (matrp2->GetMatrixArray())[3*MDIM+2];
    float disp1 = (matrp1->GetMatrixArray())[4*MDIM]*URAD;
    float disp2 = (matrp2->GetMatrixArray())[4*MDIM]*URAD;
    // angles reconstruction
    float rec_thx = (x1*disp2-x2*disp1)/(beta1*disp2-beta2*disp1)/URAD;
    float rec_thy = (y1*gamma2-y2*gamma1)/(delta1*gamma2-delta2*gamma1)/URAD;
	// Q² reconstruction
    virtuality = BE*(BE-energy)*(rec_thx*rec_thx+rec_thy*rec_thy);
	// returns ...
	return virtuality;
}

