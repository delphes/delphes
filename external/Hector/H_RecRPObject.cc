  /* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                         *
*                   --<--<--  A fast simulator --<--<--     *
*                 / --<--<--     of particle   --<--<--     *
*  ----HECTOR----<                                          *
*                 \ -->-->-- transport through -->-->--     *
*                   -->-->-- generic beamlines -->-->--     *
*                                                           *
* JINST 2:P09005 (2007)                                     *
*      X Rouby, J de Favereau, K Piotrzkowski (CP3)         *
*       http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/// \file H_RecRPObject.cc
/// \brief Class performing the reconstruction based on forward detector measurements
///
/// Units : angles [rad], distances [m], energies [GeV], c=[1].

// local #includes
#include "H_RecRPObject.h"
#include "H_RomanPot.h"
#include "H_BeamParticle.h"
using namespace std;

// reconstruction class for forward detectors.
// Featuring the brand-new reco method from the 
// louvain group !

#define MEGA 1000000.

H_RecRPObject::H_RecRPObject(): emin(0), emax(-1), x1(0), x2(0), y1(0), y2(0), s1(0), s2(0), 
				txip(NOT_YET_COMPUTED), tyip(NOT_YET_COMPUTED), energy(NOT_YET_COMPUTED), q2(NOT_YET_COMPUTED), pt(NOT_YET_COMPUTED), 
				thebeam(new H_AbstractBeamLine()),
				f_1(new TF1("f_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				f_2(new TF1("f_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				g_1(new TF1("g_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				g_2(new TF1("g_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				d_1(new TF1("d_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				d_2(new TF1("d_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				k_1(new TF1("k_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				k_2(new TF1("k_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				l_1(new TF1("l_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				l_2(new TF1("l_2","[0] + [1]*x + [2]*x*x ",emin,emax)) 
{}

H_RecRPObject::H_RecRPObject(const float ss1, const float ss2, const H_AbstractBeamLine* beam) : emin(0), emax(-1), x1(0), x2(0), y1(0), y2(0), s1(ss1), s2(ss2), 
				txip(NOT_YET_COMPUTED), tyip(NOT_YET_COMPUTED), energy(NOT_YET_COMPUTED), q2(NOT_YET_COMPUTED), pt(NOT_YET_COMPUTED), 
				thebeam(beam->clone()),
				f_1(new TF1("f_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				f_2(new TF1("f_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				g_1(new TF1("g_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				g_2(new TF1("g_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				d_1(new TF1("d_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				d_2(new TF1("d_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				k_1(new TF1("k_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				k_2(new TF1("k_2","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				l_1(new TF1("l_1","[0] + [1]*x + [2]*x*x ",emin,emax)), 
				l_2(new TF1("l_2","[0] + [1]*x + [2]*x*x ",emin,emax))
	{if(ss1==ss2) cout<<"<H_RecRPObject> WARNING : detectors are on same position"<<endl;
}


H_RecRPObject::H_RecRPObject(const float ss1, const float ss2, const H_AbstractBeamLine& beam) : emin(0), emax(-1), x1(0), x2(0), y1(0), y2(0), s1(ss1), s2(ss2),
                                txip(NOT_YET_COMPUTED), tyip(NOT_YET_COMPUTED), energy(NOT_YET_COMPUTED), q2(NOT_YET_COMPUTED), pt(NOT_YET_COMPUTED),
                                thebeam(beam.clone()),
                                f_1(new TF1("f_1","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                f_2(new TF1("f_2","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                g_1(new TF1("g_1","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                g_2(new TF1("g_2","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                d_1(new TF1("d_1","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                d_2(new TF1("d_2","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                k_1(new TF1("k_1","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                k_2(new TF1("k_2","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                l_1(new TF1("l_1","[0] + [1]*x + [2]*x*x ",emin,emax)),
                                l_2(new TF1("l_2","[0] + [1]*x + [2]*x*x ",emin,emax))
        {if(ss1==ss2) cout<<"<H_RecRPObject> WARNING : detectors are on same position"<<endl;
}


H_RecRPObject::H_RecRPObject(const H_RecRPObject& r):
   emin(r.emin), emax(r.emax), x1(r.x1), x2(r.x2), y1(r.y1), y2(r.y2), s1(r.s1), s2(r.s2),
   txip(r.txip), tyip(r.tyip), energy(r.energy), q2(r.q2),  pt(r.pt), 
   //thebeam(r.thebeam->clone()),
   thebeam(new H_AbstractBeamLine(*(r.thebeam))),
   f_1(new TF1(*(r.f_1))),  f_2(new TF1(*(r.f_2))),  g_1(new TF1(*(r.g_1))),  g_2(new TF1(*(r.g_2))),
   d_1(new TF1(*(r.d_1))),  d_2(new TF1(*(r.d_2))),  k_1(new TF1(*(r.k_1))),  k_2(new TF1(*(r.k_2))),
   l_1(new TF1(*(r.l_1))),  l_2(new TF1(*(r.l_2)))
 {}

H_RecRPObject& H_RecRPObject::operator=(const H_RecRPObject& r) {
  if (this == &r) return *this;
   emin = r.emin, emax = r.emax;
   x1 = r.x1;   x2 = r.x2; 
   y1 = r.y1;   y2 = r.y2; 
   s1 = r.s1;   s2 = r.s2;
   txip= r.txip; tyip=r.tyip;
   energy= r.energy; q2= r.q2;  pt= r.pt;
   //thebeam = r.thebeam->clone();
   thebeam = new H_AbstractBeamLine(*(r.thebeam));
   f_1 = new TF1(*(r.f_1));
   f_2 = new TF1(*(r.f_2));
   g_1 = new TF1(*(r.g_1));
   g_2 = new TF1(*(r.g_2));
   d_1 = new TF1(*(r.d_1));
   d_2 = new TF1(*(r.d_2));
   k_1 = new TF1(*(r.k_1));
   k_2 = new TF1(*(r.k_2));
   l_1 = new TF1(*(r.l_1));
   l_2 = new TF1(*(r.l_2));
  return *this;
}

void H_RecRPObject::initialize() {
	// this method sets the functions that will be used for reco later
	// it should be used only once per beamline after the energy range was fixed. 
	// copying beamline and adding detectors
	
	if(emax<0) {
		cout<<"<H_RecRPObject> ERROR : energy range has to be set first !"<<endl;
		cout<<"<H_RecRPObject> Please run setERange() or computeERange()"<<endl;
		cout<<"<H_RecRPObject> initialization aborted"<<endl;
		return;
	}

	if(emax<emin) {
		cout<<"<H_RecRPObject> ERROR : maximum energy lower than minimum !"<<endl;
		cout<<"<H_RecRPObject> Please (re-)do setERange()"<<endl;
		cout<<"<H_RecRPObject> initialization aborted"<<endl;
		return;
	}

	H_AbstractBeamLine * b1 = thebeam->clone();
	H_RomanPot * rp1 = new H_RomanPot("rp1",s1,0);
	H_RomanPot * rp2 = new H_RomanPot("rp2",s2,0);
	b1->add(rp1);
	b1->add(rp2);

	// fitting parameters
	const int N = 20;
	double e_i[N];
	double f_1i[N], f_2i[N], g_1i[N], g_2i[N], d_1i[N], d_2i[N];
	double k_1i[N], k_2i[N], l_1i[N], l_2i[N];
	for(int i = 0; i < N; i++) {
		e_i[i] = emin + i * (emax - emin)/((double)N-1);
		// 
		// the bug seems to be linked to the delete operator of the TMatrixT class in root.
		// valgrind shows memory problems at that point.
		// for unknwown reasons, copying the matrix gets around this bug.
		// valgrind (related) ouptut :
		//
		// ==13029== Invalid read of size 4
		// ==13029==    at 0x5E75DB6: H_RecRPObject::initialize() (in /home/jdf/GGamma/Hector/lib/libHector.so)
		// ==13029==    by 0x804A2D8: intelligentreco_rpo(double, double, double, double, std::string, int) (H_IntelligentReco.cpp:314)
		// ==13029==    by 0x804A937: main (H_IntelligentReco.cpp:542)
		// ==13029==  Address 0x64AC740 is 80 bytes inside a block of size 144 free'd
		// ==13029==    at 0x4021D18: operator delete[](void*) (vg_replace_malloc.c:256)
		// ==13029==    by 0x5651F57: TMatrixT<float>::Delete_m(int, float*&) (in /home/jdf/root/5.12/lib/libMatrix.so)
		// ==13029==    by 0x565CC5D: TMatrixT<float>::~TMatrixT() (in /home/jdf/root/5.12/lib/libMatrix.so)
		// ==13029==    by 0x5E75D25: H_RecRPObject::initialize() (in /home/jdf/GGamma/Hector/lib/libHector.so)
		//
		//
		TMatrix el_mattt(b1->getPartialMatrix("rp1",e_i[i],MP,QP));
		const float *el_mat1 = el_mattt.GetMatrixArray();
		// 
		// conclusion : first line of el_mat1 is completely messed-up if it is taken directly from 
		// the return of getpartialmatrix like this : 
		// const float *el_mat1 = (b1->getPartialMatrix("rp1",e_i[i],MP,QP)).GetMatrixArray()
		//
		// long-term solution (apart noticing root staff) : replacing el_mat1[i] by the equivalent el_mattt(j,k)
		// which is anyway more transparent for the reader.
		//
		f_1i[i] = el_mat1[0];
		g_1i[i] = el_mat1[1*MDIM];
		d_1i[i] = MEGA*el_mat1[4*MDIM];
		k_1i[i] = el_mat1[2*MDIM+2];
		l_1i[i] = el_mat1[3*MDIM+2];
		const float* el_mat2 = (b1->getPartialMatrix("rp2",e_i[i],MP,QP).GetMatrixArray());
		f_2i[i] = el_mat2[0];
		g_2i[i] = el_mat2[1*MDIM];
		d_2i[i] = MEGA*el_mat2[4*MDIM];
		k_2i[i] = el_mat2[2*MDIM+2];
		l_2i[i] = el_mat2[3*MDIM+2];
	}
	TGraph gf_1(N,e_i,f_1i);
	TGraph gg_1(N,e_i,g_1i);
	TGraph gd_1(N,e_i,d_1i);
	TGraph gf_2(N,e_i,f_2i);
	TGraph gg_2(N,e_i,g_2i);
	TGraph gd_2(N,e_i,d_2i);
	TGraph gk_1(N,e_i,k_1i);
	TGraph gl_1(N,e_i,l_1i);
	TGraph gk_2(N,e_i,k_2i);
	TGraph gl_2(N,e_i,l_2i);

    // functions get their final shape	
	gf_1.Fit("f_1","Q");
	gg_1.Fit("g_1","Q");
	gd_1.Fit("d_1","Q");
	gf_2.Fit("f_2","Q");
	gg_2.Fit("g_2","Q");
	gd_2.Fit("d_2","Q");
	gk_1.Fit("k_1","Q");
	gl_1.Fit("l_1","Q");
	gk_2.Fit("k_2","Q");
	gl_2.Fit("l_2","Q");

	// cleaning the rest
	delete b1;

	// the end
	return;
}

void H_RecRPObject::setDetPos(const float ss1, const float ss2) {
	energy = NOT_YET_COMPUTED;
	s1 = ss1;
	s2 = ss2;
	if(ss1==ss2) cout<<"<H_RecRPObject> WARNING : detectors are on same position"<<endl;
	return;
}

void H_RecRPObject::setPositions(const float xx1, const float xx2, const float yy1, const float yy2) {
	energy = NOT_YET_COMPUTED;
	x1 = xx1;
	x2 = xx2;
	y1 = yy1;
	y2 = yy2;
	return;
}

void H_RecRPObject::setPosition_det1(const float xx1, const float yy1) {
	energy = NOT_YET_COMPUTED;
	x1 = xx1;
	y1 = yy1;
}

void H_RecRPObject::setPosition_det2(const float xx2, const float yy2) {
	energy = NOT_YET_COMPUTED;
	x2 = xx2;
	y2 = yy2;
}

void H_RecRPObject::setERange(const float eemin, const float eemax) {
	energy = NOT_YET_COMPUTED;
	emin = eemin;
	emax = eemax;
	f_1->SetRange(emin,emax);
	f_2->SetRange(emin,emax);
	g_1->SetRange(emin,emax);
	g_2->SetRange(emin,emax);
	d_1->SetRange(emin,emax);
	d_2->SetRange(emin,emax);
	k_1->SetRange(emin,emax);
	k_2->SetRange(emin,emax);
	l_1->SetRange(emin,emax);
	l_2->SetRange(emin,emax);
    return;
}


void H_RecRPObject::computeERange() {
	// optional method to determine the energy range of the FIRST detector 
	// in order to refine the fits and get maximum precision.
	energy = NOT_YET_COMPUTED;
	H_AbstractBeamLine * b1 = thebeam->clone();
	H_RomanPot * rp1 = new H_RomanPot("rp1",s1,0);
	b1->add(rp1);
	float max = 1;
	// number of energies to check 
	const int N = 1000;
	for(int i=0; i<N; i++) {
		H_BeamParticle p;
		p.setE(BE - (emin + i*(BE-emin)/((float)N)));
		p.computePath(b1);
		if(p.stopped(b1)) {
			if(p.getStoppingElement()->getName()=="rp1") {
				max = emin + i*(BE-emin)/((float)N);
			}
		}
	}
	cout<<"<H_RecRPObject> Valid energy losses run from 0 (default) to "<<max+20.<<" GeV"<<endl;
	setERange(0,max+20.);
	delete b1;
	return;
}

void H_RecRPObject::computeAll() {
	// The big game :
	// computing E, tx, ty, Q2 and Pt and filling the variables.
	//
	// The root TF1 class features nice bugs, which explains the
	// seemingly-dumb structures happening sometimes here as
	// workarounds for these bugs. The overall thing works very
	// well but will be cleaned later anyway.

	if(energy!=NOT_YET_COMPUTED) {
		cout<<"<H_RecRPObject> already computed variables, skipping ..."<<endl;
		return;
	}

	TF1 par0("par0","[0]",emin,emax);
	par0.SetParameter(0,-x1);
	TF1 par2("par2","[0]",emin,emax);
	par2.SetParameter(0,-y1);
	TF1 par1("par1","[0]",emin,emax);
	par1.SetParameter(0,-x2);
	TF1 par3("par3","[0]",emin,emax);
	par3.SetParameter(0,-y2);

	// angle compensating method :
	TF1 xx_E("xx_E","(g_2*(par0-d_1*x)-g_1*(par1-d_2*x))/(f_2*g_1-f_1*g_2)",emin,emax);
	TF1 yy_E("yy_E","(par2*l_2 - par3*l_1) / (k_2*l_1 - k_1*l_2)",emin,emax);
	TF1 xp_E("xp_E","(f_2*(par0-d_1*x)-f_1*(par1-d_2*x))/(g_2*f_1-g_1*f_2)",emin,emax);
	TF1 yp_E("yp_E","(par2*k_2-par3*k_1)/(l_2*k_1-l_1*k_2)",emin,emax);
	// it is possible to refine study using y info, but effect was not tested.	
	// TF1 p_xy_E("p_xy_E","(-xx_E*xx_E-yy_E*yy_E)",emin,emax);
	TF1 p_xy_E("p_xy_E","(-xx_E*xx_E)",emin,emax);

	energy = p_xy_E.GetMaximumX(emin,emax);
	txip = xp_E.Eval(energy);
	tyip = yp_E.Eval(energy);
	pt = sqrt(BE*(BE-energy)*(txip*txip+tyip*tyip)/(MEGA*MEGA));

	return;
}

float H_RecRPObject::getE(int ) {
// put for backward compatibility
        if(energy==NOT_YET_COMPUTED) { computeAll(); };
        return energy;
}  // to be removed !!!!!

float H_RecRPObject::getE() {
	if(energy==NOT_YET_COMPUTED) { computeAll(); };
	return energy;
}

float H_RecRPObject::getTX() {
	if(energy==NOT_YET_COMPUTED) { computeAll(); };
	return txip;
}

float H_RecRPObject::getTY() {
	if(energy==NOT_YET_COMPUTED) { computeAll(); };
	return tyip;
}

float H_RecRPObject::getQ2() {
	if(energy==NOT_YET_COMPUTED) { computeAll(); };
	cout<<"<H_RecRPObject::getQ2> Not implemented yet"<<endl;
	return 0;
}

float H_RecRPObject::getPt() {
	if(energy==NOT_YET_COMPUTED) { computeAll(); };
	return pt;
}

std::ostream& operator<< (std::ostream& os, const H_RecRPObject& rp) {
	os << "e_min=" <<  rp.emin << "\t e_max= " << rp.emax << endl;
	os << "x1="    << rp.x1 << "\t x2= " << rp.x2 << "\t y1=" << rp.y1 << "\t y2=" << rp.y2 
	   << "\t s1=" << rp.s1 << "\t s2=" << rp.s2 << endl;
	os << "txip=" << rp.txip << "\t tyip=" << rp.tyip << "\t energy=" << rp.energy << "\t q2=" << rp.q2 << "\t pt=" << rp.pt << endl;
   return os;
}

