/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_TransportMatrices.cc
/// \brief Includes the implementation of every transport matrix.

// c++ #includes
#include <iostream>

// C #includes
#include <cmath>

// local #includes
#include "H_Parameters.h"
#include "H_TransportMatrices.h"
using namespace std;

bool relative_energy = 1;

// caution : do not change particle mass, not implemented yet.

extern double omega(const double k, const double l) {
		// [l] = [m] and [k] = [1/m�] for quadrupoles
		// [omega] = [1]
	return sqrt(fabs(k))*l;
}

extern double radius(const double k) {
		// [k] = [1/m�] for quadrupoles
		// [k] = [1/m]  for dipoles
		// [radius(k)] = [m]
	if(k==0 && VERBOSE) cout<<"ERROR : Dipole has no effect : results will be corrupted"<<endl; 
	// this is protected by the "if(k==0) -> driftmat" in every matrix below (ex vquatmat)
	return (k==0) ? 1 : 1/k;
}

extern void printMatrix(TMatrix * TMat) {
	char temp[20];
	float * el = new float[MDIM*MDIM];
	el = (TMat->GetMatrixArray());

	cout << endl << "\t";
	for(int i=0;i<MDIM*MDIM;i++) {
		if (el[i]<0)
  			 {sprintf(temp,"%.5e",el[i]);}
		else if (el[i]>0)
   			 {sprintf(temp," %.5e",el[i]);}
   		else {sprintf(temp,"       0    ");}
   		
		cout << temp << " ";
		if((i+1)%MDIM == 0) { cout << endl << "\t"; }
	}
	cout << endl;
}

extern TMatrix vquadmat(const float l, const float k, const float eloss = 0., const float p_mass=MP, const float p_charge=QP) {
		// the length l is in [m]
		// the strength k is in [1/m�] for quadrupoles
		// eloss in [GeV]
		// ke is the modified field with respect to the eloss
                //  k = e/p * dB/dx with p = mv (and m = MP)
                //  k -> ke = k * p/ (p - dp) <- chromacity
                //  ke -> ke * p_charge / QP <- if not a proton
                //  ke = 0 if charge = 0, whatever the mass
                
	const double p0 = sqrt( (BE-MP)*(BE+MP) );
	const double E  = BE - eloss;
	const double p  = sqrt( (E-p_mass)*(E+p_mass) );
	const float ke = (p_charge==0) ? 0 : k* p0/p  *p_charge/QP;
	if (ke==0) {
		TMatrix drift(driftmat(l));		
		return drift;
	}
	// else... : 
	float om = omega(ke,l);
	float * mat = new float[MDIM*MDIM];
	float tmat[MDIM*MDIM] =  {cosh(om),sqrt(ke)*sinh(om),0.,0., 0.,0.,
                           (1/sqrt(ke))*sinh(om),cosh(om),0.,0., 0.,0.,
                           0.,0.,cos(om),-sqrt(ke)*sin(om), 0.,0.,
                           0.,0.,(1/sqrt(ke))*sin(om),cos(om), 0.,0.,
                           0., 0., 0., 0., 1., 0.,
			   0., 0., 0., 0., 0., 1.
	};
	for (int i=0; i<MDIM*MDIM; i++) {mat[i] = tmat[i];}
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
    return TMat;
}

extern TMatrix hquadmat(const float l, const float k, const float eloss = 0., const float p_mass=MP, const float p_charge=QP) {
		// the length l is in [m]
		// the strength k is in [1/m�] for quadrupoles
		// ke is the modified field with respect to the eloss
                //  k = e/p * dB/dx with p = mv (and m = MP)
                //  k -> ke = k * p/ (p- dp) <- chromacity
                //  ke -> ke *p_charge/QP <- if not a proton
        const double p0 = sqrt( (BE-MP)*(BE+MP) );
        const double E  = BE - eloss;
        const double p  = sqrt( (E-p_mass)*(E+p_mass) );
        const float ke = (p_charge==0) ? 0 : fabs(k* p0/p)  *p_charge/QP;

        if (ke==0) {
		TMatrix drift(driftmat(l));
		return drift;
	}
	float om = omega(ke,l);
	float * mat = new float[MDIM*MDIM];
	float tmat[MDIM*MDIM] =  {cos(om),-sqrt(ke)*sin(om),0.,0., 0., 0.,
                           (1/sqrt(ke))*sin(om),cos(om),0.,0., 0., 0.,
                           0.,0.,cosh(om),sqrt(ke)*sinh(om), 0., 0., 
                           0.,0.,(1/sqrt(ke))*sinh(om),cosh(om), 0., 0.,
                    	    0., 0., 0., 0., 1., 0.,
                           0., 0., 0., 0., 0., 1. 
	};
	for(int i=0;i<MDIM*MDIM;i++) { mat[i] = tmat[i]; }
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
    return TMat;
}

extern TMatrix rdipmat(const float l, const float k, const float eloss = 0., const float p_mass=MP, const float p_charge=QP) {
		// the length l is in [m]
		// the strength k is in [1/m] for dipoles
		// ke is the modified field with respect to the eloss
                //  k = e/p * dB/dx with p = mv (and m = MP)
                //  k -> ke = k * p/ (p- dp) <- chromacity
                //  ke -> ke * q_mass/QP <- if not a proton

        const double p0 = sqrt( (BE-MP)*(BE+MP) );
        const double E  = BE - eloss;
        const double p  = sqrt( (E-p_mass)*(E+p_mass) );
        const float ke = (p_charge==0) ? 0 : k* p0/p  *p_charge/QP;
 
        if (ke==0) {
 		TMatrix drift(driftmat(l));
		return drift;
	}
	float r = radius(ke);
	float * mat = new float[MDIM*MDIM];
	float * efmat = new float[MDIM*MDIM];
	double simp = r*2*sin(l/(2*r))*sin(l/(2*r))/BE;
	double psy = ke*l/2.;
	float tefmat[MDIM*MDIM] = {1., (float)(tan(psy)*ke), 0., 0., 0., 0.,
	                            0., 1., 0., 0., 0., 0.,
	                            0., 0., 1., (float)(-tan(psy)*ke), 0., 0.,
	                            0., 0., 0., 1., 0., 0.,
	                            0., 0., 0., 0., 1., 0.,
	                            0., 0., 0., 0., 0., 1. };
 
	float tmat[MDIM*MDIM] =  {cos(l/r),(-1/r)*sin(l/r),0.,0., 0., 0.,
	             	           r*sin(l/r),cos(l/r),0.,0., 0., 0.,
   		           	           0.,0.,1.,0., 0., 0.,
   	            	           0.,0.,l,1., 0., 0., 
   	                        (float)simp, (float)(sin(l/r)/BE), 0., 0., 1., 0.,
   	                        0., 0., 0., 0., 0., 1. };
	for(int i=0;i<MDIM*MDIM;i++) { 
		mat[i] = tmat[i];
	    efmat[i] = tefmat[i];	
	}
	TMatrix TMat(MDIM,MDIM,mat);
	TMatrix TEfmat(MDIM,MDIM,efmat);
    if(relative_energy) {
        TMat *= TEfmat;
        TEfmat *= TMat;
    }

//	if(VERBOSE) cout<<"\t WARNING : RDipoles not implemented and replaced by SDipoles" << endl;
	delete [] mat;
	delete [] efmat;
    if(relative_energy) {
        return TEfmat;
    } else {
        return TMat;
    }
}

extern TMatrix sdipmat(const float l, const float k, const float eloss = 0., const float p_mass=MP, const float p_charge=QP) {
		// the length l is in [m]
		// the strength k is in [1/m] for dipoles
		// ke is the modified field with respect to the eloss
                //  k = e/p * dB/dx with p = mv (and m = MP)
                //  k -> ke = k * p/ (p- dp) <- chromacity
                //  ke -> ke * q_mass/QP <- if not a proton
                
        const double p0 = sqrt( (BE-MP)*(BE+MP) );
        const double E  = BE - eloss;
        const double p  = sqrt( (E-p_mass)*(E+p_mass) );
        const float ke = (p_charge==0) ? 0 : k* p0/p  *p_charge/QP;

        if (ke==0) {
		TMatrix drift(driftmat(l));
		return drift;
	}
	extern bool relative_energy;
	float r = radius(ke);
  	float * mat = new float[MDIM*MDIM];

	float simp = 2*r*sin(l/(2*r))*sin(l/(2*r))/BE;
	float tmat[MDIM*MDIM] =  {cos(l/r),(-1/r)*sin(l/r),0.,0., 0., 0.,
		                      r*sin(l/r),cos(l/r),0.,0., 0., 0.,
			                  0.,0.,1.,0., 0., 0.,
			                  0.,0.,l,1., 0., 0.,
			                  simp, (float)(sin(l/r)/BE), 0., 0., 1., 0.,
			                  0., 0., 0., 0., 0., 1.
				           };
	if(!relative_energy) {
		tmat[24] = 0;
		tmat[25] = 0;
	} 

	for(int i=0;i<MDIM*MDIM;i++) { mat[i] = tmat[i]; }
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
    return TMat;
}

extern TMatrix driftmat(const float l) {
		// the length l is in [m]
	float * mat = new float[MDIM*MDIM];
	float tmat[MDIM*MDIM] =  {1.,0.,0.,0.,0.,0.,
                           l ,1.,0.,0.,0.,0.,
                           0.,0.,1.,0.,0.,0.,
                           0.,0.,l ,1.,0.,0.,
                           0.,0.,0.,0.,1.,0.,
                           0., 0., 0., 0., 0., 1. 
	};
	for(int i=0;i<MDIM*MDIM;i++) { mat[i] = tmat[i]; }
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
	return TMat;
}


extern TMatrix hkickmat(const float l, const float k, const float eloss =0., const float p_mass=MP, const float p_charge=QP) {
                // the length l is in [m]
                // the strength k is in [rad]
        const double p0 = sqrt( (BE-MP)*(BE+MP) );
        const double E  = BE - eloss;
        const double p  = sqrt( (E-p_mass)*(E+p_mass) );
        const float ke = (p_charge==0) ? 0 : -k* p0/p  *p_charge/QP;

        if (ke==0) {
		TMatrix drift(driftmat(l));
		return drift;
	}
	float * mat = new float[MDIM*MDIM];
	float tmat[MDIM*MDIM] =  {1.,0.,0.,0.,0.,0.,
                           l ,1.,0.,0.,0.,0.,
                           0.,0.,1.,0.,0.,0.,
                           0.,0.,l ,1.,0.,0.,
                           0.,0.,0.,0.,1.,0.,
                           (float)(l*tan(ke)/2.),ke, 0., 0., 0., 1. 
	};

	for(int i=0;i<MDIM*MDIM;i++) { mat[i] = tmat[i]; }
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
	return TMat;
}

extern TMatrix vkickmat(const float l, const float k, const float eloss=0., const float p_mass=MP, const float p_charge=QP) {
                // the length l is in [m]
                // the strength k is in [rad]
        const double p0 = sqrt( (BE-MP)*(BE+MP) );
        const double E  = BE - eloss;
        const double p  = sqrt( (E-p_mass)*(E+p_mass) );
        const float ke = (p_charge==0) ? 0 : -k* p0/p  *p_charge/QP;

        if (ke==0) {
		TMatrix drift(driftmat(l));
		return drift;
	}
	float * mat = new float[MDIM*MDIM];
	float tmat[MDIM*MDIM] =  {1.,0.,0.,0.,0.,0.,
                           l ,1.,0.,0.,0.,0.,
                           0.,0.,1.,0.,0.,0.,
                           0.,0.,l ,1.,0.,0.,
                           0.,0.,0.,0.,1.,0.,
                           0.,0.,(float)(l*tan(ke)/2.),ke, 0., 1. 
	};

	for(int i=0;i<MDIM*MDIM;i++) { mat[i] = tmat[i]; }
	TMatrix TMat(MDIM,MDIM,mat);
	delete [] mat;
	return TMat;
}

