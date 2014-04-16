/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Beam.cc
/// \brief Describes a set a particles as a beam
///

// ROOT #includes
//#include "TGraph.h"
#include "TRandom.h"

// local #includes
#include "H_Beam.h"
using namespace std;

H_Beam::H_Beam() {
        setPosition(PX,PY,TX+CRANG,TY,PS);
        setE(BE);
        setDispersion(SX,SY,STX,STY,SS);
        setDE(SBE);
	Nparticles=0;
}

H_Beam::H_Beam(const H_Beam& be) {
	beamParticles = be.beamParticles;
	setPosition(be.fx_ini,be.fy_ini,tx_ini,ty_ini,be.fs_ini);
	setE(be.fe_ini);
	setDispersion(be.x_disp,be.y_disp,be.tx_disp,be.ty_disp,be.s_disp);
	setDE(be.e_disp);
	Nparticles = be.Nparticles;
}

H_Beam& H_Beam::operator=(const H_Beam& be) {
	if(this==&be) return *this;
        beamParticles = be.beamParticles;
        setPosition(be.fx_ini,be.fy_ini,tx_ini,ty_ini,be.fs_ini);
        setE(be.fe_ini);
        setDispersion(be.x_disp,be.y_disp,be.tx_disp,be.ty_disp,be.s_disp);
        setDE(be.e_disp);
        Nparticles = be.Nparticles; 
        return *this;
}

H_Beam::~H_Beam() {
	beamParticles.clear(); 
	return;
};

void H_Beam::createBeamParticles(const unsigned int Number_of_particles) {
	createBeamParticles(Number_of_particles,MP,QP);
}

void H_Beam::createBeamParticles(const unsigned int Number_of_particles, const double p_mass, const double p_charge) {
	beamParticles.clear();
	Nparticles = (Number_of_particles<1) ? 1 : Number_of_particles;
	for (unsigned int i=0; i<Nparticles; i++) {
		H_BeamParticle p(p_mass,p_charge);
		p.setPosition(fx_ini,fy_ini,tx_ini,ty_ini,fs_ini);
		p.setE(fe_ini);
		p.smearPos(x_disp,y_disp);
		p.smearAng(tx_disp,ty_disp);
		p.smearE(e_disp);
		p.smearS(s_disp);
		if (VERBOSE) {if (i==0) cout << " x_ini , tx_ini " << p.getX() << " " << p.getTX() << endl;}
		beamParticles.push_back(p);
	}
}

void H_Beam::createXScanningBeamParticles(const unsigned int Number_of_particles, const float fx_max) {
        beamParticles.clear();
	Nparticles = (Number_of_particles<2) ? 2 : Number_of_particles;
        for (unsigned int i=0; i<Nparticles; i++) {
                H_BeamParticle p;
		float fx = fx_ini + i/(float)(Nparticles-1) * (fx_max-fx_ini);
                p.setPosition(fx,fy_ini,0,0,fs_ini);
                p.setE(fe_ini);
                beamParticles.push_back(p);
        }
}

void H_Beam::createYScanningBeamParticles(const unsigned int Number_of_particles, const float fy_max) {

        beamParticles.clear();
	Nparticles = (Number_of_particles<2) ? 2 : Number_of_particles;
        for (unsigned int i=0; i<Nparticles; i++) {
                H_BeamParticle p;
                float fy = fy_ini + i/(float)(Nparticles-1) * (fy_max-fy_ini);
                p.setPosition(fx_ini,fy,0,0,fs_ini);
                p.setE(fe_ini);
                beamParticles.push_back(p);
        }
}

void H_Beam::createTXScanningBeamParticles(const unsigned int Number_of_particles, const float tx_max) {
        beamParticles.clear();
	Nparticles = (Number_of_particles<2) ? 2 : Number_of_particles;
        for (unsigned int i=0; i<Nparticles; i++) {
                H_BeamParticle p;
                float tx = tx_ini + i/(float)(Nparticles-1) * (tx_max-tx_ini);
                p.setPosition(fx_ini,fy_ini,tx,ty_ini,fs_ini);
                p.setE(fe_ini);
                beamParticles.push_back(p);
        }
}

void H_Beam::createTYScanningBeamParticles(const unsigned int Number_of_particles, const float ty_max) {
        beamParticles.clear(); 
        Nparticles = (Number_of_particles<2) ? 2 : Number_of_particles;
        for (unsigned int i=0; i<Nparticles; i++) { 
                H_BeamParticle p; 
                float ty = ty_ini + i/(float)(Nparticles-1) * (ty_max-ty_ini);
                p.setPosition(fx_ini,fy_ini,tx_ini,ty,fs_ini);
                p.setE(fe_ini); 
                beamParticles.push_back(p); 
        } 
}

const H_BeamParticle * H_Beam::getBeamParticle(const unsigned int particle_index) const {
//	const int N = (particle_index<0)?0:(( particle_index>Nparticles)?Nparticles:particle_index);
	const int N = (particle_index>Nparticles)?Nparticles:particle_index;
	return &(*(beamParticles.begin()+N));// same as "return &beamParticles[N];" but more efficient
}

H_BeamParticle * H_Beam::getBeamParticle(const unsigned int particle_index) {
//        const int N = (particle_index<0)?0:(( particle_index>Nparticles)?Nparticles:particle_index);
	const int N = (particle_index>Nparticles)?Nparticles:particle_index;
        return &(*(beamParticles.begin()+N));// same as "return &beamParticles[N];" but more efficient
}

void H_Beam::add(const H_BeamParticle &p) {
	beamParticles.push_back(p);
	Nparticles++;
}

void H_Beam::computePath(const H_AbstractBeamLine * beamline, const bool NonLinear) {
	vector<H_BeamParticle>::iterator particle_i;

	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		particle_i->computePath(beamline,NonLinear);	
	}
}

void H_Beam::computePath(const H_AbstractBeamLine * beamline) {
	computePath(beamline,false);
}

/// Propagates the beam until a given s
void H_Beam::propagate(const float position) {
	vector<H_BeamParticle>::iterator particle_i;
	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                particle_i->propagate(position);
        }
}

void H_Beam::emitGamma(const double gee, const double gq2) {
	/// @param gee = \f$ E_{\gamma} \f$ is the photon energy
	/// @param gq2 = \f$ Q^2 < 0 \f$ is virtuality of photon \f$ Q^{2} = E^{2}-\vec{k}^{2} \f$
	emitGamma(gee,gq2,0,2*PI);
}

void H_Beam::emitGamma(const double gee, const double gq2, const double phimin, const double phimax) {
	/// @param gee = \f$ E_{\gamma} \f$ is the photon energy
	/// @param gq2 = \f$ Q^2 < 0 \f$ is virtuality of photon \f$ Q^{2} = E^{2}-\vec{k}^{2} \f$
	/// @param phimin : lower bound for \f$ \phi \f$
	/// @param phimax : higher bound for \f$ \phi \f$
	vector<H_BeamParticle>::iterator particle_i;
	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++)
                particle_i->emitGamma(gee,gq2,phimin,phimax);
}

float H_Beam::getBetaX(const float s, float& error_on_betax) {
	/// @param s is the position [m] to propagate to
	/// @param error_on_betax : getBetaX(...) returns its error in this variable
	/// not a const method because does a propagate to s!
	vector<H_BeamParticle>::iterator particle_i;
	float EX2=0,dummy, mean=getX(s,dummy), temp;

	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		particle_i->propagate(s);
		temp = particle_i->getX()-mean;
		EX2 += temp*temp;
	}
	EX2 /= (float)Nparticles;
	float emitx = getEmittanceX();	
	EX2 = (emitx==0)?0:(float) (EX2 /(float) (emitx*URAD))/URAD;
	error_on_betax = EX2 / (float) sqrt((double)2*Nparticles);
	return EX2;
}

float H_Beam::getBetaY(const float s, float& error_on_betay) {
        /// @param s is the position [m] to propagate to
        /// @param error_on_betay : getBetaY(...) returns its error in this variable
	/// not a const method because does a propagate to s!
 	vector<H_BeamParticle>::iterator particle_i;
	float EY2 =0, dummy, mean=getY(s,dummy), temp;

        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                particle_i->propagate(s);
		temp = particle_i->getY() - mean;
		EY2 += temp*temp;
	}
	EY2 /= (float)Nparticles;
	float emity = getEmittanceY();
	EY2 = (emity==0)?0:(float) (EY2 / (float) (emity*URAD))/URAD;
	error_on_betay = EY2 / (float) sqrt((double)2*Nparticles);
	return EY2;
}
/*
TGraphErrors * H_Beam::getBetaX(const float length, const unsigned int number_of_points) {
	/// @param length [m]
	/// @number_of_points in the graph (typ. 200)
	const unsigned int N = number_of_points;
	float * s = new float[N], * b = new float[N], * es = new float[N], * eb = new float[N];
	for (unsigned int i=0; i<N; i++) {
		s[i] = (float) fs_ini + i/(float)(N-1) *length; 
		b[i] = getBetaX(s[i],eb[i]);
		es[i] = 0;
	}
	TGraphErrors * betax = new TGraphErrors(N,s,b,es,eb);
	betax->SetLineColor(kBlack);
	betax->SetFillColor(kYellow);
	delete [] s;
	delete [] b;
	delete [] es;
	delete [] eb;
	return betax;
}

TGraphErrors * H_Beam::getBetaY(const float length, const unsigned int number_of_points) {
        /// @param length [m] 
        /// @number_of_points in the graph (typ. 200)
        const unsigned int N = number_of_points;
        float * s = new float[N], * b = new float[N], * es = new float[N], *eb = new float[N];
        for (unsigned int i=0; i<N; i++) {
                s[i] = (float) fs_ini + i/(float)(N-1) *length;
                b[i] = getBetaY(s[i],eb[i]);
		es[i]=0;
        }
        TGraphErrors * betay = new TGraphErrors(N,s,b,es,eb);
        betay->SetLineColor(kRed);
	betay->SetFillColor(kYellow);
        delete [] s;
        delete [] b;
        delete [] es;
        delete [] eb;
	return betay;
}
*/

float H_Beam::getX(const float s, float& error_on_posx) {
	vector<H_BeamParticle>::iterator particle_i;
	float mean=0;

        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                particle_i->propagate(s);
                mean += particle_i->getX();
        }
	mean = mean / (float) Nparticles;
	error_on_posx = mean / (float) sqrt((double)Nparticles);
        return mean;
}

float H_Beam::getY(const float s, float& error_on_posy) {
	vector<H_BeamParticle>::iterator particle_i;
	float mean=0;

        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                particle_i->propagate(s);
        	mean += particle_i->getY();
	}
	mean = mean / (float) Nparticles;
	error_on_posy = mean / (float) sqrt((double)Nparticles);
        return mean;
}

unsigned int H_Beam::getStoppedNumber(const H_AbstractBeamLine * beamline) {
	int number =0;
	vector<H_BeamParticle>::iterator particle_i;
	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		if(particle_i->stopped(beamline)) number++;
	}
	return number;
}

void H_Beam::getStoppingElements(const H_AbstractBeamLine * beamline, vector<H_OpticalElement>& list, vector<int>& numb) {
        vector<H_BeamParticle>::iterator particle_i;
        vector<H_OpticalElement>::iterator element_i;
        H_OpticalElement temp_el;
        vector<int>::iterator n_i;
        int number =0;
        bool found;

        list.clear();
        numb.clear();

        // creates a list of elements where beamParticles have stopped
        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                found = false;
                if(particle_i->stopped(beamline)) {
                        temp_el = *(particle_i->getStoppingElement());
                        if(list.size()==0) {
                                number=1;
                                list.push_back(temp_el);
                                numb.push_back(number);
                        } else {
                                for (element_i = list.begin(), n_i = numb.begin(); element_i < list.end(); element_i++, n_i++) {
                                        string el_i_name = element_i->getName();
                                        string temp_el_name = temp_el.getName();
                                        if(el_i_name == temp_el_name) {
                                                number = *n_i;
                                                number++;
                                                *n_i = number;
                                                found = true;
                                        }
                                }
                                if(!found) {
                                                number=1;
                                                list.push_back(temp_el);
                                                numb.push_back(number);
                                }
                        }
                } // if particle_i->stopped
        }// for particle_i
} // H_Beam::getStoppingElements

void H_Beam::printInitialState() const {
	cout << "Initial parameters of the beam" << endl;
	cout << "(x,y,s) = (" << fx_ini << "," << fy_ini << "," << fs_ini << ") ";
	cout << "(theta_x, theta_y) = (" << tx_ini << "," << ty_ini << ") ";
	cout << "energy = " << fe_ini << endl; 
	cout << endl;
	cout << "Dispersion on these values : " << endl;
	cout << "(dx,dy,ds) = (" << x_disp << "," << y_disp << "," << s_disp << ") ";
	cout << "(dtheta_x, dtheta_y) = (" << tx_disp << "," << ty_disp << ") ";
	cout << "de = " << e_disp << endl << endl;

	float mean_ini =0;
	vector<H_BeamParticle>::const_iterator particle_i;
	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		mean_ini += particle_i->getX();
	}
	mean_ini /= (float) beamParticles.size();
	cout << "Mean ini x = " << mean_ini << endl;
}

void H_Beam::printProperties() const {
	vector<H_BeamParticle>::const_iterator particle_i;
	cout << "There are " << Nparticles << " in the beam." << endl;
	for (particle_i = beamParticles.begin();particle_i < beamParticles.end(); particle_i++) {
		particle_i->printProperties();
	}
}

void H_Beam::printStoppingElements(const vector<H_OpticalElement>& list, const vector<int>& numb) const{
	/// see also H_Beam::getStoppingElements
	vector<H_OpticalElement>::const_iterator element_i;
	vector<int>::const_iterator n_i;

	// prints the list
	for (element_i=list.begin(), n_i = numb.begin(); element_i < list.end(); element_i++, n_i++) {
		cout << *n_i << " particules in " << element_i->getName();
		cout <<  " (" << element_i->getTypeString() << ") at " << element_i->getS() << "m" << endl;
		element_i->getAperture()->printProperties();
	}
} // H_Beam::printStoppingElements
/*
TH2F *  H_Beam::drawProfile(const float s) {
	/// not a const method because does a propagate to s!
	char title[50];
	sprintf(title,"Beam profile at %.2f m",s);
	vector<H_BeamParticle>::iterator particle_i;
	float xmax, xmin, ymax, ymin;
	float xx, yy, xborder, yborder;

	particle_i=beamParticles.begin();
	xmin = particle_i->getX();
	xmax = particle_i->getX();
	ymin = particle_i->getY();
        ymax = particle_i->getY();

	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                particle_i->propagate(s);
		xx = particle_i->getX();
		yy = particle_i->getY();

                xmax = xx>xmax ? xx : xmax;
		ymax = yy>ymax ? yy : ymax;
		xmin = xx<xmin ? xx : xmin;
		ymin = yy<ymin ? yy : ymin;
        }

	// in order to avoid some drawing problems, when the beam divergence is null
	if(xmax == xmin) xmax += 0.1;
	if(ymax == ymin) ymax += 0.1;
	
	xborder = (xmax-xmin)*0.2;
	yborder = (ymax-ymin)*0.2;

	xmax += xborder;
	xmin -= xborder;
	ymax += yborder;
	ymin -= yborder;

	TH2F * profile = new TH2F("profile",title,10000,xmin,xmax,1000,ymin,ymax);
        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		profile->Fill(particle_i->getX(), particle_i->getY());
	}
	return profile;
}*/
/*
TMultiGraph * H_Beam::drawBeamX(const int color) const {
	int mycolor = color;
	vector<H_BeamParticle>::const_iterator particle_i;
	TMultiGraph * beam_profile_x = new TMultiGraph("beam_profile_x","");
	
	for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
		TGraph * ppath_x = particle_i->getPath(0,mycolor);
		beam_profile_x->Add(ppath_x);
	}
	return beam_profile_x;
}

TMultiGraph * H_Beam::drawBeamY(const int color) const {
	int mycolor = color;
        vector<H_BeamParticle>::const_iterator particle_i;
        TMultiGraph * beam_profile_y = new TMultiGraph("beam_profile_y","");

        for (particle_i = beamParticles.begin(); particle_i < beamParticles.end(); particle_i++) {
                TGraph * ppath_y = particle_i->getPath(1,mycolor);
                beam_profile_y->Add(ppath_y);
        }
        return beam_profile_y;
}
*/
