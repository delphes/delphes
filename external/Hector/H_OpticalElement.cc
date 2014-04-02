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

/// \file H_OpticalElement.cc
/// \brief Class aiming at describing any beam optical element.

// c++ #includes
#include <iostream>
#include <string>

// ROOT #includes
#include "TPaveLabel.h"

// local #includes
#include "H_Parameters.h"
#include "H_TransportMatrices.h"
#include "H_Aperture.h"
#include "H_OpticalElement.h"
using namespace std;

/// called by the constructors
void H_OpticalElement::init(const string& nameE, const int typeE, const double s, const double k, const double l) {
	// this is called by the constructors
	// must be in public section !
	// xpos and ypos are vectors with n point. They define the aperture shape of the optical element.
	name = nameE;
	fs = s;
 	fk = k;
	xpos = 0;
	ypos = 0;
	txpos = 0;
	typos = 0;
 	element_length = l;
 	type = typeE;
	element_mat.ResizeTo(MDIM,MDIM);
 	element_mat = driftmat(l);

	// do NOT use setAperture for the initialisation ! there are protections there
	
	if(element_length<0)  { if(VERBOSE) cout<<"<H_OpticalElement> ERROR : Interpenetration of elements !"<<endl; }
	if(element_length==0) { if(VERBOSE) cout<<"<H_OpticalElement> WARNING : 0-length element ! (" << name << ") " << " at " << fs << endl; }
	betax =0;
	betay =0;
}

H_OpticalElement::H_OpticalElement() : element_aperture(new H_Aperture()) {
	init("",DRIFT,0.,0.,0.1);
}

H_OpticalElement::H_OpticalElement(const string& nameE, const int typeE, const double s, const double k, const double l) : element_aperture(new H_Aperture()) {
	init(nameE,typeE,s,k,l);
}

H_OpticalElement::H_OpticalElement(const string& nameE, const int typeE, const double s, const double k, const double l, H_Aperture* the_app) : element_aperture(the_app->clone()) {
	init(nameE,typeE,s,k,l);
}

H_OpticalElement::H_OpticalElement(const int typeE, const double s, const double k, const double l, H_Aperture* the_app) : element_aperture(the_app->clone()) {
	init("",typeE,s,k,l);
}

H_OpticalElement::H_OpticalElement(const int typeE, const double s, const double k, const double l) : element_aperture(new H_Aperture()) {
	init("",typeE,s,k,l);
}

H_OpticalElement::H_OpticalElement(const H_OpticalElement& el) {
	fs = el.fs;
	element_length = el.element_length;
	fk = el.fk;
	xpos = el.xpos;
	ypos = el.ypos;
	txpos = el.txpos;
	typos = el.typos;
	betax = el.betax;
	betay = el.betay;
	type = el.type;
	name = el.name;
	typestring = el.typestring;
	element_mat.ResizeTo(MDIM,MDIM);
	element_mat = el.element_mat;
	element_aperture = el.element_aperture->clone();
}

H_OpticalElement& H_OpticalElement::operator=(const H_OpticalElement& el) {
	if(this==&el) return *this;
        fs = el.fs;
        element_length = el.element_length;
        fk = el.fk;
        xpos = el.xpos;
        ypos = el.ypos;
	txpos = el.txpos;
	typos = el.typos;
        betax = el.betax;
        betay = el.betay;
        type = el.type;
        name = el.name;
        typestring = el.typestring;
	element_mat.ResizeTo(MDIM,MDIM);
        element_mat = el.element_mat;
        element_aperture = el.element_aperture->clone();
	return *this;
}

void H_OpticalElement::setAperture(const H_Aperture* ap) {
	// do NOT use setAperture in your constructor, as element_aperture is not initialized
	// this function do not take into account ap if ap=0
	// do nothing if element_mat = ap
		if (!ap) {cout << "<H_OpticalElement> Trying to set an empty pointer for the aperture ! Nothing done.\n"; return;}
		if (element_aperture != ap) {
			delete element_aperture;
			element_aperture = ap->clone();
		}
	return;
}
/*
void H_OpticalElement::setAperture(H_Aperture* ap) {
        // do NOT use setAperture in your constructor, as element_aperture is not initialized
        // this function do not take into account ap if ap=0
        // do nothing if element_mat = ap
                if (!ap) {cout << "<H_OpticalElement> Trying to set an empty pointer for the aperture ! Nothing done.\n"; return;}
                if (element_aperture != ap) {
                        delete element_aperture;
                        element_aperture = ap; //->clone();
                }
        return;
}
*/

std::ostream& operator<< (std::ostream& os, const H_OpticalElement& el) {
        os << el.typestring << el.name << "\t at s = " << el.fs << "\t length = "<< el.element_length;
        if(el.fk!=0) os <<"\t strength = " << el.fk;
        if(el.element_aperture->getType()!=NONE) {
                os << *(el.element_aperture) << endl;
        }
        os<<endl;
        if(el.element_length<0 && VERBOSE) 
		os <<"<H_OpticalElement> ERROR : Interpenetration of elements !"<<endl; 
        else if(el.element_length==0 && VERBOSE) 
		os <<"<H_OpticalElement> WARNING : 0-length "<< el.typestring << " !" << endl;
  return os;
}

void H_OpticalElement::showMatrix() const {
	printMatrix(element_mat);
	return;
}

void H_OpticalElement::drawAperture() const {
	element_aperture->draw(1);
	return;
}

TMatrix H_OpticalElement::getMatrix() {
	setMatrix(0,MP,1);
	return element_mat;
}

TMatrix H_OpticalElement::getMatrix(const float eloss, const float p_mass, const float p_charge) {
	setMatrix(eloss,p_mass,p_charge);
	return element_mat;
}

void H_OpticalElement::draw(const float meight, const float height) const{
	/// @param meight is the minimal extend of the graph
	/// @param height is the maximal extend of the graph
	float x1 = getS();
	float x2 = getS() + getLength();
	float y1 = meight;
	float y2 = height;
	TPaveLabel* cur_box = new TPaveLabel(x1,y1,x2,y2,"");
	cur_box->SetBorderSize(1);
	cur_box->SetFillStyle(1001);
	cur_box->SetFillColor((int)getType());
	cur_box->Draw();
}

TVectorD H_OpticalElement::getHitPosition(const TVectorD& init_pos, const double energy_loss, const double mp, const double qp) {
	if(!element_length) {
		// cout<<"O-length element ("<<getName()<<"), should not appear here !"<<endl;
		return init_pos;
	}
	// some declarations
	bool inside = false;
	double vec[MDIM] = {init_pos[INDEX_X]/URAD,  tan(init_pos[INDEX_TX]/URAD),
                            init_pos[INDEX_Y]/URAD,  tan(init_pos[INDEX_TY]/URAD), 
			    -energy_loss,            1};
	TMatrixD mat_init(1,MDIM,vec);
	TMatrixD mat_min(1,MDIM,vec);
	TMatrixD mat_max(1,MDIM,vec);
	TMatrixD mat_stop(1,MDIM,vec);
	H_OpticalElement* temp_el = clone();
	// initialasing boundaries
	double min_pos = 0;
	double max_pos = element_length/2.;
	double max_old = element_length/2.;
	// number of iterations 
	// (idea : fix precision instead of number of iterations + add security)
	// (idea : interpolate between max and min and give the error)
	const int N = 10;
	// starting search loop
	for(int i = 0; i < N; i++) {
		// fixing position to be investigated 
		temp_el->setLength(max_pos);
		// initialising the vector at the initial vector + possible shift/tilt
		mat_max[0][0] = mat_init[0][0] - temp_el->getX();
		mat_max[0][1] = mat_init[0][1] - tan(temp_el->getTX());
		mat_max[0][2] = mat_init[0][2] - temp_el->getY();
		mat_max[0][3] = mat_init[0][3] - tan(temp_el->getTY());
		// propagating
		mat_max *= temp_el->getMatrix(energy_loss,mp,qp);
		// compensating for the previously stated shifts/tilts
		mat_max[0][0] = mat_max[0][0] + temp_el->getX();
		mat_max[0][1] = mat_max[0][1] + tan(temp_el->getTX());
		mat_max[0][2] = mat_max[0][2] + temp_el->getY();
		mat_max[0][3] = mat_max[0][3] + tan(temp_el->getTY());
		// fixing new boundaries 
		if(temp_el->isInside(mat_max.GetMatrixArray()[0]*URAD,mat_max.GetMatrixArray()[2]*URAD)) {
			max_old = max_pos;
			max_pos = max_pos + (max_pos - min_pos)/2.;
			min_pos = max_old;
			inside = true;
		} else {
			max_pos = min_pos + (max_pos - min_pos)/2.;
			inside = false;
		}
		// end of loop
	}
	// if it passes at the last iteration, choosing the other range²
	if(inside) min_pos = max_old;
	// here the interpolation method : now we are sure that the intercept is between min_pos and max_pos
	// getting vector at min_pos (for first boundary) :
	bool precision_estimate = false;
	if(precision_estimate) {
		temp_el->setLength(min_pos);
		mat_min[0][0] = mat_init[0][0] - temp_el->getX();
		mat_min[0][1] = mat_init[0][1] - tan(temp_el->getTX());
		mat_min[0][2] = mat_init[0][2] - temp_el->getY();
		mat_min[0][3] = mat_init[0][3] - tan(temp_el->getTY());
		mat_min *= temp_el->getMatrix(energy_loss,mp,qp);
		mat_min[0][0] = mat_min[0][0] + temp_el->getX();
   		mat_min[0][1] = mat_min[0][1] + tan(temp_el->getTX());
 	  	mat_min[0][2] = mat_min[0][2] + temp_el->getY();
		mat_min[0][3] = mat_min[0][3] + tan(temp_el->getTY());
		mat_min[0][4] = min_pos + init_pos[4];
		// getting vector at max_pos (for second boundary) :
		temp_el->setLength(max_pos);
	   	mat_max[0][0] = mat_init[0][0] - temp_el->getX();
		mat_max[0][1] = mat_init[0][1] - tan(temp_el->getTX());
		mat_max[0][2] = mat_init[0][2] - temp_el->getY();
		mat_max[0][3] = mat_init[0][3] - tan(temp_el->getTY());
		mat_max *= temp_el->getMatrix(energy_loss,mp,qp);
		mat_max[0][0] = mat_max[0][0] + temp_el->getX();
		mat_max[0][1] = mat_max[0][1] + tan(temp_el->getTX());
		mat_max[0][2] = mat_max[0][2] + temp_el->getY();
		mat_max[0][3] = mat_max[0][3] + tan(temp_el->getTY());
		mat_max[0][4] = max_pos + init_pos[4];
	}
	// getting vector in the middle (for estimate) :
	temp_el->setLength((max_pos+min_pos)/2.);
	mat_stop[0][0] = mat_init[0][0] - temp_el->getX();
	mat_stop[0][1] = mat_init[0][1] - tan(temp_el->getTX());
	mat_stop[0][2] = mat_init[0][2] - temp_el->getY();
	mat_stop[0][3] = mat_init[0][3] - tan(temp_el->getTY());
	mat_stop *= temp_el->getMatrix(energy_loss,mp,qp);
	mat_stop[0][0] = mat_stop[0][0] + temp_el->getX();
	mat_stop[0][1] = mat_stop[0][1] + tan(temp_el->getTX());
	mat_stop[0][2] = mat_stop[0][2] + temp_el->getY();
	mat_stop[0][3] = mat_stop[0][3] + tan(temp_el->getTY());
	mat_stop[0][4] = (max_pos+min_pos)/2. + init_pos[4];

	double xys[LENGTH_VEC];
	xys[INDEX_X]=  mat_stop[0][0]*URAD;
	xys[INDEX_TX]= atan(mat_stop[0][1])*URAD;
	xys[INDEX_Y]=  mat_stop[0][2]*URAD;
	xys[INDEX_TY]= atan(mat_stop[0][3])*URAD;
	xys[INDEX_S]=  mat_stop[0][4] ;
	TVectorD temp_vec(LENGTH_VEC,xys);

	if(precision_estimate) {
		cout<<"--- Results and precision estimates ---"<<endl;
		cout<<"\t Stopping element : "<<getName()<<endl;
		cout<<"\t hit point s  : "<<mat_stop[0][4]<<" m +- "<<(mat_max[0][4]-mat_stop[0][4])*1000.<<" mm"<<endl;
		cout<<"\t hit point x  : "<<mat_stop[0][0]*URAD;
		cout<<" + "<<fabs(((mat_min[0][0]<mat_max[0][0])?(mat_max[0][0]-mat_stop[0][0])*URAD:(mat_min[0][0]-mat_stop[0][0])*URAD));
		cout<<" - "<<fabs(((mat_min[0][0]<mat_max[0][0])?(mat_stop[0][0]-mat_min[0][0])*URAD:(mat_stop[0][0]-mat_max[0][0])*URAD))<<" µm"<<endl;
		cout<<"\t hit point y  : "<<mat_stop[0][2]*URAD;
		cout<<" + "<<fabs(((mat_min[0][0]<mat_max[0][2])?(mat_max[0][2]-mat_stop[0][2])*URAD:(mat_min[0][2]-mat_stop[0][2])*URAD));
		cout<<" - "<<fabs(((mat_min[0][0]<mat_max[0][2])?(mat_stop[0][2]-mat_min[0][2])*URAD:(mat_stop[0][2]-mat_max[0][2])*URAD))<<" µm"<<endl;
		cout<<"\t hit point tx : "<<mat_stop[0][1]*URAD;
		cout<<" + "<<fabs(((mat_min[0][0]<mat_max[0][1])?(mat_max[0][1]-mat_stop[0][1])*URAD:(mat_min[0][1]-mat_stop[0][1])*URAD));
		cout<<" - "<<fabs(((mat_min[0][0]<mat_max[0][1])?(mat_stop[0][1]-mat_min[0][1])*URAD:(mat_stop[0][1]-mat_max[0][1])*URAD))<<" µrad"<<endl;
		cout<<"\t hit point ty : "<<mat_stop[0][3]*URAD;
		cout<<" + "<<fabs(((mat_min[0][0]<mat_max[0][3])?(mat_max[0][3]-mat_stop[0][3])*URAD:(mat_min[0][3]-mat_stop[0][3])*URAD));
		cout<<" - "<<fabs(((mat_min[0][0]<mat_max[0][3])?(mat_stop[0][3]-mat_min[0][3])*URAD:(mat_stop[0][3]-mat_max[0][3])*URAD))<<" µrad"<<endl;

	}

	delete temp_el;
	return temp_vec;
}
