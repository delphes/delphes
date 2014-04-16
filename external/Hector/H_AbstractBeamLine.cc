/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/


/// \file H_AbstractBeamLine.cc
/// \brief Class describing ideal beamline.
///
/// Units : angles [rad], distances [m], energies [GeV], c=[1].

// c++ #includes
#include <iostream>
#include <cmath>

// ROOT #includes
//#include "TPaveLabel.h"
//#include "TLine.h"
//#include "TGaxis.h"
//#include "TLegend.h"
//#include "TF1.h"
//#include "TROOT.h"

// local #includes
#include "H_Parameters.h"
#include "H_TransportMatrices.h"
#include "H_Drift.h"
#include "H_AbstractBeamLine.h"
#include "H_BeamParticle.h"
#include "H_RomanPot.h"
using namespace std;

void H_AbstractBeamLine::init(const float length) {
	beam_mat = new TMatrix(MDIM,MDIM);
  	beam_length = length;
	H_Drift * drift0 = new H_Drift("Drift0",0.,length);
	add(drift0);
	return;
}

H_AbstractBeamLine::H_AbstractBeamLine(const H_AbstractBeamLine& beamline) {
	elements = beamline.elements;
	matrices = beamline.matrices;
	beam_mat = new TMatrix(*(beamline.beam_mat));
	beam_length = beamline.beam_length;
}

H_AbstractBeamLine& H_AbstractBeamLine::operator=(const H_AbstractBeamLine& beamline) {
	if(this== &beamline) return *this;
        elements = beamline.elements;
        matrices = beamline.matrices;
        beam_mat = new TMatrix(*(beamline.beam_mat));
	beam_length = beamline.beam_length;
	return *this;
}

H_AbstractBeamLine::~H_AbstractBeamLine() {
	vector<H_OpticalElement*>::iterator element_i;
	for (element_i = elements.begin(); element_i<elements.end(); element_i++) {
		delete (*element_i);
	}
	elements.clear(); 
	matrices.clear(); 
	delete beam_mat;
}

void H_AbstractBeamLine::add(H_OpticalElement * newElement) {
	/// @param newElement is added to the beamline
//	H_OpticalElement * el = new H_OpticalElement(*newElement);
//	H_OpticalElement * el = const_cast<H_OpticalElement*> newElement; 
//	H_OpticalElement * el = newElement;
	elements.push_back(newElement);
 	float a = newElement->getS()+newElement->getLength();
	if (a > beam_length)	{
		beam_length = a;
		if(VERBOSE) cout<<"WARNING : element ("<< newElement->getName()<<") too far away. The beam length has been extended to "<< beam_length << ". "<<endl;
	}
	calcSequence();
	calcMatrix();
}

void H_AbstractBeamLine::add(H_OpticalElement & newElement) {
	/// @param newElement is added to the beamline
//	H_OpticalElement * el = new H_OpticalElement(newElement);
//	elements.push_back(el);
	elements.push_back(&newElement);
 	float a = newElement.getS()+newElement.getLength();
	if (a > beam_length)	{
		beam_length = a;
		if(VERBOSE) cout<<"WARNING : element ("<< newElement.getName()<<") too far away. The beam length has been extended to "<< beam_length << ". "<<endl;
	}
	calcSequence();
	calcMatrix();
}

const TMatrix * H_AbstractBeamLine::getBeamMatrix() const {
	TMatrix * mat = new TMatrix(*beam_mat);
	return mat;
}

const TMatrix * H_AbstractBeamLine::getBeamMatrix(const float eloss,const float p_mass, const float p_charge) {

	vector<H_OpticalElement*>::iterator element_i;
    TMatrix calc_mat(MDIM,MDIM);

    // initialization
	calc_mat.UnitMatrix();
	
	// multiplies the matrix of each beam's element
	// and add each product matrix to the list of matrices.
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix(eloss,p_mass,p_charge);
	}
	const TMatrix* bmat = new TMatrix(calc_mat);
	return bmat;
}

const TMatrix * H_AbstractBeamLine::getPartialMatrix(const string elname, const float eloss, const float p_mass, const float p_charge) {

	vector<H_OpticalElement*>::iterator element_i;
	TMatrix calc_mat(MDIM,MDIM);

	calc_mat.UnitMatrix();

	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix(eloss,p_mass,p_charge);
		if(elname==(*element_i)->getName()) {
			const TMatrix* bmat = new TMatrix(calc_mat);
			return bmat;
		}
	}
	cout<<"Element "<<elname<<" desn't exist. Returning full beam matrix"<<endl;
	const TMatrix* bmat = new TMatrix(calc_mat);
	return bmat;
}

const TMatrix * H_AbstractBeamLine::getPartialMatrix(const unsigned int element_position) const {
 	//const int N = (element_position<0)?0:(( (element_position)>elements.size()-1)?elements.size()-1:element_position);
	const int N = (element_position>elements.size()-1)?elements.size()-1:element_position;
	return &(*(matrices.begin()+N)); // //for optimization of the code :same as return &matrices[N];
}

const TMatrix * H_AbstractBeamLine::getPartialMatrix(const H_OpticalElement * element) const{
	// returns the transport matrix to transport until the end of the specified element
	// !!! 2 elements should never have the same name in "elements" !!!

	vector<H_OpticalElement*>::const_iterator element_i;
	vector<TMatrix>::const_iterator matrix_i;
	TMatrix * calc_mat = new TMatrix(MDIM,MDIM);

	// parses the list of optical elements and find the searched one
	for(element_i = elements.begin(),matrix_i = matrices.begin(); element_i < elements.end(); element_i++, matrix_i++) {
		if(element->getName() == (*element_i)->getName()) {
			// element has been found
			calc_mat = const_cast<TMatrix*>( &(*matrix_i));
        	}
	}
	return (const TMatrix*) calc_mat;
}

H_OpticalElement * H_AbstractBeamLine::getElement(const unsigned int element_position) {
	const unsigned int N = (element_position>elements.size())?elements.size():element_position;
	return *(elements.begin()+N);//for optimization of the code :same as return &elements[N];
}

const H_OpticalElement * H_AbstractBeamLine::getElement(const unsigned int element_position) const {
	const unsigned int N = (element_position>elements.size())?elements.size():element_position;
        return *(elements.begin()+N);//for optimization of the code :same as return &elements[N];
}


H_OpticalElement * H_AbstractBeamLine::getElement(const string el_name) {
	for(unsigned int i=0; i < elements.size(); i++) {
		if( (*(elements.begin()+i))->getName() == el_name ) 
			return *(elements.begin()+i);
	} // if found -> return ; else : not found at all !	
	cout<<"Element "<<el_name<<" not found"<<endl;
	return *(elements.begin()+1);
}

const H_OpticalElement * H_AbstractBeamLine::getElement(const string el_name) const {
        for(unsigned int i=0; i < elements.size(); i++) {
		if( (*(elements.begin()+i))->getName() == el_name)
			return *(elements.begin()+i);
	} // if found -> return ; else : not found at all !
	cout<<"Element "<<el_name<<" not found"<<endl;
	return *(elements.begin()+1);
}

void H_AbstractBeamLine::printProperties() const {
	vector<H_OpticalElement*>::const_iterator element_i;
	cout << "Pointeurs des elements du faisceau" << endl;
	for (element_i = elements.begin(); element_i < elements.end(); element_i++) {
		cout << (int)(element_i-elements.begin()) << "\t" << (*element_i)->getName() << "\t" << (*element_i)->getS() << endl;	
	}
	return;
}

void H_AbstractBeamLine::showElements() const{
	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		(*element_i)->printProperties();
	}
	cout << "Beam length = " << beam_length << endl;
	cout << "Number of elements (including drifts) = " << getNumberOfElements() << endl;
	return;
}

void H_AbstractBeamLine::showElements(const int type_el) const{
	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		if ((*element_i)->getType()==type_el)
			(*element_i)->printProperties();
	}
	return;
}

void H_AbstractBeamLine::showMatrix() const {
	cout << "Transport matrix for the whole beam : " << endl;
	cout << "(x,x',...) = (x*,x'*,...) M " <<endl;
  	printMatrix(beam_mat);
	return;
}

void H_AbstractBeamLine::showMatrices() const{
	// prints the list of all transport matrices, from the whole beam.

	vector<TMatrix>::const_iterator matrix_i;
	vector<H_OpticalElement*>::const_iterator element_i;
	TMatrix temp(MDIM,MDIM);

	for(matrix_i = matrices.begin(), element_i = elements.begin(); matrix_i < matrices.end(); matrix_i++, element_i++) {
		temp = *matrix_i;
		cout << "Matrix for transport until s=" << (*element_i)->getS() + (*element_i)->getLength() << "m (" << (*element_i)->getName() << "). " << endl;
		printMatrix(&temp);
		cout << endl;
	}
	return ;
}

void H_AbstractBeamLine::calcSequence() {
		// reorders the elements, computes the drifts;

	vector<H_OpticalElement*> temp_elements;
	vector<H_OpticalElement*>::iterator element_i;
		// element_i is a pointer to elements[i]

	if(elements.size()==1) { return; }

	// getting rid of drifts before calculating
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		if((*element_i)->getType() == DRIFT) {elements.erase(element_i); }
	}

	// ordering the elements in position
	sort(elements.begin(),elements.end(),ordering());
	// inserting the drifts before the other elements
	float current_pos = 0;
	float drift_length=0;

	for(element_i=elements.begin(); element_i < elements.end(); element_i++) {
		drift_length = (*element_i)->getS() - current_pos;
		if(drift_length>0) {
  			H_Drift *dr = new H_Drift(current_pos,drift_length);
			temp_elements.push_back(dr); 
		}
		temp_elements.push_back(*element_i);
		current_pos = (*element_i)->getS() + (*element_i)->getLength();
	}
	
	//adding the last drift
	drift_length = beam_length - current_pos;
	if (drift_length>0) {
			H_Drift *dr = new H_Drift(current_pos,drift_length);
			temp_elements.push_back(dr);
	}
	elements.clear();
	for(element_i=temp_elements.begin(); element_i < temp_elements.end(); element_i++) {
		elements.push_back(*element_i);
	}
}

void H_AbstractBeamLine::calcMatrix() {
	// computes the transport matrix for the beam upto here...
	vector<H_OpticalElement*>::iterator element_i;
	TMatrix calc_mat(MDIM,MDIM);

	// initialization
	matrices.clear();
	calc_mat.UnitMatrix();

	// multiplies the matrix of each beam's element
	// and add each product matrix to the list of matrices.
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		calc_mat *= (*element_i)->getMatrix();
		matrices.push_back(calc_mat);
  	}

	*beam_mat = calc_mat;
	return;
}

float qh(float k) {
        float beta = (log((float)10.0))/0.05;
		// put (std::log((float)10.0)) instead of log(10) to avoid compilation errors
        return 0.8*(1-exp(-beta*fabs(k)));
}

float dh(float k) {
        float psi = (log((float)10.0))/0.002;
		// put (std::log((float)10.0)) instead of log(10) to avoid compilation errors
        return 0.8*(1-exp(-psi*fabs(k)));
}

void H_AbstractBeamLine::draw() const{
/*	gROOT->SetStyle("Plain");
	TLegend* leg = new TLegend(0.85,0.50,1,1,"Legend");
	leg->SetBorderSize(1);
	TBox* b1 = new TBox();
	TBox* b2 = new TBox(0,0,10,10);
	TBox* b3 = new TBox(0,0,0,0);
	TBox* b4 = new TBox(0,0,0,0);
	TBox* b5 = new TBox(0,0,0,0);
	TBox* b6 = new TBox(0,0,0,0);
	TBox* b7 = new TBox(0,0,0,0);
	b1->SetFillColor(RDIPOLE);
	b2->SetFillColor(SDIPOLE);
	b3->SetFillColor(VQUADRUPOLE);
	b4->SetFillColor(HQUADRUPOLE);
	b5->SetFillColor(HKICKER);
	b6->SetFillColor(VKICKER);
	b7->SetFillColor(RCOLLIMATOR);
	leg->AddEntry(b1,"R-Dipole");
	leg->AddEntry(b2,"S-Dipole");
	leg->AddEntry(b3,"V-Quadrupole");
	leg->AddEntry(b4,"H-Quadrupole");
	leg->AddEntry(b5,HKICKERNAME);
	leg->AddEntry(b6,VKICKERNAME);
	leg->AddEntry(b7,"RCollimator");
	leg->Draw();
*/ 
/*	TLine* l1 = new TLine(0.05,0.5,0.95,0.5);
	TLine* l2 = new TLine(0.1,0.1,0.1,0.9);
	TLine* l3 = new TLine(0.9,0.1,0.9,0.9);
	TPaveLabel* p1 = new TPaveLabel(0.05,0.5,0.1,0.6,"IP");
	TPaveLabel* p2 = new TPaveLabel(0.9,0.5,0.95,0.6,"RP");
	TGaxis* a1 = new TGaxis(0.1,0.1,0.9,0.1,0,beam_length);
	a1->SetLabelSize(0.08);
	p1->SetBorderSize(1);
	p2->SetBorderSize(1);
	p1->SetFillColor(0);
	p2->SetFillColor(0);
	l1->Draw();
	l2->Draw();
	l3->Draw();
	p1->Draw();
	p2->Draw();
	float x1,x2,y1,y2;
	vector<TPaveLabel*> boxes;
	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
		x1 = 0.1 + ((*element_i)->getS()/beam_length)*0.8;
		x2 = x1 + ((*element_i)->getLength()/beam_length)*0.8;
		if((*element_i)->getType()>5) {
			y1 = 0.3;
			y2 = 0.7;
		}
		else if((*element_i)->getType()>3) {
			y1 = 0.5 - qh((*element_i)->getK()*(*element_i)->getLength())/2.;
			y2 = 0.5 + qh((*element_i)->getK()*(*element_i)->getLength())/2.;
		} else {
			y1 = 0.5 - dh((*element_i)->getK()*(*element_i)->getLength())/2.;
			y2 = 0.5 + dh((*element_i)->getK()*(*element_i)->getLength())/2.;
		}
		TPaveLabel* cur_box = new TPaveLabel(x1,y1,x2,y2,"");
		cur_box->SetFillStyle(1);
		cur_box->SetFillColor(((int)(*element_i)->getType()));
		cur_box->SetBorderSize(1);
		if((*element_i)->getType()!=DRIFT) boxes.push_back(cur_box);
	}
	vector<TPaveLabel*>::iterator box_i;
	for(box_i = boxes.begin(); box_i < boxes.end(); box_i++) {
		(*box_i)->Draw();
	}
	a1->Draw();
*/	
	return;
}

void H_AbstractBeamLine::drawX(const float a_min, const float a_max) const{
	/// @param a_min defines the size of the drawing
	/// @param a_max defines the size of the drawing
	const int N = getNumberOfElements();
	for(int i=0;i<N;i++) {
		float height = fabs(a_max);
        	float meight = fabs(a_min);
        	float size = (height>meight)?meight:height;
		float middle = getElement(i)->getX()*URAD;
        	if(getElement(i)->getType()!=DRIFT) getElement(i)->draw(middle+size/2.,middle-size/2.);
    	}
}

void H_AbstractBeamLine::drawY(const float a_min, const float a_max) const{
        /// @param a_min defines the size of the drawing 
        /// @param a_max defines the size of the drawing
        const int N = getNumberOfElements();
        for(int i=0;i<N;i++) {
                float height = fabs(a_max);
                float meight = fabs(a_min);
                float size = (height>meight)?meight:height;
                float middle = getElement(i)->getY()*URAD;
                if(getElement(i)->getType()!=DRIFT) getElement(i)->draw(middle+size/2.,middle-size/2.);
        }
}

void H_AbstractBeamLine::moveElement(const string name, const float new_s) {
	/// @param name identifies the element to move
	/// @param new_s is where to put it
       vector<H_OpticalElement*>::iterator element_i;
       for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
	       if(name==(*element_i)->getName()) { (*element_i)->setS(new_s); }
       }

       calcSequence();
       calcMatrix();
       return;
}

void H_AbstractBeamLine::alignElement(const string name, const float disp_x, const float disp_y) {
	/// @param name identifies the element to move
	/// @param disp_x identifies the displacement to add in x [\f$ \mu m \f$]
	/// @param disp_y identifies the displacement to add in y [\f$ \mu m \f$]
		vector<H_OpticalElement*>::iterator element_i;
		for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if(name==(*element_i)->getName()) { 
				(*element_i)->setX((*element_i)->getX()+disp_x);
				(*element_i)->setY((*element_i)->getY()+disp_y);
				return ;
			}
		}
		cout<<"Element "<<name<<" not found."<<endl; 
		if(VERBOSE) cout<<"Element "<<name<<" not found."<<endl;
		return;
}

void H_AbstractBeamLine::tiltElement(const string name, const float ang_x, const float ang_y) {
	/// @param name identifies the element to move
	/// @param ang_x identifies the angle to add in x
	/// @param ang_y identifies the angle to add in y
		vector<H_OpticalElement*>::iterator element_i;
		for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if(name==(*element_i)->getName()) {
				(*element_i)->setTX((*element_i)->getTX()+ang_x);
				(*element_i)->setTY((*element_i)->getTY()+ang_y);
				return ;
			}
		}
		if(VERBOSE) cout<<"Element "<<name<<" not found."<<endl;
		return;
}

void H_AbstractBeamLine::offsetElements(const float start, const float offset) {
	/// @param start After this s [m] coordinate, all elements will be offset.
	/// @param offset In meters

	extern int relative_energy;
	if(!relative_energy) {
		vector<H_OpticalElement*>::iterator element_i;
   	     for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
			if((*element_i)->getS() > start ) {
    	       (*element_i)->setX(offset);
       	    }
		}
	}
}

/*
TGraph * H_AbstractBeamLine::getBetaX() const{
        const int N = elements.size();
        float * s = new float[N], * b = new float[N], temp;
	int i=0, n=N;

	vector<H_OpticalElement*>::const_iterator element_i;
	for(element_i = elements.begin(); element_i < elements.end(); element_i++) {	
		temp=(*element_i)->getBetaX();	
		if (temp !=0) {
			b[i] = (*element_i)->getBetaX();
			s[i] = (*element_i)->getS();
			i++;
			n=i;
		}
	}

        TGraph * betax = new TGraph(n,s,b);
        betax->SetLineColor(1);
	betax->SetLineStyle(2);
	delete [] s;
	delete [] b;
        return betax;
}

TGraph * H_AbstractBeamLine::getBetaY() const{
        const int N = elements.size();
        float * s = new float[N], * b = new float[N], temp;
        int i=0, n=N;

        vector<H_OpticalElement*>::const_iterator element_i;
        for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
                temp=(*element_i)->getBetaY();
                if (temp !=0) {
                        b[i] = (*element_i)->getBetaY();
                        s[i] = (*element_i)->getS();
                        i++;
                        n=i;
                }
        }

        TGraph * betay = new TGraph(n,s,b);
        betay->SetLineColor(2);
        betay->SetLineStyle(2);
	delete [] s;
	delete [] b;
        return betay;
}

TGraph * H_AbstractBeamLine::getDX() const{
        const int N = elements.size();
        float * s = new float[N], * d = new float[N], temp;
        int i=0, n=N;

        vector<H_OpticalElement*>::const_iterator element_i;
        for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
                temp=(*element_i)->getDX();
                if (temp !=0) {
                        d[i] = (*element_i)->getDX();
                        s[i] = (*element_i)->getS();
                        i++;
                        n=i;
                }
        }

        TGraph * dispx = new TGraph(n,s,d);
        dispx->SetLineColor(8);
        dispx->SetLineStyle(2);
        delete [] s;
        delete [] d;
        return dispx;
}

TGraph * H_AbstractBeamLine::getDY() const{
        const int N = elements.size();
        float * s = new float[N], * d = new float[N], temp;
        int i=0, n=N;
 
        vector<H_OpticalElement*>::const_iterator element_i;
        for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
                temp=(*element_i)->getDY();
                if (temp !=0) {
                        d[i] = (*element_i)->getDY();
                        s[i] = (*element_i)->getS();
                        i++;
                        n=i;
                }
        }
 
        TGraph * dispy = new TGraph(n,s,d);
        dispy->SetLineColor(kBlue);
        dispy->SetLineStyle(2);
        delete [] s;
        delete [] d;
        return dispy;
}


TGraph * H_AbstractBeamLine::getRelX() const{
        const int N = elements.size();
        float * s = new float[N], * r = new float[N], temp;
        int i=0, n=N;

        vector<H_OpticalElement*>::const_iterator element_i;
        for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
                temp=(*element_i)->getRelX();
		if((*element_i)->getType() != DRIFT) {
			r[i] = (*element_i)->getRelX();
			s[i] = (*element_i)->getS();
			i++;
			n=i;
		}
	}

        TGraph * relx = new TGraph(n,s,r);
        relx->SetLineColor(kBlack);
	relx->SetMarkerStyle(kOpenSquare);
	relx->SetMarkerSize(0.6);
        relx->SetLineStyle(2);
        delete [] s;
        delete [] r;
        return relx;
}

TGraph * H_AbstractBeamLine::getRelY() const{
        const int N = elements.size();
        float * s = new float[N], * r = new float[N], temp;
        int i=0, n=N;

        vector<H_OpticalElement*>::const_iterator element_i;
        for(element_i = elements.begin(); element_i < elements.end(); element_i++) {
                temp=(*element_i)->getRelY();
                if((*element_i)->getType() != DRIFT) {
                        r[i] = (*element_i)->getRelY();
                        s[i] = (*element_i)->getS();
                        i++;
                        n=i;
                }
        }

        TGraph * rely = new TGraph(n,s,r);
        rely->SetLineColor(kRed);
	rely->SetMarkerStyle(kOpenSquare);
	rely->SetMarkerSize(0.6);
        rely->SetLineStyle(2);
        delete [] s;
        delete [] r;
        return rely;
}
*/
