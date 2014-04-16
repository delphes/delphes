/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_BeamLine.cc
/// \brief Reads external files and retrieves features of the real beam optical elements

// c++ #includes
#include <iostream>
#include <fstream>
#include <sstream>

// local #includes
#include "H_BeamLine.h"
#include "H_BeamLineParser.h"
#include "H_RectEllipticAperture.h"
#include "H_RectangularAperture.h"
#include "H_EllipticAperture.h"
#include "H_CircularAperture.h"
#include "H_RectangularDipole.h"
#include "H_SectorDipole.h"
#include "H_HorizontalQuadrupole.h"
#include "H_VerticalQuadrupole.h"
#include "H_HorizontalKicker.h"
#include "H_VerticalKicker.h"
#include "H_RectangularCollimator.h"
#include "H_Marker.h"

using namespace std;

// Caution : mad conventions for vertically(horizontally) focusing quadrupoles
// are opposite to Wille's. See fill() for more details.

H_BeamLine::H_BeamLine(const int si, const float length) : H_AbstractBeamLine(length){
	direction = (si >= abs(si)) ? 1 : -1;
	ips=0;
	ipx=0;
	ipy=0;
	iptx=0;
	ipty=0;
}

H_BeamLine::H_BeamLine(const H_BeamLine& beam) : H_AbstractBeamLine(beam) {
	direction = beam.direction;
	ips = beam.ips;
	ipx = beam.ipx;
	ipy = beam.ipy;
	iptx = beam.iptx;
	ipty = beam.ipty;
}

H_BeamLine& H_BeamLine::operator=(const H_BeamLine& beam) {
	if(this==&beam) return *this;
	direction = beam.direction;
	ips = beam.ips;
	ipx = beam.ipx;
    ipy = beam.ipy;
    iptx = beam.iptx;
    ipty = beam.ipty;
	return *this;
}

void H_BeamLine::findIP(const string filename) {
	findIP(filename,"IP5");
	return;
}

void H_BeamLine::findIP(const string filename, const string ipname) {
	// searches for the IP position in the extended table.
	ifstream tabfile(filename.c_str());
		if (! tabfile.is_open()) cout << "\t ERROR: I Can't open \"" << filename << "\"" << endl;
	bool found = false;
	int N_col=0;
	string headers[40];  // table of all column headers
	int header_type[40]; // table of all column types

	string temp_string;
	H_BeamLineParser e;
	istringstream curstring;

	while (	getline(tabfile,temp_string)) {
		curstring.clear(); // needed when using istringstream::str(string) several times !
		curstring.str(temp_string);

		// gets the features of each element
		if (found) {
			for (int col=0; col < N_col; col++) e.setProperties(curstring,header_type[col]);
			if(strstr(e.name.c_str(),ipname.c_str())) {
				ips = e.s; //e.printProperties();
				ipx = e.x;
				ipy = e.y;
				iptx = e.px*URAD;
				ipty = e.py*URAD;
			}
		} else if (strstr(temp_string.c_str(),"K0L")) {  //searches for the beginning of the element list.
			found = true;

			// reads the title of each column
			while(curstring.good()) { curstring >> headers[N_col]; if(headers[N_col]!="*") N_col++;}
			if(VERBOSE) cout << N_col << " columns identified" << endl;

			// identifies each column
			for (int col=0; col< N_col; col++)	header_type[col] =  column_identification(headers[col]);
			
			getline(tabfile,temp_string); // skip one line
		} // if(... "K0L")
	} // while (!eof)

	if (!found) cout << "\t ERROR ! IP not found." << endl;
	tabfile.close();
}

double* H_BeamLine::getIPProperties() {
	double* temp = new double[4];
	temp[0] = ipx;
	temp[1] = ipy;
	temp[2] = iptx;
	temp[3] = ipty;
	return temp;
}

void H_BeamLine::fill(const string filename) {
	fill(filename,1,"IP5");
	return;
}


void H_BeamLine::fill(const string filename, const int dir, const string ipname) {
	string headers[40];  // table of all column headers
	int header_type[40]; // table of all column types
	findIP(filename,ipname);
	ifstream tabfile(filename.c_str());
		if (! tabfile.is_open()) cout << "\t ERROR: I Can't open \"" << filename << "\"" << endl;
	if(VERBOSE) {
		cout<<"Using file : "<< filename <<" in the "<<((direction>0)?"positive":"negative")<<" direction."<<endl;
		cout<<"IP was found at position : "<<ips<<endl<<endl;
	}
	bool found = false; // found = true when the beginning of the element list is reached
	int type, N_col=0;
	float previous_betax =0, previous_betay=0; // needed when direction ==1
	float previous_dx =0, previous_dy=0; // needed when direction ==1
	float previous_x =0, previous_y=0; // needed when direction ==1

	string temp_string;
	H_BeamLineParser e;
	istringstream curstring;
	H_OpticalElement *el = 0;
	H_Aperture * ap = 0;


	while (getline(tabfile,temp_string)) {
		curstring.clear(); // needed when using several tims istringstream::str(string)
		curstring.str(temp_string);

        // gets the features of each element
		if (found) {
			// reads each column
			for (int col=0; col < N_col; col++) {
				e.setProperties(curstring,header_type[col]);
			}
			//e.printProperties();
			
			type =0; //init
			if(e.k0l!=0) { if(strstr(e.name.c_str(),"MB.")) type = SDIPOLE; else type = RDIPOLE;} //all SBEND seems to be called MB.xxxxx.xx
			else if(e.k1l!=0) type = (e.k1l>0)?HQUADRUPOLE:VQUADRUPOLE;
			else if(e.hkick!=0) type = HKICKER;
			else if(e.vkick!=0) type = VKICKER;
			else if(strstr(e.name.c_str(),"DRIFT")) type=DRIFT;
			else if((strstr(e.name.c_str(),"DFB"))||(strstr(e.name.c_str(),"TAS"))||(strstr(e.name.c_str(),"TAN"))||(strstr(e.name.c_str(),"TCL"))) type = RCOLLIMATOR;
			else if(strstr(e.name.c_str(),ipname.c_str())) { type=IP; }
			else type=MARKER;
			
			e.s = (e.s > ips) ? (e.s -ips - e.l)*(direction) : (e.s-ips)*(direction);
			
			el=0; //init
			switch(type) {
				case RDIPOLE:     { el = new H_RectangularDipole(e.name.c_str(),e.s,dir*e.k0l/e.l,e.l); } break;
				case SDIPOLE:     { el = new H_SectorDipole(e.name.c_str(),e.s,dir*e.k0l/e.l,e.l); } break;
				case VQUADRUPOLE: { el = new H_VerticalQuadrupole(e.name.c_str(),e.s,-(e.k1l/e.l),e.l); } break;
				case HQUADRUPOLE: { el = new H_HorizontalQuadrupole(e.name.c_str(),e.s,-(e.k1l/e.l),e.l); } break;
				case VKICKER:     { el = new H_VerticalKicker(e.name.c_str(),e.s,e.vkick,e.l); } break;
				case HKICKER:     { el = new H_HorizontalKicker(e.name.c_str(),e.s,e.hkick,e.l); } break;
				case RCOLLIMATOR: { el = new H_RectangularCollimator(e.name.c_str(),e.s,e.l); } break;
				case DRIFT: {previous_betax = e.betx; previous_betay = e.bety; previous_dx = e.dx; previous_dy = e.dy; previous_x = e.x; previous_y = e.y;} break;
				case IP: { el = new H_Marker(e.name.c_str(),e.s); } break;
				default: break;
			}

			ap = 0; //init
			if(e.aper_1 !=0 || e.aper_2 !=0 || e.aper_3 !=0 || e.aper_4 != 0) {
				e.aper_1 *= URAD;
				e.aper_2 *= URAD;
				e.aper_3 *= URAD;
				e.aper_4 *= URAD; // in [m] in the tables !

				if(strstr(e.apertype.c_str(),"RECTELLIPSE")) ap = new H_RectEllipticAperture(e.aper_1,e.aper_2,e.aper_3,e.aper_4,0,0);
				else if(strstr(e.apertype.c_str(),"CIRCLE")) ap = new H_CircularAperture(e.aper_1,0,0);
				else if(strstr(e.apertype.c_str(),"RECTANGLE")) ap = new H_RectangularAperture(e.aper_1,e.aper_2,0,0);
				else if(strstr(e.apertype.c_str(),"ELLIPSE")) ap = new H_EllipticAperture(e.aper_1,e.aper_2,0,0);
			}

			if(el!=0) {
				if (direction<0) {
					el->setBetaX(e.betx);
					el->setBetaY(e.bety);
					el->setDX(e.dx);
					el->setDY(e.dy);
					el->setRelX(e.x);
					el->setRelY(e.y);
				} else {
					el->setBetaX(previous_betax);
					el->setBetaY(previous_betay);
					el->setDX(previous_dx);
					el->setDY(previous_dy);
					el->setRelX(previous_x);
					el->setRelY(previous_y);
				}
				if(ap!=0) { 
					el->setAperture(ap); 
				//		delete ap; // ap deleted in H_AbstractBeamLine::~H_AbstractBeamLine
				}
	
				/// Parses all the elements, but only keeps the ones from the IP till the desired length
				if(e.s>=0 && e.s<beam_length) add(*el);

			// delete el; // el deleted in H_AbstractBeamLine::~H_AbstractBeamLine
			} // if(el!=0)
		} // if (found)
		else if(strstr(temp_string.c_str(),"K0L")) { // if (!found)
		//searches for the beginning of the element list.
			found = true;

			// reads the title of each column
			while(curstring.good()) { curstring >> headers[N_col]; if(headers[N_col]!="*") N_col++;}
			if(VERBOSE) cout << N_col << " columns identified" << endl;

			// identifies each column
			for (int col=0; col< N_col; col++) {
				header_type[col] =  column_identification(headers[col]);
			}
			getline(tabfile,temp_string);
		} // if temp_string <-> "K0L"

	} // while (!eof)
	tabfile.close();
}

