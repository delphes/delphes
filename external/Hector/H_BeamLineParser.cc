/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_BeamLineParser.cc
/// \brief Reader for madx tables
///
/// Notes : V�ifier que tous les SBEND sont toujours appel� MB. Et seulement comme � !
/// V�ifier qu'il n'y a pas de probl�e d'inversion H/V QUADRUPOLES
/// no distinction between H and Vkickers ?
/// The identification of the element is based on the values of k1l, k2l, hkick, vkick and on their name.
/// x and y are put in m

// c++ #includes
#include <iostream>
#include <fstream>
#include <sstream>

//local #includes
#include "H_BeamLineParser.h"
#include "H_Parameters.h"
using namespace std;


/// Identifies the column content from its header
int column_identification(const string header) {
	// identifies the column type from its name

	     if (header=="NAME")    {return MADX_NAME; }
	else if (header=="KEYWORD") {return MADX_KEYWORD; }
	else if (header=="S")       {return MADX_S; }
	else if (header=="L")       {return MADX_L; }
	else if (header=="K0L")     {return MADX_K0L; }
	else if (header=="K1L")     {return MADX_K1L;}
	else if (header=="K2L")     {return MADX_K2L;}
	else if (header=="K3L")     {return MADX_K3L;}
	else if (header=="HKICK")   {return MADX_HKICK;}
	else if (header=="VKICK")   {return MADX_VKICK;}
	else if (header=="KICK")    {return MADX_KICK;}
	else if (header=="BETX")    {return MADX_BETX;}
	else if (header=="BETY")    {return MADX_BETY;}
	else if (header=="ALFX")    {return MADX_ALFX;}
	else if (header=="ALFY")    {return MADX_ALFY;}
	else if (header=="MUX")     {return MADX_MUX;}
	else if (header=="MUY")     {return MADX_MUY;}
	else if (header=="DX")      {return MADX_DX;}
	else if (header=="DY")      {return MADX_DY;}
	else if (header=="DPX")     {return MADX_DPX;}
	else if (header=="DPY")     {return MADX_DPY;}
	else if (header=="X")       {return MADX_X;}
	else if (header=="Y")       {return MADX_Y;}
	else if (header=="PX")      {return MADX_PX;}
	else if (header=="PY")      {return MADX_PY;}
	else if (header=="APERTYPE") {return MADX_APERTYPE;}
	else if (header=="APER_1")  {return MADX_APER_1;}
	else if (header=="APER_2")  {return MADX_APER_2;}
	else if (header=="APER_3")  {return MADX_APER_3;}
	else if (header=="APER_4")  {return MADX_APER_4;}
	else if (header=="PARENT")  {return MADX_PARENT;}
	return MADX_UNKNOWN;
}

void H_BeamLineParser::init() {
	name=""; 
	apertype="";
	keyword="";
	parent="";
	s=0;  l=0;  k0l=0;  k1l=0;  k2l=0;  k3l=0;  hkick=0;  vkick=0;  betx=0;
	alfx=0;  mux=0;  dx=0;  dpx=0;  x=0;  px=0;  bety=0;  alfy=0;  muy=0;
	dy=0;  dpy=0;  y=0;  py=0;  aper_1=0;  aper_2=0;  aper_3=0;  aper_4=0;
}

void H_BeamLineParser::setProperties(istream& input, const unsigned int col_type) {

	switch(col_type) {
		case   MADX_NAME:    input >> name; break;
		case   MADX_KEYWORD: input >> keyword; break;
		case   MADX_S:       input >> s; break;
		case   MADX_L:       input >> l;   break;
		case   MADX_K0L:     input >> k0l; break;
		case   MADX_K1L :    input >> k1l; break;
		case   MADX_K2L:     input >> k2l; break;
		case   MADX_K3L:     input >> k3l; break;
		case   MADX_HKICK:   input >> hkick; break;
		case   MADX_VKICK:   input >> vkick; break;
		case   MADX_BETX:    input >> betx; break;
		case   MADX_ALFX:    input >> alfx; break;
		case   MADX_MUX:     input >> mux; break;
		case   MADX_DX:      input >> dx; break;
		case   MADX_DPX:     input >> dpx; break;
		case   MADX_X:       input >> x; x *=URAD; break;
		case   MADX_PX:      input >> px; break;
		case   MADX_BETY:    input >> bety; break;
		case   MADX_ALFY:    input >> alfy; break;
		case   MADX_MUY:     input >> muy; break;
		case   MADX_DY:      input >> dy; break;
		case   MADX_DPY:     input >> dpy; break;
		case   MADX_Y:       input >> y; y*=URAD; break;
		case   MADX_PY:      input >> py; break;
		case   MADX_APERTYPE: input >> apertype; break;
		case   MADX_APER_1:  input >> aper_1; break;
		case   MADX_APER_2:  input >> aper_2; break;
		case   MADX_APER_3:  input >> aper_3; break;
		case   MADX_APER_4:  input >> aper_4; break;
 		case   MADX_PARENT:  input >> parent; break;
		default:break;
	} // switch
}


void H_BeamLineParser::printProperties() const {
//KEYWORD	NAME	PARENT	L	K0L	K1L	K2L	K3L	S	BETX	BETY	DX	DY	XC	YC	ALFX	ALFY	MUX	MUY	DPX	PXC	PYC
	cout << " keyword = " << keyword;
	cout << " name = " << name;
	cout << " l = " << l;
	cout << " k0l = " << k0l;
	cout << " k1l = " << k1l;
	cout << " k2l = " << k2l;
	cout << " k3l = " << k3l;
	cout << " s = " << s << endl;
	cout << endl;
}
