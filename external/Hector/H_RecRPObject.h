#ifndef _H_RecRPObject_
#define _H_RecRPObject_

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

/// \file H_RecRPObject.h
/// \brief Reconstructed information from objects detected by two roman pots

#include "TF1.h"
#include <cmath>

#include "H_BeamLine.h"

#define NOT_YET_COMPUTED -666

/// Reconstruction of parameters at IP using measurements with roman pots
class H_RecRPObject {
	public:
		H_RecRPObject();
		H_RecRPObject(const float, const float, const H_AbstractBeamLine* );
		H_RecRPObject(const float, const float, const H_AbstractBeamLine& ); // old-fashioned -- to be deleted
		H_RecRPObject(const H_RecRPObject&);
		H_RecRPObject& operator=(const H_RecRPObject&);
		~H_RecRPObject() {delete f_1; delete f_2; delete g_1; delete g_2; 
				  delete d_1; delete d_2; delete k_1; delete k_2; 
				  delete l_1; delete l_2; 
				  delete thebeam;};

		inline float getX1() const {return x1;};
		inline float getX2() const {return x2;};
		inline float getY1() const {return y1;};
		inline float getY2() const {return y2;};
		inline float getS1() const {return s1;};
		inline float getS2() const {return s2;};

		float getTX();
		float getTY();
		float getE();
		float getE(int );
		float getQ2();
		float getPt();

		void setPositions(const float, const float, const float, const float);
		void setPosition_det1(const float, const float);
		void setPosition_det2(const float, const float);
		void setDetPos(const float, const float);
		void setERange(const float, const float);
		void computeERange();
		void initialize();
		void computeAll();

	protected:

		float emin, emax;
		float x1, x2, y1, y2, s1, s2;
		float txip, tyip, energy, q2, pt;
		H_AbstractBeamLine* thebeam;
		TF1* f_1, *f_2;
		TF1* g_1, *g_2;
		TF1* d_1, *d_2;
		TF1* k_1, *k_2;
		TF1* l_1, *l_2;
	friend std::ostream& operator<< (std::ostream& os, const H_RecRPObject& rp);
};

#endif

