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

/// \file H_RectEllipticAperture.cc
/// \brief Defines the Rect-Elliptic aperture of beamline elements.

// C++ #includes
#include <iostream>

// C #includes
#include <cmath> // needed for fabs

// ROOT #includes
#include "TPolyLine.h"

// local #includes
#include "H_RectEllipticAperture.h"
using namespace std;

H_RectEllipticAperture::H_RectEllipticAperture(const float l, const float h, const float L, const float H, const float posx, const float posy) :H_Aperture(RECTELLIPSE,((l==0)?L:l),((h==0)?H:h),L,H,posx,posy) {
	/// @param l, h, L, H are the geometrical parameters of the rect-ellipse shape
	/// @param posx, posy defines the (x,y) of the center of the shape
}

H_RectEllipticAperture* H_RectEllipticAperture::clone() const {
	return new H_RectEllipticAperture(x1,x2,x3,x4,fx,fy);
}


TPolyLine * rectellipse(const float a_e = 2, const float b_e = 1, const float a_r = 1, const float b_r = 2, const float center_x = 0, const float center_y =0) {
        const int n = 20;  // number of points per segment
        const int N = 4*n; // there are 4 segments
        float x[N+1], y[N+1];

   if(a_e>a_r) {
        // a rectellipse has 4 segments
        // 1) upper one
        for (int i=0; i<n; i++) {
                x[i] = -a_r + i*(2*a_r)/(float)n;
                y[i] =  b_e * sqrt(1-pow(x[i]/a_e,2));
        }

        // 2) right vertical segment
        // upper right corner
        const float y2 = b_e * sqrt(1-pow(a_r/a_e,2));
        // lower right corner
        const float y3 = -b_e * sqrt(1-pow(a_r/a_e,2));
        for (int i=n; i<2*n; i++) {
                x[i] = a_r;
                y[i] = y2 - (i-n)*(2*y2)/(float)n;
        }

        // 3) lower side
        for (int i=2*n; i<3*n; i++) {
                x[i] = a_r - (i-2*n)*(2*a_r)/(float)n;
                y[i] = -b_e * sqrt(1-pow(x[i]/a_e,2));
        }

        // 4) left vertical segment
        // lower left corner
        const float y4 = y3;
        for (int i=3*n; i<4*n; i++) {
                x[i] = -a_r;
                y[i] = y4 + (i-3*n)*(2*y2)/(float)n;
        }
   } else {
	// 1) upper one : flat
	const float x1 = -a_e * sqrt(1-pow(b_r/b_e,2));
	const float x2 = -x1;
	for (int i=0; i<n; i++) {
		y[i] = b_r;
		x[i] = x1 + i * (x2-x1)/(float)n;
	}
	
	// 2) right curved border
	for (int i=n; i<2*n; i++) {
		y[i] = b_r - (i-n) * (2*b_r)/(float)n;
		x[i] = a_e * sqrt(1-pow(y[i]/b_e,2));
	}

	// 3) lower side : flat
	for (int i=2*n; i<3*n; i++) {
		y[i] = -b_r;
		x[i] = x2 - (i-2*n) * (2*x2)/(float)n;
	}
	
	// 4) left curved border
	for (int i=3*n; i<4*n; i++) {
		y[i] = -b_r + (i-3*n) * (2*b_r)/(float)n;
		x[i] = -a_e * sqrt(1-pow(y[i]/b_e,2));
	}
   }

	// closing the polyline
        x[N] = x[0];
        y[N] = y[0];

	// shifting the center
        for (int i=0; i<N+1; i++) {
                x[i] += center_x;
                y[i] += center_y;
        }

        return new TPolyLine(N+1,x,y);
}


void H_RectEllipticAperture::draw(const float scale) const {
	TPolyLine * re = rectellipse(x3*scale, x4*scale, x1*scale, x2*scale, fx*scale, fy*scale);
	re->SetLineColor(39);
	re->SetLineWidth(2);
	re->Draw("l");
	return;
}

bool H_RectEllipticAperture::isInside(const float x, const float y) const {
	return((fabs(fx-x)<x1)&&(fabs(fy-y)<x2)&&(((x-fx)/x3)*((x-fx)/x3)+((y-fy)/x4)*((y-fy)/x4)<1));
}

std::ostream& operator<< (std::ostream& os, const H_RectEllipticAperture& ap) {
	os << "Aperture shape:" << ap.aptypestring << ", parameters " << ap.x1 <<", "<< ap.x2 <<", "<< ap.x3 <<", "<< ap.x4 << endl;
	os << " \t Center : " << ap.fx <<", "<< ap.fy <<endl;
	return os;
}
