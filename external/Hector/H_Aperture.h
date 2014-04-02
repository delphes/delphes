#ifndef _H_Aperture_
#define _H_Aperture_

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

/// \file H_Aperture.h
/// \brief Defines the aperture of beamline elements.

// C++ #defines
#include <string>
#include <iostream>
using namespace std;

// local #defines
enum {NONE=0, RECTANGULAR, ELLIPTIC, CIRCULAR, RECTELLIPSE};
#define NONENAME        "None       "
#define RECTANGULARNAME "Rectangle  "
#define ELLIPTICNAME    "Ellipse    "
#define CIRCULARNAME    "Circle     "
#define RECTELLIPSENAME "Rectellipse"

/// Defines the aperture of any optics element.
class H_Aperture {

	public:
		/// Constructors, destructor and operators
		//@{
		H_Aperture();
		H_Aperture(const int, const float, const float, const float, const float, const float, const float);
		H_Aperture(const H_Aperture&);
		H_Aperture& operator=(const H_Aperture&);
		virtual ~H_Aperture() { };
		virtual H_Aperture* clone() const { return new H_Aperture(type_,x1,x2,x3,x4,fx,fy); };
		//@}

		/// Prints the aperture features
		//virtual void printProperties() const;
		void printProperties() const;
		/// Draws the aperture shape
		virtual void draw(const float) const {return;};
		/// Sets the (x,y) position in [m]
		void setPosition(const float,const float);
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const; 
		/// Returns the (int) type of aperture
		inline int getType() const {return type_;};
		/// Returns the (string) type of the aperture
		inline const string getTypeString() const { return aptypestring; }
		
	
	protected:
		/// Aperture shape (either RECTANGULAR or ELLIPTIC or ...)
		int type_;
		/// Aperture geometrical sizes (length/width or great/small radii) [m]
		//@{
		float x1, x2, x3, x4; 
		//@}
		/// Horizontal coordinate of the aperture center [m] (from the nominal beam position).
		//@{
		float fx, fy;
		//@}
                /// Aperture shape string
                string aptypestring;
		// Sets the name of the aperture from its type.
		//void setApertureString();
		/// Gets the name of the aperture from its type.
		const string getApertureString() const;

	friend std::ostream& operator<< (std::ostream& os, const H_Aperture& ap);
};


#endif
