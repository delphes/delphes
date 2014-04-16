#ifndef _H_Aperture_
#define _H_Aperture_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_Aperture.h
/// \brief Defines the aperture of beamline elements.

// C++ #defines
#include <string>
using namespace std;

// local #defines
#define NONE        0
#define RECTANGULAR 1
#define ELLIPTIC    2
#define CIRCULAR    3
#define RECTELLIPSE 4
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
		virtual ~H_Aperture() {return;};
		//@}

		/// Prints the aperture features
		virtual void printProperties() const;
		/// Draws the aperture shape
		virtual void draw() const {return;};
		/// Sets the (x,y) position in [m]
		void setPosition(const float,const float);
		/// Checks whether the point is inside the aperture or not
		virtual bool isInside(const float, const float) const; 
		/// Returns the (int) type of aperture
		inline int getType() const {return type;};
		/// Returns the (string) type of the aperture
		inline const string getTypeString() const { return aptypestring; }
		
	
	protected:
		/// Aperture shape (either RECTANGULAR or ELLIPTIC or ...)
		int type;
		/// Aperture shape string
		string aptypestring;

		/// Aperture geometrical sizes (length/width or great/small radii) [m]
		//@{
		float x1, x2, x3, x4; 
		//@}

		/// Horizontal coordinate of the aperture center [m] (from the nominal beam position).
		//@{
		float fx, fy;
		//@}
		/// Sets the name of the aperture from its type.
		void setApertureString();
};

#endif
