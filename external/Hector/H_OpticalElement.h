#ifndef _H_OpticalElement_
#define _H_OpticalElement_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/// \file H_OpticalElement.h
/// \brief Class aiming at describing any beam optical element.
///
/// This is an abstract base class
/// subclass must define setTypeString() and setMatrix() methods.
/// Quadrupole, Dipole, Kickers and Drift inherit from this class.

// c++ #includes
#include <iostream>
#include <string>

// ROOT #includes
#include "TMatrix.h"

// local #includes
#include "H_TransportMatrices.h"
#include "H_Parameters.h"
#include "H_Aperture.h"

using namespace std;

	// type #defines
#define DRIFT        1
#define RDIPOLE      2
#define SDIPOLE      3
#define VQUADRUPOLE  4
#define HQUADRUPOLE  5
#define VKICKER      6
#define HKICKER      7
#define RCOLLIMATOR  8
#define ECOLLIMATOR  9
#define CCOLLIMATOR 10
#define RP          11
#define IP          12
#define MARKER      13

	// typestring[30] #defines
#define DRIFTNAME	"Drift        "
#define RDIPOLENAME	"R-Dipole     "
#define SDIPOLENAME	"S-Dipole     "
#define VQUADRUPOLENAME	"V-Quadrupole "
#define HQUADRUPOLENAME	"H-Quadrupole "
#define VKICKERNAME     "V-Kicker     "
#define HKICKERNAME     "H-Kicker     "
#define RCOLLIMATORNAME "R-Collimator "
#define ECOLLIMATORNAME "E-Collimator "
#define CCOLLIMATORNAME "C-Collimator "
#define RPNAME          "Roman Pot    "
#define IPNAME          "IP           "
#define MARKERNAME      "Marker       "

/// Describes any beam optical element.
class H_OpticalElement {

	public:
        ///     init method for constructors
                void init(const string, const int , const double , const double , const double, H_Aperture*);
		/// Constructors and destructor
		//@{
		H_OpticalElement(const string, const int, const double, const double, const double, H_Aperture*);
		H_OpticalElement(const int, const double, const double, const double, H_Aperture*);
		H_OpticalElement(const string, const int, const double, const double, const double);
		H_OpticalElement(const int, const double, const double, const double);
		H_OpticalElement();
		H_OpticalElement(const H_OpticalElement&);
	        virtual ~H_OpticalElement() {delete element_mat; delete element_aperture;};
		//@}
		///     Prints the element features
	        virtual void printProperties() const ;
		///     Shows the element transport matrix
	        void showMatrix() const ;
		///     Draws the aperture shape
	        void drawAperture() const ;
		///     Sets the aperture of the element
		void setAperture(H_Aperture *);
		///     Ordering operator acting on the s coordinate
		inline bool operator>(const H_OpticalElement tocomp) const {if(fs>tocomp.getS()) { return true; } else { return false; }};
		///     Ordering operator acting on the s coordinate
		inline bool operator<(const H_OpticalElement tocomp) const {if(fs<tocomp.getS()) { return true; } else { return false; }};
		/// 	Copy operator
		H_OpticalElement& operator=(const H_OpticalElement&); 
		///     Sets the element longitudinal (s), and horizontal (x) and vertical (y) coordinates.
		//@{
		inline void setS(const double new_s) {fs=new_s;};
		inline void setX(const double new_pos) {
		        /// @param new_pos in [m]
	        	element_aperture->setPosition(new_pos*URAD,ypos);
		        xpos = new_pos;
		};
		inline void setY(const double new_pos) {
		        /// @param new_pos in [m]
		        element_aperture->setPosition(xpos,new_pos*URAD);
		        ypos = new_pos;
		};
		inline void setTX(const double new_ang) {
			txpos = new_ang;
		};
		inline void setTY(const double new_ang) {
			typos = new_ang;
		}
		//@}
		///     Returns the element longitudinal (s), horizontal (x) and vertical (y) coordinates, and corresponding angles.
		//@{
	        inline double getS() const {return fs;};
		inline double getX() const {return xpos;};
		inline double getY() const {return ypos;};
		inline double getTX() const {return txpos;};
		inline double getTY() const {return typos;};
		//@}
		///     Returns the element length
	        inline double getLength() const { return element_length; };
		///     Returns the element magnetic strength
	        inline double getK() const { return fk; };
		///     Returns the element type, as (int)  or (string)
		//@{
	        inline int getType() const { return type; };
		inline const string getTypeString() const { return typestring; };
		//@}
		///     Returns the element (string) name
		inline const string getName() const { return name; };
		///     Draws the element from min to max in the current pad
		void draw(const float, const float) const;
		///     Checks if the (x,y) coordinates are within the aperture acceptance
		inline bool isInside(const double x, const double y) const { return (bool) element_aperture->isInside(x,y);};
		///     Returns the element transport matrix
		//@{
		TMatrix getMatrix() const;
		TMatrix getMatrix(const float, const float, const float) const;
		//@}
		///     Returns the element aperture
		H_Aperture * getAperture() const {return element_aperture;};
		///	Sets the beta functions
		//@{
		inline void setBetaX(const double beta) { betax = beta;};
		inline void setBetaY(const double beta) { betay = beta;};
		//@}
		///	Returns the beta functions
		//@{
		inline double getBetaX() const {return betax;};
		inline double getBetaY() const {return betay;};
		//@}
                ///     Sets the dispersion functions
                //@{
                inline void setDX(const double disp) { dx = disp;};
                inline void setDY(const double disp) { dy = disp;};
                //@}
                ///     Returns the dispersion functions
                //@{
                inline double getDX() const {return dx;};
                inline double getDY() const {return dy;};
                //@}
                //@}
                ///     Sets the relative position functions (from MAD), relative to the beamline reference frame
                //@{
                inline void setRelX(const double rel_x) { relx = rel_x;};
                inline void setRelY(const double rel_y) { rely = rel_y;};
                //@}
                ///     Returns the position functions (from MAD), relative to the beamline reference frame
                //@{
                inline double getRelX() const {return relx;};
                inline double getRelY() const {return rely;};
                //@}


	protected:
		/// Optical element coordinates : fs [m], xpos [m], ypos [m], txpos [rad], typos [rad].
		//@{
        	double fs, xpos, ypos, txpos, typos;
		//@}
		/// Optical element lenght.
		double element_length;
		/// Magnetic field strength.
		double fk;
		/// Beam \f$ \beta \f$ functions.
		//@{
		double betax, betay;
		//@}
		/// Beam dispersion in position 
		//@{
		double dx, dy;
		//@}
		/// Beam position, relative to the beam ideal path, from MAD
		//@{
		double relx, rely;
		//@}
		/// Optical element type (Dipole, drift, etc).
     		int type;
		/// Optical element name and type.
		//@{
		string name, typestring;
		//@}
		virtual void setTypeString() {return;};
		/// Optical element transport matrix.
	    	virtual void setMatrix(const float, const float, const float) const {return;};
		/// Optical element transport matrix.
 		TMatrix * element_mat;
		/// Optical element aperture.
	    	H_Aperture * element_aperture;
};

#endif
