#ifndef _H_AbstractBeamLine_
#define _H_AbstractBeamLine_

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
*	http://www.fynu.ucl.ac.be/hector.html               *
*                                                           *
* Center for Cosmology, Particle Physics and Phenomenology  *
*              Universite catholique de Louvain             *
*                 Louvain-la-Neuve, Belgium                 *
 *                                                         *
   * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/// \file H_AbstractBeamLine.h
/// \brief Class aiming at simulating the LHC beamline.
///
/// Units : angles [urad], distances [um], energies [GeV], c=[1].

	/// default length of the beam line
#define LENGTH_DEF 100

// c++ #includes
#include <vector>
#include <algorithm>
#include <string>


// ROOT #includes
#include "TMatrix.h"
#include "TGraph.h"

// local #includes
#include "H_OpticalElement.h"
using namespace std;

/// Beamline made from a vector of H_OpticalElement.
class H_AbstractBeamLine {

	public:
		///     Constructors, destructor and operator
		//@{
		H_AbstractBeamLine() {init(LENGTH_DEF);};
		H_AbstractBeamLine(const float length) {init(length);};
		H_AbstractBeamLine(const H_AbstractBeamLine &);
		H_AbstractBeamLine& operator=(const H_AbstractBeamLine&);
		H_AbstractBeamLine* clone() const ;
		~H_AbstractBeamLine();
		//@}
		///     Adds an element to the beamline
		void add(H_OpticalElement *);
		///     Returns the (float) length of the beamline
  		const float getLength() const { return beam_length;};
		///     Returns the (int) number of optics element of the beamline, including drifts
  		const unsigned int getNumberOfElements() const { return elements.size();}; 
		///     Returns the transport matrix for the whole beam
		const TMatrix getBeamMatrix() const;
		///     Returns the transport matrix for the whole beam, for given energy loss/mass/charge
		const TMatrix getBeamMatrix(const float , const float, const float ) ;
		///     Returns the transport matrix for a part of the beam from the IP to a given element
		const TMatrix getPartialMatrix(const H_OpticalElement *) const;
		///     Returns the transport matrix for a part of the beam from the IP to the ith element
		const TMatrix getPartialMatrix(const unsigned int ) const;
		///     Returns the transport matrix for a part of the beam from the IP to a given element, given energy loss/mass/charge
		const TMatrix getPartialMatrix(const string&, const float, const float, const float);
		///     Returns the ith element of the beamline
		//@{
		H_OpticalElement * getElement(const unsigned int );
		H_OpticalElement * getElement(const unsigned int ) const;
		//@}
		///	Returns a given element of the beamline, choosen by name 
		//@{
		H_OpticalElement * getElement(const string& );
		H_OpticalElement * getElement(const string& ) const;
		//@}
		///	Print some info
		void printProperties() const;
		///     Prints the element list
		void showElements() const;
		///     Prints the list of elements of a give type
		void showElements(const int) const;
		///     Prints the transport matrix for the whole beam
		void showMatrix() const;
		///     Prints the transport matrix for each element of the beamline
		void showMatrices() const;
		/// 	Reorders elements and adds the drift sections
		void calcSequence();
		/// 	Computes global transport matrix
		void calcMatrix();
		/// 	Draws the legend of the elements of the beam
		void draw(const float xmin =0.85, const float ymin=0.5, const float xmax=1, const float ymax=1) const;
		/// 	Draws the elements of the beam in the (x,s) plane
		void drawX(const float, const float, const float scale=1) const;
		/// 	Draws the elements of the beam in the (y,s) plane
		void drawY(const float, const float) const;
		///	Moves an element in the list, reorders the lists and recomputes the transport matrix
		void moveElement(const string&, const float ); 
		/// Moves the given element tranversely by given amounts.
		void alignElement(const string&, const float, const float);
		/// Tilts the given element tranversely by given angles.
		void tiltElement(const string&, const float, const float);
		///	Offsets all element in X pos from the start position
		void offsetElements(const float start, const float offset);
		///	Draws the beta functions, from MAD
		//@{
		TGraph * getBetaX() const;
		TGraph * getBetaY() const;
		//@}
                ///     Draws the dispersion functions, from MAD
                //@{
                TGraph * getDX() const;
                TGraph * getDY() const;
                //@}
                ///     Draws the relative position functions, from MAD
                //@{
                TGraph * getRelX() const;
                TGraph * getRelY() const;
                //@}


	private:
        	/// list of matrices, 1 matrix = the transport till the end of each element
		vector<TMatrix> matrices; 
	        /// transport matrix for the whole beam
		TMatrix beam_mat;
                /// list of all optics elements, including drifts
                vector<H_OpticalElement*> elements;
		/// Orderting method for the vector of H_OpticalElement*
		struct ordering{ bool operator()(const H_OpticalElement* el1, const H_OpticalElement* el2) const {return (*el1 < *el2);}};
		// private method for copying the contents of "elements"
		void cloneElements(const H_AbstractBeamLine&);

	protected:
		void init(const float );
		/// total length of the beamline
		float beam_length;
	friend std::ostream& operator<< (std::ostream& os, const H_AbstractBeamLine& be);
};

#endif
