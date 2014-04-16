#ifndef _H_TransportMatrices_
#define _H_TransportMatrices_

/*
---- Hector the simulator ----
   A fast simulator of particles through generic beamlines.
   J. de Favereau, X. Rouby ~~~ hector_devel@cp3.phys.ucl.ac.be

        http://www.fynu.ucl.ac.be/hector.html

   Centre de Physique des Particules et de Phénoménologie (CP3)
   Université Catholique de Louvain (UCL)
*/

/** \file H_TransportMatrices.h
 * \brief Contains the matrices defining the propagation of the beam.
 * 
 * The matrices should have the following units :
 *	\f$
 *	\mathbf{M} =
 *	\left(
 *	\begin{array}{cccccc}
 *	    1 & 1/m & 1 & 1/m & GeV/m & 1 \\
 *	    m & 1 & m & 1 & GeV & 1 \\
 *	    1 & 1/m & 1 & 1/m & GeV/m & 1\\
 *	    m & 1 & m & 1 & GeV & 1 \\
 *	  m/GeV & 1/GeV & m/GeV & 1/GeV & 1 & 1\\
 *	    1 & 1 & 1 & 1 & 1  & 1 \\	
 *	\end{array}
 *	\right)
 *	\f$
 *
 *	Note : convention is transposed compared to ref : x0.M = x1
 * 	instead of x1 = M.x0 so the matrices should be transposed.
 */

// ROOT #includes
#include "TMatrix.h"
using namespace std;

	/// transport matrix dimension
#define MDIM 6

/// \f$ \omega(k,l) = l \sqrt{|k|} \f$ is needed for the quadrupole matrices
extern double omega(const double , const double );

/// \f$ r(k) \f$ is needed for the dipole matrices
extern double radius(const double );

/// Prints the matrix
extern void printMatrix(TMatrix * );

/// \brief Returns the matrix for a vertically focussing quadrupole (H_VerticalQuadrupole)

/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 \cosh(\omega) & \sqrt{k}\sinh(\omega) & 0 & 0 & 0 & 0\\
 (1/\sqrt{k})sinh(\omega) & \cosh(\omega) & 0 & 0 & 0 &0\\
 0 & 0 & \cos(\omega) & -\sqrt{k}sin(\omega) & 0 &0\\
 0 & 0 & (1/\sqrt{k})*sin(\omega) & \cos(\omega) & 0 &0\\
 0 & 0 & 0 & 0 & 1 &0\\
 0 & 0 & 0 & 0 & 0 &1 \\
 \end{array}
 \right)
 \f$

 assuming \f$ k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$ and \f$ \omega(k,l) = l \sqrt{|k|} \f$
*/

extern TMatrix vquadmat(const float , const float , const float , const float , const float);

/// \brief Returns the matrix for a horizontally focussing quadrupole (H_HorizontalQuadrupole)

/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 \cos(\omega) & -\sqrt{k}\sin(\omega) & 0 & 0 & 0 & 0\\
 (1/\sqrt{k})sin(\omega) & \cos(\omega) & 0 & 0 & 0 & 0\\
 0 & 0 & \cosh(\omega) & \sqrt{k}sinh(\omega) & 0 & 0\\
 0 & 0 & (1/\sqrt{k})sinh(\omega) & \cosh(\omega) & 0 & 0\\
 0 & 0 & 0 & 0 & 1 & 0\\
 0 & 0 & 0 & 0 & 0 & 1 \\
 \end{array}
 \right)
 \f$

 assuming \f$ k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$ and \f$ \omega(k,l) = l \sqrt{|k|} \f$
*/
extern TMatrix hquadmat(const float , const float , const float , const float , const float);

/// \brief Returns the matrix for a rectangle dipole (H_RectangularDipole)

/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 \cos(l/r) & \frac{-1}{r} \sin(l/r) & 0 & 0 & 0 & 0\\
 r \sin(l/r) & \cos(l/r) & 0 & 0 & 0 & 0\\
 0 & 0 & 1 & 0 & 0 &0\\
 0 & 0 & l & 1 & 0 &0\\
 2r \sin^2(l/2r)/BE & \sin(l/r)/BE & 0 & 0 & 1 &0\\
 0 & 0 & 0 & 0 & 0 & 1\\
 \end{array}
 \right)
 \f$

 assuming \f$ 1/r = k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$ and \f$ BE = 7000 GeV \f$.

Attention : numerical sensitivity with \f$ r*(1-\cos(l/r))/BE\f$. \\
Using \f$ 2\sin^2(x/2) = 1-\cos(x)\f$ instead (see the variable called "simp")
*/
extern TMatrix rdipmat(const float, const float , const float , const float , const float);

/// \brief Returns the matrix for a sector dipole (H_SectorDipole)

/*! The matrix is different if the bending is on or off.
 
 \f$
 \mathbf{M_{bending-off}} =
 \left(
 \begin{array}{cccccc}
 \cos(l/r) & \frac{-1}{r} \sin(l/r) & 0 & 0 & 0 & 0\\
 r \sin(l/r) & \cos(l/r) & 0 & 0 & 0 & 0\\
 0 & 0 & 1 & 0 & 0 & 0\\
 0 & 0 & l & 1 & 0 & 0\\
 0 & 0 & 0 & 0 & 1 & 0\\
 0 & 0 & 0 & 0 & 0 & 1\\
 \end{array}
 \right)
 \f$

 \f$
 \mathbf{M_{bending-on}} =
 \left( 
 \begin{array}{cccccc} 
 \cos(l/r) & \frac{-1}{r} \sin(l/r) & 0 & 0 & 0 & 0\\
 r \sin(l/r) & \cos(l/r) & 0 & 0 & 0 & 0\\
 0 & 0 & 1 & 0 & 0 &0\\ 
 0 & 0 & l & 1 & 0 &0\\ 
 2r \sin^2(l/2r)/BE & \sin(l/r)/BE & 0 & 0 & 1 & 0\\ 
 0 & 0 & 0 & 0 & 0 & 1\\
 \end{array} 
 \right) 
 \f$ 
 
 assuming \f$ 1/r = k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$

*/
extern TMatrix sdipmat(const float, const float , const float , const float , const float );

/// \brief Returns the matrix for a drift (H_Drift)

/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 1 & 0 & 0 & 0 & 0 & 0\\
 l & 1 & 0 & 0 & 0 & 0\\
 0 & 0 & 1 & 0 & 0 & 0\\
 0 & 0 & l & 1 & 0 & 0\\
 0 & 0 & 0 & 0 & 1 & 0\\
 0 & 0 & 0 & 0 & 0 & 1\\
 \end{array}
 \right)
 \f$
*/
extern TMatrix driftmat(const float );

/// \brief Returns the matrix for a horizontal kicker (H_HorizontalKicker)
/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 1 & 0 & 0 & 0 & 0 & 0\\
 l & 1 & 0 & 0 & 0 & 0 \\
 0 & 0 & 1 & 0 & 0 & 0\\
 0 & 0 & l & 1 & 0 & 0\\
 0 & 0 & 0 & 0 & 1 & 0\\
 l \tan(k) /2 & k & 0 & 0 & 0 & 1\\
 \end{array}
 \right) 
 \f$
 
 assuming \f$ k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$
*/
extern TMatrix hkickmat(const float, const float , const float , const float, const float);

/// \brief Returns the matrix for a vertical kicker (H_VerticalKicker)
/*! \f$
 \mathbf{M} =
 \left(
 \begin{array}{cccccc}
 1 & 0 & 0 & 0 & 0 & 0\\
 l & 1 & 0 & 0 & 0 & 0 \\
 0 & 0 & 1 & 0 & 0 & 0\\
 0 & 0 & l & 1 & 0 & 0\\
 0 & 0 & 0 & 0 & 1 & 0\\
 0 & 0 & l \tan(k) /2 & k & 0 & 1\\
 \end{array}
 \right) 
 \f$

 assuming \f$ k =  k_{0} \times \frac{p_{0}}{p_{0} - dp} \times \frac{q_{particle}}{q_{proton}} \f$
*/
extern TMatrix vkickmat(const float, const float , const float , const float , const float);



#endif
