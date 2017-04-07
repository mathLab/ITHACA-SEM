#ifndef H_PLYLIB  

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <complex>

/*
 *  LIBRARY ROUTINES FOR POLYNOMIAL CALCULUS AND INTERPOLATION
 */

/** \brief The namespace associated with the the Polylib library 
 * (\ref pagePolylib "Polylib introduction")
 */
namespace Polylib {

    /**
       \page pagePolylib The Polylib library
       \section sectionPolyLib Routines For Orthogonal Polynomial Calculus and Interpolation

       Spencer Sherwin, 
       Aeronautics, Imperial College London

       Based on codes by Einar Ronquist and Ron Henderson

       Abbreviations
       - z    -   Set of collocation/quadrature points
       - w    -   Set of quadrature weights
       - D    -   Derivative matrix
       - h    -   Lagrange Interpolant
       - I    -   Interpolation matrix
       - g    -   Gauss
	   - k	  -   Kronrod
       - gr   -   Gauss-Radau
       - gl   -   Gauss-Lobatto
       - j    -   Jacobi
       - m    -   point at minus 1 in Radau rules
       - p    -   point at plus  1 in Radau rules

       -----------------------------------------------------------------------\n
       MAIN     ROUTINES\n
       -----------------------------------------------------------------------\n

       Points and Weights:

       - zwgj        Compute Gauss-Jacobi         points and weights
       - zwgrjm      Compute Gauss-Radau-Jacobi   points and weights (z=-1)
       - zwgrjp      Compute Gauss-Radau-Jacobi   points and weights (z= 1)
       - zwglj       Compute Gauss-Lobatto-Jacobi points and weights
	   - zwgk        Compute Gauss-Kronrod-Jacobi points and weights
	   - zwrk        Compute Radau-Kronrod        points and weights
	   - zwlk        Compute Lobatto-Kronrod      points and weights
	   - JacZeros    Compute Gauss-Jacobi         points and weights

       Derivative Matrices:

       - Dgj         Compute Gauss-Jacobi         derivative matrix
       - Dgrjm       Compute Gauss-Radau-Jacobi   derivative matrix (z=-1)
       - Dgrjp       Compute Gauss-Radau-Jacobi   derivative matrix (z= 1)
       - Dglj        Compute Gauss-Lobatto-Jacobi derivative matrix

       Lagrange Interpolants:

       - hgj         Compute Gauss-Jacobi         Lagrange interpolants
       - hgrjm       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z=-1)
       - hgrjp       Compute Gauss-Radau-Jacobi   Lagrange interpolants (z= 1)
       - hglj        Compute Gauss-Lobatto-Jacobi Lagrange interpolants

       Interpolation Operators:

       - Imgj        Compute interpolation operator gj->m
       - Imgrjm      Compute interpolation operator grj->m (z=-1)
       - Imgrjp      Compute interpolation operator grj->m (z= 1)
       - Imglj       Compute interpolation operator glj->m

       Polynomial Evaluation:

       - jacobfd     Returns value and derivative of Jacobi poly. at point z
       - jacobd      Returns derivative of Jacobi poly. at point z (valid at z=-1,1)

       -----------------------------------------------------------------------\n
       LOCAL      ROUTINES\n
       -----------------------------------------------------------------------\n

       - jacobz      Returns Jacobi polynomial zeros
       - gammaf      Gamma function for integer values and halves
	   - RecCoeff    Calculates the recurrence coefficients for orthogonal poly
	   - TriQL		 QL algorithm for symmetrix tridiagonal matrix
	   - JKMatrix	 Generates the Jacobi-Kronrod matrix

       ------------------------------------------------------------------------\n

       Useful references:

       - [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
       Providence, Rhode Island, 1939.
       - [2] Abramowitz \& Stegun: Handbook of Mathematical Functions,
       Dover, New York, 1972.
       - [3] Canuto, Hussaini, Quarteroni \& Zang: Spectral Methods in Fluid
       Dynamics, Springer-Verlag, 1988.
       - [4] Ghizzetti \& Ossicini: Quadrature Formulae, Academic Press, 1970.
       - [5] Karniadakis \& Sherwin: Spectral/hp element methods for CFD, 1999


       NOTES

       -# Legendre  polynomial \f$ \alpha = \beta = 0 \f$ 
       -# Chebychev polynomial \f$ \alpha = \beta = -0.5 \f$
       -# All routines are double precision.  
       -# All array subscripts start from zero, i.e. vector[0..N-1] 
    */

    /*-----------------------------------------------------------------------
      M A I N     R O U T I N E S
      -----------------------------------------------------------------------*/

    /* Points and weights */
    LIB_UTILITIES_EXPORT void   zwgj    (double *, double *, const int , const double, const double);
    LIB_UTILITIES_EXPORT void   zwgrjm  (double *, double *, const int , const double, const double);
    LIB_UTILITIES_EXPORT void   zwgrjp  (double *, double *, const int , const double, const double);
    LIB_UTILITIES_EXPORT void   zwglj   (double *, double *, const int , const double, const double);
	LIB_UTILITIES_EXPORT void   zwgk    (double *, double *, const int , const double, const double);
	LIB_UTILITIES_EXPORT void   zwrk    (double *, double *, const int , const double, const double);
	LIB_UTILITIES_EXPORT void   zwlk    (double *, double *, const int , const double, const double);
	LIB_UTILITIES_EXPORT void   JacZeros(const int, double *, double*, const double,const double);


    /* Derivative operators */
    LIB_UTILITIES_EXPORT void   Dgj     (double *, const double *, const int, const double, 
                    const double);
    LIB_UTILITIES_EXPORT void   Dgrjm   (double *, const double *, const int, const double, 
                    const double);
    LIB_UTILITIES_EXPORT void   Dgrjp   (double *, const double *, const int, const double, 
                    const double);
    LIB_UTILITIES_EXPORT void   Dglj    (double *,const double *, const int, const double,
                    const double);

    /* Lagrangian interpolants */
    LIB_UTILITIES_EXPORT double hgj     (const int, const double, const double *, const int, 
                    const double, const double);
    LIB_UTILITIES_EXPORT double hgrjm   (const int, const double, const double *, const int, 
                    const double, const double);
    LIB_UTILITIES_EXPORT double hgrjp   (const int, const double, const double *, const int, 
                    const double, const double);
    LIB_UTILITIES_EXPORT double hglj    (const int, const double, const double *, const int, 
                    const double, const double);

    /* Interpolation operators */
    LIB_UTILITIES_EXPORT void  Imgj  (double*, const double*, const double*, const int, const int, 
                 const double, const double);
    LIB_UTILITIES_EXPORT void  Imgrjm(double*, const double*, const double*, const int, const int,
                 const double, const double);
    LIB_UTILITIES_EXPORT void  Imgrjp(double*, const double*, const double*, const int, const int, 
                 const double, const double);
    LIB_UTILITIES_EXPORT void  Imglj (double*, const double*, const double*, const int, const int, 
                 const double, const double);

    /* Polynomial functions */
    LIB_UTILITIES_EXPORT void jacobfd (const int, const double *, double *, double *, const int , 
                  const double, const double);
    LIB_UTILITIES_EXPORT void jacobd  (const int, const double *, double *,  const int , 
                  const double, const double);

    LIB_UTILITIES_EXPORT std::complex<Nektar::NekDouble> ImagBesselComp(const int,
                  std::complex<Nektar::NekDouble> );


} // end of namespace


#define H_PLYLIB
#endif          /* END OF POLYLIB.H DECLARATIONS */









