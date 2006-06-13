#ifndef H_PLYLIB  
/*
 *  LIBRARY ROUTINES FOR POLYNOMIAL CALCULUS AND INTERPOLATION
 */



namespace Polylib {


/*-----------------------------------------------------------------------
                         M A I N     R O U T I N E S
  -----------------------------------------------------------------------*/

/* Points and weights */
void   zwgj    (double *, double *, const int , const double, const double);
void   zwgrjm  (double *, double *, const int , const double, const double);
void   zwgrjp  (double *, double *, const int , const double, const double);
void   zwglj   (double *, double *, const int , const double, const double);

/* Derivative operators */
void   Dgj     (double *, const double *, const int, const double, 
		const double);
void   Dgrjm   (double *, const double *, const int, const double, 
		const double);
void   Dgrjp   (double *, const double *, const int, const double, 
		const double);
void   Dglj    (double *,const double *, const int, const double,
		const double);

/* Lagrangian interpolants */
double hgj     (const int, const double, const double *, const int, 
		const double, const double);
double hgrjm   (const int, const double, const double *, const int, 
		const double, const double);
double hgrjp   (const int, const double, const double *, const int, 
		const double, const double);
double hglj    (const int, const double, const double *, const int, 
		const double, const double);

/* Interpolation operators */
void  Imgj  (double*, const double*, const double*, const int, const int, 
	     const double, const double);
void  Imgrjm(double*, const double*, const double*, const int, const int,
	     const double, const double);
void  Imgrjp(double*, const double*, const double*, const int, const int, 
	     const double, const double);
void  Imglj (double*, const double*, const double*, const int, const int, 
	     const double, const double);

/* Polynomial functions */
void jacobfd (const int, const double *, double *, double *, const int , 
	      const double, const double);
void jacobd  (const int, const double *, double *,  const int , 
	      const double, const double);


} // end of namespace


#define H_PLYLIB
#endif          /* END OF POLYLIB.H DECLARATIONS */









