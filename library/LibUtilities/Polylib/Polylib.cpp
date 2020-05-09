#include "Polylib.h"
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <complex>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
/// Maximum number of iterations in polynomial defalation routine Jacobz
#define STOP  30
/// Precision tolerance for two points to be similar
#define EPS   100*DBL_EPSILON
/// return the sign(b)*a
#define sign(a,b) ((b)<0 ? -fabs(a) : fabs(a))

namespace Polylib {

	/// The following function is used to circumvent/reduce "Subtractive Cancellation"
	/// The expression 1/dz  is replaced by optinvsub(.,.)
	/// Added on 26 April 2017
	double optdiff(double xl, double xr)
        {
		double m_xln, m_xrn;
                int    m_expn;
                int    m_digits = static_cast<int>(fabs(floor(log10(DBL_EPSILON)))-1);

                if (fabs(xl-xr)<1.e-4){

                        m_expn = static_cast<int>(floor(log10(fabs(xl-xr))));
                        m_xln  = xl*powl(10.0L,-m_expn)-floor(xl*powl(10.0L,-m_expn)); // substract the digits overlap part
                        m_xrn  = xr*powl(10.0L,-m_expn)-floor(xl*powl(10.0L,-m_expn)); // substract the common digits overlap part
                        m_xln  = round(m_xln*powl(10.0L,m_digits+m_expn));             // git rid of rubbish
                        m_xrn  = round(m_xrn*powl(10.0L,m_digits+m_expn));

                        return powl(10.0L,-m_digits)*(m_xln-m_xrn);
                }else{
                        return (xl-xr);
                }
        }

	double laginterp(double z, int j, const double *zj, int np)
	{
        	double temp = 1.0;
        	for (int i=0; i<np; i++)
        	{
                	if (j != i)
                	{
                        	temp *=optdiff(z,zj[i])/(zj[j]-zj[i]);
                	}
        	}
        	return temp;
	}
    /// Define whether to use polynomial deflation (1)  or tridiagonal solver (0).
#define POLYNOMIAL_DEFLATION 0

#ifdef POLYNOMIAL_DEFLATION
    /// zero determination using Newton iteration with polynomial deflation
#define jacobz(n,z,alpha,beta) Jacobz(n,z,alpha,beta)
#else
    /// zero determination using eigenvalues of tridiagaonl matrix
#define jacobz(n,z,alpha,beta) JacZeros(n,z,alpha,beta)
#endif



    /* local functions */
    static void   Jacobz   (const int n, double *z, const double alpha,
        const double beta);
   // static void   JacZeros (const int n, double *a, const double alpha,
    //    const double beta);
    //static void   TriQL    (const int n, double *d, double *e);
	static void TriQL(const int, double *,double *, double **);
    double gammaF (const double);
	static void RecCoeff(const int, double *, double *,const double,
				  const double);
	void JKMatrix(int, double *, double *);
	void chri1(int,double*,double*,double*,double*,double);

    /**
    \brief  Gauss-Jacobi zeros and weights.

    \li Generate \a np Gauss Jacobi zeros, \a z, and weights,\a w,
    associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)
    \f$,

    \li Exact for polynomials of order \a 2np-1 or less
    */

    void zwgj (double *z, double *w, const int np, const double alpha,
        const double beta)
    {
        int i;
        double fac, one = 1.0, two = 2.0, apb = alpha + beta;

        jacobz (np,z,alpha,beta);
        jacobd (np,z,w,np,alpha,beta);

        fac  = pow(two,apb + one)*gammaF(alpha + np + one)*gammaF(beta + np + one);
        fac /= gammaF(np + one)*gammaF(apb + np + one);

        for(i = 0; i < np; ++i) w[i] = fac/(w[i]*w[i]*(one-z[i]*z[i]));

        return;
    }


    /**
    \brief  Gauss-Radau-Jacobi zeros and weights with end point at \a z=-1.

    \li Generate \a np Gauss-Radau-Jacobi zeros, \a z, and weights,\a w,
    associated with the  polynomial \f$(1+z) P^{\alpha,\beta+1}_{np-1}(z)
    \f$.

    \li  Exact for polynomials of order \a 2np-2 or less
    */

    void zwgrjm(double *z, double *w, const int np, const double alpha,
        const double beta)
    {

        if(np == 1){
            z[0] = 0.0;
            w[0] = 2.0;
        }
        else{
            int i;
            double fac, one = 1.0, two = 2.0, apb = alpha + beta;

            z[0] = -one;
            jacobz  (np-1,z+1,alpha,beta+1);
            jacobfd (np,z,w,NULL,np-1,alpha,beta);

            fac  = pow(two,apb)*gammaF(alpha + np)*gammaF(beta + np);
            fac /= gammaF(np)*(beta + np)*gammaF(apb + np + 1);

            for(i = 0; i < np; ++i) w[i] = fac*(1-z[i])/(w[i]*w[i]);
            w[0] *= (beta + one);
        }

        return;
    }


    /**
    \brief  Gauss-Radau-Jacobi zeros and weights with end point at \a z=1


    \li Generate \a np Gauss-Radau-Jacobi zeros, \a z, and weights,\a w,
    associated with the  polynomial \f$(1-z) P^{\alpha+1,\beta}_{np-1}(z)
    \f$.

    \li Exact for polynomials of order \a 2np-2 or less
    */

    void zwgrjp(double *z, double *w, const int np, const double alpha,
        const double beta)
    {

        if(np == 1){
            z[0] = 0.0;
            w[0] = 2.0;
        }
        else{
            int i;
            double fac, one = 1.0, two = 2.0, apb = alpha + beta;

            jacobz  (np-1,z,alpha+1,beta);
            z[np-1] = one;
            jacobfd (np,z,w,NULL,np-1,alpha,beta);

            fac  = pow(two,apb)*gammaF(alpha + np)*gammaF(beta + np);
            fac /= gammaF(np)*(alpha + np)*gammaF(apb + np + 1);

            for(i = 0; i < np; ++i) w[i] = fac*(1+z[i])/(w[i]*w[i]);
            w[np-1] *= (alpha + one);
        }

        return;
    }


    /**
    \brief  Gauss-Lobatto-Jacobi zeros and weights with end point at \a z=-1,\a 1


    \li Generate \a np Gauss-Lobatto-Jacobi points, \a z, and weights, \a w,
    associated with polynomial \f$ (1-z)(1+z) P^{\alpha+1,\beta+1}_{np-2}(z) \f$
    \li Exact for polynomials of order \a 2np-3 or less
    */

    void zwglj(double *z, double *w, const int np, const double alpha,
        const double beta)
    {

        if( np == 1 ){
            z[0] = 0.0;
            w[0] = 2.0;
        }
        else if( np == 2 ){
            z[0] = -1.0;
            z[1] =  1.0;

            w[0] =  1.0;
            w[1] =  1.0;
        }
        else{
            int i;
            double   fac, one = 1.0, apb = alpha + beta, two = 2.0;

            z[0]    = -one;
            z[np-1] =  one;
            jacobz  (np-2,z + 1,alpha + one,beta + one);
            jacobfd (np,z,w,NULL,np-1,alpha,beta);

            fac  = pow(two,apb + 1)*gammaF(alpha + np)*gammaF(beta + np);
            fac /= (np-1)*gammaF(np)*gammaF(alpha + beta + np + one);

            for(i = 0; i < np; ++i) w[i] = fac/(w[i]*w[i]);
            w[0]    *= (beta  + one);
            w[np-1] *= (alpha + one);
        }

        return;
    }

	/**
    \brief  Gauss-Kronrod-Jacobi zeros and weights.

    \li Generate \a npt=2*np+1 Gauss-Kronrod Jacobi zeros, \a z, and weights,\a w,
    associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)
    \f$,

    \li Exact for polynomials of order \a 3np+1 or less
	*/
	void zwgk(double *z, double *w, const int npt , const double alpha,
		const double beta)
	{

		int np = (npt-1)/2;

		int  i,j;

		// number of kronrod points associated with the np gauss rule
		int kpoints = 2*np + 1;

		// Define the number of required recurrence coefficents
		int ncoeffs = (int)floor(3.0*(np+1)/2);

		// Define arrays  for the recurrence coefficients
		// We will use these arrays for the Kronrod results too, hence the
		// reason for the size of the arrays
		double *a = new double[kpoints];
		double *b = new double[kpoints];

		// Initialize a and b to zero
		for(i = 0; i < kpoints; i++)
		{
			a[i] = 0.0;
			b[i] = 0.0;
		}

		// Call the routine to calculate the recurrence coefficients
		RecCoeff(ncoeffs,a,b,alpha,beta);

		// Call the routine to caluclate the jacobi-Kronrod matrix
		JKMatrix(np,a,b);

		// Set up the identity matrix
		double** zmatrix = new double*[kpoints];
		for(i = 0; i < kpoints; i++)
		{
			zmatrix[i] = new double[kpoints];
			for(j = 0; j < kpoints; j++)
			{
				zmatrix[i][j] = 0.0;
			}
		}
		for(i = 0; i < kpoints; i++)
		{
			zmatrix[i][i] = 1.0;
		}

		// Calculte the points and weights
		TriQL(kpoints, a, b, zmatrix);

		for(i = 0; i < kpoints; i++)
		{
			z[i] = a[i];
			w[i] = b[i];
		}
		delete[] a;
		delete[] b;
		for (i = 0; i < kpoints; i++)
		{
                    delete[] zmatrix[i];
		}
		delete[] zmatrix;

	}

	/**
    \brief  Gauss-Radau-Kronrod-Jacobi zeros and weights.

    \li Generate \a npt=2*np Radau-Kronrod Jacobi zeros, \a z, and weights,\a w,
    associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)
    \f$,
    */
	void zwrk(double* z, double* w, const int npt ,const double alpha,
		      const double beta)
	{

		int np = npt/2;

		if(np < 2)
		{
			fprintf(stderr,"too few points in formula\n");
			return;
		}

		double end0 = -1;

		int i,j;

		// number of kronrod points associated with the np gauss rule
		int kpoints = 2*np;

		// Define the number of required recurrence coefficents
		int ncoeffs = (int)ceil(3.0*np/2);

		// Define arrays  for the recurrence coefficients
		double *a = new double[ncoeffs+1];
		double *b = new double[ncoeffs+1];

		// Initialize a and b to zero
		for(i = 0; i < ncoeffs+1; i++)
		{
			a[i] = 0.0;
			b[i] = 0.0;
		}

		// Call the routine to calculate the recurrence coefficients
		RecCoeff(ncoeffs,a,b,alpha,beta);

		double* a0 = new double[ncoeffs];
		double* b0 = new double[ncoeffs];

		chri1(ncoeffs,a,b,a0,b0,end0);

		double s = b0[0]/fabs(b0[0]);
		b0[0] = s*b0[0];

		// Finding the 2*np-1 gauss-kronrod points
		double* z1 = new double[2*np-1];
		double* w1 = new double[2*np-1];
		for(i = 0; i < ncoeffs; i++)
		{
			z1[i] = a0[i];
			w1[i] = b0[i];
		}
		JKMatrix(np-1,z1,w1);
		// Set up the identity matrix
		double** zmatrix = new double*[2*np-1];
		for(i = 0; i < 2*np-1; i++)
		{
			zmatrix[i] = new double[2*np-1];
			for(j = 0; j < 2*np-1; j++)
			{
				zmatrix[i][j] = 0.0;
			}
		}
		for(i = 0; i < 2*np-1; i++)
		{
			zmatrix[i][i] = 1.0;
		}

		// Calculate the points and weights
		TriQL(2*np-1, z1, w1, zmatrix);

		double sumW1 = 0.0;
		for(i = 0; i < 2*np-1; i++)
		{
			w1[i] = s*w1[i]/(z1[i]-end0);
			sumW1 += w1[i];
		}

		z[0] = end0;
		w[0] = b[0]- sumW1;
		for(i = 1; i < kpoints; i++)
		{
			z[i] = z1[i-1];
			w[i] = w1[i-1];
		}


		delete[] a;
		delete[] b;
		delete[] a0;
		delete[] b0;
		delete[] z1;
		delete[] w1;
                for(i = 0; i < 2*np-1; i++)
                {
                    delete[] zmatrix[i];
                }
		delete[] zmatrix;
	}

	/**
    \brief  Gauss-Lobatto-Kronrod-Jacobi zeros and weights.

    \li Generate \a npt=2*np-1 Lobatto-Kronrod Jacobi zeros, \a z, and weights,\a w,
    associated with the Jacobi polynomial \f$ P^{\alpha,\beta}_{np}(z)
    \f$,
    */
	void zwlk(double* z, double* w, const int npt,
		      const double alpha, const double beta)
	{

		int np = (npt+1)/2;

		if(np < 4)
		{
			fprintf (stderr,"too few points in formula\n");
			return;
		}

		double endl = -1;
		double endr = 1;
		int i,j;

		// number of kronrod points associated with the np gauss rule
		int kpoints = 2*np-1;

		// Define the number of required recurrence coefficents
		int ncoeffs = (int)ceil(3.0*np/2)-1;

		// Define arrays  for the recurrence coefficients
		double *a = new double[ncoeffs+1];
		double *b = new double[ncoeffs+1];

		// Initialize a and b to zero
		for(i = 0; i < ncoeffs+1; i++)
		{
			a[i] = 0.0;
			b[i] = 0.0;
		}

		// Call the routine to calculate the recurrence coefficients
		RecCoeff(ncoeffs,a,b,alpha,beta);


		double* a0 = new double[ncoeffs];
		double* b0 = new double[ncoeffs];

		chri1(ncoeffs,a,b,a0,b0,endl);

		double* a1 = new double[ncoeffs-1];
		double* b1 = new double[ncoeffs-1];

		chri1(ncoeffs-1,a0,b0,a1,b1,endr);


		double s = b1[0]/fabs(b1[0]);
		b1[0] = s*b1[0];

		// Finding the 2*np-1 gauss-kronrod points
		double* z1 = new double[2*np-3];
		double* w1 = new double[2*np-3];
		for(i = 0; i < ncoeffs; i++)
		{
			z1[i] = a1[i];
			w1[i] = b1[i];
		}
		JKMatrix(np-2,z1,w1);
		// Set up the identity matrix
		double** zmatrix = new double*[2*np-3];
		for(i = 0; i < 2*np-3; i++)
		{
			zmatrix[i] = new double[2*np-3];
			for(j = 0; j < 2*np-3; j++)
			{
				zmatrix[i][j] = 0.0;
			}
		}
		for(i = 0; i < 2*np-3; i++)
		{
			zmatrix[i][i] = 1.0;
		}

		// Calculate the points and weights
		TriQL(2*np-3, z1, w1, zmatrix);

		double sumW1 = 0.0;
		double sumW1Z1 = 0.0;
		for(i = 0; i < 2*np-3; i++)
		{
			w1[i] = s*w1[i]/(z1[i]-endl)/(z1[i]-endr);
			sumW1 += w1[i];
			sumW1Z1 += z1[i]*w1[i];
		}

		double c0 = b[0]-sumW1;
		double c1 = a[0]*b[0]-sumW1Z1;

		z[0] = endl;
		z[2*np-2] = endr;
		w[0] = (c0*endr-c1)/(endr-endl);
		w[2*np-2] = (c1-c0*endl)/(endr-endl);

		for(i = 1; i < kpoints-1; i++)
		{
			z[i] = z1[i-1];
			w[i] = w1[i-1];
		}
		delete[] a;
		delete[] b;
		delete[] a0;
		delete[] b0;
		delete[] a1;
		delete[] b1;
		delete[] z1;
		delete[] w1;
                for(i = 0; i < 2*np-3; i++)
		{
                    delete[] zmatrix[i];
                }
		delete[] zmatrix;
	}

    /**
    \brief Compute the Derivative Matrix and its transpose associated
    with the Gauss-Jacobi zeros.

    \li Compute the derivative matrix, \a d, and its transpose, \a dt,
    associated with the n_th order Lagrangian interpolants through the
    \a np Gauss-Jacobi points \a z such that \n
    \f$  \frac{du}{dz}(z[i]) =  \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$

    */

    void Dgj(double *D,  const double *z, const int np, const double alpha,
        const double beta)
    {

        double one = 1.0, two = 2.0;

        if (np <= 0){
            D[0] = 0.0;
        }
        else{
            int i,j;
            double *pd;

            pd = (double *)malloc(np*sizeof(double));
            jacobd(np,z,pd,np,alpha,beta);

            for (i = 0; i < np; i++){
                for (j = 0; j < np; j++){

                    if (i != j)
                        D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i]));
                    else
                        D[i*np+j] = (alpha - beta + (alpha + beta + two)*z[j])/
                        (two*(one - z[j]*z[j]));
                }
            }
            free(pd);
        }
        return;
    }


    /**
    \brief Compute the Derivative Matrix and its transpose associated
    with the Gauss-Radau-Jacobi zeros with a zero at \a z=-1.

    \li Compute the derivative matrix, \a d, associated with the n_th
    order Lagrangian interpolants through the \a np Gauss-Radau-Jacobi
    points \a z such that \n \f$ \frac{du}{dz}(z[i]) =
    \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$

    */

    void Dgrjm(double *D, const double *z, const int np, const double alpha,
        const double beta)
    {

        if (np <= 0){
            D[0] = 0.0;
        }
        else{
            int i, j;
            double   one = 1.0, two = 2.0;
            double   *pd;

            pd  = (double *)malloc(np*sizeof(double));

            pd[0] = pow(-one,np-1)*gammaF(np+beta+one);
            pd[0] /= gammaF(np)*gammaF(beta+two);
            jacobd(np-1,z+1,pd+1,np-1,alpha,beta+1);
            for(i = 1; i < np; ++i) pd[i] *= (1+z[i]);

            for (i = 0; i < np; i++) {
                for (j = 0; j < np; j++){
                    if (i != j)
                        D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i]));
                    else {
                        if(j == 0)
                            D[i*np+j] = -(np + alpha + beta + one)*(np - one)/
                            (two*(beta + two));
                        else
                            D[i*np+j] = (alpha - beta + one + (alpha + beta + one)*z[j])/
                            (two*(one - z[j]*z[j]));
                    }
                }
            }
            free(pd);
        }

        return;
    }


    /**
    \brief Compute the Derivative Matrix  associated with the
    Gauss-Radau-Jacobi zeros with a zero at \a z=1.

    \li Compute the derivative matrix, \a d, associated with the n_th
    order Lagrangian interpolants through the \a np Gauss-Radau-Jacobi
    points \a z such that \n \f$ \frac{du}{dz}(z[i]) =
    \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$
    */

    void Dgrjp(double *D, const double *z, const int np, const double alpha,
        const double beta)
    {

        if (np <= 0){
            D[0] = 0.0;
        }
        else{
            int i, j;
            double   one = 1.0, two = 2.0;
            double   *pd;

            pd  = (double *)malloc(np*sizeof(double));


            jacobd(np-1,z,pd,np-1,alpha+1,beta);
            for(i = 0; i < np-1; ++i) pd[i] *= (1-z[i]);
            pd[np-1] = -gammaF(np+alpha+one);
            pd[np-1] /= gammaF(np)*gammaF(alpha+two);

            for (i = 0; i < np; i++) {
                for (j = 0; j < np; j++){
                    if (i != j)
                        D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i]));
                    else {
                        if(j == np-1)
                            D[i*np+j] = (np + alpha + beta + one)*(np - one)/
                            (two*(alpha + two));
                        else
                            D[i*np+j] = (alpha - beta - one + (alpha + beta + one)*z[j])/
                            (two*(one - z[j]*z[j]));
                    }
                }
            }
            free(pd);
        }

        return;
    }

    /**
    \brief Compute the Derivative Matrix associated with the
    Gauss-Lobatto-Jacobi zeros.

    \li Compute the derivative matrix, \a d, associated with the n_th
    order Lagrange interpolants through the \a np
    Gauss-Lobatto-Jacobi points \a z such that \n \f$
    \frac{du}{dz}(z[i]) = \sum_{j=0}^{np-1} D[i*np+j] u(z[j]) \f$

    */

    void Dglj(double *D, const double *z, const int np, const double alpha,
        const double beta)
    {
        if (np <= 1){
            D[0] = 0.0;
        }
        else{
            int i, j;
            double   one = 1.0, two = 2.0;
            double   *pd;

            pd  = (double *)malloc(np*sizeof(double));

            pd[0]  = two*pow(-one,np)*gammaF(np + beta);
            pd[0] /= gammaF(np - one)*gammaF(beta + two);
            jacobd(np-2,z+1,pd+1,np-2,alpha+1,beta+1);
            for(i = 1; i < np-1; ++i) pd[i] *= (one-z[i]*z[i]);
            pd[np-1]  = -two*gammaF(np + alpha);
            pd[np-1] /= gammaF(np - one)*gammaF(alpha + two);

            for (i = 0; i < np; i++) {
                for (j = 0; j < np; j++){
                    if (i != j)
                        D[i*np+j] = pd[j]/(pd[i]*(z[j]-z[i]));
                    else {
                        if (j == 0)
                            D[i*np+j] = (alpha - (np-1)*(np + alpha + beta))/(two*(beta+ two));
                        else if (j == np-1)
                            D[i*np+j] =-(beta - (np-1)*(np + alpha + beta))/(two*(alpha+ two));
                        else
                            D[i*np+j] = (alpha - beta + (alpha + beta)*z[j])/
                            (two*(one - z[j]*z[j]));
                    }
                }
            }
            free(pd);
        }

        return;
    }


    /**
    \brief Compute the value of the \a i th Lagrangian interpolant through
    the \a np Gauss-Jacobi points \a zgj at the arbitrary location \a z.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
    %
    \f$ \begin{array}{rcl}
    h_j(z) =  \left\{ \begin{array}{ll}
    \displaystyle \frac{P_{np}^{\alpha,\beta}(z)}
    {[P_{np}^{\alpha,\beta}(z_j)]^\prime
    (z-z_j)} & \mbox{if $z \ne z_j$}\\
    & \\
    1 & \mbox{if $z=z_j$}
    \end{array}
    \right.
    \end{array}   \f$
    */

    double hgj (const int i, const double z, const double *zgj,
        const int np, const double alpha, const double beta)
    {
        boost::ignore_unused(alpha, beta);
	double zi, dz;

        zi  = *(zgj+i);
        dz  = z-zi;
	if (fabs(dz) < EPS) return 1.0;

	return laginterp(z, i, zgj, np);

    }

    /**
    \brief Compute the value of the \a i th Lagrangian interpolant through the
    \a np Gauss-Radau-Jacobi points \a zgrj at the arbitrary location
    \a z. This routine assumes \a zgrj includes the point \a -1.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
    %
    \f$ \begin{array}{rcl}
    h_j(z) = \left\{ \begin{array}{ll}
    \displaystyle \frac{(1+z) P_{np-1}^{\alpha,\beta+1}(z)}
    {((1+z_j) [P_{np-1}^{\alpha,\beta+1}(z_j)]^\prime +
    P_{np-1}^{\alpha,\beta+1}(z_j) ) (z-z_j)} & \mbox{if $z \ne z_j$}\\
    & \\
    1 & \mbox{if $z=z_j$}
    \end{array}
    \right.
    \end{array}   \f$
    */

    double hgrjm (const int i, const double z, const double *zgrj, const int np,
        const double alpha, const double beta)
    {
        boost::ignore_unused(alpha, beta);

	double zi, dz;

        zi  = *(zgrj+i);
        dz  = z-zi;
        if (fabs(dz) < EPS) return 1.0;

	return laginterp(z, i, zgrj, np);
    }


    /**
    \brief Compute the value of the \a i th Lagrangian interpolant through the
    \a np Gauss-Radau-Jacobi points \a zgrj at the arbitrary location
    \a z. This routine assumes \a zgrj includes the point \a +1.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
    %
    \f$ \begin{array}{rcl}
    h_j(z) = \left\{ \begin{array}{ll}
    \displaystyle \frac{(1-z) P_{np-1}^{\alpha+1,\beta}(z)}
    {((1-z_j) [P_{np-1}^{\alpha+1,\beta}(z_j)]^\prime -
    P_{np-1}^{\alpha+1,\beta}(z_j) ) (z-z_j)} & \mbox{if $z \ne z_j$}\\
    & \\
    1 & \mbox{if $z=z_j$}
    \end{array}
    \right.
    \end{array}   \f$
    */

    double hgrjp (const int i, const double z, const double *zgrj, const int np,
        const double alpha, const double beta)
    {
        boost::ignore_unused(alpha, beta);
	double zi, dz;

        zi  = *(zgrj+i);
        dz  = z-zi;
        if (fabs(dz) < EPS) return 1.0;

	return laginterp(z, i, zgrj, np);
    }


    /**
    \brief Compute the value of the \a i th Lagrangian interpolant through the
    \a np Gauss-Lobatto-Jacobi points \a zgrj at the arbitrary location
    \a z.

    \li \f$ -1 \leq z \leq 1 \f$

    \li Uses the defintion of the Lagrangian interpolant:\n
    %
    \f$ \begin{array}{rcl}
    h_j(z) = \left\{ \begin{array}{ll}
    \displaystyle \frac{(1-z^2) P_{np-2}^{\alpha+1,\beta+1}(z)}
    {((1-z^2_j) [P_{np-2}^{\alpha+1,\beta+1}(z_j)]^\prime -
    2 z_j P_{np-2}^{\alpha+1,\beta+1}(z_j) ) (z-z_j)}&\mbox{if $z \ne z_j$}\\
    & \\
    1 & \mbox{if $z=z_j$}
    \end{array}
    \right.
    \end{array}   \f$
    */

    double hglj (const int i, const double z, const double *zglj, const int np,
        const double alpha, const double beta)
    {
        boost::ignore_unused(alpha, beta);

        double zi, dz;

        zi  = *(zglj+i);
        dz  = z-zi;
        if (fabs(dz) < EPS) return 1.0;

	return laginterp(z, i, zglj, np);

    }


    /**
    \brief Interpolation Operator from Gauss-Jacobi points to an
    arbitrary distribution at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Jacobi distribution of \a nz
    zeros \a zgrj to an arbitrary distribution of \a mz points \a zm, i.e.\n
    \f$
    u(zm[i]) = \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j])
    \f$

    */

    void Imgj(double *im, const double *zgj, const double *zm, const int nz,
        const int mz,const double alpha, const double beta){
            double zp;
            int i, j;

            for (i = 0; i < nz; ++i) {
                for (j = 0; j < mz; ++j)
                {
                    zp = zm[j];
                    im [i*mz+j] = hgj(i, zp, zgj, nz, alpha, beta);
                }
            }

            return;
    }

    /**
    \brief Interpolation Operator from Gauss-Radau-Jacobi points
    (including \a z=-1) to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Radau-Jacobi distribution of
    \a nz zeros \a zgrj (where \a zgrj[0]=-1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j]) \f$

    */

    void Imgrjm(double *im, const double *zgrj, const double *zm, const int nz,
        const int mz, const double alpha, const double beta)
    {
        double zp;
        int i, j;

        for (i = 0; i < nz; i++) {
            for (j = 0; j < mz; j++)
            {
                zp = zm[j];
                im [i*mz+j] = hgrjm(i, zp, zgrj, nz, alpha, beta);
            }
        }

        return;
    }

    /**
    \brief Interpolation Operator from Gauss-Radau-Jacobi points
    (including \a z=1) to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Radau-Jacobi distribution of
    \a nz zeros \a zgrj (where \a zgrj[nz-1]=1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j]) \f$

    */

    void Imgrjp(double *im, const double *zgrj, const double *zm, const int nz,
        const int mz,const double alpha, const double beta)
    {
            double zp;
            int i, j;

            for (i = 0; i < nz; i++) {
                for (j = 0; j < mz; j++)
                {
                    zp = zm[j];
                    im [i*mz+j] = hgrjp(i, zp, zgrj, nz, alpha, beta);
                }
            }

            return;
    }


    /**
    \brief Interpolation Operator from Gauss-Lobatto-Jacobi points
    to an arbitrary distrubtion at points \a zm

    \li Computes the one-dimensional interpolation matrix, \a im, to
    interpolate a function from at Gauss-Lobatto-Jacobi distribution of
    \a nz zeros \a zgrj (where \a zgrj[0]=-1) to an arbitrary
    distribution of \a mz points \a zm, i.e.
    \n
    \f$ u(zm[i]) =    \sum_{j=0}^{nz-1} im[i*nz+j] \ u(zgj[j]) \f$

    */

    void Imglj(double *im, const double *zglj, const double *zm, const int nz,
        const int mz, const double alpha, const double beta)
    {
        double zp;
        int i, j;

        for (i = 0; i < nz; i++) {
            for (j = 0; j < mz; j++)
            {
                zp = zm[j];
                im[i*mz+j] = hglj(i, zp, zglj, nz, alpha, beta);
            }
        }

        return;
    }

    /**
    \brief Routine to calculate Jacobi polynomials, \f$
    P^{\alpha,\beta}_n(z) \f$, and their first derivative, \f$
    \frac{d}{dz} P^{\alpha,\beta}_n(z) \f$.

    \li This function returns the vectors \a poly_in and \a poly_d
    containing the value of the \f$ n^th \f$ order Jacobi polynomial
    \f$ P^{\alpha,\beta}_n(z) \alpha > -1, \beta > -1 \f$ and its
    derivative at the \a np points in \a z[i]

    - If \a poly_in = NULL then only calculate derivatice

    - If \a polyd   = NULL then only calculate polynomial

    - To calculate the polynomial this routine uses the recursion
    relationship (see appendix A ref [4]) :
    \f$ \begin{array}{rcl}
    P^{\alpha,\beta}_0(z) &=& 1 \\
    P^{\alpha,\beta}_1(z) &=& \frac{1}{2} [ \alpha-\beta+(\alpha+\beta+2)z] \\
    a^1_n P^{\alpha,\beta}_{n+1}(z) &=& (a^2_n + a^3_n z)
    P^{\alpha,\beta}_n(z) - a^4_n P^{\alpha,\beta}_{n-1}(z) \\
    a^1_n &=& 2(n+1)(n+\alpha + \beta + 1)(2n + \alpha + \beta) \\
    a^2_n &=& (2n + \alpha + \beta + 1)(\alpha^2 - \beta^2)  \\
    a^3_n &=& (2n + \alpha + \beta)(2n + \alpha + \beta + 1)
    (2n + \alpha + \beta + 2)  \\
    a^4_n &=& 2(n+\alpha)(n+\beta)(2n + \alpha + \beta + 2)
    \end{array} \f$

    - To calculate the derivative of the polynomial this routine uses
    the relationship (see appendix A ref [4]) :
    \f$ \begin{array}{rcl}
    b^1_n(z)\frac{d}{dz} P^{\alpha,\beta}_n(z)&=&b^2_n(z)P^{\alpha,\beta}_n(z)
    + b^3_n(z) P^{\alpha,\beta}_{n-1}(z) \hspace{2.2cm} \\
    b^1_n(z) &=& (2n+\alpha + \beta)(1-z^2) \\
    b^2_n(z) &=& n[\alpha - \beta - (2n+\alpha + \beta)z]\\
    b^3_n(z) &=& 2(n+\alpha)(n+\beta)
    \end{array} \f$

    - Note the derivative from this routine is only valid for -1 < \a z < 1.
    */
    void jacobfd(const int np, const double *z, double *poly_in, double *polyd,
        const int n, const double alpha, const double beta){
            int i;
            double  zero = 0.0, one = 1.0, two = 2.0;

            if(!np)
                return;

            if(n == 0){
                if(poly_in)
                    for(i = 0; i < np; ++i)
                        poly_in[i] = one;
                if(polyd)
                    for(i = 0; i < np; ++i)
                        polyd[i] = zero;
            }
            else if (n == 1){
                if(poly_in)
                    for(i = 0; i < np; ++i)
                        poly_in[i] = 0.5*(alpha - beta + (alpha + beta + two)*z[i]);
                if(polyd)
                    for(i = 0; i < np; ++i)
                        polyd[i] = 0.5*(alpha + beta + two);
            }
            else{
                int k;
                double   a1,a2,a3,a4;
                double   two = 2.0, apb = alpha + beta;
                double   *poly, *polyn1,*polyn2;

                if(poly_in){ // switch for case of no poynomial function return
                    polyn1 = (double *)malloc(2*np*sizeof(double));
                    polyn2 = polyn1+np;
                    poly   = poly_in;
                }
                else{
                    polyn1 = (double *)malloc(3*np*sizeof(double));
                    polyn2 = polyn1+np;
                    poly   = polyn2+np;
                }

                for(i = 0; i < np; ++i){
                    polyn2[i] = one;
                    polyn1[i] = 0.5*(alpha - beta + (alpha + beta + two)*z[i]);
                }

                for(k = 2; k <= n; ++k){
                    a1 =  two*k*(k + apb)*(two*k + apb - two);
                    a2 = (two*k + apb - one)*(alpha*alpha - beta*beta);
                    a3 = (two*k + apb - two)*(two*k + apb - one)*(two*k + apb);
                    a4 =  two*(k + alpha - one)*(k + beta - one)*(two*k + apb);

                    a2 /= a1;
                    a3 /= a1;
                    a4 /= a1;

                    for(i = 0; i < np; ++i){
                        poly  [i] = (a2 + a3*z[i])*polyn1[i] - a4*polyn2[i];
                        polyn2[i] = polyn1[i];
                        polyn1[i] = poly  [i];
                    }
                }

                if(polyd){
                    a1 = n*(alpha - beta);
                    a2 = n*(two*n + alpha + beta);
                    a3 = two*(n + alpha)*(n + beta);
                    a4 = (two*n + alpha + beta);
                    a1 /= a4;  a2 /= a4;   a3 /= a4;

                    // note polyn2 points to polyn1 at end of poly iterations
                    for(i = 0; i < np; ++i){
                        polyd[i]  = (a1- a2*z[i])*poly[i] + a3*polyn2[i];
                        polyd[i] /= (one - z[i]*z[i]);
                    }
                }

                free(polyn1);
            }

            return;
    }


    /**
    \brief Calculate the  derivative of Jacobi polynomials

    \li Generates a vector \a poly of values of the derivative of the
    \a n th order Jacobi polynomial \f$ P^(\alpha,\beta)_n(z)\f$ at the
    \a np points \a z.

    \li To do this we have used the relation
    \n
    \f$ \frac{d}{dz} P^{\alpha,\beta}_n(z)
    = \frac{1}{2} (\alpha + \beta + n + 1)  P^{\alpha,\beta}_n(z) \f$

    \li This formulation is valid for \f$ -1 \leq z \leq 1 \f$

    */

    void jacobd(const int np, const double *z, double *polyd, const int n,
        const double alpha, const double beta)
    {
        int i;
        double one = 1.0;
        if(n == 0)
            for(i = 0; i < np; ++i) polyd[i] = 0.0;
        else{
            //jacobf(np,z,polyd,n-1,alpha+one,beta+one);
            jacobfd(np,z,polyd,NULL,n-1,alpha+one,beta+one);
            for(i = 0; i < np; ++i) polyd[i] *= 0.5*(alpha + beta + (double)n + one);
        }
        return;
    }


    /**
    \brief Calculate the Gamma function , \f$ \Gamma(n)\f$, for integer
    values and halves.

    Determine the value of \f$\Gamma(n)\f$ using:

    \f$ \Gamma(n) = (n-1)!  \mbox{ or  }  \Gamma(n+1/2) = (n-1/2)\Gamma(n-1/2)\f$

    where \f$ \Gamma(1/2) = \sqrt(\pi)\f$
    */

    double gammaF(const double x){
        double gamma = 1.0;

        if     (x == -0.5) gamma = -2.0*sqrt(M_PI);
        else if (!x) return gamma;
        else if ((x-(int)x) == 0.5){
            int n = (int) x;
            double tmp = x;

            gamma = sqrt(M_PI);
            while(n--){
                tmp   -= 1.0;
                gamma *= tmp;
            }
        }
        else if ((x-(int)x) == 0.0){
            int n = (int) x;
            double tmp = x;

            while(--n){
                tmp   -= 1.0;
                gamma *= tmp;
            }
        }
        else
            fprintf(stderr,"%lf is not of integer or half order\n",x);
        return gamma;
    }

    /**
    \brief  Calculate the \a n zeros, \a z, of the Jacobi polynomial, i.e.
    \f$ P_n^{\alpha,\beta}(z) = 0 \f$

    This routine is only value for \f$( \alpha > -1, \beta > -1)\f$
    and uses polynomial deflation in a Newton iteration
    */

    static void Jacobz(const int n, double *z, const double alpha,
        const double beta){
            int i,j,k;
            double   dth = M_PI/(2.0*(double)n);
            double   poly,pder,rlast=0.0;
            double   sum,delr,r;
            double one = 1.0, two = 2.0;

            if(!n)
                return;

            for(k = 0; k < n; ++k){
                r = -cos((two*(double)k + one) * dth);
                if(k) r = 0.5*(r + rlast);

                for(j = 1; j < STOP; ++j){
                    jacobfd(1,&r,&poly, &pder, n, alpha, beta);

                    for(i = 0, sum = 0.0; i < k; ++i) sum += one/(r - z[i]);

                    delr = -poly / (pder - sum * poly);
                    r   += delr;
                    if( fabs(delr) < EPS ) break;
                }
                z[k]  = r;
                rlast = r;
            }
            return;
    }


    /**
    \brief Zero and Weight determination through the eigenvalues and eigenvectors of a tridiagonal
    matrix from the three term recurrence relationship.

    Set up a symmetric tridiagonal matrix

    \f$ \left [  \begin{array}{ccccc}
    a[0] & b[0]   &        &        & \\
    b[0] & a[1]   & b[1]   &        & \\
    0   & \ddots & \ddots & \ddots &  \\
    &        & \ddots & \ddots & b[n-2] \\
    &        &        & b[n-2] & a[n-1] \end{array} \right ] \f$

    Where the coefficients a[n], b[n] come from the  recurrence relation

    \f$  b_j p_j(z) = (z - a_j ) p_{j-1}(z) - b_{j-1}   p_{j-2}(z) \f$

    where \f$ j=n+1\f$ and \f$p_j(z)\f$ are the Jacobi (normalized)
    orthogonal polynomials \f$ \alpha,\beta > -1\f$( integer values and
    halves). Since the polynomials are orthonormalized, the tridiagonal
    matrix is guaranteed to be symmetric. The eigenvalues of this
    matrix are the zeros of the Jacobi polynomial.
    */

    void JacZeros(const int n, double *a, double*b, const double alpha,
        const double beta){

			int i,j;
			RecCoeff(n,a,b,alpha,beta);

			double **z = new double*[n];
			for(i = 0; i < n; i++)
			{
				z[i] = new double[n];
				for(j = 0; j < n; j++)
				{
					z[i][j] = 0.0;
				}
			}
			for(i = 0; i < n; i++)
			{
				z[i][i] = 1.0;
			}

			// find eigenvalues and eigenvectors
            TriQL(n, a, b,z);

            delete[] z;
            return;
    }

	/**
    \brief  The routine finds the recurrence coefficients \a a and
	\a b of the orthogonal polynomials
	*/
	static void RecCoeff(const int n, double *a, double *b,const double alpha,
              const double beta){

        int i;
        double apb, apbi,a2b2;

        if(!n)
            return;

        // generate normalised terms
        apb  = alpha + beta;
        apbi = 2.0 + apb;

        b[0] = pow(2.0,apb+1.0)*gammaF(alpha+1.0)*gammaF(beta+1.0)/gammaF(apbi); //MuZero
		a[0]   = (beta-alpha)/apbi;
		b[1]   = (4.0*(1.0+alpha)*(1.0+beta)/((apbi+1.0)*apbi*apbi));

		a2b2 = beta*beta-alpha*alpha;

		for(i = 1; i < n-1; i++){
            apbi = 2.0*(i+1) + apb;
            a[i] = a2b2/((apbi-2.0)*apbi);
            b[i+1] = (4.0*(i+1)*(i+1+alpha)*(i+1+beta)*(i+1+apb)/
                 ((apbi*apbi-1)*apbi*apbi));
		}

		apbi   = 2.0*n + apb;
        a[n-1] = a2b2/((apbi-2.0)*apbi);

	}


    /** \brief QL algorithm for symmetric tridiagonal matrix

    This subroutine is a translation of an algol procedure,
    num. math. \b 12, 377-383(1968) by martin and wilkinson, as modified
    in num. math. \b 15, 450(1970) by dubrulle.  Handbook for
    auto. comp., vol.ii-linear algebra, 241-248(1971).  This is a
    modified version from numerical recipes.

    This subroutine finds the eigenvalues and first components of the
    eigenvectors of a symmetric tridiagonal matrix by the implicit QL
    method.

    on input:
    - n is the order of the matrix;
    - d contains the diagonal elements of the input matrix;
    - e contains the subdiagonal elements of the input matrix
    in its first n-2 positions.
	- z is the n by n identity matrix

    on output:

    - d contains the eigenvalues in ascending order.
    - e contains the weight values - modifications of the first component
	  of normalised eigenvectors
    */

    static void TriQL(const int n, double *d,double *e, double **z){
        int m,l,iter,i,k;
        double s,r,p,g,f,dd,c,b;

		double MuZero = e[0];

		// Renumber the elements of e
		for(i = 0; i < n-1; i++)
		{
			e[i] = sqrt(e[i+1]);
		}
		e[n-1] = 0.0;


        for (l=0;l<n;l++) {
            iter=0;
            do {
                for (m=l;m<n-1;m++) {
                    dd=fabs(d[m])+fabs(d[m+1]);
                    if (fabs(e[m])+dd == dd) break;
                }
                if (m != l) {
                    if (iter++ == STOP){
                        fprintf(stderr,"triQL: Too many iterations in TQLI");
                        exit(1);
                    }
                    g=(d[l+1]-d[l])/(2.0*e[l]);
                    r=sqrt((g*g)+1.0);
                    g=d[m]-d[l]+e[l]/(g+sign(r,g));
                    s=c=1.0;
                    p=0.0;
                    for (i=m-1;i>=l;i--) {
                        f=s*e[i];
                        b=c*e[i];
                        if (fabs(f) >= fabs(g)) {
                            c=g/f;
                            r=sqrt((c*c)+1.0);
                            e[i+1]=f*r;
                            c *= (s=1.0/r);
                        } else {
                            s=f/g;
                            r=sqrt((s*s)+1.0);
                            e[i+1]=g*r;
                            s *= (c=1.0/r);
                        }
                        g=d[i+1]-p;
                        r=(d[i]-g)*s+2.0*c*b;
                        p=s*r;
                        d[i+1]=g+p;
                        g=c*r-b;

						// Calculate the eigenvectors
						for(k = 0; k < n; k++)
						{
							f = z[k][i+1];
							z[k][i+1] = s*z[k][i] + c*f;
							z[k][i] = c*z[k][i] - s*f;
						}

                    }
                    d[l]=d[l]-p;
                    e[l]=g;
                    e[m]=0.0;
                }
            } while (m != l);
        }

        // order eigenvalues
		// Since we only need the first component of the eigenvectors
		// to calcualte the weight, we only swap the first components
        for(i = 0; i < n-1; ++i){
            k = i;
            p = d[i];
            for(l = i+1; l < n; ++l)
                if (d[l] < p) {
                    k = l;
                    p = d[l];
                }
            d[k] = d[i];
            d[i] = p;

			double temp = z[0][k];
			z[0][k] = z[0][i];
			z[0][i] = temp;
        }

		// Calculate the weights
		for(i =0 ; i < n; i++)
		{
			e[i] = MuZero*z[0][i]*z[0][i];
		}
    }

	/**
    \brief Calcualtes the Jacobi-kronrod matrix by determining the
	\a a and \b coefficients.

	The first \a 3n+1 coefficients are already known

	For more information refer to:
	"Dirk P. Laurie, Calcualtion of Gauss-Kronrod quadrature rules"
	*/
	void JKMatrix(int n, double *a, double *b)
	{
		int i,j,k,m;
		// Working storage
		int size = (int)floor(n/2.0)+2;
		double *s = new double[size];
		double *t = new double[size];

		// Initialize s and t to zero
		for(i = 0; i < size; i++)
		{
			s[i] = 0.0;
			t[i] = 0.0;
		}

		t[1] = b[n+1];
		for(m = 0; m <= n-2; m++)
		{
			 double u = 0.0;
			 for(k = (int)floor((m+1)/2.0); k >= 0; k--)
			 {
				int l = m-k;
				u = u+(a[k+n+1]-a[l])*t[k+1] + b[k+n+1]*s[k] - b[l]*s[k+1];
				s[k+1] = u;
			 }

			 // Swap the contents of s and t
			 double *hold = s;
			 s = t;
			 t  = hold;
		}


		for(j = (int)floor(n/2.0); j >= 0; j--)
		{
			s[j+1] = s[j];
		}

		for(m = n-1; m <= 2*n-3; m++)
		{
			double u = 0;
			for(k = m+1-n; k <= floor((m-1)/2.0); k++)
			{
				int l = m-k;
				j = n-1-l;
				u = u-(a[k+n+1]-a[l])*t[j+1] - b[k+n+1]*s[j+1] + b[l]*s[j+2];
				s[j+1] = u;
			}

			if(m%2 == 0)
			{
				k = m/2;
				a[k+n+1] = a[k] + (s[j+1]-b[k+n+1]*s[j+2])/t[j+2];

			}else
			{
				k = (m+1)/2;
				b[k+n+1] = s[j+1]/s[j+2];
			}


			// Swap the contents of s and t
			double  *hold = s;
			s = t;
			t  = hold;
		}

		a[2*n ] = a[n-1]-b[2*n]*s[1]/t[1];

	}

	/**
	\brief
	Given a weight function \f$w(t)\f$ through the first \a n+1
	coefficients \a a and \a b of its orthogonal polynomials
	this routine generates the first \a n recurrence coefficients for the orthogonal
	polynomials relative to the modified weight function \f$(t-z)w(t)\f$.

	The result will be placed in the array \a a0 and \a b0.
	*/

	void chri1(int n, double* a, double* b, double* a0,
		   double* b0,double z)
	{

		double q = ceil(3.0*n/2);
		int size = (int)q+1;
		if(size < n+1)
		{
			fprintf(stderr,"input arrays a and b are too short\n");
		}
		double* r = new double[n+1];
		r[0] = z - a[0];
		r[1] = z - a[1] - b[1]/r[0];
		a0[0] = a[1] + r[1] - r[0];
		b0[0] = -r[0]*b[0];

		if(n == 1)
		{
			delete[] r;
			return;
		}
		int k = 0;
		for(k = 1; k < n; k++)
		{
			r[k+1] = z - a[k+1] - b[k+1]/r[k];
			a0[k] = a[k+1] + r[k+1] - r[k];
			b0[k] = b[k] * r[k]/r[k-1];
		}
                delete[] r;

	}

    /**

	\brief

    Calcualte the bessel function of the first kind with complex double input y.
    Taken from Numerical Recipies in C

    Returns a complex double
    */


    std::complex<Nektar::NekDouble> ImagBesselComp(int n,std::complex<Nektar::NekDouble> y)
    {
    	std::complex<Nektar::NekDouble> z (1.0,0.0);
    	std::complex<Nektar::NekDouble> zbes (1.0,0.0);
    	std::complex<Nektar::NekDouble> zarg;
        Nektar::NekDouble tol = 1e-15;
    	int maxit = 10000;
    	int i = 1;

	    zarg = -0.25*y*y;

	    while (abs(z) > tol && i <= maxit){
		    z = z*(1.0/i/(i+n)*zarg);
	    	if  (abs(z) <= tol) break;
	    	zbes = zbes + z;
	    	i++;
    	}
        zarg = 0.5*y;
        for (i=1;i<=n;i++){
            zbes = zbes*zarg;
        }
        return zbes;

    }
} // end of namespace

