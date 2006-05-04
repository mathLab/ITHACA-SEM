///////////////////////////////////////////////////////////////////////////////
//
// File Blas.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: wrapper of functions around standard BLAS routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef BLAS_HPP
#define BLAS_HPP

#include <LibUtilities/TransF77.hpp>

/** Translations for using Fortran version of blas */
namespace Blas
{

  extern "C"
  {
    // -- BLAS Level 1:
    void   F77NAME(dcopy) (const int& n, const double *x, const int& incx,
               double *y, const int& incy);
    void   F77NAME(daxpy) (const int& n, const double& alpha, const double *x,
               const int& incx, const double *y, const int& incy);
    void   F77NAME(dswap) (const int& n, double *x, const int& incx,
               double *y, const int& incy);
    void   F77NAME(dscal) (const int& n, const double& alpha, double *x,
               const int& incx);
    void   F77NAME(drot)  (const int& n, double *x, const int& incx,
               double *y, const int& incy, const double& c,
               const double& s);
    double F77NAME(ddot)  (const int& n, const double *x,  const int& incx,
               const double *y, const int& incy);
    double F77NAME(dnrm2) (const int& n, const double *x, const int& incx);
    double F77NAME(dasum) (const int& n, const double *x, const int& incx);
    int    F77NAME(idamax)(const int& n, const double *x, const int& incx);

    // -- BLAS level 2
    void F77NAME(dgemv) (const char& trans,  const int& m,
             const int& n,       const double& alpha,
             const double* a,    const int& lda,
             const double* x,    const int& incx,
             const double& beta, double* y, const int& incy);

    void F77NAME(dspmv) (const char& trans, const int& n,    const double& alpha,
             const double* a,   const double* x, const int& incx,
             const double& beta,      double* y, const int& incy);

    // -- BLAS level 3:
    void F77NAME(dgemm) (const char& trans,   const char& transb,
             const int& m1,       const int& n,
             const int& k,        const double& alpha,
             const double* a,     const int& lda,
             const double* b,     const int& ldb,
             const double& beta,  double* c, const int& ldc);
  }

  /// \brief BLAS level 1: Copy \a x to \a y
  static void Dcopy (const int& n, const double *x, const int& incx,
             double *y, const int& incy)
  {
    F77NAME(dcopy)(n,x,incx,y,incy);
  }

  /// \brief  BLAS level 1: y = alpha \a x plus \a y
  static void Daxpy (const int& n, const double& alpha, const double *x,
             const int& incx,  const double *y, const int& incy)
  {
    F77NAME(daxpy)(n,alpha,x,incx,y,incy);
  }

  /// \brief BLAS level 1: Swap \a x with  \a y
  static void Dswap (const int& n,double *x, const int& incx,
             double *y, const int& incy)
  {
    F77NAME(dswap)(n,x,incx,y,incy);
  }

  /// \brief  BLAS level 1: x = alpha \a x
  static void Dscal (const int& n, const double& alpha, double *x,
             const int& incx)
  {
    F77NAME(dscal)(n,alpha,x,incx);
  }

  /// \brief BLAS level 1: Plane rotation by c = cos(theta), s = sin(theta)
  static void Drot (const int& n,  double *x,  const int& incx,
            double *y, const int& incy, const double& c,
            const double& s)
  {
    F77NAME(drot)(n,x,incx,y,incy,c,s);
  }

  /// \brief BLAS level 1: output =   \f$ x^T  y \f$
  static double Ddot (const int& n, const double *x, const int& incx,
              const double *y, const int& incy)
  {
    return F77NAME(ddot)(n,x,incx,y,incy);
  }

  // \brief  BLAS level 1: output = \f$ ||x||_2 \f$

  static double Dnrm2 (const int& n, const double *x, const int& incx)
  {
    return F77NAME(dnrm2)(n,x,incx);
  }

  /// \brief  BLAS level 1: output = \f$ ||x||_1 \f$
  static double Dasum (const int& n, const double *x, const int& incx)
  {
    return F77NAME(dasum)(n,x,incx);
  }

  /// \brief BLAS level 1: output = 1st value where \f$ |x[i]| = max |x|_1 \f$
  /// Note it is modified to return a value between (0,n-1) as per
  /// the standard C convention
  static int Idamax (const int& n, const double *x,  const int& incx)
  {
    return F77NAME(idamax)(n,x,incx) -1;
  }

  /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
  static void Dgemv (const char& trans,   const int& m,    const int& n,
             const double& alpha, const double* a, const int& lda,
             const double* x,     const int& incx, const double& beta,
             double* y,     const int& incy)
  {
    F77NAME(dgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
  }


  /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A
  /// is symmetric packed
  static void Dspmv (const char& trans,  const int& n,    const double& alpha,
             const double* a,    const double* x, const int& incx,
             const double& beta,       double* y, const int& incy)
  {
    F77NAME(dspmv) (trans,n,alpha,a,x,incx,beta,y,incy);
  }


  /// \brief BLAS level 3: Matrix-matrix multiply C = A x B where A[m x n],
  ///   B[n x k], C[m x k]
  static void Dgemm (const char& transa,  const char& transb, const int& m,
          const int& n,        const int& k,       const double& alpha,
          const double* a,     const int& lda,     const double* b,
          const int& ldb,      const double& beta,       double* c,
          const int& ldc)
  {
    F77NAME(dgemm) (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
  }

  // \brief Wrapper to mutliply two (row major) matrices together C =
  // a*A*B + b*C
  static void Cdgemm(const int M, const int N, const int K, const double a,
          double *A, const int ldA, double * B, const int ldB,
          const double b, double *C, const int ldC)
  {
    Dgemm('N','N',N,M,K,a,B,ldB,A,ldA,b,C,ldC) ;
  }
}
#endif //BLAS_HPP

/***
$Log: Blas.hpp,v $
Revision 1.4  2006/02/26 21:13:45  bnelson
Fixed a variety of compiler errors caused by updates to the coding standard.

Revision 1.3  2006/02/15 08:07:15  sherwin

Put codes into standard although have not yet been compiled

Revision 1.2  2006/02/12 21:51:42  sherwin

Added licence

Revision 1.1  2006/02/12 15:06:12  sherwin

Changed .h files to .hpp

**/
