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

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLAS_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLAS_HPP

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/LinearAlgebra/TransF77.hpp>
#include <LibUtilities/LinearAlgebra/BlasArray.hpp>

// Translations for using Fortran version of blas
namespace Blas
{
    extern "C"
    {
        // -- BLAS Level 1:
        void   F77NAME(dcopy) (const int& n, const double *x, const int& incx,
                   double *y, const int& incy);
        void   F77NAME(scopy) (const int& n, const float *x, const int& incx,
                   float *y, const int& incy);
        void   F77NAME(daxpy) (const int& n, const double& alpha, const double *x,
                   const int& incx, const double *y, const int& incy);
        void   F77NAME(saxpy) (const int& n, const float& alpha, const float *x,
                   const int& incx, const float *y, const int& incy);
        void   F77NAME(dswap) (const int& n, double *x, const int& incx,
                   double *y, const int& incy);
        void   F77NAME(sswap) (const int& n, float *x, const int& incx,
                   float *y, const int& incy);
        void   F77NAME(dscal) (const int& n, const double& alpha, double *x,
                   const int& incx);
        void   F77NAME(sscal) (const int& n, const float& alpha, float *x,
                   const int& incx);
        void   F77NAME(drot)  (const int& n, double *x, const int& incx,
                   double *y, const int& incy, const double& c,
                   const double& s);
        void   F77NAME(srot)  (const int& n, float *x, const int& incx,
                   float *y, const int& incy, const float& c,
                   const float& s);
        double F77NAME(ddot)  (const int& n, const double *x,  const int& incx,
                   const double *y, const int& incy);
        float F77NAME(sdot)  (const int& n, const float *x,  const int& incx,
                   const float *y, const int& incy);
        double F77NAME(dnrm2) (const int& n, const double *x, const int& incx);
        float F77NAME(snrm2) (const int& n, const float *x, const int& incx);
        double F77NAME(dasum) (const int& n, const double *x, const int& incx);
        float F77NAME(sasum) (const int& n, const float *x, const int& incx);
        int    F77NAME(idamax)(const int& n, const double *x, const int& incx);
        int    F77NAME(isamax)(const int& n, const float *x, const int& incx);

        // -- BLAS level 2
        void F77NAME(dgemv) (const char& trans,  const int& m,
                 const int& n,       const double& alpha,
                 const double* a,    const int& lda,
                 const double* x,    const int& incx,
                 const double& beta, double* y, const int& incy);
        // -- BLAS level 2
        void F77NAME(sgemv) (const char& trans,  const int& m,
                 const int& n,       const float& alpha,
                 const float* a,    const int& lda,
                 const float* x,    const int& incx,
                 const float& beta, float* y, const int& incy);

        void F77NAME(dgbmv) (const char& trans,  const int& m,
                 const int& n, const int& kl, const int& ku,
                 const double& alpha,
                 const double* a,    const int& lda,
                 const double* x,    const int& incx,
                 const double& beta, double* y, const int& incy);
        void F77NAME(sgbmv) (const char& trans,  const int& m,
                 const int& n, const int& kl, const int& ku,
                 const float& alpha,
                 const float* a,    const int& lda,
                 const float* x,    const int& incx,
                 const float& beta, float* y, const int& incy);

        void F77NAME(dtpmv) (const char& uplo, const char& trans, const char& diag,
                 const int& n, const double* ap, double* x, const int& incx);
        void F77NAME(stpmv) (const char& uplo, const char& trans, 
            const char& diag, const int& n, const float* ap, float* x, 
            const int& incx);

        void F77NAME(dspmv) (const char& uplo, const int& n,    const double& alpha,
                 const double* a,   const double* x, const int& incx,
                 const double& beta,      double* y, const int& incy);
        void F77NAME(sspmv) (const char& uplo, const int& n, const float& alpha,
                 const float* a,   const float* x, const int& incx,
                 const float& beta,      float* y, const int& incy);

        void F77NAME(dsbmv) (const char& uplo,  const int& m,
                 const int& k,       const double& alpha,
                 const double* a,    const int& lda,
                 const double* x,    const int& incx,
                 const double& beta, double* y, const int& incy);
        void F77NAME(ssbmv) (const char& uplo,  const int& m,
                 const int& k,       const float& alpha,
                 const float* a,    const int& lda,
                 const float* x,    const int& incx,
                 const float& beta, float* y, const int& incy);
        
        void F77NAME(dger) (const int& m,
                 const int& n,       const double& alpha,
                 const double* x,    const int& incx,
                 const double* y,    const int& incy,
                 double* a,    const int& lda);
        void F77NAME(sger) (const int& m,
                 const int& n,       const float& alpha,
                 const float* x,    const int& incx,
                 const float* y,    const int& incy,
                 float* a,    const int& lda);

        // -- BLAS level 3:
        void F77NAME(dgemm) (const char& trans,   const char& transb,
                 const int& m1,       const int& n,
                 const int& k,        const double& alpha,
                 const double* a,     const int& lda,
                 const double* b,     const int& ldb,
                 const double& beta,  double* c, const int& ldc);
        // -- BLAS level 3:
        void F77NAME(sgemm) (const char& trans,   const char& transb,
                 const int& m1,       const int& n,
                 const int& k,        const float& alpha,
                 const float* a,     const int& lda,
                 const float* b,     const int& ldb,
                 const float& beta,  float* c, const int& ldc);
    }

    /// \brief BLAS level 1: Copy \a x to \a y
    static inline void Dcopy (const int& n, const double *x, const int& incx,
             double *y, const int& incy)
    {
        F77NAME(dcopy)(n,x,incx,y,incy);
    }

    /// \brief  BLAS level 1: y = alpha \a x plus \a y
    static inline void Daxpy (const int& n, const double& alpha, const double *x,
             const int& incx,  const double *y, const int& incy)
    {
        F77NAME(daxpy)(n,alpha,x,incx,y,incy);
    }


    /// \brief BLAS level 1: Swap \a x with  \a y
    static inline void Dswap (const int& n,double *x, const int& incx,
             double *y, const int& incy)
    {
        F77NAME(dswap)(n,x,incx,y,incy);
    }

    /// \brief  BLAS level 1: x = alpha \a x
    static inline void Dscal (const int& n, const double& alpha, double *x,
             const int& incx)
    {
        F77NAME(dscal)(n,alpha,x,incx);
    }

    /// \brief BLAS level 1: Plane rotation by c = cos(theta), s = sin(theta)
    static inline void Drot (const int& n,  double *x,  const int& incx,
            double *y, const int& incy, const double& c,
            const double& s)
    {
        F77NAME(drot)(n,x,incx,y,incy,c,s);
    }

    /// \brief BLAS level 1: output =   \f$ x^T  y \f$
    static inline double Ddot (const int& n, const double *x, const int& incx,
              const double *y, const int& incy)
    {
        return F77NAME(ddot)(n,x,incx,y,incy);
    }

    // \brief  BLAS level 1: output = \f$ ||x||_2 \f$

    static inline double Dnrm2 (const int& n, const double *x, const int& incx)
    {
        return F77NAME(dnrm2)(n,x,incx);
    }

    /// \brief  BLAS level 1: output = \f$ ||x||_1 \f$
    static inline double Dasum (const int& n, const double *x, const int& incx)
    {
        return F77NAME(dasum)(n,x,incx);
    }

    /// \brief BLAS level 1: output = 1st value where \f$ |x[i]| = max |x|_1 \f$
    /// Note it is modified to return a value between (0,n-1) as per
    /// the standard C convention
    static inline int Idamax (const int& n, const double *x,  const int& incx)
    {
        return F77NAME(idamax)(n,x,incx) -1;
    }

    /// \brief BLAS level 2: Matrix vector multiply A = A := alpha*x*y**T + A 
    /// where A[m x n]
    static inline void Ger (const int& m,
                 const int& n,       const double& alpha,
                 const double* x,    const int& incx,
                 const double* y,    const int& incy,
                 double* a,    const int& lda)
    {
        F77NAME(dger) (m,n,alpha,x,incx,y,incy,a,lda);
    }

    /// \brief BLAS level 2: Matrix vector multiply A = A := alpha*x*y**T + A 
    /// where A[m x n]
    static inline void Ger (const int& m,
                 const int& n,       const float& alpha,
                 const float* x,    const int& incx,
                 const float* y,    const int& incy,
                 float* a,    const int& lda)
    {
        F77NAME(sger) (m,n,alpha,x,incx,y,incy,a,lda);
    }

    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
    static inline void Gemv (const char& trans, const int& m, const int& n,
             const double& alpha, const double* a, const int& lda,
             const double* x,     const int& incx, const double& beta,
             double* y,     const int& incy)
    {
        F77NAME(dgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
    }

    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
    static inline void Gemv (const char& trans, const int& m, const int& n,
             const float& alpha, const float* a, const int& lda,
             const float* x,     const int& incx, const float& beta,
             float* y,     const int& incy)
    {
        F77NAME(sgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
    }

    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A[m x n]
    static inline void Dgemv (const char& trans,   const int& m,    const int& n,
             const double& alpha, const double* a, const int& lda,
             const double* x,     const int& incx, const double& beta,
             double* y,     const int& incy)
    {
        F77NAME(dgemv) (trans,m,n,alpha,a,lda,x,incx,beta,y,incy);
    }

    static inline void Gbmv (const char& trans,  const int& m,
                 const int& n, const int& kl, const int& ku,
                 const double& alpha,
                 const double* a,    const int& lda,
                 const double* x,    const int& incx,
                 const double& beta, double* y, const int& incy)
    {
        F77NAME(dgbmv) (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, 
            incy);
    }

    static inline void Gbmv (const char& trans,  const int& m,
                 const int& n, const int& kl, const int& ku,
                 const float& alpha,
                 const float* a,    const int& lda,
                 const float* x,    const int& incx,
                 const float& beta, float* y, const int& incy)
    {
        F77NAME(sgbmv) (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, 
            incy);
    }

    static inline void Dgbmv (const char& trans,  const int& m,
                 const int& n, const int& kl, const int& ku,
                 const double& alpha,
                 const double* a,    const int& lda,
                 const double* x,    const int& incx,
                 const double& beta, double* y, const int& incy)
    {
        F77NAME(dgbmv) (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
    }

    static inline void Tpmv(const char& uplo, const char& trans, 
        const char& diag, const int& n, const double* ap, double* x, 
        const int& incx)
    {
        F77NAME(dtpmv) (uplo, trans, diag, n, ap, x, incx);
    }

    static inline void Tpmv(const char& uplo, const char& trans, 
        const char& diag, const int& n, const float* ap, float* x, 
        const int& incx)
    {
        F77NAME(stpmv) (uplo, trans, diag, n, ap, x, incx);
    }

    static inline void Dtpmv(const char& uplo, const char& trans, const char& diag,
                 const int& n, const double* ap, double* x, const int& incx)
    {
        F77NAME(dtpmv) (uplo, trans, diag, n, ap, x, incx);
    }

    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A
    /// is symmetric packed
    static inline void Spmv (const char& uplo, const int& n, 
            const double& alpha,
            const double* a,    const double* x, const int& incx,
            const double& beta,       double* y, const int& incy)
    {
        F77NAME(dspmv) (uplo,n,alpha,a,x,incx,beta,y,incy);
    } 
    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A
    /// is symmetric packed
    static inline void Spmv (const char& uplo, const int& n, 
            const float& alpha,
            const float* a,    const float* x, const int& incx,
            const float& beta,       float* y, const int& incy)
    {
        F77NAME(sspmv) (uplo,n,alpha,a,x,incx,beta,y,incy);
    }  

    /// \brief BLAS level 2: Matrix vector multiply y = A \e x where A
    /// is symmetric packed
    static inline void Dspmv (const char& uplo,  const int& n,    const double& alpha,
             const double* a,    const double* x, const int& incx,
             const double& beta,       double* y, const int& incy)
    {
        F77NAME(dspmv) (uplo,n,alpha,a,x,incx,beta,y,incy);
    }  

    static inline void Dsbmv (const char& uplo,   const int& m,    const int& k,
             const double& alpha, const double* a, const int& lda,
             const double* x,     const int& incx, const double& beta,
             double* y,     const int& incy)
    {
        F77NAME(dsbmv) (uplo,m,k,alpha,a,lda,x,incx,beta,y,incy);
    }

    /// \brief BLAS level 3: Matrix-matrix multiply C = A x B where op(A)[m x k],
    ///   op(B)[k x n], C[m x n]
    ///   DGEMM  performs one of the matrix-matrix operations:  
    ///   C := alpha*op( A )*op( B ) + beta*C,
    static inline void Gemm (const char& transa, const char& transb, 
        const int& m,
        const int& n,        const int& k,       const double& alpha,
        const double* a,     const int& lda,     const double* b,
        const int& ldb,      const double& beta,       double* c,
        const int& ldc)
    {
        F77NAME(dgemm) (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    }

    /// \brief BLAS level 3: Matrix-matrix multiply C = A x B where op(A)[m x k],
    ///   op(B)[k x n], C[m x n]
    ///   DGEMM  performs one of the matrix-matrix operations:  
    ///   C := alpha*op( A )*op( B ) + beta*C,
    static inline void Gemm (const char& transa,  const char& transb, 
        const int& m,
        const int& n,        const int& k,       const float& alpha,
        const float* a,     const int& lda,     const float* b,
        const int& ldb,      const float& beta,       float* c,
        const int& ldc)
    {
        F77NAME(sgemm) (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    }


    /// \brief BLAS level 3: Matrix-matrix multiply C = A x B 
    /// where op(A)[m x k], op(B)[k x n], C[m x n]
    /// DGEMM  performs one of the matrix-matrix operations:  
    /// C := alpha*op( A )*op( B ) + beta*C,
    static inline void Dgemm (const char& transa,  const char& transb, const int& m,
          const int& n,        const int& k,       const double& alpha,
          const double* a,     const int& lda,     const double* b,
          const int& ldb,      const double& beta,       double* c,
          const int& ldc)
    {
        F77NAME(dgemm) (transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
    }

    // \brief Wrapper to mutliply two (row major) matrices together C =
    // a*A*B + b*C
    static inline void Cdgemm(const int M, const int N, const int K, const double a,
          const double *A, const int ldA, const double * B, const int ldB,
          const double b, double *C, const int ldC)
    {
        boost::ignore_unused(ldA, ldB, ldC);
        Dgemm('N','N',N,M,K,a,B,N,A,K,b,C,N) ;
    }
}
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLAS_HPP
