///////////////////////////////////////////////////////////////////////////////
//
// File Lapack.hpp
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
// Description: wrapper of functions around standard LAPACK routines
//
///////////////////////////////////////////////////////////////////////////////
#ifndef  LAPACK_HPP
#define  LAPACK_HPP

#include <LibUtilities/BasicUtils/TransF77.hpp>

namespace Lapack
{

  extern "C"
  {
    // Matrix factorisation and solves
    void F77NAME(dsptrf) (const char& uplo, const int& n,
              double* ap, int *ipiv, int& info);
    void F77NAME(dsptrs) (const char& uplo, const int& n,
              const int& nrhs, const double* ap,
              const int  *ipiv, double* b,
              const int& ldb, int& info);
    void F77NAME(dpptrf) (const char& uplo, const int& n,
              double* ap, int& info);
    void F77NAME(dpptrs) (const char& uplo, const int& n,
              const int& nrhs, const double* ap,
              double* b, const int& ldb, int& info);
    void F77NAME(dpbtrf) (const char& uplo, const int& n, const int& kd,
              double* ab, const int& ldab, int& info);
    void F77NAME(dpbtrs) (const char& uplo, const int& n,
              const int& kd, const int& nrhs,
              const double* ab, const int& ldab,
              double* b, const int& ldb, int& info);
    void F77NAME(dgbtrf) (const int& m, const int& n, const int& kl,
              const int& ku, double* a, const int& lda,
              int* ipiv, int& info);
    void F77NAME(dgbtrs) (const char& trans, const int& n, const int& kl,
              const int &ku, const int& nrhs,   const double* a,
              const int& lda, int* ipiv, double* b,
              const int& ldb, int& info);
    void F77NAME(dgetrf) (const int& m, const int& n, double* a,
              const int& lda, int* ipiv, int& info);
    void F77NAME(dgetrs) (const char& trans, const int& n, const int& nrhs,
              const double* a,   const int& lda, int* ipiv,
              double* b, const int& ldb, int& info);
    void F77NAME(dgetri) (const int& n, double *a, const int& lda,
            const int *ipiv, double *wk,  const int& lwk,
              int& info);
    void F77NAME(dsterf) (const int& n, double *d, double *e, int& info);
    void F77NAME(dgeev)  (const char& uplo, const char& lrev, const int& n,
              double* a, const int& lda, double* wr, double* wi,
              double* rev,  const int& ldr,
              double* lev,  const int& ldv,
              double* work, const int& lwork, int& info);

    void F77NAME(dspev)  (const char& jobz, const char& uplo, const int& n,
              double* ap, double* w, double* z, const int& ldz,
              double* work, int& info);
    void F77NAME(dsbev)  (const char& jobz, const char& uplo, const int& kl,
              const int& ku,  double* ap, const int& lda,
              double* w, double* z, const int& ldz,
              double* work, int& info);
  }


  /// \brief factor a real packed-symmetric matrix using Bunch-Kaufman
  /// pivoting.
  static void Dsptrf (const char& uplo, const int& n,
              double* ap, int *ipiv, int& info)
  {
    F77NAME(dsptrf) (uplo,n,ap,ipiv,info);
  }

  /// \brief Solve a real symmetric matrix problem using Bunch-Kaufman
  /// pivoting.
  static void Dsptrs (const char& uplo, const int& n, const int& nrhs,
              const double* ap, const int  *ipiv, double* b,
              const int& ldb, int& info)
  {
    F77NAME(dsptrs) (uplo,n,nrhs,ap,ipiv,b,ldb,info);
  }

  /// \brief Cholesky factor a real Positive Definite packed-symmetric matrix.
  static void Dpptrf (const char& uplo, const int& n,
              double *ap, int& info)
  {
    F77NAME(dpptrf) (uplo,n,ap,info);
  }

  /// \brief Solve a real Positive defiinte symmetric matrix problem
  /// using Cholesky factorization.
  static void Dpptrs (const char& uplo, const int& n, const int& nrhs,
              const double *ap, double *b, const int& ldb,
              int& info)
  {
    F77NAME(dpptrs) (uplo,n,nrhs,ap,b,ldb,info);
  }

  /// \brief Cholesky factorize a real positive-definite
  /// banded-symmetric matrix
  static void Dpbtrf (const char& uplo, const int& n, const int& kd,
             double *ab, const int& ldab, int& info)
  {
    F77NAME(dpbtrf) (uplo,n,kd,ab,ldab,info);
  }


  /// \brief Solve a real, Positive definite banded symmetric matrix
  /// problem using Cholesky factorization.
  static void Dpbtrs (const char& uplo, const int& n,
              const int& kd, const int& nrhs,
              const double *ab, const int& ldab,
              double *b, const int& ldb, int& info)
  {
    F77NAME(dpbtrs) (uplo,n,kd,nrhs,ab,ldab,b,ldb,info);
  }

  /// \brief General banded matrix LU factorisation
  static void Dgbtrf (const int& m, const int& n, const int& kl,
              const int& ku, double* a, const int& lda,
              int* ipiv, int& info)
  {
    F77NAME(dgbtrf)(m,n,kl,ku,a,lda,ipiv,info);
  }

  /// \brief Solve general banded matrix using LU factorisation
  static void Dgbtrs (const char& trans, const int& n, const int& kl,
              const int &ku, const int& nrhs,   const double* a,
              const int& lda, int* ipiv, double* b,
              const int& ldb, int& info)
  {
    F77NAME(dgbtrs)(trans,n,kl,ku,nrhs,a,lda,ipiv,b,ldb,info);
  }

  /// \brief General matrix LU factorisation
  static void Dgetrf (const int& m, const int& n, double *a,
              const int& lda, int *ipiv, int& info)
  {
    F77NAME(dgetrf) (m,n,a,lda,ipiv,info);
  }

  /// \brief General matrix LU backsolve
  static void Dgetrs (const char& trans, const int& n, const int& nrhs,
              const double* a, const int& lda, int* ipiv,
              double* b, const int& ldb, int& info)
  {
    F77NAME(dgetrs) (trans,n,nrhs,a,lda,ipiv,b,ldb,info);
  }

  /// \brief Generate matrix inverse
  static void Dgetri (const int& n, double *a, const int& lda,
           const int *ipiv, double *wk,  const int& lwk, int& info)
  {
    F77NAME(dgetri) (n, a, lda, ipiv, wk, lwk,info);
  }

  /// \brief Find eigenvalues of symmetric tridiagonal matrix
  static void Dsterf(const int& n, double *d, double *e, int& info)
  {
    F77NAME(dsterf)(n,d,e,info);
  }

  /// \brief Solve general real matrix eigenproblem.
  static void Dgeev (const char& uplo, const char& lrev, const int& n,
             double* a, const int& lda, double* wr, double* wi,
             double* rev,  const int& ldr,
             double* lev,  const int& ldv,
             double* work, const int& lwork, int& info)
  {
    F77NAME(dgeev) (uplo, lrev, n, a, lda, wr, wi, rev,
            ldr, lev, ldv, work, lwork, info);
  }

  /// \brief Solve packed-symmetric real matrix eigenproblem.
  static void Dspev (const char& jobz, const char& uplo, const int& n,
             double* ap, double* w, double* z, const int& ldz,
             double* work, int& info)
  {
    F77NAME(dspev) (jobz, uplo, n, ap, w, z, ldz, work, info);
  }

  /// \brief Solve packed-banded real matrix eigenproblem.
  static void Dsbev (const char& jobz, const char& uplo, const int& kl,
             const int& ku,  double* ap, const int& lda,
             double* w, double* z, const int& ldz,
             double* work, int& info)
  {
    F77NAME(dsbev) (jobz, uplo, kl, ku, ap, lda, w, z, ldz, work, info);
  }
}
#endif //LAPACK_HPP

/***
$Log: Lapack.hpp,v $
Revision 1.1  2006/06/01 11:07:52  kirby
*** empty log message ***

Revision 1.1  2006/05/04 18:57:43  kirby
*** empty log message ***

Revision 1.4  2006/02/26 21:13:45  bnelson
Fixed a variety of compiler errors caused by updates to the coding standard.

Revision 1.3  2006/02/15 08:07:15  sherwin

Put codes into standard although have not yet been compiled

Revision 1.2  2006/02/12 21:51:42  sherwin

Added licence

Revision 1.1  2006/02/12 15:06:12  sherwin

Changed .h files to .hpp

**/
