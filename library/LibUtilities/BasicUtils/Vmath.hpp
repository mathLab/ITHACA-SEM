///////////////////////////////////////////////////////////////////////////////
//
// File Vmath.hpp
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
// Description: Collection of templated functions for vector mathematics
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP

#include <string>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Vmath 
{
    
    /***************** Math routines  ***************/
    
    /// \brief Fill a vector with a constant value
    template<class T>  LIB_UTILITIES_EXPORT void Fill( int n, const T alpha,  T *x, const int incx );



    /// \brief Generates a number from ~Normal(0,1)
    template<class T> LIB_UTILITIES_EXPORT T ran2 (long* idum);


    /// \brief Fills a vector with white noise.
    template<class T>  LIB_UTILITIES_EXPORT void FillWhiteNoise( int n, const T eps, T *x, const int incx, int seed = 0);

    /// \brief Multiply vector z = x*y
    template<class T>  LIB_UTILITIES_EXPORT void Vmul( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T*z, const int incz);

    /// \brief Scalar multiply  y = alpha*y

    template<class T>  LIB_UTILITIES_EXPORT void Smul( int n, const T alpha, const T *x, const int incx,
                                  T *y, const int incy);

    /// \brief Multiply vector z = x/y
    template<class T>  LIB_UTILITIES_EXPORT void Vdiv( int n, const T *x, const int incx, const T *y,
                  const int incy,  T*z, const int incz);
    
    /// \brief Scalar multiply  y = alpha/y
    template<class T>  LIB_UTILITIES_EXPORT void Sdiv( int n, const T alpha, const T *x, 
                                  const int incx, T *y, const int incy);

    /// \brief Add vector z = x+y
    template<class T>  LIB_UTILITIES_EXPORT void Vadd( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T *z, const int incz);
    
    /// \brief Add vector y = alpha + x
    template<class T>  LIB_UTILITIES_EXPORT void Sadd( int n, const T alpha, const T *x,
                  const int incx, T *y, const int incy);
    
    /// \brief Subtract vector z = x-y
    template<class T>  LIB_UTILITIES_EXPORT void Vsub( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T *z, const int incz);
    
    /// \brief Zero vector
    template<class T>  LIB_UTILITIES_EXPORT void Zero(int n, T *x, const int incx);
    
    /// \brief Negate x = -x
    template<class T>  LIB_UTILITIES_EXPORT void Neg( int n, T *x, const int incx);
    
    
    template<class T> void Vlog(int n, const T *x, const int incx,
                 T *y, const int incy)
    {
        while (n--)
        {
            *y = log( *x );
            x += incx;
            y += incy;
        }
    }


    template<class T> void Vexp(int n, const T *x, const int incx,
                 T *y, const int incy)
    {
        while (n--)
        {
            *y = exp( *x );
            x += incx;
            y += incy;
        }
    }

    template<class T> void Vpow(int n, const T *x, const int incx,
                const T f, T *y, const int incy)
    {
        while (n--)
        {
            *y = pow( *x, f );
            x += incx;
            y += incy;
        }
    }


    /// \brief sqrt y = sqrt(x)
    template<class T> LIB_UTILITIES_EXPORT void Vsqrt(int n, const T *x, const int incx,
                 T *y, const int incy);
    
    /// \brief vabs: y = |x|
    template<class T> LIB_UTILITIES_EXPORT void Vabs(int n, const T *x, const int incx, 
                T *y, const int incy);
    
    /********** Triad  routines  ***********************/
    
    /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvp(int n, 
                                 const T *w, const int incw, 
                                 const T *x, const int incx, 
                                 const T *y, const int incy,
                                       T *z, const int incz);
    
    /// \brief vvtvm (vector times vector plus vector): z = w*x - y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvm(int n, const T *w, const int incw, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);

    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> LIB_UTILITIES_EXPORT void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);

    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x - y
    template<class T> LIB_UTILITIES_EXPORT void Svtvm(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);

    /// \brief  vvtvvtp (vector times vector plus vector times vector): 
    // z = v*w + x*y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvvtp (int n,
                                    const T* v, int incv,
                                    const T* w, int incw,
                                    const T* x, int incx,
                                    const T* y, int incy,
                                          T* z, int incz);
    /// \brief  vvtvvtm (vector times vector minus vector times vector): 
    // z = v*w - x*y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvvtm (int n,
                                    const T* v, int incv,
                                    const T* w, int incw,
                                    const T* x, int incx,
                                    const T* y, int incy,
                                          T* z, int incz);
    /// \brief  vvtvvtp (scalar times vector plus scalar times vector): 
    // z = alpha*x + beta*y
    template<class T> LIB_UTILITIES_EXPORT void Svtsvtp (int n,
                                    const T alpha,
                                    const T* x, int incx,
                                    const T beta,
                                    const T* y, int incy,
                                          T* z, int incz);

    /// \brief  Vstvpp (scalar times vector plus vector plus vector): 
    // z = v*w + x*y
    template<class T> LIB_UTILITIES_EXPORT void Vstvpp(int n,
                                  const T alpha,
                                  const T* v, int incv,
                                  const T* w, int incw,
                                  const T* x, int incx,
                                  T* z, int incz);
    
    /************ Misc routine from Veclib (and extras)  ************/
    
    /// \brief Gather vector z[i] = x[y[i]]
    template<class T>  LIB_UTILITIES_EXPORT void Gathr(int n, const T *x, const int *y,
                  T *z);

    /// \brief Gather vector z[i] = sign[i]*x[y[i]]
    template<class T>  LIB_UTILITIES_EXPORT void Gathr(int n, const T *sign, const T *x, const int *y,
                  T *z);
    
    /// \brief Scatter vector z[y[i]] = x[i]
    template<class T>  LIB_UTILITIES_EXPORT void Scatr(int n, const T *x, const int *y,
                  T *z);

    /// \brief Scatter vector z[y[i]] = sign[i]*x[i]
    template<class T>  LIB_UTILITIES_EXPORT void Scatr(int n, const T *sign, const T *x, const int *y,
                  T *z);
    
    /// \brief Assemble z[y[i]] += x[i]; z should be zero'd first
    template<class T>  LIB_UTILITIES_EXPORT void Assmb(int n, const T *x, const int *y,
                  T *z);

    /// \brief Assemble z[y[i]] += sign[i]*x[i]; z should be zero'd first
    template<class T>  LIB_UTILITIES_EXPORT void Assmb(int n, const T *sign, const T *x, const int *y,
                  T *z);
 
    
    /************* Reduction routines  *****************/
    
    /// \brief Subtract return sum(x)
    template<class T>  LIB_UTILITIES_EXPORT T Vsum( int n, const T *x, const int incx);
    
    
    /// \brief Return the index of the maximum element in x
    template<class T>  LIB_UTILITIES_EXPORT int Imax( int n, const T *x, const int incx);
    
    /// \brief Return the maximum element in x -- called vmax to avoid
    /// conflict with max
    template<class T>  LIB_UTILITIES_EXPORT T Vmax( int n, const T *x, const int incx);
    
    /// \brief Return the index of the maximum absolute element in x
    template<class T>  LIB_UTILITIES_EXPORT int Iamax( int n, const T *x, const int incx);
    
    /// \brief Return the maximum absolute element in x
    /// called vamax to avoid conflict with max
    template<class T>  LIB_UTILITIES_EXPORT T Vamax( int n, const T *x, const int incx);
    
    
    /// \brief Return the index of the minimum element in x
    template<class T>  LIB_UTILITIES_EXPORT int Imin( int n, const T *x, const int incx);
    
    
    /// \brief Return the minimum element in x - called vmin to avoid
    /// conflict with min
    template<class T>  LIB_UTILITIES_EXPORT T Vmin( int n, const T *x, const int incx);
    
    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> LIB_UTILITIES_EXPORT T Dot(     int n,
                                 const T   *w,
                                 const T   *x);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> LIB_UTILITIES_EXPORT T Dot(     int n,
                                 const T   *w, const int incw,
                                 const T   *x, const int incx);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> LIB_UTILITIES_EXPORT T Dot2(    int n,
                                 const T   *w,
                                 const T   *x,
                                 const int *y);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> LIB_UTILITIES_EXPORT T Dot2(    int n,
                                 const T   *w, const int incw,
                                 const T   *x, const int incx,
                                 const int *y, const int incy);

    /********** Memory routines  ***********************/

#if 0     
    // \brief copy one double vector to another - This is just a wrapper
    // around Blas
    static inline void Vcopy(int n, const double *x, int incx, double *y,
                             const int incy)
    {
        Blas::Dcopy(n,x,incx,y,incy);
    }
#else
    // \brief copy one int vector to another
    LIB_UTILITIES_EXPORT void Vcopy(int n, const int *x, const int incx, int *y,
                             const int incy);

    // \brief copy one int vector to another
    LIB_UTILITIES_EXPORT void Vcopy(int n, const unsigned int *x, const int incx, unsigned int *y,
                             const int incy);

    // \brief copy one double vector to another
    LIB_UTILITIES_EXPORT void Vcopy(int n, const double *x, const int incx, double *y,
                             const int incy);

    // \brief reverse the ordering of  vector to another
    template<class T>  LIB_UTILITIES_EXPORT void  Reverse( int n, const T *x, const int incx, T *y, const int incy);

#endif

}
#endif //VECTORMATH_HPP

/***
$Log: Vmath.hpp,v $
Revision 1.21  2009/05/15 14:38:41  pvos
Changed check for regular quads so that it also includes parallellograms

Revision 1.20  2009/03/10 23:44:15  claes
Made y in z = x/y a constant in the parameter list.

Revision 1.19  2009/01/21 16:57:26  pvos
Added additional geometric factors to improve efficiency

Revision 1.18  2008/12/17 16:56:46  pvos
Performance updates

Revision 1.17  2008/09/09 14:00:55  sherwin
Fixed error in Sdiv definition

Revision 1.16  2008/08/09 19:26:08  sherwin
Corrected big in Reverse

Revision 1.15  2008/07/19 21:09:21  sherwin
Added Reverse function

Revision 1.14  2008/05/10 18:27:32  sherwin
Modifications necessary for QuadExp Unified DG Solver

Revision 1.13  2008/04/06 05:47:03  bnelson
Fixed gcc compiler warnings.

Revision 1.12  2008/03/06 04:39:55  ehan
Removed the include file <VmathArray.hpp>.

Revision 1.11  2008/02/28 09:55:57  sherwin
Added Array version of math routines

Revision 1.10  2008/01/27 09:13:04  sherwin
Added Svtvp routine

Revision 1.9  2007/12/06 22:43:57  pvos
2D Helmholtz solver updates

Revision 1.8  2007/07/26 08:32:57  sherwin
Made second vector argument to Vsub a constant

Revision 1.7  2007/07/20 00:39:37  bnelson
Replaced boost::shared_ptr with Nektar::ptr

Revision 1.6  2007/04/08 03:30:25  jfrazier
Made y in z = x*y a constant in the parameter list.

Revision 1.5  2007/04/03 03:51:44  bnelson
Moved Lapack.hpp, Blas.hpp, Transf77.hpp to LinearAlgebra

Revision 1.4  2007/01/18 20:59:26  sherwin
Before new configuration

Revision 1.3  2006/07/02 17:16:16  sherwin

Modifications to make MultiRegions work for a connected domain in 2D (Tris)

Revision 1.2  2006/06/01 13:44:28  kirby
*** empty log message ***

Revision 1.1  2006/06/01 11:07:52  kirby
*** empty log message ***

Revision 1.1  2006/05/04 18:57:44  kirby
*** empty log message ***

Revision 1.5  2006/03/01 08:25:03  sherwin

First compiling version of StdRegions

Revision 1.4  2006/02/26 21:13:45  bnelson
Fixed a variety of compiler errors caused by updates to the coding standard.

Revision 1.3  2006/02/15 08:07:15  sherwin

Put codes into standard although have not yet been compiled

Revision 1.2  2006/02/12 21:51:42  sherwin

Added licence

Revision 1.1  2006/02/12 15:06:12  sherwin

Changed .h files to .hpp

**/
