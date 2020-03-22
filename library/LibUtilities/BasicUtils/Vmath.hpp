///////////////////////////////////////////////////////////////////////////////
//
// File: Vmath.hpp
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
// Description: Collection of templated functions for vector mathematics
//
// Note: For those unfamiliar with the vector routines notation, it is
//       reverse polish notation (RPN).  For example:
//
//       In the function "Vvtvp()", it is "Vvt" means vector vector times,
//       which in infix notation is "v * v".  So "Vvtvp" is:
//
//       RPN:    vector vector times vector plus
//       Infix:  ( vector * vector ) + vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATH_HPP

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
    template<class T>  LIB_UTILITIES_EXPORT void FillWhiteNoise(
                    int n, const T eps, T *x, const int incx, int seed = 9999);

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
    template<class T> LIB_UTILITIES_EXPORT void Vvtvp( int n, 
                                                       const T *w, const int incw, 
                                                       const T *x, const int incx, 
                                                       const T *y, const int incy,
                                                       T *z,       const int incz );
    
    /// \brief vvtvm (vector times vector plus vector): z = w*x - y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvm(int n, const T *w, const int incw, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);

    /// \brief  Svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> LIB_UTILITIES_EXPORT void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);

    /// \brief  Svtvm (scalar times vector minus vector): z = alpha*x - y
    template<class T> LIB_UTILITIES_EXPORT void Svtvm ( int n, const T alpha, const T *x,
                                                        const int incx, const T *y, const int incy,
                                                        T *z, const int incz );

    /// \brief  vvtvvtp (vector times vector plus vector times vector): 
    // z = v*w + x*y
    //
    // @param n  Number of items in each vector.
    template<class T> LIB_UTILITIES_EXPORT void Vvtvvtp (       int n,
                                                          const T*  v, int incv,
                                                          const T*  w, int incw,
                                                          const T*  x, int incx,
                                                          const T*  y, int incy,
                                                                T*  z, int incz );

    /// \brief  vvtvvtm (vector times vector minus vector times vector): 
    // z = v*w - x*y
    template<class T> LIB_UTILITIES_EXPORT void Vvtvvtm (int n,
                                    const T* v, int incv,
                                    const T* w, int incw,
                                    const T* x, int incx,
                                    const T* y, int incy,
                                          T* z, int incz);
    /// \brief  Svtsvtp (scalar times vector plus scalar times vector): 
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

    /// \brief Return number of NaN elements of x
    template<class T>  LIB_UTILITIES_EXPORT int Nnan(int n, const T *x, const int incx);

    
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

    // \brief copy one vector to another
    template<class T>
    LIB_UTILITIES_EXPORT void Vcopy(int n, const T *x, const int incx,
                                                 T *y, const int incy);

    // \brief reverse the ordering of  vector to another
    template<class T>  LIB_UTILITIES_EXPORT void  Reverse( int n, const T *x, const int incx, T *y, const int incy);


}
#endif //VECTORMATH_HPP

