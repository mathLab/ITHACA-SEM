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

namespace Vmath 
{
    
    /***************** Math routines  ***************/
    
    /// \brief Fill a vector with a constant value
    template<class T>  void Fill( int n, const T alpha,  T *x, const int incx )
    {
        while( n-- )
        {
            *x = alpha;
             x += incx;
        }
    }

    /// \brief Multiply vector z = x*y
    template<class T>  void Vmul( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T*z, const int incz)
    {
        while( n-- )
        {
            *z = (*x) * (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }

    /// \brief Scalar multiply  y = alpha*y

    template<class T>  void Smul( int n, const T alpha, const T *x, const int incx,
                                  T *y, const int incy)
    {
        while( n-- )
        {
            *y = alpha * (*x);
            x += incx;
            y += incy;
        }
    }

    /// \brief Multiply vector z = x/y
    template<class T>  void Vdiv( int n, const T *x, const int incx, T *y,
                  const int incy,  T*z, const int incz)
    {
        while( n-- )
        {
            *z = (*x) / (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
    
    /// \brief Scalar multiply  y = alpha/y
    template<class T>  void Sdiv( int n, const T alpha, T*x, const int incx,
                  T *y, const int incy)
    {
        while( n-- )
        {
            *y = alpha / (*x);
            x += incx;
            y += incy;
        }
    }

    /// \brief Add vector z = x+y
    template<class T>  void Vadd( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T *z, const int incz)
    {
        while( n-- )
        {
            *z = (*x) + (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
    
    /// \brief Add vector y = alpha + x
    template<class T>  void Sadd( int n, const T alpha, const T *x,
                  const int incx, T *y, const int incy)
    {
        while( n-- )
        {
            *y = alpha + (*x);
            x += incx;
            y += incy;
        }
    }
    
    /// \brief Subtract vector z = x-y
    template<class T>  void Vsub( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T *z, const int incz)
    {
        while( n-- )
        {
            *z = (*x) - (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
    
    /// \brief Zero vector
    template<class T>  void Zero(int n, T *x, const int incx)
    {
        if(incx == 1)
        {
            std::memset(x,'\0', n*sizeof(T));
        }
        else
        {
            T zero = 0;
            while(n--)
            {
            *x = zero;
            x+=incx;
            }
        }
    }
    
    /// \brief Negate x = -x
    template<class T>  void Neg( int n, T *x, const int incx)
    {
        while( n-- )
        {
            *x = -(*x);
            x += incx;
        }
    }
    
    
    /// \brief sqrt y = sqrt(x)
    template<class T> void Vsqrt(int n, const T *x, const int incx,
                 T *y, const int incy)
    {
        while (n--)
        {
            *y  = sqrt( *x );
            x  += incx;
            y  += incy;
        }
    }
    
    /// \brief vabs: y = |x|
    template<class T> void Vabs(int n, const T *x, const int incx, 
                T *y, const int incy)
    {
        while( n-- )
        {
            *y = ( *x >0)? *x:-(*x);
            x += incx;
            y += incy;
        }
    }
    
    /********** Triad  routines  ***********************/
    
    /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
    template<class T> void Vvtvp(int n, const T *w, const int incw, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        while( n-- )
        {
            *z = (*w) * (*x) + (*y);
            w += incw;
            x += incx;
            y += incy;
            z += incz;
        }
    }

    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        while( n-- )
        {
            *z = alpha * (*x) + (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }
    
    /// \brief vvtvm (vector times vector plus vector): z = w*x - y
    template<class T> void Vvtvm(int n, const T *w, const int incw, T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        while( n-- )
        {
            *z = (*w) * (*x) - (*y);
            w += incw;
            x += incx;
            y += incy;
            z += incz;
        }
    }
    
    
    /************ Misc routine from Veclib (and extras)  ************/
    
    /// \brief Gather vector z[i] = x[y[i]]
    template<class T>  void Gathr(int n, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *z++ = *(x + *y++);
        }
        return;
    }
    
    /// \brief Scatter vector z[y[i]] = x[i]
    template<class T>  void Scatr(int n, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *(z + *(y++)) = *(x++);
        }
    }
    
    
    /// \brief Assemble z[y[i]] += x[i]; z should be zero'd first
    template<class T>  void Assmb(int n, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *(z + *(y++)) += *(x++);
        }
    }
    
    
    /************* Reduction routines  *****************/
    
    /// \brief Subtract return sum(x)
    template<class T>  T Vsum( int n, const T *x, const int incx)
    {
    
        T sum = 0;
        
        while( n-- )
        {
            sum += (*x);
            x += incx;
        }
        
        return sum;
    }
    
    
    /// \brief Return the index of the maximum element in x
    template<class T>  int Imax( int n, const T *x, const int incx)
    {
    
        int    i, indx = ( n > 0 ) ? 0 : -1;
        T      xmax = *x;
        
        for (i = 0; i < n; i++)
        {
            if (*x > xmax)
            {
            xmax = *x;
            indx = i;
            }
            x += incx;
        }
        
        return indx;
    }
    
    /// \brief Return the maximum element in x -- called vmax to avoid
    /// conflict with max
    template<class T>  T Vmax( int n, const T *x, const int incx)
    {
    
        T  xmax = *x;
        
        while( n-- )
        {
            if (*x > xmax)
            {
            xmax = *x;
            }
            x += incx;
        }
        
        return xmax;
    }
    
    /// \brief Return the index of the maximum absolute element in x
    template<class T>  int Iamax( int n, const T *x, const int incx)
    {
    
        int    i, indx = ( n > 0 ) ? 0 : -1;
        T      xmax = *x;
        T      xm;
        
        for (i = 0; i < n; i++)
        {
            xm = (*x > 0)? *x: -*x;
            if (xm > xmax)
            {
            xmax = xm;
            indx = i;
            }
            x += incx;
        }
        
        return indx;
    }
    
    /// \brief Return the maximum absolute element in x
    /// called vamax to avoid conflict with max
    template<class T>  T Vamax( int n, const T *x, const int incx)
    {
    
        T  xmax = *x;
        T  xm;
        
        while( n-- )
        {
            xm = (*x > 0)? *x: -*x;
            if (xm > xmax)
            {
            xmax = xm;
            }
            x += incx;
        }
        return xmax;
    }
    
    
    /// \brief Return the index of the minimum element in x
    template<class T>  int Imin( int n, const T *x, const int incx)
    {
    
        int    i, indx = ( n > 0 ) ? 0 : -1;
        T      xmin = *x;
        
        for(i = 0;i < n;i++)
        {
            if( *x < xmin )
            {
            xmin = *x;
            indx = i;
            }
            x += incx;
        }
        
        return indx;
    }
    
    
    /// \brief Return the minimum element in x - called vmin to avoid
    /// conflict with min
    template<class T>  T Vmin( int n, const T *x, const int incx)
    {
    
        T    xmin = *x;
        
        while( n-- )
        {
            if (*x < xmin)
            {
            xmin = *x;
            }
            x += incx;
        }
        
        return xmin;
    }
    
    /********** Memory routines  ***********************/

#if 0     
    // \brief copy one double vector to another - This is just a wrapper
    // around Blas
    static inline void Vcopy(int n, const double *x, int incx, double *y,
              int const incy)
    {
        Blas::Dcopy(n,x,incx,y,incy);
    }
#else
    // \brief copy one int vector to another
    static inline void Vcopy(int n, const int *x, const int incx, int *y,
              int const incy)
    {
        if( incx ==1 && incy == 1)
        {
            memcpy(y,x,n*sizeof(int));
        }
        else
        {
            while( n-- )
            {
                *y = *x;
                x += incx;
                y += incy;
            }
        }
    }

    // \brief copy one double vector to another
    static inline void Vcopy(int n, const double *x, const int incx, double *y,
              int const incy)
    {
        if( incx ==1 && incy == 1)
        {
            memcpy(y,x,n*sizeof(double));
        }
        else
        {
            while( n-- )
            {
                *y = *x;
                x += incx;
                y += incy;
            }
        }
    }
#endif

    
}
#endif //VECTORMATH_HPP

/***
$Log: Vmath.hpp,v $
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
