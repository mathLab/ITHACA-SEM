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
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

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

    template LIB_UTILITIES_EXPORT void Fill( int n, const Nektar::NekDouble alpha,  Nektar::NekDouble *x, const int incx );

    #define IM1   2147483563
    #define IM2   2147483399
    #define AM    (1.0/IM1)
    #define IMM1  (IM1-1)
    #define IA1   40014
    #define IA2   40692
    #define IQ1   53668
    #define IQ2   52774
    #define IR1   12211
    #define IR2   3791
    #define NTAB  32
    #define NDIV  (1+IMM1/NTAB)
    #define EPS   1.2e-7
    #define RNMX  (1.0-EPS)

    /// \brief Generates a number from ~Normal(0,1)
    template<class T> T ran2 (long* idum)
    /* ------------------------------------------------------------------------- *
     * Ran2 from NR 2e.  Returns a uniform random deviate between 0.0 &
     * 1.0 (exclusive of endpoints).  Call with idum a negative integer to
     * initialize; thereafter, do not alter idum between successive
     * deviates in a sequence.  RNMX should approximate the largest
     * floating value that is less than 1.
     * ------------------------------------------------------------------------- */
    {
      int         j;
      long        k;
      static long idum2=123456789;
      static long iy=0;
      static long iv[NTAB];
      T           temp;

      if (*idum <= 0) {
        if (-(*idum) < 1) *idum = 1;
        else              *idum = -(*idum);
        idum2 = (*idum);
        for (j=NTAB+7; j>=0; j--) {
          k = (*idum) / IQ1;
          *idum = IA1 * (*idum - k*IQ1) - k*IR1;
          if (*idum < 0) *idum += IM1;
          if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
      }

      k = (*idum) / IQ1;
      *idum = IA1*(*idum - k*IQ1) - k*IR1;
      if (*idum < 0) *idum += IM1;

      k = idum2 / IQ2;
      idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
      if (idum2 < 0) idum2 += IM2;

      j = iy / NDIV;
      iy = iv[j] - idum2;
      iv[j] = *idum;
      if (iy < 1) iy += IMM1;

      if ((temp=AM*iy) > RNMX) return RNMX;
      else                     return temp;
    }

    #undef IM1
    #undef IM2
    #undef AM
    #undef IMM1
    #undef IA1
    #undef IA2
    #undef IQ1
    #undef IQ2
    #undef IR1
    #undef IR2
    #undef NTAB
    #undef NDIV
    #undef EPS
    #undef RNMX

#ifdef NEKTAR_USE_THREAD_SAFETY
    static boost::mutex mutex;
#endif
    template LIB_UTILITIES_EXPORT Nektar::NekDouble ran2 (long* idum);

    /// \brief Fills a vector with white noise.
    template<class T>  void FillWhiteNoise( int n, const T eps, T *x,
                                      const int incx, int outseed)
    {
#ifdef NEKTAR_USE_THREAD_SAFETY
        // Protect the static vars here and in ran2
        boost::mutex::scoped_lock l(mutex);
#endif

        // Define static variables for generating random numbers
        static int     iset = 0;
        static T       gset;
        static long    seed = 0;

        // Bypass seed if outseed was specified
        if( outseed != 9999)
        {
            seed = long(outseed);
        }

        while( n-- )
        {
            T              fac, rsq, v1, v2;

            if (iset == 0)
            {
                do
                {
                    v1 = 2.0 * ran2<T> (&seed) - 1.0;
                    v2 = 2.0 * ran2<T> (&seed) - 1.0;
                    rsq = v1*v1 + v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac = sqrt(-2.0 * log (rsq) / rsq);
                gset = v1 * fac;
                iset = 1;
                *x = eps * v2 * fac;
            }
            else
            {
                iset = 0;
                *x = eps * gset;
            }
            x += incx;
        }
    }
    template  LIB_UTILITIES_EXPORT void FillWhiteNoise( int n, const Nektar::NekDouble eps, Nektar::NekDouble *x, const int incx, int outseed);

    /// \brief Multiply vector z = x*y
    template<class T>  void Vmul( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T*z, const int incz)
    {
        ++n;
        if (incx == 1 && incy == 1 && incz == 1)
        {
            while( --n )
            {
                *z = (*x) * (*y);
                ++x;
                ++y;
                ++z;
            }
        }
        else
        {
            while( --n )
            {
                *z = (*x) * (*y);
                x += incx;
                y += incy;
                z += incz;
            }
        }
    }
    template  LIB_UTILITIES_EXPORT void Vmul( int n, const Nektar::NekDouble *x, const int incx, const Nektar::NekDouble *y,
                              const int incy,  Nektar::NekDouble*z, const int incz);

    /// \brief Scalar multiply  y = alpha*x

    template<class T>  void Smul( int n, const T alpha, const T *x, const int incx,
                                  T *y, const int incy)
    {
        ++n;
        if (incx == 1 && incy == 1)
        {
            while( --n )
            {
                *y = alpha * (*x);
                ++x;
                ++y;
            }
        }
        else
        {
            while( --n )
            {
                *y = alpha * (*x);
                x += incx;
                y += incy;
            }
        }
    }

    template  LIB_UTILITIES_EXPORT void Smul( int n, const Nektar::NekDouble alpha, const Nektar::NekDouble *x, const int incx,
                              Nektar::NekDouble *y, const int incy);

    /// \brief Multiply vector z = x/y
    template<class T>  void Vdiv( int n, const T *x, const int incx, const T *y,
                  const int incy,  T*z, const int incz)
    {
        ++n;
        if (incx == 1 && incy == 1)
        {
            while( --n )
            {
                *z = (*x) / (*y);
                ++x;
                ++y;
                ++z;
            }
        }
        else
        {
            while( --n )
            {
                *z = (*x) / (*y);
                x += incx;
                y += incy;
                z += incz;
            }
        }
    }

    template LIB_UTILITIES_EXPORT void Vdiv( int n, const Nektar::NekDouble *x, const int incx, const Nektar::NekDouble *y,
              const int incy,  Nektar::NekDouble*z, const int incz);

    /// \brief Scalar multiply  y = alpha/y
    template<class T>  void Sdiv( int n, const T alpha, const T *x,
                                  const int incx, T *y, const int incy)
    {
        ++n;
        if (incx == 1 && incy == 1)
        {
            while( --n )
            {
                *y = alpha / (*x);
                ++x;
                ++y;
            }
        }
        else
        {
            while( --n )
            {
                *y = alpha / (*x);
                x += incx;
                y += incy;
            }
        }
    }

    template  LIB_UTILITIES_EXPORT void Sdiv( int n, const Nektar::NekDouble alpha, const Nektar::NekDouble *x,
                              const int incx, Nektar::NekDouble *y, const int incy);

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

    template  LIB_UTILITIES_EXPORT void Vadd( int n, const Nektar::NekDouble *x, const int incx, const Nektar::NekDouble *y,
                              const int incy,  Nektar::NekDouble *z, const int incz);

    /// \brief Add vector y = alpha + x
    template<class T>  void Sadd( int n, const T alpha, const T *x,
                  const int incx, T *y, const int incy)
    {
        ++n;
        if (incx == 1 && incy == 1)
        {
            while( --n )
            {
                *y = alpha + (*x);
                ++x;
                ++y;
            }
        }
        else
        {
            while( --n )
            {
                *y = alpha + (*x);
                x += incx;
                y += incy;
            }
        }
    }

    template LIB_UTILITIES_EXPORT void Sadd( int n, const Nektar::NekDouble alpha, const Nektar::NekDouble *x,
              const int incx, Nektar::NekDouble *y, const int incy);

    /// \brief Subtract vector z = x-y
    template<class T>  void Vsub( int n, const T *x, const int incx, const T *y,
                                  const int incy,  T *z, const int incz)
    {
        ++n;
        if (incx == 1 && incy == 1 && incz == 1)
        {
            while( --n )
            {
                *z = (*x) - (*y);
                ++x;
                ++y;
                ++z;
            }
        }
        else
        {
            while( --n )
            {
                *z = (*x) - (*y);
                x += incx;
                y += incy;
                z += incz;
            }
        }
    }

    template  LIB_UTILITIES_EXPORT void Vsub( int n, const Nektar::NekDouble *x, const int incx, const Nektar::NekDouble *y,
                              const int incy,  Nektar::NekDouble *z, const int incz);

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
            ++n;
            while(--n)
            {
            *x = zero;
            x+=incx;
            }
        }
    }

    template  LIB_UTILITIES_EXPORT void Zero(int n, Nektar::NekDouble *x, const int incx);
    template  LIB_UTILITIES_EXPORT void Zero(int n, int *x, const int incx);
    template  LIB_UTILITIES_EXPORT void Zero(int n, long *x, const int incx);

    /// \brief Negate x = -x
    template<class T>  void Neg( int n, T *x, const int incx)
    {
        while( n-- )
        {
            *x = -(*x);
            x += incx;
        }
    }

    template  LIB_UTILITIES_EXPORT void Neg( int n, Nektar::NekDouble *x, const int incx);

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

    template LIB_UTILITIES_EXPORT void Vsqrt(int n, const Nektar::NekDouble *x, const int incx,
             Nektar::NekDouble *y, const int incy);


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

    template LIB_UTILITIES_EXPORT void Vabs(int n, const Nektar::NekDouble *x, const int incx,
            Nektar::NekDouble *y, const int incy);


    /********** Triad  routines  ***********************/

    /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
    template<class T> void Vvtvp(int n,
                                 const T *w, const int incw,
                                 const T *x, const int incx,
                                 const T *y, const int incy,
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

    template LIB_UTILITIES_EXPORT void Vvtvp(int n,
                             const Nektar::NekDouble *w, const int incw,
                             const Nektar::NekDouble *x, const int incx,
                             const Nektar::NekDouble *y, const int incy,
                                   Nektar::NekDouble *z, const int incz);

    /// \brief vvtvm (vector times vector plus vector): z = w*x - y
    template<class T> void Vvtvm(int n, const T *w, const int incw, const T *x,
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

    template LIB_UTILITIES_EXPORT void Vvtvm(int n, const Nektar::NekDouble *w, const int incw, const Nektar::NekDouble *x,
             const int incx, const Nektar::NekDouble *y, const int incy,
             Nektar::NekDouble *z, const int incz);


    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> LIB_UTILITIES_EXPORT void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        ++n;
        if (incx == 1 && incy == 1 && incz == 1)
        {
            while( --n )
            {
                *z = alpha * (*x) + (*y);
                ++x;
                ++y;
                ++z;
            }
        }
        else
        {
            while( --n )
            {
                *z = alpha * (*x) + (*y);
                x += incx;
                y += incy;
                z += incz;
            }
        }
    }

    template LIB_UTILITIES_EXPORT void Svtvp(int n, const Nektar::NekDouble alpha, const Nektar::NekDouble *x,
                 const int incx, const Nektar::NekDouble *y, const int incy,
                 Nektar::NekDouble *z, const int incz);


    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x - y
    template<class T> void Svtvm(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        while( n-- )
        {
            *z = alpha * (*x) - (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }

    template LIB_UTILITIES_EXPORT void Svtvm(int n, const Nektar::NekDouble alpha, const Nektar::NekDouble *x,
             const int incx, const Nektar::NekDouble *y, const int incy,
             Nektar::NekDouble *z, const int incz);

    /// \brief  vvtvvtp (vector times vector plus vector times vector):
    // z = v*w + x*y
    template<class T> void Vvtvvtp (int n,
                                    const T* v, int incv,
                                    const T* w, int incw,
                                    const T* x, int incx,
                                    const T* y, int incy,
                                          T* z, int incz)
    {
        while( n-- )
        {
            *z = (*v) * (*w) + (*x) * (*y);
            v += incv;
            w += incw;
            x += incx;
            y += incy;
            z += incz;
        }
    }
    template LIB_UTILITIES_EXPORT void Vvtvvtp (int n,
                                const Nektar::NekDouble* v, int incv,
                                const Nektar::NekDouble* w, int incw,
                                const Nektar::NekDouble* x, int incx,
                                const Nektar::NekDouble* y, int incy,
                                      Nektar::NekDouble* z, int incz);


    /// \brief  vvtvvtm (vector times vector minus vector times vector):
    // z = v*w - x*y
    template<class T> void Vvtvvtm (int n,
                                    const T* v, int incv,
                                    const T* w, int incw,
                                    const T* x, int incx,
                                    const T* y, int incy,
                                          T* z, int incz)
    {
        while( n-- )
        {
            *z = (*v) * (*w) - (*x) * (*y);
            v += incv;
            w += incw;
            x += incx;
            y += incy;
            z += incz;
        }
    }

    template LIB_UTILITIES_EXPORT void Vvtvvtm (int n,
                                const Nektar::NekDouble* v, int incv,
                                const Nektar::NekDouble* w, int incw,
                                const Nektar::NekDouble* x, int incx,
                                const Nektar::NekDouble* y, int incy,
                                      Nektar::NekDouble* z, int incz);

    /// \brief  vvtvvtp (scalar times vector plus scalar times vector):
    // z = alpha*x + beta*y
    template<class T> void Svtsvtp (int n,
                                    const T alpha,
                                    const T* x, int incx,
                                    const T beta,
                                    const T* y, int incy,
                                          T* z, int incz)
    {
        while( n-- )
        {
            *z = alpha * (*x) + beta * (*y);
            x += incx;
            y += incy;
            z += incz;
        }
    }

    template LIB_UTILITIES_EXPORT void Svtsvtp (int n,
                                const Nektar::NekDouble alpha,
                                const Nektar::NekDouble* x, int incx,
                                const Nektar::NekDouble beta,
                                const Nektar::NekDouble* y, int incy,
                                      Nektar::NekDouble* z, int incz);


    /// \brief  Vstvpp (scalar times vector plus vector plus vector):
    // z = v*w + x*y
    template<class T> void Vstvpp(int n,
                                  const T alpha,
                                  const T* v, int incv,
                                  const T* w, int incw,
                                  const T* x, int incx,
                                  T* z, int incz)
    {
        while( n-- )
        {
            *z = alpha * (*v) + (*w) + (*x);
            v += incv;
            w += incw;
            x += incx;
            z += incz;
        }
    }

    template LIB_UTILITIES_EXPORT void Vstvpp(int n,
                              const Nektar::NekDouble alpha,
                              const Nektar::NekDouble* v, int incv,
                              const Nektar::NekDouble* w, int incw,
                              const Nektar::NekDouble* x, int incx,
                              Nektar::NekDouble* z, int incz);

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

    template LIB_UTILITIES_EXPORT  void Gathr(int n, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);


    /// \brief Gather vector z[i] = sign[i]*x[y[i]]
    template<class T>  void Gathr(int n, const T *sign, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *z++ = *(sign++) * (*(x + *y++));
        }
        return;
    }

    template LIB_UTILITIES_EXPORT  void Gathr(int n, const Nektar::NekDouble *sign, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);

    /// \brief Scatter vector z[y[i]] = x[i]
    template<class T>  void Scatr(int n, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *(z + *(y++)) = *(x++);
        }
    }

    template LIB_UTILITIES_EXPORT  void Scatr(int n, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);

    /// \brief Scatter vector z[y[i]] = sign[i]*x[i]
    template<class T>  void Scatr(int n, const T *sign, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            if(*sign)
            {
                *(z + *(y++)) = *(sign++) * (*(x++));
            }
            else
            {
                x++;
                y++;
                sign++;
            }
        }
    }

    template LIB_UTILITIES_EXPORT  void Scatr(int n, const Nektar::NekDouble *sign, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);


    /// \brief Assemble z[y[i]] += x[i]; z should be zero'd first
    template<class T>  void Assmb(int n, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *(z + *(y++)) += *(x++);
        }
    }

    template LIB_UTILITIES_EXPORT  void Assmb(int n, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);

    /// \brief Assemble z[y[i]] += sign[i]*x[i]; z should be zero'd first
    template<class T>  void Assmb(int n, const T *sign, const T *x, const int *y,
                  T *z)
    {
        while (n--)
        {
            *(z + *(y++)) += *(sign++) * (*(x++));
        }
    }

    template LIB_UTILITIES_EXPORT  void Assmb(int n, const Nektar::NekDouble *sign, const Nektar::NekDouble *x, const int *y,
              Nektar::NekDouble *z);

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

    template LIB_UTILITIES_EXPORT  Nektar::NekDouble Vsum( int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT  int Vsum( int n, const int *x, const int incx);

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

    template LIB_UTILITIES_EXPORT  int Imax( int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT  int Imax( int n, const int *x, const int incx);

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

    template LIB_UTILITIES_EXPORT  Nektar::NekDouble Vmax( int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT  int Vmax( int n, const int *x, const int incx);

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

    template LIB_UTILITIES_EXPORT int Iamax( int n, const Nektar::NekDouble *x, const int incx);

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

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Vamax( int n, const Nektar::NekDouble *x, const int incx);


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

    template LIB_UTILITIES_EXPORT int Imin( int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT  int Imin( int n, const int *x, const int incx);

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

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Vmin( int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT int Vmin( int n, const int *x, const int incx);

   /// \brief Return number of NaN elements of x
    template<class T>  int Nnan(int n, const T *x, const int incx)
    {

        int nNan = 0;

        while (n--)
        {
            if (*x != *x)
            {
                nNan++;
            }
            x += incx;
        }

        return nNan;
    }

    template LIB_UTILITIES_EXPORT int Nnan(int n, const Nektar::NekDouble *x, const int incx);
    template LIB_UTILITIES_EXPORT int Nnan(int n, const float *x, const int incx);
    template LIB_UTILITIES_EXPORT int Nnan(int n, const int *x, const int incx);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Dot(     int n,
                                 const T   *w,
                                 const T   *x)
    {
        T sum = 0;

        while( n-- )
        {
            sum += (*w) * (*x);
            ++w;
            ++x;
        }
        return sum;
    }

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Dot(     int n,
                             const Nektar::NekDouble   *w,
                             const Nektar::NekDouble   *x);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Dot(     int n,
                                 const T   *w, const int incw,
                                 const T   *x, const int incx)
    {
        T sum = 0;

        while( n-- )
        {
            sum += (*w) * (*x);
            w += incw;
            x += incx;
        }
        return sum;
    }

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Dot(     int n,
                             const Nektar::NekDouble   *w, const int incw,
                             const Nektar::NekDouble   *x, const int incx);

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Dot2(    int n,
                                 const T   *w,
                                 const T   *x,
                                 const int *y)
    {
        T sum = 0;

        while( n-- )
        {
            sum += (*y == 1 ? (*w) * (*x) : 0 );
            ++w;
            ++x;
            ++y;
        }
        return sum;
    }

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Dot2(    int n,
                             const Nektar::NekDouble   *w,
                             const Nektar::NekDouble   *x,
                             const int *y);


    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Dot2(    int n,
                                 const T   *w, const int incw,
                                 const T   *x, const int incx,
                                 const int *y, const int incy)
    {
        T sum = 0;

        while( n-- )
        {
            sum += (*y == 1 ? (*w) * (*x) : 0.0 );
            w += incw;
            x += incx;
            y += incy;
        }
        return sum;
    }

    template LIB_UTILITIES_EXPORT Nektar::NekDouble Dot2(    int n,
                             const Nektar::NekDouble   *w, const int incw,
                             const Nektar::NekDouble   *x, const int incx,
                             const int *y, const int incy);

/*
    // \brief copy one int vector to another
    void Vcopy(int n, const int *x, const int incx, int *y,
                             const int incy)
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

    // \brief copy one unsigned int vector to another
    void Vcopy(int n, const unsigned int *x, const int incx, unsigned int *y,
                             const int incy)
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
    void Vcopy(int n, const double *x, const int incx, double *y,
                             const int incy)
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
*/

    // \brief copy one vector to another
    template<typename T>
    void Vcopy(int n, const T *x, const int incx,
                            T *y, const int incy)
    {
        if( incx ==1 && incy == 1)
        {
            memcpy(y,x,n*sizeof(T));
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

    template  LIB_UTILITIES_EXPORT void  Vcopy( int n, const int *x, const int incx, int *y, const int incy);
    template  LIB_UTILITIES_EXPORT void  Vcopy( int n, const unsigned int *x, const int incx, unsigned int *y, const int incy);
    template  LIB_UTILITIES_EXPORT void  Vcopy( int n, const Nektar::NekDouble *x, const int incx, Nektar::NekDouble *y, const int incy);


    // \brief reverse the ordering of  vector to another
    template<class T>  void  Reverse( int n, const T *x, const int incx, T *y, const int incy)
    {
        int i;
        T store;

        // Perform element by element swaps in case x and y reference the same
        // array.
        int nloop = n/2;

        // copy value in case of n is odd number
        y[nloop] = x[nloop];

        const T* x_end = x + (n-1)*incx;
        T*       y_end = y + (n-1)*incy;
        for (i = 0; i < nloop; ++i) {
            store  = *x_end;
            *y_end = *x;
            *y     = store;
            x     += incx;
            y     += incy;
            x_end -= incx;
            y_end -= incy;
        }
    }

    template  LIB_UTILITIES_EXPORT void  Reverse( int n, const Nektar::NekDouble *x, const int incx, Nektar::NekDouble *y, const int incy);
}
