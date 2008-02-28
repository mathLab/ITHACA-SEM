///////////////////////////////////////////////////////////////////////////////
//
// File VmathArray.hpp
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
// Description: Wrappers around Vmath routines using Array<OneD,T> as arugments
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATHARRAY_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORMATHARRAY_HPP

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp> 

using namespace Nektar;
    namespace Vmath 
    {
    
        /***************** Math routines  ***************/
        /// \brief Fill a vector with a constant value
        template<class T>  void Fill( int n, const T alpha,  Array<OneD, T> &x, const int incx )
        {
            
            ASSERTL1(n*incx <= x.num_elements(),"Out of bounds");
            
            Fill(n,alpha,&x[0],incx);
        }
        
        
        /// \brief Multiply vector z = x*y
        template<class T>  void Vmul( int n, const ConstArray<OneD,T> &x, const int incx, const ConstArray<OneD,T> &y, const int incy,  Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");
            
            Vmul(n,&x[0],incx,&y[0],incy,&z[0],incz);
        }

        /// \brief Scalar multiply  y = alpha*y
        
        template<class T>  void Smul( int n, const T alpha, const ConstArray<OneD,T> &x,  const int incx,  Array<OneD,T>  &y, const int incy)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            
            Smul(n,alpha, &x[0],incx,&y[0],incy);
        }
        
        /// \brief Multiply vector z = x/y
        template<class T>  void Vdiv( int n, const ConstArray<OneD,T> &x, const int incx, Array<OneD,T> &y, const int incy,  Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");
            
            Vdiv(n,&x[0],incx,&y[0],incy,&z[0],incz);
            
        }
        
        /// \brief Scalar multiply  y = alpha/y
        template<class T>  void Sdiv( int n, const T alpha, Array<OneD,T> &x, const int incx,  Array<OneD,T> &y, const int incy)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            
            Sdiv(n,alpha,&x[0],incx,&y[0],incy);
        }
        
        /// \brief Add vector z = x+y
        template<class T>  void Vadd( int n, const ConstArray<OneD,T> &x, const int incx, const ConstArray<OneD,T> &y,  const int incy,  Array<OneD,T> &z, const int incz)
        {

            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");
         
            Vadd(n,&x[0],incx,&y[0],incy,&z[0],incz);
        }
    
        /// \brief Add vector y = alpha + x
        template<class T>  void Sadd( int n, const T alpha, const ConstArray<OneD,T> &x,const int incx, Array<OneD,T> &y, const int incy)
        {

            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");

            Sadd(n,alpha,&x[0],incx,&y[0],incy);
        }
    
        /// \brief Subtract vector z = x-y
        template<class T>  void Vsub( int n, const ConstArray<OneD,T> &x, const int incx, const ConstArray<OneD,T> &y, const int incy,  Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");

            Vsub(n,&x[0],incx,&y[0],incy,&z[0],incz);

        }
    
        /// \brief Zero vector
        template<class T>  void Zero(int n, Array<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");

            Zero(n,&x[0],incx);

        }
        
        /// \brief Negate x = -x
        template<class T>  void Neg( int n, Array<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");

            Neg(n,&x[0],incx);
            
        }
    
    
        /// \brief sqrt y = sqrt(x)
        template<class T> void Vsqrt(int n, const ConstArray<OneD,T> &x, const int incx, Array<OneD,T> &y, const int incy)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");

            Vsqrt(n,&x[0],incx,&y[0],incy);
        }
    
        /// \brief vabs: y = |x|
        template<class T> void Vabs(int n, const ConstArray<OneD,T> &x, const int incx, Array<OneD,T> &y, const int incy)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");

            Vabs(n,&x[0],incx,&y[0],incy);
        }
    
        /********** Triad  routines  ***********************/
        
        /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
        template<class T> void Vvtvp(int n, const ConstArray<OneD,T> &w, const int incw, const ConstArray<OneD,T> &x, const int incx, const ConstArray<OneD,T> &y, const int incy, Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incw <= w.num_elements(),"Array out of bounds");
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");

            Vvtvp(n,&w[0],incw,&x[0],incx,&y[0],incy,&z[0],incz);
        }

        /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
        template<class T> void Svtvp(int n, const T alpha, const ConstArray<OneD,T> &x,  const int incx, const ConstArray<OneD,T> &y, const int incy, Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");
            
            Svtvp(n,alpha,&x[0],incx,&y[0],incy,&z[0],incz);
            
        }
    
        /// \brief vvtvm (vector times vector plus vector): z = w*x - y
        template<class T> void Vvtvm(int n, const ConstArray<OneD,T> &w, const int incw, Array<OneD,T> &x, const int incx, const ConstArray<OneD,T> &y, const int incy,  Array<OneD,T> &z, const int incz)
        {
            ASSERTL1(n*incw <= w.num_elements(),"Array out of bounds");
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n*incz <= z.num_elements(),"Array out of bounds");

            Vvtvm(n,&w[0],incw,&x[0],incx,&y[0],incy,&z[0],incz);
            
        }
    
    
        /************ Misc routine from Veclib (and extras)  ************/
        
        /// \brief Gather vector z[i] = x[y[i]]
        template<class T>  void Gathr(int n, const ConstArray<OneD,T> &x, const ConstArray<OneD,int> &y,  Array<OneD,T> &z)
        {
            
            ASSERTL1(n <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n <= z.num_elements(),"Array out of bounds");

            Gathr(n,&x[0],&y[0],&z[0]);

        }
    
        /// \brief Scatter vector z[y[i]] = x[i]
        template<class T>  void Scatr(int n, const ConstArray<OneD,T> &x, const ConstArray<OneD,int> &y,  Array<OneD,T> &z)
        {
            ASSERTL1(n <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n <= z.num_elements(),"Array out of bounds");

            Scatr(n,&x[0],&y[0],&z[0]);
        }
    
    
        /// \brief Assemble z[y[i]] += x[i]; z should be zero'd first
        template<class T>  void Assmb(int n, const ConstArray<OneD,T> &x, const ConstArray<OneD,int> &y, Array<OneD,T> &z)
        {
            ASSERTL1(n <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n <= y.num_elements(),"Array out of bounds");
            ASSERTL1(n <= z.num_elements(),"Array out of bounds");

            Assmb(n,&x[0],&y[0],&z[0]);
        }
    
    
        /************* Reduction routines  *****************/
        
        /// \brief Subtract return sum(x)
        template<class T>  T Vsum( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");

            return Vsum(n,&x[0],incx);
        }
    
    
        /// \brief Return the index of the maximum element in x
        template<class T>  int Imax( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
    
            return Imax(n,&x[0],incx);
        }
    
        /// \brief Return the maximum element in x -- called vmax to avoid
        /// conflict with max
        template<class T>  T Vmax( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
    
            return Vmax(n,&x[0],incx);
        }
    
        /// \brief Return the index of the maximum absolute element in x
        template<class T>  int Iamax( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
    
            return Iamax(n,&x[0],incx);

        }
    
        /// \brief Return the maximum absolute element in x
        /// called vamax to avoid conflict with max
        template<class T>  T Vamax( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
    
            return Vamax(n,&x[0],incx);            
        }
    
    
        /// \brief Return the index of the minimum element in x
        template<class T>  int Imin( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            
            return Imin(n,&x[0],incx);
        }
    
    
        /// \brief Return the minimum element in x - called vmin to avoid
        /// conflict with min
        template<class T>  T Vmin( int n, const ConstArray<OneD,T> &x, const int incx)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            
            return Vmin(n,&x[0],incx);
        }

        /********** Memory routines  ***********************/
        
        template<class T> void Vcopy(int n, const ConstArray<OneD, T> &x, int incx, Array<OneD,T> &y, int const incy)
        {
            ASSERTL1(n*incx <= x.num_elements(),"Array out of bounds");
            ASSERTL1(n*incy <= y.num_elements(),"Array out of bounds");

            Vcopy(n,&x[0],incx,&y[0],incy);
        }
        
    }
#endif //VECTORMATHARRAY_HPP

/***
$Log:$

**/
