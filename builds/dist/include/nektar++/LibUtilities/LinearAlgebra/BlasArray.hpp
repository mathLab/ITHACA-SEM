///////////////////////////////////////////////////////////////////////////////
//
// File BlasArray.hpp
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
// using Array's as calling arguments
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP

#include <LibUtilities/LinearAlgebra/TransF77.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

// Translations for using Fortran version of blas
namespace Blas
{

    extern "C"
    {
        // -- BLAS Level 1:
        void   F77NAME(dcopy) (const int& n, const double *x, const int& incx,
                   double *y, const int& incy);
        void   F77NAME(daxpy) (const int& n, const double& alpha, const double *x,
                   const int& incx, const double *y, const int& incy);
        double F77NAME(ddot)  (const int& n, const double *x,  const int& incx,
                   const double *y, const int& incy);
    }

#ifdef NEKTAR_USING_BLAS
    static inline void Dcopy (const int& n, const Nektar::Array <Nektar::OneD, const double> &x, const int& incx,  Nektar::Array<Nektar::OneD,double> &y, const int& incy)
    {
        ASSERTL1(static_cast<unsigned int>(n*incx) <= x.num_elements()+x.GetOffset(),"Array out of bounds");
        ASSERTL1(static_cast<unsigned int>(n*incy) <= y.num_elements()+y.GetOffset(),"Array out of bounds");
        
        F77NAME(dcopy)(n,&x[0],incx,&y[0],incy);
    }

    /// \brief  BLAS level 1: y = alpha \a x plus \a y
    static inline void Daxpy (const int& n, const double& alpha, const Nektar::Array <Nektar::OneD,const double> &x, const int& incx,  Nektar::Array<Nektar::OneD,double> &y, const int& incy)
    {
        ASSERTL1(static_cast<unsigned int>(n*incx) <= x.num_elements()+x.GetOffset(),"Array out of bounds");
        ASSERTL1(static_cast<unsigned int>(n*incy) <= y.num_elements()+y.GetOffset(),"Array out of bounds");
        
        F77NAME(daxpy)(n,alpha,&x[0],incx,&y[0],incy);
    }

    /// \brief BLAS level 1: output =   \f$ x^T  y \f$
    static inline double Ddot (const int& n, const Nektar::Array <Nektar::OneD,const double> &x, const int& incx,
              const Nektar::Array <Nektar::OneD,const double> &y, const int& incy)
    {
        return F77NAME(ddot)(n,&x[0],incx,&y[0],incy);
    }
#endif // NEKTAR_USING BLAS
}
#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_BLASARRAY_HPP

/***
$Log: BlasArray.hpp,v $
Revision 1.5  2008/05/10 18:27:32  sherwin
Modifications necessary for QuadExp Unified DG Solver

Revision 1.4  2008/04/30 02:55:51  bnelson
Fixed gcc compiler warning.

Revision 1.3  2008/04/06 05:55:11  bnelson
Changed ConstArray to Array<const>

Revision 1.2  2008/03/12 15:22:45  pvos
Clean up of the code

Revision 1.1  2008/02/28 09:57:08  sherwin
Added array version of some routines

**/
