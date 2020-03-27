///////////////////////////////////////////////////////////////////////////////
//
// File DVmath.hpp
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
// Description: Collection of templated functions for distributed vector
//              mathematics
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORDISTMATH_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORDISTMATH_HPP

#include <LibUtilities/Communication/Comm.h>

namespace VDmath
{
    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Ddot2(   Nektar::LibUtilities::CommSharedPtr& pComm,
                                 int n,
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
        pComm->AllReduce(sum, Nektar::LibUtilities::ReduceSum);
        return sum;
    }

    /// \brief  vvtvp (vector times vector times vector): z = w*x*y
    template<class T> T Ddot2(   Nektar::LibUtilities::CommSharedPtr& pComm,
                                 int n,
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
        pComm->AllReduce(sum, Nektar::LibUtilities::ReduceSum);
        return sum;
    }
}
#endif /* DVMATH_HPP_ */
