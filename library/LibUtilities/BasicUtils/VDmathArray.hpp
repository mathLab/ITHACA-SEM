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

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORDISTMATHARRAY_HPP
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_VECTORDISTMATHARRAY_HPP

#include <LibUtilities/BasicUtils/VDmath.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>


namespace VDmath
{
    using namespace Nektar;

    /// \brief
    template<class T> T Ddot2(LibUtilities::CommSharedPtr& pComm, int n,
                              const Array<OneD, const T> &w,
                              const Array<OneD, const T> &x,
                              const Array<OneD, const int> &y)
    {
        ASSERTL1(n <= w.size()+w.GetOffset(),"Array out of bounds");
        ASSERTL1(n <= x.size()+x.GetOffset(),"Array out of bounds");
        ASSERTL1(n <= y.size()+y.GetOffset(),"Array out of bounds");

        return Ddot2(pComm, n,&w[0],&x[0],&y[0]);
    }

    /// \brief
    template<class T> T Ddot2(LibUtilities::CommSharedPtr& pComm, int n,
                              const Array<OneD, const T> &w, const int incw,
                              const Array<OneD, const T> &x, const int incx,
                              const Array<OneD, const int> &y, const int incy)
    {
        ASSERTL1(n*incw <= w.size()+w.GetOffset(),"Array out of bounds");
        ASSERTL1(n*incx <= x.size()+x.GetOffset(),"Array out of bounds");
        ASSERTL1(n*incy <= y.size()+y.GetOffset(),"Array out of bounds");

        return Ddot2(pComm, n,&w[0],incw,&x[0],incx,&y[0],incy);
    }

}

#endif /* VDMATHARRAY_HPP_ */
