///////////////////////////////////////////////////////////////////////////////
//
// File: LibSMV.hpp
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
// Description: wrapper of functions around SMV routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SMV_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SMV_HPP

#ifdef NEKTAR_USING_SMV

#include <LibUtilities/LinearAlgebra/TransF77.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <boost/preprocessor/iteration/local.hpp>

#include <smv.h>


namespace Smv
{
    // Translations for using Fortran version of SMV
    extern "C"
    {
        /// \brief  Matrix-vector multiply C = C + A*B where
        ///         A is [m x m], B is [m] and C is [m].
        /// Expected: no matrix transpose, row-major ordering,
        /// unit increments and no matrix views (lda != m|n|k),
        /// type double, no constant factors.
        void F77NAME(smv) (
                 const int& m,
                 const double* a, const double* b, double* c);

        /// \brief  Rank-specific matrix-vector LibSMV multiply
        ///         kernels. Row-major ordering, unit increments,
        ///         type double. May eventually call dgemv
        ///         implementation that LibSMV is linked against.
        #define BOOST_PP_LOCAL_MACRO(n)         \
                       void F77NAME(smv_##n)    \
                             (const double* a,  \
                              const double* b,  \
                                    double* c);
        #define BOOST_PP_LOCAL_LIMITS     (1, LIBSMV_MAX_RANK)
        #include BOOST_PP_LOCAL_ITERATE()
    }


    /// \brief LibSmv matrix-vector multiply: Y = Y + A*X where A is [m x m]-matrix
    template <typename T>
    LIB_UTILITIES_EXPORT   void Smvn (const int& m,
                             const T* a,    const T* x,    T* y);

}
#endif // NEKTAR_USING_SMV
#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SMV_HPP
