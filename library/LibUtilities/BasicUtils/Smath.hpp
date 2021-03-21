///////////////////////////////////////////////////////////////////////////////
//
// File: Smath.hpp
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
// Description: Collection of templated functions for scalar mathematics
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_SCALARMATH_H
#define NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_SCALARMATH_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <algorithm>
#include <cstdlib>
#include <math.h>
using namespace std;
using namespace Nektar;

namespace Smath
{

/***************** Math routines  ***************/

#if 1 // try removing LIB_UTILITIES_EXPORT 
/// \brief Return the soft max of between two scalars
template <class T> T Smax(const T a, const T b, const T k)
{
    T maxi = std::max(a, b) * k;
    T mini = std::min(a, b) * k;
    T xmax = (maxi + log(1.0 + exp(mini - maxi))) / k;
    return xmax;
}

template NekDouble Smax(const NekDouble a, const NekDouble b,
                        const NekDouble k);
    
template int Smax(const int a, const int b, const int k);
#else
/// \brief Return the soft max of between two scalars
template <class T, typename = typename std::enable_if
          < std::is_floating_point<T>::value ||
            std::is_integral<T> >::type > 
LIB_UTILITIES_EXPORT T Smax(const T a, const T b, const T k)
{
    T maxi = std::max(a, b) * k;
    T mini = std::min(a, b) * k;
    T xmax = (maxi + log(1.0 + exp(mini - maxi))) / k;
    return xmax;
}
#endif
} // namespace Smath
#endif // NEKTAR_LIB_LIBUTILITIES_BASSICUTILS_SCALARMATH_H
