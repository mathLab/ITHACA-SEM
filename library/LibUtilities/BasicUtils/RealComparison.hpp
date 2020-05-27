///////////////////////////////////////////////////////////////////////////////
//
// File RealComparison.hpp
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
// Description: simple routines to compare 2 real
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_REALCOMPARISON_H
#define NEKTAR_LIB_UTILITIES_REALCOMPARISON_H

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <limits>
#include <type_traits>
#include <cmath>

namespace Nektar
{
namespace LibUtilities
{
/// compare reals of same type with relative tolerance
template
<
    class T1, class T2,
    class = typename std::enable_if
    <
        std::is_floating_point
        <
            typename std::remove_cv<
                typename std::remove_reference<T1>::type>::type
        >::value &&
        std::is_same
        <
            typename std::remove_cv<
                typename std::remove_reference<T1>::type>::type,
            typename std::remove_cv<
                typename std::remove_reference<T2>::type>::type
        >::value
    >::type
>
inline bool IsRealEqual(T1&& lhs, T2&& rhs,
    const unsigned int factor = NekConstants::kNekFloatCompFact)
{
    // Check precondition in debug mode
    ASSERTL1(factor >= 1, "real comparison factor needs to be >= 1");
    // Get base type
    typedef typename std::remove_reference<T1>::type Tbase;
    // Tolerance
    Tbase tol = factor * std::numeric_limits<Tbase>::epsilon();
    // Distance
    Tbase dist = std::abs(lhs-rhs);
    // Reference
    Tbase ref = std::max(std::abs(lhs), std::abs(rhs));
    return dist < ref * tol || dist == 0;
}

}
}

#endif
