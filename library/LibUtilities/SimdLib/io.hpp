///////////////////////////////////////////////////////////////////////////////
//
// File: io.hpp
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_LIBUTILITES_SIMDLIB_IO_HPP
#define NEKTAR_LIB_LIBUTILITES_SIMDLIB_IO_HPP

#include "traits.hpp"
#include "tinysimd.hpp"
#include <ostream>
#include <memory>

namespace tinysimd
{

template <class T, typename = typename std::enable_if
        <tinysimd::is_vector<T>::value>::type
    >
std::ostream& operator<<(std::ostream& os, const T& avec)
{
    // Note the type cast to 'unsigned int' is only necessary to
    // overcome a bug in Centos 7 gcc 4.8.5 package
    alignas((unsigned int) T::alignment) typename T::scalarArray tmp;
    avec.store(tmp);
    for (unsigned short i = 0; i < T::width; ++i)
    {
        os << tmp[i] << '\t';
    }
    return os;
}

} // namespace tinysimd
#endif
