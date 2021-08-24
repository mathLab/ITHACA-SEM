///////////////////////////////////////////////////////////////////////////////
//
// File: RawType.hpp
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
// Description: Template metafunction which removes pointers, const, and 
// volatile modifiers to a type.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASICUTILS_RAW_TYPE_HPP
#define NEKTAR_LIB_UTILITIES_BASICUTILS_RAW_TYPE_HPP

#include <type_traits>
#include <memory>

namespace Nektar
{
    template<typename T>
    struct RawType
    {
        typedef typename
                  std::decay<typename std::remove_pointer<T>::type >::type type;
    };

    template<typename T>
    struct RawType<std::shared_ptr<T>>
    {
        typedef typename RawType<T>::type type;
    };

    template<typename T>
    struct RawType<const std::shared_ptr<T>>
    {
        typedef typename RawType<T>::type type;
    };

    template<typename T>
    struct RawType<volatile std::shared_ptr<T>>
    {
        typedef typename RawType<T>::type type;
    };

    template<typename T>
    struct RawType<const volatile std::shared_ptr<T>>
    {
        typedef typename RawType<T>::type type;
    };

    template<typename T>
    using RawType_t = typename RawType<T>::type;
}

#endif //NEKTAR_LIB_UTILITIES_BASICUTILS_RAW_TYPE_HPP
