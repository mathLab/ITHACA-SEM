///////////////////////////////////////////////////////////////////////////////
//
// File: Deprecated.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Defines a macro for deprecated functions.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_DEPRECATED_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_DEPRECATED_HPP

/*
 * Defines a deprecated macro. Note that gcc 4.5+ and it seems virtually all
 * clang versions support supplying a macro, as does MSVC from VS 2015 onwards.
 */
#if defined(__GNUC__) || defined(__clang__)
#  if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR >= 5)) || defined(__clang__)
#    define DEPRECATED(since, alt) \
    __attribute__ ((deprecated(\
                        "since version " #since "; use '" #alt "' instead")))
#  else
#    define DEPRECATED(since, alt) __attribute__ ((deprecated))
#  endif
#elif defined(_MSC_VER) && _MSC_VER >= 1900
#  define DEPRECATED(since, alt) \
    __declspec(deprecated, "since version " #since "; use '" #alt "' instead")
#else
#  define DEPRECATED(since, alt)
#endif

#endif
