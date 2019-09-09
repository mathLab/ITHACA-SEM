///////////////////////////////////////////////////////////////////////////////
//
// File: TestRawType.cpp
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
// Description: Test the RawType metafunction.  These are compile tests, so if
// this file compiles successfully, then all tests have passed.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/RawType.hpp>
#include <type_traits>

namespace Nektar
{
static_assert(std::is_same<double, RawType<double>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const double>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile double>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile double>::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<double*>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const double*>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile double*>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile double*>::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<double* const>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const double* const>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile double* const>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile double* const>::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<double* volatile>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const double* volatile>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile double* volatile >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile double* volatile>::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<double* const volatile>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const double* const volatile>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile double* const volatile>::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile double* const volatile>::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<std::shared_ptr<double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const std::shared_ptr<double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile std::shared_ptr<double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile std::shared_ptr<double> >::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<std::shared_ptr<const double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const std::shared_ptr<const double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile std::shared_ptr<const double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile std::shared_ptr<const double> >::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<std::shared_ptr<volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const std::shared_ptr<volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile std::shared_ptr<volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile std::shared_ptr<volatile double> >::type>::value,
              "RawType error");

static_assert(std::is_same<double, RawType<std::shared_ptr<const volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const std::shared_ptr<const volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<volatile std::shared_ptr<const volatile double> >::type>::value,
              "RawType error");
static_assert(std::is_same<double, RawType<const volatile std::shared_ptr<const volatile double> >::type>::value,
              "RawType error");

}
