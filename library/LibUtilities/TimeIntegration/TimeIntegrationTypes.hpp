///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationTypes.hpp
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
// Description: implementation of time integration types that are
// common to all schemes.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_TYPES
#define NEKTAR_LIB_UTILITIES_TIME_INTEGRATION_TIME_INTEGRATION_TYPES

#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace LibUtilities
{

// Typedefs double arrays
template <class T> using AT = Array<OneD, T>;

// clang-format off
typedef const AT<const AT<AT<NekDouble>>> ConstTripleArray;
typedef       AT<      AT<AT<NekDouble>>>      TripleArray;
typedef const AT<const AT<   NekDouble>>  ConstDoubleArray;
typedef       AT<      AT<   NekDouble>>       DoubleArray;
typedef const AT<      const NekDouble>   ConstSingleArray;
typedef       AT<            NekDouble>        SingleArray;

typedef       AT<      AT<AT<std::complex<NekDouble>>>>      ComplexTripleArray;
typedef const AT<const AT<   std::complex<NekDouble>>>  ConstComplexDoubleArray;
typedef       AT<      AT<   std::complex<NekDouble>>>       ComplexDoubleArray;
typedef const AT<      const std::complex<NekDouble>>   ConstComplexSingleArray;
typedef       AT<            std::complex<NekDouble>>        ComplexSingleArray;
// clang-format on

// Functors
typedef std::function<void(ConstDoubleArray &, DoubleArray &, const NekDouble)>
    FunctorType1;
typedef std::function<void(ConstDoubleArray &, DoubleArray &, const NekDouble,
                           const NekDouble)>
    FunctorType2;

// Shared pointers
class TimeIntegrationScheme;

typedef std::shared_ptr<TimeIntegrationScheme> TimeIntegrationSchemeSharedPtr;

typedef std::vector<TimeIntegrationSchemeSharedPtr> TimeIntegrationSchemeVector;

//
class FractionalInTimeIntegrationScheme;

typedef std::shared_ptr<FractionalInTimeIntegrationScheme>
    FractionalInTimeIntegrationSchemeSharedPtr;

typedef std::vector<FractionalInTimeIntegrationSchemeSharedPtr>
    FractionalInTimeIntegrationSchemeVector;

//
class TimeIntegrationSchemeGLM;

typedef std::shared_ptr<TimeIntegrationSchemeGLM>
    TimeIntegrationSchemeGLMSharedPtr;

typedef std::vector<TimeIntegrationSchemeGLMSharedPtr>
    TimeIntegrationSchemeGLMVector;

//
class TimeIntegrationAlgorithmGLM;

typedef std::shared_ptr<TimeIntegrationAlgorithmGLM>
    TimeIntegrationAlgorithmGLMSharedPtr;

typedef std::vector<TimeIntegrationAlgorithmGLMSharedPtr>
    TimeIntegrationAlgorithmGLMVector;

//
class TimeIntegrationSolutionGLM;

typedef std::shared_ptr<TimeIntegrationSolutionGLM>
    TimeIntegrationSolutionGLMSharedPtr;

typedef std::vector<TimeIntegrationSolutionGLMSharedPtr>
    TimeIntegrationSolutionGLMVector;

//
enum TimeIntegrationSchemeType
{
    eNoTimeIntegrationSchemeType,
    eExplicit,           //!< Formally explicit scheme
    eDiagonallyImplicit, //!< Diagonally implicit scheme (e.g. the DIRK schemes)
    eIMEX,               //!< Implicit Explicit General Linear Method
    eImplicit,           //!< Fully implicit scheme
    eExponential,        //!< Exponential scheme
    eFractionalInTime,   //!< Fractional in Time scheme
};

} // end namespace LibUtilities
} // end namespace Nektar

#endif
