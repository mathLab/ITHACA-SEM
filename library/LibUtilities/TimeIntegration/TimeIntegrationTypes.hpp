///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationTypes.cpp
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
// Description: implementation of time integration key class
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace LibUtilities
{
// Typedefs double arrays
typedef const Array<OneD, const Array<OneD, Array<OneD, NekDouble>>> ConstTripleArray;
typedef       Array<OneD,       Array<OneD, Array<OneD, NekDouble>>>      TripleArray;
typedef const Array<OneD, const Array<OneD,             NekDouble>>  ConstDoubleArray;
typedef       Array<OneD,       Array<OneD,             NekDouble>>       DoubleArray;
typedef const Array<OneD, const                         NekDouble>   ConstSingleArray;
typedef       Array<OneD,                               NekDouble>        SingleArray;

// Typedefs complex double arrays
typedef       Array<OneD,       Array<OneD, Array<OneD, std::complex<NekDouble>>>>      ComplexTripleArray;
typedef const Array<OneD, const Array<OneD,             std::complex<NekDouble>>>  ConstComplexDoubleArray;
typedef       Array<OneD,       Array<OneD,             std::complex<NekDouble>>>       ComplexDoubleArray;
typedef const Array<OneD, const                         std::complex<NekDouble>>   ConstComplexSingleArray;
typedef       Array<OneD,                               std::complex<NekDouble>>        ComplexSingleArray;

// Functors
typedef std::function<void(ConstDoubleArray &, DoubleArray &,
                           const NekDouble)>
    FunctorType1;
typedef std::function<void(ConstDoubleArray &, DoubleArray &,
                           const NekDouble, const NekDouble)>
    FunctorType2;

// Shared pointers
class TimeIntegrationScheme;

typedef std::shared_ptr<TimeIntegrationScheme>
    TimeIntegrationSchemeSharedPtr;

typedef std::vector<TimeIntegrationSchemeSharedPtr>
    TimeIntegrationSchemeVector;

typedef std::vector<TimeIntegrationSchemeSharedPtr>::iterator
    TimeIntegrationSchemeVectorIter;
//
class FractionalInTimeIntegrationScheme;

typedef std::shared_ptr<FractionalInTimeIntegrationScheme>
    FractionalInTimeIntegrationSchemeSharedPtr;

typedef std::vector<FractionalInTimeIntegrationSchemeSharedPtr>
    FractionalInTimeIntegrationSchemeVector;

typedef std::vector<FractionalInTimeIntegrationSchemeSharedPtr>::iterator
    FractionalInTimeIntegrationSchemeVectorIter;
//
class TimeIntegrationSchemeGLM;
  
typedef std::shared_ptr<TimeIntegrationSchemeGLM>
    TimeIntegrationSchemeGLMSharedPtr;

typedef std::vector<TimeIntegrationSchemeGLMSharedPtr>
    TimeIntegrationSchemeGLMVector;

typedef std::vector<TimeIntegrationSchemeGLMSharedPtr>::iterator
    TimeIntegrationSchemeGLMVectorIter;
//
class TimeIntegrationAlgorithmGLM;

typedef std::shared_ptr<TimeIntegrationAlgorithmGLM>
    TimeIntegrationAlgorithmGLMSharedPtr;

typedef std::vector<TimeIntegrationAlgorithmGLMSharedPtr>
    TimeIntegrationAlgorithmGLMVector;

typedef std::vector<TimeIntegrationAlgorithmGLMSharedPtr>::iterator
    TimeIntegrationAlgorithmGLMVectorIter;  
//
class TimeIntegrationSolutionGLM;

typedef std::shared_ptr<TimeIntegrationSolutionGLM>
    TimeIntegrationSolutionGLMSharedPtr;

typedef std::vector<TimeIntegrationSolutionGLMSharedPtr>
    TimeIntegrationSolutionGLMVector;

typedef std::vector<TimeIntegrationSolutionGLMSharedPtr>::iterator
    TimeIntegrationSolutionGLMVectorIter;

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

}
}
