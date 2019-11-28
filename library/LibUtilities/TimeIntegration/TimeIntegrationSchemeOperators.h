///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegrationOperators.h
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
// Description: Header file of time integration operators class
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

///////////////////////////////////////////////////////////////////////////////

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{
/**
 * @brief Binds a set of functions for use by time integration schemes.
 */
class TimeIntegrationSchemeOperators
{
public:
    typedef const Array<OneD, const Array<OneD, NekDouble>> InArrayType;
    typedef Array<OneD, Array<OneD, NekDouble>> OutArrayType;

    typedef std::function<void(InArrayType &, OutArrayType &, const NekDouble)>
        FunctorType1;
    typedef std::function<void(InArrayType &, OutArrayType &, const NekDouble,
                               const NekDouble)>
        FunctorType2;

    typedef const FunctorType1 &ConstFunctorType1Ref;
    typedef const FunctorType2 &ConstFunctorType2Ref;

    typedef Array<OneD, FunctorType1> FunctorType1Array;
    typedef Array<OneD, FunctorType2> FunctorType2Array;

    TimeIntegrationSchemeOperators(void) : m_functors1(4), m_functors2(1)
    {
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineOdeRhs(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[0] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineOdeExplicitRhs(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[1] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineOdeImplicitRhs(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[2] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineProjection(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors1[3] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3);
    }

    template <typename FuncPointerT, typename ObjectPointerT>
    void DefineImplicitSolve(FuncPointerT func, ObjectPointerT obj)
    {
        m_functors2[0] =
            std::bind(func, obj, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4);
    }

    inline void DoOdeRhs(InArrayType &inarray, OutArrayType &outarray,
                         const NekDouble time) const
    {
        ASSERTL1(m_functors1[0], "OdeRhs should be defined for this time "
                                 "integration scheme");
        m_functors1[0](inarray, outarray, time);
    }

    inline void DoOdeExplicitRhs(InArrayType &inarray, OutArrayType &outarray,
                                 const NekDouble time) const
    {
        ASSERTL1(m_functors1[1], "OdeExplicitRhs should be defined for this "
                                 "time integration scheme");
        m_functors1[1](inarray, outarray, time);
    }

    inline void DoOdeImplicitRhs(InArrayType &inarray, OutArrayType &outarray,
                                 const NekDouble time) const
    {
        ASSERTL1(m_functors1[2], "OdeImplictRhs should be defined for this "
                                 "time integration scheme");
        m_functors1[2](inarray, outarray, time);
    }

    inline void DoProjection(InArrayType &inarray, OutArrayType &outarray,
                             const NekDouble time) const
    {
        ASSERTL1(m_functors1[3], "Projection operation should be defined for "
                                 "this time integration scheme");
        m_functors1[3](inarray, outarray, time);
    }

    inline void DoImplicitSolve(InArrayType &inarray, OutArrayType &outarray,
                                const NekDouble time,
                                const NekDouble lambda) const
    {
        ASSERTL1(m_functors2[0], "ImplicitSolve should be defined for this "
                                 "time integration scheme");
        m_functors2[0](inarray, outarray, time, lambda);
    }

protected:
    FunctorType1Array m_functors1;
    FunctorType2Array m_functors2;

private:

};

}
}

