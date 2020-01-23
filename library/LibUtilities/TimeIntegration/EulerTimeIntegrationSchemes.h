///////////////////////////////////////////////////////////////////////////////
//
// File: EulerTimeIntegrationSchemes.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Division of Applied Mathematics, Brown University (USA),
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
// Description: Combined header file for all basic Euler based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Backward Euler

class BackwardEulerTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    BackwardEulerTimeIntegrationScheme(std::string variant, unsigned int order,
				       std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        BackwardEulerTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~BackwardEulerTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
						 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          BackwardEulerTimeIntegrationScheme>::AllocateSharedPtr("", 1, freeParams);
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("BackwardEuler");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order = 1;
        phase->m_name = std::string("BackwardEulerOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 1.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 1.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

}; // end class BackwardEulerTimeIntegrator

///////////////////////////////////////////////////////////////////////////////
// Forward Euler

class ForwardEulerTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    ForwardEulerTimeIntegrationScheme(std::string variant, unsigned int order,
				      std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme("", 1, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        ForwardEulerTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~ForwardEulerTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
						 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          ForwardEulerTimeIntegrationScheme>::AllocateSharedPtr("", 1, freeParams);
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("ForwardEuler");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.784;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;
        phase->m_order = 1;
        phase->m_name = std::string("ForwardEulerOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 1.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

}; // end class ForwardEulerTimeIntegrator

} // end namespace LibUtilities
} // end namespace Nektar
