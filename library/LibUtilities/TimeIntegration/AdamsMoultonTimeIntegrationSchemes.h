///////////////////////////////////////////////////////////////////////////////
//
// File: AdamsMoultonTimeIntegrationSchemes.h
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
// Description: Combined header file for all Adams Moulton based time
// integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/EulerTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Adams Bashforth Order 2

class AdamsMoultonOrder1TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    AdamsMoultonOrder1TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        AdamsMoultonOrder1TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~AdamsMoultonOrder1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsMoultonOrder1TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("AdamsMoultonOrder1");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        BackwardEulerTimeIntegrationScheme::SetupSchemeData(phase);
    }

}; // end class AdamsMoultonOrder1TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// Adams Bashforth Order 2

class AdamsMoultonOrder2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    AdamsMoultonOrder2TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(2);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[1] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        ForwardEulerTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
        AdamsMoultonOrder2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[1]);
    }

    virtual ~AdamsMoultonOrder2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            AdamsMoultonOrder2TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("AdamsMoultonOrder2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eDiagonallyImplicit;

        phase->m_numsteps  = 2;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.5);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = 0.5;
        phase->m_B[0][1][0] = 1.0;

        phase->m_U[0][0] = 1.0;
        phase->m_U[0][1] = 0.5;

        phase->m_V[0][0] = 1.0;
        phase->m_V[0][1] = 0.5;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 1;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 0;

        phase->m_firstStageEqualsOldSolution =
            phase->CheckIfFirstStageEqualsOldSolution(phase->m_A, phase->m_B,
                                                      phase->m_U, phase->m_V);
        phase->m_lastStageEqualsNewSolution =
            phase->CheckIfLastStageEqualsNewSolution(phase->m_A, phase->m_B,
                                                     phase->m_U, phase->m_V);

        ASSERTL1(phase->VerifyIntegrationSchemeType(phase->m_schemeType,
                                                    phase->m_A, phase->m_B,
                                                    phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class AdamsMoultonOrder2TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
