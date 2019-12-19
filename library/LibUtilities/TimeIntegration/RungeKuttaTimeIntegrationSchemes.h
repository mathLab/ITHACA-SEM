///////////////////////////////////////////////////////////////////////////////
//
// File: RungeKuttaTimeIntegrationSchemes.h
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
// Description: Combined header file for all Runge Kutta based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

namespace Nektar
{
namespace LibUtilities
{

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta 2

class RungeKutta2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    RungeKutta2TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKutta2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~RungeKutta2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            RungeKutta2TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 2;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 1.0);

        phase->m_A[0][1][0] = 0.5;
        phase->m_B[0][0][1] = 1.0;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

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

}; // end class RungeKutta2TimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta 2 Improved Euler

class RungeKutta2_ImprovedEulerTimeIntegrationScheme
    : public TimeIntegrationScheme
{
public:
    RungeKutta2_ImprovedEulerTimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKutta2_ImprovedEulerTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~RungeKutta2_ImprovedEulerTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<RungeKutta2_ImprovedEulerTimeIntegrationScheme>::
                AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta2_ImprovedEuler");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 2;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 1.0);

        phase->m_A[0][1][0] = 1.0;

        phase->m_B[0][0][0] = 0.5;
        phase->m_B[0][0][1] = 0.5;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->m_firstStageEqualsOldSolution =
            phase->CheckIfFirstStageEqualsOldSolution(phase->m_A, phase->m_B,
                                                      phase->m_U, phase->m_V);
        phase->m_lastStageEqualsNewSolution =
            phase->CheckIfLastStageEqualsNewSolution(phase->m_A, phase->m_B,
                                                     phase->m_U, phase->m_V);

        ASSERTL1(phase->VerifyIntegrationSchemeType(phase->m_schemeType,
                                                    phase->m_A, phase->m_B,
                                                    phase->m_U, phase->m_V),
                 "Time integration phase coefficients do not match its type");
    }

}; // end class RungeKutta2_ImprovedEulerTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta 2 SSP

class RungeKutta2_SSPTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    RungeKutta2_SSPTimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKutta2_SSPTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~RungeKutta2_SSPTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            RungeKutta2_SSPTimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta2_SSP");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        RungeKutta2_ImprovedEulerTimeIntegrationScheme::SetupSchemeData(phase);
    }

}; // end class RungeKutta2_SSPTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta 3 SSP

class RungeKutta3_SSPTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    RungeKutta3_SSPTimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKutta3_SSPTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~RungeKutta3_SSPTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            RungeKutta3_SSPTimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta3_SSP");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 3;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 1.0);

        phase->m_A[0][1][0] = 1.0;
        phase->m_A[0][2][0] = 0.25;
        phase->m_A[0][2][1] = 0.25;

        phase->m_B[0][0][0] = 1.0 / 6.0;
        phase->m_B[0][0][1] = 1.0 / 6.0;
        phase->m_B[0][0][2] = 2.0 / 3.0;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->m_firstStageEqualsOldSolution =
            phase->CheckIfFirstStageEqualsOldSolution(phase->m_A, phase->m_B,
                                                      phase->m_U, phase->m_V);
        phase->m_lastStageEqualsNewSolution =
            phase->CheckIfLastStageEqualsNewSolution(phase->m_A, phase->m_B,
                                                     phase->m_U, phase->m_V);

        ASSERTL1(phase->VerifyIntegrationSchemeType(phase->m_schemeType,
                                                    phase->m_A, phase->m_B,
                                                    phase->m_U, phase->m_V),
                 "Time integration phase coefficients do not match its type");
    }

}; // end class RungeKutta3_SSPTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta 5

class RungeKutta5TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    RungeKutta5TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKutta5TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~RungeKutta5TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            RungeKutta5TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta5");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.784;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 6;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 1.0);

        phase->m_A[0][1][0] = 1.0 / 4.0;
        phase->m_A[0][2][0] = 1.0 / 8.0;
        phase->m_A[0][2][1] = 1.0 / 8.0;
        phase->m_A[0][3][2] = 1.0 / 2.0;
        phase->m_A[0][4][0] = 3.0 / 16.0;
        phase->m_A[0][4][1] = -3.0 / 8.0;
        phase->m_A[0][4][2] = 3.0 / 8.0;
        phase->m_A[0][4][3] = 9.0 / 16.0;
        phase->m_A[0][5][0] = -3.0 / 7.0;
        phase->m_A[0][5][1] = 8.0 / 7.0;
        phase->m_A[0][5][2] = 6.0 / 7.0;
        phase->m_A[0][5][3] = -12.0 / 7.0;
        phase->m_A[0][5][4] = 8.0 / 7.0;

        phase->m_B[0][0][0] = 7.0 / 90.0;
        phase->m_B[0][0][1] = 32.0 / 90.0;
        phase->m_B[0][0][2] = 12.0 / 90.0;
        phase->m_B[0][0][3] = 32.0 / 90.0;
        phase->m_B[0][0][4] = 7.0 / 90.0;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->m_firstStageEqualsOldSolution =
            phase->CheckIfFirstStageEqualsOldSolution(phase->m_A, phase->m_B,
                                                      phase->m_U, phase->m_V);
        phase->m_lastStageEqualsNewSolution =
            phase->CheckIfLastStageEqualsNewSolution(phase->m_A, phase->m_B,
                                                     phase->m_U, phase->m_V);

        ASSERTL1(phase->VerifyIntegrationSchemeType(phase->m_schemeType,
                                                    phase->m_A, phase->m_B,
                                                    phase->m_U, phase->m_V),
                 "Time integration phase coefficients do not match its type");
    }

}; // end class RungeKutta5TimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Classic RungeKutta 4

class ClassicalRungeKutta4TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    ClassicalRungeKutta4TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        ClassicalRungeKutta4TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~ClassicalRungeKutta4TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            ClassicalRungeKutta4TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("ClassicalRungeKutta4");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 2.784;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExplicit;

        phase->m_numsteps  = 1;
        phase->m_numstages = 4;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 1.0);

        phase->m_A[0][1][0] = 0.5;
        phase->m_A[0][2][1] = 0.5;
        phase->m_A[0][3][2] = 1.0;

        phase->m_B[0][0][0] = 1.0 / 6.0;
        phase->m_B[0][0][1] = 1.0 / 3.0;
        phase->m_B[0][0][2] = 1.0 / 3.0;
        phase->m_B[0][0][3] = 1.0 / 6.0;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->m_firstStageEqualsOldSolution =
            phase->CheckIfFirstStageEqualsOldSolution(phase->m_A, phase->m_B,
                                                      phase->m_U, phase->m_V);
        phase->m_lastStageEqualsNewSolution =
            phase->CheckIfLastStageEqualsNewSolution(phase->m_A, phase->m_B,
                                                     phase->m_U, phase->m_V);

        ASSERTL1(phase->VerifyIntegrationSchemeType(phase->m_schemeType,
                                                    phase->m_A, phase->m_B,
                                                    phase->m_U, phase->m_V),
                 "Time integration phase coefficients do not match its type");
    }

}; // end class ClassicalRungeKutta4TimeIntegrator

} // end namespace LibUtilities
} // end namespace Nektar
