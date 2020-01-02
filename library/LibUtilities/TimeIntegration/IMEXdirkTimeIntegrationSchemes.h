///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXdirkTimeIntegrationSchemes.h
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
// Description: Combined header file for all IMEX Dirk based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#define LUE LIB_UTILITIES_EXPORT

namespace Nektar
{
namespace LibUtilities
{

//  IMEXDirk-(s, sigma, p), where s is the number of implicit stage
//  schemes, sigma is the number of explicit stage scheme and p is the
//  combined order of the scheme.

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 1 1 : Forward - Backward Euler IMEX

class IMEXdirk_1_1_1TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_1_1_1TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(1, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_1_1_1TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);

        // for( unsigned int n=0; n<m_order; ++n )
        //   std::cout << m_integration_phases[n] << std::endl;
    }

    virtual ~IMEXdirk_1_1_1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          IMEXdirk_1_1_1TimeIntegrationScheme>::AllocateSharedPtr(1, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_1_1_1");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 1;
        phase->m_name = std::string("IMEXdirk_1_1_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 2;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  0.0);

        phase->m_A[0][0][0] = 1.0;

        phase->m_B[0][0][0] = 1.0;
        phase->m_B[1][1][0] = 1.0;

        phase->m_U[0][0]    = 1.0;
        phase->m_U[0][1]    = 1.0;

        phase->m_V[0][0]    = 1.0;
        phase->m_V[0][1]    = 1.0;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_1_1_1TimeIntegrator

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 2 1 : Forward - Backward Euler IMEX

class IMEXdirk_1_2_1TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_1_2_1TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(1, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_1_2_1TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);

        // for( unsigned int n=0; n<m_order; ++n )
        //   std::cout << m_integration_phases[n] << std::endl;
    }

    virtual ~IMEXdirk_1_2_1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_1_2_1TimeIntegrationScheme>::AllocateSharedPtr(1, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_1_2_1");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 1;
        phase->m_name = std::string("IMEXdirk_1_2_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 2;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_A[0][1][1] = 1.0;
        phase->m_B[0][0][1] = 1.0;

        phase->m_A[1][1][0] = 1.0;
        phase->m_B[1][0][1] = 1.0;

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

}; // end class IMEXdirk_1_2_1TimeIntegrator

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 2 2 : Implict-Explicit Midpoint IMEX

class IMEXdirk_1_2_2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_1_2_2TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(1, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_1_2_2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_1_2_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_1_2_2TimeIntegrationScheme>::AllocateSharedPtr(1, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_1_2_2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 1;
        phase->m_name = std::string("IMEXdirk_1_2_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 2;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_A[0][1][1] = 1.0 / 2.0;
        phase->m_B[0][0][1] = 1.0;

        phase->m_A[1][1][0] = 1.0 / 2.0;
        phase->m_B[1][0][1] = 1.0;

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

}; // end class IMEXdirk_1_2_2TimeIntegrator

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 2 2 : L Stable, two implict stages, two explict stages, second order IMEX

class IMEXdirk_2_2_2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_2_2_2TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(2, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_2_2_2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_2_2_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_2_2_2TimeIntegrationScheme>::AllocateSharedPtr(2, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_2_2_2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 2;
        phase->m_name = std::string("IMEXdirk_2_2_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 3;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        NekDouble glambda = 0.2928932188134524756;
        NekDouble gdelta  = -0.7071067811865475244;

        phase->m_A[0][1][1] = glambda;
        phase->m_A[0][2][1] = 1.0 - glambda;
        phase->m_A[0][2][2] = glambda;

        phase->m_B[0][0][1] = 1.0 - glambda;
        phase->m_B[0][0][2] = glambda;

        phase->m_A[1][1][0] = glambda;
        phase->m_A[1][2][0] = gdelta;
        phase->m_A[1][2][1] = 1.0 - gdelta;

        phase->m_B[1][0][0] = gdelta;
        phase->m_B[1][0][1] = 1.0 - gdelta;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_2_2_2TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 3 2 : L Stable, two implict stages, three explict stages, second order IMEX

class IMEXdirk_2_3_2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_2_3_2TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(2, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_2_3_2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_2_3_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_2_3_2TimeIntegrationScheme>::AllocateSharedPtr(2, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_2_3_2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 2;
        phase->m_name = std::string("IMEXdirk_2_3_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 3;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        NekDouble lambda = (2.0 - sqrt(2.0)) / 2.0;
        NekDouble delta  = -2.0 * sqrt(2.0) / 3.0;

        phase->m_A[0][1][1] = lambda;
        phase->m_A[0][2][1] = 1.0 - lambda;
        phase->m_A[0][2][2] = lambda;

        phase->m_B[0][0][1] = 1.0 - lambda;
        phase->m_B[0][0][2] = lambda;

        phase->m_A[1][1][0] = lambda;
        phase->m_A[1][2][0] = delta;
        phase->m_A[1][2][1] = 1.0 - delta;

        phase->m_B[1][0][1] = 1.0 - lambda;
        phase->m_B[1][0][2] = lambda;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_2_3_2TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 3 3 : L Stable, two implict stages, three explict stages, third order IMEX

class IMEXdirk_2_3_3TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_2_3_3TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(3, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_2_3_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_2_3_3TimeIntegrationScheme>::AllocateSharedPtr(3, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_2_3_3");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 2;
        phase->m_name = std::string("IMEXdirk_2_2_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 3;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        NekDouble glambda = 0.788675134594813;

        phase->m_A[0][1][1] = glambda;
        phase->m_A[0][2][1] = 1.0 - 2.0 * glambda;
        phase->m_A[0][2][2] = glambda;

        phase->m_B[0][0][1] = 0.5;
        phase->m_B[0][0][2] = 0.5;

        phase->m_A[1][1][0] = glambda;
        phase->m_A[1][2][0] = glambda - 1.0;
        phase->m_A[1][2][1] = 2.0 * (1 - glambda);

        phase->m_B[1][0][1] = 0.5;
        phase->m_B[1][0][2] = 0.5;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_2_3_3TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 3 4 3 : L Stable, three implict stages, four explict stages, third order IMEX

class IMEXdirk_3_4_3TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_3_4_3TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(3, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_3_4_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_3_4_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_3_4_3TimeIntegrationScheme>::AllocateSharedPtr(3, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_3_4_3");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 3;
        phase->m_name = std::string("IMEXdirk_3_4_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 4;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        NekDouble lambda = 0.4358665215;

        phase->m_A[0][1][1] = lambda;
        phase->m_A[0][2][1] = 0.5 * (1.0 - lambda);
        phase->m_A[0][3][1] =
            0.25 * (-6.0 * lambda * lambda + 16.0 * lambda - 1.0);
        phase->m_A[0][2][2] = lambda;
        phase->m_A[0][3][2] =
            0.25 * (6.0 * lambda * lambda - 20.0 * lambda + 5.0);
        phase->m_A[0][3][3] = lambda;

        phase->m_B[0][0][1] =
            0.25 * (-6.0 * lambda * lambda + 16.0 * lambda - 1.0);
        phase->m_B[0][0][2] =
            0.25 * (6.0 * lambda * lambda - 20.0 * lambda + 5.0);
        phase->m_B[0][0][3] = lambda;

        phase->m_A[1][1][0] = 0.4358665215;
        phase->m_A[1][2][0] = 0.3212788860;
        phase->m_A[1][2][1] = 0.3966543747;
        phase->m_A[1][3][0] = -0.105858296;
        phase->m_A[1][3][1] = 0.5529291479;
        phase->m_A[1][3][2] = 0.5529291479;

        phase->m_B[1][0][1] =
            0.25 * (-6.0 * lambda * lambda + 16.0 * lambda - 1.0);
        phase->m_B[1][0][2] =
            0.25 * (6.0 * lambda * lambda - 20.0 * lambda + 5.0);
        phase->m_B[1][0][3] = lambda;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_3_4_3TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 4 4 3 : L Stable, four implict stages, four explict stages, third order IMEX

class IMEXdirk_4_4_3TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirk_4_4_3TimeIntegrationScheme(int order, std::string type) :
        TimeIntegrationScheme(3, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_4_4_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXdirk_4_4_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_4_4_3TimeIntegrationScheme>::AllocateSharedPtr(3, "");
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXdirk_4_4_3");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;
        phase->m_order = 3;
        phase->m_name = std::string("IMEXdirk_4_4_" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = 5;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        phase->m_A[0][1][1] = 1.0 / 2.0;
        phase->m_A[0][2][1] = 1.0 / 6.0;
        phase->m_A[0][2][2] = 1.0 / 2.0;
        phase->m_A[0][3][1] = -1.0 / 2.0;
        phase->m_A[0][3][2] = 1.0 / 2.0;
        phase->m_A[0][3][3] = 1.0 / 2.0;
        phase->m_A[0][4][1] = 3.0 / 2.0;
        phase->m_A[0][4][2] = -3.0 / 2.0;
        phase->m_A[0][4][3] = 1.0 / 2.0;
        phase->m_A[0][4][4] = 1.0 / 2.0;

        phase->m_B[0][0][1] = 3.0 / 2.0;
        phase->m_B[0][0][2] = -3.0 / 2.0;
        phase->m_B[0][0][3] = 1.0 / 2.0;
        phase->m_B[0][0][4] = 1.0 / 2.0;

        phase->m_A[1][1][0] = 1.0 / 2.0;
        phase->m_A[1][2][0] = 11.0 / 18.0;
        phase->m_A[1][2][1] = 1.0 / 18.0;
        phase->m_A[1][3][0] = 5.0 / 6.0;
        phase->m_A[1][3][1] = -5.0 / 6.0;
        phase->m_A[1][3][2] = 1.0 / 2.0;
        phase->m_A[1][4][0] = 1.0 / 4.0;
        phase->m_A[1][4][1] = 7.0 / 4.0;
        phase->m_A[1][4][2] = 3.0 / 4.0;
        phase->m_A[1][4][3] = -7.0 / 4.0;

        phase->m_B[1][0][0] = 1.0 / 4.0;
        phase->m_B[1][0][1] = 7.0 / 4.0;
        phase->m_B[1][0][2] = 3.0 / 4.0;
        phase->m_B[1][0][3] = -7.0 / 4.0;

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

        ASSERTL1(phase->VerifyIntegrationSchemeType(
                     phase->GetIntegrationSchemeType(), phase->m_A, phase->m_B,
                     phase->m_U, phase->m_V),
                 "Time integration scheme coefficients do not match its type");
    }

}; // end class IMEXdirk_4_4_3TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
