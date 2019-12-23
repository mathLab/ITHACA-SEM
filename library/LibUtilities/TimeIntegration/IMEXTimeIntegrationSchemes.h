///////////////////////////////////////////////////////////////////////////////
//
// File: IMEXTimeIntegrationSchemes.h
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
// Description: Combined header file for all basic IMEX time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// IMEX Order 1

class IMEXOrder1TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXOrder1TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXOrder1TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
    }

    virtual ~IMEXOrder1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXOrder1TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXOrder1");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;

        phase->m_numstages = 1;
        phase->m_numsteps  = 2;

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 1;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 0;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 1.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = 1.0;
        phase->m_B[1][1][0] = 1.0;
        phase->m_V[0][0]    = 1.0;
        phase->m_V[0][1]    = 1.0;

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

}; // end class IMEXOrder1TimeIntegrator

///////////////////////////////////////////////////////////////////////////////
// IMEX Order 2

class IMEXOrder2TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXOrder2TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(2);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[1] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXOrder1TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
        IMEXOrder2TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[1]);
    }

    virtual ~IMEXOrder2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXOrder2TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXOrder2");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;

        phase->m_numsteps  = 4;
        phase->m_numstages = 1;

        phase->m_numMultiStepValues = 2;
        phase->m_numMultiStepDerivs = 2;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);

        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 0;
        phase->m_timeLevelOffset[3] = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        NekDouble third = 1.0 / 3.0;

        phase->m_A[0] = Array<TwoD, NekDouble>(phase->m_numstages,
                                               phase->m_numstages, 2 * third);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U = Array<TwoD, NekDouble>(phase->m_numstages,
                                            phase->m_numsteps, 4 * third);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = 2 * third;
        phase->m_B[1][2][0] = 1.0;
        phase->m_U[0][1]    = -third;
        phase->m_U[0][3]    = -2 * third;

        phase->m_V[0][0] = 4 * third;
        phase->m_V[0][1] = -third;
        phase->m_V[0][2] = 4 * third;
        phase->m_V[0][3] = -2 * third;
        phase->m_V[1][0] = 1.0;
        phase->m_V[3][2] = 1.0;

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

}; // end class IMEXOrder2TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Order 3

class IMEXOrder3TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXOrder3TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(3);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[1] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[2] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_3_4_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
        IMEXdirk_3_4_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[1]);
        IMEXOrder3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[2]);
    }

    virtual ~IMEXOrder3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXOrder3TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXOrder3");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;

        phase->m_numsteps  = 6;
        phase->m_numstages = 1;

        phase->m_numMultiStepValues = 3;
        phase->m_numMultiStepDerivs = 3;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);

        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 2;
        phase->m_timeLevelOffset[3] = 0;
        phase->m_timeLevelOffset[4] = 1;
        phase->m_timeLevelOffset[5] = 2;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        NekDouble eleventh = 1.0 / 11.0;

        phase->m_A[0] = Array<TwoD, NekDouble>(
            phase->m_numstages, phase->m_numstages, 6 * eleventh);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U = Array<TwoD, NekDouble>(phase->m_numstages,
                                            phase->m_numsteps, 18 * eleventh);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = 6 * eleventh;
        phase->m_B[1][3][0] = 1.0;
        phase->m_U[0][1]    = -9 * eleventh;
        phase->m_U[0][2]    = 2 * eleventh;
        phase->m_U[0][4]    = -18 * eleventh;
        phase->m_U[0][5]    = 6 * eleventh;

        phase->m_V[0][0] = 18 * eleventh;
        phase->m_V[0][1] = -9 * eleventh;
        phase->m_V[0][2] = 2 * eleventh;
        phase->m_V[0][3] = 18 * eleventh;
        phase->m_V[0][4] = -18 * eleventh;
        phase->m_V[0][5] = 6 * eleventh;
        phase->m_V[1][0] = 1.0;
        phase->m_V[2][1] = 1.0;
        phase->m_V[4][3] = 1.0;
        phase->m_V[5][4] = 1.0;

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

}; // end class IMEXOrder3TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Order 3

class IMEXOrder4TimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXOrder4TimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_integration_phases    = TimeIntegrationSchemeDataVector(4);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[1] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[2] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));
        m_integration_phases[3] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0]);
        IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[1]);
        IMEXOrder3TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[2]);
        IMEXOrder4TimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[3]);
    }

    virtual ~IMEXOrder4TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p =
            MemoryManager<IMEXOrder4TimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXOrder4");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eIMEX;

        phase->m_numsteps  = 8;
        phase->m_numstages = 1;

        phase->m_numMultiStepValues = 4;
        phase->m_numMultiStepDerivs = 4;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);

        phase->m_timeLevelOffset[0] = 0;
        phase->m_timeLevelOffset[1] = 1;
        phase->m_timeLevelOffset[2] = 2;
        phase->m_timeLevelOffset[3] = 3;
        phase->m_timeLevelOffset[4] = 0;
        phase->m_timeLevelOffset[5] = 1;
        phase->m_timeLevelOffset[6] = 2;
        phase->m_timeLevelOffset[7] = 3;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(2);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(2);

        NekDouble twentyfifth = 1.0 / 25.0;

        phase->m_A[0] = Array<TwoD, NekDouble>(
            phase->m_numstages, phase->m_numstages, 12 * twentyfifth);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_A[1] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[1] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
        phase->m_U = Array<TwoD, NekDouble>(
            phase->m_numstages, phase->m_numsteps, 48 * twentyfifth);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        phase->m_B[0][0][0] = 12 * twentyfifth;
        phase->m_B[1][4][0] = 1.0;
        phase->m_U[0][1]    = -36 * twentyfifth;
        phase->m_U[0][2]    = 16 * twentyfifth;
        phase->m_U[0][3]    = -3 * twentyfifth;
        phase->m_U[0][5]    = -72 * twentyfifth;
        phase->m_U[0][7]    = -12 * twentyfifth;

        phase->m_V[0][0] = 48 * twentyfifth;
        phase->m_V[0][1] = -36 * twentyfifth;
        phase->m_V[0][2] = 16 * twentyfifth;
        phase->m_V[0][3] = -3 * twentyfifth;
        phase->m_V[0][4] = 48 * twentyfifth;
        phase->m_V[0][5] = -72 * twentyfifth;
        phase->m_V[0][6] = 48 * twentyfifth;
        phase->m_V[0][7] = -12 * twentyfifth;
        phase->m_V[1][0] = 1.0;
        phase->m_V[2][1] = 1.0;
        phase->m_V[3][2] = 1.0;
        phase->m_V[5][4] = 1.0;
        phase->m_V[6][5] = 1.0;
        phase->m_V[7][6] = 1.0;

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

}; // end class IMEXOrder4TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
