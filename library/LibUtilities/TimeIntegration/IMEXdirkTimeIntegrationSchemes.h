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

///////////////////////////////////////////////////////////////////////////////
//  IMEXDirk-(s, sigma, p), where s is the number of implicit stage
//  schemes, sigma is the number of explicit stage scheme and p is the
//  combined order of the scheme.

class IMEXdirkTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    IMEXdirkTimeIntegrationScheme(std::string variant, unsigned int order,
                                  std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme(variant, order, freeParams)
    {
        ASSERTL1(freeParams.size() == 2,
                 "IMEX DIRK Time integration scheme invalid number "
                 "of free parameters, expected two, received  " +
                 std::to_string(freeParams.size()));

        std::cerr << __LINE__ << std::endl;

        int s     = freeParams[0];
        int sigma = freeParams[1];

        std::cerr << __LINE__ << "  " << s << "  " << sigma << std::endl;

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        if( order == 1 && s == 1 && sigma == 1 )
        {
            // This phase is Forward Backward Euler which has two steps.
            IMEXdirkTimeIntegrationScheme::SetupSchemeData_1_1_1(
                m_integration_phases[0]);
        }
        else
        {
            IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[0], order, freeParams);
        }
    }

    virtual ~IMEXdirkTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetFullName() const
    {
        return m_integration_phases[m_integration_phases.size()-1]->m_name;
    }

    LUE virtual std::string GetName() const
    {
        return std::string("IMEX");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    unsigned int order,
                                    std::vector<NekDouble> freeParams)
    {
        // This parsing is a hack of stuffing two values into a string.
        ASSERTL1(freeParams.size() == 2,
                 "IMEX DIRK Time integration scheme invalid number "
                 "of free parameters, expected two, received  " +
                 std::to_string(freeParams.size()));

        int s     = freeParams[0];
        int sigma = freeParams[1];

        phase->m_schemeType = eIMEX;
        phase->m_variant = "dirk";
        phase->m_order = order;
        phase->m_name = "IMEXdirk"
            "_" + std::to_string(s) +
            "_" + std::to_string(sigma) +
            "_" + std::to_string(phase->m_order);

        phase->m_numsteps  = 1;
        phase->m_numstages = s + 1;

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

        if( s == 1 && sigma == 2 && order == 1)
        {
            SetupSchemeData_1_2_1(phase);
        }
        else if( s == 1 && sigma == 2 && order == 2)
        {
            SetupSchemeData_1_2_2(phase);
        }
        else if( s == 2 && sigma == 2 && order == 2)
        {
            SetupSchemeData_2_2_2(phase);
        }
        else if( s == 2 && sigma == 3 && order == 2)
        {
            SetupSchemeData_2_3_2(phase);
        }
        else if( s == 2 && sigma == 3 && order == 3)
        {
            SetupSchemeData_2_3_3(phase);
        }
        else if( s == 3 && sigma == 4 && order == 3)
        {
            SetupSchemeData_3_4_3(phase);
        }
        else if( s == 4 && sigma == 4 && order == 3)
        {
            SetupSchemeData_4_4_3(phase);
        }
        else
        {
            ASSERTL1(false,
                     "IMEX DIRK Time integration scheme bad type. ");
        }

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 1 1 : Forward - Backward Euler IMEX
    LUE static void SetupSchemeData_1_1_1(TimeIntegrationSchemeDataSharedPtr &phase)
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

        phase->CheckAndVerify();
    }

 ///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 2 1 : Forward - Backward Euler IMEX w/B implicit == B explicit
  LUE static void SetupSchemeData_1_2_1(TimeIntegrationSchemeDataSharedPtr &phase)
  {
      phase->m_A[0][1][1] = 1.0;
      phase->m_B[0][0][1] = 1.0;

      phase->m_A[1][1][0] = 1.0;
      phase->m_B[1][0][1] = 1.0;

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 2 2 : Implict-Explicit Midpoint IMEX
  LUE static void SetupSchemeData_1_2_2(TimeIntegrationSchemeDataSharedPtr &phase)
  {
      phase->m_A[0][1][1] = 1.0 / 2.0;
      phase->m_B[0][0][1] = 1.0;

      phase->m_A[1][1][0] = 1.0 / 2.0;
      phase->m_B[1][0][1] = 1.0;

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 2 2 : L Stable, two stage, second order IMEX
  LUE static void SetupSchemeData_2_2_2(TimeIntegrationSchemeDataSharedPtr &phase)
  {
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

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 3 2 : L Stable, two stage, second order IMEX
  LUE static void SetupSchemeData_2_3_2(TimeIntegrationSchemeDataSharedPtr &phase)
  {
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

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 3 3 : L Stable, two stage, third order IMEX
  LUE static void SetupSchemeData_2_3_3(TimeIntegrationSchemeDataSharedPtr &phase)
  {
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

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 3 4 3 : L Stable, three stage, third order IMEX
  LUE static void SetupSchemeData_3_4_3(TimeIntegrationSchemeDataSharedPtr &phase)
  {
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

      // U and V set to 1 when allocated.
  }

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 4 4 3 : L Stable, four stage, third order IMEX
  LUE static void SetupSchemeData_4_4_3(TimeIntegrationSchemeDataSharedPtr &phase)
  {
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

      // U and V set to 1 when allocated.
  }

}; // end class IMEXdirkTimeIntegrationScheme


///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 1 1 1 : Forward - Backward Euler IMEX
class IMEXdirk_1_1_1TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_1_1_1TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 1, std::vector<NekDouble>{1, 1})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_1_1_1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_1_1_1TimeIntegrationScheme>::AllocateSharedPtr("", 1, std::vector<NekDouble>{1, 1});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_1_1_1TimeIntegrator


// IMEX Dirk 1 2 1 : Forward - Backward Euler IMEX w/B implicit == B explicit
class IMEXdirk_1_2_1TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_1_2_1TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 1, std::vector<NekDouble>{1, 2})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_1_2_1TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 1, std::vector<NekDouble>{1, 2});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_1_2_1TimeIntegrator


// IMEX Dirk 1 2 2 : Implict-Explicit Midpoint IMEX
class IMEXdirk_1_2_2TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_1_2_2TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 2, std::vector<NekDouble>{1, 2})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_1_2_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 2, std::vector<NekDouble>{1, 2});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_1_2_2TimeIntegrator


// IMEX Dirk 2 2 2 : L Stable, two stage, second order IMEX
class IMEXdirk_2_2_2TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_2_2_2TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 2, std::vector<NekDouble>{2, 2})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_2_2_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirk_2_2_2TimeIntegrationScheme>::AllocateSharedPtr("", 2, std::vector<NekDouble>{2, 2});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_2_2_2TimeIntegrationScheme


// IMEX Dirk 2 3 2 : L Stable, two stage, second order IMEX
class IMEXdirk_2_3_2TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_2_3_2TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 2, std::vector<NekDouble>{2, 3})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_2_3_2TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 2, std::vector<NekDouble>{2, 3});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_2_3_2TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 2 3 3 : L Stable, two stage, third order IMEX
class IMEXdirk_2_3_3TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_2_3_3TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 3, std::vector<NekDouble>{2, 3})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_2_3_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 3, std::vector<NekDouble>{2, 3});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_2_3_3TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 3 4 3 : L Stable, three stage, third order IMEX
class IMEXdirk_3_4_3TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_3_4_3TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 3, std::vector<NekDouble>{3, 4})
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_3_4_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 3, std::vector<NekDouble>{3, 4});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_3_4_3TimeIntegrationScheme

///////////////////////////////////////////////////////////////////////////////
// IMEX Dirk 4 4 3 : L Stable, four stage, third order IMEX
class IMEXdirk_4_4_3TimeIntegrationScheme : public IMEXdirkTimeIntegrationScheme
{
public:
    IMEXdirk_4_4_3TimeIntegrationScheme(std::string variant, unsigned int order,
                                        std::vector<NekDouble> freeParams) :
        IMEXdirkTimeIntegrationScheme("", 3, std::vector<NekDouble>{4, 4})
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);
        boost::ignore_unused(freeParams);
    }

    virtual ~IMEXdirk_4_4_3TimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(variant);
        boost::ignore_unused(freeParams);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            IMEXdirkTimeIntegrationScheme>::AllocateSharedPtr("", 3, std::vector<NekDouble>{4, 4});

        return p;
    }

    static std::string className;

}; // end class IMEXdirk_4_4_3TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
