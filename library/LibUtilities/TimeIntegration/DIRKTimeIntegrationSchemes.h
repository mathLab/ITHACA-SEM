///////////////////////////////////////////////////////////////////////////////
//
// File: DIRKTimeIntegrationSchemes.h
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
// Description: Combined header file for all basic DIRK based time integration
// schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>

#define LUE LIB_UTILITIES_EXPORT
using namespace std;

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// DIRK Order N
class DIRKTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
  DIRKTimeIntegrationScheme(std::string variant, unsigned int order,
                            std::vector<NekDouble> freeParams) :
    TimeIntegrationScheme(variant, order, freeParams)
    {
        // Currently 2nd and 3rd order are implemented.
        ASSERTL1(2 <= order && order <= 4,
                 "Runge Kutta Time Diagonally Implicit integration scheme bad order (2-4): " +
                 std::to_string(order));

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        ASSERTL1((variant == "" || variant == "ES5" || variant == "ES6"),
                 "DIRK Time integration scheme bad variant: "
		 + variant + ". "
                  "Must blank, 'ES5', or 'ES6'");
        
        if ("" == variant)
        {
            DIRKTimeIntegrationScheme::SetupSchemeData( m_integration_phases[0],
                                                    order);
        }
        else
        {
            DIRKTimeIntegrationScheme::SetupSchemeDataESDIRK( 
                m_integration_phases[0], variant, order, freeParams);
        }
    }

    virtual ~DIRKTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, 
                                                 unsigned int order,
                                   std::vector<NekDouble> freeParams)
    {
         TimeIntegrationSchemeSharedPtr p =
           MemoryManager<DIRKTimeIntegrationScheme>::AllocateSharedPtr(
                                             variant, order, freeParams);
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("DIRK");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    unsigned int order)
    {
        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order = order;
        phase->m_name = std::string("DIRKOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 1;
        phase->m_numstages = phase->m_order;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);

        switch( phase->m_order )
        {
            case 2:
            {
                // Two-stage, 2nd order Diagonally Implicit Runge
                // Kutta method. It is A-stable if and only if x ≥ 1/4
                // and is L-stable if and only if x equals one of the
                // roots of the polynomial x^2 - 2x + 1/2, i.e. if
                // x = 1 ± sqrt(2)/2.
                NekDouble lambda = (2.0 - sqrt(2.0)) / 2.0;

                phase->m_A[0][0][0] = lambda;
                phase->m_A[0][1][0] = 1.0 - lambda;
                phase->m_A[0][1][1] = lambda;

                phase->m_B[0][0][0] = 1.0 - lambda;
                phase->m_B[0][0][1] = lambda;
            }
            break;

            case 3:
            {
                // Three-stage, 3rd order, L-stable Diagonally Implicit
                // Runge Kutta method:
                NekDouble lambda = 0.4358665215;

                phase->m_A[0][0][0] = lambda;
                phase->m_A[0][1][0] = 0.5 * (1.0 - lambda);
                phase->m_A[0][2][0] =
                    0.25 * (-6.0 * lambda * lambda + 16.0 * lambda - 1.0);

                phase->m_A[0][1][1] = lambda;

                phase->m_A[0][2][1] =
                    0.25 * (6.0 * lambda * lambda - 20.0 * lambda + 5.0);
                phase->m_A[0][2][2] = lambda;

                phase->m_B[0][0][0] =
                    0.25 * (-6.0 * lambda * lambda + 16.0 * lambda - 1.0);
                phase->m_B[0][0][1] =
                    0.25 * (6.0 * lambda * lambda - 20.0 * lambda + 5.0);
                phase->m_B[0][0][2] = lambda;
            }
            break;
        }

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

    LUE static void SetupSchemeDataESDIRK(
        TimeIntegrationSchemeDataSharedPtr &phase,
        std::string                         variant,
        unsigned int                        order,
        std::vector<NekDouble>              freeParams)
    {
        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order = order;
        phase->m_name = std::string("DIRKOrder" +
                                    std::to_string(phase->m_order)
                                    + variant);

        phase->m_numsteps  = 1;
        char const &stage = variant.back();
        phase->m_numstages =  std::atoi(&stage);

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  1.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  1.0);
        
        const NekDouble ConstSqrt2  = sqrt(2.0);
        switch( phase->m_order )
        {
            case 3:
            {
                ASSERTL0(5==phase->m_numstages,
                    std::string("DIRKOrder3_ES" +
                    std::to_string(phase->m_numstages)+
                    " not defined"));
                NekDouble lambda;
                if(freeParams.size())
                {
                    lambda = freeParams[0];
                }
                else
                {
                    lambda = 9.0/40.0;
                }

                phase->m_A[0][0][0] = 0.0;
                phase->m_A[0][1][0] = lambda;
                phase->m_A[0][2][0] = 9.0 * (1.0 + ConstSqrt2) / 80.0;
                phase->m_A[0][3][0] = (22.0 + 15.0 * ConstSqrt2) / 
                                      (80.0 * (1 + ConstSqrt2));
                phase->m_A[0][4][0] = (2398.0 + 1205.0 * ConstSqrt2) / 
                                      (2835.0 * (4 + 3.0 * ConstSqrt2));

                phase->m_A[0][1][1] = phase->m_A[0][1][0];
                phase->m_A[0][2][1] = phase->m_A[0][2][0];
                phase->m_A[0][3][1] = phase->m_A[0][3][0];
                phase->m_A[0][4][1] = phase->m_A[0][4][0];

                phase->m_A[0][2][2] = lambda;
                phase->m_A[0][3][2] = -7.0 / (40.0 * (1.0 + ConstSqrt2));
                phase->m_A[0][4][2] = -2374 * (2.0 + ConstSqrt2)/ 
                                      (2835.0 * (4.0 + 3.0 * ConstSqrt2));

                phase->m_A[0][3][3] = lambda;
                phase->m_A[0][4][3] = 5827.0 / 7560.0;

                phase->m_A[0][4][4] = lambda;

                phase->m_B[0][0][0] = phase->m_A[0][4][0];
                phase->m_B[0][0][1] = phase->m_A[0][4][1];
                phase->m_B[0][0][2] = phase->m_A[0][4][2];
                phase->m_B[0][0][3] = phase->m_A[0][4][3];
                phase->m_B[0][0][4] = phase->m_A[0][4][4];
            }
            break;

            case 4:
            {
                ASSERTL0(6 == phase->m_numstages,
                    std::string("DIRKOrder4_ES" +
                    std::to_string(phase->m_numstages)+
                    " not defined"));
                NekDouble lambda;
                if (freeParams.size())
                {
                    lambda = freeParams[0];
                }
                else
                {
                    lambda = 1.0 / 4.0;
                }

                const NekDouble ConstSqrt2  = sqrt(2.0);

                phase->m_A[0][0][0] = 0.0;
                phase->m_A[0][1][0] = lambda;
                phase->m_A[0][2][0] = (1.0 - ConstSqrt2) / 8.0;
                phase->m_A[0][3][0] = (5.0 - 7.0 * ConstSqrt2) / 64.0;
                phase->m_A[0][4][0] = (-13796.0 - 54539 * ConstSqrt2) / 
                                      125000.0;
                phase->m_A[0][5][0] = (1181.0 - 987.0 * ConstSqrt2) /
                                       13782.0;

                phase->m_A[0][1][1] = phase->m_A[0][1][0];
                phase->m_A[0][2][1] = phase->m_A[0][2][0];
                phase->m_A[0][3][1] = phase->m_A[0][3][0];
                phase->m_A[0][4][1] = phase->m_A[0][4][0];
                phase->m_A[0][5][1] = phase->m_A[0][5][0];

                phase->m_A[0][2][2] = lambda;
                phase->m_A[0][3][2] = 7.0 * (1.0 + ConstSqrt2) /
                                      32.0;
                phase->m_A[0][4][2] = (506605.0 + 132109.0 * ConstSqrt2) /
                                      437500.0;
                phase->m_A[0][5][2] = 47.0 * (-267.0 + 1783.0 * ConstSqrt2) /
                                      273343.0;

                phase->m_A[0][3][3] = lambda;
                phase->m_A[0][4][3] = 166.0 * (-97.0 + 376.0 * ConstSqrt2) /
                                      109375.0;
                phase->m_A[0][5][3] = -16.0 * (-22922.0 + 3525.0 * ConstSqrt2) /
                                      571953.0;

                phase->m_A[0][4][4] = lambda;
                phase->m_A[0][5][4] = -15625.0 * (97.0 + 376.0 * ConstSqrt2) /
                                      90749876.0;

                phase->m_A[0][5][5] = lambda;
            }
            break;
            default:
            {
                ASSERTL0(false, std::string("ESDIRK of order" +
                                    std::to_string(phase->m_order)+
                                    " not defined"));
                break;
            }
        }

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

}; // end class DIRKTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class DIRKOrder2TimeIntegrationScheme :
    public DIRKTimeIntegrationScheme
{
public:
    DIRKOrder2TimeIntegrationScheme(std::string variant, unsigned int order,
                                    std::vector<NekDouble> freeParams) :
        DIRKTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            DIRKTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);

        return p;
    }

    static std::string className;

}; // end class DIRKOrder2TimeIntegrationScheme


class DIRKOrder3TimeIntegrationScheme :
    public DIRKTimeIntegrationScheme
{
public:
    DIRKOrder3TimeIntegrationScheme(std::string variant, unsigned int order,
                                    std::vector<NekDouble> freeParams) :
        DIRKTimeIntegrationScheme("", 3, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            DIRKTimeIntegrationScheme>::AllocateSharedPtr("", 3, freeParams);

        return p;
    }

    static std::string className;

}; // end class DIRKOrder3TimeIntegrationScheme

class DIRKOrder3_ES5TimeIntegrationScheme :
    public DIRKTimeIntegrationScheme
{
public:
    DIRKOrder3_ES5TimeIntegrationScheme(std::string variant, unsigned int order,
                                    std::vector<NekDouble> freeParams) :
        DIRKTimeIntegrationScheme("ES5", 3, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            DIRKTimeIntegrationScheme>::AllocateSharedPtr("ES5", 3, freeParams);

        return p;
    }

    static std::string className;

}; // end class DIRKOrder3_ES5TimeIntegrationScheme

class DIRKOrder4_ES6TimeIntegrationScheme :
    public DIRKTimeIntegrationScheme
{
public:
    DIRKOrder4_ES6TimeIntegrationScheme(std::string variant, unsigned int order,
                                    std::vector<NekDouble> freeParams) :
        DIRKTimeIntegrationScheme("ES6", 4, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order,
        std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            DIRKTimeIntegrationScheme>::AllocateSharedPtr("ES6", 4, freeParams);

        return p;
    }

    static std::string className;

}; 


} // end namespace LibUtilities
} // end namespace Nektar
