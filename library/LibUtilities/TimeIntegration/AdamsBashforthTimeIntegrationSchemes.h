///////////////////////////////////////////////////////////////////////////////
//
// File: AdamsBashforthTimeIntegrationSchemes.h
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
// Description: Combined header file for all Adams Bashforth based time
// integration schemes.
//
///////////////////////////////////////////////////////////////////////////////

// Note : If adding a new integrator be sure to register the
// integrator with the Time Integration Scheme Facatory in
// SchemeInitializor.cpp.

#pragma once

#define LUE LIB_UTILITIES_EXPORT

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>
#include <LibUtilities/TimeIntegration/RungeKuttaTimeIntegrationSchemes.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Adams Bashforth Order N

class AdamsBashforthTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
  AdamsBashforthTimeIntegrationScheme(std::string variant, unsigned int order,
				      std::vector<NekDouble> freeParams) :
    TimeIntegrationScheme(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(0 < order && order <= 4,
                 "AdamsBashforth Time integration scheme bad order (1-4): " +
                 std::to_string(order));

        m_integration_phases = TimeIntegrationSchemeDataVector(order);

        for( unsigned int n=0; n<order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));
        }

        // Next to last phase
        if( order > 1 )
            AdamsBashforthTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order-2], order-1);

        // Last phase
        AdamsBashforthTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[order-1], order);

        // Initial phases
        switch( order )
        {
            case 1:
                // No intial phases.
                break;

            case 2:
                // Done above.
                break;

            case 3:
	        // Order 2
                RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], "", 2, std::vector<NekDouble>());
                break;

            case 4:
	        // SSP Order 3
	        RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], "SSP", 3, std::vector<NekDouble>());
		// SSP Order 3
                RungeKuttaTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], "SSP", 3, std::vector<NekDouble>());
                break;

            default:
              ASSERTL1(false,
                       "AdamsBashforth Time integration scheme bad order: " +
                       std::to_string(order));
        }
    }

    virtual ~AdamsBashforthTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("AdamsBashforth");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    int order)
    {
        const NekDouble coefficients[5][4] =
            { {      0.,       0.,      0.,      0. },
              // 1st Order
              {      1.,       0.,      0.,      0. },
              // 2nd Order
              {  3./ 2.,  -1./ 2.,      0.,      0. },
              // 3rd Order
              { 23./12., -16./12.,  5./12.,      0. },
              // 4th Order
              { 55./24., -59./24., 37./24., -9./24.} };

        phase->m_schemeType = eExplicit;
        phase->m_order = order;
        phase->m_name = std::string("AdamsBashforthOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = phase->m_order;
        phase->m_numstages = 1;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numstages, 0.0);
        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps,  0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps,  phase->m_numsteps,  0.0);

        // Coefficients

        // When multiple steps are taken B[0][0] and V[0][1...s] must be
        // weighted so the time contribution is correct.

        // B Coefficient for first row first column
        phase->m_B[0][0][0] = coefficients[phase->m_order][0];

        // B evaluation value shuffling second row first column
        if( phase->m_order > 1 )
        {
            phase->m_B[0][1][0] = 1.0; // constant 1
        }

        // U Curent time step evaluation first row first column
        phase->m_U[0][0] = 1.0;
        phase->m_V[0][0] = 1.0;

        // V Coefficients for first row additional columns
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_V[0][n] = coefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for( int n=2; n<phase->m_order; ++n )
        {
            phase->m_V[n][n-1] = 1.0;
        }
        
        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = phase->m_order-1;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        // For order > 1 derivatives are needed.
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_timeLevelOffset[n] = n;
        }

        phase->CheckAndVerify();
    }

}; // end class AdamsBashforthTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class AdamsBashforthOrder1TimeIntegrationScheme :
    public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder1TimeIntegrationScheme(std::string variant, unsigned int order,
					      std::vector<NekDouble> freeParams) :
      AdamsBashforthTimeIntegrationScheme("", 1, freeParams)
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
            AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 1, freeParams);
        return p;
    }

    static std::string className;

}; // end class AdamsBashforthOrder1TimeIntegrationScheme


class AdamsBashforthOrder2TimeIntegrationScheme :
    public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder2TimeIntegrationScheme(std::string variant, unsigned int order, std::vector<NekDouble> freeParams) :
      AdamsBashforthTimeIntegrationScheme("", 2, freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);
    }

    static TimeIntegrationSchemeSharedPtr create(
        std::string variant, unsigned int order, std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(variant);
        boost::ignore_unused(order);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);
        return p;
    }

    static std::string className;

}; // end class AdamsBashforthOrder2TimeIntegrationScheme


class AdamsBashforthOrder3TimeIntegrationScheme :
    public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder3TimeIntegrationScheme(std::string variant, unsigned int order,
					      std::vector<NekDouble> freeParams) :
      AdamsBashforthTimeIntegrationScheme("", 3, freeParams)
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
          AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 3, freeParams);
        return p;
    }

    static std::string className;

}; // end class AdamsBashforthOrder3TimeIntegrationScheme


class AdamsBashforthOrder4TimeIntegrationScheme :
    public AdamsBashforthTimeIntegrationScheme
{
public:
    AdamsBashforthOrder4TimeIntegrationScheme(std::string variant, unsigned int order,
					      std::vector<NekDouble> freeParams) :
      AdamsBashforthTimeIntegrationScheme("", 4, freeParams)
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
          AdamsBashforthTimeIntegrationScheme>::AllocateSharedPtr("", 4, freeParams);
        return p;
    }

    static std::string className;

}; // end class AdamsBashforthOrder4TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
