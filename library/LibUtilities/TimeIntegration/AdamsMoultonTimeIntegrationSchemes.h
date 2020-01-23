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

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

#include <LibUtilities/TimeIntegration/EulerTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>

namespace Nektar
{
namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
// Adams Moulton Order N

class AdamsMoultonTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    AdamsMoultonTimeIntegrationScheme(std::string variant, unsigned int order,
                                      std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(0 < order && order <= 4,
                 "AdamsMoulton Time integration scheme bad order (1-4): " +
                 std::to_string(order));

        m_integration_phases = TimeIntegrationSchemeDataVector(order);

        for( unsigned int n=0; n<order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));
        }

        // Next to last phase
        if( order > 1 )
            AdamsMoultonTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order-2], order-1);

        // Last phase
        AdamsMoultonTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[order-1], order);

        // Initial phases
        switch( order )
        {
            case 1:
                // No intial phases.
                break;

            case 2:
                // Why forward euler and not backward euler???
                ForwardEulerTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0]);
                break;

            case 3:
                // The first and second phases needed to be set correctly
                ForwardEulerTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0]);
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], 3, std::vector<NekDouble>{3, 4});
                break;

            case 4:
                // The first and second phases needed to be set correctly
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], 3, std::vector<NekDouble>{2, 3});
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], 3, std::vector<NekDouble>{2, 3});
                break;

            default:
              ASSERTL1(false,
                       "AdamsMoulton Time integration scheme bad order: " +
                       std::to_string(order));
        }
    }

    virtual ~AdamsMoultonTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("AdamsMoulton");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    int order)
    {
        // The 3rd and 4th order tableaus have not been validated!!!!!
        // The 3rd and 4th order tableaus have not been validated!!!!!
        const NekDouble coefficients[5][4] =
            { {      0.,       0.,      0.,     0. },
              // 1st Order
              {      1.,       0.,      0.,     0. },
              // 2nd Order
              {  1./ 2.,   1./ 2.,      0.,     0. },
              // 3rd Order
              {  5./12.,   8./12., -1./12.,     0. },
              // 4th Order
              {  9./24.,  19./24., -5./24.,  1./24. } };

        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order = order;
        phase->m_name = std::string("AdamsMoultonOrder" +
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

        // When multiple steps are taken A/B[0][0] and U/V[0][1...s]
        // must be weighted so the time contribution is correct.

        // A/B Coefficient for first row first column
        phase->m_A[0][0][0] = coefficients[phase->m_order][0];
        phase->m_B[0][0][0] = coefficients[phase->m_order][0];

        // B evaluation value shuffling second row first column
        if( phase->m_order > 1 )
        {
            phase->m_B[0][1][0] = 1.0; // constant 1
        }

        // U/V Coefficient for first row first column
        phase->m_U[0][0] = 1.0;
        phase->m_V[0][0] = 1.0;

        // U/V Coefficients for first row additional columns
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_U[0][n] = coefficients[phase->m_order][n];
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
            phase->m_timeLevelOffset[n] = n-1;
        }

        phase->CheckAndVerify();
    }

}; // end class AdamsMoultonTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class AdamsMoultonOrder1TimeIntegrationScheme :
    public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder1TimeIntegrationScheme(std::string variant, unsigned int order,
                                            std::vector<NekDouble> freeParams) :
        AdamsMoultonTimeIntegrationScheme("", 1, freeParams)
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
            AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr("", 1, freeParams);

        return p;
    }

    static std::string className;

}; // end class AdamsMoultonOrder1TimeIntegrationScheme


class AdamsMoultonOrder2TimeIntegrationScheme :
    public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder2TimeIntegrationScheme(std::string variant, unsigned int order,
                                            std::vector<NekDouble> freeParams) :
        AdamsMoultonTimeIntegrationScheme("", 2, freeParams)
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
            AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);

        return p;
    }

    static std::string className;

}; // end class AdamsMoultonOrder2TimeIntegrationScheme


class AdamsMoultonOrder3TimeIntegrationScheme :
    public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder3TimeIntegrationScheme(std::string variant, unsigned int order,
                                            std::vector<NekDouble> freeParams) :
        AdamsMoultonTimeIntegrationScheme("", 3, freeParams)
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
            AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr("", 3, freeParams);

        return p;
    }

    static std::string className;

}; // end class AdamsMoultonOrder3TimeIntegrationScheme


class AdamsMoultonOrder4TimeIntegrationScheme :
    public AdamsMoultonTimeIntegrationScheme
{
public:
    AdamsMoultonOrder4TimeIntegrationScheme(std::string variant, unsigned int order,
                                            std::vector<NekDouble> freeParams) :
        AdamsMoultonTimeIntegrationScheme("", 4, freeParams)
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
            AdamsMoultonTimeIntegrationScheme>::AllocateSharedPtr("", 4, freeParams);

        return p;
    }

    static std::string className;

}; // end class AdamsMoultonOrder4TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
