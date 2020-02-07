///////////////////////////////////////////////////////////////////////////////
//
// File: BDFImplicitTimeIntegrationSchemes.h
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
// Description: Combined header file for all BDF Implicit based time integration
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
// BDF Implicit Order N

class BDFImplicitTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    BDFImplicitTimeIntegrationScheme(std::string variant, unsigned int order,
                                     std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented.

        // Methods with order > 6 are not zero-stable.
        ASSERTL1(0 < order && order <= 4,
                 "BDFImplicit Time integration scheme bad order (1-4): " +
                 std::to_string(order));

        m_integration_phases = TimeIntegrationSchemeDataVector(order);

        for( unsigned int n=0; n<order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));
        }

        // Next to last phase
        if( order > 1 )
            BDFImplicitTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[order-2], order-1);

        // Last phase
        BDFImplicitTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[order-1], order);

        // Initial phases
        switch( order )
        {
            case 1:
                // No intial phase.
                break;

            case 2:
                // Intial phase set above
                break;

            case 3:
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], 3, std::vector<NekDouble>{3, 4});
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], 3, std::vector<NekDouble>{3, 4});
                break;

            case 4:
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0], 3, std::vector<NekDouble>{2, 3});
                IMEXdirkTimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1], 3, std::vector<NekDouble>{2, 3});
                break;

            default:
                ASSERTL1(false,
                         "BDFImplicit Time integration scheme bad order: " +
                         std::to_string(order));
        }
    }

    virtual ~BDFImplicitTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
						 std::vector<NekDouble> freeParams)
  {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
          BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("BDFImplicit");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    unsigned int order)
    {
        const NekDouble ABcoefficients[5] = {      0.,
                                                   1.,    // 1st Order
                                              2. / 3.,    // 2nd Order
                                              6. / 11.,   // 3rd Order
                                             12. / 25. }; // 4th Order

        // The 3rd and 4th order tableaus have not been validated!!!!!
        // The 3rd and 4th order tableaus have not been validated!!!!!
        const NekDouble UVcoefficients[5][4] =
            { {       0.,      0.,     0.,       0. },
              // 1st Order
              {       1.,      0.,     0.,       0. },
              // 2nd Order
              {  4./ 3.,  -1./ 3.,     0.,       0. },
              // 3rd Order
              { 18./11.,  -9./11.,  2./11.,      0. },
              // 4th Order
              { 48./25., -36./25., 16./25., -3./25. } };

        phase->m_schemeType = eDiagonallyImplicit;
        phase->m_order = order;
        phase->m_name = std::string("BDFImplicitOrder" +
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
        phase->m_A[0][0][0] = ABcoefficients[phase->m_order];
        phase->m_B[0][0][0] = ABcoefficients[phase->m_order];

        // U/V Coefficients for first row additional columns
        for( int n=0; n<phase->m_order; ++n )
        {
            phase->m_U[0][n] = UVcoefficients[phase->m_order][n];
            phase->m_V[0][n] = UVcoefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_V[n][n-1] = 1.0;
        }

        phase->m_numMultiStepValues = phase->m_order;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);

        // For order >= 1 values are needed.
        for( int n=0; n<phase->m_order; ++n )
        {
            phase->m_timeLevelOffset[n] = n;
        }

        phase->CheckAndVerify();
    }

}; // end class BDFImplicitTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class BDFImplicitOrder1TimeIntegrationScheme :
    public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder1TimeIntegrationScheme(std::string variant, unsigned int order,
                                           std::vector<NekDouble> freeParams) :
      BDFImplicitTimeIntegrationScheme("", 1, freeParams)
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
            BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr("", 1, freeParams);
        return p;
    }

    static std::string className;

}; // end class BDFImplicitOrder1TimeIntegrationScheme


class BDFImplicitOrder2TimeIntegrationScheme :
    public BDFImplicitTimeIntegrationScheme
{
public:
    BDFImplicitOrder2TimeIntegrationScheme(std::string variant, unsigned int order,
                                           std::vector<NekDouble> freeParams) :
        BDFImplicitTimeIntegrationScheme("", 2, freeParams)
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
            BDFImplicitTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);
        return p;
    }

    static std::string className;

}; // end class BDFImplicitOrder2TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
