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
// IMEX Order N

class IMEXTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
  IMEXTimeIntegrationScheme(int order, std::string type) :
    TimeIntegrationScheme(order, type)
    {
        // Currently up to 4th order is implemented.
        ASSERTL1(0 < order && order <= 4,
                 "IMEX Time integration scheme bad order (1-4): " +
                 std::to_string(order));

        m_integration_phases = TimeIntegrationSchemeDataVector(m_order);

        for( unsigned int n=0; n<m_order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));
        }

        // Next to last phase
        if( m_order > 1 )
            IMEXTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[m_order-2], m_order-1);

        // Last phase
        IMEXTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[m_order-1], m_order);

        // Initial phases
        switch( m_order )
        {
            case 1:
                // No intial phase.
                break;

            case 2:
                // Intial phase set above
                break;

            case 3:
                IMEXdirk_3_4_3TimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0]);
                IMEXdirk_3_4_3TimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1]);
                break;

            case 4:
                IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[0]);
                IMEXdirk_2_3_3TimeIntegrationScheme::SetupSchemeData(
                    m_integration_phases[1]);
                break;

            default:
                ASSERTL1(false,
                         "IMEX Time integration scheme bad order: " +
                         std::to_string(order));
        }
    }

    virtual ~IMEXTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        TimeIntegrationSchemeSharedPtr p =
	  MemoryManager<IMEXTimeIntegrationScheme>::AllocateSharedPtr(order, type);
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("IMEXOrder" + std::to_string(m_order));
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    int order)
    {
        const NekDouble ABcoefficients[5] = {      0.,
                                                   1.,    // 1st Order
                                               2./ 3.,    // 2nd Order
                                               6./11.,    // 3rd Order
                                              12./25. };  // 4th Order

	// Nsteps = 2 * order
        const NekDouble UVcoefficients[5][8] =
            { {         0.,    0.,     0.,        0.,
                        0.,    0.,     0.,        0. },
              // 1st Order
              {         1.,    1.,     0.,        0.,
                        0.,    0.,     0.,        0. },
              // 2nd Order
              {  4./ 3.,  -1./ 3.,  4./3.,   -2./ 3.,
                 0.,           0.,     0.,        0. },
              // 3rd Order
              {  18./11.,  -9./11.,  2./11.,  18./11.,
                -18./11.,   6./11.,      0.,       0. },
              // 4th Order
              { 48./25., -36./25., 16./25.,  -3./25.,
                48./25., -72./25., 48./25., -12./25. } };

        phase->m_schemeType = eIMEX;
        phase->m_order = order;
        phase->m_name = std::string("IMEXOrder" +
                                    std::to_string(phase->m_order));

        phase->m_numsteps  = 2 * phase->m_order;
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

        // Coefficients
        phase->m_B[1][phase->m_order][0] =  1.0;

        phase->m_A[0][0][0] = ABcoefficients[phase->m_order];
        phase->m_B[0][0][0] = ABcoefficients[phase->m_order];

        for( int n=0; n<2*phase->m_order; ++n )
        {
            phase->m_U[0][n] = UVcoefficients[phase->m_order][n];
            phase->m_V[0][n] = UVcoefficients[phase->m_order][n];
        }

        // V evaluation value shuffling row n column n-1
        for( int n=1; n<2*phase->m_order; ++n )
        {
            if( n != phase->m_order )
                phase->m_V[n][n-1] = 1.0; // constant 1
        }

        phase->m_numMultiStepValues = phase->m_order;
        phase->m_numMultiStepDerivs = phase->m_order;

        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);

        // Values and derivatives needed.
        for( int n=0; n<phase->m_order; ++n )
        {
            phase->m_timeLevelOffset[n] = n;
            phase->m_timeLevelOffset[phase->m_order+n] = n;
        }

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

}; // end class IMEXTimeIntegrationScheme

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class IMEXOrder1TimeIntegrationScheme :
    public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder1TimeIntegrationScheme(int order, std::string type) :
        IMEXTimeIntegrationScheme(1, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
	  IMEXTimeIntegrationScheme>::AllocateSharedPtr(1, "");
        return p;
    }

    static std::string className;

}; // end class IMEXOrder1TimeIntegrationScheme


class IMEXOrder2TimeIntegrationScheme :
    public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder2TimeIntegrationScheme(int order, std::string type) :
        IMEXTimeIntegrationScheme(2, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
	  IMEXTimeIntegrationScheme>::AllocateSharedPtr(2, "");
        return p;
    }

    static std::string className;

}; // end class IMEXOrder2TimeIntegrationScheme


class IMEXOrder3TimeIntegrationScheme :
    public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder3TimeIntegrationScheme(int order, std::string type) :
        IMEXTimeIntegrationScheme(3, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
	  IMEXTimeIntegrationScheme>::AllocateSharedPtr(3, "");
        return p;
    }

    static std::string className;

}; // end class IMEXOrder3TimeIntegrationScheme


class IMEXOrder4TimeIntegrationScheme :
    public IMEXTimeIntegrationScheme
{
public:
    IMEXOrder4TimeIntegrationScheme(int order, std::string type) :
        IMEXTimeIntegrationScheme(4, "")
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);
    }

    static TimeIntegrationSchemeSharedPtr create(int order, std::string type)
    {
        boost::ignore_unused(order);
        boost::ignore_unused(type);

        TimeIntegrationSchemeSharedPtr p = MemoryManager<
	  IMEXTimeIntegrationScheme>::AllocateSharedPtr(4, "");
        return p;
    }

    static std::string className;

}; // end class IMEXOrder4TimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
