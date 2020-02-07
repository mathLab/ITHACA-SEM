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
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

namespace Nektar
{
namespace LibUtilities
{

////////////////////////////////////////////////////////////////////////////////
// Runge Kutta Order N where the number of stages == order

class RungeKuttaTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    RungeKuttaTimeIntegrationScheme(std::string variant, unsigned int order,
				    std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme(variant, order, freeParams)
    {
        ASSERTL1((variant == "" || variant == "Classic" ||
		  variant == "SSP" || variant == "ImprovedEuler"),
                 "Runge Kutta Time integration scheme bad variant: "
		 + variant + ". "
                  "Must blank, 'Classic', 'SSP', or 'ImprovedEuler'");

        // Std - Currently up to 5th order is implemented.
        // SSP - Currently 1st through 3rd order is implemented.
        ASSERTL1(((variant == "SSP" || variant == "ImprovedEuler") == 0 &&
                  1 <= order && order <= 5) ||
                 ((variant == "SSP" || variant == "ImprovedEuler") == 1 &&
                  1 <= order && order <= 3),
                 "Runge Kutta Time integration scheme bad order "
                 "Std (1-5) or SSP (1-3): " + std::to_string(order));

        m_integration_phases    = TimeIntegrationSchemeDataVector(1);
        m_integration_phases[0] = TimeIntegrationSchemeDataSharedPtr(
            new TimeIntegrationSchemeData(this));

        RungeKuttaTimeIntegrationScheme::SetupSchemeData(
            m_integration_phases[0], variant, order, freeParams);
    }

    virtual ~RungeKuttaTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
						 std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
	    RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("RungeKutta");
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        if( GetOrder() == 4 || GetOrder() == 5 )
        {
            return 2.784;
        }
        else
        {
            return 2.0;
        }
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    std::string variant, unsigned int order,
				    std::vector<NekDouble> freeParams)
    {
        boost::ignore_unused(freeParams);

        const unsigned int nStages[6] = { 0, 1, 2, 3, 4, 6 };

        // A Coefficients for the lower diagonal quadrant stored in a
        // contiguous fashion. For the fourth order, six coefficients
        // from the Butcher tableau would be stored as the following.
        //
        //                0 0 0 0 
        //    Butcher     a 0 0 0   Stored as   a
        //    Tableau     b c 0 0               b c
        //                d e f 0               d e f 0 ... 0
        const NekDouble Acoefficients[2][6][15] =
            { { {     0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 1st Order
                {     0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 2nd Order - midpoint
                {   1./2,      // Last entry
                      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 3rd Order - Ralston's
                {  1./2.,
                      0.,   3./4.,      // Last entry
                      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 4th Order - Classic
                {  1./2.,
                      0.,   1./2.,
                      0.,      0.,      1.,      // Last entry
                      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 5th Order - 6 stages
                {  1./4.,
                   1./8.,   1./8.,
                      0.,      0.,    1./2.,
                  3./16.,  -3./8.,    3./8.,  9./16.,
                  -3./7.,   8./7.,    6./7., -12./7.,   8./7. } },
              // Strong Stability Preserving
              { {     0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 1st Order
                {     0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 2nd Order - strong scaling - improved
                {     1.,      // Last entry
                      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 3rd Order - strong scaling
                {     1.,
                      1./4.,   1./4.,      // Last entry
                      0,       0.,
                      0.,      0.,       0.,     0.,      0.,
                      0.,      0.,       0.,     0.,      0. },
                // 4th Order - Classic - not used
                {  1./2.,
                      0.,   1./2.,
                      0.,      0.,      1.,      // Last entry
                      0.,      0.,      0.,      0.,
                      0.,      0.,      0.,      0.,      0. },
                // 5th Order - 6 stages - not used
                {  1./4.,
                   1./8.,   1./8.,
                      0.,      0.,    1./2.,
                  3./16.,  -3./8.,    3./8.,  9./16.,
                  -3./7.,   8./7.,    6./7., -12./7.,   8./7. } } };

        // B Coefficients for the finial summing.
        const NekDouble Bcoefficients[2][6][5] =
            { { {    0.,       0.,    0.,       0.,      0. },
                // 1st Order
                {    1.,       0.,    0.,       0.,      0. },
                // 2nd Order - midpoint
                {    0.,       1.,     0.,      0.,      0. },
                // 3rd Order - Ralston's
                { 2./9.,    3./9.,  4./9.,      0.,      0. },
                // 4th Order - Classic
                { 1./6.,    2./6.,  2./6.,   1./6.,      0. },
                // 5th Order - 6 stages
                { 7./90., 32./90., 12./90., 32./90., 7./90. } },
              // Strong Stability Preserving
              { {    0.,       0.,     0.,      0.,      0. },
                // 1st Order
                {    1.,       0.,     0.,      0.,      0. },
                // 2nd Order - improved
                { 1./2.,    1./2.,     0.,      0.,      0. },
                // 3rd Order - strong scaling
                { 1./6.,    1./6.,  4./6.,      0.,      0. },
                // 4th Order - Classic
                { 1./6.,    2./6.,  2./6.,   1./6.,      0. },
                // 5th Order - 6 stages
                { 7./90., 32./90., 12./90., 32./90., 7./90. } } };

        unsigned int index = (variant == "SSP" || variant == "ImprovedEuler");

        phase->m_schemeType = eExplicit;
        phase->m_variant = variant;
        phase->m_order = order;
        phase->m_name =
          std::string("RungeKutta") + phase->m_variant +
          std::string("Order") + std::to_string(phase->m_order);

        phase->m_numsteps  = 1;
        phase->m_numstages = nStages[phase->m_order];

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

        // Coefficients

        // A Coefficients for each stages along the lower diagonal quadrant.
        unsigned int cc = 0;

        for( int s=1; s<phase->m_numstages; ++s )
        {
            for( int i=0; i<s; ++i )
            {
                phase->m_A[0][s][i] =
                  Acoefficients[index][phase->m_order][cc++];
            }
        }

        // B Coefficients for the finial summing.
        for( int n=0; n<phase->m_order; ++n )
        {
            phase->m_B[0][0][n] = Bcoefficients[index][phase->m_order][n];
        }

        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = 0;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        phase->CheckAndVerify();
    }

}; // end class RungeKuttaTimeIntegrator

////////////////////////////////////////////////////////////////////////////////
// Backwards compatibility
class RungeKutta2TimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    RungeKutta2TimeIntegrationScheme(std::string variant, unsigned int order,
				     std::vector<NekDouble> freeParams) :
        RungeKuttaTimeIntegrationScheme("", 2, freeParams)
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
          RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("", 2, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta2TimeIntegrationScheme

class ClassicalRungeKutta4TimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    ClassicalRungeKutta4TimeIntegrationScheme(std::string variant, unsigned int order,
					      std::vector<NekDouble> freeParams) :
        RungeKuttaTimeIntegrationScheme("Classic", 4, freeParams)
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
          RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("Classic", 4, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta2TimeIntegrationScheme

class RungeKutta5TimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    RungeKutta5TimeIntegrationScheme(std::string variant, unsigned int order,
				     std::vector<NekDouble> freeParams) :
      RungeKuttaTimeIntegrationScheme("", 5, freeParams)
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
          RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("", 5, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta2TimeIntegrationScheme

class RungeKutta2_ImprovedEulerTimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    RungeKutta2_ImprovedEulerTimeIntegrationScheme(std::string variant, unsigned int order,
						   std::vector<NekDouble> freeParams) :
        RungeKuttaTimeIntegrationScheme("ImprovedEuler", 2, freeParams)
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
            RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("ImprovedEuler", 2, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta2_ImprovedEulerTimeIntegrationScheme

class RungeKutta2_SSPTimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    RungeKutta2_SSPTimeIntegrationScheme(std::string variant, unsigned int order,
					 std::vector<NekDouble> freeParams) :
      RungeKuttaTimeIntegrationScheme("SSP", 2, freeParams)
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
          RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("SSP", 2, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta2_SSPTimeIntegrationScheme

class RungeKutta3_SSPTimeIntegrationScheme :
    public RungeKuttaTimeIntegrationScheme
{
public:
    RungeKutta3_SSPTimeIntegrationScheme(std::string variant, unsigned int order,
					 std::vector<NekDouble> freeParams) :
      RungeKuttaTimeIntegrationScheme("SSP", 3, freeParams)
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
          RungeKuttaTimeIntegrationScheme>::AllocateSharedPtr("SSP", 3, freeParams);
        return p;
    }

    static std::string className;

}; // end class RungeKutta3_SSPTimeIntegrationScheme

} // end namespace LibUtilities
} // end namespace Nektar
