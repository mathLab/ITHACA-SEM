///////////////////////////////////////////////////////////////////////////////
//
// File: EulerExponentialTimeIntegrationSchemes.h
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
// Description: Combined header file for all basic Euler exponential
// based time integration schemes.
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

///////////////////////////////////////////////////////////////////////////////
//  EulerExponential

class EulerExponentialTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    EulerExponentialTimeIntegrationScheme() : TimeIntegrationScheme()
    {
        unsigned int order = 2;

        m_integration_phases = TimeIntegrationSchemeDataVector(order);

        for( unsigned int n=0; n<order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));

            m_integration_phases[n]->m_order = n+1;

            EulerExponentialTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[n]);
        }
    }

    virtual ~EulerExponentialTimeIntegrationScheme()
    {
    }

    LUE virtual std::string GetName() const = 0;

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase)
    {
        phase->m_schemeType = eExponential;

        // Parameters for the compact 1 step implementation.
        phase->m_numstages = 1;
        phase->m_numsteps  = phase->m_order;

        phase->m_A = Array<OneD, Array<TwoD, NekDouble>>(1);
        phase->m_B = Array<OneD, Array<TwoD, NekDouble>>(1);

        phase->m_A[0] =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
        phase->m_B[0] =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);

        phase->m_U =
            Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
        phase->m_V =
            Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

        // Coefficients

        // B Phi function for first row first column
        phase->m_B[0][0][0] = 1.0 / phase->m_order; // phi(0) or phi(1)

        // B evaluation value shuffling first row second column
        if( phase->m_order > 1 )
          phase->m_B[0][0][1] = 1.0; // constant 1

        // U Curent time step evaluation first row first column
        phase->m_U[0][0] = 1.0; // constant 1

        // V Phi function for first row first column
        phase->m_V[0][0] = 1.0; // phi(0)

        // V Phi function for first row additional columns
        for( int n=2; n<=phase->m_order; ++n )
        {
            phase->m_V[0][n-1] = 1.0 / phase->m_order; // phi(n+1)
        }

        // V evaluation value shuffling row n column n-1
        for( int n=2; n<phase->m_order; ++n )
        {
            phase->m_V[n][n-1] = 1.0; // constant 1
        }
        
        phase->m_numMultiStepValues = 1;
        phase->m_numMultiStepDerivs = phase->m_order-1;
        phase->m_timeLevelOffset = Array<OneD, unsigned int>(phase->m_numsteps);
        phase->m_timeLevelOffset[0] = 0;

        // For order > 1 then derivatives are needed.
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_timeLevelOffset[n] = n;
        }

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

    virtual void SetupSchemeExponentialData(TimeIntegrationSchemeData *phase,
                                            NekDouble deltaT) const
    {
        // Assumptions the two-dimensional Lambda matrix is a diagonal
        // matrix thus values are non zero if and only i=j. As such,
        // the diagonal Lambda values are stored as two vectors so to
        // accomodate complex numbers m_L[0] real, m_L[1] imaginary.
        ASSERTL1(phase->m_nvars == phase->m_L.GetColumns(),
                 "The number of variables does not match "
                 "the number of exponential coefficents.");

        phase->m_A_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_B_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_U_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_V_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        
        for( unsigned int k=0; k<phase->m_nvars; ++k )
        {
            phase->m_A_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
            phase->m_B_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
            phase->m_U_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
            phase->m_V_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

            // B Phi function for first row first column
            if( GetName() == "LawsonEuler" )
            {
                phase->m_B_phi[k][0][0] =
                    phi_function(0, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
            }
            else if( GetName() == "NorsettEuler" )
            {
                phase->m_B_phi[k][0][0] =
                    phi_function(1, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
                // FIXME - Hack for order == 2 ONLY
                if( phase->m_order == 2 )
                  phase->m_B_phi[k][0][0] -=
                    phi_function(2, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
            }
            else
            {
              ASSERTL1(false,
                       "Cannot call EulerExponential directly "
                       "use LawsonEuler or NorsettEuler.");
            }

            // B evaluation value shuffling first row second column
            if( phase->m_order > 1 )
                phase->m_B_phi[k][0][1] = 1.0; // constant 1

            // U Curent time step evaluation first row first column
            phase->m_U_phi[k][0][0] = 1.0; // constant 1

            // V Phi function for first row first column
            phase->m_V_phi[k][0][0] =
                phi_function(0, deltaT, phase->m_L[0][k], phase->m_L[1][k]);

            // V Phi function for first row additional columns
            for( int n=2; n<=phase->m_order; ++n )
            {
                phase->m_V_phi[k][0][n-1] =
                    phi_function(n, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
            }

            // V evaluation value shuffling row n column n-1
            for( int n=2; n<phase->m_order; ++n )
            {
                phase->m_V_phi[k][n][n-1] = 1.0; // constant 1
            }
        }

        std::cout << *phase << std::endl;
    }
}; // end class EulerExponentialTimeIntegrator

class LawsonEulerTimeIntegrationScheme : public EulerExponentialTimeIntegrationScheme
{
public:
    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            LawsonEulerTimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("LawsonEuler");
    }
}; // end class LawsonEulerTimeIntegrationScheme
  
class NorsettEulerTimeIntegrationScheme : public EulerExponentialTimeIntegrationScheme
{
public:
    static TimeIntegrationSchemeSharedPtr create()
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            NorsettEulerTimeIntegrationScheme>::AllocateSharedPtr();
        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("NorsettEuler");
    }
}; // end class NorsettEulerTimeIntegrationScheme
  
} // end namespace LibUtilities
} // end namespace Nektar
