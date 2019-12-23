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

#include <utilities/NekMesh/ProcessModules/ProcessVarOpti/Evaluator.hxx>

namespace Nektar
{

using namespace Utilities;

namespace LibUtilities
{

///////////////////////////////////////////////////////////////////////////////
//  EulerExponential

class EulerExponentialTimeIntegrationScheme : public TimeIntegrationScheme
{
public:
    EulerExponentialTimeIntegrationScheme() : TimeIntegrationScheme()
    {
        m_order = 3;

        m_integration_phases = TimeIntegrationSchemeDataVector(m_order);

        for( unsigned int n=0; n<m_order; ++n )
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

        // When multiple steps are taken B[0][0] and V[0][1...s] must be
        // weighted so the time contribution is correct.

        // B Phi function for first row first column
        phase->m_B[0][0][0] = 1.0 / phase->m_order; // phi_func(0 or 1)

        // B evaluation value shuffling first row second column
        if( phase->m_order > 1 )
          phase->m_B[0][0][1] = 1.0; // constant 1

        // U Curent time step evaluation first row first column
        phase->m_U[0][0] = 1.0; // constant 1

        // V Phi function for first row first column
        phase->m_V[0][0] = 1.0; // phi_func(0)

        // V Phi function for first row additional columns
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_V[0][n] = 1.0 / phase->m_order; // phi_func(n+1)
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
            NekDouble phi[phase->m_order];

            // B Phi function for first row first column
            if( GetName().find( "LawsonEuler" ) == 0 )
            {
                phi[0] =
                  phi_function(0, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
            }
            else if( GetName().find( "NorsettEuler" ) == 0 )
            {
                phi[0] =
                  phi_function(1, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
            }
            else
            {
              ASSERTL1(false,
                       "Cannot call EulerExponential directly "
                       "use LawsonEuler or NorsettEuler.");
            }

            // Set up for multiple steps. For multiple steps one needs
            // to weight the phi functions in much the same way there
            // are weights for multi-step Adams-Bashfort methods.

            // For order N the weights are an N x N matrix with
            // values: W[j][i] = std::pow(i, j) and phi_func = W phi.
            // There are other possible wieghting schemes
            if( phase->m_order == 1 )
            {
                // Nothing to do as the value is set above and the
                // wieght is just 1.
            }
            else if( phase->m_order == 2 )
            {
                // For Order 2 the weights are simply : 1 1
                //                                      0 1

                // If one were to solve the system of equations it
                // simply results in subtracting the second order
                // value from the first order value.
                phi[1] =
                    phi_function(2, deltaT, phase->m_L[0][k], phase->m_L[1][k]);

                phi[0] -= phi[1];
            }
            else if( phase->m_order == 3 )
            {
                NekDouble phi_func[phase->m_order];

                phi_func[0] = phi[0];
                phi_func[1] =
                    phi_function(2, deltaT, phase->m_L[0][k], phase->m_L[1][k]);
                phi_func[2] =
                    phi_function(3, deltaT, phase->m_L[0][k], phase->m_L[1][k]);

                NekDouble W[3][3];

                // Set up teh wieghts and calculate the determinant.
                // for( unsigned int j=0; j<phase->m_order; ++j )
                // {
                //     for( unsigned int i=0; i<phase->m_order; ++i )
                //     {
                //      W[j][i] = std::pow(i, j);
                //     }
                // }
                // NekDouble W_det = Determinant<3>(W);

                // No need to calculate the determinant value as it is fixed;
                NekDouble W_det = 2.0;

                // Solve the series of equations using Cramer's rule.
                for( unsigned int m=0; m<phase->m_order; ++m )
                {
                    // Assemble the working matrix for this solution.
                    for( unsigned int j=0; j<phase->m_order; ++j )
                    {
                        for( unsigned int i=0; i<phase->m_order; ++i )
                        {
                            // Fill in the mth column for the mth
                            // solution using the phi function value
                            // otherwise utilize the weights.
                            W[i][j] = (j == m) ? phi_func[i] : std::pow(j, i);
                        }
                    }

                    // Get the mth solutiion.
                    phi[m] = Determinant<3>(W) / W_det;
                }
            }
            else
            {
                ASSERTL1(false, "Not set up for more than 3rd Order.");
            }

            // Create the phi based Butcher tableau matrices. Note
            // these matrices are set up using a general formational
            // based on the number of steps (i.e. the order).
            phase->m_A_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numstages, 0.0);
            phase->m_B_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numstages, 0.0);
            phase->m_U_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numstages, phase->m_numsteps, 0.0);
            phase->m_V_phi[k] =
              Array<TwoD, NekDouble>(phase->m_numsteps, phase->m_numsteps, 0.0);

            // B Phi function for first row first column.
            phase->m_B_phi[k][0][0] = phi[0];

            // B evaluation value shuffling first row second column.
            if( phase->m_order > 1 )
                phase->m_B_phi[k][0][1] = 1.0; // constant 1

            // U Curent time step evaluation first row first column.
            phase->m_U_phi[k][0][0] = 1.0; // constant 1

            // V Phi function for first row first column.
            phase->m_V_phi[k][0][0] =
                phi_function(0, deltaT, phase->m_L[0][k], phase->m_L[1][k]);

            // V Phi function for first row additional columns.
            for( int n=1; n<phase->m_order; ++n )
            {
                phase->m_V_phi[k][0][n] = phi[n];
            }

            // V evaluation value shuffling row n column n-1.
            for( int n=2; n<phase->m_order; ++n )
            {
                phase->m_V_phi[k][n][n-1] = 1.0; // constant 1
            }
        }

        std::cout << *phase << std::endl;
    }

    protected:
        unsigned int m_order;

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
        return std::string("LawsonEuler" + std::to_string(m_order));
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
        return std::string("NorsettEuler" + std::to_string(m_order));
    }
}; // end class NorsettEulerTimeIntegrationScheme
  
} // end namespace LibUtilities
} // end namespace Nektar
