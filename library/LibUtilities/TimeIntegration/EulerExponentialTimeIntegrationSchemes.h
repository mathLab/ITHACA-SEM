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
    EulerExponentialTimeIntegrationScheme(std::string variant, unsigned int order,
                                          std::vector<NekDouble> freeParams) :
        TimeIntegrationScheme(variant, order, freeParams)
    {
        // Currently up to 4th order is implemented because the number
        // of steps is the same as the order.
        // Currently up to 4th order is implemented.
        ASSERTL1(variant == "Lawson" || variant == "Norsett",
                 "EulerExponential Time integration scheme bad variant: " +
                 variant);

        ASSERTL1(0 < order && order <= 4,
                 "EulerExponential Time integration scheme bad order: " +
                 std::to_string(order));

        m_integration_phases = TimeIntegrationSchemeDataVector(order);

        // Currently the next lowest order is used to seed the current
        // order. This is not correct but is an okay approximation.
        for( unsigned int n=0; n<order; ++n )
        {
            m_integration_phases[n] = TimeIntegrationSchemeDataSharedPtr(
                new TimeIntegrationSchemeData(this));

            EulerExponentialTimeIntegrationScheme::SetupSchemeData(
                m_integration_phases[n], variant, n+1);
        }

        // for( unsigned int n=0; n<order; ++n )
        //   std::cout << m_integration_phases[n] << std::endl;
    }

    virtual ~EulerExponentialTimeIntegrationScheme()
    {
    }

    static TimeIntegrationSchemeSharedPtr create(std::string variant, unsigned int order,
                                                 std::vector<NekDouble> freeParams)
    {
        TimeIntegrationSchemeSharedPtr p = MemoryManager<
            EulerExponentialTimeIntegrationScheme>::AllocateSharedPtr(variant, order, freeParams);

        return p;
    }

    static std::string className;

    LUE virtual std::string GetName() const
    {
        return std::string("EulerExponential");
    }

    LUE virtual std::string GetFullName () const
    {
        return GetVariant() + GetName() + "Order" + std::to_string(GetOrder());
    }

    LUE virtual NekDouble GetTimeStability() const
    {
        return 1.0;
    }

    LUE static void SetupSchemeData(TimeIntegrationSchemeDataSharedPtr &phase,
                                    std::string variant, int order)
    {
        phase->m_schemeType = eExponential;
        phase->m_variant  = variant;
        phase->m_order = order;
        phase->m_name = variant + std::string("EulerExponentialOrder") +
          std::to_string(phase->m_order);

        // Parameters for the compact 1 step implementation.
        phase->m_numstages = 1;
        phase->m_numsteps  = phase->m_order;  // Okay up to 4th order.

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

        // B Phi function for first row first column
        phase->m_B[0][0][0] = 1.0 / phase->m_order; // phi_func(0 or 1)

        // B evaluation value shuffling second row first column
        if( phase->m_order > 1 )
          phase->m_B[0][1][0] = 1.0; // constant 1

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

        // For order > 1 derivatives are needed.
        for( int n=1; n<phase->m_order; ++n )
        {
            phase->m_timeLevelOffset[n] = n;
        }

        phase->CheckAndVerify();
    }

    virtual void SetupSchemeExponentialData(TimeIntegrationSchemeData *phase,
                                            NekDouble deltaT) const
    {
        // Assumptions the two-dimensional Lambda matrix is a diagonal
        // matrix thus values are non zero if and only i=j. As such,
        // the diagonal Lambda values are stored as two vectors so to
        // accomodate complex numbers m_L[0] real, m_L[1] imaginary.
        ASSERTL1(phase->m_nvars == phase->m_L.size(),
                 "The number of variables does not match "
                 "the number of exponential coefficents.");

        phase->m_A_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_B_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_U_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);
        phase->m_V_phi = Array<OneD, Array<TwoD, NekDouble>>(phase->m_nvars);

        for( unsigned int k=0; k<phase->m_nvars; ++k )
        {
            Array<OneD, NekDouble> phi = Array<OneD, NekDouble>(phase->m_nvars);

            // B Phi function for first row first column
            if( phase->m_variant == "Lawson" )
            {
	        phi[0] = phi_function(0, deltaT * phase->m_L[k]).real();
            }
            else if( phase->m_variant == "Norsett" )
            {
                phi[0] = phi_function(1, deltaT * phase->m_L[k]).real();
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
            // There are other possible wieghting schemes.
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
                phi[1] = phi_function(2, deltaT * phase->m_L[k]).real();

                phi[0] -= phi[1];
            }
            else if( phase->m_order == 3 )
            {
                Array<OneD, NekDouble> phi_func =
                  Array<OneD, NekDouble>(phase->m_nvars);

                phi_func[0] = phi[0];

                for( unsigned int m=1; m<phase->m_order; ++m )
                {
                    phi_func[m] = phi_function(m+1, deltaT * phase->m_L[k]).real();
                }

                NekDouble W[3][3];

                // Set up the wieghts and calculate the determinant.
                for( unsigned int j=0; j<phase->m_order; ++j )
                {
                    for( unsigned int i=0; i<phase->m_order; ++i )
                    {
                        W[j][i] = std::pow(i, j);
                    }
                }

                NekDouble W_det = Determinant<3>(W);

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
            else if( phase->m_order == 4 )
            {
                Array<OneD, NekDouble> phi_func =
                  Array<OneD, NekDouble>(phase->m_nvars);

                phi_func[0] = phi[0];

                for( unsigned int m=1; m<phase->m_order; ++m )
                {
                    phi_func[m] = phi_function(m+1, deltaT * phase->m_L[k]).real();
                }

                NekDouble W[4][4];

                // Set up the wieghts and calculate the determinant.
                for( unsigned int j=0; j<phase->m_order; ++j )
                {
                    for( unsigned int i=0; i<phase->m_order; ++i )
                    {
                        W[j][i] = std::pow(i, j);
                    }
                }

                NekDouble W_det = Determinant<4>(W);

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
                    phi[m] = Determinant<4>(W) / W_det;
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

            // B evaluation value shuffling second row first column.
            if( phase->m_order > 1 )
                phase->m_B_phi[k][1][0] = 1.0; // constant 1

            // U Curent time step evaluation first row first column.
            phase->m_U_phi[k][0][0] = 1.0; // constant 1

            // V Phi function for first row first column.
            phase->m_V_phi[k][0][0] = phi_function(0, deltaT * phase->m_L[k]).real();

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

	// std::cout << *phase << std::endl;
    }

}; // end class EulerExponentialTimeIntegrator

} // end namespace LibUtilities
} // end namespace Nektar
