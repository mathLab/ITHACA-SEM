///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationSchemeData.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
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
// Description: Implementation of time integration scheme data class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeData.h>

#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>

#include <iostream>

#include <boost/core/ignore_unused.hpp>

#include <cmath>

namespace Nektar
{
namespace LibUtilities
{

////////////////////////////////////////////////////////////////////////////////
TimeIntegrationSchemeData::TimeIntegrationSolutionSharedPtr
    TimeIntegrationSchemeData::InitializeData(
        const NekDouble deltaT, ConstDoubleArray &y_0, const NekDouble time,
        const TimeIntegrationSchemeOperators &op)
{
    // create a TimeIntegrationSolution object based upon the
    // initial value. Initialise all other multi-step values
    // and derivatives to zero

    TimeIntegrationSolutionSharedPtr y_out =
        MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(this, y_0,
                                                                  time, deltaT);

    if (m_schemeType == eExplicit || m_schemeType == eExponential)
    {
        // ensure initial solution is in correct space
        op.DoProjection(y_0, y_out->UpdateSolution(), time);
    }

    // calculate the initial derivative, if is part of the
    // solution vector of the current scheme
    if (GetNmultiStepDerivs() > 0)
    {
        const Array<OneD, const unsigned int> offsets = GetTimeLevelOffset();

        if (offsets[GetNmultiStepValues()] == 0)
        {
            int i;
            int nvar    = y_0.size();
            int npoints = y_0[0].size();
            DoubleArray f_y_0(nvar);
            for (i = 0; i < nvar; i++)
            {
                f_y_0[i] = Array<OneD, NekDouble>(npoints);
            }
            // calculate the derivative of the initial value
            op.DoOdeRhs(y_0, f_y_0, time);

            // multiply by the step size
            for (i = 0; i < nvar; i++)
            {
                Blas::Dscal(npoints, deltaT, f_y_0[i].get(), 1);
            }
            y_out->SetDerivative(0, f_y_0, deltaT);
        }
    }

    return y_out;
}

TimeIntegrationScheme::ConstDoubleArray &TimeIntegrationSchemeData::
    TimeIntegrate(const NekDouble deltaT,
                  TimeIntegrationSolutionSharedPtr &solvector,
                  const TimeIntegrationSchemeOperators &op)
{
    // ASSERTL1( !(m_parent->GetIntegrationSchemeType() == eImplicit), "Fully
    // Implicit integration scheme cannot be handled by this routine." );

    int nvar    = solvector->GetFirstDim();
    int npoints = solvector->GetSecondDim();

    if (solvector->GetIntegrationSchemeData() !=
        this) // FIMXE ... verify that this is correct...
    {
        // This branch will be taken when the solution vector
        // (solvector) is set up for a different scheme than
        // the object this method is called from.  (typically
        // needed to calculate the first time-levels of a
        // multi-step scheme)

        // To do this kind of 'non-matching' integration, we
        // perform the following three steps:
        //
        // 1: copy the required input information from the
        //    solution vector of the master scheme to the
        //    input solution vector of the current scheme
        //
        // 2: time-integrate for one step using the current
        //    scheme
        //
        // 3: copy the information contained in the output
        //    vector of the current scheme to the solution
        //    vector of the master scheme

        // STEP 1: copy the required input information from
        //          the solution vector of the master scheme
        //          to the input solution vector of the
        //          current scheme

        // 1.1 Determine which information is required for the
        // current scheme

        DoubleArray y_n;
        NekDouble t_n = 0;
        DoubleArray dtFy_n;

        unsigned int nCurSchemeVals =
            GetNmultiStepValues(); // number of required values of the current
                                   // scheme
        unsigned int nCurSchemeDers =
            GetNmultiStepDerivs(); // number of required derivs of the current
                                   // scheme
        unsigned int nCurSchemeSteps =
            GetNsteps(); // number of steps in the current scheme    // FIMXE...
                         // what does it mean by "current scheme"... is this now
                         // a SchemeData issue?

        unsigned int nMasterSchemeVals =
            solvector->GetNvalues(); // number of values of the master scheme
        unsigned int nMasterSchemeDers =
            solvector->GetNderivs(); // number of derivs of the master scheme

        // The arrays below contains information to which time-level
        // the values and derivatives of the schemes belong
        const Array<OneD, const unsigned int> &curTimeLevels =
            GetTimeLevelOffset();
        const Array<OneD, const unsigned int> &masterTimeLevels =
            solvector->GetTimeLevelOffset();

        // 1.2 Copy the required information from the master
        //     solution vector to the input solution vector of
        //     the current scheme
        TimeIntegrationSolutionSharedPtr solvector_in =
            MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(
                this); // input solution vector of the current scheme

        for (int n = 0; n < nCurSchemeVals; n++)
        {
            // Get the required value out of the master solution vector
            // DoubleArray& y_n = solvector->GetValue    ( curTimeLevels[n] );
            // NekDouble    t_n = solvector->GetValueTime( curTimeLevels[n] );

            y_n = solvector->GetValue(curTimeLevels[n]);
            t_n = solvector->GetValueTime(curTimeLevels[n]);

            // Set the required value in the input solution
            // vector of the current scheme
            solvector_in->SetValue(curTimeLevels[n], y_n, t_n);
        }

        for (int n = nCurSchemeVals; n < nCurSchemeSteps; n++)
        {
            // Get the required derivative out of the master
            // solution vector
            // DoubleArray& dtFy_n = solvector->GetDerivative    (
            // curTimeLevels[n] );
            dtFy_n = solvector->GetDerivative(curTimeLevels[n]);

            // Set the required derivative in the input
            // solution vector of the current scheme
            solvector_in->SetDerivative(curTimeLevels[n], dtFy_n, deltaT);
        }

        // STEP 2: time-integrate for one step using the
        // current scheme
        TimeIntegrationSolutionSharedPtr solvector_out =
            MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(
                this, nvar,
                npoints); // output solution vector of the current scheme

        // Integrate one step
        TimeIntegrate(deltaT,
                      solvector_in->GetSolutionVector(),
                      solvector_in->GetTimeVector(),
                      solvector_out->UpdateSolutionVector(),
                      solvector_out->UpdateTimeVector(), op);

        // STEP 3: copy the information contained in the
        //         output vector of the current scheme to the
        //         solution vector of the master scheme

        // 3.1 Check whether the current time scheme updates
        //     the most recent derivative that should be
        //     updated in the master scheme.  If not,
        //     calculate the derivative. This can be done
        //     based upon the corresponding value and the
        //     DoOdeRhs operator.
        bool CalcNewDeriv = false; // This flag inidicates whether
                                   // the new derivative is availble
                                   // in the output of the current
                                   // scheme or whether it should be
                                   // calculated.

        if (nMasterSchemeDers > 0)
        {
            if (nCurSchemeDers == 0)
            {
                CalcNewDeriv = true;
            }
            else
            {
                if (masterTimeLevels[nMasterSchemeVals] <
                    curTimeLevels[nCurSchemeVals])
                {
                    CalcNewDeriv = true;
                }
            }
        }

        if (CalcNewDeriv)
        {
            int newDerivTimeLevel =
                masterTimeLevels[nMasterSchemeVals]; // Contains the
                                                     // time level at
                                                     // which the
                                                     // derivative of
                                                     // the master
                                                     // scheme is
                                                     // known.
            // DoubleArray  y_n;
            // NekDouble    t_n;
            // If the time level corresponds to 0, calculate the
            // derivative based upon the solution value at the new
            // time-level
            if (newDerivTimeLevel == 0)
            {
                y_n = solvector_out->GetValue(0);
                t_n = solvector_out->GetValueTime(0);
            }
            // If the time level corresponds to 1, calculate the
            // derivative based upon the solution value at the new
            // old-level.
            else if (newDerivTimeLevel == 1)
            {
                y_n = solvector->GetValue(0);
                t_n = solvector->GetValueTime(0);
            }
            else
            {
                ASSERTL1(false, "Problems with initialising scheme");
            }

            DoubleArray f_n(nvar);
            for (int j = 0; j < nvar; j++)
            {
                f_n[j] = Array<OneD, NekDouble>(npoints);
            }

            // Calculate the derivative
            op.DoOdeRhs(y_n, f_n, t_n);

            // Multiply by dt (as required by the General Linear Method
            // framework)
            for (int j = 0; j < nvar; j++)
            {
                Vmath::Smul(npoints, deltaT, f_n[j], 1, f_n[j], 1);
            }

            // Rotate the solution vector (i.e. updating without
            // calculating/inserting new values).
            solvector->RotateSolutionVector();

            // Set the calculated derivative in the master solution
            // vector.
            solvector->SetDerivative(newDerivTimeLevel, f_n, deltaT);
        }
        else
        {
            // Rotate the solution vector (i.e. updating without
            // calculating/inserting new values).
            solvector->RotateSolutionVector();
        }

        // 1.2 Copy the information calculated using the
        //     current scheme from the output solution vector
        //     to the master solution vector
        for (int n = 0; n < nCurSchemeVals; n++)
        {
            // Get the calculated value out of the output
            // solution vector of the current scheme
            // DoubleArray& y_n = solvector_out->GetValue    ( curTimeLevels[n]
            // );
            // NekDouble    t_n = solvector_out->GetValueTime( curTimeLevels[n]
            // );
            y_n = solvector_out->GetValue(curTimeLevels[n]);
            t_n = solvector_out->GetValueTime(curTimeLevels[n]);

            // Set the calculated value in the master solution vector
            solvector->SetValue(curTimeLevels[n], y_n, t_n);
        }

        for (int n = nCurSchemeVals; n < nCurSchemeSteps; n++)
        {
            // Get the calculated derivative out of the output
            // solution vector of the current scheme.
            // DoubleArray& dtFy_n =
            // solvector_out->GetDerivative (curTimeLevels[n]);
            dtFy_n = solvector_out->GetDerivative(curTimeLevels[n]);

            // Set the calculated derivative in the master
            // solution vector.
            solvector->SetDerivative(curTimeLevels[n], dtFy_n, deltaT);
        }
    }
    else
    {
        TimeIntegrationSolutionSharedPtr solvector_new =
            MemoryManager<TimeIntegrationSolution>::AllocateSharedPtr(
                this, nvar, npoints);

        TimeIntegrate(deltaT,
                      solvector->GetSolutionVector(),
                      solvector->GetTimeVector(),
                      solvector_new->UpdateSolutionVector(),
                      solvector_new->UpdateTimeVector(), op);

        solvector = solvector_new;
    }

    return solvector->GetSolution();

} // end TimeIntegrate()


// Does the actual multi-stage multi-step integration.
void
TimeIntegrationSchemeData::TimeIntegrate( const NekDouble deltaT,
                                          ConstTripleArray &y_old,
                                          ConstSingleArray &t_old,
                                          TripleArray &y_new,
                                          SingleArray &t_new,
                                          const TimeIntegrationSchemeOperators &op )
{
    ASSERTL1(
        CheckTimeIntegrateArguments(/*deltaT,*/ y_old, t_old, y_new, t_new, op),
        "Arguments not well defined");

    TimeIntegrationSchemeType type = m_schemeType;

    // Check if storage has already been initialised.
    // If so, we just zero the temporary storage.
    if (m_initialised && m_nvars == GetFirstDim(y_old) &&
        m_npoints == GetSecondDim(y_old))
    {
        for (int j = 0; j < m_nvars; j++)
        {
            Vmath::Zero(m_npoints, m_tmp[j], 1);
        }
    }
    else
    {
        m_nvars   = GetFirstDim(y_old);
        m_npoints = GetSecondDim(y_old);

        // First, calculate the various stage values and stage
        // derivatives (this is the multi-stage part of the method)
        // - m_Y   corresponds to the stage values
        // - m_F   corresponds to the stage derivatives
        // - m_T   corresponds to the time at the different stages
        // - m_tmp corresponds to the explicit right hand side of
        //   each stage equation
        //   (for explicit schemes, this correspond to m_Y)

        // Allocate memory for the arrays m_Y and m_F and m_tmp. The same
        // storage will be used for every stage -> m_Y is a
        // DoubleArray
        m_tmp = DoubleArray(m_nvars);
        for (int j = 0; j < m_nvars; j++)
        {
            m_tmp[j] = Array<OneD, NekDouble>(m_npoints, 0.0);
        }

        // The same storage will be used for every stage -> m_tmp is a
        // DoubleArray
        if (type == eExplicit || m_schemeType == eExponential)
        {
            m_Y = m_tmp;
        }
        else
        {
            m_Y = DoubleArray(m_nvars);
            for (int j = 0; j < m_nvars; j++)
            {
                m_Y[j] = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }

        // Different storage for every stage derivative as the data
        // will be re-used to update the solution -> m_F is a TripleArray
        m_F = TripleArray(m_numstages);
        for (int i = 0; i < m_numstages; ++i)
        {
            m_F[i] = DoubleArray(m_nvars);
            for (int j = 0; j < m_nvars; j++)
            {
                m_F[i][j] = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }

        if (type == eIMEX)
        {
            m_F_IMEX = TripleArray(m_numstages);
            for (int i = 0; i < m_numstages; ++i)
            {
                m_F_IMEX[i] = DoubleArray(m_nvars);
                for (int j = 0; j < m_nvars; j++)
                {
                    m_F_IMEX[i][j] = Array<OneD, NekDouble>(m_npoints, 0.0);
                }
            }
        }

        // Finally, flag that the memory has been initialised.
        m_initialised = true;
    } // end else

    // For an exponential integrator if the time increment or the
    // number of variables has changed then the exponenial matrices
    // must be recomputed.
    if (type == eExponential)
    {
        if (m_lastDeltaT != deltaT || m_lastNVars != GetFirstDim(y_old) )
        {
            m_parent->SetupSchemeExponentialData( this, deltaT );

            m_lastDeltaT = deltaT;
            m_lastNVars  = GetFirstDim(y_old);
        }
    }

    // The loop below calculates the stage values and derivatives
    for (int stage = 0; stage < m_numstages; stage++)
    {
        if ((stage == 0) && m_firstStageEqualsOldSolution)
        {
            for (int k = 0; k < m_nvars; k++)
            {
                Vmath::Vcopy(m_npoints, y_old[0][k], 1, m_Y[k], 1);
            }

            m_T = t_old[0];
        }
        else
        {
            // The stage values m_Y are a linear combination of:
            // 1: The stage derivatives:
            if (stage != 0)
            {
                for (int k = 0; k < m_nvars; k++)
                {
                    if (type == eExponential)
                    {
                        Vmath::Smul(m_npoints, deltaT * m_A_phi[k][stage][0],
                                    m_F[0][k], 1, m_tmp[k], 1);
                    }
                    else
                    {
                        Vmath::Smul(m_npoints, deltaT * A(stage, 0),
                                    m_F[0][k], 1, m_tmp[k], 1);
                    }

                    if (type == eIMEX)
                    {
                        Vmath::Svtvp(m_npoints, deltaT * A_IMEX(stage, 0),
                                     m_F_IMEX[0][k], 1, m_tmp[k], 1,
                                     m_tmp[k], 1);
                    }
                }
            }

            m_T = A(stage, 0) * deltaT;

            for (int j = 1; j < stage; j++)
            {
                for (int k = 0; k < m_nvars; k++)
                {
                    if (type == eExponential)
                    {
                        Vmath::Svtvp(m_npoints, deltaT * m_A_phi[k][stage][j],
                                     m_F[j][k], 1, m_tmp[k], 1, m_tmp[k], 1);
                    }
                    else
                    {
                        Vmath::Svtvp(m_npoints, deltaT * A(stage, j),
                                     m_F[j][k], 1, m_tmp[k], 1, m_tmp[k], 1);
                    }

                    if (type == eIMEX)
                    {
                        Vmath::Svtvp(m_npoints, deltaT * A_IMEX(stage, j),
                                     m_F_IMEX[j][k], 1, m_tmp[k], 1,
                                     m_tmp[k], 1);
                    }
                }

                m_T += A(stage, j) * deltaT;
            }

            // 2: The imported multi-step solution of the previous time level:
            for (int j = 0; j < m_numsteps; j++)
            {
                for (int k = 0; k < m_nvars; k++)
                {
                    if (type == eExponential)
                    {
                        Vmath::Svtvp(m_npoints, m_U_phi[k][stage][j],
                                     y_old[j][k], 1, m_tmp[k], 1, m_tmp[k], 1);
                    }
                    else
                    {
                        Vmath::Svtvp(m_npoints, U(stage, j),
                                     y_old[j][k], 1, m_tmp[k], 1, m_tmp[k], 1);
                    }
                }

                m_T += U(stage, j) * t_old[j];
            }
        } // end else

        // Calculate the stage derivative based upon the stage value
        if (type == eDiagonallyImplicit)
        {
            if (m_numstages == 1)
            {
                m_T = t_old[0] + deltaT;
            }
            else
            {
                m_T = t_old[0];
                for (int j = 0; j <= stage; ++j)
                {
                    m_T += A(stage, j) * deltaT;
                }
            }

            op.DoImplicitSolve(m_tmp, m_Y, m_T, A(stage, stage) * deltaT);

            for (int k = 0; k < m_nvars; k++)
            {
                Vmath::Vsub(m_npoints, m_Y[k], 1, m_tmp[k], 1, m_F[stage][k],
                            1);
                Vmath::Smul(m_npoints, 1.0 / (A(stage, stage) * deltaT),
                            m_F[stage][k], 1, m_F[stage][k], 1);
            }
        }
        else if (type == eIMEX)
        {
            if (m_numstages == 1)
            {
                m_T = t_old[0] + deltaT;
            }
            else
            {
                m_T = t_old[0];
                for (int j = 0; j <= stage; ++j)
                {
                    m_T += A(stage, j) * deltaT;
                }
            }

            if (fabs(A(stage, stage)) > NekConstants::kNekZeroTol)
            {
                op.DoImplicitSolve(m_tmp, m_Y, m_T, A(stage, stage) * deltaT);

                for (int k = 0; k < m_nvars; k++)
                {
                    Vmath::Vsub(m_npoints, m_Y[k], 1, m_tmp[k], 1,
                                m_F[stage][k], 1);
                    Vmath::Smul(m_npoints, 1.0 / (A(stage, stage) * deltaT),
                                m_F[stage][k], 1, m_F[stage][k], 1);
                }
            }
            op.DoOdeRhs(m_Y, m_F_IMEX[stage], m_T);
        }
        else if (type == eExplicit || m_schemeType == eExponential)
        {
            // Avoid projecting the same solution twice
            if (!((stage == 0) && m_firstStageEqualsOldSolution))
            {
                // ensure solution is in correct space
                op.DoProjection(m_Y, m_Y, m_T);
            }

            op.DoOdeRhs(m_Y, m_F[stage], m_T);
        }
    }

    // Next, the solution vector y at the new time level will
    // be calculated.
    //
    // For multi-step methods, this includes updating the
    // values of the auxiliary parameters
    //
    // The loop below calculates the solution at the new time
    // level
    //
    // If last stage equals the new solution, the new solution
    // needs not be calculated explicitly but can simply be
    // copied. This saves a solve.
    int i_start = 0;
    if (m_lastStageEqualsNewSolution)
    {
        for (int k = 0; k < m_nvars; k++)
        {
            Vmath::Vcopy(m_npoints, m_Y[k], 1, y_new[0][k], 1);
        }

        if (m_numstages == 1 && type == eIMEX)
        {
            t_new[0] = t_old[0] + deltaT;
        }
        else
        {
            t_new[0] = B(0, 0) * deltaT;

            for (int j = 1; j < m_numstages; j++)
            {
                t_new[0] += B(0, j) * deltaT;
            }

            for (int j = 0; j < m_numsteps; j++)
            {
                t_new[0] += V(0, j) * t_old[j];
            }
        }

        i_start = 1;
    }

    for (int i = i_start; i < m_numsteps; i++)
    {
        // The solution at the new time level is a linear
        // combination of:
        // 1: the stage derivatives
        for (int k = 0; k < m_nvars; k++)
        {
            if (type == eExponential)
            {
                Vmath::Smul(m_npoints, deltaT * m_B_phi[k][i][0],
                            m_F[0][k], 1, y_new[i][k], 1);
            }
            else
            {
                Vmath::Smul(m_npoints, deltaT * B(i, 0),
                            m_F[0][k], 1, y_new[i][k], 1);
            }

            if (type == eIMEX)
            {
                Vmath::Svtvp(m_npoints, deltaT * B_IMEX(i, 0),
                             m_F_IMEX[0][k], 1, y_new[i][k], 1, y_new[i][k], 1);
            }
        }

        if (m_numstages != 1 || type != eIMEX)
        {
            t_new[i] = B(i, 0) * deltaT;
        }

        for (int j = 1; j < m_numstages; j++)
        {
            for (int k = 0; k < m_nvars; k++)
            {
                if (type == eExponential)
                {
                    Vmath::Svtvp(m_npoints, deltaT * m_B_phi[k][i][j],
                                 m_F[j][k], 1, y_new[i][k], 1, y_new[i][k], 1);
                }
                else
                {
                    Vmath::Svtvp(m_npoints, deltaT * B(i, j),
                                 m_F[j][k], 1, y_new[i][k], 1, y_new[i][k], 1);
                }

                if (type == eIMEX)
                {
                    Vmath::Svtvp(m_npoints, deltaT * B_IMEX(i, j),
                                 m_F_IMEX[j][k], 1, y_new[i][k], 1,
                                 y_new[i][k], 1);
                }
            }

            if (m_numstages != 1 || type != eIMEX)
            {
                t_new[i] += B(i, j) * deltaT;
            }
        }

        // 2: the imported multi-step solution of the previous
        // time level
        for (int j = 0; j < m_numsteps; j++)
        {
            for (int k = 0; k < m_nvars; k++)
            {
                if (type == eExponential)
                {
                    Vmath::Svtvp(m_npoints, m_V_phi[k][i][j], y_old[j][k], 1,
                                 y_new[i][k], 1, y_new[i][k], 1);
                }
                else
                {
                    Vmath::Svtvp(m_npoints, V(i, j), y_old[j][k], 1,
                                 y_new[i][k], 1, y_new[i][k], 1);
                }
            }

            if (m_numstages != 1 || type != eIMEX)
            {
                t_new[i] += V(i, j) * t_old[j];
            }
        }
    }

    // Ensure that the new solution is projected if necessary
    if (type == eExplicit || m_schemeType == eExponential)
    {
        op.DoProjection(y_new[0], y_new[0], t_new[0]);
    }

} // end TimeIntegrate()

void TimeIntegrationSchemeData::CheckIfFirstStageEqualsOldSolution()
{
    // First stage equals old solution if:
    // 1. the first row of the coefficient matrix A consists of zeros
    // 2. U[0][0] is equal to one and all other first row entries of U are zero

    // 1. Check first condition
    for (int m = 0; m < m_A.size(); m++)
    {
        for (int i = 0; i < m_numstages; i++)
        {
            if (fabs(m_A[m][0][i]) > NekConstants::kNekZeroTol)
            {
                m_firstStageEqualsOldSolution = false;
                return;
            }
        }
    }

    // 2. Check second condition
    if (fabs(m_U[0][0] - 1.0) > NekConstants::kNekZeroTol)
    {
        m_firstStageEqualsOldSolution = false;
        return;
    }

    for (int i = 1; i < m_numsteps; i++)
    {
        if (fabs(m_U[0][i]) > NekConstants::kNekZeroTol)
        {
            m_firstStageEqualsOldSolution = false;
            return;
        }
    }

    m_firstStageEqualsOldSolution = true;
}

void TimeIntegrationSchemeData::CheckIfLastStageEqualsNewSolution()
{
    // Last stage equals new solution if:
    // 1. the last row of the coefficient matrix A is equal to the first row of
    // matrix B
    // 2. the last row of the coefficient matrix U is equal to the first row of
    // matrix V

    // 1. Check first condition
    for (int m = 0; m < m_A.size(); m++)
    {
        for (int i = 0; i < m_numstages; i++)
        {
            if (fabs(m_A[m][m_numstages-1][i] - m_B[m][0][i]) >
                NekConstants::kNekZeroTol)
            {
                m_lastStageEqualsNewSolution = false;
                return;
            }
        }
    }

    // 2. Check second condition
    for (int i = 0; i < m_numsteps; i++)
    {
        if (fabs(m_U[m_numstages-1][i] - m_V[0][i]) > NekConstants::kNekZeroTol)
        {
            m_lastStageEqualsNewSolution = false;
            return;
        }
    }

    m_lastStageEqualsNewSolution = true;
}

void TimeIntegrationSchemeData::VerifyIntegrationSchemeType()
{
#ifdef DEBUG
    int IMEXdim = m_A.size();
    int dim     = m_A[0].GetRows();

    Array<OneD, TimeIntegrationSchemeType> vertype(IMEXdim, eExplicit);

    if (m_schemeType == eExponential)
        vertype[0] = eExponential;

    for (int m = 0; m < IMEXdim; m++)
    {
        for (int i = 0; i < dim; i++)
        {
            if (fabs(m_A[m][i][i]) > NekConstants::kNekZeroTol)
            {
                vertype[m] = eDiagonallyImplicit;
            }
        }

        for (int i = 0; i < dim; i++)
        {
            for (int j = i + 1; j < dim; j++)
            {
                if (fabs(m_A[m][i][j]) > NekConstants::kNekZeroTol)
                {
                    vertype[m] = eImplicit;
                    ASSERTL1(false, "Fully Implicit schemes cannnot be handled "
                                    "by the TimeIntegrationScheme class");
                }
            }
        }
    }

    if (IMEXdim == 2)
    {
        ASSERTL1(m_B.size() == 2, "Coefficient Matrix B should have an "
                 "implicit and explicit part for IMEX schemes");

        if ((vertype[0] == eDiagonallyImplicit) && (vertype[1] == eExplicit))
        {
            vertype[0] = eIMEX;
        }
        else
        {
            ASSERTL1(false, "This is not a proper IMEX scheme");
        }
    }

    ASSERTL1(vertype[0] == m_schemeType,
                 "Time integration scheme coefficients do not match its type");
#endif
}

bool TimeIntegrationSchemeData::CheckTimeIntegrateArguments(
          ConstTripleArray &y_old,
          ConstSingleArray &t_old,
          TripleArray &y_new,
          SingleArray &t_new,
    const TimeIntegrationSchemeOperators &op) const
{
#ifdef DEBUG
    boost::ignore_unused(op);
#else
    boost::ignore_unused(y_old, t_old, y_new, t_new, op);
#endif

    // Check if arrays are all of consistent size

    ASSERTL1(y_old.size() == m_numsteps,
             "Non-matching number of steps.");
    ASSERTL1(y_new.size() == m_numsteps,
             "Non-matching number of steps.");

    ASSERTL1(y_old[0].size() == y_new[0].size(),
             "Non-matching number of variables.");
    ASSERTL1(y_old[0][0].size() == y_new[0][0].size(),
             "Non-matching number of coefficients.");

    ASSERTL1(t_old.size() == m_numsteps,
             "Non-matching number of steps.");
    ASSERTL1(t_new.size() == m_numsteps,
             "Non-matching number of steps.");

    return true;
}

std::ostream &operator<<(
    std::ostream &os,
    const TimeIntegrationSchemeData::TimeIntegrationSchemeDataSharedPtr &rhs)
{
    return operator<<(os, *rhs);
}

std::ostream &operator<<(std::ostream &os, const TimeIntegrationSchemeData &rhs)
{
    int r                          = rhs.m_numsteps;
    int s                          = rhs.m_numstages;
    TimeIntegrationSchemeType type = rhs.m_schemeType;

    int oswidth     = 9;
    int osprecision = 6;

    os << "Time Integration Scheme (Master): "
       << rhs.m_parent->GetFullName() << "\n"
       << "Time Integration Phase  : " << rhs.m_name << "\n"
       << "- number of steps:  " << r << "\n"
       << "- number of stages: " << s << "\n"
       << "General linear method tableau:\n";

    for (int i = 0; i < s; i++)
    {
        for (int j = 0; j < s; j++)
        {
            os.width(oswidth);
            os.precision(osprecision);
            os << std::right << rhs.A(i, j) << " ";
        }
        if (type == eIMEX)
        {
            os << " '";
            for (int j = 0; j < s; j++)
            {
                os.width(oswidth);
                os.precision(osprecision);
                os << std::right << rhs.A_IMEX(i, j) << " ";
            }
        }

        os << " |";

        for (int j = 0; j < r; j++)
        {
            os.width(oswidth);
            os.precision(osprecision);
            os << std::right << rhs.U(i, j);
        }
        os << std::endl;
    }

    int imexflag = (type == eIMEX) ? 2 : 1;
    for (int i = 0; i < (r + imexflag * s) * (oswidth + 1) + imexflag * 2 - 1;
         i++)
    {
        os << "-";
    }
    os << std::endl;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < s; j++)
        {
            os.width(oswidth);
            os.precision(osprecision);
            os << std::right << rhs.B(i, j) << " ";
        }
        if (type == eIMEX)
        {
            os << " '";
            for (int j = 0; j < s; j++)
            {
                os.width(oswidth);
                os.precision(osprecision);
                os << std::right << rhs.B_IMEX(i, j) << " ";
            }
        }

        os << " |";

        for (int j = 0; j < r; j++)
        {
            os.width(oswidth);
            os.precision(osprecision);
            os << std::right << rhs.V(i, j);
        }

        os << "  |";

        os.width(oswidth);
        os.precision(osprecision);
        os << std::right << rhs.m_timeLevelOffset[i];

        if( i < rhs.m_numMultiStepValues )
        {
          os << std::right << " value";
        }
        else
        {
          os << std::right << " derivative";
        }

        os << std::endl;
    }

    if( type == eExponential )
    {
        for (int k = 0; k < rhs.m_nvars; k++)
        {
            os << std::endl
               << "General linear method exponential tableau for variable "
               << k << ":\n";

            for (int i = 0; i < s; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << std::right << rhs.m_A_phi[k][i][j] << " ";
                }

                os << " |";

                for (int j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << std::right << rhs.m_U_phi[k][i][j];
                }
                os << std::endl;
            }

            int imexflag = (type == eIMEX) ? 2 : 1;
            for (int i = 0; i < (r + imexflag * s) * (oswidth + 1) + imexflag * 2 - 1;
                 i++)
            {
                os << "-";
            }
            os << std::endl;

            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << std::right << rhs.m_B_phi[k][i][j] << " ";
                }

                os << " |";

                for (int j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os.precision(osprecision);
                    os << std::right << rhs.m_V_phi[k][i][j];
                }
                os << std::endl;
            }
        }
    }

    return os;
} // end function operator<<

} // end namespace LibUtilities
} // end namespace NekTar
