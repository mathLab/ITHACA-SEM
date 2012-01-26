///////////////////////////////////////////////////////////////////////////////
//
// File Monodomain.cpp
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
// Description: Monodomain cardiac electrophysiology homogenised model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <CardiacEPSolver/EquationSystems/Monodomain.h>

namespace Nektar
{
    /**
     * @class Monodomain
     *
     * Base model of cardiac electrophysiology of the form
     * \f{align*}{
     *     \frac{\partial u}{\partial t} = \nabla^2 u + J_{ion},
     * \f}
     * where the reaction term, \f$J_{ion}\f$ is defined by a specific cell
     * model.
     *
     * This implementation, at present, treats the reaction terms explicitly
     * and the diffusive element implicitly.
     */

    /**
     * Registers the class with the Factory.
     */
    string Monodomain::className
            = GetEquationSystemFactory().RegisterCreatorFunction(
                "Monodomain",
                Monodomain::create,
                "Phenomological model of canine cardiac electrophysiology.");


    /**
     *
     */
    Monodomain::Monodomain(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    void Monodomain::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("epsilon",    m_epsilon,   1.0);

        std::string vCellModel;
        m_session->LoadSolverInfo("CELLMODEL", vCellModel, "");

        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session, m_fields[0]->GetNpoints());

        // Load variable coefficients
        /// @todo: move to an EvaluateFunction routine.
        std::string varCoeffs[3] = {"d00", "d11", "d22"};
        StdRegions::VarCoeffType varCoeffEnum[3] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD22
        };
        std::string varName = "intensity";

        for (int i = 0; i < 3; ++i)
        {
            if (m_session->DefinesFunction(varCoeffs[i]))
            {
                // Load from FLD file.
                if (m_session->GetFunctionType(varCoeffs[i]) == LibUtilities::eFunctionTypeFile)
                {
                    ImportFld(m_session->GetFunctionFilename(varCoeffs[i]),
                              m_fields[0],
                              varName);

                    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(),
                                          m_fields[0]->UpdatePhys());

                    // Normalise and invert (assuming image intensity data)
                    int nq = m_fields[0]->GetNpoints();
                    NekDouble max = Vmath::Vmax(nq, m_fields[0]->GetPhys(), 1);
                    Vmath::Smul(nq, 1.0/max, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
                    Vmath::Neg(nq, m_fields[0]->UpdatePhys(), 1);
                    Vmath::Sadd(nq, 1.0, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);

                    m_vardiff[varCoeffEnum[i]] = m_fields[0]->GetPhys();
                }
                // Evaluate expression
                else
                {
                    int nq = m_fields[0]->GetNpoints();
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);

                    // get the coordinates
                    m_fields[0]->GetCoords(x0,x1,x2);

                    Array<OneD, NekDouble> tmp(nq);

                    LibUtilities::EquationSharedPtr ifunc
                            = m_session->GetFunction(varCoeffs[i], varName);
                    for(int j = 0; j < nq; j++)
                    {
                        tmp[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],0.0);
                    }
                    m_vardiff[varCoeffEnum[i]] = tmp;
                }

                // Dump actual variable coefficients for verification.
                m_fields[0]->FwdTrans_IterPerExp(m_fields[0]->GetPhys(),
                                                 m_fields[0]->UpdateCoeffs());
                WriteFld(varCoeffs[i] + ".fld");
            }
        }

        if (m_session->DefinesParameter("StimulusDuration"))
        {
            ASSERTL0(m_session->DefinesFunction("Stimulus", "u"),
                    "Stimulus function not defined.");
            m_session->LoadParameter("StimulusDuration", m_stimDuration);
        }
        else
        {
            m_stimDuration = 0;
        }

        if (!m_explicitDiffusion)
        {
            m_ode.DefineImplicitSolve (&Monodomain::DoImplicitSolve, this);
        }
        m_ode.DefineOdeRhs(&Monodomain::DoOdeRhs, this);
    }


    /**
     *
     */
    Monodomain::~Monodomain()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array.
     * @param   time            Current simulation time.
     * @param   lambda          Timestep.
     */
    void Monodomain::DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables  = inarray.num_elements();
        int ncoeffs     = inarray[0].num_elements();
        int nq          = m_fields[0]->GetNpoints();
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 1.0/lambda/m_epsilon;

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Only apply diffusion to first variable.
            if (i > 0) {
                Vmath::Vcopy(nq, &inarray[i][0], 1, &outarray[i][0], 1);
                continue;
            }

            // Multiply 1.0/timestep/lambda
            Vmath::Smul(nq, -1.0/lambda/m_epsilon, inarray[i], 1,
                                            m_fields[i]->UpdatePhys(), 1);

            // Solve a system of equations with Helmholtz solver and transform
            // back into physical space.
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(), NullFlagList,
                                   factors, m_vardiff);

            m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(),
                                   m_fields[i]->UpdatePhys());
            m_fields[i]->SetPhysState(true);

            // Copy the solution vector (required as m_fields must be set).
            outarray[i] = m_fields[i]->GetPhys();
        }
    }


    void Monodomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        if (m_stimDuration > 0 && time < m_stimDuration)
        {
            int nq = m_fields[0]->GetNpoints();
            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);

            // get the coordinates
            m_fields[0]->GetCoords(x0,x1,x2);

            LibUtilities::EquationSharedPtr ifunc
                    = m_session->GetFunction("Stimulus", "u");
            for(int j = 0; j < nq; j++)
            {
                outarray[0][j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
            }
        }
        m_cell->Update(inarray, outarray, time);
    }


    void Monodomain::v_SetInitialConditions(NekDouble initialtime,
                        bool dumpInitialConditions)
    {
        EquationSystem::v_SetInitialConditions(initialtime, dumpInitialConditions);
    }


    /**
     *
     */
    void Monodomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        out << "\tEpsilon         : " << m_epsilon << endl;
        if (m_session->DefinesParameter("d00"))
        {
            out << "\tDiffusivity-x   : " << m_session->GetParameter("d00") << endl;
        }
        m_cell->PrintSummary(out);
    }

}
