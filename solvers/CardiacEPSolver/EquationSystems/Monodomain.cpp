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
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

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
                "Monodomain model of cardiac electrophysiology.");


    /**
     *
     */
    Monodomain::Monodomain(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }


    /**
     *
     */
    void Monodomain::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("Chi",        m_chi);
        m_session->LoadParameter("Cm",         m_capMembrane);

        std::string vCellModel;
        m_session->LoadSolverInfo("CELLMODEL", vCellModel, "");

        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = GetCellModelFactory().CreateInstance(
                                        vCellModel, m_session, m_fields[0]);

        m_intVariables.push_back(0);

        // Load variable coefficients
        StdRegions::VarCoeffType varCoeffEnum[3] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD22
        };
        std::string varName = "intensity";
        std::string varCoeffs[3] = {
                "AnisotropicConductivityX",
                "AnisotropicConductivityY",
                "AnisotropicConductivityZ"
        };
        int nq = m_fields[0]->GetNpoints();
        Array<OneD, NekDouble> vTemp;

        // Allocate storage for variable coeffs and initialize to 1.
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vardiff[varCoeffEnum[i]] = Array<OneD, NekDouble>(nq, 1.0);
        }

        // Apply intensity map (range d_min -> d_max)
        if (m_session->DefinesFunction("IsotropicConductivity"))
        {
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                cout << "Loading Isotropic Conductivity map." << endl;
            }
            EvaluateFunction(varName, vTemp, "IsotropicConductivity");
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, vTemp, 1, m_vardiff[varCoeffEnum[i]], 1, m_vardiff[varCoeffEnum[i]], 1);
            }
        }

        // Apply fibre map (range 0 -> 1)
        if (m_session->DefinesFunction(varCoeffs[0]))
        {
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                cout << "Loading Anisotropic Fibre map." << endl;
            }
            for (int i = 0; i < m_spacedim; ++i)
            {
                ASSERTL0(m_session->DefinesFunction(varCoeffs[i], varName),
                    "Function '" + varCoeffs[i] + "' not correctly defined.");
                EvaluateFunction(varName, vTemp, varCoeffs[i]);
                Vmath::Vmul(nq, vTemp, 1, m_vardiff[varCoeffEnum[i]], 1, m_vardiff[varCoeffEnum[i]], 1);
            }
        }

        if (!m_vardiff.empty())
        {
            // Process each of the defined variable coefficients
            for (int i = 0; i < m_spacedim; ++i)
            {
                // If scaling parameters defined, do scaling
                if (m_session->DefinesParameter("d_min"))
                {
                    // Normalise and invert
                    NekDouble f_min = m_session->GetParameter("d_min");
                    NekDouble f_max = m_session->GetParameter("d_max");
                    NekDouble f_range = f_max - f_min;
                    NekDouble o_min = m_session->GetParameter("o_min");
                    NekDouble o_max = m_session->GetParameter("o_max");
                    Vmath::Sadd(nq, -f_min, 
                                    m_vardiff[varCoeffEnum[i]], 1,
                                    m_vardiff[varCoeffEnum[i]], 1);
                    for (int j = 0; j < nq; ++j)
                    {
                        if (m_vardiff[varCoeffEnum[i]][j] < 0)
                        {
                            m_vardiff[varCoeffEnum[i]][j] = 0.0;
                        }
                        if (m_vardiff[varCoeffEnum[i]][j] > f_range)
                        {
                            m_vardiff[varCoeffEnum[i]][j] = f_range;
                        }
                    }
                    Vmath::Smul(nq, -1.0/f_range, 
                                    m_vardiff[varCoeffEnum[i]], 1,
                                    m_vardiff[varCoeffEnum[i]], 1);
                    Vmath::Sadd(nq, 1.0, 
                                    m_vardiff[varCoeffEnum[i]], 1,
                                    m_vardiff[varCoeffEnum[i]], 1);
                    Vmath::Smul(nq, o_max-o_min, 
                                    m_vardiff[varCoeffEnum[i]], 1,
                                    m_vardiff[varCoeffEnum[i]], 1);
                    Vmath::Sadd(nq, o_min, 
                                    m_vardiff[varCoeffEnum[i]], 1,
                                    m_vardiff[varCoeffEnum[i]], 1);
                }

                // Transform variable coefficient and write out to file.
                m_fields[0]->FwdTrans_IterPerExp(m_vardiff[varCoeffEnum[i]],
                                                 m_fields[0]->UpdateCoeffs());
                std::stringstream filename;
                filename << varCoeffs[i];
                if (m_comm->GetSize() > 1)
                {
                    filename << "_P" << m_comm->GetRank();
                }
                filename << ".fld";
                WriteFld(filename.str());
            }
        }

        // Search through the loaded filters and pass the cell model to any
        // CheckpointCellModel filters loaded.
        int k = 0;
        const LibUtilities::FilterMap& f = m_session->GetFilters();
        LibUtilities::FilterMap::const_iterator x;
        for (x = f.begin(); x != f.end(); ++x, ++k)
        {
            if (x->first == "CheckpointCellModel")
            {
                boost::shared_ptr<FilterCheckpointCellModel> c
                    = boost::shared_dynamic_cast<FilterCheckpointCellModel>(
                                                                m_filters[k]);
                c->SetCellModel(m_cell);
            }
        }

        // Load stimuli
        m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);

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
        int nq          = m_fields[0]->GetNpoints();
        StdRegions::ConstFactorMap factors;
        // lambda = \Delta t
        factors[StdRegions::eFactorLambda] = 1.0/lambda*m_chi*m_capMembrane;

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep
            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
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


    /**
     *
     */
    void Monodomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int nq = m_fields[0]->GetNpoints();

        // Compute I_ion
        m_cell->TimeIntegrate(inarray, outarray, time);

        // Compute I_stim
        for (unsigned int i = 0; i < m_stimulus.size(); ++i)
        {   
            m_stimulus[i]->Update(outarray, time);
        }
    }


    /**
     *
     */
    void Monodomain::v_SetInitialConditions(NekDouble initialtime,
                        bool dumpInitialConditions)
    {
        EquationSystem::v_SetInitialConditions(initialtime,
                                               dumpInitialConditions);
        m_cell->Initialise();
    }


    /**
     *
     */
    void Monodomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        if (m_session->DefinesFunction("d00") &&
            m_session->GetFunctionType("d00", "intensity") 
                    == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d00", "intensity")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("d11") &&
            m_session->GetFunctionType("d11", "intensity") 
                    == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d11", "intensity")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("d22") &&
            m_session->GetFunctionType("d22", "intensity") 
                    == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d22", "intensity")->GetExpression()
                << endl;
        }
        m_cell->PrintSummary(out);
    }
}
