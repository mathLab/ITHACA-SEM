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
#include <CardiacEPSolver/Filters/FilterCellHistoryPoints.h>

using namespace std;

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
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph)
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
        StdRegions::VarCoeffType varCoeffEnum[6] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD01,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD02,
                StdRegions::eVarCoeffD12,
                StdRegions::eVarCoeffD22
        };
        std::string varCoeffString[6] = {"xx","xy","yy","xz","yz","zz"};
        std::string aniso_var[3] = {"fx", "fy", "fz"};

        const int nq            = m_fields[0]->GetNpoints();
        const int nVarDiffCmpts = m_spacedim * (m_spacedim + 1) / 2;

        // Allocate storage for variable coeffs and initialize to 1.
        for (int i = 0, k = 0; i < m_spacedim; ++i)
        {
            for (int j = 0; j < i+1; ++j)
            {
                if (i == j)
                {
                    m_vardiff[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 1.0);
                }
                else
                {
                    m_vardiff[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 0.0);
                }
                ++k;
            }
        }

        // Apply fibre map f \in [0,1], scale to conductivity range
        // [o_min,o_max], specified by the session parameters o_min and o_max
        if (m_session->DefinesFunction("AnisotropicConductivity"))
        {
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                cout << "Loading Anisotropic Fibre map." << endl;
            }

            NekDouble   o_min        = m_session->GetParameter("o_min");
            NekDouble   o_max        = m_session->GetParameter("o_max");
            int         k            = 0;

            Array<OneD, NekDouble> vTemp_i;
            Array<OneD, NekDouble> vTemp_j;

            /*
             * Diffusivity matrix D is upper triangular and defined as
             *    d_00   d_01  d_02
             *           d_11  d_12
             *                 d_22
             *
             * Given a principle fibre direction _f_ the diffusivity is given
             * by
             *    d_ij = { D_2 + (D_1 - D_2) f_i f_j   if i==j
             *           {       (D_1 - D_2) f_i f_j   if i!=j
             *
             * The vector _f_ is given in terms of the variables fx,fy,fz in the
             * function AnisotropicConductivity. The values of D_1 and D_2 are
             * the parameters o_max and o_min, respectively.
             */

            // Loop through columns of D
            for (int j = 0; j < m_spacedim; ++j)
            {
                ASSERTL0(m_session->DefinesFunction("AnisotropicConductivity",
                                                    aniso_var[j]),
                         "Function 'AnisotropicConductivity' not correctly "
                         "defined.");
                GetFunction("AnisotropicConductivity")->Evaluate(aniso_var[j], vTemp_j);

                // Loop through rows of D
                for (int i = 0; i < j + 1; ++i)
                {
                    ASSERTL0(m_session->DefinesFunction(
                                        "AnisotropicConductivity",aniso_var[i]),
                             "Function 'AnisotropicConductivity' not correctly "
                             "defined.");
                    GetFunction("AnisotropicConductivity")->Evaluate(aniso_var[i], vTemp_i);

                    Vmath::Vmul(nq, vTemp_i, 1, vTemp_j, 1,
                                    m_vardiff[varCoeffEnum[k]], 1);

                    Vmath::Smul(nq, o_max-o_min,
                                    m_vardiff[varCoeffEnum[k]], 1,
                                    m_vardiff[varCoeffEnum[k]], 1);

                    if (i == j)
                    {
                        Vmath::Sadd(nq, o_min,
                                        m_vardiff[varCoeffEnum[k]], 1,
                                        m_vardiff[varCoeffEnum[k]], 1);
                    }

                    ++k;
                }
            }
        }
        else
        {
            // Otherwise apply isotropic conductivity value (o_max) to
            // diagonal components of tensor
            NekDouble o_max = m_session->GetParameter("o_max");
            for (int i = 0; i < nVarDiffCmpts; ++i)
            {
                Vmath::Smul(nq,o_max,
                                m_vardiff[varCoeffEnum[i]], 1,
                                m_vardiff[varCoeffEnum[i]], 1);
            }
        }

        // Scale by scar map (range 0->1) derived from intensity map
        // (range d_min -> d_max)
        if (m_session->DefinesFunction("IsotropicConductivity"))
        {
            if (m_session->DefinesCmdLineArgument("verbose"))
            {
                cout << "Loading Isotropic Conductivity map." << endl;
            }

            const std::string varName  = "intensity";
            Array<OneD, NekDouble> vTemp;
            GetFunction( "IsotropicConductivity")->Evaluate(varName,  vTemp);

            // If the d_min and d_max parameters are defined, then we need to
            // rescale the isotropic conductivity to convert from the source
            // domain (e.g. late-gad intensity) to conductivity
            if ( m_session->DefinesParameter("d_min") ||
                 m_session->DefinesParameter("d_max") ) {
                const NekDouble   f_min    = m_session->GetParameter("d_min");
                const NekDouble   f_max    = m_session->GetParameter("d_max");
                const NekDouble   scar_min = 0.0;
                const NekDouble   scar_max = 1.0;

                // Threshold based on d_min, d_max
                for (int j = 0; j < nq; ++j)
                {
                    vTemp[j] = (vTemp[j] < f_min ? f_min : vTemp[j]);
                    vTemp[j] = (vTemp[j] > f_max ? f_max : vTemp[j]);
                }

                // Rescale to s \in [0,1] (0 maps to d_max, 1 maps to d_min)
                Vmath::Sadd(nq, -f_min, vTemp, 1, vTemp, 1);
                Vmath::Smul(nq, -1.0/(f_max-f_min), vTemp, 1, vTemp, 1);
                Vmath::Sadd(nq, 1.0, vTemp, 1, vTemp, 1);
                Vmath::Smul(nq, scar_max - scar_min, vTemp, 1, vTemp, 1);
                Vmath::Sadd(nq, scar_min, vTemp, 1, vTemp, 1);
            }

            // Scale anisotropic conductivity values
            for (int i = 0; i < nVarDiffCmpts; ++i)
            {
                Vmath::Vmul(nq, vTemp, 1,
                                m_vardiff[varCoeffEnum[i]], 1,
                                m_vardiff[varCoeffEnum[i]], 1);
            }
        }


        // Write out conductivity values
        for (int j = 0, k = 0; j < m_spacedim; ++j)
        {
            // Loop through rows of D
            for (int i = 0; i < j + 1; ++i)
            {
                // Transform variable coefficient and write out to file.
                m_fields[0]->FwdTrans_IterPerExp(m_vardiff[varCoeffEnum[k]],
                                                 m_fields[0]->UpdateCoeffs());
                std::stringstream filename;
                filename << "Conductivity_" << varCoeffString[k] << ".fld";
                WriteFld(filename.str());

                ++k;
            }
        }

        // Search through the loaded filters and pass the cell model to any
        // CheckpointCellModel filters loaded.
        for (auto &x : m_filters)
        {
            if (x.first == "CheckpointCellModel")
            {
                std::shared_ptr<FilterCheckpointCellModel> c
                    = std::dynamic_pointer_cast<FilterCheckpointCellModel>(
                        x.second);
                c->SetCellModel(m_cell);
            }
            if (x.first == "CellHistoryPoints")
            {
                std::shared_ptr<FilterCellHistoryPoints> c
                    = std::dynamic_pointer_cast<FilterCellHistoryPoints>(
                        x.second);
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
        int nvariables  = inarray.size();
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
                                   m_fields[i]->UpdateCoeffs(), 
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
                        bool dumpInitialConditions,
                        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime,
                                               dumpInitialConditions,
                                               domain);
        m_cell->Initialise();
    }


    /**
     *
     */
    void Monodomain::v_GenerateSummary(SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
        if (m_session->DefinesFunction("d00") &&
            m_session->GetFunctionType("d00", "intensity")
                    == LibUtilities::eFunctionTypeExpression)
        {
            AddSummaryItem(s, "Diffusivity-x",
                m_session->GetFunction("d00", "intensity")->GetExpression());
        }
        if (m_session->DefinesFunction("d11") &&
            m_session->GetFunctionType("d11", "intensity")
                    == LibUtilities::eFunctionTypeExpression)
        {
            AddSummaryItem(s, "Diffusivity-y",
                m_session->GetFunction("d11", "intensity")->GetExpression());
        }
        if (m_session->DefinesFunction("d22") &&
            m_session->GetFunctionType("d22", "intensity")
                    == LibUtilities::eFunctionTypeExpression)
        {
            AddSummaryItem(s, "Diffusivity-z",
                m_session->GetFunction("d22", "intensity")->GetExpression());
        }
        m_cell->GenerateSummary(s);
    }
}
