///////////////////////////////////////////////////////////////////////////////
//
// File BidomainRoth.cpp
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
// Description: Bidomain cardiac electrophysiology model - Roth formulation.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <CardiacEPSolver/EquationSystems/BidomainRoth.h>
#include <CardiacEPSolver/Filters/FilterCheckpointCellModel.h>

using namespace std;

namespace Nektar
{

/**
 * Registers the class with the Factory.
 */
string BidomainRoth::className
        = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "BidomainRoth",
            BidomainRoth::create,
            "Bidomain Roth model of cardiac electrophysiology.");


/**
 *
 */
BidomainRoth::BidomainRoth(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}


/**
 *
 */
void BidomainRoth::v_InitObject()
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

    // Allocate storage for variable coeffs and initialize to 1.
    for (int i = 0, k = 0; i < m_spacedim; ++i)
    {
        for (int j = 0; j < i+1; ++j)
        {
            if (i == j)
            {
                m_vardiffi[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 1.0);
                m_vardiffe[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 1.0);
                m_vardiffie[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 1.0);
            }
            else
            {
                m_vardiffi[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 0.0);
                m_vardiffe[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 0.0);
                m_vardiffie[varCoeffEnum[k]] = Array<OneD, NekDouble>(nq, 0.0);
            }
            ++k;
        }
    }

    // Apply fibre map f \in [0,1], scale to conductivity range
    // [o_min,o_max], specified by the session parameters o_min and o_max
    if (m_session->DefinesFunction("ExtracellularAnisotropicConductivity"))
    {
        if (m_session->DefinesCmdLineArgument("verbose"))
        {
            cout << "Loading Extracellular Anisotropic Fibre map." << endl;
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
            ASSERTL0(m_session->DefinesFunction(
                                        "ExtracellularAnisotropicConductivity",
                                        aniso_var[j]),
                     "Function 'AnisotropicConductivity' not correctly "
                     "defined.");

            GetFunction("ExtracellularAnisotropicConductivity")->Evaluate(aniso_var[j], vTemp_j);

            // Loop through rows of D
            for (int i = 0; i < j + 1; ++i)
            {
                ASSERTL0(m_session->DefinesFunction(
                                    "ExtracellularAnisotropicConductivity",
                                    aniso_var[i]),
                         "Function 'ExtracellularAnisotropicConductivity' not "
                         "correctly defined.");

                GetFunction("ExtracellularAnisotropicConductivity")->Evaluate(aniso_var[i], vTemp_i);

                Vmath::Vmul(nq, vTemp_i, 1, vTemp_j, 1,
                                m_vardiffe[varCoeffEnum[k]], 1);

                Vmath::Smul(nq, o_max-o_min,
                                m_vardiffe[varCoeffEnum[k]], 1,
                                m_vardiffe[varCoeffEnum[k]], 1);

                if (i == j)
                {
                    Vmath::Sadd(nq, o_min,
                                    m_vardiffe[varCoeffEnum[k]], 1,
                                    m_vardiffe[varCoeffEnum[k]], 1);
                }

            }
        }
    }

    // Apply fibre map f \in [0,1], scale to conductivity range
    // [o_min,o_max], specified by the session parameters o_min and o_max
    if (m_session->DefinesFunction("IntracellularAnisotropicConductivity"))
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
            ASSERTL0(m_session->DefinesFunction(
                                        "IntracellularAnisotropicConductivity",
                                        aniso_var[j]),
                     "Function 'IntracellularAnisotropicConductivity' not "
                     "correctly defined.");

            GetFunction("IntracellularAnisotropicConductivity")->Evaluate(aniso_var[j], vTemp_j);

            // Loop through rows of D
            for (int i = 0; i < j + 1; ++i)
            {
                ASSERTL0(m_session->DefinesFunction(
                                    "IntracellularAnisotropicConductivity",
                                    aniso_var[i]),
                         "Function 'IntracellularAnisotropicConductivity' not "
                         "correctly defined.");
                GetFunction("IntracellularAnisotropicConductivity")->Evaluate(aniso_var[i], vTemp_i);

                Vmath::Vmul(nq, vTemp_i, 1, vTemp_j, 1,
                                m_vardiffi[varCoeffEnum[k]], 1);

                Vmath::Smul(nq, o_max-o_min,
                                m_vardiffi[varCoeffEnum[k]], 1,
                                m_vardiffi[varCoeffEnum[k]], 1);

                if (i == j)
                {
                    Vmath::Sadd(nq, o_min,
                                    m_vardiffi[varCoeffEnum[k]], 1,
                                    m_vardiffi[varCoeffEnum[k]], 1);
                }

                Vmath::Vadd(nq, m_vardiffe[varCoeffEnum[k]], 1,
                                m_vardiffi[varCoeffEnum[k]], 1,
                                m_vardiffie[varCoeffEnum[k]], 1);

                ++k;
            }
        }
    }


    // Write out conductivity values
    for (int j = 0, k = 0; j < m_spacedim; ++j)
    {
        // Loop through rows of D
        for (int i = 0; i < j + 1; ++i)
        {
            // Transform variable coefficient and write out to file.
            m_fields[0]->FwdTrans_IterPerExp(m_vardiffi[varCoeffEnum[k]],
                                             m_fields[0]->UpdateCoeffs());
            std::stringstream filenamei;
            filenamei << "IConductivity_" << varCoeffString[k] << ".fld";
            WriteFld(filenamei.str());

            // Transform variable coefficient and write out to file.
            m_fields[0]->FwdTrans_IterPerExp(m_vardiffe[varCoeffEnum[k]],
                                             m_fields[0]->UpdateCoeffs());
            std::stringstream filenamee;
            filenamee << "EConductivity_" << varCoeffString[k] << ".fld";
            WriteFld(filenamee.str());

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
    }
    // Load stimuli
    m_stimulus = Stimulus::LoadStimuli(m_session, m_fields[0]);

    if (!m_explicitDiffusion)
    {
        m_ode.DefineImplicitSolve (&BidomainRoth::DoImplicitSolve, this);
    }
    m_ode.DefineOdeRhs(&BidomainRoth::DoOdeRhs, this);
}


/**
 *
 */
BidomainRoth::~BidomainRoth()
{

}


/**
 * @param   inarray         Input array.
 * @param   outarray        Output array.
 * @param   time            Current simulation time.
 * @param   lambda          Timestep.
 */
void BidomainRoth::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
              Array<OneD, Array<OneD, NekDouble> >&outarray,
        const NekDouble time,
        const NekDouble lambda)
{
    int nq          = m_fields[0]->GetNpoints();

    StdRegions::ConstFactorMap factorsHelmholtz;
    // lambda = \Delta t
    factorsHelmholtz[StdRegions::eFactorLambda]
                    = 1.0/lambda*m_chi*m_capMembrane;

    // ------------------------------
    // Solve Helmholtz problem for Vm
    // ------------------------------
    // Multiply 1.0/timestep
    //Vmath::Vadd(nq, inarray[0], 1, ggrad, 1, m_fields[0]->UpdatePhys(), 1);
    Vmath::Smul(nq, -factorsHelmholtz[StdRegions::eFactorLambda], inarray[0], 1,
                                    m_fields[0]->UpdatePhys(), 1);

    // Solve a system of equations with Helmholtz solver and transform
    // back into physical space.
    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),
                           m_fields[0]->UpdateCoeffs(), 
                           factorsHelmholtz, m_vardiffe);

    m_fields[0]->BwdTrans( m_fields[0]->GetCoeffs(),
                           m_fields[0]->UpdatePhys());
    m_fields[0]->SetPhysState(true);

    // Copy the solution vector (required as m_fields must be set).
    outarray[0] = m_fields[0]->GetPhys();
}


/**
 *
 */
void BidomainRoth::DoOdeRhs(
        const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
              Array<OneD,        Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
{
    int nq          = m_fields[0]->GetNpoints();

    // Compute I_ion
    m_cell->TimeIntegrate(inarray, outarray, time);

    // Compute I_stim
    for (unsigned int i = 0; i < m_stimulus.size(); ++i)
    {
        m_stimulus[i]->Update(outarray, time);
    }

    Array<OneD, NekDouble> ggrad0(nq), ggrad1(nq), ggrad2(nq), ggrad(nq);
    StdRegions::ConstFactorMap factorsPoisson;
    factorsPoisson[StdRegions::eFactorLambda] = 0.0;

    // ----------------------------
    // Compute \nabla g_i \nabla Vm
    // ----------------------------
    m_fields[0]->PhysDeriv(inarray[0], ggrad0, ggrad1, ggrad2);
    m_fields[0]->PhysDeriv(0, ggrad0, ggrad0);
    m_fields[0]->PhysDeriv(1, ggrad1, ggrad1);
    m_fields[0]->PhysDeriv(2, ggrad2, ggrad2);
    if (m_session->DefinesFunction("IntracellularAnisotropicConductivity") &&
        m_session->DefinesFunction("ExtracellularAnisotropicConductivity"))
    {
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD00], 1, ggrad0,
                1, ggrad0, 1);
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD11], 1, ggrad1,
                1, ggrad1, 1);
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD22], 1, ggrad2,
                1, ggrad2, 1);
    }
    // Add partial derivatives together
    Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
    Vmath::Vadd(nq, ggrad2, 1, ggrad, 1, ggrad, 1);

    Vmath::Smul(nq, -1.0, ggrad, 1, m_fields[1]->UpdatePhys(), 1);

    // ----------------------------
    // Solve Poisson problem for Ve
    // ----------------------------
    m_fields[1]->HelmSolve(m_fields[1]->GetPhys(),
            m_fields[1]->UpdateCoeffs(), factorsPoisson, m_vardiffie);
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(),
            m_fields[1]->UpdatePhys());
    m_fields[1]->SetPhysState(true);

    // ------------------------------
    // Compute Laplacian of Ve (forcing term)
    // ------------------------------
    m_fields[1]->PhysDeriv(m_fields[1]->GetPhys(), ggrad0, ggrad1, ggrad2);
    m_fields[1]->PhysDeriv(0, ggrad0, ggrad0);
    m_fields[1]->PhysDeriv(1, ggrad1, ggrad1);
    m_fields[1]->PhysDeriv(2, ggrad2, ggrad2);
    if (m_session->DefinesFunction("IntracellularAnisotropicConductivity") &&
        m_session->DefinesFunction("ExtracellularAnisotropicConductivity"))
    {
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD00], 1, ggrad0,
                1, ggrad0, 1);
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD11], 1, ggrad1,
                1, ggrad1, 1);
        Vmath::Vmul(nq, m_vardiffi[StdRegions::eVarCoeffD22], 1, ggrad2,
                1, ggrad2, 1);
    }
    // Add partial derivatives together
    Vmath::Vadd(nq, ggrad0, 1, ggrad1, 1, ggrad, 1);
    Vmath::Vadd(nq, ggrad2, 1, ggrad, 1, ggrad, 1);

    Vmath::Vadd(nq, ggrad, 1, outarray[0], 1, outarray[0], 1);
}


/**
 *
 */
void BidomainRoth::v_SetInitialConditions(NekDouble initialtime,
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
void BidomainRoth::v_GenerateSummary(SummaryList& s)
{
    UnsteadySystem::v_GenerateSummary(s);
    m_cell->GenerateSummary(s);
}

}
