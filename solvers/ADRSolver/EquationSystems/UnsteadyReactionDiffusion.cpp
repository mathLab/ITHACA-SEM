///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyReactionDiffusion.cpp
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
// Description: Unsteady reaction-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <boost/core/ignore_unused.hpp>

#include <ADRSolver/EquationSystems/UnsteadyReactionDiffusion.h>

using namespace std;

namespace Nektar
{
string UnsteadyReactionDiffusion::className = GetEquationSystemFactory().
    RegisterCreatorFunction("UnsteadyReactionDiffusion",
                            UnsteadyReactionDiffusion::create);

UnsteadyReactionDiffusion::UnsteadyReactionDiffusion(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

/**
 * @brief Initialisation object for the unsteady reaction-diffusion problem.
 */
void UnsteadyReactionDiffusion::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    ASSERTL0(m_intScheme->GetIntegrationSchemeType() == LibUtilities::eIMEX,
             "Reaction-diffusion requires an implicit-explicit timestepping"
             " (e.g. IMEXOrder2)");
    ASSERTL0(m_projectionType == MultiRegions::eGalerkin,
             "Reaction-diffusion requires use of continuous Galerkin"
             "projection.");

    // Load diffusion parameter
    m_session->LoadParameter("epsilon", m_epsilon,  0.0);

    // Forcing terms
    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    m_ode.DefineOdeRhs       (&UnsteadyReactionDiffusion::DoOdeRhs,        this);
    m_ode.DefineProjection   (&UnsteadyReactionDiffusion::DoOdeProjection, this);
    m_ode.DefineImplicitSolve(&UnsteadyReactionDiffusion::DoImplicitSolve, this);
}

/**
 * @brief Unsteady diffusion problem destructor.
 */
UnsteadyReactionDiffusion::~UnsteadyReactionDiffusion()
{
}

/**
 * @brief Compute the right-hand side for the unsteady reaction diffusion
 * problem.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyReactionDiffusion::DoOdeRhs(
    const Array<OneD, const  Array<OneD, NekDouble> > &inarray,
          Array<OneD,        Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
{
    // RHS should be set to zero.
    for (int i = 0; i < outarray.size(); ++i)
    {
        Vmath::Zero(outarray[i].size(), &outarray[i][0], 1);
    }

    // Add forcing terms for reaction.
    for (auto &x : m_forcing)
    {
        // set up non-linear terms
        x->Apply(m_fields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the projection for the unsteady diffusion problem.
 *
 * @param inarray    Given fields.
 * @param outarray   Calculated solution.
 * @param time       Time.
 */
void UnsteadyReactionDiffusion::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
          Array<OneD,       Array<OneD, NekDouble> > &outarray,
    const NekDouble time)
{
    int i;
    int nvariables = inarray.size();
    SetBoundaryConditions(time);

    Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

    for(i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(inarray[i], coeffs);
        m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
    }
}

/**
 * @brief Implicit solution of the unsteady diffusion problem.
 */
void UnsteadyReactionDiffusion::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
          Array<OneD,       Array<OneD, NekDouble> > &outarray,
    const NekDouble time,
    const NekDouble lambda)
{
    boost::ignore_unused(time);

    StdRegions::ConstFactorMap factors;

    int nvariables = inarray.size();
    int npoints    = m_fields[0]->GetNpoints();
    factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i] inarray = input:
    // \hat{rhs} -> output: \hat{Y} outarray = output: nabla^2 \hat{Y} where
    // \hat = modal coeffs
    for (int i = 0; i < nvariables; ++i)
    {
        // Multiply 1.0/timestep/lambda
        Vmath::Smul(npoints, -factors[StdRegions::eFactorLambda],
                    inarray[i], 1, outarray[i], 1);

        // Solve a system of equations with Helmholtz solver
        m_fields[i]->HelmSolve(
            outarray[i], m_fields[i]->UpdateCoeffs(), factors);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);
        m_fields[i]->SetPhysState(false);
    }
}

}
