///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePropagation.cpp
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
// Description: Pulse Wave Propagation solve routines based on the weak
// formulation (1):
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <PulseWaveSolver/EquationSystems/PulseWavePropagation.h>

using namespace std;

namespace Nektar
{
string PulseWavePropagation::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "PulseWavePropagation", PulseWavePropagation::create,
        "Pulse Wave Propagation equation.");
/**
 *  @class PulseWavePropagation
 *
 *  Set up the routines based on the weak formulation from
 *  "Computational Modelling of 1D blood flow with variable
 *  mechanical properties" by S. J. Sherwin et al. The weak
 *  formulation (1) reads:
 *  \f$ \sum_{e=1}^{N_{el}} \left[ \left( \frac{\partial \mathbf{U}^{\delta}
 * }{\partial t} , \mathbf{\psi}^{\delta} \right)_{\Omega_e} - \left(
 * \frac{\partial \mathbf{F(\mathbf{U})}^{\delta} }
 *    {\partial x}, \mathbf{\psi}^{\delta}  \right)_{\Omega_e} + \left[
 * \mathbf{\psi}^{\delta} \cdot \{ \mathbf{F}^u -
 * \mathbf{F}(\mathbf{U}^{\delta}) \} \right]_{x_e^l}^{x_eû} \right] = 0 \f$
 */
PulseWavePropagation::PulseWavePropagation(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : PulseWaveSystem(pSession, pGraph)
{
}

void PulseWavePropagation::v_InitObject()
{
    PulseWaveSystem::v_InitObject();

    m_pressureArea = GetPressureAreaFactory().CreateInstance(
        "Lymphatic", m_vessels, m_session);
    m_pressureArea->DoPressure();

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&PulseWavePropagation::DoOdeRhs, this);
        m_ode.DefineProjection(&PulseWavePropagation::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit Pulse Wave Propagation not set up.");
    }

    // Create advection object
    string advName;
    string riemName;
    switch (m_upwindTypePulse)
    {
        case eUpwindPulse:
        {
            advName  = "WeakDG";
            riemName = "UpwindPulse";
        }
        break;
        default:
        {
            ASSERTL0(false, "populate switch statement for upwind flux");
        }
        break;
    }
    m_advObject =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
    m_advObject->SetFluxVector(&PulseWavePropagation::GetFluxVector, this);
    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
        riemName, m_session);
    m_riemannSolver->SetScalar("A0", &PulseWavePropagation::GetA0, this);
    m_riemannSolver->SetScalar("beta", &PulseWavePropagation::GetBeta, this);
    m_riemannSolver->SetScalar("N", &PulseWavePropagation::GetN, this);
    m_riemannSolver->SetParam("rho", &PulseWavePropagation::GetRho, this);
    m_riemannSolver->SetParam("pext", &PulseWavePropagation::GetPext, this);

    m_advObject->SetRiemannSolver(m_riemannSolver);
    m_advObject->InitObject(m_session, m_fields);
}

PulseWavePropagation::~PulseWavePropagation()
{
}

/**
 *  Computes the right hand side of (1). The RHS is everything
 *  except the term that contains the time derivative
 *  \f$\frac{\partial \mathbf{U}}{\partial t}\f$. In case of a
 *  Discontinuous Galerkin projection, m_advObject->Advect
 *  will be called
 *
 */
void PulseWavePropagation::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int i;

    Array<OneD, Array<OneD, NekDouble>> physarray(m_nVariables);

    // Dummy array for WeakDG advection
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    // Output array for advection
    Array<OneD, Array<OneD, NekDouble>> out(m_nVariables);

    int cnt = 0;

    // Set up Inflow and Outflow boundary conditions.
    SetPulseWaveBoundaryConditions(inarray, outarray, time);

    // Set up any interface conditions and write into boundary condition
    EnforceInterfaceConditions(inarray);

    // do advection evauation in all domains
    for (int omega = 0; omega < m_nDomains; ++omega)
    {
        m_currentDomain = omega;
        int nq          = m_vessels[omega * m_nVariables]->GetTotPoints();

        for (i = 0; i < m_nVariables; ++i)
        {
            physarray[i] = inarray[i] + cnt;
            out[i]       = outarray[i] + cnt;
        }

        for (i = 0; i < m_nVariables; ++i)
        {
            m_fields[i] = m_vessels[omega * m_nVariables + i];
        }

        m_advObject->Advect(m_nVariables, m_fields, advVel, physarray, out,
                            time);
        for (i = 0; i < m_nVariables; ++i)
        {
            Vmath::Neg(nq, out[i], 1);
        }

        cnt += nq;
    }
}

void PulseWavePropagation::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Just copy over array
    for (int i = 0; i < m_nVariables; ++i)
    {
        Vmath::Vcopy(inarray[i].size(), inarray[i], 1, outarray[i], 1);
    }
}

/**
 *	Does the projection between ... space and the ... space. Also checks for
 *Q-inflow boundary conditions at the inflow of the current arterial segment and
 *applies the Q-inflow if specified
 */
void PulseWavePropagation::SetPulseWaveBoundaryConditions(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)

{
    int omega;

    Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

    int offset = 0;

    // This will be moved to the RCR boundary condition once factory is setup
    if (time == 0)
    {
        m_Boundary = Array<OneD, PulseWaveBoundarySharedPtr>(2 * m_nDomains);

        for (omega = 0; omega < m_nDomains; ++omega)
        {
            vessel[0] = m_vessels[2 * omega];
            vessel[1] = m_vessels[2 * omega + 1];

            for (int j = 0; j < 2; ++j)
            {
                std::string BCType =
                    vessel[0]->GetBndConditions()[j]->GetUserDefined();
                if (BCType.empty()) // if not condition given define it to be
                                    // NoUserDefined
                {
                    BCType = "NoUserDefined";
                }

                m_Boundary[2 * omega + j] = GetBoundaryFactory().CreateInstance(
                    BCType, m_vessels, m_session, m_pressureArea);

                // turn on time depedent BCs
                if (BCType == "Q-inflow")
                {
                    vessel[0]->GetBndConditions()[j]->SetIsTimeDependent(true);
                }
                if (BCType == "A-inflow")
                {
                    vessel[0]->GetBndConditions()[j]->SetIsTimeDependent(true);
                }
                if (BCType == "U-inflow")
                {
                    vessel[1]->GetBndConditions()[j]->SetIsTimeDependent(true);
                }
                else if (BCType == "RCR-terminal")
                {
                    vessel[0]->GetBndConditions()[j]->SetIsTimeDependent(true);
                }
            }
        }
    }

    SetBoundaryConditions(time);

    // Loop over all vessesls and set boundary conditions
    for (omega = 0; omega < m_nDomains; ++omega)
    {
        for (int n = 0; n < 2; ++n)
        {
            m_Boundary[2 * omega + n]->DoBoundary(inarray, m_A_0, m_beta, time,
                                                  omega, offset, n);
        }
        offset += m_vessels[2 * omega]->GetTotPoints();
    }
}

/**
 *  Calculates the second term of the weak form (1): \f$
 *  \left( \frac{\partial \mathbf{F(\mathbf{U})}^{\delta}
 *  }{\partial x}, \mathbf{\psi}^{\delta} \right)_{\Omega_e}
 *  \f$
 *  The variables of the system are $\mathbf{U} = [A,u]^T$
 *  physfield[0] = A        physfield[1] = u
 *  flux[0] = F[0] = A*u    flux[1] = F[1] = u^2/2 + p/rho
 *  p-A-relationship: p = p_ext + beta*(sqrt(A)-sqrt(A_0))
 */
void PulseWavePropagation::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int nq        = m_vessels[m_currentDomain * m_nVariables]->GetTotPoints();
    NekDouble p   = 0.0;
    NekDouble p_t = 0.0;

    for (int j = 0; j < nq; j++)
    {
        flux[0][0][j] = physfield[0][j] * physfield[1][j];

        ASSERTL0(physfield[0][j] >= 0, "Negative A not allowed.");

        p = m_pext +
            m_beta[m_currentDomain][j] *
                (sqrt(physfield[0][j]) - sqrt(m_A_0[m_currentDomain][j]));

        p_t           = (physfield[1][j] * physfield[1][j]) / 2 + p / m_rho;
        flux[1][0][j] = p_t;
    }
}

Array<OneD, NekDouble> &PulseWavePropagation::GetA0()
{
    return m_A_0_trace[m_currentDomain];
}

Array<OneD, NekDouble> &PulseWavePropagation::GetBeta()
{
    return m_beta_trace[m_currentDomain];
}

Array<OneD, NekDouble> &PulseWavePropagation::GetN()
{
    return m_trace_fwd_normal[m_currentDomain];
}

NekDouble PulseWavePropagation::GetRho()
{
    return m_rho;
}

NekDouble PulseWavePropagation::GetPext()
{
    return m_pext;
}

/**
 *  Print summary routine, calls virtual routine reimplemented in
 *  UnsteadySystem
 */
void PulseWavePropagation::v_GenerateSummary(SolverUtils::SummaryList &s)
{
    PulseWaveSystem::v_GenerateSummary(s);
}
} // namespace Nektar
