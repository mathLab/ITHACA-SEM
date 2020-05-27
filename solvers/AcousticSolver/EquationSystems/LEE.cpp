///////////////////////////////////////////////////////////////////////////////
//
// File LEE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Linearized Euler Equations
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <AcousticSolver/EquationSystems/LEE.h>

using namespace std;

namespace Nektar
{
string LEE::className = GetEquationSystemFactory().RegisterCreatorFunction(
    "LEE", LEE::create, "Linearized Euler Equations");

LEE::LEE(const LibUtilities::SessionReaderSharedPtr &pSession,
         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AcousticSystem(pSession, pGraph)
{
    m_ip   = 0;
    m_irho = 1;
    m_iu   = 2;

    m_conservative = true;
}

/**
 * @brief Initialization object for the LEE class.
 */
void LEE::v_InitObject()
{
    AcousticSystem::v_InitObject();

    m_bfNames.push_back("gamma");

    // Initialize basefield again
    m_bf = Array<OneD, Array<OneD, NekDouble>>(m_bfNames.size());
    for (int i = 0; i < m_bf.size(); ++i)
    {
        m_bf[i] = Array<OneD, NekDouble>(GetTotPoints());
    }
    GetFunction("Baseflow", m_fields[0], true)
        ->Evaluate(m_bfNames, m_bf, m_time);

    // Define the normal velocity fields
    m_bfFwdBwd = Array<OneD, Array<OneD, NekDouble>>(2 * m_bfNames.size());
    for (int i = 0; i < m_bfFwdBwd.size(); i++)
    {
        m_bfFwdBwd[i] = Array<OneD, NekDouble>(GetTraceNpoints(), 0.0);
    }

    string riemName;
    m_session->LoadSolverInfo("UpwindType", riemName, "Upwind");
    if (boost::to_lower_copy(riemName) == "characteristics" ||
        boost::to_lower_copy(riemName) == "leeupwind" ||
        boost::to_lower_copy(riemName) == "upwind")
    {
        riemName = "LEEUpwind";
    }
    if (boost::to_lower_copy(riemName) == "laxfriedrichs")
    {
        riemName = "LEELaxFriedrichs";
    }
    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
        riemName, m_session);
    m_riemannSolver->SetVector("N", &LEE::GetNormals, this);
    m_riemannSolver->SetVector("basefieldFwdBwd", &LEE::GetBasefieldFwdBwd,
                               this);
    m_riemannSolver->SetAuxVec("vecLocs", &LEE::GetVecLocs, this);

    // Set up advection operator
    string advName;
    m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    m_advection =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
    m_advection->SetFluxVector(&LEE::v_GetFluxVector, this);
    m_advection->SetRiemannSolver(m_riemannSolver);
    m_advection->InitObject(m_session, m_fields);

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&LEE::DoOdeRhs, this);
        m_ode.DefineProjection(&LEE::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit LEE not set up.");
    }
}

/**
 * @brief Destructor for LEE class.
 */
LEE::~LEE()
{
}

/**
 * @brief Return the flux vector for the LEE equations.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux. flux[eq][dir][pt]
 */
void LEE::v_GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int nq = physfield[0].size();

    ASSERTL1(flux[0].size() == m_spacedim,
             "Dimension of flux array and velocity array do not match");

    Array<OneD, const NekDouble> c0sq = m_bf[0];
    Array<OneD, Array<OneD, const NekDouble>> u0(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        u0[i] = m_bf[2 + i];
    }

    Array<OneD, const NekDouble> p   = physfield[m_ip];
    Array<OneD, const NekDouble> rho = physfield[m_irho];
    Array<OneD, Array<OneD, const NekDouble>> ru(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        ru[i] = physfield[m_iu + i];
    }

    // F_{adv,p',j} = c0^2 * ru_j + u0_j * p
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 0;
        Vmath::Vvtvvtp(nq, c0sq, 1, ru[j], 1, u0[j], 1, p, 1, flux[i][j], 1);
    }

    // F_{adv,rho',j} = u0_j * rho' + ru_j
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 1;
        // u0_j * rho' + ru_j
        Vmath::Vvtvp(nq, u0[j], 1, rho, 1, ru[j], 1, flux[i][j], 1);
    }

    for (int i = 0; i < m_spacedim; ++i)
    {
        // F_{adv,u'_i,j} = ru_i * u0_j + delta_ij * p
        for (int j = 0; j < m_spacedim; ++j)
        {
            // ru_i * u0_j
            Vmath::Vmul(nq, ru[i], 1, u0[j], 1, flux[m_iu + i][j], 1);

            // kronecker delta
            if (i == j)
            {
                // delta_ij + p
                Vmath::Vadd(nq, p, 1, flux[m_iu + i][j], 1, flux[m_iu + i][j],
                            1);
            }
        }
    }
}

void LEE::v_AddLinTerm(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                       Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int nq = GetTotPoints();

    Array<OneD, Array<OneD, NekDouble>> linTerm(m_spacedim + 2);
    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        if (i == 1)
        {
            // skip rho
            continue;
        }

        linTerm[i] = Array<OneD, NekDouble>(nq);
    }

    Array<OneD, const NekDouble> c0sq  = m_bf[0];
    Array<OneD, const NekDouble> rho0  = m_bf[1];
    Array<OneD, const NekDouble> gamma = m_bf[m_iu + m_spacedim];

    Array<OneD, NekDouble> gammaMinOne(nq);
    Vmath::Sadd(nq, -1.0, gamma, 1, gammaMinOne, 1);

    Array<OneD, NekDouble> p0(nq);
    Vmath::Vmul(nq, c0sq, 1, rho0, 1, p0, 1);
    Vmath::Vdiv(nq, p0, 1, gamma, 1, p0, 1);

    Array<OneD, Array<OneD, const NekDouble>> u0(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        u0[i] = m_bf[2 + i];
    }

    Array<OneD, const NekDouble> p   = inarray[0];
    Array<OneD, const NekDouble> rho = inarray[1];

    Array<OneD, Array<OneD, const NekDouble>> ru(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        ru[i] = inarray[2 + i];
    }

    Array<OneD, NekDouble> grad(nq);
    Array<OneD, NekDouble> tmp1(nq);

    // p
    {
        Vmath::Zero(nq, linTerm[m_ip], 1);
        // (1-gamma) ( ru_j / rho0 * dp0/dx_j - p * du0_j/dx_j )
        for (int j = 0; j < m_spacedim; ++j)
        {
            // ru_j / rho0 * dp0/dx_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], p0, grad);
            Vmath::Vmul(nq, grad, 1, ru[j], 1, tmp1, 1);
            Vmath::Vdiv(nq, tmp1, 1, rho0, 1, tmp1, 1);
            // p * du0_j/dx_j - ru_j / rho0 * dp0/dx_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], u0[j],
                                   grad);
            Vmath::Vvtvm(nq, grad, 1, p, 1, tmp1, 1, tmp1, 1);
            // (gamma-1) (p * du0_j/dx_j - ru_j / rho0 * dp0/dx_j)
            Vmath::Vvtvp(nq, gammaMinOne, 1, tmp1, 1, linTerm[m_ip], 1,
                         linTerm[m_ip], 1);
        }
    }

    // rho has no linTerm

    // ru_i
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Zero(nq, linTerm[m_iu + i], 1);
        // du0_i/dx_j * (u0_j * rho + ru_j)
        for (int j = 0; j < m_spacedim; ++j)
        {
            // d u0_i / d x_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], u0[i],
                                   grad);
            // u0_j * rho + ru_j
            Vmath::Vvtvp(nq, u0[j], 1, rho, 1, ru[j], 1, tmp1, 1);
            // du0_i/dx_j * (u0_j * rho + ru_j)
            Vmath::Vvtvp(nq, grad, 1, tmp1, 1, linTerm[m_iu + i], 1,
                         linTerm[m_iu + i], 1);
        }
    }

    Array<OneD, NekDouble> tmpC(GetNcoeffs());
    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        if (i == 1)
        {
            // skip rho
            continue;
        }

        m_fields[0]->FwdTrans(linTerm[i], tmpC);
        m_fields[0]->BwdTrans(tmpC, linTerm[i]);

        Vmath::Vadd(nq, outarray[i], 1, linTerm[i], 1, outarray[i], 1);
    }
}

/**
 * @brief Outflow characteristic boundary conditions for compressible
 * flow problems.
 */
void LEE::v_RiemannInvariantBC(int bcRegion, int cnt,
                               Array<OneD, Array<OneD, NekDouble>> &Fwd,
                               Array<OneD, Array<OneD, NekDouble>> &BfFwd,
                               Array<OneD, Array<OneD, NekDouble>> &physarray)
{
    int id1, id2, nBCEdgePts;
    int nVariables = physarray.size();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    for (int e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]
                         ->GetBndCondExpansions()[bcRegion]
                         ->GetExp(e)
                         ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

        // Calculate (v.n)
        Array<OneD, NekDouble> RVn(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &Fwd[m_iu + i][id2], 1,
                         &m_traceNormals[i][id2], 1, &RVn[0], 1, &RVn[0], 1);
        }

        // Calculate (v0.n)
        Array<OneD, NekDouble> RVn0(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &BfFwd[2 + i][id2], 1,
                         &m_traceNormals[i][id2], 1, &RVn0[0], 1, &RVn0[0], 1);
        }

        for (int i = 0; i < nBCEdgePts; ++i)
        {
            NekDouble c = sqrt(BfFwd[0][id2 + i]);

            NekDouble h1, h4, h5;

            if (RVn0[i] > 0)
            {
                // rho - p / c^2
                h1 = Fwd[m_irho][id2 + i] - Fwd[m_ip][id2 + i] / (c * c);
            }
            else
            {
                h1 = 0.0;
            }

            if (RVn0[i] - c > 0)
            {
                // ru / 2 - p / (2*c)
                h4 = RVn[i] / 2 - Fwd[m_ip][id2 + i] / (2 * c);
            }
            else
            {
                h4 = 0.0;
            }

            if (RVn0[i] + c > 0)
            {
                // ru / 2 + p / (2*c)
                h5 = RVn[i] / 2 + Fwd[m_ip][id2 + i] / (2 * c);
            }
            else
            {
                h5 = 0.0;
            }

            // compute conservative variables
            // p = c0*(h5-h4)
            // rho = h1 + (h5-h4)/c0
            // ru = h4+h5
            Fwd[m_ip][id2 + i]   = c * (h5 - h4);
            Fwd[m_irho][id2 + i] = h1 + (h5 - h4) / c;
            NekDouble RVnNew     = h4 + h5;

            // adjust velocity pert. according to new value
            // here we just omit the wall parallel velocity components, i.e
            // setting them to zero. This is equivalent to setting the two
            // vorticity characteristics h2 and h3 to zero. Mathematically,
            // this is only legitimate for incoming characteristics. However,
            // as h2 and h3 are convected by the flow, the value we precribe at
            // an the boundary for putgoing characteristics does not matter.
            // This implementation saves a few operations and is more robust
            // for mixed in/outflow boundaries and at the boundaries edges.
            for (int j = 0; j < m_spacedim; ++j)
            {
                Fwd[m_iu + j][id2 + i] = RVnNew * m_traceNormals[j][id2 + i];
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

} // namespace Nektar
