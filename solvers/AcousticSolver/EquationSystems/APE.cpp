///////////////////////////////////////////////////////////////////////////////
//
// File APE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Kilian Lackhove
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
// Description: APE1/APE4 (Acoustic Perturbation Equations)
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <AcousticSolver/EquationSystems/APE.h>

using namespace std;

namespace Nektar
{
string APE::className = GetEquationSystemFactory().RegisterCreatorFunction(
    "APE", APE::create, "APE1/APE4 (Acoustic Perturbation Equations)");

APE::APE(const LibUtilities::SessionReaderSharedPtr &pSession,
         const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AcousticSystem(pSession, pGraph)
{
    m_ip   = 0;
    m_irho = -1;
    m_iu   = 1;

    m_conservative = false;
}

/**
 * @brief Initialization object for the APE class.
 */
void APE::v_InitObject()
{
    AcousticSystem::v_InitObject();

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
        boost::to_lower_copy(riemName) == "apeupwind" ||
        boost::to_lower_copy(riemName) == "upwind")
    {
        riemName = "APEUpwind";
    }
    if (boost::to_lower_copy(riemName) == "laxfriedrichs")
    {
        riemName = "APELaxFriedrichs";
    }
    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
        riemName, m_session);
    m_riemannSolver->SetVector("N", &APE::GetNormals, this);
    m_riemannSolver->SetVector("basefieldFwdBwd", &APE::GetBasefieldFwdBwd,
                               this);
    m_riemannSolver->SetAuxVec("vecLocs", &APE::GetVecLocs, this);

    // Set up advection operator
    string advName;
    m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    m_advection =
        SolverUtils::GetAdvectionFactory().CreateInstance(advName, advName);
    m_advection->SetFluxVector(&APE::v_GetFluxVector, this);
    m_advection->SetRiemannSolver(m_riemannSolver);
    m_advection->InitObject(m_session, m_fields);

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&APE::DoOdeRhs, this);
        m_ode.DefineProjection(&APE::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit APE not set up.");
    }
}

/**
 * @brief Destructor for APE class.
 */
APE::~APE()
{
}

/**
 * @brief Return the flux vector for the APE equations.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux. flux[eq][dir][pt]
 */
void APE::v_GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    int nq = physfield[0].size();
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    ASSERTL1(flux[0].size() == m_spacedim,
             "Dimension of flux array and velocity array do not match");

    // F_{adv,p',j} = \bar{rho}  \bar{c^2} u'_j + p' \bar{u}_j
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Zero(nq, flux[0][j], 1);

        // construct \bar{rho}  \bar{c^2} u'_j
        Vmath::Vmul(nq, m_bf[0], 1, m_bf[1], 1, tmp1, 1);
        Vmath::Vmul(nq, tmp1, 1, physfield[j + 1], 1, tmp1, 1);

        // construct p' \bar{u}_j term
        Vmath::Vmul(nq, physfield[0], 1, m_bf[j + 2], 1, tmp2, 1);

        // add both terms
        Vmath::Vadd(nq, tmp1, 1, tmp2, 1, flux[0][j], 1);
    }

    for (int i = 1; i < flux.size(); ++i)
    {
        ASSERTL1(flux[i].size() == m_spacedim,
                 "Dimension of flux array and velocity array do not match");

        // F_{adv,u'_i,j} = (p'/ \bar{rho} + \bar{u}_k u'_k) \delta_{ij}
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Zero(nq, flux[i][j], 1);

            if (i - 1 == j)
            {
                // contruct p'/ \bar{rho} term
                Vmath::Vdiv(nq, physfield[0], 1, m_bf[1], 1, flux[i][j], 1);

                // construct \bar{u}_k u'_k term
                Vmath::Zero(nq, tmp1, 1);
                for (int k = 0; k < m_spacedim; ++k)
                {
                    Vmath::Vvtvp(nq, physfield[k + 1], 1, m_bf[k + 2], 1, tmp1,
                                 1, tmp1, 1);
                }

                // add terms
                Vmath::Vadd(nq, flux[i][j], 1, tmp1, 1, flux[i][j], 1);
            }
        }
    }
}

/**
 * @brief Outflow characteristic boundary conditions for compressible
 * flow problems.
 */
void APE::v_RiemannInvariantBC(int bcRegion, int cnt,
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
        Array<OneD, NekDouble> Vn(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &Fwd[m_iu + i][id2], 1,
                         &m_traceNormals[i][id2], 1, &Vn[0], 1, &Vn[0], 1);
        }

        // Calculate (v0.n)
        Array<OneD, NekDouble> Vn0(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &BfFwd[2 + i][id2], 1,
                         &m_traceNormals[i][id2], 1, &Vn0[0], 1, &Vn0[0], 1);
        }

        for (int i = 0; i < nBCEdgePts; ++i)
        {
            NekDouble c = sqrt(BfFwd[0][id2 + i]);

            // LODI
            NekDouble h1, h2;

            // outgoing
            if (Vn0[i] - c > 0)
            {
                // u/2 - p/(2*rho0*sqr(c0sq))
                h1 = Vn[i] / 2 -
                     Fwd[m_ip][id2 + i] / (2 * BfFwd[1][id2 + i] * c);
            }
            // incoming
            else
            {
                h1 = 0.0;
            }

            // outgoing
            if (Vn0[i] + c > 0)
            {
                // u/2 + p/(2*rho0*sqr(c0sq))
                h2 = Vn[i] / 2 +
                     Fwd[m_ip][id2 + i] / (2 * BfFwd[1][id2 + i] * c);
            }
            // incoming
            else
            {
                h2 = 0.0;
            }

            // compute primitive variables
            // p = rho0*sqr(c0sq) * (h2 - h1)
            Fwd[m_ip][id2 + i] = BfFwd[1][id2 + i] * c * (h2 - h1);
            // u = h1 + h2
            NekDouble VnNew = h1 + h2;

            // adjust velocity pert. according to new value
            for (int j = 0; j < m_spacedim; ++j)
            {
                Fwd[m_iu + j][id2 + i] =
                    Fwd[m_iu + j][id2 + i] +
                    (VnNew - Vn[i]) * m_traceNormals[j][id2 + i];
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
