///////////////////////////////////////////////////////////////////////////////
//
// File IsentropicVortex.cpp
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
// Description: Euler equations for isentropic vortex
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/EquationSystems/IsentropicVortex.h>

using namespace std;

namespace Nektar
{
    string IsentropicVortex::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "IsentropicVortex", IsentropicVortex::create,
        "Euler equations for the isentropic vortex test case.");

    IsentropicVortex::IsentropicVortex(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          EulerCFE(pSession, pGraph)
    {
    }

    /**
     * @brief Destructor for EulerCFE class.
     */
    IsentropicVortex::~IsentropicVortex()
    {
    }

    /**
     * @brief Print out a summary with some relevant information.
     */
    void IsentropicVortex::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type", "Isentropic Vortex");
    }

    /**
     * @brief Set the initial conditions.
     */
    void IsentropicVortex::v_SetInitialConditions(
        NekDouble   initialtime,
        bool        dumpInitialConditions,
        const int   domain)
    {
        boost::ignore_unused(domain);

        int nTotQuadPoints  = GetTotPoints();
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);
        Array<OneD, Array<OneD, NekDouble> > u(m_spacedim+2);

        m_fields[0]->GetCoords(x, y, z);

        for (int i = 0; i < m_spacedim + 2; ++i)
        {
            u[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }

        EvaluateIsentropicVortex(x, y, z, u, initialtime);

        // Forward transform to fill the coefficient space
        for(int i = 0; i < m_fields.size(); ++i)
        {
            Vmath::Vcopy(nTotQuadPoints, u[i], 1, m_fields[i]->UpdatePhys(), 1);
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                                  m_fields[i]->UpdateCoeffs());
        }

        if (dumpInitialConditions && m_checksteps)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }

    /**
     * @brief Get the exact solutions for isentropic vortex and Ringleb
     * flow problems.
     */
    void IsentropicVortex::v_EvaluateExactSolution(
        unsigned int                         field,
        Array<OneD, NekDouble>              &outfield,
        const NekDouble                      time)
    {
        int nTotQuadPoints  = GetTotPoints();
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);
        Array<OneD, Array<OneD, NekDouble> > u(m_spacedim+2);

        m_fields[0]->GetCoords(x, y, z);

        for (int i = 0; i < m_spacedim + 2; ++i)
        {
            u[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }

        EvaluateIsentropicVortex(x, y, z, u, time);

        Vmath::Vcopy(nTotQuadPoints, u[field], 1, outfield, 1);
    }

    void IsentropicVortex::EvaluateIsentropicVortex(
        const Array<OneD, NekDouble>               &x,
        const Array<OneD, NekDouble>               &y,
        const Array<OneD, NekDouble>               &z,
        Array<OneD, Array<OneD, NekDouble> > &u,
        NekDouble                             time,
        const int                                   o)
    {
        boost::ignore_unused(z);

        int nq = x.size();

        // Flow parameters
        const NekDouble x0    = 5.0;
        const NekDouble y0    = 0.0;
        const NekDouble beta  = 5.0;
        const NekDouble u0    = 1.0;
        const NekDouble v0    = 0.5;
        const NekDouble gamma = m_gamma;
        NekDouble r, xbar, ybar, tmp;
        NekDouble fac = 1.0/(16.0*gamma*M_PI*M_PI);

        // In 3D zero rhow field.
        if (m_spacedim == 3)
        {
            Vmath::Zero(nq, &u[3][o], 1);
        }

        // Fill storage
        for (int i = 0; i < nq; ++i)
        {
            xbar      = x[i] - u0*time - x0;
            ybar      = y[i] - v0*time - y0;
            r         = sqrt(xbar*xbar + ybar*ybar);
            tmp       = beta*exp(1-r*r);
            u[0][i+o] = pow(1.0 - (gamma-1.0)*tmp*tmp*fac, 1.0/(gamma-1.0));
            u[1][i+o] = u[0][i+o]*(u0 - tmp*ybar/(2*M_PI));
            u[2][i+o] = u[0][i+o]*(v0 + tmp*xbar/(2*M_PI));
            u[m_spacedim+1][i+o] = pow(u[0][i+o], gamma)/(gamma-1.0) +
            0.5*(u[1][i+o]*u[1][i+o] + u[2][i+o]*u[2][i+o]) / u[0][i+o];
        }
    }

}
