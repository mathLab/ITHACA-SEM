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
// Description: Euler equations for isentropic vortex
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <CompressibleFlowSolver/EquationSystems/IsentropicVortex.h>

using namespace std;

namespace Nektar
{
    string IsentropicVortex::className = 
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "IsentropicVortex", IsentropicVortex::create, 
        "Euler equations for the isentropic vortex test case.");
    
    IsentropicVortex::IsentropicVortex(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : EulerCFE(pSession)
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
        CompressibleFlowSystem::v_SetInitialConditions();

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
        for(int i = 0; i < m_fields.num_elements(); ++i)
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

    /**
     * @brief Set the boundary conditions for the isentropic vortex problem.
     */
    void IsentropicVortex::v_SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {
        int cnt        = 0;
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.num_elements();

        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        const Array<OneD, const int> &bndTraceMap = m_fields[0]->GetTraceBndMap();
        // loop over Boundary Regions
        int id2, e_max;
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            e_max = m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();

            for(int e = 0; e < e_max; ++e)
            {
                int npoints = m_fields[0]->
                    GetBndCondExpansions()[n]->GetExp(e)->GetTotPoints();
                int id1  = m_fields[0]->
                    GetBndCondExpansions()[n]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(bndTraceMap[cnt++]);

                Array<OneD,NekDouble> x(npoints, 0.0);
                Array<OneD,NekDouble> y(npoints, 0.0);
                Array<OneD,NekDouble> z(npoints, 0.0);

                m_fields[0]->GetBndCondExpansions()[n]->
                GetExp(e)->GetCoords(x, y, z);

                EvaluateIsentropicVortex(x, y, z, Fwd, time, id2);

                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, &Fwd[i][id2], 1, 
                                 &(m_fields[i]->GetBndCondExpansions()[n]->
                                 UpdatePhys())[id1], 1);
                }
            }
        }
    }

    void IsentropicVortex::EvaluateIsentropicVortex(
        const Array<OneD, NekDouble>               &x,
        const Array<OneD, NekDouble>               &y,
        const Array<OneD, NekDouble>               &z,
        Array<OneD, Array<OneD, NekDouble> > &u,
        NekDouble                             time,
        const int                                   o)
    {
        int nq = x.num_elements();
        
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
