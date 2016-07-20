///////////////////////////////////////////////////////////////////////////////
//
// File RinglebFlow.cpp
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
// Description: Euler equations for Ringleb flow
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <CompressibleFlowSolver/EquationSystems/RinglebFlow.h>

using namespace std;

namespace Nektar
{
    string RinglebFlow::className = 
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "RinglebFlow", RinglebFlow::create, 
        "Euler equations for Ringleb flow.");
    
    RinglebFlow::RinglebFlow(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : EulerCFE(pSession)
    {
    }

    /**
     * @brief Destructor for EulerCFE class.
     */
    RinglebFlow::~RinglebFlow()
    {
    }

    /**
     * @brief Print out a summary with some relevant information.
     */
    void RinglebFlow::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type", "RinglebFlow");
    }

    /**
     * @brief Get the exact solutions for isentropic vortex and Ringleb
     * flow problems.
     */
    void RinglebFlow::v_EvaluateExactSolution(
        unsigned int                         field,
        Array<OneD, NekDouble>              &outfield,
        const NekDouble                      time)
    {
        GetExactRinglebFlow( field, outfield);
    }

    /**
     * @brief Compute the exact solution for the Ringleb flow problem.
     */
    void RinglebFlow::GetExactRinglebFlow(
        int                              field, 
        Array<OneD, NekDouble>          &outarray)
    {
        int nTotQuadPoints  = GetTotPoints();
        
        Array<OneD, NekDouble> rho(nTotQuadPoints, 100.0);
        Array<OneD, NekDouble> rhou(nTotQuadPoints);
        Array<OneD, NekDouble> rhov(nTotQuadPoints);
        Array<OneD, NekDouble> E(nTotQuadPoints);
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);
        
        m_fields[0]->GetCoords(x, y, z);
        
        // Flow parameters
        NekDouble c, k, phi, r, J, VV, pp, sint, P, ss;
        NekDouble J11, J12, J21, J22, det;
        NekDouble Fx, Fy;
        NekDouble xi, yi;
        NekDouble dV;
        NekDouble dtheta;
        NekDouble par1;
        NekDouble theta     = M_PI / 4.0;
        NekDouble kExt      = 0.7;
        NekDouble V         = kExt * sin(theta);
        NekDouble toll      = 1.0e-8;
        NekDouble errV      = 1.0;
        NekDouble errTheta  = 1.0;
        NekDouble gamma     = m_gamma;
        NekDouble gamma_1_2 = (gamma - 1.0) / 2.0;
        
        for (int i = 0; i < nTotQuadPoints; ++i)
        {
            while ((abs(errV) > toll) || (abs(errTheta) > toll))
            {
                VV   = V * V;
                sint = sin(theta);
                c    = sqrt(1.0 - gamma_1_2 * VV);
                k    = V / sint;
                phi  = 1.0 / k;
                pp   = phi * phi;
                J    = 1.0 / c + 1.0 / (3.0 * c * c * c) + 
                1.0 / (5.0 * c * c * c * c * c) - 
                0.5 * log((1.0 + c) / (1.0 - c));
                
                r    = pow(c, 1.0 / gamma_1_2);
                xi   = 1.0 / (2.0 * r) * (1.0 / VV - 2.0 * pp) + J / 2.0;
                yi   = phi / (r * V) * sqrt(1.0 - VV * pp);
                par1 = 25.0 - 5.0 * VV;
                ss   = sint * sint;
                
                Fx   = xi - x[i];
                Fy   = yi - y[i];
                
                J11  = 39062.5 / pow(par1, 3.5) * (1.0 / VV - 2.0 / VV * ss) * 
                V + 1562.5 / pow(par1, 2.5) * (-2.0 / (VV * V) + 4.0 / 
                (VV * V) * ss) + 12.5 / pow(par1, 1.5) * V + 312.5 / 
                pow(par1, 2.5) * V + 7812.5 / pow(par1, 3.5) * V - 
                0.25 * (-1.0 / pow(par1, 0.5) * V/(1.0 - 0.2 * 
                pow(par1, 0.5)) - (1.0 + 0.2 * pow(par1, 0.5)) / 
                pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) / 
                pow(par1, 0.5) * V) / (1.0 + 0.2 * pow(par1, 0.5)) * 
                (1.0 - 0.2 * pow(par1, 0.5));
                
                J12  = -6250.0 / pow(par1, 2.5) / VV * sint * cos(theta);
                J21  = -6250.0 / (VV * V) * sint / 
                pow(par1, 2.5) * pow((1.0 - ss), 0.5) + 
                78125.0 / V * sint / pow(par1, 3.5) * 
                pow((1.0 - ss), 0.5);
                
                // the matrix is singular when theta = pi/2
                if(abs(y[i])<toll && abs(cos(theta))<toll)
                {
                    J22 = -39062.5 / pow(par1, 3.5) / V + 3125 / 
                    pow(par1, 2.5) / (VV * V) + 12.5 / pow(par1, 1.5) * 
                    V + 312.5 / pow(par1, 2.5) * V + 7812.5 / 
                    pow(par1, 3.5) * V - 0.25 * (-1.0 / pow(par1, 0.5) *
                    V / (1.0 - 0.2 * pow(par1, 0.5)) - (1.0 + 0.2 * 
                    pow(par1, 0.5)) / pow((1.0 - 0.2 * 
                    pow(par1, 0.5)), 2.0) / pow(par1, 0.5) * V) / 
                    (1.0 + 0.2 * pow(par1, 0.5)) * (1.0 - 0.2 * 
                    pow(par1,0.5));
                    
                    // dV = -dV/dx * Fx
                    dV      = -1.0 / J22 * Fx;
                    dtheta  = 0.0;
                    theta   = M_PI / 2.0;
                }
                else
                {
                    J22 = 3125.0 / VV * cos(theta) / pow(par1, 2.5) * 
                    pow((1.0 - ss), 0.5) - 3125.0 / VV * ss / 
                    pow(par1, 2.5) / pow((1.0 - ss), 0.5) * cos(theta);
                    
                    det = -1.0 / (J11 * J22 - J12 * J21);
                    
                    // [dV dtheta]' = -[invJ]*[Fx Fy]'
                    dV     = det * ( J22 * Fx - J12 * Fy);
                    dtheta = det * (-J21 * Fx + J11 * Fy);
                }
                
                V     = V + dV;
                theta = theta + dtheta;
                
                errV     = abs(dV);
                errTheta = abs(dtheta);
                
            }
            
            c = sqrt(1.0 - gamma_1_2 * VV);
            r = pow(c, 1.0 / gamma_1_2);
            
            rho[i]  = r;
            rhou[i] = rho[i] * V * cos(theta);
            rhov[i] = rho[i] * V * sin(theta);
            P       = (c * c) * rho[i] / gamma;
            E[i]    = P / (gamma - 1.0) + 0.5 * 
            (rhou[i] * rhou[i] / rho[i] + rhov[i] * rhov[i] / rho[i]);
            
            // Resetting the guess value
            errV     = 1.0;
            errTheta = 1.0;
            theta    = M_PI/4.0;
            V        = kExt*sin(theta);
        }
        
        switch (field)
        {
            case 0:
                outarray = rho;
                break;
            case 1:
                outarray = rhou;
                break;
            case 2:
                outarray = rhov;
                break;
            case 3:
                outarray = E;
                break;
            default:
                ASSERTL0(false, "Error in variable number!");
                break;
        }
    }

    /**
     * @brief Set the boundary conditions for the ringleb flow.
     */
    void RinglebFlow::v_SetBoundaryConditions(
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

        std::string userDefStr;
        int nreg = m_fields[0]->GetBndConditions().num_elements();
        // Loop over Boundary Regions
        for (int n = 0; n < nreg; ++n)
        {
            userDefStr = m_fields[0]->GetBndConditions()[n]->GetUserDefined();
            if (boost::iequals(userDefStr,"RinglebFlow"))
            {
                SetBoundaryRinglebFlow(n, time, cnt, Fwd, physarray);
            }
            else
            {
                // set up userdefined BC
                SetCommonBC(userDefStr, n, time, cnt, Fwd, physarray);
            }
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }

    /**
     * @brief Set the boundary conditions for the Ringleb flow problem.
     */
    void RinglebFlow::SetBoundaryRinglebFlow(
        int                                      bcRegion, 
        NekDouble                                time, 
        int                                      cnt, 
        Array<OneD, Array<OneD, NekDouble> >    &Fwd,
        Array<OneD, Array<OneD, NekDouble> >    &physarray)
    {
        int nvariables      = physarray.num_elements();
        
        // For 3DHomogenoeus1D
        int n_planes = 1;
        if (m_expdim == 2 &&  m_HomogeneousType == eHomogeneous1D)
        {
            int nPointsTot = m_fields[0]->GetTotPoints();
            int nPointsTot_plane = m_fields[0]->GetPlane(0)->GetTotPoints();
            n_planes = nPointsTot/nPointsTot_plane;
        }
        
        int id2, id2_plane, e_max;
        
        e_max = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
        
        for(int e = 0; e < e_max; ++e)
        {
            int npoints = m_fields[0]->
            GetBndCondExpansions()[bcRegion]->GetExp(e)->GetTotPoints();
            int id1  = m_fields[0]->
            GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
            
            // For 3DHomogenoeus1D
            if (m_expdim == 2 &&  m_HomogeneousType == eHomogeneous1D)
            {
                int cnt_plane = cnt/n_planes;
                int e_plane;
                int e_max_plane = e_max/n_planes;
                int nTracePts_plane = GetTraceTotPoints();
                
                int planeID = floor((e + 0.5 )/ e_max_plane );
                e_plane = e - e_max_plane*planeID;
                
                id2_plane  = m_fields[0]->GetTrace()->GetPhys_Offset(
                    m_fields[0]->GetTraceMap()->
                    GetBndCondCoeffsToGlobalCoeffsMap(
                        cnt_plane + e_plane));
                id2 = id2_plane + planeID*nTracePts_plane;
            }
            else // For general case
            {
                id2 = m_fields[0]->
                GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->
                GetBndCondTraceToGlobalTraceMap(cnt+e));
            }
            
            Array<OneD,NekDouble> x0(npoints, 0.0);
            Array<OneD,NekDouble> x1(npoints, 0.0);
            Array<OneD,NekDouble> x2(npoints, 0.0);
            
            m_fields[0]->GetBndCondExpansions()[bcRegion]->
            GetExp(e)->GetCoords(x0, x1, x2);
            
            // Flow parameters
            NekDouble c, k, phi, r, J, VV, pp, sint, P, ss;
            NekDouble J11, J12, J21, J22, det;
            NekDouble Fx, Fy;
            NekDouble xi, yi;
            NekDouble dV;
            NekDouble dtheta;
            NekDouble par1;
            NekDouble theta     = M_PI / 4.0;
            NekDouble kExt      = 0.7;
            NekDouble V         = kExt * sin(theta);
            NekDouble toll      = 1.0e-8;
            NekDouble errV      = 1.0;
            NekDouble errTheta  = 1.0;
            NekDouble gamma     = m_gamma;
            NekDouble gamma_1_2 = (gamma - 1.0) / 2.0;
            
            // Loop on all the points of that edge
            for (int j = 0; j < npoints; j++)
            {
                
                while ((abs(errV) > toll) || (abs(errTheta) > toll))
                {
                    VV   = V * V;
                    sint = sin(theta);
                    c    = sqrt(1.0 - gamma_1_2 * VV);
                    k    = V / sint;
                    phi  = 1.0 / k;
                    pp   = phi * phi;
                    J    = 1.0 / c + 1.0 / (3.0 * c * c * c) + 
                    1.0 / (5.0 * c * c * c * c * c) - 
                    0.5 * log((1.0 + c) / (1.0 - c));
                    
                    r    = pow(c, 1.0 / gamma_1_2);
                    xi   = 1.0 / (2.0 * r) * (1.0 / VV - 2.0 * pp) + J / 2.0;
                    yi   = phi / (r * V) * sqrt(1.0 - VV * pp);
                    par1 = 25.0 - 5.0 * VV;
                    ss   = sint * sint;
                    
                    Fx   = xi - x0[j];
                    Fy   = yi - x1[j];
                    
                    J11 = 39062.5 / pow(par1, 3.5) * 
                    (1.0 / VV - 2.0 / VV * ss) * V + 1562.5 / 
                    pow(par1, 2.5) * (-2.0 / (VV * V) + 4.0 / 
                    (VV * V) * ss) + 12.5 / pow(par1, 1.5) * V + 
                    312.5 / pow(par1, 2.5) * V + 7812.5 / 
                    pow(par1, 3.5) * V - 0.25 * 
                    (-1.0 / pow(par1, 0.5) * V / (1.0 - 0.2 * 
                    pow(par1, 0.5)) - (1.0 + 0.2 * pow(par1, 0.5)) / 
                    pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) / 
                    pow(par1, 0.5) * V) / (1.0 + 0.2 * pow(par1, 0.5)) *
                    (1.0 - 0.2 * pow(par1, 0.5));
                    
                    J12 = -6250.0 / pow(par1, 2.5) / VV * sint * cos(theta);
                    J21 = -6250.0 / (VV * V) * sint / pow(par1, 2.5) * 
                    pow((1.0 - ss), 0.5) + 78125.0 / V * sint / 
                    pow(par1, 3.5) * pow((1.0 - ss), 0.5);
                    
                    // the matrix is singular when theta = pi/2
                    if (abs(x1[j]) < toll && abs(cos(theta)) < toll)
                    {
                        J22 = -39062.5 / pow(par1, 3.5) / V + 3125 / 
                        pow(par1, 2.5) / (VV * V) + 12.5 / 
                        pow(par1, 1.5) * V + 312.5 / pow(par1, 2.5) * 
                        V + 7812.5 / pow(par1, 3.5) * V - 0.25 * 
                        (-1.0 / pow(par1, 0.5) * V / (1.0 - 0.2 * 
                        pow(par1, 0.5)) - (1.0 + 0.2 * pow(par1, 0.5)) /
                        pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) / 
                        pow(par1, 0.5) * V) / (1.0 + 0.2 * 
                        pow(par1, 0.5)) * (1.0 - 0.2 * pow(par1, 0.5));
                        
                        // dV = -dV/dx * Fx
                        dV      = -1.0 / J22 * Fx;
                        dtheta  = 0.0;
                        theta   = M_PI / 2.0;
                    }
                    else
                    {
                        J22 = 3125.0 / VV * cos(theta) / pow(par1, 2.5) * 
                        pow((1.0 - ss), 0.5) - 3125.0 / VV * ss / 
                        pow(par1, 2.5) / pow((1.0 - ss), 0.5) * 
                        cos(theta);
                        
                        det = -1.0 / (J11 * J22 - J12 * J21);
                        
                        // [dV dtheta]' = -[invJ]*[Fx Fy]'
                        dV     = det * ( J22 * Fx - J12 * Fy);
                        dtheta = det * (-J21 * Fx + J11 * Fy);
                    }
                    
                    V     = V + dV;
                    theta = theta + dtheta;
                    
                    errV     = abs(dV);
                    errTheta = abs(dtheta);
                }
                
                c                      = sqrt(1.0 - gamma_1_2 * VV);
                int kk                 = id2 + j;
                NekDouble timeramp     = 200.0;
                std::string restartstr = "RESTART";
                if (time<timeramp &&
                    !(m_session->DefinesFunction("InitialConditions") &&
                    m_session->GetFunctionType("InitialConditions", 0) ==
                    LibUtilities::eFunctionTypeFile))
                {
                    Fwd[0][kk] = pow(c, 1.0 / gamma_1_2) * 
                    exp(-1.0 + time /timeramp);
                    
                    Fwd[1][kk] = Fwd[0][kk] * V * cos(theta) * 
                    exp(-1 + time / timeramp);
                    
                    Fwd[2][kk] = Fwd[0][kk] * V * sin(theta) * 
                    exp(-1 + time / timeramp);
                }
                else
                {
                    Fwd[0][kk]  = pow(c, 1.0 / gamma_1_2);
                    Fwd[1][kk] = Fwd[0][kk] * V * cos(theta);
                    Fwd[2][kk] = Fwd[0][kk] * V * sin(theta);
                }
                
                P  = (c * c) * Fwd[0][kk] / gamma;
                Fwd[3][kk]  = P / (gamma - 1.0) + 0.5 * 
                (Fwd[1][kk] * Fwd[1][kk] / Fwd[0][kk] + 
                Fwd[2][kk] * Fwd[2][kk] / Fwd[0][kk]);
                
                errV     = 1.0;
                errTheta = 1.0;
                theta    = M_PI / 4.0;
                V        = kExt * sin(theta);
            }
            
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, &Fwd[i][id2], 1, 
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1],1);
            }
        }
    }
}
