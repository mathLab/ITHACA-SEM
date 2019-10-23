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

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/EquationSystems/RinglebFlow.h>

using namespace std;

namespace Nektar
{
    string RinglebFlow::className = 
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "RinglebFlow", RinglebFlow::create, 
        "Euler equations for Ringleb flow.");
    
    RinglebFlow::RinglebFlow(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          EulerCFE(pSession, pGraph)
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
        boost::ignore_unused(time);
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
}
