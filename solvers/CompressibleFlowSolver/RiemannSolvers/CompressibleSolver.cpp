///////////////////////////////////////////////////////////////////////////////
//
// File: CompressibleSolver.cpp
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
// Description: Compressible Riemann solver.
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/RiemannSolvers/CompressibleSolver.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
//#include <LocalRegions/MatrixKey.h>


namespace Nektar
{
    CompressibleSolver::CompressibleSolver(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : RiemannSolver(pSession), m_pointSolve(true)
    {
        m_requiresRotation = true;

        // Create equation of state object
        std::string eosType;
        pSession->LoadSolverInfo("EquationOfState",
                                  eosType, "IdealGas");
        m_eos = GetEquationOfStateFactory()
                                .CreateInstance(eosType, pSession);
        // Check if using ideal gas
        m_idealGas = boost::iequals(eosType,"IdealGas");
    }
    
    void CompressibleSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
    {
        if (m_pointSolve)
        {
            int expDim      = nDim;
            int nvariables  = Fwd.num_elements();
            
            NekDouble rhouf, rhovf;
            
            // Check if PDE-based SC is used
            if (expDim == 1)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], 0.0,   0.0,   Fwd [2][i],
                        Bwd [0][i], Bwd [1][i], 0.0,   0.0,   Bwd [2][i],
                        flux[0][i], flux[1][i], rhouf, rhovf, flux[2][i]);
                }
            }
            else if (expDim == 2)
            {
                if (nvariables == expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolve(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0,   Fwd [3][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0,   Bwd [3][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i]);
                    }
                }
                
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], 0.0, Fwd [3][i], Fwd [4][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], 0.0, Bwd [3][i], Bwd [4][i],
                            flux[0][i], flux[1][i], flux[2][i], rhovf, flux[3][i], flux[4][i]);
                    }
                }
                
            }
            else if (expDim == 3)
            {
                for (int i = 0; i < Fwd[0].num_elements(); ++i)
                {
                    v_PointSolve(
                        Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i],
                        Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i],
                        flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i]);
                }
                if (nvariables > expDim+2)
                {
                    for (int i = 0; i < Fwd[0].num_elements(); ++i)
                    {
                        v_PointSolveVisc(
                            Fwd [0][i], Fwd [1][i], Fwd [2][i], Fwd [3][i], Fwd [4][i], Fwd [5][i],
                            Bwd [0][i], Bwd [1][i], Bwd [2][i], Bwd [3][i], Bwd [4][i], Bwd [5][i],
                            flux[0][i], flux[1][i], flux[2][i], flux[3][i], flux[4][i], flux[5][i]);
                    }
                }
            }
        }
        else
        {
            v_ArraySolve(Fwd, Bwd, flux);
        }
    }

    void CompressibleSolver::v_CalcFluxJacobian(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
        const Array<OneD, const Array<OneD, NekDouble> > &normals,
              DNekBlkMatSharedPtr                         FJac,
              DNekBlkMatSharedPtr                         BJac)
    {
        int expDim      = nDim;
        int nvariables  = Fwd.num_elements();
        int nnomals     = normals.num_elements();
        // int nvariables3D = nvariables+2;
        int nvariables3D = 5;

        if (nvariables > expDim+2)
        {
            ASSERTL0(false,"nvariables > expDim+2 case not coded")
        }

        DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables3D, nvariables3D);
        DNekMatSharedPtr PointBJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables3D, nvariables3D);
        
        Array<OneD, NekDouble> PointFwd(nvariables3D,0.0),PointBwd(nvariables3D,0.0);
        Array<OneD, NekDouble> PointNormal(3,0.0);

        Array<OneD, unsigned int> index(nvariables);

        index[nvariables-1] = 4;
        for(int i=0;i<nvariables-1;i++)
        {
            index[i] = i;
        }

        int nj=0;
        int nk=0;
        for (int i = 0; i < Fwd[0].num_elements(); ++i)
        {
            for(int j=0; j< nnomals; j++)
            {
                PointNormal[j] = normals [j][i];
            }

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                PointFwd[nj] = Fwd [j][i];
                PointBwd[nj] = Bwd [j][i];
            }
            
            v_PointFluxJacobian(PointFwd,PointBwd,PointNormal,PointFJac3D,PointBJac3D);

            DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables, nvariables);
            DNekMatSharedPtr PointBJac = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables, nvariables);

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                for(int k=0; k< nvariables; k++)
                {
                    nk = index[k];
                    (*PointFJac)(j,k) = (*PointFJac3D)(nj,nk); 
                    (*PointBJac)(j,k) = (*PointBJac3D)(nj,nk); 
                }
            }

            FJac->SetBlock(i, i, PointFJac);
            BJac->SetBlock(i, i, PointBJac);
        }
    }
   

    NekDouble CompressibleSolver::GetRoeSoundSpeed(
        NekDouble rhoL, NekDouble pL, NekDouble eL, NekDouble HL, NekDouble srL,
        NekDouble rhoR, NekDouble pR, NekDouble eR, NekDouble HR, NekDouble srR,
        NekDouble HRoe, NekDouble URoe2, NekDouble srLR)
    {
        static NekDouble gamma = m_params["gamma"]();
        NekDouble cRoe;
        if(m_idealGas)
        {
            cRoe = sqrt((gamma - 1.0)*(HRoe - 0.5 * URoe2));
        }
        else
        {
            // Calculate static enthalpy of left and right states
            NekDouble hL = eL + pL/rhoL;
            NekDouble hR = eR + pR/rhoR;

            // Get partial derivatives of P(rho,e)
            NekDouble dpdeL   = m_eos->GetDPDe_rho(rhoL,eL);
            NekDouble dpdeR   = m_eos->GetDPDe_rho(rhoR,eR);
            NekDouble dpdrhoL = m_eos->GetDPDrho_e(rhoL,eL);
            NekDouble dpdrhoR = m_eos->GetDPDrho_e(rhoR,eR);

            // Evaluate chi and kappa parameters
            NekDouble chiL    = dpdrhoL - eL / rhoL * dpdeL;
            NekDouble kappaL  = dpdeL / rhoL;
            NekDouble chiR    = dpdrhoR - eR / rhoR * dpdeR;
            NekDouble kappaR  = dpdeR / rhoR;

            //
            // Calculate interface speed of sound using procedure from
            //    Vinokur, M.; Montagné, J.-L. "Generalized Flux-Vector
            //    Splitting and Roe Average for an Equilibrium Real Gas",
            //    JCP (1990).
            //

            // Calculate averages
            NekDouble avgChi    = 0.5 * (chiL      + chiR);
            NekDouble avgKappa  = 0.5 * (kappaL    + kappaR);
            NekDouble avgKappaH = 0.5 * (kappaL*hL + kappaR*hR);

            // Calculate jumps
            NekDouble deltaP    = pR      - pL;
            NekDouble deltaRho  = rhoR    - rhoL;
            NekDouble deltaRhoe = rhoR*eR - rhoL*eL;

            // Evaluate dP: equation (64) from Vinokur-Montagné
            NekDouble dP = deltaP - avgChi * deltaRho - avgKappa * deltaRhoe;
            // s (eq 66)
            NekDouble s  = avgChi + avgKappaH;
            // D (eq 65)
            NekDouble D  = (s*deltaRho)*(s*deltaRho) + deltaP*deltaP;
            // chiRoe and kappaRoe (eq 66)
            NekDouble chiRoe, kappaRoe;
            NekDouble fac = D - deltaP*deltaRho;
            if( abs(fac) > NekConstants::kNekZeroTol)
            {
                chiRoe   = (D*avgChi + s*s*deltaRho*dP) / fac;
                kappaRoe = D*avgKappa / fac;
            }
            else
            {
                chiRoe = avgChi;
                kappaRoe = avgKappa;
            }
            // Speed of sound (eq 53)
            cRoe = sqrt( chiRoe + kappaRoe*(HRoe - 0.5 * URoe2));
        }
        return cRoe;
    }

}
