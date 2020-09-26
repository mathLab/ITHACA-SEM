///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFEAxisym.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFEAxisym.h>

using namespace std;

namespace Nektar
{
    string NavierStokesCFEAxisym::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFEAxisym", NavierStokesCFEAxisym::create,
            "Axisymmetric NavierStokes equations in conservative variables.");

    NavierStokesCFEAxisym::NavierStokesCFEAxisym(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          NavierStokesCFE(pSession, pGraph)
    {
    }

    NavierStokesCFEAxisym::~NavierStokesCFEAxisym()
    {

    }

    void NavierStokesCFEAxisym::v_InitObject()
    {
        NavierStokesCFE::v_InitObject();

        int nVariables = m_fields.size();
        int npoints    = GetNpoints();
        m_viscousForcing = Array<OneD, Array<OneD, NekDouble>> (nVariables);
        for (int i = 0; i < nVariables; ++i)
        {
            m_viscousForcing[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }
    }

    void NavierStokesCFEAxisym::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        int npoints    = GetNpoints();
        int nvariables = inarray.size();

        NavierStokesCFE::v_DoDiffusion(inarray, outarray, pFwd, pBwd);

        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(npoints,
                        m_viscousForcing[i], 1,
                        outarray[i], 1,
                        outarray[i], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFEAxisym::v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.size();
        int nPts       = physfield[0].size();

        // 1/r
        Array<OneD, Array<OneD, NekDouble> > coords(3);
        Array<OneD, NekDouble> invR (nPts,0.0);
        for (int i = 0; i < 3; i++)
        {
            coords[i] = Array<OneD, NekDouble> (nPts);
        }
        m_fields[0]->GetCoords(coords[0], coords[1], coords[2]);
        for (int i = 0; i < nPts; ++i)
        {
            if (coords[0][i] < NekConstants::kNekZeroTol)
            {
                invR[i] = 0;
            }
            else
            {
                invR[i] = 1.0/coords[0][i];
            }
        }

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > divVel             (nPts, 0.0);
        Array<OneD, NekDouble > tmp                (nPts, 0.0);

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nVariables-2], m_mu,
            m_thermalConductivity);

        // Velocity divergence = d(u_r)/dr + d(u_z)/dz + u_r/r
        Vmath::Vadd(nPts, derivativesO1[0][0], 1, derivativesO1[1][1], 1,
                        divVel, 1);
        Vmath::Vvtvp(nPts, physfield[0], 1 , invR, 1, divVel, 1, divVel, 1);

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, m_mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)

        for (i = 0; i < 2; ++i)
        {
            for (j = i; j < 2; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, m_mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }
        // Swirl case
        if(m_spacedim == 3)
        {
            // Tau_theta_theta = mu ( 2*u_r/r - 2/3*div(u))
            Vmath::Vmul(nPts, physfield[0], 1 , invR, 1,
                              viscousTensor[2][3], 1);
            Vmath::Smul(nPts, 2.0, viscousTensor[2][3], 1,
                              viscousTensor[2][3], 1);
            Vmath::Vmul(nPts, m_mu, 1, viscousTensor[2][3], 1,
                              viscousTensor[2][3], 1);
            Vmath::Vadd(nPts, viscousTensor[2][3], 1,
                              divVel, 1,
                              viscousTensor[2][3], 1);

            // Tau_r_theta = mu (-u_theta/r + d(u_theta)/dr )
            Vmath::Vmul(nPts, physfield[2], 1 , invR, 1,
                              viscousTensor[2][1], 1);
            Vmath::Smul(nPts, -1.0, viscousTensor[2][1], 1,
                              viscousTensor[2][1], 1);
            Vmath::Vadd(nPts, derivativesO1[0][2], 1 , viscousTensor[2][1], 1,
                              viscousTensor[2][1], 1);
            Vmath::Vmul(nPts, m_mu, 1, viscousTensor[2][1], 1,
                              viscousTensor[2][1], 1);
            Vmath::Vcopy(nPts, viscousTensor[2][1], 1,
                                       viscousTensor[0][3], 1);

            // Tau_z_theta = mu (d(u_theta)/dz )
            Vmath::Vmul(nPts, m_mu, 1, derivativesO1[1][2], 1,
                              viscousTensor[2][2], 1);
            Vmath::Vcopy(nPts, viscousTensor[2][2], 1,
                                       viscousTensor[1][3], 1);
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            if (i != 2)
            {
                Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                                   derivativesO1[i][m_spacedim], 1,
                                   viscousTensor[i][m_spacedim+1], 1,
                                   viscousTensor[i][m_spacedim+1], 1);
            }
            else
            {
                Vmath::Vmul(nPts, derivativesO1[i][m_spacedim], 1 ,
                                    invR, 1, tmp, 1);
                Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                                   tmp, 1,
                                   viscousTensor[i][m_spacedim+1], 1,
                                   viscousTensor[i][m_spacedim+1], 1);
            }
        }

        // Update viscous forcing
        //   r-momentum: F = 1/r * (tau_rr - tau_theta_theta)
        if(m_spacedim == 3)
        {
            Vmath::Vsub(nPts, viscousTensor[0][1], 1, viscousTensor[2][3], 1,
                                m_viscousForcing[1], 1);
            Vmath::Vmul(nPts, m_viscousForcing[1], 1 ,
                                    invR, 1, m_viscousForcing[1], 1);
        }
        else
        {
            Vmath::Vmul(nPts, viscousTensor[0][1], 1 ,
                                    invR, 1, m_viscousForcing[1], 1);
        }

        //   z-momentum: F = 1/r * tau_r_z
        Vmath::Vmul(nPts, viscousTensor[0][2], 1 ,
                                    invR, 1, m_viscousForcing[2], 1);

        //  Theta_momentum: F = 2* tau_r_theta
        if(m_spacedim == 3)
        {
            Vmath::Vmul(nPts, viscousTensor[0][3], 1 ,
                                    invR, 1, m_viscousForcing[3], 1);
            Vmath::Smul(nPts, 2.0, m_viscousForcing[3], 1,
                                    m_viscousForcing[3], 1);
        }

        // Energy: F = 1/r* viscousTensor_T_r
        Vmath::Vmul(nPts, viscousTensor[0][m_spacedim+1], 1 ,
                                    invR, 1, m_viscousForcing[m_spacedim+1], 1);
    }
}
