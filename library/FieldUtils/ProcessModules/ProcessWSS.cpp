////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWSS.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes wss field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include "ProcessWSS.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessWSS::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "wss"),
    ProcessWSS::create,
    "Computes wall shear stress field.");

ProcessWSS::ProcessWSS(FieldSharedPtr f) : ProcessBoundaryExtract(f)
{
}

ProcessWSS::~ProcessWSS()
{
}

void ProcessWSS::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);

    int i, j;
    int nfields = m_f->m_variables.size();
    int expdim  = m_f->m_graph->GetSpaceDimension();
    m_spacedim  = expdim + m_f->m_numHomogeneousDir;


    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    if (m_spacedim == 1)
    {
        ASSERTL0(false, "Error: wss for a 1D problem cannot "
                        "be computed");
    }

    // Declare arrays
    int nshear    = m_spacedim + 1;
    int nstress   = m_spacedim * m_spacedim;
    int ngrad     = m_spacedim * m_spacedim;

    Array<OneD, Array<OneD, NekDouble> > velocity(nfields);
    Array<OneD, Array<OneD, NekDouble> > grad(ngrad);
    Array<OneD, Array<OneD, NekDouble> > stress(nstress), fstress(nstress);
    Array<OneD, Array<OneD, NekDouble> > fshear(nshear);

    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(nshear);
    Array<OneD, MultiRegions::ExpListSharedPtr> BndElmtExp(nfields);

    // will resuse nfields expansions to write shear components.
    if(nshear > nfields)
    {
        m_f->m_exp.resize(nshear);
        for (i = nfields; i < nshear; ++i)
        {
            m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        }
    }

    // Create map of boundary ids for partitioned domains
    SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                           m_f->m_exp[0]->GetGraph());
    const SpatialDomains::BoundaryRegionCollection bregions =
        bcs.GetBoundaryRegions();
    map<int, int> BndRegionMap;
    int cnt = 0;
    for (auto &breg_it : bregions)
    {
        BndRegionMap[breg_it.first] = cnt++;
    }

    // Loop over boundaries to Write
    for (int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
    {
        if (BndRegionMap.count(m_f->m_bndRegionsToWrite[b]) == 1)
        {
            int bnd = BndRegionMap[m_f->m_bndRegionsToWrite[b]];
            // Get expansion list for boundary and for elements containing this
            // bnd
            for (i = 0; i < nshear; i++)
            {
                BndExp[i] = m_f->m_exp[i]->UpdateBndCondExpansion(bnd);
            }
            for (i = 0; i < nfields; i++)
            {
                m_f->m_exp[i]->GetBndElmtExpansion(bnd, BndElmtExp[i]);
            }

            // Get number of points in expansions
            int nqb = BndExp[0]->GetTotPoints();
            int nqe = BndElmtExp[0]->GetTotPoints();

            // Initialise local arrays for the velocity gradients, and stress
            // components
            // size of total number of quadrature points for elements in this
            // bnd
            for (i = 0; i < ngrad; ++i)
            {
                grad[i] = Array<OneD, NekDouble>(nqe);
            }

            for (i = 0; i < nstress; ++i)
            {
                stress[i] = Array<OneD, NekDouble>(nqe);
            }

            Array<OneD, NekDouble> div(nqe, 0.0);

            // initialise arrays in the boundary
            for (i = 0; i < nstress; ++i)
            {
                fstress[i] = Array<OneD, NekDouble>(nqb);
            }

            for (i = 0; i < nshear; ++i)
            {
                fshear[i] = Array<OneD, NekDouble>(nqb, 0.0);
            }

            // Extract Velocities
            GetVelocity( BndElmtExp, velocity);

            // Extract viscosity coefficients
            NekDouble lambda;
            Array<OneD, NekDouble> mu(nqe, 0.0);
            GetViscosity( BndElmtExp, mu, lambda);

            // Compute gradients
            for (i = 0; i < m_spacedim; ++i)
            {
                if (m_spacedim == 2)
                {
                    BndElmtExp[i]->PhysDeriv(velocity[i],
                                             grad[i * m_spacedim + 0],
                                             grad[i * m_spacedim + 1]);
                }
                else
                {
                    BndElmtExp[i]->PhysDeriv(velocity[i],
                                             grad[i * m_spacedim + 0],
                                             grad[i * m_spacedim + 1],
                                             grad[i * m_spacedim + 2]);
                }
                // Add contribution to div(u)
                Vmath::Vadd(nqe, grad[i * m_spacedim + i], 1, div, 1, div, 1);
            }

            // Velocity divergence scaled by lambda * mu
            Vmath::Smul(nqe, lambda, div, 1, div, 1);
            Vmath::Vmul(nqe, mu, 1, div, 1, div, 1);

            // Compute stress component terms
            //    tau_ij = mu*(u_i,j + u_j,i) + mu*lambda*delta_ij*div(u)
            for (i = 0; i < m_spacedim; ++i)
            {
                for (j = i; j < m_spacedim; ++j)
                {
                    Vmath::Vadd(nqe, grad[i * m_spacedim + j], 1,
                                     grad[j * m_spacedim + i], 1,
                                     stress[i * m_spacedim + j], 1);

                    Vmath::Vmul(nqe, mu, 1,
                                     stress[i * m_spacedim + j], 1,
                                     stress[i * m_spacedim + j], 1);

                    if (i == j)
                    {
                        // Add divergence term to diagonal
                        Vmath::Vadd(nqe, stress[i * m_spacedim + j], 1,
                                         div, 1,
                                         stress[i * m_spacedim + j], 1);
                    }
                    else
                    {
                        // Copy to make symmetric
                        Vmath::Vcopy(nqe, stress[i * m_spacedim + j], 1,
                                          stress[j * m_spacedim + i], 1);
                    }
                }
            }

            // Get boundary stress values.
            for (j = 0; j < nstress; ++j)
            {
                m_f->m_exp[0]->ExtractElmtToBndPhys(bnd, stress[j], fstress[j]);
            }

            // Get normals
            Array<OneD, Array<OneD, NekDouble> > normals;
            m_f->m_exp[0]->GetBoundaryNormals(bnd, normals);
            // Reverse normals, to get correct orientation for the body
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Neg(nqb, normals[i], 1);
            }

            // calculate wss, and update coeffs in the boundary expansion
            // S = tau_ij * n_j
            for (i = 0; i < m_spacedim; ++i)
            {
                for (j = 0; j < m_spacedim; ++j)
                {
                    Vmath::Vvtvp(nqb, normals[j], 1,
                                      fstress[i * m_spacedim + j], 1,
                                      fshear[i], 1, fshear[i], 1);
                }
            }

            // T = S - (S.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nqb, normals[i], 1, fshear[i], 1,
                             fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            }
            Vmath::Smul(nqb, -1.0, fshear[nshear - 1], 1,
                                   fshear[nshear - 1], 1);

            for (i = 0; i < m_spacedim; i++)
            {
                Vmath::Vvtvp(nqb, normals[i], 1, fshear[nshear - 1], 1,
                             fshear[i], 1, fshear[i], 1);
                BndExp[i]->FwdTrans_IterPerExp(fshear[i],
                                        BndExp[i]->UpdateCoeffs());
            }

            // Tw
            Vmath::Zero(nqb, fshear[nshear - 1], 1);
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(nqb, fshear[i], 1, fshear[i], 1,
                             fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            }
            Vmath::Vsqrt(nqb, fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            BndExp[nshear - 1]->FwdTrans_IterPerExp(fshear[nshear - 1],
                                         BndExp[nshear - 1]->UpdateCoeffs());
        }
    }

    if (m_spacedim == 2)
    {
        m_f->m_variables[0] = "Shear_x";
        m_f->m_variables[1] = "Shear_y";
        m_f->m_variables[2] = "Shear_mag";
    }
    else
    {
        m_f->m_variables[0] = "Shear_x";
        m_f->m_variables[1] = "Shear_y";
        m_f->m_variables[2] = "Shear_z";
        m_f->m_variables[3] = "Shear_mag";
    }
}

void ProcessWSS::GetViscosity(
        const Array<OneD, MultiRegions::ExpListSharedPtr> exp,
              Array<OneD, NekDouble> &mu,
              NekDouble &lambda)
{
    NekDouble m_mu;
    int npoints = exp[0]->GetNpoints();

    if(boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        m_mu = m_f->m_session->GetParameter("Kinvis");
        Vmath::Fill(npoints, m_mu, mu, 1);
        lambda = 0;
    }
    else if(boost::iequals(m_f->m_variables[0], "rho") &&
            boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        std::string m_ViscosityType;
        m_f->m_session->LoadParameter ("mu",     m_mu, 1.78e-05);
        m_f->m_session->LoadParameter ("lambda", lambda, -2.0/3.0);
        m_f->m_session->LoadSolverInfo("ViscosityType", m_ViscosityType
                                      , "Constant");

        if (m_ViscosityType == "Variable")
        {
            // Check equation of state
            std::string eosType;
            bool m_idealGas;
            m_f->m_session->LoadSolverInfo("EquationOfState", eosType,
                "IdealGas");
            m_idealGas = boost::iequals(eosType,"IdealGas");
            ASSERTL0(m_idealGas,
                "Only IdealGas EOS implemented for Variable ViscosityType");

            // Get relevant parameters
            NekDouble  m_gamma;
            NekDouble  m_pInf;
            NekDouble  m_rhoInf;
            NekDouble  m_gasConstant;
            NekDouble  cv_inv;
            m_f->m_session->LoadParameter("Gamma", m_gamma, 1.4);
            m_f->m_session->LoadParameter("pInf", m_pInf, 101325);
            m_f->m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);
            m_f->m_session->LoadParameter("GasConstant", m_gasConstant
                                         , 287.058);

            // Get temperature from flowfield
            cv_inv = (m_gamma - 1.0) / m_gasConstant;
            // e = 1/rho ( E - 1/2 ( rhou^2/rho + ... ) )
            Array<OneD, NekDouble> tmp(npoints, 0.0);
            Array<OneD, NekDouble> energy(npoints, 0.0);
            Array<OneD, NekDouble> temperature(npoints, 0.0);
            Vmath::Vcopy(npoints, exp[m_spacedim+1]->GetPhys(), 1, energy, 1);
            for (int i = 0; i < m_spacedim; i++)
            {
                // rhou^2
                Vmath::Vmul(npoints, exp[i + 1]->GetPhys(), 1
                           , exp[i + 1]->GetPhys(), 1, tmp, 1);
                // rhou^2/rho
                Vmath::Vdiv(npoints, tmp, 1, exp[0]->GetPhys(), 1, tmp, 1);
                // 0.5 rhou^2/rho
                Vmath::Smul(npoints, 0.5, tmp, 1, tmp, 1);
                // E - 0.5 rhou^2/rho - ...
                Vmath::Vsub(npoints, energy, 1, tmp, 1, energy, 1);
            }
            // rhoe/rho
            Vmath::Vdiv(npoints, energy, 1, exp[0]->GetPhys(), 1, energy, 1);
            // T = e/Cv
            Vmath::Smul(npoints, cv_inv, energy, 1, temperature, 1 );

            // Variable viscosity through the Sutherland's law
            //
            // WARNING, if this routine is modified the same must be done in the
            // CompressibleFlowSolver function in VariableConverter.cpp
            // (this class should be restructured).

            const NekDouble C = .38175;
            NekDouble mu_star = m_mu;
            NekDouble T_star  = m_pInf / (m_rhoInf * m_gasConstant);
            NekDouble ratio;
            for (int i = 0; i < npoints; ++i)
            {
                ratio = temperature[i] / T_star;
                mu[i] = mu_star * ratio * sqrt(ratio) * (1 + C) / (ratio + C);
            }
        }
        else
        {
            Vmath::Fill(npoints, m_mu, mu, 1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Invalid variables for WSS");
    }
}

void ProcessWSS::GetVelocity(
        const Array<OneD, MultiRegions::ExpListSharedPtr> exp,
              Array<OneD, Array<OneD, NekDouble> > &vel)
{
    int npoints = exp[0]->GetNpoints();
    if(boost::iequals(m_f->m_variables[0], "u"))
    {
        // IncNavierStokesSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vcopy(npoints,
                         exp[i]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else if(boost::iequals(m_f->m_variables[0], "rho") &&
            boost::iequals(m_f->m_variables[1], "rhou"))
    {
        // CompressibleFlowSolver
        for (int i = 0; i < m_spacedim; ++i)
        {
            vel[i] = Array<OneD, NekDouble>(npoints);
            Vmath::Vdiv(npoints,
                         exp[i + 1]->GetPhys(), 1,
                         exp[0]->GetPhys(), 1,
                         vel[i], 1);
        }
    }
    else
    {
        // Unknown
        ASSERTL0(false, "Could not identify velocity for WSS");
    }
}

}
}
