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
//  License for the specific language governing rights and limitations under
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
using namespace std;

#include "ProcessWSS.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

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

    NekDouble kinvis = m_f->m_session->GetParameter("Kinvis");

    int i, j;
    int nfields = m_f->m_variables.size();
    int spacedim = m_f->m_graph->GetSpaceDimension();
    if ((m_f->m_numHomogeneousDir) == 1 || (m_f->m_numHomogeneousDir) == 2)
    {
        spacedim += m_f->m_numHomogeneousDir;
    }

    if (spacedim == 2)
    {
        m_f->m_variables.push_back("Shear_x");
        m_f->m_variables.push_back("Shear_y");
        m_f->m_variables.push_back("Shear_mag");
    }
    else
    {
        m_f->m_variables.push_back("Shear_x");
        m_f->m_variables.push_back("Shear_y");
        m_f->m_variables.push_back("Shear_z");
        m_f->m_variables.push_back("Shear_mag");
    }
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    ASSERTL0(m_f->m_variables[0] == "u",
             "Implicit assumption that input is in incompressible format of "
             "(u,v,p) or (u,v,w,p)");

    if (spacedim == 1)
    {
        ASSERTL0(false, "Error: wss for a 1D problem cannot "
                        "be computed");
    }

    int newfields = spacedim + 1;
    int nshear    = spacedim + 1;
    int nstress   = spacedim * spacedim;
    int ngrad     = spacedim * spacedim;

    Array<OneD, Array<OneD, NekDouble> > velocity(nfields), grad(ngrad),
        fgrad(ngrad);
    Array<OneD, Array<OneD, NekDouble> > stress(nstress), fstress(nstress);
    Array<OneD, Array<OneD, NekDouble> > fshear(nshear);

    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp(newfields);
    Array<OneD, MultiRegions::ExpListSharedPtr> BndElmtExp(spacedim);

    m_f->m_exp.resize(nfields + newfields);
    for (i = 0; i < newfields; ++i)
    {
        m_f->m_exp[nfields + i] =
            m_f->AppendExpList(m_f->m_numHomogeneousDir);
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
            for (i = 0; i < newfields; i++)
            {
                BndExp[i] =
                    m_f->m_exp[nfields + i]->UpdateBndCondExpansion(bnd);
            }
            for (i = 0; i < spacedim; i++)
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

            // initialise arrays in the boundary
            for (i = 0; i < nstress; ++i)
            {
                fstress[i] = Array<OneD, NekDouble>(nqb);
            }

            for (i = 0; i < ngrad; ++i)
            {
                fgrad[i] = Array<OneD, NekDouble>(nqb);
            }

            for (i = 0; i < nshear; ++i)
            {
                fshear[i] = Array<OneD, NekDouble>(nqb, 0.0);
            }

            // Extract Velocities
            for (i = 0; i < spacedim; ++i)
            {
                velocity[i] = BndElmtExp[i]->GetPhys();
            }

            // Compute gradients (velocity correction scheme method)
            for (i = 0; i < spacedim; ++i)
            {
                if (spacedim == 2)
                {
                    BndElmtExp[i]->PhysDeriv(velocity[i],
                                             grad[i * spacedim + 0],
                                             grad[i * spacedim + 1]);
                }
                else
                {
                    BndElmtExp[i]->PhysDeriv(
                        velocity[i], grad[i * spacedim + 0],
                        grad[i * spacedim + 1], grad[i * spacedim + 2]);
                }
            }

            // Compute stress component terms  tau_ij = mu*(u_i,j + u_j,i)
            for (i = 0; i < spacedim; ++i)
            {
                for (j = 0; j < spacedim; ++j)
                {
                    Vmath::Vadd(nqe, grad[i * spacedim + j], 1,
                                grad[j * spacedim + i], 1,
                                stress[i * spacedim + j], 1);

                    Vmath::Smul(nqe, kinvis, stress[i * spacedim + j], 1,
                                stress[i * spacedim + j], 1);
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
            for (i = 0; i < spacedim; ++i)
            {
                Vmath::Neg(nqb, normals[i], 1);
            }

            // calculate wss, and update coeffs in the boundary expansion
            // S = tau_ij * n_j
            for (i = 0; i < spacedim; ++i)
            {
                for (j = 0; j < spacedim; ++j)
                {
                    Vmath::Vvtvp(nqb, normals[j], 1, fstress[i * spacedim + j],
                                 1, fshear[i], 1, fshear[i], 1);
                }
            }

            // T = S - (S.n)n
            for (i = 0; i < spacedim; ++i)
            {
                Vmath::Vvtvp(nqb, normals[i], 1, fshear[i], 1,
                             fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            }
            Vmath::Smul(nqb, -1.0, fshear[nshear - 1], 1, fshear[nshear - 1],
                        1);

            for (i = 0; i < spacedim; i++)
            {
                Vmath::Vvtvp(nqb, normals[i], 1, fshear[nshear - 1], 1,
                             fshear[i], 1, fshear[i], 1);
                BndExp[i]->FwdTrans_IterPerExp(fshear[i],
                                        BndExp[i]->UpdateCoeffs());
            }

            // Tw
            Vmath::Zero(nqb, fshear[nshear - 1], 1);
            for (i = 0; i < spacedim; ++i)
            {
                Vmath::Vvtvp(nqb, fshear[i], 1, fshear[i], 1,
                             fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            }
            Vmath::Vsqrt(nqb, fshear[nshear - 1], 1, fshear[nshear - 1], 1);
            BndExp[nshear - 1]->FwdTrans_IterPerExp(fshear[nshear - 1],
                                         BndExp[nshear - 1]->UpdateCoeffs());
        }
    }
}
}
}
