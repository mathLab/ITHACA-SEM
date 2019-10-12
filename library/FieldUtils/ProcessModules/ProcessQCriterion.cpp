////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessQCriterion.cpp
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
//  Description: Computes Q Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessQCriterion.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessQCriterion::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "QCriterion"),
        ProcessQCriterion::create,
        "Computes Q-Criterion.");

ProcessQCriterion::ProcessQCriterion(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessQCriterion::~ProcessQCriterion()
{
}

void ProcessQCriterion::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nfields = m_f->m_variables.size();
    m_f->m_variables.push_back("Q");
    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int i, s;
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim + (m_f->m_numHomogeneousDir);

    ASSERTL0(spacedim == 3,
        "ProcessQCriterion must be computed for a 3D (or quasi-3D) case.");

    int npoints = m_f->m_exp[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble> > grad(spacedim * spacedim);

    Array<OneD, NekDouble> omega(npoints);
    Array<OneD, NekDouble> S(npoints);

    // Will store the Q-Criterion
    Array<OneD, NekDouble> outfield (npoints);
    Array<OneD, NekDouble> outfield1(npoints);
    Array<OneD, NekDouble> outfield2(npoints);
    Array<OneD, NekDouble> outfield3(npoints);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    for (i = 0; i < spacedim * spacedim; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    MultiRegions::ExpListSharedPtr Exp;

    for (s = 0; s < nstrips; ++s) // homogeneous strip varient
    {
        for (i = 0; i < spacedim; ++i)
        {
            m_f->m_exp[s * nfields + i]->PhysDeriv(
                m_f->m_exp[s * nfields + i]->GetPhys(), grad[i * spacedim],
                grad[i * spacedim + 1], grad[i * spacedim + 2]);
        }

        // W_x = Wy - Vz
        Vmath::Vsub(npoints, grad[2 * spacedim + 1], 1,
                             grad[1 * spacedim + 2], 1,
                             outfield1, 1);
        // W_x^2
        Vmath::Vmul(npoints, outfield1, 1, outfield1, 1, outfield1, 1);

        // W_y = Uz - Wx
        Vmath::Vsub(npoints, grad[0 * spacedim + 2], 1,
                             grad[2 * spacedim + 0], 1,
                             outfield2, 1);
        // W_y^2
        Vmath::Vmul(npoints, outfield2, 1, outfield2, 1, outfield2, 1);

        // W_z = Vx - Uy
        Vmath::Vsub(npoints, grad[1 * spacedim + 0], 1,
                             grad[0 * spacedim + 1], 1,
                             outfield3, 1);
        // W_z^2
        Vmath::Vmul(npoints, outfield3, 1, outfield3, 1, outfield3, 1);

        // Omega = 0.5*(W_x^2 + W_y^2 + W_z^2)
        NekDouble fac = 0.5;
        Vmath::Vadd(npoints, outfield1, 1, outfield2, 1, omega, 1);
        Vmath::Vadd(npoints, omega,     1, outfield3, 1, omega, 1);
        Vmath::Smul(npoints, fac, omega, 1, omega, 1);

        // Ux^2
        Vmath::Vmul(npoints, grad[0 * spacedim + 0], 1,
                             grad[0 * spacedim + 0], 1,
                             outfield1, 1);
        // Vy^2
        Vmath::Vmul(npoints, grad[1 * spacedim + 1], 1,
                             grad[1 * spacedim + 1], 1,
                             outfield2, 1);
        // Wz^2
        Vmath::Vmul(npoints, grad[2 * spacedim + 2], 1,
                             grad[2 * spacedim + 2], 1,
                             outfield3, 1);

        //
        Vmath::Vadd(npoints, outfield1, 1, outfield2, 1, S, 1);
        Vmath::Vadd(npoints, S,         1, outfield3, 1, S, 1);

        // Wy + Vz
        Vmath::Vadd(npoints, grad[2 * spacedim + 1], 1,
                             grad[1 * spacedim + 2], 1,
                             outfield1, 1);
        Vmath::Vmul(npoints, outfield1, 1, outfield1, 1, outfield1, 1);

        // Uz + Wx
        Vmath::Vadd(npoints, grad[0 * spacedim + 2], 1,
                             grad[2 * spacedim + 0], 1,
                             outfield2, 1);
        Vmath::Vmul(npoints, outfield2, 1, outfield2, 1, outfield2, 1);

        // Vx + Uy
        Vmath::Vadd(npoints, grad[1 * spacedim + 0], 1,
                             grad[0 * spacedim + 1], 1,
                             outfield3, 1);
        Vmath::Vmul(npoints, outfield3, 1, outfield3, 1, outfield3, 1);

        Vmath::Vadd(npoints, outfield1, 1, outfield2, 1, outfield2, 1);
        Vmath::Vadd(npoints, outfield2, 1, outfield3, 1, outfield3, 1);

        Vmath::Smul(npoints, fac, outfield3, 1, outfield3, 1);

        Vmath::Vadd(npoints, outfield3, 1, S, 1, S, 1);
        Vmath::Vsub(npoints, omega, 1, S, 1, outfield, 1);

        Vmath::Smul(npoints, fac, outfield, 1, outfield, 1);

        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, outfield, 1, Exp->UpdatePhys(), 1);
        Exp->FwdTrans_IterPerExp(outfield, Exp->UpdateCoeffs());

        auto it = m_f->m_exp.begin() + s * (nfields + 1) + nfields;
        m_f->m_exp.insert(it, Exp);
    }
}

}
}
