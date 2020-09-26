////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessL2Criterion.cpp
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
//  Description: Computes Lambda 2 Criterion field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessL2Criterion.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessL2Criterion::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "L2Criterion"), ProcessL2Criterion::create,
        "Computes Lambda 2 Criterion.");

ProcessL2Criterion::ProcessL2Criterion(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessL2Criterion::~ProcessL2Criterion()
{
}

/**
 * @brief Calculates eigenvalues of a 3x3 Symmetric matrix.
 *
 * @param d1, d2, d3 - matrix diagonal entries at [0,0], [1,1] and [2,2]
 * @param a - matrix value at [0,1] and [1,0]
 * @param b - matrix value at [0,2] and [2,0]
 * @param c - matrix value at [1,2] and [2,1]
 * @param l1, l2, l3 the computed eigenvalues, ordered l3 >= l2 >= l1
 */
void MatSymEVals(NekDouble d1, NekDouble d2, NekDouble d3, NekDouble a,
                 NekDouble b, NekDouble c, NekDouble &l1, NekDouble &l2,
                 NekDouble &l3)
{
    NekDouble p = a * a + b * b + c * c;
    if (p == 0)
    {
        l1 = d1;
        l2 = d2;
        l3 = d3;
        if (l1 > l3)
        {
            swap(l1, l3);
        }
        if (l1 > l2)
        {
            swap(l1, l2);
        }
        if (l2 > l3)
        {
            swap(l2, l3);
        }
    }
    else
    {
        NekDouble q = (d1 + d2 + d3) / 3.0;
        p = (d1 - q) * (d1 - q) + (d2 - q) * (d2 - q) + (d3 - q) * (d3 - q) +
            2.0 * p;
        p = sqrt(p / 6.0);
        NekDouble r =
            -0.5 *
            (a * a * d3 - a * a * q - 2.0 * a * b * c + b * b * d2 - b * b * q +
             c * c * d1 - c * c * q - d1 * d2 * d3 + d1 * d2 * q + d1 * d3 * q -
             d1 * q * q + d2 * d3 * q - d2 * q * q - d3 * q * q + q * q * q) /
            (p * p * p);

        NekDouble phi = 0;
        if (r <= -1)
        {
            phi = M_PI / 3.0;
        }
        else if (r >= 1)
        {
            phi = 0.0;
        }
        else
        {
            phi = acos(r) / 3.0;
        }

        // the eigenvalues satisfy eig3 >= eig2 >= eig1
        l3 = q + 2.0 * p * cos(phi);
        l1 = q + 2.0 * p * cos(phi + (2.0 * M_PI / 3.0));
        // since trace(A) = eig1 + eig2 + eig3
        l2 = 3.0 * q - l1 - l3;
    }
}

void ProcessL2Criterion::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    auto nfields = m_f->m_variables.size();
    m_f->m_variables.push_back("L2");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int i, s;
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim + (m_f->m_numHomogeneousDir);

    ASSERTL0(
        spacedim == 3,
        "ProcessL2Criterion must be computed for a 3D (or quasi-3D) case.");

    int npoints = m_f->m_exp[0]->GetNpoints();

    Array<OneD, Array<OneD, NekDouble>> grad(spacedim * spacedim);

    // Will store the Lambdas
    NekDouble a00, a11, a22, a01, a02, a12;
    NekDouble t1, t2, t3, t4, t5, t6, t7, t8, t10, t11, t13, t14, t15;
    NekDouble outfield1, outfield3;
    Array<OneD, NekDouble> outfield2(npoints);

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

        /*
         * For each node calculate the S^2+W^2 tensor
         * where S and W are the symmetric and the skew-symmetric
         * parts of the velocity gradient tensor D=grad(u),
         * S=0.5(D+transpose(D)) and W=0.5((D-transpose(D)))
         */
        for (int j = 0; j < npoints; ++j)
        {
            // diff(u,y) + diff(v,x);
            t1 = grad[0 * spacedim + 1][j] + grad[1 * spacedim + 0][j];
            // diff(u,z) + diff(w,x);
            t2 = grad[0 * spacedim + 2][j] + grad[2 * spacedim + 0][j];
            // diff(u,y) - diff(v,x);
            t3 = grad[0 * spacedim + 1][j] - grad[1 * spacedim + 0][j];
            // diff(u,z) - diff(w,x);
            t4 = grad[0 * spacedim + 2][j] - grad[2 * spacedim + 0][j];

            t5 = t2 * t2;
            t6 = t4 * t4;
            t7 = t3 * t3;
            t8 = t1 * t1;

            // diff(w,y) + diff(v,z);
            t10 = grad[2 * spacedim + 1][j] + grad[1 * spacedim + 2][j];
            // diff(w,y) - diff(v,z);
            t11 = grad[2 * spacedim + 1][j] - grad[1 * spacedim + 2][j];

            t13 = 0.25 * (t10 * t2 + t11 * t4) +
                  0.5 * t1 *
                      (grad[0 * spacedim + 0][j] + grad[1 * spacedim + 1][j]);
            t14 = 0.5 * t2 *
                      (grad[0 * spacedim + 0][j] + grad[2 * spacedim + 2][j]) +
                  0.25 * (t1 * t10 - t11 * t3);
            t15 = t10 * t10;
            t11 = t11 * t11;
            t1  = 0.5 * t10 *
                     (grad[1 * spacedim + 1][j] + grad[2 * spacedim + 2][j]) -
                 0.25 * (-t1 * t2 + t3 * t4);

            a00 = 0.25 * (t5 - t6 - t7 + t8) +
                  grad[0 * spacedim + 0][j] * grad[0 * spacedim + 0][j];
            a01 = t13;
            a02 = t14;
            a11 = 0.25 * (-t7 + t8 + t15 - t11) +
                  grad[1 * spacedim + 1][j] * grad[1 * spacedim + 1][j];
            a12 = t1;
            a22 = 0.25 * (t5 - t6 + t15 - t11) +
                  grad[2 * spacedim + 2][j] * grad[2 * spacedim + 2][j];

            // Compute the eigenvalues of a symmetric 3x3 matrix
            MatSymEVals(a00, a11, a22, a01, a02, a12, outfield1, outfield2[j],
                        outfield3);
        }

        Exp = m_f->AppendExpList(m_f->m_numHomogeneousDir);
        Vmath::Vcopy(npoints, outfield2, 1, Exp->UpdatePhys(), 1);
        Exp->FwdTrans_IterPerExp(outfield2, Exp->UpdateCoeffs());
        auto it = m_f->m_exp.begin() + s * (nfields + 1) + nfields;
        m_f->m_exp.insert(it, Exp);
    }
}
}
}
