////////////////////////////////////////////////////////////////////////////////
//
//  File: NodeOpti.cpp
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "NodeOptiCAD.h"
#include "Evaluator.hxx"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

std::mutex mtx;

int NodeOpti1D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    13, NodeOpti1D3D::create, "1D3D");

void NodeOpti1D3D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<3>(minJacNew);
    NekDouble newVal   = currentW;

    if (m_grad[0] * m_grad[0] + m_grad[1] * m_grad[1] + m_grad[2] * m_grad[2] >
        gradTol())
    {
        // modify the gradient to be on the cad system
        ProcessGradient();

        // needs to optimise
        NekDouble tc = m_node->GetCADCurveInfo(curve->GetId());
        NekDouble xc = m_node->m_x;
        NekDouble yc = m_node->m_y;
        NekDouble zc = m_node->m_z;
        NekDouble nt;

        vector<NekDouble> sk(1);

        if (m_grad[1] < 1e-6)
        {
            m_grad[1] = 1e-6 - m_grad[1];
        }

        sk[0] = m_grad[0] / m_grad[1] * -1.0;

        vector<NekDouble> bd(2);
        curve->GetBounds(bd[0], bd[1]);

        bool found = false;

        NekDouble pg = m_grad[0] * sk[0];

        // normal gradient line Search
        NekDouble alpha = 1.0;

        while (alpha > alphaTol())
        {
            // Update node
            nt = tc + alpha * sk[0];
            if (nt < bd[0] || nt > bd[1])
            {
                alpha /= 2.0;
                continue;
            }

            curve->P(nt, m_node->m_x, m_node->m_y, m_node->m_z);

            newVal = GetFunctional<3>(minJacNew, false);

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * pg)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if (!found)
        {
            // reset the node
            nt = tc;
            curve->P(nt, m_node->m_x, m_node->m_y, m_node->m_z);

            mtx.lock();
            m_res->nReset[0]++;
            mtx.unlock();
        }
        else
        {
            m_minJac = minJacNew;
            m_node->Move(m_node->m_x, m_node->m_y, m_node->m_z, curve->GetId(),
                         nt);
        }
        mtx.lock();
        m_res->val = max(sqrt((m_node->m_x - xc) * (m_node->m_x - xc) +
                              (m_node->m_y - yc) * (m_node->m_y - yc) +
                              (m_node->m_z - zc) * (m_node->m_z - zc)),
                         m_res->val);
        mtx.unlock();
    }

    mtx.lock();
    m_res->func += newVal;
    mtx.unlock();
}

int NodeOpti2D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    23, NodeOpti2D3D::create, "2D3D");

void NodeOpti2D3D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<3>(minJacNew);
    NekDouble newVal   = currentW;

    if (m_grad[0] * m_grad[0] + m_grad[1] * m_grad[1] + m_grad[2] * m_grad[2] >
        gradTol())
    {
        // modify the gradient to be on the cad system
        ProcessGradient();

        // needs to optimise
        Array<OneD, NekDouble> uvc = m_node->GetCADSurfInfo(surf->GetId());
        NekDouble xc = m_node->m_x;
        NekDouble yc = m_node->m_y;
        NekDouble zc = m_node->m_z;
        Array<OneD, NekDouble> uvt(2);
        vector<NekDouble> bd(4);
        surf->GetBounds(bd[0], bd[1], bd[2], bd[3]);

        vector<NekDouble> sk(2);
        NekDouble val;

        // Calculate minimum eigenvalue
        MinEigen<2>(val);

        if (val < 1e-6)
        {
            // Add constant identity to Hessian matrix.
            m_grad[2] += 1e-6 - val;
            m_grad[4] += 1e-6 - val;
        }

        sk[0] = -1.0 / (m_grad[2] * m_grad[4] - m_grad[3] * m_grad[3]) *
                (m_grad[4] * m_grad[0] - m_grad[3] * m_grad[1]);
        sk[1] = -1.0 / (m_grad[2] * m_grad[4] - m_grad[3] * m_grad[3]) *
                (m_grad[2] * m_grad[1] - m_grad[3] * m_grad[0]);

        bool found   = false;
        NekDouble pg = (m_grad[0] * sk[0] + m_grad[1] * sk[1]);

        // normal gradient line Search
        NekDouble alpha = 1.0;
        while (alpha > alphaTol())
        {
            uvt[0] = uvc[0] + alpha * sk[0];
            uvt[1] = uvc[1] + alpha * sk[1];

            if (uvt[0] < bd[0] || uvt[0] > bd[1] || uvt[1] < bd[2] ||
                uvt[1] > bd[3])
            {
                alpha /= 2.0;
                continue;
            }

            surf->P(uvt, m_node->m_x, m_node->m_y, m_node->m_z);

            newVal = GetFunctional<3>(minJacNew, false);

            // Wolfe conditions
            if (newVal <= currentW + c1() * (alpha * pg))
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if (!found)
        {
            // reset the node
            surf->P(uvc, m_node->m_x, m_node->m_y, m_node->m_z);

            mtx.lock();
            m_res->nReset[1]++;
            mtx.unlock();
        }
        else
        {
            m_minJac = minJacNew;
            m_node->Move(m_node->m_x, m_node->m_y, m_node->m_z, surf->GetId(),
                         uvt);
        }

        mtx.lock();
        m_res->val = max(sqrt((m_node->m_x - xc) * (m_node->m_x - xc) +
                              (m_node->m_y - yc) * (m_node->m_y - yc) +
                              (m_node->m_z - zc) * (m_node->m_z - zc)),
                         m_res->val);
        mtx.unlock();

        Array<OneD, NekDouble> uva = m_node->GetCADSurfInfo(surf->GetId());
        if (uva[0] < bd[0] || uva[0] > bd[1] || uva[1] < bd[2] ||
            uva[1] > bd[3])
        {
            ASSERTL0(false, "something when very wrong and the node finished "
                            "out of its parameter plane");
        }
    }

    mtx.lock();
    m_res->func += newVal;
    mtx.unlock();
}

int NodeOpti1D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    12, NodeOpti1D2D::create, "1D2D");

void NodeOpti1D2D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<2>(minJacNew);
    NekDouble newVal   = currentW;

    if (m_grad[0] * m_grad[0] + m_grad[1] * m_grad[1] > gradTol())
    {
        // modify the gradient to be on the cad system
        ProcessGradient();

        // needs to optimise
        NekDouble tc = m_node->GetCADCurveInfo(curve->GetId());
        NekDouble xc = m_node->m_x;
        NekDouble yc = m_node->m_y;
        NekDouble zc = m_node->m_z;
        NekDouble nt;
        Array<OneD, NekDouble> p;

        Array<OneD, NekDouble> sk(1);

        if (m_grad[1] < 1e-6)
        {
            m_grad[1] = 1e-6 - m_grad[1];
        }

        sk[0] = m_grad[0] / m_grad[1] * -1.0;

        bool found = false;

        NekDouble pg = m_grad[0] * sk[0];

        // normal gradient line Search
        NekDouble alpha = 1.0;

        while (alpha > alphaTol())
        {
            // Update node
            nt = tc + alpha * sk[0];
            if (nt < m_bd[0] || nt > m_bd[1])
            {
                alpha /= 2.0;
                continue;
            }

            p           = curve->P(nt);
            m_node->m_x = p[0];
            m_node->m_y = p[1];
            m_node->m_z = p[2];

            newVal = GetFunctional<2>(minJacNew, false);

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * pg)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if (!found)
        {
            // reset the node
            nt          = tc;
            p           = curve->P(nt);
            m_node->m_x = p[0];
            m_node->m_y = p[1];
            m_node->m_z = p[2];

            mtx.lock();
            m_res->nReset[0]++;
            mtx.unlock();
        }
        else
        {
            m_minJac = minJacNew;
            m_node->Move(p, curve->GetId(), nt);
        }

        mtx.lock();
        m_res->val = max(sqrt((m_node->m_x - xc) * (m_node->m_x - xc) +
                              (m_node->m_y - yc) * (m_node->m_y - yc) +
                              (m_node->m_z - zc) * (m_node->m_z - zc)),
                         m_res->val);
        mtx.unlock();
    }

    mtx.lock();
    m_res->func += newVal;
    mtx.unlock();
}

void NodeOpti1D3D::ProcessGradient()
{
    NekDouble tc           = m_node->GetCADCurveInfo(curve->GetId());
    vector<NekDouble> grad = m_grad;
    m_grad                 = vector<NekDouble>(2, 0.0);

    // Grab first and second order CAD derivatives
    Array<OneD, NekDouble> d2 = curve->D2(tc);

    // Multiply gradient by derivative of CAD
    m_grad[0] = grad[0] * d2[3] + grad[1] * d2[4] + grad[2] * d2[5];

    // Second order: product rule of above, so multiply gradient by second
    // order CAD derivatives and Hessian by gradient of CAD
    m_grad[1] = grad[0] * d2[6] + grad[1] * d2[7] + grad[2] * d2[8] +
                d2[3] * (grad[3] * d2[3] + grad[4] * d2[4] + grad[5] * d2[5]) +
                d2[4] * (grad[4] * d2[3] + grad[6] * d2[4] + grad[7] * d2[5]) +
                d2[5] * (grad[5] * d2[3] + grad[7] * d2[4] + grad[8] * d2[5]);
}

void NodeOpti2D3D::ProcessGradient()
{
    Array<OneD, NekDouble> uvc = m_node->GetCADSurfInfo(surf->GetId());

    vector<NekDouble> grad = m_grad;
    m_grad                 = vector<NekDouble>(5, 0.0);

    Array<OneD, NekDouble> d2 = surf->D2(uvc);
    // r[0]   x
    // r[1]   y
    // r[2]   z
    // r[3]   dx/du a
    // r[4]   dy/du b
    // r[5]   dz/du c
    // r[6]   dx/dv d
    // r[7]   dy/dv e
    // r[8]   dz/dv f
    // r[9]   d2x/du2
    // r[10]  d2y/du2
    // r[11]  d2z/du2
    // r[12]  d2x/dv2
    // r[13]  d2y/dv2
    // r[14]  d2z/dv2
    // r[15]  d2x/dudv
    // r[16]  d2y/dudv
    // r[17]  d2z/dudv
    //
    // grad[0] d/dx
    // grad[1] d/dy
    // grad[2] d/dz
    // grad[3] d2/dx2
    // grad[4] d2/dxdy
    // grad[5] d2/dxdz
    // grad[6] d2/dy2
    // grad[7] d2/dydz
    // grad[8] d2/dz2

    // Gradients
    m_grad[0] = d2[3] * grad[0] + d2[4] * grad[1] + d2[5] * grad[2];
    m_grad[1] = d2[6] * grad[0] + d2[7] * grad[1] + d2[8] * grad[2];

    m_grad[2] = grad[0] * d2[9] + grad[1] * d2[10] + grad[2] * d2[11] +
                grad[3] * d2[3] * d2[3] + grad[6] * d2[4] * d2[4] +
                grad[8] * d2[5] * d2[5] + 2.0 * grad[4] * d2[3] * d2[4] +
                2.0 * grad[5] * d2[3] * d2[5] + 2.0 * grad[7] * d2[4] * d2[5];

    m_grad[4] = grad[0] * d2[12] + grad[1] * d2[13] + grad[2] * d2[14] +
                grad[3] * d2[6] * d2[6] + grad[6] * d2[7] * d2[7] +
                grad[8] * d2[8] * d2[8] + 2.0 * grad[4] * d2[6] * d2[7] +
                2.0 * grad[5] * d2[6] * d2[8] + 2.0 * grad[7] * d2[7] * d2[8];

    m_grad[3] = grad[0] * d2[15] + grad[1] * d2[16] + grad[2] * d2[17] +
                grad[3] * d2[3] * d2[6] + grad[6] * d2[4] * d2[7] +
                grad[8] * d2[5] * d2[8] +
                grad[4] * (d2[3] * d2[7] + d2[4] * d2[6]) +
                grad[5] * (d2[3] * d2[8] + d2[5] * d2[6]) +
                grad[7] * (d2[4] * d2[8] + d2[5] * d2[7]);
}

void NodeOpti1D2D::ProcessGradient()
{
    NekDouble tc = m_node->GetCADCurveInfo(curve->GetId());
    vector<NekDouble> grad = m_grad;
    m_grad = vector<NekDouble>(2, 0.0);

    // Grab first and second order CAD derivatives
    Array<OneD, NekDouble> d2 = curve->D2(tc);

    // Multiply gradient by derivative of CAD
    m_grad[0] = grad[0] * d2[3] + grad[1] * d2[4];

    // Second order: product rule of above, so multiply gradient by second
    // order CAD derivatives and Hessian by gradient of CAD
    m_grad[1] = grad[0] * d2[6] + grad[1] * d2[7] +
                d2[3] * (grad[2] * d2[3] + grad[3] * d2[4]) +
                d2[4] * (grad[3] * d2[3] + grad[4] * d2[4]);
}
}
}
