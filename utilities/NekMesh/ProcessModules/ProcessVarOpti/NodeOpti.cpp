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

#include <boost/multi_array.hpp>

#include <limits>

#include "NodeOpti.h"
#include "Evaluator.hxx"
#include "Hessian.hxx"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

NodeOptiFactory &GetNodeOptiFactory()
{
    static NodeOptiFactory asd;
    return asd;
}

void NodeOpti::CalcMinJac()
{
    m_minJac = numeric_limits<double>::max();
    for (auto &typeIt : m_data)
    {
        for (int i = 0; i < typeIt.second.size(); i++)
        {
            m_minJac = min(m_minJac, typeIt.second[i]->GetMinJac());
        }
    }
}

int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<2>(minJacNew);
    NekDouble newVal   = currentW;

    // Gradient already zero
    if (m_grad[0] * m_grad[0] + m_grad[1] * m_grad[1] > gradTol())
    {
        // needs to optimise
        NekDouble xc = m_node->m_x;
        NekDouble yc = m_node->m_y;

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
            // Update node
            m_node->m_x = xc + alpha * sk[0];
            m_node->m_y = yc + alpha * sk[1];

            newVal = GetFunctional<2>(minJacNew, false);

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
            m_node->m_x = xc;
            m_node->m_y = yc;

            mtx.lock();
            m_res->nReset[2]++;
            mtx.unlock();
        }
        else
        {
            m_minJac = minJacNew;

            mtx.lock();
            if (alpha < 1.0)
            {
                m_res->alphaI++;
            }
            mtx.unlock();
        }

        mtx.lock();
        m_res->val = max(sqrt((m_node->m_x - xc) * (m_node->m_x - xc) +
                              (m_node->m_y - yc) * (m_node->m_y - yc)),
                         m_res->val);
        mtx.unlock();
    }

    mtx.lock();
    m_res->func += newVal;
    mtx.unlock();
}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");

void NodeOpti3D3D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<3>(minJacNew);
    NekDouble newVal   = currentW;

    if (m_grad[0] * m_grad[0] + m_grad[1] * m_grad[1] + m_grad[2] * m_grad[2] >
        gradTol())
    {
        // needs to optimise
        NekDouble xc = m_node->m_x;
        NekDouble yc = m_node->m_y;
        NekDouble zc = m_node->m_z;

        vector<NekDouble> sk(3);
        NekDouble val;

        // Calculate minimum eigenvalue
        MinEigen<3>(val);

        if (val < 1e-6)
        {
            // Add constant identity to Hessian matrix.
            m_grad[3] += 1e-6 - val;
            m_grad[6] += 1e-6 - val;
            m_grad[8] += 1e-6 - val;
        }

        // calculate sk
        NekDouble det =
            m_grad[3] * (m_grad[6] * m_grad[8] - m_grad[7] * m_grad[7]) -
            m_grad[4] * (m_grad[4] * m_grad[8] - m_grad[5] * m_grad[7]) +
            m_grad[5] * (m_grad[4] * m_grad[7] - m_grad[5] * m_grad[6]);

        sk[0] = m_grad[0] * (m_grad[6] * m_grad[8] - m_grad[7] * m_grad[7]) +
                m_grad[1] * (m_grad[5] * m_grad[7] - m_grad[4] * m_grad[8]) +
                m_grad[2] * (m_grad[4] * m_grad[7] - m_grad[3] * m_grad[7]);
        sk[1] = m_grad[0] * (m_grad[7] * m_grad[5] - m_grad[4] * m_grad[5]) +
                m_grad[1] * (m_grad[3] * m_grad[8] - m_grad[5] * m_grad[5]) +
                m_grad[2] * (m_grad[4] * m_grad[5] - m_grad[3] * m_grad[7]);
        sk[2] = m_grad[0] * (m_grad[4] * m_grad[7] - m_grad[6] * m_grad[5]) +
                m_grad[1] * (m_grad[4] * m_grad[5] - m_grad[3] * m_grad[7]) +
                m_grad[2] * (m_grad[3] * m_grad[6] - m_grad[4] * m_grad[4]);

        sk[0] /= det * -1.0;
        sk[1] /= det * -1.0;
        sk[2] /= det * -1.0;

        bool found = false;

        NekDouble pg =
            (m_grad[0] * sk[0] + m_grad[1] * sk[1] + m_grad[2] * sk[2]);

        // normal gradient line Search
        NekDouble alpha = 1.0;

        while (alpha > alphaTol())
        {
            // Update node
            m_node->m_x = xc + alpha * sk[0];
            m_node->m_y = yc + alpha * sk[1];
            m_node->m_z = zc + alpha * sk[2];

            newVal = GetFunctional<3>(minJacNew, false);
            // dont need the hessian again this function updates G to be the new
            // location
            //
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
            m_node->m_x = xc;
            m_node->m_y = yc;
            m_node->m_z = zc;

            mtx.lock();
            m_res->nReset[2]++;
            mtx.unlock();
        }
        else
        {
            m_minJac = minJacNew;

            mtx.lock();
            if (alpha < 1.0)
            {
                m_res->alphaI++;
            }
            mtx.unlock();
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

NodeOptiJob *NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}
}
}
