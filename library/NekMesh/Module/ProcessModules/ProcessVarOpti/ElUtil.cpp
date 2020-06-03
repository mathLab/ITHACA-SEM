////////////////////////////////////////////////////////////////////////////////
//
//  File: ElUtil.cpp
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

#include "ElUtil.h"
#include "ProcessVarOpti.h"

#include <mutex>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace Utilities
{

std::mutex mtx2;

ElUtil::ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d, ResidualSharedPtr r,
               int n, int o)
{
    m_el        = e;
    m_derivUtil = d;
    m_res       = r;
    m_mode      = n;
    m_order     = o;
    m_dim       = m_el->GetDim();
    vector<NodeSharedPtr> ns;
    m_el->GetCurvedNodes(ns);
    nodes.resize(ns.size());
    for (int i = 0; i < ns.size(); ++i)
    {
        nodes[i].resize(m_dim);
        nodes[i][0] = &ns[i]->m_x;

        if (m_dim >= 2)
        {
            nodes[i][1] = &ns[i]->m_y;
        }

        if (m_dim >= 3)
        {
            nodes[i][2] = &ns[i]->m_z;
        }

        m_idmap[ns[i]->m_id] = i;
    }
    MappingIdealToRef();
}

void ElUtil::MappingIdealToRef()
{
    if (m_el->GetConf().m_e == LibUtilities::eQuadrilateral)
    {
        LibUtilities::PointsKey pkey1(m_mode, LibUtilities::eNodalQuadElec);
        LibUtilities::PointsKey pkey2(m_mode + m_order,
                                      LibUtilities::eNodalQuadElec);

        Array<OneD, NekDouble> u1(m_derivUtil->ptsStd), v1(m_derivUtil->ptsStd),
            u2(m_derivUtil->pts), v2(m_derivUtil->pts);

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
        LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);

        vector<vector<NekDouble> > xyz(4);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for (int i = 0; i < 4; i++)
        {
            vector<NekDouble> x(3);
            x[0]   = ns[i]->m_x;
            x[1]   = ns[i]->m_y;
            x[2]   = ns[i]->m_z;
            xyz[i] = x;
        }

        for (int i = 0; i < m_derivUtil->ptsStd; ++i)
        {
            NekDouble a1 = 0.5 * (1 - u1[i]);
            NekDouble a2 = 0.5 * (1 + u1[i]);
            NekDouble b1 = 0.5 * (1 - v1[i]);
            NekDouble b2 = 0.5 * (1 + v1[i]);

            DNekMat J(2, 2, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0] +
                      0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1, 0) = -0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1] +
                      0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];

            J(0, 1) = -0.5 * a1 * xyz[0][0] - 0.5 * a2 * xyz[1][0] +
                      0.5 * a2 * xyz[2][0] + 0.5 * a1 * xyz[3][0];
            J(1, 1) = -0.5 * a1 * xyz[0][1] - 0.5 * a2 * xyz[1][1] +
                      0.5 * a2 * xyz[2][1] + 0.5 * a1 * xyz[3][1];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = 0.0;
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            m_mapsStd.push_back(r);
        }

        for (int i = 0; i < m_derivUtil->pts; ++i)
        {
            NekDouble a1 = 0.5 * (1 - u2[i]);
            NekDouble a2 = 0.5 * (1 + u2[i]);
            NekDouble b1 = 0.5 * (1 - v2[i]);
            NekDouble b2 = 0.5 * (1 + v2[i]);

            DNekMat J(2, 2, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0] +
                      0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1, 0) = -0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1] +
                      0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];

            J(0, 1) = -0.5 * a1 * xyz[0][0] - 0.5 * a2 * xyz[1][0] +
                      0.5 * a2 * xyz[2][0] + 0.5 * a1 * xyz[3][0];
            J(1, 1) = -0.5 * a1 * xyz[0][1] - 0.5 * a2 * xyz[1][1] +
                      0.5 * a2 * xyz[2][1] + 0.5 * a1 * xyz[3][1];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = 0.0;
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            m_maps.push_back(r);
        }
    }
    else if (m_el->GetConf().m_e == LibUtilities::eTriangle)
    {
        DNekMat J(2, 2, 0.0);
        J(0, 0) = (*nodes[1][0] - *nodes[0][0]);
        J(1, 0) = (*nodes[1][1] - *nodes[0][1]);
        J(0, 1) = (*nodes[2][0] - *nodes[0][0]);
        J(1, 1) = (*nodes[2][1] - *nodes[0][1]);

        J.Invert();

        DNekMat R(2, 2, 0.0);
        R(0, 0) = 2.0;
        R(1, 1) = 2.0;

        J = J * R;

        for (int i = 0; i < m_derivUtil->pts; i++)
        {
            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0));
            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = 0.0;
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            m_maps.push_back(r);
            m_mapsStd.push_back(r);
        }
    }
    else if (m_el->GetConf().m_e == LibUtilities::eTetrahedron)
    {
        DNekMat J(3, 3, 0.0);
        J(0, 0) = (*nodes[1][0] - *nodes[0][0]);
        J(1, 0) = (*nodes[1][1] - *nodes[0][1]);
        J(2, 0) = (*nodes[1][2] - *nodes[0][2]);
        J(0, 1) = (*nodes[2][0] - *nodes[0][0]);
        J(1, 1) = (*nodes[2][1] - *nodes[0][1]);
        J(2, 1) = (*nodes[2][2] - *nodes[0][2]);
        J(0, 2) = (*nodes[3][0] - *nodes[0][0]);
        J(1, 2) = (*nodes[3][1] - *nodes[0][1]);
        J(2, 2) = (*nodes[3][2] - *nodes[0][2]);

        J.Invert();

        DNekMat R(3, 3, 0.0);
        R(0, 0) = 2.0;
        R(1, 1) = 2.0;
        R(2, 2) = 2.0;

        J = J * R;

        for (int i = 0; i < m_derivUtil->pts; i++)
        {
            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * (J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)) -
                          J(0, 1) * (J(1, 0) * J(2, 2) - J(2, 0) * J(1, 2)) +
                          J(0, 2) * (J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1)));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = J(2, 0);
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = J(2, 1);
            r[6] = J(0, 2);
            r[7] = J(1, 2);
            r[8] = J(2, 2);
            m_maps.push_back(r);
            m_mapsStd.push_back(r);
        }
    }
    else if (m_el->GetConf().m_e == LibUtilities::ePrism)
    {
        LibUtilities::PointsKey pkey1(m_mode, LibUtilities::eNodalPrismElec);
        LibUtilities::PointsKey pkey2(m_mode + m_order,
                                      LibUtilities::eNodalPrismSPI);
        Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
        LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

        vector<vector<NekDouble> > xyz(6);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for (int i = 0; i < 6; i++)
        {
            vector<NekDouble> x(3);
            x[0]   = ns[i]->m_x;
            x[1]   = ns[i]->m_y;
            x[2]   = ns[i]->m_z;
            xyz[i] = x;
        }

        for (int i = 0; i < m_derivUtil->ptsStd; ++i)
        {
            NekDouble a2 = 0.5 * (1 + u1[i]);
            NekDouble b1 = 0.5 * (1 - v1[i]);
            NekDouble b2 = 0.5 * (1 + v1[i]);
            NekDouble c2 = 0.5 * (1 + w1[i]);
            NekDouble d  = 0.5 * (u1[i] + w1[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0] +
                      0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1, 0) = -0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1] +
                      0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];
            J(2, 0) = -0.5 * b1 * xyz[0][2] + 0.5 * b1 * xyz[1][2] +
                      0.5 * b2 * xyz[2][2] - 0.5 * b2 * xyz[3][2];

            J(0, 1) = 0.5 * d * xyz[0][0] - 0.5 * a2 * xyz[1][0] +
                      0.5 * a2 * xyz[2][0] - 0.5 * d * xyz[3][0] -
                      0.5 * c2 * xyz[4][0] + 0.5 * c2 * xyz[5][0];
            J(1, 1) = 0.5 * d * xyz[0][1] - 0.5 * a2 * xyz[1][1] +
                      0.5 * a2 * xyz[2][1] - 0.5 * d * xyz[3][1] -
                      0.5 * c2 * xyz[4][1] + 0.5 * c2 * xyz[5][1];
            J(2, 1) = 0.5 * d * xyz[0][2] - 0.5 * a2 * xyz[1][2] +
                      0.5 * a2 * xyz[2][2] - 0.5 * d * xyz[3][2] -
                      0.5 * c2 * xyz[4][2] + 0.5 * c2 * xyz[5][2];

            J(0, 2) = -0.5 * b1 * xyz[0][0] - 0.5 * b2 * xyz[3][0] +
                      0.5 * b1 * xyz[4][0] + 0.5 * b2 * xyz[5][0];
            J(1, 2) = -0.5 * b1 * xyz[0][1] - 0.5 * b2 * xyz[3][1] +
                      0.5 * b1 * xyz[4][1] + 0.5 * b2 * xyz[5][1];
            J(2, 2) = -0.5 * b1 * xyz[0][2] - 0.5 * b2 * xyz[3][2] +
                      0.5 * b1 * xyz[4][2] + 0.5 * b2 * xyz[5][2];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * (J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)) -
                          J(0, 1) * (J(1, 0) * J(2, 2) - J(2, 0) * J(1, 2)) +
                          J(0, 2) * (J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1)));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = J(2, 0);
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = J(2, 1);
            r[6] = J(0, 2);
            r[7] = J(1, 2);
            r[8] = J(2, 2);
            m_mapsStd.push_back(r);
        }
        for (int i = 0; i < m_derivUtil->pts; ++i)
        {
            NekDouble a2 = 0.5 * (1 + u2[i]);
            NekDouble b1 = 0.5 * (1 - v2[i]);
            NekDouble b2 = 0.5 * (1 + v2[i]);
            NekDouble c2 = 0.5 * (1 + w2[i]);
            NekDouble d  = 0.5 * (u2[i] + w2[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0] +
                      0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1, 0) = -0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1] +
                      0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];
            J(2, 0) = -0.5 * b1 * xyz[0][2] + 0.5 * b1 * xyz[1][2] +
                      0.5 * b2 * xyz[2][2] - 0.5 * b2 * xyz[3][2];

            J(0, 1) = 0.5 * d * xyz[0][0] - 0.5 * a2 * xyz[1][0] +
                      0.5 * a2 * xyz[2][0] - 0.5 * d * xyz[3][0] -
                      0.5 * c2 * xyz[4][0] + 0.5 * c2 * xyz[5][0];
            J(1, 1) = 0.5 * d * xyz[0][1] - 0.5 * a2 * xyz[1][1] +
                      0.5 * a2 * xyz[2][1] - 0.5 * d * xyz[3][1] -
                      0.5 * c2 * xyz[4][1] + 0.5 * c2 * xyz[5][1];
            J(2, 1) = 0.5 * d * xyz[0][2] - 0.5 * a2 * xyz[1][2] +
                      0.5 * a2 * xyz[2][2] - 0.5 * d * xyz[3][2] -
                      0.5 * c2 * xyz[4][2] + 0.5 * c2 * xyz[5][2];

            J(0, 2) = -0.5 * b1 * xyz[0][0] - 0.5 * b2 * xyz[3][0] +
                      0.5 * b1 * xyz[4][0] + 0.5 * b2 * xyz[5][0];
            J(1, 2) = -0.5 * b1 * xyz[0][1] - 0.5 * b2 * xyz[3][1] +
                      0.5 * b1 * xyz[4][1] + 0.5 * b2 * xyz[5][1];
            J(2, 2) = -0.5 * b1 * xyz[0][2] - 0.5 * b2 * xyz[3][2] +
                      0.5 * b1 * xyz[4][2] + 0.5 * b2 * xyz[5][2];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * (J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)) -
                          J(0, 1) * (J(1, 0) * J(2, 2) - J(2, 0) * J(1, 2)) +
                          J(0, 2) * (J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1)));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = J(2, 0);
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = J(2, 1);
            r[6] = J(0, 2);
            r[7] = J(1, 2);
            r[8] = J(2, 2);
            m_maps.push_back(r);
        }
    }
    else if (m_el->GetConf().m_e == LibUtilities::eHexahedron)
    {
        LibUtilities::PointsKey pkey1(m_mode, LibUtilities::eNodalHexElec);
        LibUtilities::PointsKey pkey2(m_mode + m_order,
                                      LibUtilities::eNodalHexElec);
        Array<OneD, NekDouble> u1(m_derivUtil->ptsStd), v1(m_derivUtil->ptsStd),
            w1(m_derivUtil->ptsStd), u2(m_derivUtil->pts), v2(m_derivUtil->pts),
            w2(m_derivUtil->pts);

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
        LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);

        vector<vector<NekDouble> > xyz(8);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for (int i = 0; i < 8; i++)
        {
            vector<NekDouble> x(3);
            x[0]   = ns[i]->m_x;
            x[1]   = ns[i]->m_y;
            x[2]   = ns[i]->m_z;
            xyz[i] = x;
        }

        for (int i = 0; i < m_derivUtil->ptsStd; ++i)
        {
            NekDouble a1 = 0.5 * (1 - u1[i]);
            NekDouble a2 = 0.5 * (1 + u1[i]);
            NekDouble b1 = 0.5 * (1 - v1[i]);
            NekDouble b2 = 0.5 * (1 + v1[i]);
            NekDouble c1 = 0.5 * (1 - w1[i]);
            NekDouble c2 = 0.5 * (1 + w1[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * c1 * xyz[0][0] + 0.5 * b1 * c1 * xyz[1][0] +
                      0.5 * b2 * c1 * xyz[2][0] - 0.5 * b2 * c1 * xyz[3][0] -
                      0.5 * b1 * c2 * xyz[5][0] + 0.5 * b1 * c2 * xyz[5][0] +
                      0.5 * b2 * c2 * xyz[6][0] - 0.5 * b2 * c2 * xyz[7][0];
            J(1, 0) = -0.5 * b1 * c1 * xyz[0][1] + 0.5 * b1 * c1 * xyz[1][1] +
                      0.5 * b2 * c1 * xyz[2][1] - 0.5 * b2 * c1 * xyz[3][1] -
                      0.5 * b1 * c2 * xyz[5][1] + 0.5 * b1 * c2 * xyz[5][1] +
                      0.5 * b2 * c2 * xyz[6][1] - 0.5 * b2 * c2 * xyz[7][1];
            J(2, 0) = -0.5 * b1 * c1 * xyz[0][2] + 0.5 * b1 * c1 * xyz[1][2] +
                      0.5 * b2 * c1 * xyz[2][2] - 0.5 * b2 * c1 * xyz[3][2] -
                      0.5 * b1 * c2 * xyz[5][2] + 0.5 * b1 * c2 * xyz[5][2] +
                      0.5 * b2 * c2 * xyz[6][2] - 0.5 * b2 * c2 * xyz[7][2];

            J(0, 1) = -0.5 * a1 * c1 * xyz[0][0] - 0.5 * a2 * c1 * xyz[1][0] +
                      0.5 * a2 * c1 * xyz[2][0] + 0.5 * a1 * c1 * xyz[3][0] -
                      0.5 * a1 * c2 * xyz[5][0] - 0.5 * a2 * c2 * xyz[5][0] +
                      0.5 * a2 * c2 * xyz[6][0] + 0.5 * a1 * c2 * xyz[7][0];
            J(1, 1) = -0.5 * a1 * c1 * xyz[0][1] - 0.5 * a2 * c1 * xyz[1][1] +
                      0.5 * a2 * c1 * xyz[2][1] + 0.5 * a1 * c1 * xyz[3][1] -
                      0.5 * a1 * c2 * xyz[5][1] - 0.5 * a2 * c2 * xyz[5][1] +
                      0.5 * a2 * c2 * xyz[6][1] + 0.5 * a1 * c2 * xyz[7][1];
            J(2, 1) = -0.5 * a1 * c1 * xyz[0][2] - 0.5 * a2 * c1 * xyz[1][2] +
                      0.5 * a2 * c1 * xyz[2][2] + 0.5 * a1 * c1 * xyz[3][2] -
                      0.5 * a1 * c2 * xyz[5][2] - 0.5 * a2 * c2 * xyz[5][2] +
                      0.5 * a2 * c2 * xyz[6][2] + 0.5 * a1 * c2 * xyz[7][2];

            J(0, 0) = -0.5 * b1 * a1 * xyz[0][0] - 0.5 * b1 * a2 * xyz[1][0] -
                      0.5 * b2 * a2 * xyz[2][0] - 0.5 * b2 * a1 * xyz[3][0] +
                      0.5 * b1 * a1 * xyz[5][0] + 0.5 * b1 * a2 * xyz[5][0] +
                      0.5 * b2 * a2 * xyz[6][0] + 0.5 * b2 * a1 * xyz[7][0];
            J(1, 0) = -0.5 * b1 * a1 * xyz[0][1] - 0.5 * b1 * a2 * xyz[1][1] -
                      0.5 * b2 * a2 * xyz[2][1] - 0.5 * b2 * a1 * xyz[3][1] +
                      0.5 * b1 * a1 * xyz[5][1] + 0.5 * b1 * a2 * xyz[5][1] +
                      0.5 * b2 * a2 * xyz[6][1] + 0.5 * b2 * a1 * xyz[7][1];
            J(2, 0) = -0.5 * b1 * a1 * xyz[0][2] - 0.5 * b1 * a2 * xyz[1][2] -
                      0.5 * b2 * a2 * xyz[2][2] - 0.5 * b2 * a1 * xyz[3][2] +
                      0.5 * b1 * a1 * xyz[5][2] + 0.5 * b1 * a2 * xyz[5][2] +
                      0.5 * b2 * a2 * xyz[6][2] + 0.5 * b2 * a1 * xyz[7][2];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * (J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)) -
                          J(0, 1) * (J(1, 0) * J(2, 2) - J(2, 0) * J(1, 2)) +
                          J(0, 2) * (J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1)));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = J(2, 0);
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = J(2, 1);
            r[6] = J(0, 2);
            r[7] = J(1, 2);
            r[8] = J(2, 2);
            m_mapsStd.push_back(r);
        }

        for (int i = 0; i < m_derivUtil->pts; ++i)
        {
            NekDouble a1 = 0.5 * (1 - u2[i]);
            NekDouble a2 = 0.5 * (1 + u2[i]);
            NekDouble b1 = 0.5 * (1 - v2[i]);
            NekDouble b2 = 0.5 * (1 + v2[i]);
            NekDouble c1 = 0.5 * (1 - w2[i]);
            NekDouble c2 = 0.5 * (1 + w2[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0, 0) = -0.5 * b1 * c1 * xyz[0][0] + 0.5 * b1 * c1 * xyz[1][0] +
                      0.5 * b2 * c1 * xyz[2][0] - 0.5 * b2 * c1 * xyz[3][0] -
                      0.5 * b1 * c2 * xyz[5][0] + 0.5 * b1 * c2 * xyz[5][0] +
                      0.5 * b2 * c2 * xyz[6][0] - 0.5 * b2 * c2 * xyz[7][0];
            J(1, 0) = -0.5 * b1 * c1 * xyz[0][1] + 0.5 * b1 * c1 * xyz[1][1] +
                      0.5 * b2 * c1 * xyz[2][1] - 0.5 * b2 * c1 * xyz[3][1] -
                      0.5 * b1 * c2 * xyz[5][1] + 0.5 * b1 * c2 * xyz[5][1] +
                      0.5 * b2 * c2 * xyz[6][1] - 0.5 * b2 * c2 * xyz[7][1];
            J(2, 0) = -0.5 * b1 * c1 * xyz[0][2] + 0.5 * b1 * c1 * xyz[1][2] +
                      0.5 * b2 * c1 * xyz[2][2] - 0.5 * b2 * c1 * xyz[3][2] -
                      0.5 * b1 * c2 * xyz[5][2] + 0.5 * b1 * c2 * xyz[5][2] +
                      0.5 * b2 * c2 * xyz[6][2] - 0.5 * b2 * c2 * xyz[7][2];

            J(0, 1) = -0.5 * a1 * c1 * xyz[0][0] - 0.5 * a2 * c1 * xyz[1][0] +
                      0.5 * a2 * c1 * xyz[2][0] + 0.5 * a1 * c1 * xyz[3][0] -
                      0.5 * a1 * c2 * xyz[5][0] - 0.5 * a2 * c2 * xyz[5][0] +
                      0.5 * a2 * c2 * xyz[6][0] + 0.5 * a1 * c2 * xyz[7][0];
            J(1, 1) = -0.5 * a1 * c1 * xyz[0][1] - 0.5 * a2 * c1 * xyz[1][1] +
                      0.5 * a2 * c1 * xyz[2][1] + 0.5 * a1 * c1 * xyz[3][1] -
                      0.5 * a1 * c2 * xyz[5][1] - 0.5 * a2 * c2 * xyz[5][1] +
                      0.5 * a2 * c2 * xyz[6][1] + 0.5 * a1 * c2 * xyz[7][1];
            J(2, 1) = -0.5 * a1 * c1 * xyz[0][2] - 0.5 * a2 * c1 * xyz[1][2] +
                      0.5 * a2 * c1 * xyz[2][2] + 0.5 * a1 * c1 * xyz[3][2] -
                      0.5 * a1 * c2 * xyz[5][2] - 0.5 * a2 * c2 * xyz[5][2] +
                      0.5 * a2 * c2 * xyz[6][2] + 0.5 * a1 * c2 * xyz[7][2];

            J(0, 0) = -0.5 * b1 * a1 * xyz[0][0] - 0.5 * b1 * a2 * xyz[1][0] -
                      0.5 * b2 * a2 * xyz[2][0] - 0.5 * b2 * a1 * xyz[3][0] +
                      0.5 * b1 * a1 * xyz[5][0] + 0.5 * b1 * a2 * xyz[5][0] +
                      0.5 * b2 * a2 * xyz[6][0] + 0.5 * b2 * a1 * xyz[7][0];
            J(1, 0) = -0.5 * b1 * a1 * xyz[0][1] - 0.5 * b1 * a2 * xyz[1][1] -
                      0.5 * b2 * a2 * xyz[2][1] - 0.5 * b2 * a1 * xyz[3][1] +
                      0.5 * b1 * a1 * xyz[5][1] + 0.5 * b1 * a2 * xyz[5][1] +
                      0.5 * b2 * a2 * xyz[6][1] + 0.5 * b2 * a1 * xyz[7][1];
            J(2, 0) = -0.5 * b1 * a1 * xyz[0][2] - 0.5 * b1 * a2 * xyz[1][2] -
                      0.5 * b2 * a2 * xyz[2][2] - 0.5 * b2 * a1 * xyz[3][2] +
                      0.5 * b1 * a1 * xyz[5][2] + 0.5 * b1 * a2 * xyz[5][2] +
                      0.5 * b2 * a2 * xyz[6][2] + 0.5 * b2 * a1 * xyz[7][2];

            J.Invert();

            vector<NekDouble> r(10, 0.0); // store det in 10th entry

            r[9] = 1.0 / (J(0, 0) * (J(1, 1) * J(2, 2) - J(2, 1) * J(1, 2)) -
                          J(0, 1) * (J(1, 0) * J(2, 2) - J(2, 0) * J(1, 2)) +
                          J(0, 2) * (J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1)));

            r[0] = J(0, 0);
            r[1] = J(1, 0);
            r[2] = J(2, 0);
            r[3] = J(0, 1);
            r[4] = J(1, 1);
            r[5] = J(2, 1);
            r[6] = J(0, 2);
            r[7] = J(1, 2);
            r[8] = J(2, 2);
            m_maps.push_back(r);
        }
    }
    else
    {
        ASSERTL0(false, "not coded");
    }

    for (int i = 0; i < m_maps.size(); ++i)
    {
        maps.emplace_back(10, 0.0);
        mapsStd.emplace_back(10, 0.0);
    }

    UpdateMapping();
}

void ElUtil::Evaluate()
{
    NekDouble mx2 = -1.0 * numeric_limits<double>::max();
    NekDouble mn2 = numeric_limits<double>::max();

    ASSERTL0(nodes.size() == m_derivUtil->ptsStd, "node count wrong");
    const int nNodes = nodes.size();

    if (m_dim == 2)
    {
        std::vector<NekDouble> X(nNodes), Y(nNodes);
        for (int j = 0; j < nNodes; j++)
        {
            X[j] = *nodes[j][0];
            Y[j] = *nodes[j][1];
        }

        std::vector<NekDouble> x1i(nNodes), y1i(nNodes);
        std::vector<NekDouble> x2i(nNodes), y2i(nNodes);

        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[0].GetRawPtr(),
            nNodes, &X[0], 1, 0.0, &x1i[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[0].GetRawPtr(),
            nNodes, &Y[0], 1, 0.0, &y1i[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[1].GetRawPtr(),
            nNodes, &X[0], 1, 0.0, &x2i[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[1].GetRawPtr(),
            nNodes, &Y[0], 1, 0.0, &y2i[0], 1.0);

        for (int j = 0; j < nNodes; j++)
        {
            NekDouble jacDet = x1i[j] * y2i[j] - x2i[j] * y1i[j];
            mx2              = max(mx2, jacDet / mapsStd[j][9]);
            mn2              = min(mn2, jacDet / mapsStd[j][9]);
        }
    }
    else if (m_dim == 3)
    {
        std::vector<NekDouble> X(nNodes), Y(nNodes), Z(nNodes);
        for (int j = 0; j < nNodes; j++)
        {
            X[j] = *nodes[j][0];
            Y[j] = *nodes[j][1];
            Z[j] = *nodes[j][2];
        }
        std::vector<NekDouble> x1i2(nNodes), y1i2(nNodes), z1i2(nNodes);
        std::vector<NekDouble> x2i2(nNodes), y2i2(nNodes), z2i2(nNodes);
        std::vector<NekDouble> x3i2(nNodes), y3i2(nNodes), z3i2(nNodes);

        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[0].GetRawPtr(),
            nNodes, &X[0], 1, 0.0, &x1i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[0].GetRawPtr(),
            nNodes, &Y[0], 1, 0.0, &y1i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[0].GetRawPtr(),
            nNodes, &Z[0], 1, 0.0, &z1i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[1].GetRawPtr(),
            nNodes, &X[0], 1, 0.0, &x2i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[1].GetRawPtr(),
            nNodes, &Y[0], 1, 0.0, &y2i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[1].GetRawPtr(),
            nNodes, &Z[0], 1, 0.0, &z2i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[2].GetRawPtr(),
            nNodes, &X[0], 1, 0.0, &x3i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[2].GetRawPtr(),
            nNodes, &Y[0], 1, 0.0, &y3i2[0], 1.0);
        Blas::Dgemv(
            'N', nNodes, nNodes, 1.0, m_derivUtil->VdmDStd[2].GetRawPtr(),
            nNodes, &Z[0], 1, 0.0, &z3i2[0], 1.0);

        for (int j = 0; j < nNodes; j++)
        {
            NekDouble jacDet =
                x1i2[j] * (y2i2[j] * z3i2[j] - z2i2[j] * y3i2[j]) -
                x2i2[j] * (y1i2[j] * z3i2[j] - z1i2[j] * y3i2[j]) +
                x3i2[j] * (y1i2[j] * z2i2[j] - z1i2[j] * y2i2[j]);

            mx2 = max(mx2, jacDet / mapsStd[j][9]);
            mn2 = min(mn2, jacDet / mapsStd[j][9]);
        }
    }

    mtx2.lock();
    if (mn2 < 0)
    {
        m_res->startInv++;
    }
    m_res->worstJac = min(m_res->worstJac, (mn2 / mx2));
    mtx2.unlock();

    m_scaledJac = (mn2 / mx2);
}

void ElUtil::InitialMinJac()
{
    NekDouble mx = -1.0 * numeric_limits<double>::max();
    NekDouble mn = numeric_limits<double>::max();

    ASSERTL0(nodes.size() == m_derivUtil->ptsStd, "node count wrong");

    if (m_dim == 2)
    {
        NekVector<NekDouble> X(nodes.size()), Y(nodes.size());
        for (int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
        }

        NekVector<NekDouble> x1i2(m_derivUtil->pts), y1i2(m_derivUtil->pts),
            x2i2(m_derivUtil->pts), y2i2(m_derivUtil->pts);

        x1i2 = m_derivUtil->VdmD[0] * X;
        y1i2 = m_derivUtil->VdmD[0] * Y;
        x2i2 = m_derivUtil->VdmD[1] * X;
        y2i2 = m_derivUtil->VdmD[1] * Y;

        for (int j = 0; j < m_derivUtil->pts; j++)
        {
            NekDouble jacDet = x1i2(j) * y2i2(j) - x2i2(j) * y1i2(j);
            jacDet /= maps[j][9];
            mx = max(mx, jacDet);
            mn = min(mn, jacDet);
        }
    }
    else if (m_dim == 3)
    {
        NekVector<NekDouble> X(nodes.size()), Y(nodes.size()), Z(nodes.size());
        for (int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
            Z(j) = *nodes[j][2];
        }

        NekVector<NekDouble> x1i(m_derivUtil->pts), y1i(m_derivUtil->pts),
            z1i(m_derivUtil->pts), x2i(m_derivUtil->pts), y2i(m_derivUtil->pts),
            z2i(m_derivUtil->pts), x3i(m_derivUtil->pts), y3i(m_derivUtil->pts),
            z3i(m_derivUtil->pts);

        x1i = m_derivUtil->VdmD[0] * X;
        y1i = m_derivUtil->VdmD[0] * Y;
        z1i = m_derivUtil->VdmD[0] * Z;
        x2i = m_derivUtil->VdmD[1] * X;
        y2i = m_derivUtil->VdmD[1] * Y;
        z2i = m_derivUtil->VdmD[1] * Z;
        x3i = m_derivUtil->VdmD[2] * X;
        y3i = m_derivUtil->VdmD[2] * Y;
        z3i = m_derivUtil->VdmD[2] * Z;

        for (int j = 0; j < m_derivUtil->pts; j++)
        {
            DNekMat dxdz(3, 3, 1.0, eFULL);
            dxdz(0, 0) = x1i(j);
            dxdz(0, 1) = x2i(j);
            dxdz(0, 2) = x3i(j);
            dxdz(1, 0) = y1i(j);
            dxdz(1, 1) = y2i(j);
            dxdz(1, 2) = y3i(j);
            dxdz(2, 0) = z1i(j);
            dxdz(2, 1) = z2i(j);
            dxdz(2, 2) = z3i(j);

            NekDouble jacDet =
                dxdz(0, 0) *
                    (dxdz(1, 1) * dxdz(2, 2) - dxdz(2, 1) * dxdz(1, 2)) -
                dxdz(0, 1) *
                    (dxdz(1, 0) * dxdz(2, 2) - dxdz(2, 0) * dxdz(1, 2)) +
                dxdz(0, 2) *
                    (dxdz(1, 0) * dxdz(2, 1) - dxdz(2, 0) * dxdz(1, 1));

            mx = max(mx, jacDet);
            mn = min(mn, jacDet);
        }
    }

    m_minJac = mn;
}

void ElUtil::UpdateMapping()
{
    if (m_interp.GetInField())
    {
        if (!m_interpField)
        {
            Array<OneD, Array<OneD, NekDouble> > centre(m_dim + 1);
            for (int i = 0; i < m_dim + 1; ++i)
            {
                centre[i] = Array<OneD, NekDouble>(1, 0.0);
            }

            vector<string> fieldNames;
            fieldNames.push_back("");

            map<LibUtilities::PtsInfo, int> ptsInfo =
                LibUtilities::NullPtsInfoMap;

            m_interpField = MemoryManager<LibUtilities::PtsField>
                ::AllocateSharedPtr(m_dim, fieldNames, centre, ptsInfo);
        }

        vector<NodeSharedPtr> nodes = m_el->GetVertexList();

        vector<NekDouble> centre(m_dim, 0.0);
        for (int i = 0; i < nodes.size(); ++i)
        {
            centre[0] += nodes[i]->m_x;
            centre[1] += nodes[i]->m_y;
            if (m_dim > 2)
            {
                centre[2] += nodes[i]->m_z;
            }
        }

        m_interpField->SetPointVal(0, 0, centre[0] / nodes.size());
        m_interpField->SetPointVal(1, 0, centre[1] / nodes.size());
        if (m_dim > 2)
        {
            m_interpField->SetPointVal(2, 0, centre[2] / nodes.size());
        }

        m_interp.CalcWeights(m_interp.GetInField(), m_interpField, true);
        m_interp.Interpolate(m_interp.GetInField(), m_interpField);
    }

    NekDouble scaling = 1.0;

    if (m_interp.GetInField())
    {
        scaling = m_interpField->GetPointVal(m_dim + 0, 0);
    }

    for (int i = 0; i < m_maps.size(); ++i)
    {
        for (int j = 0; j < 9; ++j)
        {
            maps[i][j]    = m_maps[i][j] / scaling;
            mapsStd[i][j] = m_mapsStd[i][j] / scaling;
        }

        if (m_dim == 2)
        {
            maps[i][9]    = m_maps[i][9] * scaling * scaling;
            mapsStd[i][9] = m_mapsStd[i][9] * scaling * scaling;
        }
        else if (m_dim == 3)
        {
            maps[i][9]    = m_maps[i][9] * scaling * scaling * scaling;
            mapsStd[i][9] = m_mapsStd[i][9] * scaling * scaling * scaling;
        }
        else
        {
            ASSERTL0(false, "not coded");
        }
    }
}

ElUtilJob *ElUtil::GetJob(bool update)
{
    return new ElUtilJob(this, update);
}
}
}
