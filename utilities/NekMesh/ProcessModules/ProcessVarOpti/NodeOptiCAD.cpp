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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "NodeOptiCAD.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx;

int NodeOpti1D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    13, NodeOpti1D3D::create, "1D3D");

void NodeOpti1D3D::Optimise()
{
    CalcDX();

    CalcMinJac();

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0]) > 1e-10)
    {
        //needs to optimise
        NekDouble tc = node->GetCADCurveInfo(curve->GetId());
        NekDouble currentW = GetFunctional<3>();
        NekDouble functional;
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delT;
        NekDouble nt;
        Array<OneD, NekDouble> p;

        delT = G[0] / G[1];

        Array<OneD, NekDouble> bd = curve->Bounds();

        bool found = false;
        while(alpha > 1e-10)
        {
            nt = tc - alpha * delT;
            if(nt < bd[0] || nt > bd[1])
            {
                alpha /= 2.0;
                continue;
            }
            p = curve->P(nt);
            node->m_x = p[0];
            node->m_y = p[1];
            node->m_z = p[2];
            functional = GetFunctional<3>();
            if(functional < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            nt = tc;
            p = curve->P(nt);
            node->m_x = p[0];
            node->m_y = p[1];
            node->m_z = p[2];
            functional = currentW;
            // cout << "warning: had to reset node" << endl;
        }
        else
        {
            node->MoveCurve(p,curve->GetId(),nt);
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
}

int NodeOpti2D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    23, NodeOpti2D3D::create, "1D3D");

void NodeOpti2D3D::Optimise()
{
    CalcDX();

    CalcMinJac();

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1]) > 1e-10)
    {
        //needs to optimise
        Array<OneD, NekDouble> uvc = node->GetCADSurfInfo(surf->GetId());
        NekDouble currentW = GetFunctional<3>();
        NekDouble functional;
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        Array<OneD, NekDouble> uvt(2);
        Array<OneD, NekDouble> p;
        NekDouble delU = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[3]*G[0] - G[4]*G[1]);
        NekDouble delV = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[2]*G[1] - G[4]*G[0]);

        Array<OneD, NekDouble> bd = surf->GetBounds();

        bool found = false;
        while(alpha > 1e-10)
        {
            uvt[0] = uvc[0] - alpha * delU;
            uvt[1] = uvc[1] - alpha * delV;

            if(uvt[0] < bd[0] || uvt[0] > bd[1] ||
               uvt[1] < bd[2] || uvt[1] > bd[3])
            {
                alpha /= 2.0;
                continue;
            }

            p = surf->P(uvt);
            node->m_x = p[0];
            node->m_y = p[1];
            node->m_z = p[2];
            functional = GetFunctional<3>();
            if(functional < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            p = surf->P(uvc);
            node->m_x = p[0];
            node->m_y = p[1];
            node->m_z = p[2];
            functional = currentW;
            // cout << "warning: had to reset node" << endl;
        }
        else
        {
            node->Move(p,surf->GetId(),uvt);
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
}

Array<OneD, NekDouble> NodeOpti1D3D::GetGrad()
{
    NekDouble tc = node->GetCADCurveInfo(curve->GetId());
    Array<OneD, NekDouble> d1 = curve->D1(tc);

    NekDouble dr = sqrt(d1[3]*d1[3] + d1[4]*d1[4] + d1[5]*d1[5]);

    NekDouble dt = dx / dr;

    vector<NekDouble> w(3);

    w[0] = GetFunctional<3>();

    NekDouble nt = tc + dt;
    Array<OneD, NekDouble> p = curve->P(nt);
    node->m_x = p[0];
    node->m_y = p[1];
    node->m_z = p[2];
    w[1] = GetFunctional<3>();

    nt = tc - dt;
    p = curve->P(nt);
    node->m_x = p[0];
    node->m_y = p[1];
    node->m_z = p[2];
    w[2] = GetFunctional<3>();

    nt = tc;
    p = curve->P(nt);
    node->m_x = p[0];
    node->m_y = p[1];
    node->m_z = p[2];

    Array<OneD, NekDouble> ret(2,0.0);

    ret[0] = (w[1] - w[2]) / 2.0 / dt;
    ret[1] = (w[1] + w[2] - 2.0*w[0]) / dt / dt;

    return ret;
}

Array<OneD, NekDouble> NodeOpti2D3D::GetGrad()
{
    Array<OneD, NekDouble> uvc = node->GetCADSurfInfo(surf->GetId());
    Array<OneD, NekDouble> d1 = surf->D1(uvc);

    NekDouble dru = sqrt(d1[3]*d1[3] + d1[4]*d1[4] + d1[5]*d1[5]);
    NekDouble drv = sqrt(d1[6]*d1[6] + d1[7]*d1[7] + d1[8]*d1[8]);

    NekDouble du = dx / dru;
    NekDouble dv = dx / drv;

    vector<NekDouble> w(7);

    for(int i = 0; i < 7; i++)
    {
        Array<OneD, NekDouble> uvt(2);
        uvt[0] = uvc[0] + dir[i][0] * du;
        uvt[1] = uvc[1] + dir[i][1] * dv;
        Array<OneD, NekDouble> p = surf->P(uvt);
        node->m_x = p[0];
        node->m_y = p[1];
        node->m_z = p[2];
        w[i] = GetFunctional<3>();
    }

    Array<OneD, NekDouble> p = surf->P(uvc);
    node->m_x = p[0];
    node->m_y = p[1];
    node->m_z = p[2];

    Array<OneD, NekDouble> ret(5,0.0);

    //ret[0] d/dx
    //ret[1] d/dy

    //ret[2] d2/dx2
    //ret[3] d2/dy2
    //ret[4] d2/dxdy


    ret[0] = (w[1] - w[4]) / 2.0 / du;
    ret[1] = (w[3] - w[6]) / 2.0 / dv;
    ret[2] = (w[1] + w[4] - 2.0*w[0]) / du / du;
    ret[3] = (w[3] + w[6] - 2.0*w[0]) / dv / dv;
    ret[4] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / du / dv;

    return ret;
}

}
}
