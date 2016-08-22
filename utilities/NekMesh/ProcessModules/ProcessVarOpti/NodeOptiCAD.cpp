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
#include "Evaluator.hxx"

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
    CalcMinJac();

    NekDouble currentW = GetFunctional<3>();

    if (G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
    {
        //modify the gradient to be on the cad system
        ProcessGradient();

        //needs to optimise
        NekDouble tc = node->GetCADCurveInfo(curve->GetId());
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delT;
        NekDouble nt;
        Array<OneD, NekDouble> p;

        delT = G[0] / G[1];

        Array<OneD, NekDouble> bd = curve->Bounds();

        // Dot product of p_k with gradient
        NekDouble tmp = (G[0] * delT) * -1.0;

        bool found = false;
        NekDouble newVal;

        while (alpha > alphaTol())
        {
            // Update node
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
            node->MoveCurve(p,curve->GetId(),nt);

            newVal = GetFunctional<3>(true,false);
            ProcessGradient();

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * tmp &&
                -1.0 * (G[0] * delT) >= c2() * tmp)
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
            node->MoveCurve(p,curve->GetId(),nt);

            mtx.lock();
            res->nReset++;
            mtx.unlock();
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        res->func += newVal;
        mtx.unlock();
    }
}

int NodeOpti2D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    23, NodeOpti2D3D::create, "1D3D");

void NodeOpti2D3D::Optimise()
{
    CalcMinJac();

    NekDouble currentW = GetFunctional<3>();

    if (G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
    {
        //cout << "***********************************" << endl;
        //cout << G[0] << " " << G[1] << " " << G[2] << endl
        //     << G[3] << " " << G[4] << " " << G[5] << " " << G[6] << " " << G[7] << " " << G[8] << endl;
        //cout << endl;

        //modify the gradient to be on the cad system
        ProcessGradient();

        //cout << endl << G[0] << " " << G[1] << endl << G[2] << " " << G[3] << " " << G[4] << endl;
        //cout << endl;

        //needs to optimise
        Array<OneD, NekDouble> uvc = node->GetCADSurfInfo(surf->GetId());
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        Array<OneD, NekDouble> uvt(2);
        Array<OneD, NekDouble> p;

        NekDouble delU = 1.0/(G[2]*G[4]-G[3]*G[3])*(G[4]*G[0] - G[3]*G[1]);
        NekDouble delV = 1.0/(G[2]*G[4]-G[3]*G[3])*(G[2]*G[1] - G[3]*G[0]);

        //cout << delU << " " << delV << endl << endl;

        // Dot product of p_k with gradient
        NekDouble tmp = (G[0] * delU + G[1] * delV) * -1.0;

        Array<OneD, NekDouble> bd = surf->GetBounds();

        bool found = false;
        NekDouble newVal;
        while(alpha > alphaTol())
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
            node->Move(p,surf->GetId(),uvt);

            newVal = GetFunctional<3>(true,false);
            ProcessGradient();

            //cout << currentW << " " <<  newVal << endl;

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * tmp &&
                -1.0 * (G[0] * delU + G[1] * delV) >= c2() * tmp)
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
            node->Move(p,surf->GetId(),uvc);

            mtx.lock();

            /*NekDouble newVal = GetFunctional<3>();
            cout << node->m_id << " " << delU << " " << delV << endl
                 << "Gradient 3D\t" << G[0] << " " << G[1] << " " << G[2] << endl
                 << "Hessian 3D\t" << G[3] << " " << G[4] << " " << G[5] << " " << G[6] << " " << G[7] << " " << G[8] << endl;
            ProcessGradient();
            cout << "Gradient surf\t" << G[0] << " " << G[1] << endl
                 << "Hessian surf\t" << G[2] << " " << G[3] << " " << G[4] << endl
                 << (G[2]*G[4]-G[3]*G[3]) << " " << minJac << endl << endl;*/
            res->nReset++;
            mtx.unlock();
        }

        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        res->func += newVal;
        mtx.unlock();
    }
}

void NodeOpti1D3D::ProcessGradient()
{
    NekDouble tc = node->GetCADCurveInfo(curve->GetId());
    Array<OneD, NekDouble> grad = G;
    G = Array<OneD, NekDouble>(2,0.0);

    // Grab first and second order CAD derivatives
    Array<OneD, NekDouble> d2 = curve->D2(tc);

    // Multiply gradient by derivative of CAD
    G[0] = grad[0] * d2[3] + grad[1] * d2[4] + grad[2] * d2[5];

    // Second order: product rule of above, so multiply gradient by second
    // order CAD derivatives and Hessian by gradient of CAD
    G[1] = grad[0] * d2[6] + grad[1] * d2[7] + grad[2] * d2[8]
        + d2[3] * (grad[3] * d2[3] + grad[4] * d2[4] + grad[5] * d2[5])
        + d2[4] * (grad[4] * d2[3] + grad[6] * d2[4] + grad[7] * d2[5])
        + d2[5] * (grad[5] * d2[3] + grad[7] * d2[4] + grad[8] * d2[5]);
}

void NodeOpti2D3D::ProcessGradient()
{
    Array<OneD, NekDouble> uvc = node->GetCADSurfInfo(surf->GetId());

    Array<OneD, NekDouble> grad = G;
    G = Array<OneD, NekDouble>(5,0.0);

    Array<OneD, NekDouble> d2 = surf->D2(uvc);
    //r[0]   x
    //r[1]   y
    //r[2]   z
    //r[3]   dx/du a
    //r[4]   dy/du b
    //r[5]   dz/du c
    //r[6]   dx/dv d
    //r[7]   dy/dv e
    //r[8]   dz/dv f
    //r[9]   d2x/du2
    //r[10]  d2y/du2
    //r[11]  d2z/du2
    //r[12]  d2x/dv2
    //r[13]  d2y/dv2
    //r[14]  d2z/dv2
    //r[15]  d2x/dudv
    //r[16]  d2y/dudv
    //r[17]  d2z/dudv
    //
    //grad[0] d/dx
    //grad[1] d/dy
    //grad[2] d/dz
    //grad[3] d2/dx2
    //grad[4] d2/dxdy
    //grad[5] d2/dxdz
    //grad[6] d2/dy2
    //grad[7] d2/dydz
    //grad[8] d2/dz2

    // Gradients
    G[0] = d2[3] * grad[0] + d2[4] * grad[1] + d2[5] * grad[2];
    G[1] = d2[6] * grad[0] + d2[7] * grad[1] + d2[8] * grad[2];

    G[2]  = grad[0] * d2[9] + grad[1] * d2[10] + grad[2] * d2[11]
          + grad[3]*d2[3]*d2[3] + grad[6]*d2[4]*d2[4] + grad[8]*d2[5]*d2[5]
          + 2.0*grad[4]*d2[3]*d2[4] + 2.0*grad[5]*d2[3]*d2[5] + 2.0*grad[7]*d2[4]*d2[5];

    G[4]  = grad[0] * d2[12] + grad[1] * d2[13] + grad[2] * d2[14]
          + grad[3]*d2[6]*d2[6] + grad[6]*d2[7]*d2[7] + grad[8]*d2[8]*d2[8]
          + 2.0*grad[4]*d2[6]*d2[7] + 2.0*grad[5]*d2[6]*d2[8] + 2.0*grad[7]*d2[7]*d2[8];

    G[3]  = grad[0] * d2[15] + grad[1] * d2[16] + grad[2] * d2[17]
          + grad[3]*d2[3]*d2[6] + grad[6]*d2[4]*d2[7] + grad[8]*d2[5]*d2[8]
          + grad[4]*(d2[3]*d2[7] + d2[4]*d2[6])
          + grad[5]*(d2[3]*d2[8] + d2[5]*d2[6])
          + grad[7]*(d2[4]*d2[8] + d2[5]*d2[7]);


    //G[2] = 1.0;
    ///G[3] = 0.0;
    //G[4] = 1.0;
}

}
}
