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
    /*
    typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
                                  Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    return Type::Instance();
    */
    //typedef Loki::SingletonHolder<NodeOptiFactory, Loki::CreateUsingNew,
    //                              Loki::NoDestroy, Loki::ClassLevelLockable> Type;
    //return Type::Instance();

    static NodeOptiFactory asd;
    return asd;
}

const NekDouble NodeOpti::gam = numeric_limits<float>::epsilon();

void NodeOpti::CalcMinJac()
{
    minJac = numeric_limits<double>::max();
    map<LibUtilities::ShapeType,vector<ElUtilSharedPtr> >::iterator typeIt;
    for(typeIt = data.begin(); typeIt != data.end(); typeIt++)
    {
        for(int i = 0; i < typeIt->second.size(); i++)
        {
            minJac = min(minJac,typeIt->second[i]->GetMinJac());
        }
    }
}

NekDouble dir[12][3] = {{1,0,0},
                       {-1,0,0},
                       {0,1,0},
                       {0,-1,0},
                       {0,0,1},
                       {0,0,-1},
                       {1,1,0},
                       {-1,-1,0},
                       {1,0,1},
                       {-1,0,-1},
                       {0,1,1},
                       {0,-1,-1}};

int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<2>(minJacNew);
    NekDouble newVal = currentW;

    // Gradient already zero
    if (G[0]*G[0] + G[1]*G[1] > gradTol())
    {
        //needs to optimise
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;

        /*Array<OneD, NekDouble> tmp = G;
        //cout << endl;
        //cout << G[0] << " " << G[1] << " " << G[2] << " " << G[3] << " " << G[4] << endl;

        Array<OneD, NekDouble> GW(6);
        NekDouble dx = 1e-8;
        for(int i = 0; i < 4; i++)
        {
            node->m_x = xc + dir[i][0] * dx;
            node->m_y = yc + dir[i][1] * dx;
            //node->m_z = zc + dir[i][2] * dx;
            GW[i] = GetFunctional<2>(minJacNew,false);
        }
        node->m_x = xc + dir[6][0] * dx;
        node->m_y = yc + dir[6][1] * dx;
        //node->m_z = zc + dir[i][2] * dx;
        GW[4] = GetFunctional<2>(minJacNew,false);
        node->m_x = xc + dir[7][0] * dx;
        node->m_y = yc + dir[7][1] * dx;
        //node->m_z = zc + dir[i][2] * dx;
        GW[5] = GetFunctional<2>(minJacNew,false);

        Array<OneD, NekDouble> grad(5);
        grad[0] = (GW[0] - GW[1]) /dx/2.0;
        grad[1] = (GW[2] - GW[3]) /dx/2.0;
        grad[2] = (GW[0] + GW[1] - 2*currentW) / dx / dx;
        grad[3] = (GW[4] + GW[5] - GW[0] - GW[1] - GW[2] - GW[3] + 2*currentW) / 2 / dx / dx;
        grad[4] = (GW[2] + GW[3] - 2*currentW) / dx / dx;

        //cout << grad[0] << " " << grad[1] << " " << grad[2] << " " << grad[3] << " " << grad[4] << endl;
        node->m_x = xc;
        node->m_y = yc;
        //node->m_z = zc;
        G = tmp;*/

        Array<OneD, NekDouble> sk(2), dk(2);
        NekDouble val;

        // Calculate minimum eigenvalue
        MinEigen<2>(val, dk);

        if (val < 1e-6)
        {
            // Add constant identity to Hessian matrix.
            G[2] += 1e-6 - val;
            G[4] += 1e-6 - val;
        }

        sk[0] = -1.0/(G[2]*G[4]-G[3]*G[3])*(G[4]*G[0] - G[3]*G[1]);
        sk[1] = -1.0/(G[2]*G[4]-G[3]*G[3])*(G[2]*G[1] - G[3]*G[0]);

        bool found  = false;
        NekDouble pg = (G[0]*sk[0]+G[1]*sk[1]);
        //normal gradient line Search
        NekDouble alpha    = 1.0;
        NekDouble hes = sk[0] * (sk[0]*G[2] + sk[1]*G[3]) +
                        sk[1] * (sk[0]*G[3] + sk[1]*G[4]);
        hes = min(hes,0.0);

        while (alpha > alphaTol())
        {
            // Update node
            node->m_x = xc + alpha * sk[0];
            node->m_y = yc + alpha * sk[1];

            newVal = GetFunctional<2>(minJacNew,false);
            //dont need the hessian again this function updates G to be the new
            //location
            //
            // Wolfe conditions
            if (newVal <= currentW + c1() * (alpha*pg+ 0.5*alpha*alpha*hes))
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if (!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;

            mtx.lock();
            res->nReset[2]++;
            mtx.unlock();
        }
        else
        {
            minJac = minJacNew;

            mtx.lock();
            if(alpha < 1.0)
            {
                res->alphaI++;
            }
            mtx.unlock();
        }

        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+
                            (node->m_y-yc)*(node->m_y-yc)),
                       res->val);
        mtx.unlock();
    }

    mtx.lock();
    res->func += newVal;
    mtx.unlock();
}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");

void NodeOpti3D3D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<3>(minJacNew);
    NekDouble newVal = currentW;

    if(G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
    {
        //needs to optimise
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;

        Array<OneD, NekDouble> sk(3), dk(3);
        NekDouble val;

        // Calculate minimum eigenvalue
        MinEigen<3>(val, dk);

        if (val < 1e-4)
        {
            // Add constant identity to Hessian matrix.
            G[3] += 1e-4 - val;
            G[6] += 1e-4 - val;
            G[8] += 1e-4 - val;
        }

        //calculate sk
        NekDouble det = G[3]*(G[6]*G[8]-G[7]*G[7])
                       -G[4]*(G[4]*G[8]-G[5]*G[7])
                       +G[5]*(G[4]*G[7]-G[5]*G[6]);

        sk[0] = G[0]*(G[6]*G[8]-G[7]*G[7]) +
                G[1]*(G[5]*G[7]-G[4]*G[8]) +
                G[2]*(G[4]*G[7]-G[3]*G[7]);
        sk[1] = G[0]*(G[7]*G[5]-G[4]*G[5]) +
                G[1]*(G[3]*G[8]-G[5]*G[5]) +
                G[2]*(G[4]*G[5]-G[3]*G[7]);
        sk[2] = G[0]*(G[4]*G[7]-G[6]*G[5]) +
                G[1]*(G[4]*G[5]-G[3]*G[7]) +
                G[2]*(G[3]*G[6]-G[4]*G[4]);

        sk[0] /= det * -1.0;
        sk[1] /= det * -1.0;
        sk[2] /= det * -1.0;

        bool found  = false;

        NekDouble pg = (G[0]*sk[0]+G[1]*sk[1]+G[2]*sk[2]);

        //normal gradient line Search
        NekDouble alpha    = 1.0;
        NekDouble hes = sk[0] * (sk[0]*G[3] + sk[1]*G[4] + sk[2]*G[5]) +
                        sk[1] * (sk[0]*G[4] + sk[1]*G[6] + sk[2]*G[7]) +
                        sk[2] * (sk[0]*G[5] + sk[1]*G[7] + sk[2]*G[8]);
        hes = min(hes,0.0);

        while (alpha > alphaTol())
        {
            // Update node
            node->m_x = xc + alpha * sk[0];
            node->m_y = yc + alpha * sk[1];
            node->m_z = zc + alpha * sk[2];

            newVal = GetFunctional<3>(minJacNew,false);
            //dont need the hessian again this function updates G to be the new
            //location
            //
            // Wolfe conditions
            if (newVal <= currentW + c1() * (
                    alpha*pg+ 0.5*alpha*alpha*hes))
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            node->m_x = xc;
            node->m_y = yc;
            node->m_z = zc;

            mtx.lock();
            res->nReset[2]++;
            mtx.unlock();
        }
        else
        {
            minJac = minJacNew;

            mtx.lock();
            if(alpha < 1.0)
            {
                res->alphaI++;
            }
            mtx.unlock();
        }

        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
    mtx.lock();
    res->func += newVal;
    mtx.unlock();
}

NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}

}
}
