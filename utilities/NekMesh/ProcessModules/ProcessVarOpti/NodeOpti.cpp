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

#include "NodeOpti.h"
#include "Evaluator.hxx"

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

    for(int i = 0; i < data.size(); i++)
    {
        minJac = min(minJac, data[i]->minJac);
    }
}

int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise()
{
    CalcMinJac();

    NekDouble currentW = GetFunctional<2>();

    // Gradient already zero
    if (G[0]*G[0] + G[1]*G[1] > gradTol())
    {
        NekDouble alpha    = 1.0;
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;

        bool      found    = false;

        // Search direction
        NekDouble delX = 1.0/(G[2]*G[4]-G[3]*G[3])*(G[4]*G[0] - G[3]*G[1]);
        NekDouble delY = 1.0/(G[2]*G[4]-G[3]*G[3])*(G[2]*G[1] - G[3]*G[0]);
        // Dot product of p_k with gradient
        NekDouble tmp = (G[0] * delX + G[1] * delY) * -1.0;

        while (alpha > alphaTol())
        {
            // Update node
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;

            NekDouble newVal = GetFunctional<2>(true,false);
            //dont need the hessian again this function updates G to be the new
            //location

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * tmp &&
                -1.0 * (G[0] * delX + G[1] * delY) >= c2() * tmp)
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
            res->nReset++;
            mtx.unlock();
        }

        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+
                            (node->m_y-yc)*(node->m_y-yc)),
                       res->val);
        mtx.unlock();
    }
}

int NodeOpti3D3D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    33, NodeOpti3D3D::create, "3D3D");

void NodeOpti3D3D::Optimise()
{
    CalcMinJac();

    NekDouble currentW = GetFunctional<3>();

    if(G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
    {
        //needs to optimise
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delX, delY, delZ;

        NekDouble det = G[3]*(G[6]*G[8]-G[7]*G[7])
                       -G[4]*(G[4]*G[8]-G[5]*G[7])
                       +G[5]*(G[4]*G[7]-G[5]*G[6]);

        delX = G[0]*(G[6]*G[8]-G[7]*G[7]) +
               G[1]*(G[5]*G[7]-G[4]*G[8]) +
               G[2]*(G[4]*G[7]-G[3]*G[7]);
        delY = G[0]*(G[7]*G[5]-G[4]*G[5]) +
               G[1]*(G[3]*G[8]-G[5]*G[5]) +
               G[2]*(G[4]*G[5]-G[3]*G[7]);
        delZ = G[0]*(G[4]*G[7]-G[6]*G[5]) +
               G[1]*(G[4]*G[5]-G[3]*G[7]) +
               G[2]*(G[3]*G[6]-G[4]*G[4]);

        delX /= det;
        delY /= det;
        delZ /= det;

        bool found = false;
        // Dot product of p_k with gradient
        NekDouble tmp = (G[0] * delX + G[1] * delY + G[2] * delZ) * -1.0;

        while (alpha > alphaTol())
        {
            // Update node
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;

            NekDouble newVal = GetFunctional<3>(true,false);
            //dont need the hessian again this function updates G to be the new
            //location

            // Wolfe conditions
            if (newVal <= currentW + c1() * alpha * tmp &&
                -1.0 * (G[0] * delX + G[1] * delY + G[2] * delZ) >= c2() * tmp)
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
            res->nReset++;
            mtx.unlock();
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
}

NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}

}
}
