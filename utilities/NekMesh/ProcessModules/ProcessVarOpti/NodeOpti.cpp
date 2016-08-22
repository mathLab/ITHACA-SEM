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

NekDouble NodeOpti::ModifyHessian()
{
    Array<OneD, NekDouble> eigR(3);
    Array<OneD, NekDouble> eigI(3);
    Array<OneD, NekDouble> eigV(9);
    NekMatrix<NekDouble> H(3,3);
    H(0,0) = G[3];
    H(1,0) = G[4];
    H(0,1) = H(1,0);
    H(2,0) = G[5];
    H(0,2) = H(2,0);
    H(1,1) = G[6];
    H(2,1) = G[7];
    H(1,2) = H(2,1);
    H(2,2) = G[8];

    int nVel = 3;
    char jobvl = 'N', jobvr = 'V';
    int worklen = 8*nVel, info;

    DNekMat eval   (nVel, nVel, 0.0, eDIAGONAL);
    DNekMat evec   (nVel, nVel, 0.0, eFULL);
    Array<OneD, NekDouble> vl  (nVel*nVel);
    Array<OneD, NekDouble> work(worklen);
    Array<OneD, NekDouble> wi  (nVel);

    Lapack::Dgeev(jobvl, jobvr, nVel, H.GetRawPtr(), nVel,
                  &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                  &(evec.GetPtr())[0], nVel,
                  &work[0], worklen, info);

    DNekMat evecT = evec;
    evecT.Transpose();

    NekDouble delta = 1e-6;

    if(eval(0,0) < 0.0 || eval(1,1) < 0.0 || eval(2,2) < 0.0)
    {
        for(int i = 0; i < 3; i++)
        {
            eval(i,i) = fabs(eval(i,i));
            if(eval(i,i) < delta)
            {
                eval(i,i) = delta;
            }
        }

        DNekMat Hb = evec * eval * evecT;

        G[3] = Hb(0,0);
        G[4] = Hb(1,0);
        G[5] = Hb(2,0);
        G[6] = Hb(1,1);
        G[7] = Hb(2,1);
        G[8] = Hb(2,2);

        Lapack::Dgeev(jobvl, jobvr, nVel, Hb.GetRawPtr(), nVel,
                      &(eval.GetPtr())[0], &wi[0], &vl[0], nVel,
                      &(evec.GetPtr())[0], nVel,
                      &work[0], worklen, info);
        if(eval(0,0) < 0.0 || eval(1,1) < 0.0 || eval(2,2) < 0.0)
        {
            cout << "still negative" << endl;
        }
    }

    return min(min(eval(0,0),eval(1,1)),eval(2,2));
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
        NekDouble newVal;

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

            newVal = GetFunctional<2>(true,false);
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
        res->func += newVal;
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

    NekDouble minHes = ModifyHessian();

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
        NekDouble newVal;
        // Dot product of p_k with gradient
        NekDouble tmp = (G[0] * delX + G[1] * delY + G[2] * delZ) * -1.0;

        while (alpha > alphaTol())
        {
            // Update node
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;

            newVal = GetFunctional<3>(true,false);
            //dont need the hessian again this function updates G to be the new
            //location
            //
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
        res->func +=newVal;
        mtx.unlock();
    }
}

NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}

}
}
