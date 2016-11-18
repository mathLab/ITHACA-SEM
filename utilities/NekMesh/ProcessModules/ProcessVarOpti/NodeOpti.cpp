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
            minJac = min(minJac,typeIt->second[i]->minJac);
        }
    }
}

int NodeOpti2D2D::m_type = GetNodeOptiFactory().RegisterCreatorFunction(
    22, NodeOpti2D2D::create, "2D2D");

void NodeOpti2D2D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<2>(minJacNew);
    NekDouble newVal = currentW;

    /*cout << endl;
    cout << G[0] << " " << G[1] << " " << G[2] << " " << G[3] << " " << G[4] << endl;

    Array<OneD, NekDouble> Gt = G;

    NekDouble xc       = node->m_x;
    NekDouble yc       = node->m_y;

    vector<NekDouble> d;
    for(int i = 0; i < 6; i++)
    {
        node->m_x = xc + 1e-8*dir[i][0];
        node->m_y = yc + 1e-8*dir[i][1];

        d.push_back(GetFunctional<2>(false,false));
    }

    G = Array<OneD, NekDouble>(5);
    G[0] = (d[1] - d[3]) / 2e-8;
    G[1] = (d[2] - d[4]) / 2e-8;
    G[2] = (d[1] + d[3] - 2*currentW) / 1e-16;
    G[3] = (d[0] - d[1] - d[2] + 2*currentW - d[3] - d[4] + d[5]) / 2e-16;
    G[4] = (d[2] + d[4] - 2*currentW) / 1e-16;

    cout << G[0] << " " << G[1] << " " << G[2] << " " << G[3] << " " << G[4] << endl;

    G = Gt;*/

    // Gradient already zero
    if (G[0]*G[0] + G[1]*G[1] > gradTol())
    {
        //needs to optimise
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;

        Array<OneD, NekDouble> sk(2), dk(2);
        bool DNC = false;
        NekDouble lhs;

        int def = IsIndefinite<2>();
        if(def)
        {
            //the dk vector needs calculating
            NekDouble val;
            MinEigen<2>(val,dk);

            if(dk[0]*G[0] + dk[1]*G[1]> 0.0)
            {
                for(int i = 0; i < 2; i++)
                {
                    dk[i] *= -1.0;
                }
            }

            lhs = dk[0] * (dk[0]*G[2] + dk[1]*G[3]) +
                  dk[1] * (dk[0]*G[3] + dk[1]*G[4]);

            ASSERTL0(lhs < 0.0 , "weirdness");

            DNC = true;
        }

        sk[0] = -1.0/(G[2]*G[4]-G[3]*G[3])*(G[4]*G[0] - G[3]*G[1]);
        sk[1] = -1.0/(G[2]*G[4]-G[3]*G[3])*(G[2]*G[1] - G[3]*G[0]);

        bool runDNC = false; //so want to make this varible runDMC
        bool found  = false;

        NekDouble skmag = sqrt(sk[0]*sk[0] + sk[1]*sk[1]);

        if(DNC)
        {
            runDNC = !((G[0]*sk[0]+G[1]*sk[1])/skmag <=
                        2.0*(0.5*lhs + G[0]*dk[0]+G[1]*dk[1]));
        }

        NekDouble pg = (G[0]*sk[0]+G[1]*sk[1]);

        if(!runDNC)
        {
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

                newVal = GetFunctional<2>(minJacNew,false,false);
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
        }
        else
        {
            NekDouble sig = 1.0;
            NekDouble beta = 0.5;
            int l = 0;
            NekDouble alpha = pow(beta,l);

            NekDouble hes = lhs;

            pg = (G[0]*dk[0]+G[1]*dk[1]);

            //choose whether to do forward or reverse line search
            node->m_x = xc + dk[0];
            node->m_y = yc + dk[1];

            newVal = GetFunctional<2>(minJacNew,false,false);

            if(newVal <= currentW + c1() * (
                pg + 0.5*hes))
            {
                //this is a minimser so see if we can extend further
                while (l > -10)
                {
                    // Update node
                    node->m_x = xc + alpha * dk[0];
                    node->m_y = yc + alpha * dk[1];

                    newVal = GetFunctional<2>(minJacNew,false,false);

                    node->m_x = xc + alpha/beta * dk[0];
                    node->m_y = yc + alpha/beta * dk[1];

                    NekDouble dbVal = GetFunctional<2>(minJacNew,false,false);

                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes) &&
                        dbVal > currentW + c1() *(
                        alpha/beta*pg + 0.5*alpha*alpha*hes/beta/beta))
                    {
                        found = true;
                        break;
                    }

                    l--;
                    alpha = pow(beta,l);
                }
            }
            else
            {
                //this is not a minimser so reverse line search
                while (alpha > alphaTol())
                {
                    // Update node
                    node->m_x = xc + alpha * dk[0];
                    node->m_y = yc + alpha * dk[1];

                    newVal = GetFunctional<2>(minJacNew,false,false);

                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes))
                    {
                        found = true;
                        break;
                    }

                    l++;
                    alpha = pow(beta,l);
                }
            }
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


void NodeOpti3D3D::Optimise()
{
    NekDouble minJacNew;
    NekDouble currentW = GetFunctional<3>(minJacNew);
    NekDouble newVal = currentW;

    /*//cout << endl;
    //cout << G[0] << " " << G[1] << " " << G[2] << " " << G[3] << " " << G[4] << " " << G[5] << " " << G[6] << " " << G[7] << " " << G[8] << endl;

    Array<OneD, NekDouble> Gt = G;

    NekDouble xc       = node->m_x;
    NekDouble yc       = node->m_y;
    NekDouble zc       = node->m_z;

    vector<NekDouble> d;
    for(int i = 0; i < 12; i++)
    {
        node->m_x = xc + 1e-6*dir[i][0];
        node->m_y = yc + 1e-6*dir[i][1];
        node->m_z = zc + 1e-6*dir[i][2];

        d.push_back(GetFunctional<3>(false,false));
    }

    G = Array<OneD, NekDouble>(9);
    G[0] = (d[0] - d[1]) / 2e-6;
    G[1] = (d[2] - d[3]) / 2e-6;
    G[2] = (d[4] - d[5]) / 2e-6;

    G[3] = (d[0] + d[1] - 2*currentW) / 1e-12;
    G[6] = (d[2] + d[3] - 2*currentW) / 1e-12;
    G[8] = (d[4] + d[5] - 2*currentW) / 1e-12;

    G[4] = (d[6] - d[0] - d[2] + 2*currentW - d[1] - d[3] + d[7]) / 2e-12;
    G[5] = (d[8] - d[0] - d[4] + 2*currentW - d[1] - d[5] + d[9]) / 2e-12;
    G[7] = (d[10] - d[2] - d[4] + 2*currentW - d[3] - d[5] + d[11]) / 2e-12;

    //cout << G[0] << " " << G[1] << " " << G[2] << " " << G[3] << " " << G[4] << " " << G[5] << " " << G[6] << " " << G[7] << " " << G[8] << endl;

    //G = Gt;*/

    if(G[0]*G[0] + G[1]*G[1] + G[2]*G[2] > gradTol())
    {
        //needs to optimise
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;

        Array<OneD, NekDouble> sk(3), dk(3);
        bool DNC = false;
        NekDouble lhs;

        int def = IsIndefinite<3>();
        if(def)
        {
            //the dk vector needs calculating
            NekDouble val;
            MinEigen<3>(val,dk);

            if(dk[0]*G[0] + dk[1]*G[1] + dk[2]*G[2] > 0.0)
            {
                for(int i = 0; i < 3; i++)
                {
                    dk[i] *= -1.0;
                }
            }

            lhs = dk[0] * (dk[0]*G[3] + dk[1]*G[4] + dk[2]*G[5]) +
                  dk[1] * (dk[0]*G[4] + dk[1]*G[6] + dk[2]*G[7]) +
                  dk[2] * (dk[0]*G[5] + dk[1]*G[7] + dk[2]*G[8]);

            ASSERTL0(lhs < 0.0 , "weirdness");

            DNC = true;
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

        bool runDNC = false; //so want to make this varible runDMC
        bool found  = false;

        NekDouble skmag = sqrt(sk[0]*sk[0] + sk[1]*sk[1] + sk[2]*sk[2]);
        if(DNC)
        {
            runDNC = !((G[0]*sk[0]+G[1]*sk[1]+G[2]*sk[2])/skmag <=
                        2.0*(0.5*lhs + G[0]*dk[0]+G[1]*dk[1]+G[2]*dk[2]));
        }

        NekDouble pg = (G[0]*sk[0]+G[1]*sk[1]+G[2]*sk[2]);

        if(!runDNC)
        {
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

                newVal = GetFunctional<3>(minJacNew,false,false);
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
        }
        else
        {
            NekDouble sig = 1.0;
            NekDouble beta = 0.5;
            int l = 0;
            NekDouble alpha = pow(beta,l);

            NekDouble hes = lhs;

            pg = (G[0]*dk[0]+G[1]*dk[1]+G[2]*dk[2]);

            //choose whether to do forward or reverse line search
            node->m_x = xc + dk[0];
            node->m_y = yc + dk[1];
            node->m_z = zc + dk[2];
            newVal = GetFunctional<3>(minJacNew,false,false);

            if(newVal <= currentW + c1() * (
                pg + 0.5*hes))
            {
                //this is a minimser so see if we can extend further
                while (l > -10)
                {
                    // Update node
                    node->m_x = xc + alpha * dk[0];
                    node->m_y = yc + alpha * dk[1];
                    node->m_z = zc + alpha * dk[2];

                    newVal = GetFunctional<3>(minJacNew,false,false);

                    node->m_x = xc + alpha/beta * dk[0];
                    node->m_y = yc + alpha/beta * dk[1];
                    node->m_z = zc + alpha/beta * dk[2];

                    NekDouble dbVal = GetFunctional<3>(minJacNew,false,false);

                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes) &&
                        dbVal > currentW + c1() *(
                        alpha/beta*pg + 0.5*alpha*alpha*hes/beta/beta))
                    {
                        found = true;
                        break;
                    }

                    l--;
                    alpha = pow(beta,l);
                }
            }
            else
            {
                //this is not a minimser so reverse line search
                while (alpha > alphaTol())
                {
                    // Update node
                    node->m_x = xc + alpha * dk[0];
                    node->m_y = yc + alpha * dk[1];
                    node->m_z = zc + alpha * dk[2];

                    newVal = GetFunctional<3>(minJacNew,false,false);

                    if (newVal <= currentW + c1() * (
                        alpha*pg + 0.5*alpha*alpha*hes))
                    {
                        found = true;
                        break;
                    }

                    l++;
                    alpha = pow(beta,l);
                }
            }
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
