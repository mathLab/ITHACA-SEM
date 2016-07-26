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

#include "NodeOpti.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx;

void NodeOpti::CalcDX()
{
    dx = numeric_limits<double>::max();

    for(int i = 0; i < data.size(); i++)
    {
        dx = min(dx, data[i]->delta);
    }
}

void NodeOpti2D2D::Optimise()
{
    CalcDX();

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1]) > 1e-10)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional<2>();
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble alpha    = 1.0;
        NekDouble delX = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[3]*G[0] - G[4]*G[1]);
        NekDouble delY = 1.0/(G[2]*G[3]-G[4]*G[4])*(G[2]*G[1] - G[4]*G[0]);

        bool found = false;
        while(alpha > 1e-10)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            if(GetFunctional<2>() < currentW)
            {
                found = true;
                break;
            }

            alpha /= 2.0;
        }

        if(!found)
        {
            //reset the node
            node->m_x = xc;
            node->m_y = yc;
            // cout << "warning: had to reset node" << endl;
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)),res->val);
        mtx.unlock();
    }
}


void NodeOpti3D3D::Optimise()
{
    CalcDX();

    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]) > 1e-10)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional<3>();
        NekDouble functional;
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delX;
        NekDouble delY;
        NekDouble delZ;

        NekDouble det = G[3]*(G[4]*G[5]-G[8]*G[8])
                       -G[6]*(G[6]*G[5]-G[7]*G[8])
                       +G[7]*(G[6]*G[8]-G[7]*G[4]);

        delX = G[0]*(G[4]*G[5]-G[8]*G[8]) + G[1]*(G[7]*G[8]-G[6]*G[5]) + G[2]*(G[6]*G[8]-G[7]*G[4]);
        delX /= det;
        delY = G[0]*(G[8]*G[7]-G[6]*G[5]) + G[1]*(G[3]*G[5]-G[7]*G[7]) + G[2]*(G[6]*G[7]-G[3]*G[8]);
        delY /= det;
        delZ = G[0]*(G[6]*G[8]-G[4]*G[7]) + G[1]*(G[6]*G[7]-G[3]*G[8]) + G[2]*(G[3]*G[4]-G[6]*G[6]);
        delZ /= det;


        bool found = false;
        while(alpha > 1e-10)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;
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
            node->m_x = xc;
            node->m_y = yc;
            node->m_z = zc;
            functional = currentW;
            // cout << "warning: had to reset node" << endl;
        }
        mtx.lock();
        res->val = max(sqrt((node->m_x-xc)*(node->m_x-xc)+(node->m_y-yc)*(node->m_y-yc)+
                            (node->m_z-zc)*(node->m_z-zc)),res->val);
        mtx.unlock();
    }
}

void NodeOpti1D3D::Optimise()
{
    CalcDX();

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

void NodeOpti2D3D::Optimise()
{
    CalcDX();

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

NekDouble dir[13][3] = {{  0.0,  0.0,  0.0 },  // 0  (x   , y   , z   )
                        {  1.0,  0.0,  0.0 },  // 1  (x+dx, y   , z   )
                        {  1.0,  1.0,  0.0 },  // 2  (x+dx, y+dy, z   )
                        {  0.0,  1.0,  0.0 },  // 3  (x   , y+dy, z   )
                        { -1.0,  0.0,  0.0 },  // 4  (x-dx, y   , z   )
                        { -1.0, -1.0,  0.0 },  // 5  (x-dx, y-dy, z   )
                        {  0.0, -1.0,  0.0 },  // 6  (x   , y-dy, z   )
                        { -1.0,  0.0, -1.0 },  // 7  (x-dx, y   , z-dz)
                        {  0.0,  0.0, -1.0 },  // 8 (x   , y   , z-dz)
                        {  0.0,  0.0,  1.0 },  // 9 (x   , y   , z+dz)
                        {  1.0,  0.0,  1.0 },  // 10 (x+dx, y   , z+dz)
                        {  0.0,  1.0,  1.0 },  // 11 (x   , y+dy, z+dz)
                        {  0.0, -1.0, -1.0 }}; // 12 (x   , y-dy, z-dz)

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

Array<OneD, NekDouble> NodeOpti2D2D::GetGrad()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    vector<NekDouble> w(9);

    for(int i = 0; i < 7; i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        w[i] = GetFunctional<2>();
    }
    node->m_x = xc;
    node->m_y = yc;

    Array<OneD, NekDouble> ret(5,0.0);

    //ret[0] d/dx
    //ret[1] d/dy

    //ret[3] d2/dx2
    //ret[4] d2/dy2
    //ret[5] d2/dxdy


    ret[0] = (w[1] - w[4]) / 2.0 / dx;
    ret[1] = (w[3] - w[6]) / 2.0 / dx;
    ret[2] = (w[1] + w[4] - 2.0*w[0]) / dx / dx;
    ret[3] = (w[3] + w[6] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / dx / dx;

    return ret;
}

Array<OneD, NekDouble> NodeOpti3D3D::GetGrad()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble zc = node->m_z;

    vector<NekDouble> w;

    for(int i = 0; i < 13; i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        node->m_z = zc + dir[i][2] * dx;
        w.push_back(GetFunctional<3>());
    }
    node->m_x = xc;
    node->m_y = yc;
    node->m_z = zc;

    Array<OneD, NekDouble> ret(9,0.0);

    //ret[0] d/dx
    //ret[1] d/dy
    //ret[2] d/dz

    //ret[3] d2/dx2
    //ret[4] d2/dy2
    //ret[5] d2/dz2
    //ret[6] d2/dxdy
    //ret[7] d2/dxdz
    //ret[8] d2/dydz

    ret[0] = (w[1] - w[4]) / 2.0 / dx;
    ret[1] = (w[3] - w[6]) / 2.0 / dx;
    ret[2] = (w[9] - w[8]) / 2.0 / dx;

    ret[3] = (w[1] + w[4] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[3] + w[6] - 2.0*w[0]) / dx / dx;
    ret[5] = (w[9] + w[8] - 2.0*w[0]) / dx / dx;

    ret[6] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / dx / dx;
    ret[7] = (w[10] - w[1] - w[9] + 2.0*w[0] - w[4] - w[8] + w[7]) / 2.0 / dx / dx;
    ret[8] = (w[11] - w[3] - w[9] + 2.0*w[0] - w[6] - w[8] + w[12]) / 2.0 / dx / dx;

    return ret;
}

template<int DIM> inline NekDouble JacDet(NekDouble *jac)
{
    return 0.0;
}

template<> inline NekDouble JacDet<2>(NekDouble *jac)
{
    return jac[0] * jac[3] - jac[2] * jac[1];
}

template<> inline NekDouble JacDet<3>(NekDouble *jac)
{
    return jac[0]*(jac[4]*jac[8]-jac[5]*jac[7])
          -jac[3]*(jac[1]*jac[8]-jac[2]*jac[7])
          +jac[6]*(jac[1]*jac[5]-jac[2]*jac[4]);
}

template<int DIM> inline NekDouble LinElasTrace(NekDouble *jac)
{
    return 0.0;
}

template<> inline NekDouble LinElasTrace<2>(NekDouble *jac)
{
    return 0.25 * (
        (jac[0]*jac[0]+jac[1]*jac[1]-1.0)*(jac[0]*jac[0]+jac[1]*jac[1]-1.0) +
        (jac[2]*jac[2]+jac[3]*jac[3]-1.0)*(jac[2]*jac[2]+jac[3]*jac[3]-1.0))
        + 0.5 * (
            (jac[0]*jac[2]+jac[1]*jac[3])*(jac[0]*jac[2]+jac[1]*jac[3]));
}

template<> inline NekDouble LinElasTrace<3>(NekDouble *jac)
{
    return 0.25 *(
        (jac[0]*jac[0]+jac[1]*jac[1]+jac[2]*jac[2]-1.0)*
        (jac[0]*jac[0]+jac[1]*jac[1]+jac[2]*jac[2]-1.0) +
        (jac[3]*jac[3]+jac[4]*jac[4]+jac[5]*jac[5]-1.0)*
        (jac[3]*jac[3]+jac[4]*jac[4]+jac[5]*jac[5]-1.0) +
        (jac[6]*jac[6]+jac[7]*jac[7]+jac[8]*jac[8]-1.0)*
        (jac[6]*jac[6]+jac[7]*jac[7]+jac[8]*jac[8]-1.0))
        + 0.5 * (
            (jac[0]*jac[6]+jac[1]*jac[7]+jac[2]*jac[8])*
            (jac[0]*jac[6]+jac[1]*jac[7]+jac[2]*jac[8])+
            (jac[3]*jac[6]+jac[4]*jac[7]+jac[5]*jac[8])*
            (jac[3]*jac[6]+jac[4]*jac[7]+jac[5]*jac[8])+
            (jac[0]*jac[3]+jac[1]*jac[4]+jac[3]*jac[5])*
            (jac[0]*jac[3]+jac[1]*jac[4]+jac[3]*jac[5]));
}

template<int DIM>
NekDouble NodeOpti::GetFunctional()
{
    const int nElmt      = data.size();
    const int totpts = derivUtil->ptsLow * nElmt;
    NekDouble X[DIM * totpts];

    // Store x/y components of each element sequentially in memory
    for (int i = 0, cnt = 0; i < nElmt; ++i)
    {
        for (int j = 0; j < derivUtil->ptsLow; ++j)
        {
            for (int d = 0; d < DIM; ++d)
            {
                X[cnt + d*derivUtil->ptsLow + j] = *(data[i]->nodes[j][d]);
            }
        }

        cnt += DIM*derivUtil->ptsLow;
    }

    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - number of elements
    //   - cartesian coordinate direction
    //   - quadrature points
    NekDouble deriv[DIM][nElmt][DIM][derivUtil->ptsHigh];

    // Calculate x- and y-gradients
    for (int d = 0; d < DIM; ++d)
    {
        Blas::Dgemm('N', 'N', derivUtil->ptsHigh, DIM * nElmt, derivUtil->ptsLow, 1.0,
                    derivUtil->VdmD[d].GetRawPtr(), derivUtil->ptsHigh, X, derivUtil->ptsLow, 0.0,
                    &deriv[d][0][0][0], derivUtil->ptsHigh);
    }

    NekDouble integral = 0.0;

    switch(opti)
    {
        case eLinEl:
        {
            const NekDouble nu = 0.4;
            const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);
            NekDouble jacDet[nElmt][derivUtil->ptsHigh],
                       trEtE[nElmt][derivUtil->ptsHigh];
            NekDouble jacMin = 0.0;

            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    NekDouble jacIdeal[DIM*DIM];
                    int cnt = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n, ++cnt)
                        {
                            jacIdeal[cnt] = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacIdeal[cnt] +=
                                    deriv[l][i][n][k] * data[i]->maps[k][m * 3 + l];
                            }
                        }
                    }
                    jacDet[i][k] = JacDet<DIM>(jacIdeal);
                    trEtE[i][k]  = LinElasTrace<DIM>(jacIdeal);
                    jacMin = min(jacMin,jacDet[i][k]);
                }
            }
            NekDouble gam = numeric_limits<float>::epsilon();
            NekDouble ep;
            if(jacMin < gam)
            {
                ep = sqrt(gam*(gam-jacMin));
            }
            else
            {
                ep = 0.0;
            }

            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    NekDouble sigma = 0.5*(jacDet[i][k] + sqrt(jacDet[i][k]*jacDet[i][k] + 4.0*ep*ep));
                    NekDouble lsigma = log(sigma);
                    integral += derivUtil->quadW[k] * fabs(data[i]->maps[k][9]) * (K * 0.5 * lsigma * lsigma + mu * trEtE[i][k]);
                }
            }
            break;
        }


        case eHypEl:
        {
            const NekDouble nu = 0.4;
            const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
            const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

            for (int i = 0; i < nElmt; ++i)
            {
                bool valid = true;
                NekDouble jacDet[derivUtil->ptsHigh], I1[derivUtil->ptsHigh];

                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    NekDouble jacIdeal[DIM*DIM];
                    int cnt = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n, ++cnt)
                        {
                            jacIdeal[cnt] = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacIdeal[cnt] +=
                                    deriv[l][i][n][k] * data[i]->maps[k][m * 3 + l];
                            }
                        }
                    }

                    jacDet[k] = JacDet<DIM>(jacIdeal);
                    I1[k] = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        I1[k] += jacIdeal[m]*jacIdeal[m];
                    }
                    valid = valid ? jacDet[k] > 0.0 : false;
                }

                if (!valid)
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        NekDouble de = 1e-1;
                        NekDouble sigma = 0.5*(jacDet[k] + sqrt(jacDet[k]*jacDet[k] + 4.0*de*de));
                        NekDouble lsigma = log(sigma);
                        integral += derivUtil->quadW[k]*fabs(data[i]->maps[k][9]) * (0.5 * mu * (I1[k] - 3.0 - 2.0*lsigma) + 0.5 * K * lsigma * lsigma);
                    }
                }
                else
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        NekDouble lsigma = log(jacDet[k]);
                        integral += derivUtil->quadW[k]*fabs(data[i]->maps[k][9]) * (0.5 * mu * (I1[k] - 3.0 - 2.0*lsigma) + 0.5 * K * lsigma * lsigma);
                    }
                }
            }
            break;
        }

        case eRoca:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                bool valid = true;
                NekDouble jacDet[derivUtil->ptsHigh], frob[derivUtil->ptsHigh];

                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    NekDouble jacIdeal[DIM*DIM];
                    int cnt = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n, ++cnt)
                        {
                            jacIdeal[cnt] = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacIdeal[cnt] +=
                                    deriv[l][i][n][k] * data[i]->maps[k][m * 3 + l];
                            }
                        }
                    }

                    frob[k] = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        frob[k] += jacIdeal[m] * jacIdeal[m];
                    }
                    jacDet[k] = JacDet<DIM>(jacIdeal);
                    valid = valid ? jacDet[k] > 0.0 : false;
                }

                if (!valid)
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        NekDouble de = 1e-2;
                        NekDouble sigma = 0.5*(jacDet[k] + sqrt(jacDet[k]*jacDet[k] + 4.0*de*de));
                        integral += derivUtil->quadW[k] * fabs(data[i]->maps[k][9]) * (frob[k] / DIM / pow(fabs(sigma), 2.0/DIM) -1.0);
                    }
                }
                else
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        integral += derivUtil->quadW[k] * fabs(data[i]->maps[k][9]) * (frob[k] / DIM / pow(fabs(jacDet[k]), 2.0/DIM) -1.0);
                    }
                }
            }
            break;
        }

        case eWins:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                bool valid = true;
                NekDouble jacDet[derivUtil->ptsHigh], frob[derivUtil->ptsHigh];

                for(int k = 0; k < derivUtil->ptsHigh; ++k)
                {
                    NekDouble jacIdeal[DIM*DIM];
                    int cnt = 0;
                    for (int m = 0; m < DIM; ++m)
                    {
                        for (int n = 0; n < DIM; ++n, ++cnt)
                        {
                            jacIdeal[cnt] = 0.0;
                            for (int l = 0; l < DIM; ++l)
                            {
                                jacIdeal[cnt] +=
                                    deriv[l][i][n][k] * data[i]->maps[k][m * 3 + l];
                            }
                        }
                    }

                    frob[k] = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        frob[k] += jacIdeal[m] * jacIdeal[m];
                    }
                    jacDet[k] = JacDet<DIM>(jacIdeal);
                    valid = valid ? jacDet[k] > 0.0 : false;
                }

                if (!valid)
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        NekDouble de = 1e-2;
                        NekDouble sigma = 0.5*(jacDet[k] + sqrt(jacDet[k]*jacDet[k] + 4.0*de*de));
                        integral += derivUtil->quadW[k]*fabs(data[i]->maps[k][9])*(frob[k] / sigma);
                    }
                }
                else
                {
                    for(int k = 0; k < derivUtil->ptsHigh; ++k)
                    {
                        integral += derivUtil->quadW[k]*fabs(data[i]->maps[k][9])*(frob[k] / jacDet[k]);
                    }
                }
            }
            break;
        }
    }

    return integral;
}

NodeOptiJob* NodeOpti::GetJob()
{
    return new NodeOptiJob(this);
}

}
}
