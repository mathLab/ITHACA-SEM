////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJac.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessVarOpti.h"

#include <boost/thread/mutex.hpp>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

#include <LibUtilities/Foundations/NodalUtil.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx;
int ptsLow;
int ptsHigh;
NekMatrix<NekDouble> VdmD[3];
NekVector<NekDouble> quadW;
optimiser opti;

ModuleKey ProcessVarOpti::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varopti"),
    ProcessVarOpti::create,
    "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["linearelastic"] =
        ConfigOption(true, "", "Optimise for linear elasticity");
    m_config["winslow"] =
        ConfigOption(true, "", "Optimise for winslow");
    m_config["roca"] =
        ConfigOption(true, "", "Optimise for roca method");
    m_config["hyperelastic"] =
        ConfigOption(true, "", "Optimise for hyper elasticity");
    m_config["numthreads"] =
        ConfigOption(false, "1", "Number of threads");
    m_config["nq"] =
        ConfigOption(false, "0", "Number of quad points");
    m_config["stats"] =
        ConfigOption(false, "", "Write a file with list of scaled jacobians");
    m_config["restol"] =
        ConfigOption(false, "1e-6", "Tolerance criterion");
    m_config["maxiter"] =
        ConfigOption(false, "500", "Maximum number of iterations");
}

ProcessVarOpti::~ProcessVarOpti()
{
}

void ProcessVarOpti::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessVarOpti: Optimising... " << endl;
    }

    if(m_config["linearelastic"].beenSet)
    {
        opti = eLinEl;
    }
    else if(m_config["winslow"].beenSet)
    {
        opti = eWins;
    }
    else if(m_config["roca"].beenSet)
    {
        opti = eRoca;
    }
    else if(m_config["hyperelastic"].beenSet)
    {
        opti = eHypEl;
    }
    else
    {
        ASSERTL0(false,"not opti type set");
    }

    const int maxIter = m_config["maxiter"].as<int>();
    const NekDouble restol = m_config["restol"].as<NekDouble>();

    m_mesh->m_nummode = m_config["nq"].as<int>();

    ASSERTL0(m_mesh->m_nummode > 2,"not specified high-order");

    if(m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"cannot deal with manifolds");
    }

    res = boost::shared_ptr<Residual>(new Residual);
    res->val = 1.0;

    FillQuadPoints();

    //build Vandermonde information
    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;

            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTriElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTriSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2;

            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);
            NekVector<NekDouble> U1(u1), V1(v1);
            NekVector<NekDouble> U2(u2), V2(v2);
            ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            NekMatrix<NekDouble> interp = LibUtilities::GetInterpolationMatrix(U1, V1, U2, V2);

            NekMatrix<NekDouble> Vandermonde = LibUtilities::GetVandermonde(U1,V1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            VdmD[0] = interp * (
              LibUtilities::GetVandermondeForXDerivative(U1,V1) * VandermondeI);
            VdmD[1] = interp * (
              LibUtilities::GetVandermondeForYDerivative(U1,V1) * VandermondeI);
            //quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            quadW = quadWi;
        }
        break;
        case 3:
        {
            ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)*(m_mesh->m_nummode+2)/6;
            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTetElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTetSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);
            NekVector<NekDouble> U1(u1), V1(v1), W1(w1);
            NekVector<NekDouble> U2(u2), V2(v2), W2(w2);
            ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            NekMatrix<NekDouble> interp =
                        LibUtilities::GetTetInterpolationMatrix(U1, V1, W1,
                                                                U2, V2, W2);

            NekMatrix<NekDouble> Vandermonde =
                                LibUtilities::GetTetVandermonde(U1,V1,W1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            VdmD[0] = interp * (
              LibUtilities::GetVandermondeForTetXDerivative(U1,V1,W1) * VandermondeI);
            VdmD[1] = interp * (
              LibUtilities::GetVandermondeForTetYDerivative(U1,V1,W1) * VandermondeI);
            VdmD[2] = interp * (
              LibUtilities::GetVandermondeForTetZDerivative(U1,V1,W1) * VandermondeI);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            quadW = quadWi;
        }
    }

    GetElementMap();

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes();
    vector<vector<NodeOpti> > optiNodes;

    for(int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOpti> ns;
        for(int j = 0; j < freenodes[i].size(); j++)
        {
            NodeElMap::iterator it = nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != nodeElMap.end(), "could not find");
            ns.push_back(NodeOpti(freenodes[i][j],it->second,res));
        }
        optiNodes.push_back(ns);
    }

    int nset = optiNodes.size();
    int p = 0;
    int mn = numeric_limits<int>::max();
    int mx = 0;
    for(int i = 0; i < nset; i++)
    {
        p += optiNodes[i].size();
        mn = min(mn, int(optiNodes[i].size()));
        mx = max(mx, int(optiNodes[i].size()));
    }

    cout << scientific << endl;
    cout << "N elements:\t\t" << m_mesh->m_element[m_mesh->m_expDim].size() << endl
         << "N elements invalid:\t" << res->startInv << endl
         << "Worst jacobian:\t\t" << res->worstJac << endl
         << "N free nodes:\t\t" << res->n << endl
         << "N Dof:\t\t\t" << res->nDoF << endl
         << "N Dirclet:\t\t" << res->nDirc << endl
         << "N color sets:\t\t" << nset << endl
         << "Avg set colors:\t\t" << p/nset << endl
         << "Min set:\t\t" << mn << endl
         << "Max set:\t\t" << mx << endl
         << "Residual tolerance:\t\t" << restol << endl;

    int nThreads = m_config["numthreads"].as<int>();

    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm =
                tms.CreateInstance(Thread::ThreadMaster::SessionJob, nThreads);

    while (res->val > restol)
    {
        ctr++;
        res->val = 0.0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob*> jobs(optiNodes[i].size());
            if (m_mesh->m_spaceDim == 2)
            {
                for(int j = 0; j < optiNodes[i].size(); j++)
                {
                    jobs[j] = optiNodes[i][j].GetJob<2>();
                }
            }
            else
            {
                for(int j = 0; j < optiNodes[i].size(); j++)
                {
                    jobs[j] = optiNodes[i][j].GetJob<3>();
                }
            }

            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(nThreads);
            tm->Wait();
        }

        cout << ctr <<  "\tResidual: " << res->val << endl;
        if(ctr > maxIter)
            break;
    }

    EvaluateMesh();

    cout << "Invalid at end:\t\t" << res->startInv << endl;
    cout << "Worst at end:\t\t" << res->worstJac << endl;

    if(m_config["stats"].beenSet)
    {
        string file = m_config["stats"].as<string>();
        cout << "writing stats to " << file.c_str() << endl;
        WriteStats(file);
    }
}

template<int DIM>
void ProcessVarOpti::NodeOpti::Optimise()
{
    //it doesnt matter at this point what the dim is so long that
    //in the 2d case z is left as zero

    Array<OneD, NekDouble> G = GetGrad<DIM>();

    if(sqrt(G[0]*G[0] + G[1]*G[1] + G[2]*G[2]) > 1e-10)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional<DIM>();
        NekDouble xc       = node->m_x;
        NekDouble yc       = node->m_y;
        NekDouble zc       = node->m_z;
        NekDouble alpha    = 1.0;
        NekDouble delX     = 0.0;
        NekDouble delY     = 0.0;
        NekDouble delZ     = 0.0;

        if(DIM == 2)
        {
             delX = 1.0/(G[3]*G[4]-G[6]*G[6])*(G[4]*G[0] - G[6]*G[1]);
             delY = 1.0/(G[3]*G[4]-G[6]*G[6])*(G[3]*G[1] - G[6]*G[0]);
        }
        else
        {
            DNekMat H(3,3,0.0);
            H(0,0) = G[3];
            H(1,1) = G[4];
            H(2,2) = G[5];
            H(0,1) = G[6];
            H(1,0) = G[6];
            H(0,2) = G[7];
            H(2,0) = G[7];
            H(2,1) = G[8];
            H(1,2) = G[8];
            H.Invert();
            NekVector<NekDouble> g(3);
            g[0] = G[0];
            g[1] = G[1];
            g[2] = G[2];
            NekVector<NekDouble> del = H * g;
            delX = del[0];
            delY = del[1];
            delZ = del[2];
        }

        bool found = false;
        while(alpha > 1e-10)
        {
            node->m_x = xc - alpha * delX;
            node->m_y = yc - alpha * delY;
            node->m_z = zc - alpha * delZ;
            if(GetFunctional<DIM>() < currentW)
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
            // cout << "warning: had to reset node" << endl;
            // cout << G[0] << " " << G[1] << " " << G[2] << " " << node->m_id << endl;
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

template<int DIM>
Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad()
{
    return Array<OneD, NekDouble>();
}

template<>
Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad<2>()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble dx = 1e-4;
    vector<NekDouble> w(9);

    for(int i = 0; i < 7; i++)
    {
        node->m_x = xc + dir[i][0] * dx;
        node->m_y = yc + dir[i][1] * dx;
        w[i] = GetFunctional<2>();
    }
    node->m_x = xc;
    node->m_y = yc;

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
    ret[3] = (w[1] + w[4] - 2.0*w[0]) / dx / dx;
    ret[4] = (w[3] + w[6] - 2.0*w[0]) / dx / dx;
    ret[6] = (w[2] - w[1] - w[3] + 2.0*w[0] - w[4] - w[6] + w[5]) / 2.0 / dx / dx;

    return ret;
}

template<>
Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad<3>()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble zc = node->m_z;
    NekDouble dx = 1e-6;

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
NekDouble ProcessVarOpti::NodeOpti::GetFunctional()
{
    const int nElmt      = data.size();
    const int totPtsLow  = ptsLow * nElmt;
    NekDouble X[DIM * totPtsLow];

    // Store x/y components of each element sequentially in memory
    for (int i = 0, cnt = 0; i < nElmt; ++i)
    {
        for (int j = 0; j < ptsLow; ++j)
        {
            for (int d = 0; d < DIM; ++d)
            {
                X[cnt + d*ptsLow + j] = *(data[i]->nodes[j][d]);
            }
        }

        cnt += DIM*ptsLow;
    }

    // Storage for derivatives, ordered by:
    //   - standard coordinate direction
    //   - number of elements
    //   - cartesian coordinate direction
    //   - quadrature points
    NekDouble deriv[DIM][nElmt][DIM][ptsHigh];

    // Calculate x- and y-gradients
    for (int d = 0; d < DIM; ++d)
    {
        Blas::Dgemm('N', 'N', ptsHigh, DIM * nElmt, ptsLow, 1.0,
                    VdmD[d].GetRawPtr(), ptsHigh, X, ptsLow, 0.0,
                    &deriv[d][0][0][0], ptsHigh);
    }

    NekDouble integral = 0.0;

    switch(opti)
    {
        case eLinEl:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < ptsHigh; ++k)
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
                    NekDouble jacDet = JacDet<DIM>(jacIdeal);
                    NekDouble trEtE = LinElasTrace<DIM>(jacIdeal);

                    const NekDouble nu = 0.45;
                    const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
                    const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

                    NekDouble de = fabs(data[i]->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    NekDouble lsigma = log(sigma);
                    integral += quadW[k] * (K * 0.5 * lsigma * lsigma + mu * trEtE);
                }
            }
            break;
        }

        case eHypEl:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < ptsHigh; ++k)
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

                    NekDouble I1 = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        I1 += jacIdeal[m]*jacIdeal[m];
                    }

                    const NekDouble nu = 0.45;
                    const NekDouble mu = 1.0 / 2.0 / (1.0+nu);
                    const NekDouble K  = 1.0 / 3.0 / (1.0 - 2.0 * nu);

                    NekDouble jacDet = JacDet<DIM>(jacIdeal);
                    NekDouble de = fabs(data[i]->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    NekDouble lsigma = log(sigma);
                    integral += quadW[k]*(0.5 * mu * (I1 - 3.0 - 2.0*lsigma) + 0.5 * K * lsigma * lsigma);
                }
            }
            break;
        }

        case eRoca:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < ptsHigh; ++k)
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

                    NekDouble frob = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        frob += jacIdeal[m] * jacIdeal[m];
                    }
                    NekDouble jacDet = JacDet<DIM>(jacIdeal);

                    NekDouble de = fabs(data[i]->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    integral += quadW[k]*(frob / DIM / pow(fabs(sigma), 2.0/DIM) -1.0);
                }
            }
            break;
        }

        case eWins:
        {
            for (int i = 0; i < nElmt; ++i)
            {
                for(int k = 0; k < ptsHigh; ++k)
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

                    NekDouble frob = 0.0;
                    for (int m = 0; m < DIM*DIM; ++m)
                    {
                        frob += jacIdeal[m] * jacIdeal[m];
                    }
                    NekDouble jacDet = JacDet<DIM>(jacIdeal);

                    NekDouble de = fabs(data[i]->maps[k][9]) * sqrt(1E-3*1E-3 + 1E-3);
                    NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet + 4.0*de*de));
                    integral += quadW[k]*(frob / sigma);
                }
            }
            break;
        }
    }

    return integral;
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes()
{
    //this figures out the dirclet nodes and colors the others into paralell sets
    NodeSet boundaryNodes;

    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            EdgeSet::iterator it;
            for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                boundaryNodes.insert((*it)->m_n1);
                boundaryNodes.insert((*it)->m_n2);
                for(int i = 0; i < (*it)->m_edgeNodes.size(); i++)
                {
                    boundaryNodes.insert((*it)->m_edgeNodes[i]);
                }
            }
            break;
        }
        case 3:
        {
            FaceSet::iterator it;
            for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
            {
                if((*it)->m_elLink.size() == 2)
                {
                    continue;
                }

                vector<NodeSharedPtr> vs = (*it)->m_vertexList;
                for(int j = 0; j < vs.size(); j++)
                {
                    boundaryNodes.insert(vs[j]);
                }

                vector<EdgeSharedPtr> es = (*it)->m_edgeList;
                for(int j = 0; j < es.size(); j++)
                {
                    for(int k = 0; k < es[j]->m_edgeNodes.size(); k++)
                    {
                        boundaryNodes.insert(es[j]->m_edgeNodes[k]);
                    }
                }

                for(int i = 0; i < (*it)->m_faceNodes.size(); i++)
                {
                    boundaryNodes.insert((*it)->m_faceNodes[i]);
                }
            }
            break;
        }
        default:
            ASSERTL0(false,"space dim issue");
    }

    res->nDirc = boundaryNodes.size();

    vector<NodeSharedPtr> remain;

    EdgeSet::iterator eit;
    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        vector<NodeSharedPtr> n = (*eit)->m_edgeNodes;
        n.push_back((*eit)->m_n1);
        n.push_back((*eit)->m_n2);
        for(int j = 0; j < n.size(); j++)
        {
            NodeSet::iterator nit = boundaryNodes.find(n[j]);
            if(nit == boundaryNodes.end())
            {
                remain.push_back(n[j]);
            }
        }
    }

    if(m_mesh->m_expDim == 2)
    {

        for(int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVolumeNodes();
            for(int j = 0; j < ns.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(ns[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(ns[j]);
                }
            }
        }
    }
    else
    {
        FaceSet::iterator fit;
        for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
        {
            for(int j = 0; j < (*fit)->m_faceNodes.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find((*fit)->m_faceNodes[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back((*fit)->m_faceNodes[j]);
                }
            }
        }

        for(int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVolumeNodes();
            for(int j = 0; j < ns.size(); j++)
            {
                NodeSet::iterator nit = boundaryNodes.find(ns[j]);
                if(nit == boundaryNodes.end())
                {
                    remain.push_back(ns[j]);
                }
            }
        }
    }

    res->n = remain.size();
    res->nDoF = res->n * m_mesh->m_spaceDim;

    vector<vector<NodeSharedPtr> > ret;

    while (remain.size() > 0)
    {
        vector<NodeSharedPtr> layer;
        set<int> locked;
        set<int> completed;
        for(int i = 0; i < remain.size(); i++)
        {
            NodeElMap::iterator it = nodeElMap.find(remain[i]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            bool islocked = false;
            for(int j = 0; j < it->second.size(); j++)
            {
                set<int>::iterator sit = locked.find(it->second[j]->el->GetId());
                if(sit != locked.end())
                {
                    islocked = true;
                    break;
                }
            }
            if(!islocked)
            {
                layer.push_back(remain[i]);
                completed.insert(remain[i]->m_id);
                for(int j = 0; j < it->second.size(); j++)
                {
                    locked.insert(it->second[j]->el->GetId());
                }
            }
        }

        vector<NodeSharedPtr> tmp = remain;
        remain.clear();
        for(int i = 0; i < tmp.size(); i++)
        {
            set<int>::iterator sit = completed.find(tmp[i]->m_id);
            if(sit == completed.end())
            {
                remain.push_back(tmp[i]);
            }
        }
        ret.push_back(layer);
    }
    return ret;
}

void ProcessVarOpti::GetElementMap()
{
    //build ideal maps and structs;
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);
        ElDataSharedPtr d = boost::shared_ptr<ElData>(new ElData(ns, m_mesh->m_spaceDim));
        d->el   = el;
        d->maps = MappingIdealToRef(el);
        dataSet.push_back(d);
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
        vector<NodeSharedPtr> ns;
        el->GetCurvedNodes(ns);

        for(int j = 0; j < ns.size(); j++)
        {
            nodeElMap[ns[j]->m_id].push_back(dataSet[i]);
        }
    }
}

vector<Array<OneD, NekDouble> > ProcessVarOpti::MappingIdealToRef(ElementSharedPtr el)
{
    //need to make ideal element out of old element
    ElmtConfig ec = el->GetConf();
    ec.m_order  = 1;
    ec.m_faceNodes = false;
    ec.m_volumeNodes = false;
    ec.m_reorient = false;

    ElementSharedPtr E = GetElementFactory().CreateInstance(
                            ec.m_e, ec, el->GetVertexList(),
                            el->GetTagList());

    SpatialDomains::GeometrySharedPtr    geom = E->GetGeom(el->GetDim());
    geom->FillGeom();
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

    vector<Array<OneD, NekDouble> > ret;

    if(geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        ASSERTL0(false,"Not coded");
        /*vector<Array<OneD, NekDouble> > xy;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(2);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xy.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();

        for(int j = 0; j < b[1]->GetNumPoints(); j++)
        {
            for(int i = 0; i < b[0]->GetNumPoints(); i++)
            {
                NekDouble a1 = 0.5*(1.0-u[i]), a2 = 0.5*(1.0+u[i]);
                NekDouble b1 = 0.5*(1.0-v[j]), b2 = 0.5*(1.0+v[j]);
                DNekMat dxdz(2,2,1.0,eFULL);

                dxdz(0,0) = 0.5*(-b1*xy[0][0] + b1*xy[1][0] + b2*xy[2][0] - b2*xy[3][0]);
                dxdz(1,0) = 0.5*(-b1*xy[0][1] + b1*xy[1][1] + b2*xy[2][1] - b2*xy[3][1]);

                dxdz(0,1) = 0.5*(-a1*xy[0][0] - a2*xy[1][0] + a2*xy[2][0] + a1*xy[3][0]);
                dxdz(1,1) = 0.5*(-a1*xy[0][1] - a2*xy[1][1] + a2*xy[2][1] + a1*xy[3][1]);

                NekDouble det = 1.0/(dxdz(0,0)*dxdz(1,1) - dxdz(1,0)*dxdz(0,1));

                dxdz.Invert();
                Array<OneD, NekDouble> r(9,0.0);
                r[0] = dxdz(0,0);
                r[1] = dxdz(1,0);
                r[3] = dxdz(0,1);
                r[4] = dxdz(1,1);
                ret.push_back(r);
            }
        }*/
    }
    else if(geom->GetShapeType() == LibUtilities::eTriangle)
    {
        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u, v;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);

        NekVector<NekDouble> X(ptsLow),Y(ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(2);
            xp[0] = u[j];
            xp[1] = v[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
        }

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh);

        x1i = VdmD[0]*X;
        y1i = VdmD[0]*Y;
        x2i = VdmD[1]*X;
        y2i = VdmD[1]*Y;

        for(int i = 0 ; i < ptsHigh; i++)
        {
            DNekMat dxdz(2,2,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);

            Array<OneD, NekDouble> r(10,0.0);
            r[9] = dxdz(0,0)*dxdz(1,1)-dxdz(1,0)*dxdz(0,1);

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                     LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u, v, w;
        LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

        Array<OneD, NekDouble> xc(chi->GetTotPoints());
        Array<OneD, NekDouble> yc(chi->GetTotPoints());
        Array<OneD, NekDouble> zc(chi->GetTotPoints());

        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        chi->BwdTrans(coeffs0,xc);
        chi->BwdTrans(coeffs1,yc);
        chi->BwdTrans(coeffs2,zc);

        NekVector<NekDouble> X(ptsLow),Y(ptsLow),Z(ptsLow);
        for(int j = 0; j < u.num_elements(); j++)
        {
            Array<OneD, NekDouble> xp(3);
            xp[0] = u[j];
            xp[1] = v[j];
            xp[2] = w[j];

            X(j) = chi->PhysEvaluate(xp, xc);
            Y(j) = chi->PhysEvaluate(xp, yc);
            Z(j) = chi->PhysEvaluate(xp, zc);
        }

        NekVector<NekDouble> x1i(ptsHigh),y1i(ptsHigh),z1i(ptsHigh),
                             x2i(ptsHigh),y2i(ptsHigh),z2i(ptsHigh),
                             x3i(ptsHigh),y3i(ptsHigh),z3i(ptsHigh);

        x1i = VdmD[0]*X;
        y1i = VdmD[0]*Y;
        x2i = VdmD[1]*X;
        y2i = VdmD[1]*Y;
        z1i = VdmD[0]*Z;
        z2i = VdmD[1]*Z;
        x3i = VdmD[2]*X;
        y3i = VdmD[2]*Y;
        z3i = VdmD[2]*Z;

        for(int i = 0 ; i < ptsHigh; i++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i(i);
            dxdz(0,1) = x2i(i);
            dxdz(0,2) = x3i(i);
            dxdz(1,0) = y1i(i);
            dxdz(1,1) = y2i(i);
            dxdz(1,2) = y3i(i);
            dxdz(2,0) = z1i(i);
            dxdz(2,1) = z2i(i);
            dxdz(2,2) = z3i(i);

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                   -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                   +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            dxdz.Invert();

            r[0] = dxdz(0,0);
            r[1] = dxdz(1,0);
            r[2] = dxdz(2,0);
            r[3] = dxdz(0,1);
            r[4] = dxdz(1,1);
            r[5] = dxdz(2,1);
            r[6] = dxdz(0,2);
            r[7] = dxdz(1,2);
            r[8] = dxdz(2,2);
            ret.push_back(r);
        }
    }
    else if(geom->GetShapeType() == LibUtilities::ePrism)
    {
        ASSERTL0(false, "not coded");
        /*vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> eta1 = b[0]->GetZ();
        Array<OneD, NekDouble> eta2 = b[1]->GetZ();
        Array<OneD, NekDouble> eta3 = b[2]->GetZ();

        for(int k = 0; k < b[2]->GetNumPoints(); k++)
        {

            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int i = 0; i < b[0]->GetNumPoints(); i++)
                {
                    NekDouble xi1 = 0.5*(1+eta1[i])*(1-eta3[k])-1.0;
                    NekDouble a1 = 0.5*(1-xi1),     a2 = 0.5*(1+xi1);
                    NekDouble b1 = 0.5*(1-eta2[j]), b2 = 0.5*(1+eta2[j]);
                    NekDouble c1 = 0.5*(1-eta3[k]), c2 = 0.5*(1+eta3[k]);

                    DNekMat dxdz(3,3,1.0,eFULL);

                    dxdz(0,0) = 0.5*(-b1*xyz[0][0] + b1*xyz[1][0] + b2*xyz[2][0] - b2*xyz[3][0]);
                    dxdz(1,0) = 0.5*(-b1*xyz[0][1] + b1*xyz[1][1] + b2*xyz[2][1] - b2*xyz[3][1]);
                    dxdz(2,0) = 0.5*(-b1*xyz[0][2] + b1*xyz[1][2] + b2*xyz[2][2] - b2*xyz[3][2]);

                    dxdz(0,1) = 0.5*((a2-c1)*xyz[0][0] - a2*xyz[1][0] + a2*xyz[2][0] + (c1-a2)*xyz[3][0] - c2*xyz[4][0] + c2*xyz[5][0]);
                    dxdz(1,1) = 0.5*((a2-c1)*xyz[0][1] - a2*xyz[1][1] + a2*xyz[2][1] + (c1-a2)*xyz[3][1] - c2*xyz[4][1] + c2*xyz[5][1]);
                    dxdz(2,1) = 0.5*((a2-c1)*xyz[0][2] - a2*xyz[1][2] + a2*xyz[2][2] + (c1-a2)*xyz[3][2] - c2*xyz[4][2] + c2*xyz[5][2]);

                    dxdz(0,2) = 0.5*(-b1*xyz[0][0] - b2*xyz[3][0] + b1*xyz[4][0] + b2*xyz[5][0]);
                    dxdz(1,2) = 0.5*(-b1*xyz[0][1] - b2*xyz[3][1] + b1*xyz[4][1] + b2*xyz[5][1]);
                    dxdz(2,2) = 0.5*(-b1*xyz[0][2] - b2*xyz[3][2] + b1*xyz[4][2] + b2*xyz[5][2]);

                    dxdz.Invert();
                    Array<OneD, NekDouble> r(9,0.0);
                    r[0] = dxdz(0,0);
                    r[1] = dxdz(1,0);
                    r[3] = dxdz(0,1);
                    r[4] = dxdz(1,1);
                    r[2] = dxdz(2,0);
                    r[5] = dxdz(2,1);
                    r[6] = dxdz(0,2);
                    r[7] = dxdz(1,2);
                    r[8] = dxdz(2,2);
                    ret.push_back(r);
                }
            }
        }*/
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

void ProcessVarOpti::FillQuadPoints()
{
    int nq = m_mesh->m_nummode;
    int id = m_mesh->m_vertexSet.size();

    LibUtilities::PointsKey ekey(m_mesh->m_nummode,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    EdgeSet::iterator eit;

    for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
    {
        SpatialDomains::Geometry1DSharedPtr geom =
                                            (*eit)->GetGeom(m_mesh->m_spaceDim);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

        vector<NodeSharedPtr> hons;

        switch (m_mesh->m_spaceDim)
        {
            case 2:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),0.0)));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
            case 3:
            {
                Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
                Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
                Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

                Array<OneD, NekDouble> xc(xmap->GetTotPoints());
                Array<OneD, NekDouble> yc(xmap->GetTotPoints());
                Array<OneD, NekDouble> zc(xmap->GetTotPoints());

                xmap->BwdTrans(coeffs0,xc);
                xmap->BwdTrans(coeffs1,yc);
                xmap->BwdTrans(coeffs2,zc);

                for(int j = 1; j < m_mesh->m_nummode - 1; j++)
                {
                    Array<OneD, NekDouble> xp(2);
                    xp[0] = gll[j];

                    hons.push_back(boost::shared_ptr<Node>(new Node(
                            id++,xmap->PhysEvaluate(xp,xc),
                                 xmap->PhysEvaluate(xp,yc),
                                 xmap->PhysEvaluate(xp,zc))));
                }
                (*eit)->m_edgeNodes.clear();
                (*eit)->m_edgeNodes = hons;
                (*eit)->m_curveType = LibUtilities::eGaussLobattoLegendre;
            }
            break;
        }
    }

    if(m_mesh->m_expDim == 2)
    {
        //for faces need to do volume nodes
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom =
                                            el->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);

            vector<NodeSharedPtr> hons;

            for(int j = nq * (nq + 1) / 2 - (nq-2)*(nq-3) / 2;
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),0.0)));
            }

            el->SetVolumeNodes(hons);
            el->SetCurveType(LibUtilities::eNodalTriElec);
        }
    }
    else
    {
        FaceSet::iterator it;
        for(it = m_mesh->m_faceSet.begin(); it != m_mesh->m_faceSet.end(); it++)
        {
            //this is a hack and needs to be fixed
            //it really should take the get geom of the whole element and
            //then pick the correct parts
            /////////
            (*it)->m_faceNodes.clear();
            ////////
            SpatialDomains::Geometry2DSharedPtr geom =
                                            (*it)->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTriElec);
            Array<OneD, NekDouble> u, v;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v);

            vector<NodeSharedPtr> hons;

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            for(int j = nq * (nq + 1) / 2 - (nq-2)*(nq-3) / 2;
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(2);
                xp[0] = u[j];
                xp[1] = v[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }
            (*it)->m_faceNodes.clear();
            (*it)->m_faceNodes = hons;
            (*it)->m_curveType = LibUtilities::eNodalTriElec;
        }
        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            SpatialDomains::GeometrySharedPtr geom =
                                            el->GetGeom(m_mesh->m_spaceDim);
            geom->FillGeom();
            StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();

            LibUtilities::PointsKey pkey(m_mesh->m_nummode,
                                         LibUtilities::eNodalTetElec);
            Array<OneD, NekDouble> u, v, w;
            LibUtilities::PointsManager()[pkey]->GetPoints(u, v, w);

            Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
            Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
            Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

            Array<OneD, NekDouble> xc(xmap->GetTotPoints());
            Array<OneD, NekDouble> yc(xmap->GetTotPoints());
            Array<OneD, NekDouble> zc(xmap->GetTotPoints());

            xmap->BwdTrans(coeffs0,xc);
            xmap->BwdTrans(coeffs1,yc);
            xmap->BwdTrans(coeffs2,zc);

            vector<NodeSharedPtr> hons;

            //need to finish for tet
            for(int j = 4 + 6*(nq-2) + 4 * ((nq-2)*(nq-3) / 2);
                                                    j < u.num_elements(); j++)
            {
                Array<OneD, NekDouble> xp(3);
                xp[0] = u[j];
                xp[1] = v[j];
                xp[2] = w[j];

                hons.push_back(boost::shared_ptr<Node>(new Node(
                        id++,xmap->PhysEvaluate(xp,xc),
                             xmap->PhysEvaluate(xp,yc),
                             xmap->PhysEvaluate(xp,zc))));
            }

            el->SetVolumeNodes(hons);
            el->SetCurveType(LibUtilities::eNodalTetElec);
        }
    }

    EvaluateMesh();
}

void ProcessVarOpti::EvaluateMesh()
{
    res->startInv =0;
    res->worstJac = numeric_limits<double>::max();

    if(m_mesh->m_expDim == 2)
    {
        LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                      LibUtilities::eNodalTriElec);
        Array<OneD, NekDouble> u1, v1;

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);

        NekVector<NekDouble> U1(u1), V1(v1);

        NekMatrix<NekDouble> Vandermonde = LibUtilities::GetVandermonde(U1,V1);
        NekMatrix<NekDouble> VandermondeI = Vandermonde;
        VandermondeI.Invert();
        NekMatrix<NekDouble> VdmDxt =  (
          LibUtilities::GetVandermondeForXDerivative(U1,V1) * VandermondeI);
        NekMatrix<NekDouble> VdmDyt =  (
          LibUtilities::GetVandermondeForYDerivative(U1,V1) * VandermondeI);



        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            vector<NodeSharedPtr> ns;
            el->GetCurvedNodes(ns);

            NekDouble mx = -1.0 * numeric_limits<double>::max();
            NekDouble mn =  numeric_limits<double>::max();


            NekVector<NekDouble> X(u1.num_elements()),Y(u1.num_elements());
            for(int j = 0; j < u1.num_elements(); j++)
            {
                X(j) = ns[j]->m_x;
                Y(j) = ns[j]->m_y;
            }

            NekVector<NekDouble> x1i(u1.num_elements()),y1i(u1.num_elements()),
                                 x2i(u1.num_elements()),y2i(u1.num_elements());

            x1i = VdmDxt*X;
            y1i = VdmDxt*Y;
            x2i = VdmDyt*X;
            y2i = VdmDyt*Y;

            for(int j = 0; j < u1.num_elements(); j++)
            {
                NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
                mx = max(mx,jacDet);
                mn = min(mn,jacDet);
            }

            if(mn < 0)
            {
                res->startInv++;
            }
            res->worstJac = min(res->worstJac,mn/mx);
        }
    }
    else
    {
        LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                      LibUtilities::eNodalTetElec);
        Array<OneD, NekDouble> u1, v1,w1;

        LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1,w1);

        NekVector<NekDouble> U1(u1), V1(v1), W1(w1);

        NekMatrix<NekDouble> Vandermonde = LibUtilities::GetTetVandermonde(U1,V1,W1);
        NekMatrix<NekDouble> VandermondeI = Vandermonde;
        VandermondeI.Invert();
        NekMatrix<NekDouble> VdmDxt =  (
          LibUtilities::GetVandermondeForTetXDerivative(U1,V1,W1) * VandermondeI);
        NekMatrix<NekDouble> VdmDyt =  (
          LibUtilities::GetVandermondeForTetYDerivative(U1,V1,W1) * VandermondeI);
        NekMatrix<NekDouble> VdmDzt =  (
          LibUtilities::GetVandermondeForTetZDerivative(U1,V1,W1) * VandermondeI);

        for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

            vector<NodeSharedPtr> ns;
            el->GetCurvedNodes(ns);

            NekDouble mx = -1.0 * numeric_limits<double>::max();
            NekDouble mn =  numeric_limits<double>::max();

            NekVector<NekDouble> X(u1.num_elements()),Y(u1.num_elements()),Z(u1.num_elements());
            for(int j = 0; j < u1.num_elements(); j++)
            {
                X(j) = ns[j]->m_x;
                Y(j) = ns[j]->m_y;
                Z(j) = ns[j]->m_z;
            }

            NekVector<NekDouble> x1i(u1.num_elements()),y1i(u1.num_elements()),z1i(u1.num_elements()),
                                 x2i(u1.num_elements()),y2i(u1.num_elements()),z2i(u1.num_elements()),
                                 x3i(u1.num_elements()),y3i(u1.num_elements()),z3i(u1.num_elements());

            x1i = VdmDxt*X;
            y1i = VdmDxt*Y;
            z1i = VdmDxt*Z;
            x2i = VdmDyt*X;
            y2i = VdmDyt*Y;
            z2i = VdmDyt*Z;
            x3i = VdmDzt*X;
            y3i = VdmDzt*Y;
            z3i = VdmDzt*Z;

            for(int j = 0; j < u1.num_elements(); j++)
            {

                DNekMat dxdz(3,3,1.0,eFULL);
                dxdz(0,0) = x1i(j);
                dxdz(0,1) = x2i(j);
                dxdz(0,2) = x3i(j);
                dxdz(1,0) = y1i(j);
                dxdz(1,1) = y2i(j);
                dxdz(1,2) = y3i(j);
                dxdz(2,0) = z1i(j);
                dxdz(2,1) = z2i(j);
                dxdz(2,2) = z3i(j);

                NekDouble jacDet = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                       -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                       +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

                mx = max(mx,jacDet);
                mn = min(mn,jacDet);
            }



            if(mn < 0)
            {
                res->startInv++;
            }

            res->worstJac = min(res->worstJac,mn/mx);
        }
    }
}

void ProcessVarOpti::WriteStats(string file)
{
    ASSERTL0(file != "", "no file name given");

    ofstream out;
    out.open(file.c_str());
    out << scientific;

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> ns;
        m_mesh->m_element[m_mesh->m_expDim][i]->GetCurvedNodes(ns);
        int pts = ns.size();
        int dim = m_mesh->m_element[m_mesh->m_expDim][i]->GetDim();

        NekDouble mx = -1.0 * numeric_limits<double>::max();
        NekDouble mn = numeric_limits<double>::max();

        if(dim == 2)
        {
            NekVector<NekDouble> X(pts),Y(pts),Z(pts),
                                 x1(pts),y1(pts),
                                 x2(pts),y2(pts);
            for(int i = 0; i < pts; i++)
            {
                X(i) = ns[i]->m_x;
                Y(i) = ns[i]->m_y;
            }

            x1 = VdmD[0]*X;
            y1 = VdmD[0]*Y;
            x2 = VdmD[0]*X;
            y2 = VdmD[0]*Y;

            for(int i = 0; i < pts; i++)
            {
                mx = max(mx, x1(i)*y2(i)-y1(i)*x2(i));
                mn = min(mn, x1(i)*y2(i)-y1(i)*x2(i));
            }

        }
        else
        {
            ASSERTL0(false,"not coded");
        }
        out << mn / mx << endl;
    }
}

}
}
