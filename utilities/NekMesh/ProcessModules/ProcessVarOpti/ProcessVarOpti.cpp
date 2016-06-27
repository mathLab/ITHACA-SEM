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
#include "NodeOpti.h"
#include "NodeOptiJob.h"

#include <boost/thread/mutex.hpp>

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/NodalUtil.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

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

    derivUtil = boost::shared_ptr<DerivUtil>(new DerivUtil);
    ptsHelp = boost::shared_ptr<PtsHelper>(new PtsHelper);

    FillQuadPoints();

    //build Vandermonde information
    switch (m_mesh->m_spaceDim)
    {
        case 2:
        {
            ptsHelp->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)/2;

            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTriElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTriSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2;

            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2);
            NekVector<NekDouble> U1(u1), V1(v1);
            NekVector<NekDouble> U2(u2), V2(v2);
            ptsHelp->ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            NekMatrix<NekDouble> interp = LibUtilities::GetInterpolationMatrix(U1, V1, U2, V2);

            NekMatrix<NekDouble> Vandermonde = LibUtilities::GetVandermonde(U1,V1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            derivUtil->VdmD[0] = interp * (
              LibUtilities::GetVandermondeForXDerivative(U1,V1) * VandermondeI);
            derivUtil->VdmD[1] = interp * (
              LibUtilities::GetVandermondeForYDerivative(U1,V1) * VandermondeI);
            //derivUtil->quadW = LibUtilities::MakeQuadratureWeights(U2,V1);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            derivUtil->quadW = quadWi;
        }
        break;
        case 3:
        {
            ptsHelp->ptsLow  = m_mesh->m_nummode*(m_mesh->m_nummode+1)*(m_mesh->m_nummode+2)/6;
            LibUtilities::PointsKey pkey1(m_mesh->m_nummode,
                                          LibUtilities::eNodalTetElec);
            LibUtilities::PointsKey pkey2(m_mesh->m_nummode+2,
                                          LibUtilities::eNodalTetSPI);
            Array<OneD, NekDouble> u1, v1, u2, v2, w1, w2;
            LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);
            LibUtilities::PointsManager()[pkey2]->GetPoints(u2, v2, w2);
            NekVector<NekDouble> U1(u1), V1(v1), W1(w1);
            NekVector<NekDouble> U2(u2), V2(v2), W2(w2);
            ptsHelp->ptsHigh = LibUtilities::PointsManager()[pkey2]->GetNumPointsAlt();

            NekMatrix<NekDouble> interp =
                        LibUtilities::GetTetInterpolationMatrix(U1, V1, W1,
                                                                U2, V2, W2);

            NekMatrix<NekDouble> Vandermonde =
                                LibUtilities::GetTetVandermonde(U1,V1,W1);
            NekMatrix<NekDouble> VandermondeI = Vandermonde;
            VandermondeI.Invert();
            derivUtil->VdmD[0] = interp * (
              LibUtilities::GetVandermondeForTetXDerivative(U1,V1,W1) * VandermondeI);
            derivUtil->VdmD[1] = interp * (
              LibUtilities::GetVandermondeForTetYDerivative(U1,V1,W1) * VandermondeI);
            derivUtil->VdmD[2] = interp * (
              LibUtilities::GetVandermondeForTetZDerivative(U1,V1,W1) * VandermondeI);
            Array<OneD, NekDouble> qds = LibUtilities::PointsManager()[pkey2]->GetW();
            NekVector<NekDouble> quadWi(qds);
            derivUtil->quadW = quadWi;
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
            ns.push_back(NodeOpti(freenodes[i][j],it->second,res,derivUtil,ptsHelp,opti));
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
         << "Residual tolerance:\t" << restol << endl;

    int nThreads = m_config["numthreads"].as<int>();

    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm =
                tms.CreateInstance(Thread::ThreadMaster::SessionJob, nThreads);

    Timer t;
    t.Start();

    while (res->val > restol)
    {
        ctr++;
        res->val = 0.0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob*> jobs(optiNodes[i].size());
            for(int j = 0; j < optiNodes[i].size(); j++)
            {
                jobs[j] = optiNodes[i][j].GetJob();
            }

            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(nThreads);
            tm->Wait();
        }

        EvaluateMesh();

        cout << ctr <<  "\tResidual: " << res->val << " " << res->worstJac << endl;
        if(ctr >= maxIter)
            break;
    }

    t.Stop();
    cout << "Time to compute: " << t.TimePerTest(1) << endl;

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

}
}
