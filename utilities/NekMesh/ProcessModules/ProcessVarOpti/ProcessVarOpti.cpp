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

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include "ElUtil.h"
#include "NodeOpti.h"
#include "ProcessVarOpti.h"
#include <NekMeshUtils/MeshElements/Element.h>

#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdTriExp.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

// Including Timer.h includes Windows.h, which causes GetJob to be set as a
// macro for some reason.
#if _WIN32
#undef GetJob
#endif

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessVarOpti::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "varopti"), ProcessVarOpti::create,
        "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{
    // clang-format off
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
    m_config["restol"] =
        ConfigOption(false, "1e-6", "Tolerance criterion");
    m_config["maxiter"] =
        ConfigOption(false, "500", "Maximum number of iterations");
    m_config["subiter"] =
        ConfigOption(false, "0", "Number of iterations between updates for r-adaptation");
    m_config["nq"] =
        ConfigOption(false, "-1", "Order of mesh");
    m_config["region"] =
        ConfigOption(false, "0.0", "create regions based on target");
    m_config["resfile"] =
        ConfigOption(false, "", "writes residual values to file");
    m_config["histfile"] =
        ConfigOption(false, "", "histogram of scaled jac");
    m_config["overint"] =
        ConfigOption(false, "6", "over integration order");
    m_config["analytics"] =
        ConfigOption(false, "", "basic analytics module");
    m_config["scalingfile"] =
        ConfigOption(false, "", "Read scaling field from file for r-adaptation");
    // clang-format on
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

    if (m_config["linearelastic"].beenSet)
    {
        m_opti = eLinEl;
    }
    else if (m_config["winslow"].beenSet)
    {
        m_opti = eWins;
    }
    else if (m_config["roca"].beenSet)
    {
        m_opti = eRoca;
    }
    else if (m_config["hyperelastic"].beenSet)
    {
        m_opti = eHypEl;
    }
    else
    {
        ASSERTL0(false, "not opti type set");
    }

    const int maxIter      = m_config["maxiter"].as<int>();
    const NekDouble restol = m_config["restol"].as<NekDouble>();
    int subIter            = m_config["subiter"].as<int>();
    if (subIter == 0)
    {
        subIter = maxIter;
    }

    // m_mesh->m_nummode = m_config["nq"].as<int>();

    bool fd = false;

    if (m_config["nq"].beenSet)
    {
        m_mesh->m_nummode = m_config["nq"].as<int>();
        fd                = true;
    }

    if (!fd)
    {
        for (auto &edge : m_mesh->m_edgeSet)
        {
            if (edge->m_edgeNodes.size() > 0)
            {
                m_mesh->m_nummode = edge->m_edgeNodes.size() + 2;
                fd                = true;
                break;
            }
        }
    }
    ASSERTL0(fd, "failed to find order of mesh");

    // Safety feature: limit over-integration order for high-order triangles
    // over order 5.
    int intOrder = m_config["overint"].as<int>();
    intOrder = m_mesh->m_nummode + intOrder <= 11 ?
        intOrder : 11 - m_mesh->m_nummode;

    if (m_mesh->m_verbose)
    {
        cout << "Identified mesh order as: " << m_mesh->m_nummode - 1 << endl;
    }

    if (m_mesh->m_expDim == 2 && m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false, "cannot deal with manifolds");
    }

    m_res      = std::shared_ptr<Residual>(new Residual);
    m_res->val = 1.0;

    
    m_mesh->MakeOrder(m_mesh->m_nummode - 1,
                      LibUtilities::eGaussLobattoLegendre);

    if (m_config["analytics"].beenSet)
    {
        Analytics();
        return;
    }

    map<LibUtilities::ShapeType, DerivUtilSharedPtr> derivUtils =
        BuildDerivUtil(intOrder);

    GetElementMap(intOrder, derivUtils);

    m_res->startInv = 0;
    m_res->worstJac = numeric_limits<double>::max();
    for (int i = 0; i < m_dataSet.size(); i++)
    {
        m_dataSet[i]->Evaluate();
        m_dataSet[i]->InitialMinJac();
    }

    vector<ElUtilSharedPtr> elLock;

    if (m_config["region"].beenSet)
    {
        elLock = GetLockedElements(m_config["region"].as<NekDouble>());
    }

    
    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes(elLock);
    vector<vector<NodeOptiSharedPtr> > optiNodes;
    
    // turn the free nodes into optimisable objects with all required data
    set<int> check;
    for (int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOptiSharedPtr> ns;
        for (int j = 0; j < freenodes[i].size(); j++)
        {
            auto it = m_nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != m_nodeElMap.end(), "could not find");

            int optiKind = m_mesh->m_spaceDim;

            if (freenodes[i][j]->GetNumCadCurve() == 1)
            {
                optiKind += 10;
            }
            else if (freenodes[i][j]->GetNumCADSurf() == 1)
            {
                optiKind += 20;
            }
            else
            {
                optiKind += 10 * m_mesh->m_expDim;
            }

            auto c = check.find(freenodes[i][j]->m_id);
            ASSERTL0(c == check.end(), "duplicate node");
            check.insert(freenodes[i][j]->m_id);

            ns.push_back(GetNodeOptiFactory().CreateInstance(
                optiKind, freenodes[i][j], it->second, m_res, derivUtils,
                m_opti));
        }
        optiNodes.push_back(ns);
    }

    int nset = optiNodes.size();
    int p    = 0;
    int mn   = numeric_limits<int>::max();
    int mx   = 0;
    for (int i = 0; i < nset; i++)
    {
        p += optiNodes[i].size();
        mn = min(mn, int(optiNodes[i].size()));
        mx = max(mx, int(optiNodes[i].size()));
    }

    if (m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_start.txt";
        histFile.open(name.c_str());

        for (int i = 0; i < m_dataSet.size(); i++)
        {
            histFile << m_dataSet[i]->GetScaledJac() << endl;
        }
        histFile.close();
    }

    if(m_mesh->m_verbose)
    {
        cout << scientific << endl;
        cout << "N elements:\t\t" << m_mesh->m_element[m_mesh->m_expDim].size() - elLock.size() << endl
             << "N elements invalid:\t" << m_res->startInv << endl
             << "Worst jacobian:\t\t" << m_res->worstJac << endl
             << "N free nodes:\t\t" << m_res->n << endl
             << "N Dof:\t\t\t" << m_res->nDoF << endl
             << "N color sets:\t\t" << nset << endl
             << "Avg set colors:\t\t" << p/nset << endl
             << "Min set:\t\t" << mn << endl
             << "Max set:\t\t" << mx << endl
             << "Residual tolerance:\t" << restol << endl;
    }

    int nThreads = m_config["numthreads"].as<int>();

    if (m_mesh->m_cad)
    {
        if (boost::equals(m_mesh->m_cad->GetEngine(), "cfi"))
        {
            WARNINGL0(false,
                      "CFI is not thread-safe; forcing to 'numthreads=1'.");
            nThreads = 1;
        }
    }

    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm =
        tms.CreateInstance(Thread::ThreadMaster::SessionJob, nThreads);

    LibUtilities::Timer t;
    t.Start();

    ofstream resFile;
    if (m_config["resfile"].beenSet)
    {
        resFile.open(m_config["resfile"].as<string>().c_str());
    }

    for (int i = 0; i < optiNodes.size(); i++)
    {
        vector<Thread::ThreadJob *> jobs(optiNodes[i].size());
        for (int j = 0; j < optiNodes[i].size(); j++)
        {
            optiNodes[i][j]->CalcMinJac();
        }
    }

    while (m_res->val > restol && ctr < maxIter)
    {
        ctr++;
        m_res->val       = 0.0;
        m_res->func      = 0.0;
        m_res->nReset[0] = 0;
        m_res->nReset[1] = 0;
        m_res->nReset[2] = 0;
        m_res->alphaI    = 0;
        for (int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob *> jobs(optiNodes[i].size());
            for (int j = 0; j < optiNodes[i].size(); j++)
            {
                jobs[j] = optiNodes[i][j]->GetJob();
            }

            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(nThreads);
            tm->Wait();
        }

        m_res->startInv = 0;
        m_res->worstJac = numeric_limits<double>::max();

        bool update =
            m_config["scalingfile"].beenSet && (ctr % subIter) == 0;

        vector<Thread::ThreadJob *> elJobs(m_dataSet.size());
        for (int i = 0; i < m_dataSet.size(); i++)
        {
            elJobs[i] = m_dataSet[i]->GetJob(update);
        }

        tm->SetNumWorkers(0);
        tm->QueueJobs(elJobs);
        tm->SetNumWorkers(nThreads);
        tm->Wait();

        if (m_config["resfile"].beenSet)
        {
            resFile << m_res->val << " " << m_res->worstJac << " "
                    << m_res->func << endl;
        }

        if(m_mesh->m_verbose)
        {
            cout << ctr
                 << "\tResidual: " << m_res->val
                 << "\tMin Jac: " << m_res->worstJac
                 << "\tInvalid: " << m_res->startInv
                 << "\tReset nodes: " << m_res->nReset[0] << "/" << m_res->nReset[1]
                 << "/" << m_res->nReset[2] << "\tFunctional: " << m_res->func
                 << endl;
        }

        if(ctr >= maxIter)
        {
            break;
        }

        if (m_mesh->m_verbose && update)
        {
            cout << "Mapping updated!" << endl;
        }
    }

    if (m_config["histfile"].beenSet)
    {
        ofstream histFile;
        string name = m_config["histfile"].as<string>() + "_end.txt";
        histFile.open(name.c_str());

        for (int i = 0; i < m_dataSet.size(); i++)
        {
            histFile << m_dataSet[i]->GetScaledJac() << endl;
        }
        histFile.close();
    }
    if (m_config["resfile"].beenSet)
    {
        resFile.close();
    }

    t.Stop();

    //RemoveLinearCurvature();

    if(m_mesh->m_verbose)
    {
        cout << "Time to compute: " << t.TimePerTest(1) << endl;
        cout << "Invalid at end: " << m_res->startInv << endl;
        cout << "Worst at end: " << m_res->worstJac << endl;
    }
}

class NodalUtilTriMonomial : public LibUtilities::NodalUtilTriangle
{
public:
    NodalUtilTriMonomial(int degree, Array<OneD, NekDouble> r,
                         Array<OneD, NekDouble> s)
        : NodalUtilTriangle(degree, r, s)
    {
    }

    virtual ~NodalUtilTriMonomial()
    {
    }

protected:
    virtual NekVector<NekDouble> v_OrthoBasis(const int mode)
    {
        // Monomial basis.
        std::pair<int, int> modes = m_ordering[mode];
        NekVector<NekDouble> ret(m_numPoints);

        for (int i = 0; i < m_numPoints; ++i)
        {
            ret(i) =
                pow(m_xi[0][i], modes.first) * pow(m_xi[1][i], modes.second);
        }

        return ret;
    }

    virtual NekVector<NekDouble> v_OrthoBasisDeriv(const int dir,
                                                   const int mode)
    {
        boost::ignore_unused(dir, mode);
        NEKERROR(ErrorUtil::efatal, "OrthoBasisDeriv: not supported");
        return NekVector<NekDouble>();
    }

    virtual std::shared_ptr<NodalUtil> v_CreateUtil(
        Array<OneD, Array<OneD, NekDouble> > &xi)
    {
        return MemoryManager<NodalUtilTriMonomial>::AllocateSharedPtr(
            m_degree, xi[0], xi[1]);
    }
};

void ProcessVarOpti::Analytics()
{
    // Grab the first element from the list
    ElementSharedPtr elmt = m_mesh->m_element[m_mesh->m_expDim][0];

    // Get curved nodes
    vector<NodeSharedPtr> nodes;
    elmt->GetCurvedNodes(nodes);

    // We're going to investigate only the first node (corner node)
    NodeSharedPtr node = nodes[4];

    // Loop over overintegration orders
    const int nPoints       = 200;
    const int overInt       = 40;
    const NekDouble originX = -1.0;
    const NekDouble originY = -1.0;
    const NekDouble length  = 2.0;
    const NekDouble dx      = length / (nPoints - 1);

    cout << "# overint = " << overInt << endl;
    cout << "# Columns: x, y, over-integration orders (0 -> " << overInt - 1
         << "), "
         << " min(scaledJac)" << endl;

    // Loop over square defined by (originX, originY), length
    for (int k = 0; k < nPoints; ++k)
    {
        node->m_y = originY + k * dx;
        for (int j = 0; j < nPoints; ++j)
        {
            node->m_x = originX + j * dx;
            cout << node->m_x << " " << node->m_y << " ";

            NekDouble minJacNew;

            for (int i = 0; i < overInt; ++i)
            {
                // Clear any existing node to element mapping.
                m_dataSet.clear();
                m_nodeElMap.clear();

                // Build deriv utils and element map.
                map<LibUtilities::ShapeType, DerivUtilSharedPtr> derivUtils =
                    BuildDerivUtil(i);

                // Reconstruct element map
                GetElementMap(i, derivUtils);

                for (int j = 0; j < m_dataSet.size(); j++)
                {
                    m_dataSet[j]->Evaluate();
                    m_dataSet[j]->InitialMinJac();
                }

                // Create NodeOpti object.
                NodeOptiSharedPtr nodeOpti =
                    GetNodeOptiFactory().CreateInstance(
                        m_mesh->m_spaceDim * 11, node,
                        m_nodeElMap.find(node->m_id)->second, m_res, derivUtils,
                        m_opti);

                minJacNew = 0.0;

                // Evaluate functional.
                nodeOpti->CalcMinJac();
                cout << nodeOpti->GetFunctional<2>(minJacNew) << " ";
                // NekDouble eigen;
                // nodeOpti->GetFunctional<2>(minJacNew);
                // nodeOpti->MinEigen<2>(eigen);
                // cout << eigen << " ";
            }

            cout << minJacNew << endl;
        }
    }
}
}
}
