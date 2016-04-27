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

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPrismExp.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{
inline NekDouble GetElFunctional(ElDataSharedPtr d, optimiser opti)
{
    SpatialDomains::GeometrySharedPtr    geom = d->el->GetGeom(2);
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();
    LibUtilities::PointsKeyVector        p    = chi->GetPointsKeys();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
    const int expDim = 2;

    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

    const int pts = deriv[0][0].num_elements();

    ASSERTL0(pts == d->maps.size(), "what");

    Array<OneD, NekDouble> dW(pts);

    for (int k = 0; k < pts; ++k)
    {
        Array<TwoD, NekDouble> jacIdeal(2,2);
        Array<TwoD, NekDouble> jac(2,2);

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                jac[j][i] = deriv[i][j][k];
            }
        }

        //jacIdeal = jac * d->maps[k];
        jacIdeal[0][0] = jac[0][0] * d->maps[k](0,0) + jac[0][1] * d->maps[k](1,0);
        jacIdeal[1][1] = jac[1][0] * d->maps[k](0,1) + jac[1][1] * d->maps[k](1,1);
        jacIdeal[0][1] = jac[0][0] * d->maps[k](0,1) + jac[0][1] * d->maps[k](1,1);
        jacIdeal[1][0] = jac[1][0] * d->maps[k](0,0) + jac[1][1] * d->maps[k](1,0);

        switch (opti)
        {
            case eLinEl:
            {
                NekDouble trEtE = 0.25*(
                                  (jacIdeal[0][0] * jacIdeal[0][0] + jacIdeal[1][0] * jacIdeal[1][0] - 1.0)*
                                  (jacIdeal[0][0] * jacIdeal[0][0] + jacIdeal[1][0] * jacIdeal[1][0] - 1.0) +
                                  (jacIdeal[0][0] * jacIdeal[0][1] + jacIdeal[1][0] * jacIdeal[1][1])*
                                  (jacIdeal[0][0] * jacIdeal[0][1] + jacIdeal[1][0] * jacIdeal[1][1]) +
                                  (jacIdeal[0][0] * jacIdeal[0][1] + jacIdeal[1][0] * jacIdeal[1][1])*
                                  (jacIdeal[0][0] * jacIdeal[0][1] + jacIdeal[1][0] * jacIdeal[1][1]) +
                                  (jacIdeal[0][1] * jacIdeal[0][1] + jacIdeal[1][1] * jacIdeal[1][1] - 1.0)*
                                  (jacIdeal[0][1] * jacIdeal[0][1] + jacIdeal[1][1] * jacIdeal[1][1] - 1.0));
                NekDouble ljacDet = log(jacIdeal[0][0] * jacIdeal[1][1] - jacIdeal[0][1]*jacIdeal[1][0]);

                NekDouble nu = 0.45;
                NekDouble mu = 1.0/2.0/(1.0+nu);
                NekDouble K = 1.0 / 3.0 / (1.0 - 2.0 * nu);
                dW[k] = K *0.5 * ljacDet * ljacDet + mu * trEtE;
                break;
            }

            case eWins:
            {
                NekDouble jacDet = jacIdeal[0][0] * jacIdeal[1][1] - jacIdeal[0][1]*jacIdeal[1][0];
                NekDouble frob = 0.0;

                for (int i = 0; i < expDim; ++i)
                {
                    for (int j = 0; j < expDim; ++j)
                    {
                        frob += jacIdeal[i][j] * jacIdeal[i][j];
                    }
                }

                if(jacDet < 1E-6)
                {
                    jacDet = 1E-6;
                }
                dW[k] = frob / jacDet;
                break;
            }

            case eRoca:
            {
                NekDouble jacDet = jacIdeal[0][0] * jacIdeal[1][1] - jacIdeal[0][1]*jacIdeal[1][0];
                NekDouble frob = 0.0;

                for (int i = 0; i < expDim; ++i)
                {
                    for (int j = 0; j < expDim; ++j)
                    {
                        frob += jacIdeal[i][j] * jacIdeal[i][j];
                    }
                }

                NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet));
                dW[k] = frob / expDim * pow(sigma, 2.0/expDim);
                break;
            }
        }
    }

    d->lastEval = chi->Integral(dW);
    return d->lastEval;
}

inline vector<DNekMat> MappingIdealToRef(SpatialDomains::GeometrySharedPtr geom,
                                 StdRegions::StdExpansionSharedPtr chi)
{
    int dim = geom->GetShapeDim();
    vector<DNekMat> ret;

    if(geom->GetShapeType() == LibUtilities::eQuadrilateral)
    {
        vector<Array<OneD, NekDouble> > xy;
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
                ret.push_back(dxdz);
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTriangle)
    {
        vector<Array<OneD, NekDouble> > xy;
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

        for(int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                DNekMat dxdz(2,2,1.0,eFULL);
                dxdz(0,0) = -xy[0][0]/2.0 + xy[1][0]/2.0;

                dxdz(0,1) = -xy[0][0]/2.0 + xy[2][0]/2.0;

                dxdz(1,0) = -xy[0][1]/2.0 + xy[1][1]/2.0;

                dxdz(1,1) = -xy[0][1]/2.0 + xy[2][1]/2.0;

                dxdz.Invert();
                ret.push_back(dxdz);
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::eTetrahedron)
    {
        vector<Array<OneD, NekDouble> > xyz;
        for(int i = 0; i < geom->GetNumVerts(); i++)
        {
            Array<OneD, NekDouble> loc(3);
            SpatialDomains::PointGeomSharedPtr p = geom->GetVertex(i);
            p->GetCoords(loc);
            xyz.push_back(loc);
        }

        Array<OneD, const LibUtilities::BasisSharedPtr> b = chi->GetBase();
        Array<OneD, NekDouble> u = b[0]->GetZ();
        Array<OneD, NekDouble> v = b[1]->GetZ();
        Array<OneD, NekDouble> z = b[2]->GetZ();

        for(int i = 0; i < b[0]->GetNumPoints(); i++)
        {
            for(int j = 0; j < b[1]->GetNumPoints(); j++)
            {
                for(int k = 0; k < b[2]->GetNumPoints(); k++)
                {
                    DNekMat dxdz(3,3,1.0,eFULL);
                    dxdz(0,0) = -xyz[0][0]/2.0 + xyz[1][0]/2.0;

                    dxdz(0,1) = -xyz[0][0]/2.0 + xyz[2][0]/2.0;

                    dxdz(0,2) = -xyz[0][0]/2.0 + xyz[3][0]/2.0;


                    dxdz(1,0) = -xyz[0][1]/2.0 + xyz[1][1]/2.0;

                    dxdz(1,1) = -xyz[0][1]/2.0 + xyz[2][1]/2.0;

                    dxdz(1,2) = -xyz[0][1]/2.0 + xyz[3][1]/2.0;


                    dxdz(2,0) = -xyz[0][2]/2.0 + xyz[1][2]/2.0;

                    dxdz(2,1) = -xyz[0][2]/2.0 + xyz[2][2]/2.0;

                    dxdz(2,2) = -xyz[0][2]/2.0 + xyz[3][2]/2.0;

                    dxdz.Invert();
                    ret.push_back(dxdz);
                }
            }
        }
    }
    else if(geom->GetShapeType() == LibUtilities::ePrism)
    {
        vector<Array<OneD, NekDouble> > xyz;
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
                    ret.push_back(dxdz);
                }
            }
        }
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

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
    else
    {
        ASSERTL0(false,"not opti type set");
    }

    if(m_mesh->m_expDim == 3 || m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"wrong mesh dim");
    }

    FillQuadPoints();

    GetElementMap();

    vector<vector<NodeSharedPtr> > freenodes = GetColouredNodes();

    vector<vector<NodeOpti> > optiNodes;
    for(int i = 0; i < freenodes.size(); i++)
    {
        vector<NodeOpti> ns;
        for(int j = 0; j < freenodes[i].size(); j++)
        {
            NodeElMap::iterator it = nodeElMap.find(freenodes[i][j]->m_id);
            ASSERTL0(it != nodeElMap.end(),"could not find");
            ns.push_back(NodeOpti(freenodes[i][j],it->second,opti));
        }
        optiNodes.push_back(ns);
    }

    NekDouble functionalStart = 0.0;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        functionalStart += GetElFunctional(dataSet[i], opti);
    }

    cout << scientific << endl;

    NekDouble functionalEnd = functionalStart;
    NekDouble functionalLast = 0.0;
    int ctr = 0;
    Thread::ThreadMaster tms;
    tms.SetThreadingType("ThreadManagerBoost");
    Thread::ThreadManagerSharedPtr tm = tms.CreateInstance(Thread::ThreadMaster::SessionJob, 4);
    while (fabs(functionalLast - functionalEnd) > 1e-5)
    {
        ctr++;
        functionalLast = functionalEnd;
        int c = 0;
        for(int i = 0; i < optiNodes.size(); i++)
        {
            vector<Thread::ThreadJob*> jobs;
            for(int j = 0; j < optiNodes[i].size(); j++)
            {
                jobs.push_back(optiNodes[i][j].GetJob());
            }
            cout << jobs.size() << endl;
            tm->SetNumWorkers(0);
            tm->QueueJobs(jobs);
            tm->SetNumWorkers(1);
            tm->Wait();

            /*for(int j = 0; j < optiNodes[i].size(); j++)
            {
                optiNodes[i][j].Run();
            }*/
        }
        functionalEnd = 0.0;
        for(int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            functionalEnd += dataSet[i]->lastEval;
        }
        cout << ctr << "  " << functionalStart << " " <<
                functionalEnd << endl;
        if(ctr > 1000)
            break;
    }
}

NekDouble dir[4][2] = {{1.0,0},{0,1.0},{-1.0,0},{0,-1.0}};

void ProcessVarOpti::NodeOpti::Optimise()
{
    Array<OneD, NekDouble> G = GetGrad();

    if(sqrt(G[0]*G[0] + G[1]*G[1]) > 1e-3)
    {
        //needs to optimise
        NekDouble currentW = GetFunctional();
        NekDouble xc = node->m_x;
        NekDouble yc = node->m_y;
        NekDouble alpha = 1.0;
        bool found = false;
        while(alpha > 1e-6)
        {
            node->m_x = xc - alpha * G[0];
            node->m_y = yc - alpha * G[1];

            if(GetFunctional() < currentW)
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
            cout << "warning: had to reset node" << endl;
        }
    }
    cout << "done" << endl;
}

Array<OneD, NekDouble> ProcessVarOpti::NodeOpti::GetGrad()
{
    NekDouble xc = node->m_x;
    NekDouble yc = node->m_y;
    NekDouble dx = 1e-6;

    vector<NekDouble> w;

    for(int i = 0; i < 4; i++)
    {
        node->m_x += dir[i][0] * dx;
        node->m_y += dir[i][1] * dx;
        w.push_back(GetFunctional());
        node->m_x = xc;
        node->m_y = yc;
    }

    Array<OneD, NekDouble> ret(2);

    ret[0] = (w[0] - w[2]) / 2.0 / dx;
    ret[1] = (w[1] - w[3]) / 2.0 / dx;

    return ret;
}

NekDouble ProcessVarOpti::NodeOpti::GetFunctional()
{
    NekDouble r = 0.0;
    for(int i = 0; i < data.size(); i++)
    {
        r += GetElFunctional(data[i], opti);
    }
    return r;
}

vector<vector<NodeSharedPtr> > ProcessVarOpti::GetColouredNodes()
{
    //loop over the composites to build a set of nodes which lie in
    //boundary composites
    //then iterate over all nodes, if the node is not in the set its free
    //add it to the vector
    NodeSet boundaryNodes;

    EdgeSet::iterator it;
    int ct = 0;
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        if((*it)->m_elLink.size() == 2)
        {
            continue;
        }

        ct++;

        vector<NodeSharedPtr> n;
        (*it)->GetCurvedNodes(n);
        for(int i = 0; i < n.size(); i++)
        {
            boundaryNodes.insert(n[i]);
        }
    }

    vector<NodeSharedPtr> remain;
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        vector<NodeSharedPtr> n;
        (*it)->GetCurvedNodes(n);
        for(int j = 0; j < n.size(); j++)
        {
            NodeSet::iterator nit = boundaryNodes.find(n[j]);
            if(nit == boundaryNodes.end())
            {
                remain.push_back(n[j]);
            }
        }
    }

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
        ElDataSharedPtr d = boost::shared_ptr<ElData>(new ElData);
        d->el = el;

        SpatialDomains::GeometrySharedPtr    geom = el->GetGeom(m_mesh->m_spaceDim);
        StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();

        vector<LibUtilities::BasisKey> basisKeys;

        for (int i = 0; i < m_mesh->m_expDim; ++i)
        {
            basisKeys.push_back(chi->GetBasis(i)->GetBasisKey());
        }

        StdRegions::StdExpansionSharedPtr chiMod;
        switch(chi->DetShapeType())
        {
            case LibUtilities::eTriangle:
                chiMod = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1]);
                break;
            case LibUtilities::eQuadrilateral:
                chiMod = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1]);
                break;
            case LibUtilities::eTetrahedron:
                chiMod = MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1], basisKeys[2]);
                break;
            case LibUtilities::ePrism:
                chiMod = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(
                    basisKeys[0], basisKeys[1], basisKeys[2]);
                break;
            default:
                ASSERTL0(false, "nope");
        }

        d->maps = MappingIdealToRef(geom, chiMod);

        dataSet.push_back(d);
    }

    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        vector<NodeSharedPtr> n;
        el->GetCurvedNodes(n);
        for(int j = 0; j < 3 * (5 - 1); j++)
        {
            //data set and elements have same index in vector
            nodeElMap[n[j]->m_id].push_back(dataSet[i]);
        }
    }
}

void ProcessVarOpti::FillQuadPoints()
{
    //not all quadrature points are there
    //this function adds GLL points to linear edges
    EdgeSet::iterator it;
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        if((*it)->m_edgeNodes.size() > 0)
        {
            //already high-order ignore
            continue;
        }

        //need to fix nummode at some point, set to 5 to match sqcr.xml
        LibUtilities::PointsKey ekey(5,LibUtilities::eGaussLobattoLegendre);
        Array<OneD, NekDouble> gll;
        LibUtilities::PointsManager()[ekey]->GetPoints(gll);

        Array<OneD, Array<OneD, NekDouble> > xyi(5 - 2);
        for (int k = 1; k < 5 - 1; k++)
        {
            Array<OneD, NekDouble> xy(2);
            xy[0] = (*it)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                    (*it)->m_n2->m_x * (1.0 + gll[k]) / 2.0;
            xy[1] = (*it)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                    (*it)->m_n2->m_y * (1.0 + gll[k]) / 2.0;
            xyi[k-1] = xy;
        }

        vector<NodeSharedPtr> ns;
        for(int i = 0; i < xyi.num_elements(); i++)
        {
            ns.push_back(boost::shared_ptr<Node>(new Node(0,xyi[i][0],xyi[i][1],0.0)));
        }

        (*it)->m_edgeNodes = ns;
        (*it)->m_curveType = LibUtilities::eGaussLobattoLegendre;
    }

    int id = m_mesh->m_vertexSet.size();
    //enumerate the curved nodes for mapping purposes
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        vector<NodeSharedPtr> n;
        (*it)->GetCurvedNodes(n);
        for(int j = 1; j < n.size()-1; j++)
        {
            n[j]->m_id = id++;
        }
    }

}

}
}
