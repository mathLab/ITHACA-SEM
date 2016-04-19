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

ModuleKey ProcessVarOpti::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "varopti"),
    ProcessVarOpti::create,
    "Optimise mesh locations.");

ProcessVarOpti::ProcessVarOpti(MeshSharedPtr m) : ProcessModule(m)
{

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

    if(m_mesh->m_expDim == 3 || m_mesh->m_spaceDim == 3)
    {
        ASSERTL0(false,"wrong mesh dim");
    }

    FillQuadPoints();

    vector<NodeSharedPtr> optiNodes = GetFreeNodes();

    GetElementMap();

    for(int i = 0; i < optiNodes.size(); i++)
    {
        NekDouble w = GetFunctional(optiNodes[i]);
    }

}

NekDouble ProcessVarOpti::GetFunctional(NodeSharedPtr n)
{
    NodeElMap::iterator it = nodeElMap.find(n);
    ASSERTL0(it != nodeElMap.end(),"could not find");
    vector<ElementSharedPtr> els = it->second;

    NekDouble r = 0.0;
    for(int i = 0; i < els.size(); i++)
    {
        r += GetElFunctional(els[i]);
    }
}

NekDouble ProcessVarOpti::GetElFunctional(ElementSharedPtr el)
{
    SpatialDomains::GeometrySharedPtr    geom = el->GetGeom(m_mesh->m_spaceDim);
    StdRegions::StdExpansionSharedPtr    chi  = geom->GetXmap();
    LibUtilities::PointsKeyVector        p    = chi->GetPointsKeys();
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();
    const int expDim = chi->GetNumBases();
    int nElemPts = 1;

    vector<LibUtilities::BasisKey> basisKeys;

    for (int i = 0; i < expDim; ++i)
    {
        basisKeys.push_back(chi->GetBasis(i)->GetBasisKey());
    }

    for(int i = 0; i < basisKeys.size(); i++)
    {
        cout << basisKeys[i].GetTotNumPoints() << endl;
    }
    exit(-1);

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

    SpatialDomains::DerivStorage deriv = gfac->GetDeriv(p);

    const int pts = deriv[0][0].num_elements();
    cout << pts << endl;
    exit(-1);
    /*
    const int nq  = chiMod->GetTotPoints();

    ASSERTL0(pts == nq, "what");

    vector<DNekMat> i2rm = MappingIdealToRef(geom, chiMod);
    Array<OneD, NekDouble> eta(nq);

    for (int k = 0; k < pts; ++k)
    {
        DNekMat jac     (expDim, expDim, 0.0, eFULL);
        DNekMat jacIdeal(expDim, expDim, 0.0, eFULL);

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                jac(j,i) = deriv[i][j][k];
            }
        }

        jacIdeal = jac * i2rm[k];
        NekDouble jacDet;

        if(expDim == 2)
        {
            jacDet = jacIdeal(0,0) * jacIdeal(1,1) - jacIdeal(0,1)*jacIdeal(1,0);
        }
        else if(expDim == 3)
        {
            jacDet = jacIdeal(0,0) * (jacIdeal(1,1)*jacIdeal(2,2) - jacIdeal(2,1)*jacIdeal(1,2)) -
                     jacIdeal(0,1) * (jacIdeal(1,0)*jacIdeal(2,2) - jacIdeal(2,0)*jacIdeal(1,2)) +
                     jacIdeal(0,2) * (jacIdeal(1,0)*jacIdeal(2,1) - jacIdeal(2,0)*jacIdeal(1,1));
        }
        else
        {
            ASSERTL0(false,"silly exp dim");
        }

        NekDouble frob = 0.0;

        for (int i = 0; i < expDim; ++i)
        {
            for (int j = 0; j < expDim; ++j)
            {
                frob += jacIdeal(i,j) * jacIdeal(i,j);
            }
        }

        NekDouble sigma = 0.5*(jacDet + sqrt(jacDet*jacDet));
        eta[k] = expDim * pow(sigma, 2.0/expDim) / frob;
    }*/
}

void ProcessVarOpti::GetElementMap()
{
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        vector<NodeSharedPtr> n;
        el->GetCurvedNodes(n);
        for(int j = 0; j < 3 * (5 - 1); j++)
        {
            nodeElMap[n[j]].push_back(el);
        }
    }
}

vector<NodeSharedPtr> ProcessVarOpti::GetFreeNodes()
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

    vector<NodeSharedPtr> ret;
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        vector<NodeSharedPtr> n;
        (*it)->GetCurvedNodes(n);
        for(int j = 0; j < n.size(); j++)
        {
            NodeSet::iterator nit = boundaryNodes.find(n[j]);
            if(nit == boundaryNodes.end())
            {
                ret.push_back(n[j]);
            }
        }
    }

    return ret;
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
}

}
}
