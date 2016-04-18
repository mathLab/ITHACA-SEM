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

    map<NodeSharedPtr, vector<ElementSharedPtr> > nodeElMap = GetElementMap();

    for(int i = 0; i < optiNodes.size(); i++)
    {
        map<NodeSharedPtr, vector<ElementSharedPtr> >::iterator it;
        it = nodeElMap.find(optiNodes[i]);
        ASSERTL0(it != nodeElMap.end(), "not found");
        cout << i << " " << it->second.size() << endl;
    }

}

map<NodeSharedPtr, vector<ElementSharedPtr> > ProcessVarOpti::GetElementMap()
{
    map<NodeSharedPtr, vector<ElementSharedPtr> > ret;
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

        vector<NodeSharedPtr> n;
        el->GetCurvedNodes(n);
        for(int j = 0; j < 3 * (5 - 1); j++)
        {
            ret[n[j]].push_back(el);
        }
    }
    return ret;
}

vector<NodeSharedPtr> ProcessVarOpti::GetFreeNodes()
{
    //loop over the composites to build a set of nodes which lie in
    //boundary composites
    //then iterate over all nodes, if the node is not in the set its free
    //add it to the vector
    NodeSet boundaryNodes;

    CompositeMap cm = m_mesh->m_composite;
    CompositeMap::iterator it;
    for(it = cm.begin(); it != cm.end(); it++)
    {
        if(it->second->m_tag != "E")
        {
            continue;
        }

        for(int i = 0; i < it->second->m_items.size(); i++)
        {
            EdgeSharedPtr e = it->second->m_items[i]->GetEdgeLink();
            vector<NodeSharedPtr> n;
            e->GetCurvedNodes(n);
            for(int j = 0; j < n.size(); j++)
            {
                boundaryNodes.insert(n[j]);
            }
        }
    }

    vector<NodeSharedPtr> ret;
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); i++)
    {
        vector<NodeSharedPtr> n;
        m_mesh->m_element[m_mesh->m_expDim][i]->GetCurvedNodes(n);
        for(int j = 0; j < 3 * (5 - 1); j++)
        {
            if(!n[j])
            {
                cout << "error in node " << j << endl;
                exit(-1);
            }
            NodeSet::iterator it = boundaryNodes.find(n[j]);
            if(it == boundaryNodes.end())
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
