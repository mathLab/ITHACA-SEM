///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessExtrude.cpp
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
//  Description: Extrude a two-dimensional mesh to a three-dimensional mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessExtrude.h"

using namespace std;
using namespace Nektar::NekMeshUtils;


namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessExtrude::className = GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extrude"), ProcessExtrude::create);

ProcessExtrude::ProcessExtrude(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["layers"] = ConfigOption(
        false, "5", "Number of layers to extrude");
    m_config["length"] = ConfigOption(
        false, "1.0", "Length of extrusion");
}

ProcessExtrude::~ProcessExtrude()
{
}

void ProcessExtrude::Process()
{
    int nLayers = m_config["layers"].as<int>();
    NekDouble length = m_config["length"].as<NekDouble>();

    NekDouble dz = length / (nLayers - 1);

    ASSERTL0(m_mesh->m_spaceDim == 2,
             "Extrude should only be called for a two dimensional mesh");

    // Increment space and expansion dimensions.
    m_mesh->m_spaceDim++;
    m_mesh->m_expDim++;

    // Grab a copy of the existing two-dimensional elements.
    vector<ElementSharedPtr> el = m_mesh->m_element[2];

    // Reset mesh.
    for (int d = 0; d <= 3; ++d)
    {
        m_mesh->m_element[d].clear();
    }

    NodeSet nodes = m_mesh->m_vertexSet;

    map<int, NodeSharedPtr> id2node;

    NodeSet::iterator it;
    for (it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end(); ++it)
    {
        id2node[(*it)->m_id] = *it;
    }

    // Create vertices for subsequent layers.
    for (int i = 1; i < nLayers; ++i)
    {
        for (it = nodes.begin(); it != nodes.end(); ++it)
        {
            NodeSharedPtr n = *it;
            NodeSharedPtr newNode(new Node(i*nodes.size()+n->m_id, n->m_x, n->m_y, i*dz));
            m_mesh->m_vertexSet.insert(newNode);
            id2node[i*nodes.size()+n->m_id] = newNode;
        }
    }

    for (int j = 0; j < nLayers-1; ++j)
    {
        for (int i = 0; i < el.size(); ++i)
        {
            vector<NodeSharedPtr> quadVerts = el[i]->GetVertexList();
            vector<NodeSharedPtr> nodeList(8);
            nodeList[0] = id2node[quadVerts[0]->m_id + j*nodes.size()];
            nodeList[1] = id2node[quadVerts[1]->m_id + j*nodes.size()];
            nodeList[2] = id2node[quadVerts[1]->m_id + (j+1)*nodes.size()];
            nodeList[3] = id2node[quadVerts[0]->m_id + (j+1)*nodes.size()];
            nodeList[4] = id2node[quadVerts[3]->m_id + j*nodes.size()];
            nodeList[5] = id2node[quadVerts[2]->m_id + j*nodes.size()];
            nodeList[6] = id2node[quadVerts[2]->m_id + (j+1)*nodes.size()];
            nodeList[7] = id2node[quadVerts[3]->m_id + (j+1)*nodes.size()];

            vector<int> tags(1);
            tags[0] = 0;

            ElmtConfig conf(LibUtilities::eHexahedron,1,false,false);
            ElementSharedPtr E = GetElementFactory().
                CreateInstance(LibUtilities::eHexahedron,conf,nodeList,tags);

            m_mesh->m_element[3].push_back(E);
        }
    }

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}
}
}
