////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputGmsh.cpp
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
//  Description: Gmsh file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include <NekMeshUtils/MeshElements/Triangle.h>
#include <NekMeshUtils/MeshElements/Quadrilateral.h>

#include "OutputGmsh.h"
#include "../InputModules/InputGmsh.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputGmsh::className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "msh"),
                                               OutputGmsh::create,
                                               "Writes Gmsh msh file.");

OutputGmsh::OutputGmsh(MeshSharedPtr m) : OutputModule(m)
{
    // Populate #InputGmsh::elmMap and use this to construct an
    // inverse mapping from %ElmtConfig to Gmsh ID.
    for (auto &it : InputGmsh::GenElmMap())
    {
        elmMap[it.second] = it.first;
    }

    m_config["order"] = ConfigOption(false, "-1", "Enforce a polynomial order");
}

OutputGmsh::~OutputGmsh()
{
}

/**
 * @brief Process a mesh to output to Gmsh MSH format.
 *
 * Gmsh output is fairly straightforward. The file first contains a
 * list of nodes, followed by a list of elements. Since
 * Mesh::vertexSet only contains vertices of the linear elements, we
 * first loop over the elements so that any high-order vertices can be
 * enumerated and then added to the node list. We then print out the
 * list of nodes and finally print the element list.
 */
void OutputGmsh::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "OutputGmsh: Writing file..." << endl;
    }

    std::unordered_map<int, vector<int> > orderingMap;

    // Open the file stream.
    OpenStream();

    // Write MSH header
    m_mshFile << "$MeshFormat" << endl
              << "2.2 0 8" << endl
              << "$EndMeshFormat" << endl;

    int id = m_mesh->m_vertexSet.size();
    vector<ElementSharedPtr> toComplete;

    int order = m_config["order"].as<int>();

    if (order != -1)
    {
        if (m_mesh->m_verbose)
        {
            cout << "Making mesh of order " << order << endl;
        }
    }
    else
    {
        // Do first pass over elements of expansion dimension to determine
        // which elements need completion.
        for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
        {
            ElementSharedPtr e = m_mesh->m_element[m_mesh->m_expDim][i];
            if (e->GetMaxOrder() > order)
            {
                order = e->GetMaxOrder();
            }
        }
    }

    // Convert this mesh into a high-order mesh of uniform order.
    if (m_mesh->m_verbose)
    {
        cout << "Mesh order of " << order << " detected" << endl;
    }

    m_mesh->MakeOrder(order, LibUtilities::ePolyEvenlySpaced);

    // Add edge- and face-interior nodes to vertex set.
    for (auto &eIt : m_mesh->m_edgeSet)
    {
        m_mesh->m_vertexSet.insert(eIt->m_edgeNodes.begin(),
                                   eIt->m_edgeNodes.end());
    }

    for (auto &fIt : m_mesh->m_faceSet)
    {
        m_mesh->m_vertexSet.insert(fIt->m_faceNodes.begin(),
                                   fIt->m_faceNodes.end());
    }

    // Do second pass over elements for volume nodes.
    for (int d = 1; d <= 3; ++d)
    {
        for (int i = 0; i < m_mesh->m_element[d].size(); ++i)
        {
            ElementSharedPtr e = m_mesh->m_element[d][i];
            vector<NodeSharedPtr> volList = e->GetVolumeNodes();
            m_mesh->m_vertexSet.insert(volList.begin(), volList.end());
        }
    }

    // Create ordered set of nodes - not required but looks nicer.
    map<int, NodeSharedPtr> tmp;
    for (const auto &it : m_mesh->m_vertexSet)
    {
        tmp[it->GetID() + 1] = it;
    }

    // Write out nodes section.
    m_mshFile << "$Nodes" << endl << m_mesh->m_vertexSet.size() << endl;

    for (auto &it : tmp)
    {
        m_mshFile << it.first << " " << scientific << setprecision(10)
                  << it.second->m_x << " " << it.second->m_y << " "
                  << it.second->m_z << endl;
    }

    m_mshFile << "$EndNodes" << endl;

    // Write elements section. All other sections are not currently
    // supported (physical names etc).
    m_mshFile << "$Elements" << endl;
    m_mshFile << m_mesh->GetNumEntities() << endl;

    id = 1;

    for (int d = 1; d <= 3; ++d)
    {
        for (int i = 0; i < m_mesh->m_element[d].size(); ++i, ++id)
        {
            ElementSharedPtr e = m_mesh->m_element[d][i];

            // First output element ID and type.
            int elmtType = elmMap[e->GetConf()];
            m_mshFile << id << " " << elmtType << " ";

            // Write out number of element tags and then the tags
            // themselves.
            vector<int> tags = e->GetTagList();

            if (tags.size() == 1)
            {
                tags.push_back(tags[0]);
                tags.push_back(0);
            }

            m_mshFile << tags.size() << " ";

            for (int j = 0; j < tags.size(); ++j)
            {
                m_mshFile << tags[j] << " ";
            }

            // Finally write out node list. First write vertices, then
            // internal edge nodes, then face nodes.
            vector<NodeSharedPtr> nodeList = e->GetVertexList();
            vector<EdgeSharedPtr> edgeList = e->GetEdgeList();
            vector<FaceSharedPtr> faceList = e->GetFaceList();
            vector<NodeSharedPtr> volList  = e->GetVolumeNodes();

            tags.clear();

            for (int j = 0; j < nodeList.size(); ++j)
            {
                tags.push_back(nodeList[j]->m_id);
            }

            // Process edge-interior points
            for (int j = 0; j < edgeList.size(); ++j)
            {
                nodeList = edgeList[j]->m_edgeNodes;

                if (e->GetEdgeOrient(j, edgeList[j]) == StdRegions::eForwards)
                {
                    for (int k = 0; k < nodeList.size(); ++k)
                    {
                        tags.push_back(nodeList[k]->m_id);
                    }
                }
                else
                {
                    for (int k = nodeList.size() - 1; k >= 0; --k)
                    {
                        tags.push_back(nodeList[k]->m_id);
                    }
                }
            }

            // Process face-interior points
            for (int j = 0; j < faceList.size(); ++j)
            {
                nodeList = faceList[j]->m_faceNodes;
                int nFaceVerts = faceList[j]->m_vertexList.size();
                vector<int> faceIds(nFaceVerts), volFaceIds(nFaceVerts);

                for (int k = 0; k < nFaceVerts; ++k)
                {
                    faceIds   [k] = faceList[j]->m_vertexList[k]->m_id;
                    volFaceIds[k] =
                        e->GetVertexList()[e->GetFaceVertex(j, k)]->m_id;
                }

                if (nFaceVerts == 3)
                {
                    HOTriangle<NodeSharedPtr> hoTri(faceIds, nodeList);
                    hoTri.Align(volFaceIds);
                    for (int k = 0; k < hoTri.surfVerts.size(); ++k)
                    {
                        tags.push_back(hoTri.surfVerts[k]->m_id);
                    }
                }
                else
                {
                    HOQuadrilateral<NodeSharedPtr> hoQuad(faceIds, nodeList);
                    hoQuad.Align(volFaceIds);

                    for (int k = 0; k < hoQuad.surfVerts.size(); ++k)
                    {
                        tags.push_back(hoQuad.surfVerts[k]->m_id);
                    }
                }
            }

            // Process volume nodes
            for (int j = 0; j < volList.size(); ++j)
            {
                tags.push_back(volList[j]->m_id);
            }

            // Construct inverse of input reordering. First try to find it in
            // our cache.
            auto oIt = orderingMap.find(elmtType);

            // If it's not created, then create it.
            if (oIt == orderingMap.end())
            {
                vector<int> reordering = InputGmsh::CreateReordering(elmtType);
                vector<int> inv(tags.size());

                ASSERTL1(tags.size() == reordering.size(),
                         "Reordering map size not equal to element tags 1.");

                for (int j = 0; j < tags.size(); ++j)
                {
                    inv[reordering[j]] = j;
                }

                oIt = orderingMap.insert(make_pair(elmtType, inv)).first;
            }

            ASSERTL1(tags.size() == oIt->second.size(),
                     "Reordering map size not equal to element tags 2.");

            // Finally write element nodes.
            for (int j = 0; j < tags.size(); ++j)
            {
                m_mshFile << tags[oIt->second[j]] + 1 << " ";
            }

            m_mshFile << endl;
        }
    }
    m_mshFile << "$EndElements" << endl;
}
}
}
