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
    map<unsigned int, ElmtConfig> igelmap = InputGmsh::GenElmMap();
    map<unsigned int, ElmtConfig>::iterator it;

    // Populate #InputGmsh::elmMap and use this to construct an
    // inverse mapping from %ElmtConfig to Gmsh ID.
    for (it = igelmap.begin(); it != igelmap.end(); ++it)
    {
        elmMap[it->second] = it->first;
    }
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

    // Open the file stream.
    OpenStream();

    // Write MSH header
    m_mshFile << "$MeshFormat" << endl
              << "2.2 0 8" << endl
              << "$EndMeshFormat" << endl;

    int id = m_mesh->m_vertexSet.size();
    vector<ElementSharedPtr> toComplete;

    int maxOrder = -1;

    // Do first pass over elements of expansion dimension to determine
    // which elements need completion.
    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        ElementSharedPtr e = m_mesh->m_element[m_mesh->m_expDim][i];
        if (e->GetMaxOrder() > maxOrder)
        {
            maxOrder = e->GetMaxOrder();
        }
    }

    // Convert this mesh into a high-order mesh of uniform order.
    cout << "Mesh order of " << maxOrder << " detected" << endl;
    maxOrder = 4;
    m_mesh->MakeOrder(maxOrder, LibUtilities::ePolyEvenlySpaced);

    // Add edge- and face-interior nodes to vertex set.
    EdgeSet::iterator eIt;
    FaceSet::iterator fIt;

    for (eIt = m_mesh->m_edgeSet.begin(); eIt != m_mesh->m_edgeSet.end(); ++eIt)
    {
        m_mesh->m_vertexSet.insert((*eIt)->m_edgeNodes.begin(),
                                   (*eIt)->m_edgeNodes.end());
    }

    for (fIt = m_mesh->m_faceSet.begin(); fIt != m_mesh->m_faceSet.end(); ++fIt)
    {
        m_mesh->m_vertexSet.insert((*fIt)->m_faceNodes.begin(),
                                   (*fIt)->m_faceNodes.end());
    }

    // Do second pass over elements for volume nodes.
    for (int d = 1; d <= 3; ++d)
    {
        for (int i = 0; i < m_mesh->m_element[d].size(); ++i)
        {
            ElementSharedPtr e = m_mesh->m_element[d][i];
            vector<NodeSharedPtr> volList = e->GetVolumeNodes();

            for (int j = 0; j < volList.size(); ++j)
            {
                m_mesh->m_vertexSet.insert(volList[j]);
            }
        }
    }

    // Create ordered set of nodes - not required but looks nicer.
    std::set<NodeSharedPtr>::iterator it;
    std::set<NodeSharedPtr> tmp(m_mesh->m_vertexSet.begin(),
                                m_mesh->m_vertexSet.end());

    // Write out nodes section.
    m_mshFile << "$Nodes" << endl << m_mesh->m_vertexSet.size() << endl;

    for (it = tmp.begin(); it != tmp.end(); ++it)
    {
        m_mshFile << (*it)->m_id+1 << " " << scientific << setprecision(10)
                  << (*it)->m_x << " " << (*it)->m_y << " " << (*it)->m_z
                  << endl;
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

            // Construct inverse of input reordering
            vector<int> reordering = InputGmsh::CreateReordering(elmtType);
            vector<int> inv(tags.size());

            ASSERTL1(tags.size() == reordering.size(),
                     "Reordering map size not equal to element tags.");

            for (int j = 0; j < tags.size(); ++j)
            {
                inv[reordering[j]] = j;
            }

            // Finally write element nodes.
            for (int j = 0; j < tags.size(); ++j)
            {
                m_mshFile << tags[inv[j]] + 1 << " ";
            }

            m_mshFile << endl;
        }
    }
    m_mshFile << "$EndElements" << endl;
}
}
}
