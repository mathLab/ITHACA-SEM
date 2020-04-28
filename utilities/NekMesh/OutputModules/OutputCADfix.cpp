///////////////////////////////////////////////////////////////////////////////
//
//  File: OutputCADfix.cpp
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
//  Description: CADfix file format output.
//
///////////////////////////////////////////////////////////////////////////////

#include "OutputCADfix.h"

#include <NekMeshUtils/CADSystem/CFI/CADCurveCFI.h>
#include <NekMeshUtils/CADSystem/CFI/CADSurfCFI.h>
#include <NekMeshUtils/CADSystem/CFI/CADElementCFI.h>

using namespace std;
using namespace Nektar::NekMeshUtils;
using namespace Nektar::SpatialDomains;

namespace Nektar
{
namespace Utilities
{
ModuleKey OutputCADfix::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "fbm"), OutputCADfix::create,
    "Appends to a CADfix database file.");

OutputCADfix::OutputCADfix(MeshSharedPtr m) : OutputModule(m)
{
    m_config["order"] = ConfigOption(false, "1", "Enforce a polynomial order");
}

OutputCADfix::~OutputCADfix()
{
}

bool compareT(NodeSharedPtr n1, NodeSharedPtr n2)
{
    int id = n1->GetCADCurves()[0]->GetId();
    return n1->GetCADCurveInfo(id) < n2->GetCADCurveInfo(id);
}

void OutputCADfix::Process()
{
    m_cad = std::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    ASSERTL0(m_cad, "CFI system must be kept in memory")

    if (m_mesh->m_verbose)
    {
        cout << "OutputCADfix: Writing file..." << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    m_model        = m_cad->GetCFIModel();
    NekDouble scal = m_cad->GetScaling();

    int order = m_config["order"].as<int>();

    if (m_mesh->m_verbose)
    {
        cout << "Making mesh of order " << order << endl;
    }

    m_mesh->MakeOrder(order, LibUtilities::ePolyEvenlySpaced);

    // Delete old nodes
    vector<cfi::NodeDefinition> *oldNodes = m_model->getFenodes();
    vector<cfi::Node *> nodesToDel;
    for (int i = 0; i < oldNodes->size(); ++i)
    {
        nodesToDel.push_back((*oldNodes)[i].node);
    }
    m_model->deleteNodes(nodesToDel);

    // Delete old elements
    vector<cfi::ElementDefinition> *oldEls =
        m_model->getElements(cfi::SUBTYPE_ALL, 8);
    vector<cfi::Element *> elsToDel;
    for (int i = 0; i < oldEls->size(); ++i)
    {
        elsToDel.push_back((*oldEls)[i].element);
    }
    m_model->deleteElements(elsToDel);

    // map of new nodes
    map<NodeSharedPtr, cfi::Node *> newMap;

    // Make list of nodes to write out
    // Nodes must be created contiguously for each parent entity
    map<cfi::MeshableEntity *, vector<NodeSharedPtr>> mapParentNode;
    for (auto &el : m_mesh->m_element[3])
    {
        vector<NodeSharedPtr> nodes = el->GetVertexList();

        for (auto &edge : el->GetEdgeList())
        {
            edge->GetCurvedNodes(nodes);
        }

        for (auto &face : el->GetFaceList())
        {
            vector<NodeSharedPtr> fnodes;
            face->GetCurvedNodes(fnodes);
            nodes.insert(nodes.end(), fnodes.begin(), fnodes.end());
        }

        vector<NodeSharedPtr> vnodes = el->GetVolumeNodes();
        nodes.insert(nodes.end(), vnodes.begin(), vnodes.end());

        for (auto &node : nodes)
        {
            if (newMap.count(node))
            {
                continue;
            }

            // Default: volume parent
            CADElementCFISharedPtr cadParent = std::dynamic_pointer_cast<
                CADElementCFI>(el->m_parentCAD);
            ASSERTL0(cadParent, "Expected a CFI parent.");
            cfi::MeshableEntity *parent = cadParent->GetCfiPointer();

            // Point parent
            if (node->GetNumCadCurve() > 1)
            {
                map<cfi::Point *, int> allVerts;

                for (auto &curve : node->GetCADCurves())
                {
                    vector<cfi::Oriented<cfi::TopoEntity *>> *vertList =
                        std::dynamic_pointer_cast<CADCurveCFI>(curve)
                            ->GetCfiPointer()
                            ->getChildList();

                    for (auto &vert : *vertList)
                    {
                        cfi::Point *v = static_cast<cfi::Point *>(vert.entity);
                        if (allVerts.count(v))
                        {
                            allVerts[v]++;
                        }
                        else
                        {
                            allVerts[v] = 1;
                        }
                    }
                }

                // Search for most likely parent vertex
                map<cfi::Point *, int>::iterator maxIt = allVerts.begin();
                for (auto it = allVerts.begin(); it != allVerts.end(); ++it)
                {
                    if (it->second > maxIt->second)
                    {
                        maxIt = it;
                    }
                }

                parent = maxIt->first;
            }
            // Line parent
            else if (node->GetNumCadCurve())
            {
                vector<CADCurveSharedPtr> curves = node->GetCADCurves();
                parent = std::dynamic_pointer_cast<CADCurveCFI>(curves[0])
                             ->GetCfiPointer();
            }
            // Face parent
            else if (node->GetNumCADSurf())
            {
                vector<CADSurfSharedPtr> surfs = node->GetCADSurfs();
                parent = std::dynamic_pointer_cast<CADSurfCFI>(surfs[0])
                             ->GetCfiPointer();
            }

            newMap[node] = NULL;
            mapParentNode[parent].push_back(node);
        }
    }

    // Write out nodes
    for (auto &parent : mapParentNode)
    {
        // Order nodes by parametric coordinate on lines
        if (dynamic_cast<cfi::Line *>(parent.first))
        {
            sort(parent.second.begin(), parent.second.end(), compareT);
        }

        for (auto &node : parent.second)
        {
            newMap[node] = parent.first->createFenode(
                0, node->m_x / scal, node->m_y / scal, node->m_z / scal);
        }
    }

    // Write out elements
    for (auto &el : m_mesh->m_element[3])
    {
        vector<cfi::Node *> cfiNodes;
        vector<NodeSharedPtr> nekNodes  = el->GetVertexList();
        vector<EdgeSharedPtr> nekEdges  = el->GetEdgeList();
        vector<FaceSharedPtr> nekFaces  = el->GetFaceList();
        vector<NodeSharedPtr> nekVNodes = el->GetVolumeNodes();
        int type;

        if (el->GetTag() == "H")
        {
            type = order % 2 ? CFI_SUBTYPE_HE8 : CFI_SUBTYPE_HE27;

            if (!(order % 2))
            {
                // Edges need re-ordering
                // Swapping edges 4->7 with 8->11
                swap_ranges(nekEdges.begin() + 4, nekEdges.begin() + 8,
                            nekEdges.begin() + 8);

                // Faces need re-ordering
                // Moving face 0 to second to last
                swap_ranges(nekFaces.begin(), nekFaces.begin() + 4,
                            nekFaces.begin() + 1);
            }
        }
        else if (el->GetTag() == "R")
        {
            type = order % 2 ? CFI_SUBTYPE_PE6 : CFI_SUBTYPE_PE18;

            // Nodes need re-ordering
            vector<NodeSharedPtr> newNekNodes;
            newNekNodes.push_back(nekNodes[0]);
            newNekNodes.push_back(nekNodes[4]);
            newNekNodes.push_back(nekNodes[1]);
            newNekNodes.push_back(nekNodes[3]);
            newNekNodes.push_back(nekNodes[5]);
            newNekNodes.push_back(nekNodes[2]);
            nekNodes.swap(newNekNodes);

            // Edges need re-ordering
            vector<EdgeSharedPtr> newNekEdges;
            newNekEdges.push_back(nekEdges[4]);
            newNekEdges.push_back(nekEdges[5]);
            newNekEdges.push_back(nekEdges[0]);
            newNekEdges.push_back(nekEdges[7]);
            newNekEdges.push_back(nekEdges[6]);
            newNekEdges.push_back(nekEdges[2]);
            newNekEdges.push_back(nekEdges[3]);
            newNekEdges.push_back(nekEdges[8]);
            newNekEdges.push_back(nekEdges[1]);
            nekEdges.swap(newNekEdges);

            // Faces need re-ordering
            vector<FaceSharedPtr> newNekFaces;
            newNekFaces.push_back(nekFaces[4]);
            newNekFaces.push_back(nekFaces[2]);
            newNekFaces.push_back(nekFaces[0]);
            nekFaces.swap(newNekFaces);
        }
        else if (el->GetTag() == "A")
        {
            type = order % 2 ? CFI_SUBTYPE_TE4 : CFI_SUBTYPE_TE10;
        }
        else
        {
            WARNINGL0(false, "Element type not supported");
            continue;
        }

        for (auto &node : nekNodes)
        {
            cfiNodes.push_back(newMap.at(node));
        }

        if (!(order % 2))
        {
            for (auto &edge : nekEdges)
            {
                cfiNodes.push_back(
                    newMap.at(edge->m_edgeNodes[(order - 1) / 2]));
            }
            for (auto &face : nekFaces)
            {
                // Could be a triangular face without a face node
                if (face->m_edgeList.size() == 4)
                {
                    cfiNodes.push_back(
                        newMap.at(face->m_faceNodes[pow(order - 1, 2) / 2]));
                }
            }
            // Could be an element without a volume node
            if (el->GetTag() == "H")
            {
                cfiNodes.push_back(newMap.at(nekVNodes[pow(order - 1, 3) / 2]));
            }
        }

        CADElementCFISharedPtr cadParent = std::dynamic_pointer_cast<
            CADElementCFI>(el->m_parentCAD);
        ASSERTL0(cadParent, "Expected a CFI parent.");
        cadParent->GetCfiPointer()->createElement(
            0, cfi::EntitySubtype(type), cfiNodes);
    }

    m_model->saveCopy(m_config["outfile"].as<string>());
}
}
}
