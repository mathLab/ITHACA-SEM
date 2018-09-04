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

#include "NekMeshUtils/CADSystem/CFI/CADCurveCFI.h"
#include "NekMeshUtils/CADSystem/CFI/CADSurfCFI.h"

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
    m_config["from"]  = ConfigOption(false, "", "FBM file to load");
    m_config["order"] = ConfigOption(false, "1", "Enforce a polynomial order");
}

OutputCADfix::~OutputCADfix()
{
}

void OutputCADfix::Process()
{
    if (!m_mesh->m_cad)
    {
        ModuleSharedPtr module = GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "loadcad"), m_mesh);
        module->RegisterConfig("CFIMesh", "");
        if (m_mesh->m_verbose)
        {
            module->RegisterConfig("verbose", "");
        }

        // If no input file specified, load output file
        module->RegisterConfig("filename",
                               m_config["from"].beenSet
                                   ? m_config["from"].as<string>()
                                   : m_config["outfile"].as<string>());

        module->SetDefaults();
        module->Process();
    }

    if (m_mesh->m_verbose)
    {
        cout << "OutputCADfix: Writing file..." << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    m_cad          = std::dynamic_pointer_cast<CADSystemCFI>(m_mesh->m_cad);
    m_model        = m_cad->GetCFIModel();
    NekDouble scal = m_cad->GetScaling();

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
            ElementSharedPtr e            = m_mesh->m_element[d][i];
            vector<NodeSharedPtr> volList = e->GetVolumeNodes();
            m_mesh->m_vertexSet.insert(volList.begin(), volList.end());
        }
    }

    // Find main body
    cfi::Body *body;
    vector<cfi::Entity *> *bds =
        m_model->getEntityList(cfi::TYPE_BODY, cfi::SUBTYPE_ALL);

    for (auto &i : *bds)
    {
        cfi::Body *b = static_cast<cfi::Body *>(i);
        if (b->getTopoSubtype() != cfi::SUBTYPE_COMBINED)
        {
            body = b;
            break;
        }
    }

    // map of new nodes
    map<NodeSharedPtr, cfi::Node *> newMap;

    // Write out nodes
    for (auto &it : m_mesh->m_vertexSet)
    {
        newMap[it] = m_model->createOrphanFenode(
            0, it->m_x / scal, it->m_y / scal, it->m_z / scal);

        // Point parent
        if (it->GetNumCadCurve() > 1)
        {
            map<cfi::Point *, int> allVerts;

            for (auto &curve : it->GetCADCurves())
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

            newMap[it]->setParent(maxIt->first);
        }
        // Line parent
        else if (it->GetNumCadCurve())
        {
            vector<CADCurveSharedPtr> curves = it->GetCADCurves();
            cfi::Line *c = std::dynamic_pointer_cast<CADCurveCFI>(curves[0])
                               ->GetCfiPointer();
            newMap[it]->setParent(c);
        }
        // Face parent
        else if (it->GetNumCADSurf())
        {
            vector<CADSurfSharedPtr> surfs = it->GetCADSurfs();
            cfi::Face *s = std::dynamic_pointer_cast<CADSurfCFI>(surfs[0])
                               ->GetCfiPointer();
            newMap[it]->setParent(s);
        }
        // Body parent
        else
        {
            newMap[it]->setParent(body);
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
            cfiNodes.push_back(newMap[node]);
        }

        if (!(order % 2))
        {
            for (auto &edge : nekEdges)
            {
                cfiNodes.push_back(newMap[edge->m_edgeNodes[(order - 1) / 2]]);
            }
            for (auto &face : nekFaces)
            {
                // Could be a triangular face without a face node
                if (face->m_edgeList.size() == 4)
                {
                    cfiNodes.push_back(
                        newMap[face->m_faceNodes[pow(order - 1, 2) / 2]]);
                }
            }
            // Could be an element without a volume node
            if (el->GetTag() == "H")
            {
                cfiNodes.push_back(newMap[nekVNodes[pow(order - 1, 3) / 2]]);
            }
        }

        m_model->createOrphanElement(0, cfi::EntitySubtype(type), cfiNodes);
    }

    m_model->saveCopy(m_config["outfile"].as<string>());
}
}
}
