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
ModuleKey ProcessExtrude::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extrude"), ProcessExtrude::create);

ProcessExtrude::ProcessExtrude(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["layers"] =
        ConfigOption(false, "5", "Number of layers to extrude");
    m_config["length"] = ConfigOption(false, "1.0", "Length of extrusion");
}

ProcessExtrude::~ProcessExtrude()
{
}

void ProcessExtrude::Process()
{
    int nLayers      = m_config["layers"].as<int>();
    NekDouble length = m_config["length"].as<NekDouble>();

    NekDouble dz = length / nLayers;

    ASSERTL0(m_mesh->m_spaceDim == 2,
             "Extrude should only be called for a two dimensional mesh");

    // Increment space and expansion dimensions.
    m_mesh->m_spaceDim++;
    m_mesh->m_expDim++;

    // Grab a copy of the existing two-dimensional elements.
    vector<ElementSharedPtr> el = m_mesh->m_element[2];

    // Grab a copy of existing composites.
    CompositeMap oldComp = m_mesh->m_composite;
    if (m_mesh->m_verbose)
    {
        cout << "Boundary composites" << endl;
        for (auto &it : oldComp)
        {
            if (it.second->m_tag != "E")
            {
                continue;
            }
            cout << it.first << "\t" << it.second->m_tag;
            for (int i = 0; i < it.second->m_items.size(); ++i)
            {
                cout << "\t" << it.second->m_items[i]->GetId()
                     << " (" << it.second->m_items[i]->GetVertex(0)
                     << ", " << it.second->m_items[i]->GetVertex(1) << ")";
                vector<NodeSharedPtr> vv = it.second->m_items[i]->GetVertexList();
                cout << "\t(" << vv[0]->GetID()<< ", " << vv[1]->GetID() <<")";
            }
            cout << endl;
        }
    }

    // Reset mesh.
    for (int d = 0; d <= 3; ++d)
    {
        m_mesh->m_element[d].clear();
    }

    NodeSet nodes = m_mesh->m_vertexSet;

    map<int, NodeSharedPtr> id2node;

    for (auto &n : nodes)
    {
        id2node[n->m_id] = n;
    }
    // Save z plane coordinate
    NekDouble z0 = 0;

    // Create vertices for subsequent layers.
    for (int i = 1; i < nLayers + 1; ++i)
    {
        for (auto &n : nodes)
        {
            z0 = n->m_z;
            NodeSharedPtr newNode(
                new Node(i * nodes.size() + n->m_id, n->m_x, n->m_y,
                    n->m_z + i * dz));
            m_mesh->m_vertexSet.insert(newNode);
            id2node[i * nodes.size() + n->m_id] = newNode;
        }
    }

    EdgeSet es = m_mesh->m_edgeSet; // copy edges for curvature

    for (int j = 0; j < nLayers; ++j)
    {
        for (int i = 0; i < el.size(); ++i)
        {
            vector<NodeSharedPtr> verts = el[i]->GetVertexList();
            if (verts.size() == 4)
            {
                vector<NodeSharedPtr> nodeList(8);
                nodeList[0] = id2node[verts[0]->m_id + j * nodes.size()];
                nodeList[1] = id2node[verts[1]->m_id + j * nodes.size()];
                nodeList[2] = id2node[verts[2]->m_id + j * nodes.size()];
                nodeList[3] = id2node[verts[3]->m_id + j * nodes.size()];
                nodeList[4] = id2node[verts[0]->m_id + (j + 1) * nodes.size()];
                nodeList[5] = id2node[verts[1]->m_id + (j + 1) * nodes.size()];
                nodeList[6] = id2node[verts[2]->m_id + (j + 1) * nodes.size()];
                nodeList[7] = id2node[verts[3]->m_id + (j + 1) * nodes.size()];


                vector<int> tags(1);
                tags[0] = 0;

                ElmtConfig conf(LibUtilities::eHexahedron, 1, false, false,
                    false);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eHexahedron, conf, nodeList, tags);

                m_mesh->m_element[3].push_back(E);
            }
            else
            {
                vector<NodeSharedPtr> nodeList(6);
                nodeList[0] = id2node[verts[0]->m_id + (j + 1) * nodes.size()];
                nodeList[1] = id2node[verts[1]->m_id + (j + 1) * nodes.size()];
                nodeList[2] = id2node[verts[1]->m_id + j * nodes.size()];
                nodeList[3] = id2node[verts[0]->m_id + j * nodes.size()];
                nodeList[4] = id2node[verts[2]->m_id + (j + 1) * nodes.size()];
                nodeList[5] = id2node[verts[2]->m_id + j * nodes.size()];

                vector<int> tags(1);
                tags[0] = 1;

                ElmtConfig conf(LibUtilities::ePrism, 1, false, false, false);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::ePrism, conf, nodeList, tags);

                m_mesh->m_element[3].push_back(E);
            }
        }
    }

    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    // Copy edge information
    for (auto &edge : es)
    {
        if (edge->m_edgeNodes.size() > 0)
        {
            for (int j = 0; j < nLayers + 1; ++j)
            {
                vector<NodeSharedPtr> ns(edge->m_edgeNodes.size());
                for (int i = 0; i < ns.size(); i++)
                {
                    NodeSharedPtr n = edge->m_edgeNodes[i];
                    ns[i]           = std::shared_ptr<Node>(
                        new Node(0, n->m_x, n->m_y, n->m_z + j * dz));
                }

                EdgeSharedPtr e = std::shared_ptr<Edge>(
                    new Edge(id2node[edge->m_n1->m_id + j * nodes.size()],
                             id2node[edge->m_n2->m_id + j * nodes.size()]));

                auto f = m_mesh->m_edgeSet.find(e);
                ASSERTL0(f != m_mesh->m_edgeSet.end(), "could not find edge");

                // Copy edge type
                (*f)->m_curveType = edge->m_curveType;
                // Copy points
                if ((*f)->m_n1 == e->m_n1)
                {
                    (*f)->m_edgeNodes = ns;
                }
                else
                {
                    reverse(ns.begin(), ns.end());
                    (*f)->m_edgeNodes = ns;
                }
            }
        }
    }

    // Get composites max id
    unsigned int maxCompId = 0;
    for (auto &it : oldComp)
    {
        if(it.second->m_id >= maxCompId)
        {
            maxCompId = it.second->m_id;
        }
    }

    // First rename surface to volume composites to out of range
    int outCompId = maxCompId + 1;
    std::vector<int> toErase;
    for (auto &it2 : m_mesh->m_composite)
    {
        if (it2.second->m_id > maxCompId)
        {
            // done!
            break;
        }
        if (it2.second->m_tag == "H" || it2.second->m_tag == "R")
        {
            it2.second->m_id = outCompId;
            m_mesh->m_composite.insert(std::make_pair(outCompId,
                it2.second));
            toErase.push_back(it2.first);
            outCompId += 1;
        }
    }

    for (auto &e : toErase)
    {
        m_mesh->m_composite.erase(e);
    }

    toErase.clear();

    // Then copy surface to volume composites names
    for (auto &it2 : m_mesh->m_composite)
    {
        if (it2.second->m_tag == "H" || it2.second->m_tag == "R")
        {
            for (auto &it1 : oldComp)
            {
                if (it2.second->m_tag == "H" && it1.second->m_tag == "Q")
                {
                    it2.second->m_id = it1.second->m_id;
                    m_mesh->m_composite.insert(std::make_pair(it1.second->m_id,
                        it2.second));
                    toErase.push_back(it2.first);
                    oldComp.erase(it1.first);
                    break;
                }
                else if (it2.second->m_tag == "R" && it1.second->m_tag == "T")
                {
                    it2.second->m_id = it1.second->m_id;
                    m_mesh->m_composite.insert(std::make_pair(it1.second->m_id,
                        it2.second));
                    toErase.push_back(it2.first);
                    oldComp.erase(it1.first);
                    break;
                }
            }
        }
    }

    for (auto &e : toErase)
    {
        m_mesh->m_composite.erase(e);
    }

    // Add new composite to be filled with all boundary faces
    CompositeSharedPtr comp(new Composite());
    comp->m_id  = ++maxCompId;
    unsigned int compAllFaceId = maxCompId; // save it so we can remove it later on
    comp->m_tag = "F";
    m_mesh->m_composite.insert(std::make_pair(maxCompId, comp));

    // Add all boundary faces to the composite
    auto allFaceC =  m_mesh->m_composite.find(maxCompId);
    if (m_mesh->m_verbose)
    {
        cout << "Faces boundary list" << endl;
    }
    for (auto &it : m_mesh->m_faceSet)
    {
        // Add to composite if boundary face
        if ( it->m_elLink.size() < 2 )
        {
            if(it->m_vertexList.size() == 3)
            {
                // Triangle
                ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
                vector<int> tags(1);
                tags[0] = 1;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                LibUtilities::eTriangle, conf, it->m_vertexList, tags);
                E->SetId(it->m_id);
                allFaceC->second->m_items.push_back(E);
            }
            else if(it->m_vertexList.size() == 4)
            {
                // Quad
                ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
                vector<int> tags(1);
                tags[0] = 0;
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                LibUtilities::eQuadrilateral, conf, it->m_vertexList, tags);
                E->SetId(it->m_id);
                allFaceC->second->m_items.push_back(E);
            }
        }
    }

    // Create boundary composites
    for (auto &itOc : oldComp)
    {
        CompositeSharedPtr comp(new Composite());
        comp->m_id  = itOc.second->m_id;
        comp->m_tag = "F";
        m_mesh->m_composite.insert(std::make_pair(itOc.second->m_id, comp));
    }
    // Create periodic composites
    for (int i = 0; i < 2; i++)
    {
        CompositeSharedPtr comp(new Composite());
        comp->m_id  = ++maxCompId;
        comp->m_tag = "F";
        m_mesh->m_composite.insert(std::make_pair(maxCompId, comp));
    }

    // Populates boundary composites
    for (auto &itQ : allFaceC->second->m_items)
    {
        // Check if this quad belongs to previous boundary
        for (auto &itOc : oldComp)
        {
            for (int iEd = 0; iEd < itOc.second->m_items.size(); ++iEd)
            {
                int inCommon = 0 ;
                for (int iV = 0; iV < itQ->GetVertexList().size(); iV++)
                {
                    for( int j = 0; j < 2; j++)
                    {
                        if (LibUtilities::IsRealEqual(itQ->GetVertex(iV)->m_x,
                                itOc.second->m_items[iEd]->GetVertex(j)->m_x) &&
                            LibUtilities::IsRealEqual(itQ->GetVertex(iV)->m_y,
                                itOc.second->m_items[iEd]->GetVertex(j)->m_y))
                        {
                            ++inCommon;
                        }
                    }
                }
                // If the face contains 4 xy pairs in common with 1 edge it
                // must be an extruded edge and it should be added to the
                // corresponding composite
                if(inCommon == 4)
                {
                    if (m_mesh->m_verbose)
                    {
                        cout << "Face " << itQ->GetId() << "\t"
                            << "Composite " << itOc.second->m_id << "\t"
                            << endl;
                    }
                    auto newC =  m_mesh->m_composite.find(itOc.second->m_id);
                    // Quad
                    ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false,
                        false);
                    vector<int> tags(1);
                    tags[0] = 0;
                    ElementSharedPtr E = GetElementFactory().CreateInstance(
                    LibUtilities::eQuadrilateral, conf, itQ->GetVertexList(),
                        tags);
                    E->SetId(itQ->GetId());
                    newC->second->m_items.push_back(E);
                }
            }


        }

        // Populates periodic composites
        NekDouble zdist = 0.0;
        for (int iV = 0; iV < itQ->GetVertexList().size(); iV++)
        {
            zdist += itQ->GetVertex(iV)->m_z;
        }
        zdist = zdist / itQ->GetVertexList().size();
        unsigned int compPerId = 0;
        if(LibUtilities::IsRealEqual(zdist, z0))
        {
            compPerId = maxCompId-1;
        }
        else if(LibUtilities::IsRealEqual(zdist - z0, length))
        {
           compPerId = maxCompId;
        }
        if(compPerId > 0 && itQ->GetVertexList().size() == 3)
        {
            // Triangle
            auto perC =  m_mesh->m_composite.find(compPerId);
            ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);
            vector<int> tags(1);
            tags[0] = 1;
            ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, conf, itQ->GetVertexList(), tags);
            E->SetId(itQ->GetId());
            perC->second->m_items.push_back(E);
        }
        else if(compPerId > 0 && itQ->GetVertexList().size() == 4)
        {
            // Quad
            auto perC =  m_mesh->m_composite.find(compPerId);
            ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false, false);
            vector<int> tags(1);
            tags[0] = 0;
            ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::eQuadrilateral, conf, itQ->GetVertexList(), tags);
            E->SetId(itQ->GetId());
            perC->second->m_items.push_back(E);
        }
    }
    // Remove all faces composite
    m_mesh->m_composite.erase(compAllFaceId);



}
}
}
