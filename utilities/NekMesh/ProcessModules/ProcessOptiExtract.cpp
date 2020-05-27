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

#include <boost/algorithm/string.hpp>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessOptiExtract.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessOptiExtract::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "opti"),
        ProcessOptiExtract::create,
        "Pulls out blobs for linear elastic solver.");

ProcessOptiExtract::ProcessOptiExtract(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["insert"] =
        ConfigOption(false, "-1", "Name of mesh file to be combined.");
}

ProcessOptiExtract::~ProcessOptiExtract()
{
}

void ProcessOptiExtract::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessOptiExtract: ... " << endl;
    }

    string ins = m_config["insert"].as<string>();

    bool extract = boost::iequals(ins, "-1");

    if (extract)
    {
        vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];

        m_mesh->m_element[m_mesh->m_expDim].clear();
        m_mesh->m_element[m_mesh->m_expDim - 1].clear();

        vector<ElementSharedPtr> invalid;

        // get invalid elements
        for (int i = 0; i < el.size(); ++i)
        {
            // Create elemental geometry.
            SpatialDomains::GeometrySharedPtr geom =
                el[i]->GetGeom(m_mesh->m_spaceDim);

            // Generate geometric factors.
            SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

            // Get the Jacobian and, if it is negative, print a warning
            // message.
            if (!gfac->IsValid())
            {
                invalid.push_back(el[i]);
            }
        }

        std::unordered_set<int> inmesh;
        vector<ElementSharedPtr> totest;

        for (int i = 0; i < invalid.size(); i++)
        {
            auto t = inmesh.insert(invalid[i]->GetId());
            if (t.second)
            {
                m_mesh->m_element[m_mesh->m_expDim].push_back(invalid[i]);
            }

            vector<FaceSharedPtr> f = invalid[i]->GetFaceList();
            for (int j = 0; j < f.size(); j++)
            {
                for (int k = 0; k < f[j]->m_elLink.size(); k++)
                {
                    if (f[j]->m_elLink[k].first.lock()->GetId() ==
                        invalid[i]->GetId())
                    {
                        continue;
                    }

                    t = inmesh.insert(f[j]->m_elLink[k].first.lock()->GetId());
                    if (t.second)
                    {
                        m_mesh->m_element[m_mesh->m_expDim].push_back(
                            f[j]->m_elLink[k].first.lock());
                        totest.push_back(f[j]->m_elLink[k].first.lock());
                    }
                }
            }
        }

        for (int i = 0; i < 12; i++)
        {
            vector<ElementSharedPtr> tmp = totest;
            totest.clear();
            for (int j = 0; j < tmp.size(); j++)
            {
                vector<FaceSharedPtr> f = tmp[j]->GetFaceList();
                for (int k = 0; k < f.size(); k++)
                {
                    for (int l = 0; l < f[k]->m_elLink.size(); l++)
                    {
                        if (f[k]->m_elLink[l].first.lock()->GetId() ==
                            tmp[j]->GetId())
                        {
                            continue;
                        }

                        auto t = inmesh.insert(f[k]->m_elLink[l].first.lock()->GetId());
                        if (t.second)
                        {
                            m_mesh->m_element[m_mesh->m_expDim].push_back(
                                f[k]->m_elLink[l].first.lock());
                            totest.push_back(f[k]->m_elLink[l].first.lock());
                        }
                    }
                }
            }
        }

        ClearElementLinks();
        m_mesh->m_vertexSet.clear();
        m_mesh->m_edgeSet.clear();
        m_mesh->m_faceSet.clear();

        el = m_mesh->m_element[m_mesh->m_expDim];

        if (m_mesh->m_verbose)
        {
            cout << el.size() << " elements in blobs" << endl;
        }

        m_mesh->m_faceSet.clear();

        // re build face links
        for (int i = 0; i < el.size(); ++i)
        {
            for (int j = 0; j < el[i]->GetFaceCount(); ++j)
            {
                auto testIns = m_mesh->m_faceSet.insert(el[i]->GetFace(j));

                if (testIns.second)
                {
                    (*(testIns.first))
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(el[i], j));
                }
                else
                {
                    el[i]->SetFace(j, *testIns.first);
                    // Update face to element map.
                    (*(testIns.first))
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(el[i], j));
                }
            }
        }

        // build surface composite from faces
        for (int i = 0; i < el.size(); i++)
        {
            vector<FaceSharedPtr> f = el[i]->GetFaceList();
            for (int j = 0; j < f.size(); j++)
            {
                if (f[j]->m_elLink.size() == 1)
                {
                    // boundary element make new composite
                    ElmtConfig conf(LibUtilities::eTriangle, 1, false, false);

                    vector<int> tags;
                    tags.push_back(1);
                    ElementSharedPtr E = GetElementFactory().CreateInstance(
                        LibUtilities::eTriangle,
                        conf,
                        f[j]->m_vertexList,
                        tags);
                    m_mesh->m_element[m_mesh->m_expDim - 1].push_back(E);
                }
            }
        }

        ClearElementLinks();
        for (int i = 0; i < el.size(); ++i)
        {
            for (int j = 0; j < el[i]->GetVertexCount(); ++j)
            {
                auto testIns = m_mesh->m_vertexSet.insert(el[i]->GetVertex(j));

                if (!testIns.second)
                {
                    el[i]->SetVertex(j, *testIns.first);
                }
            }
        }

        for (int i = 0; i < el.size(); ++i)
        {
            for (int j = 0; j < el[i]->GetEdgeCount(); ++j)
            {
                EdgeSharedPtr ed = el[i]->GetEdge(j);
                auto testIns     = m_mesh->m_edgeSet.insert(ed);

                if (testIns.second)
                {
                    EdgeSharedPtr ed2 = *testIns.first;
                    ed2->m_elLink.push_back(
                        pair<ElementSharedPtr, int>(el[i], j));
                }
                else
                {
                    EdgeSharedPtr e2 = *(testIns.first);
                    el[i]->SetEdge(j, e2);
                    if (e2->m_edgeNodes.size() == 0 &&
                        ed->m_edgeNodes.size() > 0)
                    {
                        e2->m_curveType = ed->m_curveType;
                        e2->m_edgeNodes = ed->m_edgeNodes;

                        // Reverse nodes if appropriate.
                        if (e2->m_n1->m_id != ed->m_n1->m_id)
                        {
                            reverse(e2->m_edgeNodes.begin(),
                                    e2->m_edgeNodes.end());
                        }
                    }

                    // Update edge to element map.
                    e2->m_elLink.push_back(
                        pair<ElementSharedPtr, int>(el[i], j));
                }
            }
        }
        for (int i = 0; i < el.size(); ++i)
        {
            for (int j = 0; j < el[i]->GetFaceCount(); ++j)
            {
                auto testIns = m_mesh->m_faceSet.insert(el[i]->GetFace(j));

                if (testIns.second)
                {
                    (*(testIns.first))
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(el[i], j));
                }
                else
                {
                    el[i]->SetFace(j, *testIns.first);
                    // Update face to element map.
                    (*(testIns.first))
                        ->m_elLink.push_back(
                            pair<ElementSharedPtr, int>(el[i], j));
                }
            }
        }
        ProcessFaces(false);
        ProcessComposites();
    }
    else
    {
        // insert other mesh
        cout << ins << endl;
        MeshSharedPtr inp_mesh = std::shared_ptr<Mesh>(new Mesh());
        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
            ModuleKey(eInputModule, "xml"), inp_mesh);
        mod->RegisterConfig("infile", ins);
        mod->Process();

        // need to update the vertices manually then the edges and faces can
        // be updated using more simple means.
        map<int, NodeSharedPtr> nmap;

        for (auto &node : inp_mesh->m_vertexSet)
        {
            nmap[node->m_id] = node;
        }
        // for all the nodes in the main mesh see if they are in nmap, if so
        // update the node
        for (auto &node : m_mesh->m_vertexSet)
        {
            if (nmap.count(node->m_id) == 1)
            {
                auto s          = nmap.find(node->m_id);
                NodeSharedPtr n = s->second;
                node->m_x = n->m_x;
                node->m_y = n->m_y;
                node->m_z = n->m_z;
            }
        }

        for (auto &eit : inp_mesh->m_edgeSet)
        {
            auto et = m_mesh->m_edgeSet.find(eit);
            if (et != m_mesh->m_edgeSet.end())
            {
                (*et)->m_edgeNodes = eit->m_edgeNodes;
                (*et)->m_curveType = eit->m_curveType;
            }
        }

        for (auto &fit : inp_mesh->m_faceSet)
        {
            auto ft = m_mesh->m_faceSet.find(fit);
            if (ft != m_mesh->m_faceSet.end())
            {
                (*ft)->m_faceNodes = fit->m_faceNodes;
                (*ft)->m_curveType = fit->m_curveType;
            }
        }
    }
}
}
}
