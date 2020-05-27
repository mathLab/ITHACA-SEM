////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessExtractSurf.cpp
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
//  Description: Extract one or more surfaces from mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessExtractSurf.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessExtractSurf::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extract"),
        ProcessExtractSurf::create,
        "Process elements to extract a specified surface(s) or composites(s).");

ProcessExtractSurf::ProcessExtractSurf(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["surf"] = ConfigOption(
        false, "NotSet", "Tag identifying surface/composite to process.");
    m_config["detectbnd"] =
        ConfigOption(false, "-1", "Tag to detect on boundary composites");
}

ProcessExtractSurf::~ProcessExtractSurf()
{
}

void ProcessExtractSurf::Process()
{
    int i, j;
    string surf    = m_config["surf"].as<string>();
    bool detectbnd = m_config["detectbnd"].beenSet;

    // Obtain vector of surface IDs from string.
    vector<unsigned int> surfs;
    ASSERTL0(ParseUtils::GenerateSeqVector(surf, surfs),
             "Failed to interp surf string. Have you specified this string?");
    sort(surfs.begin(), surfs.end());

    // If we're running in verbose mode print out a list of surfaces.
    if (m_mesh->m_verbose)
    {
        cout << "ProcessExtractSurf: extracting surface"
             << (surfs.size() > 1 ? "s" : "") << " " << surf << endl;
    }

    // Make a copy of all existing elements of one dimension lower.
    vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim - 1];

    // Clear all elements.
    m_mesh->m_element[m_mesh->m_expDim].clear();
    m_mesh->m_element[m_mesh->m_expDim - 1].clear();

    // Clear existing vertices, edges and faces.
    m_mesh->m_vertexSet.clear();

    m_mesh->m_edgeSet.clear();
    m_mesh->m_faceSet.clear();

    // Clear all edge -> element links.
    for (i = 0; i < el.size(); ++i)
    {
        vector<EdgeSharedPtr> edges = el[i]->GetEdgeList();
        for (j = 0; j < edges.size(); ++j)
        {
            edges[j]->m_elLink.clear();
        }

        FaceSharedPtr f = el[i]->GetFaceLink();
        if (f)
        {
            for (j = 0; j < f->m_edgeList.size(); ++j)
            {
                f->m_edgeList[j]->m_elLink.clear();
            }
        }
    }

    // keptIds stores IDs of elements we processed earlier.
    std::unordered_set<int> keptIds;

    EdgeSet bndEdgeSet;

    // Iterate over list of surface elements.
    for (i = 0; i < el.size(); ++i)
    {
        // Work out whether this lies on our surface of interest.
        vector<int> inter, tags = el[i]->GetTagList();

        sort(tags.begin(), tags.end());
        set_intersection(surfs.begin(),
                         surfs.end(),
                         tags.begin(),
                         tags.end(),
                         back_inserter(inter));

        // It doesn't continue to next element.
        if (inter.size() != 1)
        {
            continue;
        }

        // Get list of element vertices and edges.
        ElementSharedPtr elmt       = el[i];
        vector<NodeSharedPtr> verts = elmt->GetVertexList();
        vector<EdgeSharedPtr> edges = elmt->GetEdgeList();

        // Insert surface vertices.
        for (j = 0; j < verts.size(); ++j)
        {
            m_mesh->m_vertexSet.insert(verts[j]);
        }

        // Problem: edges and element IDs aren't enumerated with
        // geometry IDs by some input modules/the Module ProcessEdges
        // function. Get around this by replacing everything in the
        // edge/face with information from edge/face link.
        EdgeSharedPtr e = elmt->GetEdgeLink();
        FaceSharedPtr f = elmt->GetFaceLink();
        if (e)
        {
            elmt->SetId(e->m_id);
        }
        else if (f)
        {
            for (j = 0; j < f->m_vertexList.size(); j++)
            {
                elmt->SetVertex(j, f->m_vertexList[j]);
            }

            for (j = 0; j < edges.size(); ++j)
            {
                m_mesh->m_edgeSet.insert(f->m_edgeList[j]);
                elmt->SetEdge(j, f->m_edgeList[j]);
                f->m_edgeList[j]->m_elLink.push_back(std::make_pair(elmt, j));

                // generate a list of edges on boundary of surfaces being
                // extracted
                auto edit = bndEdgeSet.find(f->m_edgeList[j]);
                if (edit != bndEdgeSet.end())
                {
                    // remove since visited more than once
                    bndEdgeSet.erase(edit);
                }
                else
                {
                    bndEdgeSet.insert(f->m_edgeList[j]);
                }
            }
            elmt->SetVolumeNodes(f->m_faceNodes);
            elmt->SetId(f->m_id);
            elmt->SetCurveType(f->m_curveType);
        }
        else
        {
            for (j = 0; j < edges.size(); ++j)
            {
                m_mesh->m_edgeSet.insert(edges[j]);
            }
        }

        // Nullify edge/face links to get correct tag
        elmt->SetFaceLink(FaceSharedPtr());
        elmt->SetEdgeLink(EdgeSharedPtr());
        keptIds.insert(elmt->GetId());

        // Push element back into the list.
        m_mesh->m_element[m_mesh->m_expDim - 1].push_back(elmt);
    }

    // Decrement the expansion dimension to get manifold embedding.
    m_mesh->m_expDim--;

    // Now process composites. This is necessary because 2D surfaces may
    // contain both quadrilaterals and triangles and so need to be split
    // up.
    CompositeMap tmp = m_mesh->m_composite;

    m_mesh->m_composite.clear();
    int maxId = -1;

    // Loop over composites for first time to determine any composites
    // which don't have elements of the correct dimension.
    for (auto &it : tmp)
    {
        if (it.second->m_items[0]->GetDim() != m_mesh->m_expDim)
        {
            continue;
        }

        vector<ElementSharedPtr> el = it.second->m_items;
        it.second->m_items.clear();

        for (i = 0; i < el.size(); ++i)
        {
            if (keptIds.count(el[i]->GetId()) > 0)
            {
                it.second->m_items.push_back(el[i]);
            }
        }

        if (it.second->m_items.size() == 0)
        {
            continue;
        }

        m_mesh->m_composite.insert(it);

        // Figure out the maximum ID so if we need to create new
        // composites we can give them a unique ID.
        maxId = (std::max)(maxId, (int)it.second->m_id) + 1;
    }

    tmp = m_mesh->m_composite;
    m_mesh->m_composite.clear();

    // Now do another loop over the composites to remove composites
    // which don't contain any elements in the new mesh.
    for (auto &it : tmp)
    {
        CompositeSharedPtr c        = it.second;
        vector<ElementSharedPtr> el = c->m_items;

        // Remove all but the first element from this composite.
        string initialTag = el[0]->GetTag();
        c->m_items.resize(1);
        c->m_tag = initialTag;

        // newComps stores the new composites. The key is the composite
        // type (e.g. Q for quad) and value is the composite.
        map<string, CompositeSharedPtr> newComps;
        newComps[initialTag] = c;

        // Loop over remaining elements in composite and figure out
        // whether it needs to be split up.
        for (i = 1; i < el.size(); ++i)
        {
            // See if tag exists. If it does, we append this to the
            // composite, otherwise we create a new composite and store
            // it in newComps.
            string tag = el[i]->GetTag();
            auto it2 = newComps.find(tag);
            if (it2 == newComps.end())
            {
                CompositeSharedPtr newComp(new Composite());
                newComp->m_id  = maxId++;
                newComp->m_tag = tag;
                newComp->m_items.push_back(el[i]);
                newComps[tag] = newComp;
            }
            else
            {
                it2->second->m_items.push_back(el[i]);
            }
        }

        // Print out mapping information if we remapped composite IDs.
        if (m_mesh->m_verbose && newComps.size() > 1)
        {
            cout << "- Mapping composite " << it.first << " ->";
        }

        // Insert new composites.
        i = 0;
        for (auto &it2 : newComps)
        {
            if (m_mesh->m_verbose && newComps.size() > 1)
            {
                cout << (i > 0 ? ", " : " ") << it2.second->m_id << "("
                     << it2.second->m_tag << ")";
            }
            m_mesh->m_composite[it2.second->m_id] = it2.second;
            ++i;
        }

        if (m_mesh->m_verbose && newComps.size() > 1)
        {
            cout << endl;
        }
    }

    // Detect composites for boundaries. This is done by looping over all
    // elements identifiying if they are not part of required surfaces and if
    // not setting up a list of boundary edges (identified by only being visited
    // once). This list is then compared against an earlier identification of
    // boundary edges on the required surfaces and if the two overlap add a
    // segment element and put segment element in composite as well
    if (detectbnd)
    {
        if (m_mesh->m_expDim != 2)
        {
            cerr << "Surface boundary detection only implemented for 2D meshes"
                 << endl;
            return;
        }

        map<int, EdgeSet> surfBndEdgeSet;
        map<int, string> surfLabels;

        // Iterate over list of surface elements.
        for (i = 0; i < el.size(); ++i)
        {
            // Work out whether this lies on our surface of interest.
            vector<int> inter, tags = el[i]->GetTagList();

            ASSERTL0(tags.size() == 1, "Not sure what mutliple tags implies");

            sort(tags.begin(), tags.end());
            set_intersection(surfs.begin(),
                             surfs.end(),
                             tags.begin(),
                             tags.end(),
                             back_inserter(inter));

            // It does so continue to next element.
            if (inter.size() == 1)
            {
                continue;
            }

            int surf = tags[0];

            // gather surface labels if they exist.
            if (m_mesh->m_faceLabels.count(surf))
            {
                surfLabels[surf] = m_mesh->m_faceLabels[surf];
            }

            // Get list of element vertices and edges.
            ElementSharedPtr elmt       = el[i];
            vector<EdgeSharedPtr> edges = elmt->GetEdgeList();

            FaceSharedPtr f = elmt->GetFaceLink();
            if (f)
            {
                for (j = 0; j < edges.size(); ++j)
                {
                    // generate a list of edges on boundary of surfaces being
                    // extracted
                    if (surfBndEdgeSet.count(surf))
                    {
                        auto edit = surfBndEdgeSet[surf].find(f->m_edgeList[j]);
                        if (edit != surfBndEdgeSet[surf].end())
                        {
                            // remove since visited more than once
                            surfBndEdgeSet[surf].erase(edit);
                        }
                        else
                        {
                            surfBndEdgeSet[surf].insert(f->m_edgeList[j]);
                        }
                    }
                    else
                    {
                        EdgeSet newEdgeSet;
                        surfBndEdgeSet[surf] = newEdgeSet;
                        surfBndEdgeSet[surf].insert(f->m_edgeList[j]);
                    }
                }
            }
        }

        m_mesh->m_faceLabels.clear();

        // iteratve over surfBndEdgeSet and see if they are in BndEdgeSet
        for (auto &esetit : surfBndEdgeSet)
        {
            CompositeSharedPtr newComp(new Composite());
            newComp->m_id  = maxId;
            newComp->m_tag = "E";
            // set up labels if they exist
            if (surfLabels.count(esetit.first))
            {
                newComp->m_label = surfLabels[esetit.first];
            }

            for (auto &edit : esetit.second)
            {
                auto locit = bndEdgeSet.find(edit);
                if (locit != bndEdgeSet.end())
                {
                    // make 1D segment element
                    LibUtilities::ShapeType elType = LibUtilities::eSegment;

                    vector<int> tags;
                    tags.push_back(maxId);

                    // make unique node list
                    vector<NodeSharedPtr> nodeList;
                    nodeList.push_back((*locit)->m_n1);
                    nodeList.push_back((*locit)->m_n2);

                    ElmtConfig conf(elType, 1, true, true);
                    ElementSharedPtr E = GetElementFactory().CreateInstance(
                        elType, conf, nodeList, tags);
                    E->SetId((*locit)->m_id);
                    m_mesh->m_element[E->GetDim()].push_back(E);
                    newComp->m_items.push_back(E);
                }
            }

            if (newComp->m_items.size())
            {
                m_mesh->m_composite[maxId++] = newComp;
            }
        }
    }
}
}
}
