////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessDetectSurf.cpp
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
#include "ProcessDetectSurf.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessDetectSurf::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "detect"),
        ProcessDetectSurf::create,
        "Process elements to detect a surface.");

ProcessDetectSurf::ProcessDetectSurf(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["vol"] =
        ConfigOption(false, "-1", "Tag identifying surface to process.");
}

ProcessDetectSurf::~ProcessDetectSurf()
{
}

struct EdgeInfo
{
    EdgeInfo() : count(0)
    {
    }
    int count;
    EdgeSharedPtr edge;
    unsigned int group;
};

void ProcessDetectSurf::Process()
{
    if (m_mesh->m_expDim > 2)
    {
        cerr << "Surface detection only implemented for 2D meshes" << endl;
        return;
    }

    int i, j;
    string surf = m_config["vol"].as<string>();

    // Obtain vector of surface IDs from string.
    vector<unsigned int> surfs;
    if (surf != "-1")
    {
        ParseUtils::GenerateSeqVector(surf, surfs);
        sort(surfs.begin(), surfs.end());
    }

    // If we're running in verbose mode print out a list of surfaces.
    if (m_mesh->m_verbose)
    {
        cout << "ProcessDetectSurf: detecting surfaces";
        if (surfs.size() > 0)
        {
            cout << " for surface" << (surfs.size() == 1 ? "" : "s") << " "
                 << surf << endl;
        }
    }

    vector<ElementSharedPtr> &el = m_mesh->m_element[m_mesh->m_expDim];
    map<int, EdgeInfo> edgeCount;
    set<int> doneIds;
    map<int, int> idMap;

    // Iterate over list of surface elements.
    for (i = 0; i < el.size(); ++i)
    {
        // Work out whether this lies on our surface of interest.
        if (surfs.size() > 0)
        {
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
        }

        // List all edges.
        ElementSharedPtr elmt = el[i];
        for (j = 0; j < elmt->GetEdgeCount(); ++j)
        {
            EdgeSharedPtr e = elmt->GetEdge(j);
            int eId         = e->m_id;
            edgeCount[eId].count++;
            edgeCount[eId].edge = e;
        }

        doneIds.insert(elmt->GetId());
        ASSERTL0(idMap.count(elmt->GetId()) == 0, "Shouldn't happen");
        idMap[elmt->GetId()] = i;
    }

    unsigned int maxId = 0;

    for (auto &cIt : m_mesh->m_composite)
    {
        maxId = (std::max)(cIt.first, maxId);
    }

    ++maxId;

    while (doneIds.size() > 0)
    {
        ElementSharedPtr start =
            m_mesh->m_element[m_mesh->m_expDim][idMap[*(doneIds.begin())]];

        vector<ElementSharedPtr> block;
        FindContiguousSurface(start, doneIds, block);
        ASSERTL0(block.size() > 0, "Contiguous block not found");

        // Loop over all edges in block.
        for (i = 0; i < block.size(); ++i)
        {
            // Find edge info.
            ElementSharedPtr elmt = block[i];

            for (j = 0; j < elmt->GetEdgeCount(); ++j)
            {
                auto eIt = edgeCount.find(elmt->GetEdge(j)->m_id);
                ASSERTL0(eIt != edgeCount.end(), "Couldn't find edge");
                eIt->second.group = maxId;
            }
        }

        ++maxId;
    }

    for (auto &eIt : edgeCount)
    {
        if (eIt.second.count > 1)
        {
            continue;
        }

        unsigned int compId = eIt.second.group;
        auto cIt            = m_mesh->m_composite.find(compId);

        if (cIt == m_mesh->m_composite.end())
        {
            CompositeSharedPtr comp(new Composite());
            comp->m_id  = compId;
            comp->m_tag = "E";
            cIt =
                m_mesh->m_composite.insert(std::make_pair(compId, comp)).first;
        }

        vector<int> tags(1);
        tags[0] = compId;
        vector<NodeSharedPtr> nodeList(2);
        nodeList[0] = eIt.second.edge->m_n1;
        nodeList[1] = eIt.second.edge->m_n2;

        ElmtConfig conf(LibUtilities::eSegment, 1, false, false);
        ElementSharedPtr elmt = GetElementFactory().CreateInstance(
            LibUtilities::eSegment, conf, nodeList, tags);
        elmt->SetEdgeLink(eIt.second.edge);

        cIt->second->m_items.push_back(elmt);
    }
}

void ProcessDetectSurf::FindContiguousSurface(ElementSharedPtr start,
                                              set<int> &doneIds,
                                              vector<ElementSharedPtr> &block)
{
    block.push_back(start);
    doneIds.erase(start->GetId());

    vector<EdgeSharedPtr> edges = start->GetEdgeList();

    for (int i = 0; i < edges.size(); ++i)
    {
        for (int j = 0; j < edges[i]->m_elLink.size(); ++j)
        {
            ElementSharedPtr elmt = (edges[i]->m_elLink[j].first).lock();
            if (elmt == start)
            {
                continue;
            }

            if (doneIds.count(elmt->GetId()) == 0)
            {
                continue;
            }

            FindContiguousSurface(elmt, doneIds, block);
        }
    }
}
}
}
