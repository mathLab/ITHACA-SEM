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
//  Description: Extract one or more surfaces from mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include "MeshElements.h"
#include "ProcessExtractSurf.h"

#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

#include <vector>
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessExtractSurf::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "extract"), ProcessExtractSurf::create,
                "Process elements based on values of Jacobian.");

        ProcessExtractSurf::ProcessExtractSurf(MeshSharedPtr m) : ProcessModule(m)
        {
            config["surf"] = ConfigOption(false, "-1",
                "Tag identifying surface to process.");
        }

        ProcessExtractSurf::~ProcessExtractSurf()
        {
        }
        
        void ProcessExtractSurf::Process()
        {
            if (m->verbose)
            {
                cout << "ProcessExtractSurf: Calculating Jacobians..." << endl;
            }

            int i, j, k;
            string surf = config["surf"].as<string>();

            vector<unsigned int> surfs;
            ParseUtils::GenerateSeqVector(surf.c_str(), surfs);
            sort(surfs.begin(), surfs.end());

            vector<ElementSharedPtr> el = m->element[m->expDim-1];

            // Clear list of 3D elements
            m->element[m->expDim]  .clear();
            m->element[m->expDim-1].clear();

            // Clear vertices, edges and faces.
            m->vertexSet.clear();
            m->edgeSet.clear();
            m->faceSet.clear();

            boost::unordered_set<int> keptIds;

            // Iterate over list of elements of expansion dimension.
            for (i = 0; i < el.size(); ++i)
            {
                vector<int> inter;
                vector<int> tags = el[i]->GetTagList();

                sort(tags.begin(), tags.end());
                set_intersection(surfs.begin(), surfs.end(),
                                 tags .begin(), tags .end(),
                                 back_inserter(inter));
                
                if (inter.size() != 1)
                {
                    continue;
                }

                ElementSharedPtr      elmt  = el[i];
                vector<NodeSharedPtr> verts = elmt->GetVertexList();
                vector<EdgeSharedPtr> edges = elmt->GetEdgeList();

                for (j = 0; j < verts.size(); ++j)
                {
                    m->vertexSet.insert(verts[j]);
                }

                FaceSharedPtr f = elmt->GetFaceLink();
                for (j = 0; j < edges.size(); ++j)
                {
                    for (k = 0; k < edges.size(); ++k)
                    {
                        if (f->edgeList[k] == edges[j])
                        {
                            break;
                        }
                    }

                    ASSERTL0(k < edges.size(), "Unable to identify face edges");
                    edges[j]->id = f->edgeList[k]->id;
                    elmt->SetId(f->id);
                    m->edgeSet.insert(edges[j]);
                }

                elmt->SetFaceLink(FaceSharedPtr());
                elmt->SetEdgeLink(EdgeSharedPtr());
                keptIds.insert(elmt->GetId());

                m->element[m->expDim-1].push_back(elmt);
            }

            m->expDim--;

            CompositeMap tmp = m->composite;
            CompositeMap::iterator it;
            
            m->composite.clear();
            int maxId = -1;

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                if (it->second->items[0]->GetDim() != m->expDim)
                {
                    continue;
                }

                vector<ElementSharedPtr> el = it->second->items;
                it->second->items.clear();

                for (i = 0; i < el.size(); ++i)
                {
                    if (keptIds.count(el[i]->GetId()) > 0)
                    {
                        it->second->items.push_back(el[i]);
                    }
                }

                if (it->second->items.size() == 0)
                {
                    continue;
                }

                m->composite.insert(*it);
                maxId = std::max(maxId, (int)it->second->id);
            }

            tmp = m->composite;
            m->composite.clear();
            map<string, CompositeSharedPtr>::iterator it2;

            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                vector<ElementSharedPtr> el = it->second->items;

                string initialTag = el[0]->GetTag();
                it->second->items.resize(1);
                it->second->tag = initialTag;

                map<string, CompositeSharedPtr> newComps;
                newComps[initialTag] = it->second;

                for (i = 1; i < el.size(); ++i)
                {
                    string tag = el[i]->GetTag();
                    if (tag == initialTag)
                    {
                        it->second->items.push_back(el[i]);
                        continue;
                    }

                    it2 = newComps.find(tag);
                    if (it2 == newComps.end())
                    {
                        CompositeSharedPtr newComp(new Composite());
                        newComp->id  = maxId++;
                        newComp->tag = tag;
                        newComp->items.push_back(el[i]);
                    }
                    else
                    {
                        it2->second->items.push_back(el[i]);
                    }
                }

                for (it2 = newComps.begin(); it2 != newComps.end(); ++it2)
                {
                    m->composite[it2->second->id] = it2->second;
                }
            }
        }
    }
}
