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
            int i, j;
            string surf = config["surf"].as<string>();

            // Obtain vector of surface IDs from string.
            vector<unsigned int> surfs;
            ParseUtils::GenerateSeqVector(surf.c_str(), surfs);
            sort(surfs.begin(), surfs.end());

            // If we're running in verbose mode print out a list of surfaces.
            if (m->verbose)
            {
                cout << "ProcessExtractSurf: extracting surface"
                     << (surfs.size() > 1 ? "s" : "") << " " << surf << endl;
            }

            // Make a copy of all existing elements of one dimension lower.
            vector<ElementSharedPtr> el = m->element[m->expDim-1];

            // Clear all elements.
            m->element[m->expDim]  .clear();
            m->element[m->expDim-1].clear();

            // Clear existing vertices, edges and faces.
            m->vertexSet.clear();
            m->edgeSet.clear();
            m->faceSet.clear();

            // Clear all edge -> element links.
            for (i = 0; i < el.size(); ++i)
            {
                vector<EdgeSharedPtr> edges = el[i]->GetEdgeList();
                for (j = 0; j < edges.size(); ++j)
                {
                    edges[j]->elLink.clear();
                }

                FaceSharedPtr f = el[i]->GetFaceLink();
                if (f)
                {
                    for (j = 0; j < f->edgeList.size(); ++j)
                    {
                        f->edgeList[j]->elLink.clear();
                    }
                }
            }

            // keptIds stores IDs of elements we processed earlier.
            boost::unordered_set<int> keptIds;

            // Iterate over list of surface elements.
            for (i = 0; i < el.size(); ++i)
            {
                // Work out whether this lies on our surface of interest.
                vector<int> inter, tags = el[i]->GetTagList();

                sort(tags.begin(), tags.end());
                set_intersection(surfs.begin(), surfs.end(),
                                 tags .begin(), tags .end(),
                                 back_inserter(inter));
                
                // It doesn't continue to next element.
                if (inter.size() != 1)
                {
                    continue;
                }

                // Get list of element vertices and edges.
                ElementSharedPtr      elmt  = el[i];
                vector<NodeSharedPtr> verts = elmt->GetVertexList();
                vector<EdgeSharedPtr> edges = elmt->GetEdgeList();

                // Insert surface vertices.
                for (j = 0; j < verts.size(); ++j)
                {
                    m->vertexSet.insert(verts[j]);
                }

                // Problem: edges and element IDs aren't enumerated with
                // geometry IDs by some input modules/the Module ProcessEdges
                // function. Get around this by replacing everything in the
                // edge/face with information from edge/face link.
                EdgeSharedPtr e = elmt->GetEdgeLink();
                FaceSharedPtr f = elmt->GetFaceLink();
                if (e)
                {
                    elmt->SetId(e->id);
                }
                else if (f)
                {
                    for (j = 0; j < edges.size(); ++j)
                    {
                        m->edgeSet.insert(f->edgeList[j]);
                        elmt->SetEdge(j, f->edgeList[j]);
                        f->edgeList[j]->elLink.push_back(
                            std::make_pair(elmt, j));
                    }
                    elmt->SetId(f->id);
                }
                else
                {
                    for (j = 0; j < edges.size(); ++j)
                    {
                        m->edgeSet.insert(edges[j]);
                    }
                }

                // Nullify edge/face links to get correct tag
                elmt->SetFaceLink(FaceSharedPtr());
                elmt->SetEdgeLink(EdgeSharedPtr());
                keptIds.insert(elmt->GetId());

                // Push element back into the list.
                m->element[m->expDim-1].push_back(elmt);
            }

            // Decrement the expansion dimension to get manifold embedding.
            m->expDim--;

            // Now process composites. This is necessary because 2D surfaces may
            // contain both quadrilaterals and triangles and so need to be split
            // up.
            CompositeMap tmp = m->composite;
            CompositeMap::iterator it;
            
            m->composite.clear();
            int maxId = -1;

            // Loop over composites for first time to determine any composites
            // which don't have elements of the correct dimension.
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

                // Figure out the maximum ID so if we need to create new
                // composites we can give them a unique ID.
                maxId = std::max(maxId, (int)it->second->id) + 1;
            }

            tmp = m->composite;
            m->composite.clear();
            map<string, CompositeSharedPtr>::iterator it2;

            // Now do another loop over the composites to remove composites
            // which don't contain any elements in the new mesh.
            for (it = tmp.begin(); it != tmp.end(); ++it)
            {
                CompositeSharedPtr c = it->second;
                vector<ElementSharedPtr> el = c->items;

                // Remove all but the first element from this composite.
                string initialTag = el[0]->GetTag();
                c->items.resize(1);
                c->tag = initialTag;

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
                    it2 = newComps.find(tag);
                    if (it2 == newComps.end())
                    {
                        CompositeSharedPtr newComp(new Composite());
                        newComp->id  = maxId++;
                        newComp->tag = tag;
                        newComp->items.push_back(el[i]);
                        newComps[tag] = newComp;
                    }
                    else
                    {
                        it2->second->items.push_back(el[i]);
                    }
                }

                // Print out mapping information if we remapped composite IDs.
                if (m->verbose && newComps.size() > 1)
                {
                    cout << "- Mapping composite " << it->first << " ->";
                }

                // Insert new composites.
                for (i = 0, it2 = newComps.begin(); it2 != newComps.end();
                     ++it2, ++i)
                {
                    if (m->verbose && newComps.size() > 1)
                    {
                        cout << (i > 0 ? ", " : " ") << it2->second->id << "("
                             << it2->second->tag << ")";
                    }
                    m->composite[it2->second->id] = it2->second;
                }

                if (m->verbose && newComps.size() > 1)
                {
                    cout << endl;
                }
            }
        }
    }
}
