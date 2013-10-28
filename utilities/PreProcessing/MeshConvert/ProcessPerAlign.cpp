////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPerAlign.cpp
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
//  Description: Reorder composites to align periodic boundaries.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessPerAlign.h"

#include <LocalRegions/SegExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/NodalTriExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessPerAlign::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "peralign"),
                ProcessPerAlign::create);

        /**
         * @class ProcessPerAlign
         */

        /**
         * @brief Default constructor.
         */
        ProcessPerAlign::ProcessPerAlign(MeshSharedPtr m) : ProcessModule(m)
        {
            config["surf1"]  = ConfigOption(false, "-1",
                "Tag identifying first surface.");
            config["surf2"]  = ConfigOption(false, "-1",
                "Tag identifying first surface.");
            config["dir"]    = ConfigOption(false, "",
                "Direction in which to align (either x, y, or z)");
            config["orient"] = ConfigOption(true,  "0",
                "Attempt to reorient tets and prisms");
        }

        /**
         * @brief Destructor.
         */
        ProcessPerAlign::~ProcessPerAlign()
        {

        }

        void ProcessPerAlign::Process()
        {
            static int tetFaceNodes[4][3] = {
                {0,1,2},{0,1,3},{1,2,3},{0,2,3}};

            int    surf1  = config["surf1"]. as<int>   ();
            int    surf2  = config["surf2"]. as<int>   ();
            string dir    = config["dir"].   as<string>();
            bool   orient = config["orient"].as<bool>  ();

            if (surf1 == -1)
            {
                cerr << "WARNING: surf1 must be set to a positive integer. "
                     << "Skipping periodic alignment." << endl;
                return;
            }

            if (surf2 == -1)
            {
                cerr << "WARNING: surf2 must be set to a positive integer. "
                     << "Skipping periodic alignment." << endl;
                return;
            }

            if (dir != "x" && dir != "y" && dir != "z")
            {
                cerr << "WARNING: dir must be set to either x, y or z. "
                     << "Skipping periodic alignment." << endl;
                return;
            }

            NekDouble vec[3];
            vec[0] = dir == "x" ? 1.0 : 0.0;
            vec[1] = dir == "y" ? 1.0 : 0.0;
            vec[2] = dir == "z" ? 1.0 : 0.0;

            CompositeMap::iterator it1 = m->composite.find(surf1);
            CompositeMap::iterator it2 = m->composite.find(surf2);

            if (it1 == m->composite.end())
            {
                cerr << "WARNING: Couldn't find surface " << surf1
                     << ". Skipping periodic alignment." << endl;
                return;
            }

            if (it2 == m->composite.end())
            {
                cerr << "WARNING: Couldn't find surface " << surf2 << ", "
                     << "skipping periodic alignment." << endl;
                return;
            }

            CompositeSharedPtr c1 = it1->second;
            CompositeSharedPtr c2 = it2->second;

            if (c1->items.size() != c2->items.size())
            {
                cerr << "WARNING: Surfaces " << surf1 << " and " << surf2
                     << " have different numbers of elements. Skipping periodic"
                     << " alignment." << endl;
                return;
            }

            c1->reorder = false;
            c2->reorder = false;

            map<int, pair<FaceSharedPtr, vector<int> > > perFaces;

            // Loop over elements, calculate centroids of elements in c2.
            map<int, Node> centroidMap;
            map<int, Node>::iterator it;
            for (int i = 0; i < c2->items.size(); ++i)
            {
                Node centroid;
                for (int j = 0; j < c2->items[i]->GetVertexCount(); ++j)
                {
                    centroid += *(c2->items[i]->GetVertex(j));
                }
                centroid /= (NekDouble)c2->items[i]->GetVertexCount();
                centroidMap[i] = centroid;
            }

            boost::unordered_set<int> elmtDone;
            map<int, int> elmtPairs;
            map<int, int> vertCheck;

            for (int i = 0; i < c1->items.size(); ++i)
            {
                Node centroid;
                for (int j = 0; j < c1->items[i]->GetVertexCount(); ++j)
                {
                    centroid += *(c1->items[i]->GetVertex(j));
                }
                centroid /= (NekDouble)c1->items[i]->GetVertexCount();

                for (it = centroidMap.begin(); it != centroidMap.end(); ++it)
                {
                    if (elmtDone.count(it->first) > 0)
                    {
                        continue;
                    }

                    Node dx = it->second - centroid;
                    if (fabs(fabs(dx.x*vec[0] + dx.y*vec[1] + dx.z*vec[2])/
                             sqrt(dx.abs2()) - 1.0) < 1e-6)
                    {
                        // Found match
                        int id1 = c1->items[i]        ->GetFaceLink()->id;
                        int id2 = c2->items[it->first]->GetFaceLink()->id;

                        elmtDone.insert(it->first);
                        elmtPairs[i] = it->first;

                        // Identify periodic vertices
                        int nVerts = c1->items[i]->GetVertexCount();
                        vector<int> perVerts(nVerts, 0), perVertsInv(nVerts, 0);

                        if (orient)
                        {
                            for (int k = 0; k < nVerts; ++k)
                            {
                                NodeSharedPtr n1 = c1->items[i]->GetFaceLink()->vertexList[k];
                                int l;
                                
                                for (l = 0; l < nVerts; ++l)
                                {
                                    NodeSharedPtr n2 =
                                        c2->items[it->first]->GetFaceLink()->vertexList[l];
                                    
                                    Node dn = *n2 - *n1;
                                    if (fabs(fabs(dn.x*vec[0] + dn.y*vec[1] +
                                                  dn.z*vec[2])/
                                             sqrt(dn.abs2()) - 1.0) < 1e-6)
                                    {
                                        perVerts   [k] = l;
                                        perVertsInv[l] = k;
                                        
                                        int id1 = n1->id;
                                        int id2 = n2->id;
                                        if (vertCheck.count(id1) == 0)
                                        {
                                            vertCheck[id1] = id2;
                                        }
                                        else
                                        {
                                            ASSERTL0(vertCheck[id1] == id2,
                                                     "Periodic vertex already "
                                                     "identified!");
                                        }
                                        break;
                                    }
                                }
                                ASSERTL1(l < nVerts,
                                         "Could not identify periodic vertices.");
                            }
                            
                            int tot1 = 0, tot2 = 0;
                            for (int k = 0; k < nVerts; ++k)
                            {
                                tot1 += perVerts   [k];
                                tot2 += perVertsInv[k];
                            }
                            ASSERTL0(tot1 == nVerts*(nVerts-1)/2 &&
                                     tot2 == nVerts*(nVerts-1)/2,
                                     "Error identifying periodic vertices");
                        }

                        perFaces[id1] = make_pair(
                            c2->items[it->first]->GetFaceLink(), perVerts);
                        perFaces[id2] = make_pair(
                            c1->items[i]        ->GetFaceLink(), perVertsInv);
                        break;
                    }
                }

                if (it == centroidMap.end())
                {
                    cerr << "WARNING: Could not find matching edge for surface "
                         << "element " << c1->items[i]->GetId() << ". "
                         << "Skipping periodic alignment." << endl;
                    return;
                }
            }

            // Reorder vectors.
            vector<ElementSharedPtr> tmp = c2->items;

            map<int, int>::iterator mIt;

            for (int i = 0; i < tmp.size(); ++i)
            {
                c2->items[i] = tmp[elmtPairs[i]];
            }

            if (orient)
            {
                ReorderPrisms(perFaces);
            }
        }
    }
}
