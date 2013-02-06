////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessTetSplit.cpp
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
//  Description: Split prisms -> tets
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessTetSplit.h"

#include <LocalRegions/PrismExp.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessTetSplit::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "tetsplit"), ProcessTetSplit::create,
                "Split prismatic elements to tetrahedra");
        
        ProcessTetSplit::ProcessTetSplit(MeshSharedPtr m) : ProcessModule(m)
        {
	    config["nq"]   = ConfigOption(false, "3",
                "Number of points in high order elements.");
        }
      
        ProcessTetSplit::~ProcessTetSplit()
        {
          
        }
        
        void ProcessTetSplit::Process()
        {
            int nodeId = m->vertexSet.size();
            
            // Set up map which identifies edges (as pairs of vertex ids)
            // including their vertices to the offset/stride in the 3d array of
            // collapsed co-ordinates. Note that this map also includes the
            // diagonal edges along quadrilateral faces which will be used to
            // add high-order information to the split tetrahedra.
            map<pair<int,int>, pair<int,int> > edgeMap;
            map<pair<int,int>, pair<int,int> >::iterator it;
            
            // Standard prismatic edges (0->8)
            int nq = config["nq"].as<int>();
            edgeMap[pair<int,int>(0,1)] = pair<int,int>(0,            1);
            edgeMap[pair<int,int>(3,2)] = pair<int,int>(nq*(nq-1),    1);
            edgeMap[pair<int,int>(0,3)] = pair<int,int>(0,            nq);
            edgeMap[pair<int,int>(1,2)] = pair<int,int>(nq-1,         nq);
	    edgeMap[pair<int,int>(4,5)] = pair<int,int>((nq-1)*nq*nq, nq);
            edgeMap[pair<int,int>(0,4)] = pair<int,int>(0,            nq*nq);
            edgeMap[pair<int,int>(1,4)] = pair<int,int>(nq-1,         nq*nq);
            edgeMap[pair<int,int>(2,5)] = pair<int,int>(nq*nq-1,      nq*nq);
            edgeMap[pair<int,int>(3,5)] = pair<int,int>(nq*(nq-1),    nq*nq);
            // Face 0 diagonals
            edgeMap[pair<int,int>(0,2)] = pair<int,int>(0,            nq+1);
            edgeMap[pair<int,int>(1,3)] = pair<int,int>(nq-1,         nq-1);
            // Face 2 diagonals
            edgeMap[pair<int,int>(1,5)] = pair<int,int>(nq-1,         nq*nq+nq);
            edgeMap[pair<int,int>(2,4)] = pair<int,int>(nq*nq-1,      nq*nq-nq);
            // Face 4 diagonals
            edgeMap[pair<int,int>(0,5)] = pair<int,int>(0,            nq*nq+nq);
            edgeMap[pair<int,int>(3,4)] = pair<int,int>(nq*(nq-1),    nq*nq-nq);

            // Default PointsType.
            LibUtilities::PointsType pt = LibUtilities::eGaussLobattoLegendre;
            
            /*
             * Split all element types into tetrahedra. This is based on the
             * algorithm found in:
             * 
             * "How to Subdivide Pyramids, Prisms and Hexahedra into
             * Tetrahedra", J. Dompierre et al.
             */

            // Denotes a set of indices inside m->element[m->expDim-1] which are
            // to be removed. These are precisely the quadrilateral boundary
            // faces which will be replaced by two triangular faces.
            set<int> toRemove;
                
            // Represents table 2 of paper; each row i represents a rotation of
            // the prism nodes such that vertex i is placed at position 0.
            static int indir[6][6] = {
                {0,1,2,3,4,5},
                {1,2,0,4,5,3},
                {2,0,1,5,3,4},
                {3,5,4,0,2,1},
                {4,3,5,1,0,2},
                {5,4,3,2,1,0}
            };
            
            // Represents table 3 of paper; the first three rows are the three
            // tetrahedra if the first condition is met; the latter three rows
            // are the three tetrahedra if the second condition is met.
            static int prismTet[6][4] = {
                {0,1,2,5},
                {0,1,5,4},
                {0,4,5,3},
                {0,1,2,4},
                {0,4,2,5},
                {0,4,5,3}
            };
                
            // Represents the order of tetrahedral edges (in Nektar++ ordering).
            static int tetEdges[6][2] = {
                {0,1}, {1,2},
                {0,2}, {0,3},
                {1,3}, {2,3}};
            
            // A tetrahedron nodes -> faces map.
            static int tetFaceNodes[4][3] = {
                {0,1,2},{0,1,3},{1,2,3},{0,2,3}};
            
            // Make a copy of the element list.
            vector<ElementSharedPtr> el = m->element[m->expDim];
            m->element[m->expDim].clear();
            
            for (int i = 0; i < el.size(); ++i)
            {
                if (el[i]->GetConf().e != ePrism)
                {
                    m->element[m->expDim].push_back(el[i]);
                    continue;
                }                    
                
                vector<NodeSharedPtr> nodeList(6);
                
                // Map Nektar++ ordering (vertices 0,1,2,3 are base quad) to
                // paper ordering (vertices 0,1,2 are first triangular face).
                int mapPrism[6] = {0,1,4,3,2,5};
                for (int j = 0; j < 6; ++j)
                {
                    nodeList[j] = el[i]->GetVertex(mapPrism[j]);
                }
                    
                // Determine minimum ID of the nodes in this prism.
                int minElId = nodeList[0]->id;
                int minId   = 0;
                for (int j = 1; j < 6; ++j)
                {
                    int curId = nodeList[j]->id;
                    if (curId < minElId)
                    {
                        minElId = curId;
                        minId   = j;
                    }
                }
                
                int offset;
                
                // Split prism using paper criterion.
                if (min(nodeList[indir[minId][1]]->id, nodeList[indir[minId][5]]->id) <
                    min(nodeList[indir[minId][2]]->id, nodeList[indir[minId][4]]->id))
                {
                    offset = 0;
                }
                else if (min(nodeList[indir[minId][1]]->id, nodeList[indir[minId][5]]->id) >
                         min(nodeList[indir[minId][2]]->id, nodeList[indir[minId][4]]->id))
                {
                    offset = 3;
                }
                else
                {
                    cerr << "Connectivity issue with prism->tet splitting." << endl;
                    abort();
                }
                
                // Create local prismatic region so that co-ordinates of the
                // mapped element can be read from.
                SpatialDomains::PrismGeomSharedPtr geomLayer = 
                    boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                        el[i]->GetGeom(m->spaceDim));
                LibUtilities::BasisKey B0(LibUtilities::eModified_A, nq,
                                          LibUtilities::PointsKey(nq,pt));
                LibUtilities::BasisKey B1(LibUtilities::eModified_B, nq,
                                          LibUtilities::PointsKey(nq,pt));
                LocalRegions::PrismExpSharedPtr qs = 
                    MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                        B0,B0,B1,geomLayer);
                    
                // Get the coordiantes of the high order prismatic element.
                Array<OneD, NekDouble> xs(nq*nq*nq);
                Array<OneD, NekDouble> ys(nq*nq*nq);
                Array<OneD, NekDouble> zs(nq*nq*nq);
                qs->GetCoords(xs,ys,zs);
                
                for (int j = 0; j < 3; ++j)
                {
                    vector<NodeSharedPtr> tetNodes(4);
                    
                    for (int k = 0; k < 4; ++k)
                    {
                        tetNodes[k] = nodeList[indir[minId][prismTet[j+offset][k]]];
                    }
                    
                    // Add high order information to tetrahedral edges.
                    for (int k = 0; k < 6; ++k)
                    {
                        // Determine prismatic nodes which correspond with this
                        // edge. Apply prism map to this to get Nektar++
                        // ordering (this works since as a permutation,
                        // prismMap^{-1} = prismMap).
                        int n1 = mapPrism[indir[minId][prismTet[j+offset][tetEdges[k][0]]]];
                        int n2 = mapPrism[indir[minId][prismTet[j+offset][tetEdges[k][1]]]];
                        
                        // Find offset/stride
                        it = edgeMap.find(pair<int,int>(n1,n2));
                        if (it == edgeMap.end())
                        {
                            it = edgeMap.find(pair<int,int>(n2,n1));
                            if (it == edgeMap.end())
                            {
                                cerr << "Couldn't find prism edges " << n1 
                                     << " " << n2 << endl;
                                abort();
                            }
                            // Extract vertices -- reverse order.
                            for (int l = nq-2; l >= 1; --l)
                            {
                                int pos = it->second.first + l*it->second.second;
                                tetNodes.push_back(
                                    NodeSharedPtr(
                                        new Node(nodeId++, xs[pos], ys[pos], zs[pos])));
                            }
                        }
                        else
                        {
                            // Extract vertices -- forwards order.
                            for (int l = 1; l < nq-1; ++l)
                            {
                                int pos = it->second.first + l*it->second.second;
                                tetNodes.push_back(
                                    NodeSharedPtr(
                                        new Node(nodeId++, xs[pos], ys[pos], zs[pos])));
                            }
                        }
                    }
                    
                    // Create new tetrahedron with edge curvature. TODO: face
                    // nodes.
                    vector<int> tags = el[i]->GetTagList();
                    ElmtConfig conf(eTetrahedron, nq-1, false, false);
                    ElementSharedPtr elmt = GetElementFactory().
                        CreateInstance(eTetrahedron,conf,tetNodes,tags);
                    
                    m->element[m->expDim].push_back(elmt);
                }
                
                // Now check to see if this one of the quadrilateral faces is
                // associated with a boundary condition. If it is, we split the
                // face into two triangles and mark the existing face for
                // removal.
                //
                // Note that this algorithm has significant room for improvement
                // and is likely one of the least optimal approachs - however
                // implementation is simple.
                for (int fid = 0; fid < 5; fid += 2)
                {
                    int bl = el[i]->GetBoundaryLink(fid);
                    
                    if (bl == -1)
                    {
                        continue;
                    }
                    
                    vector<NodeSharedPtr> triNodeList(3);
                    vector<int>           faceNodes  (3);
                    vector<int>           tmp;
                    vector<int>           tagBE;
                    ElmtConfig            bconf(eTriangle, 1, true, true);
                    ElementSharedPtr      elmt;
                    
                    // Mark existing boundary face for removal.
                    toRemove.insert(bl);
                    tagBE =  m->element[m->expDim-1][bl]->GetTagList();
                    
                    // First loop over tets.
                    for (int j = 0; j < 3; ++j)
                    {
                        // Now loop over faces.
                        for (int k = 0; k < 4; ++k)
                        {
                            // Finally determine Nektar++ local node
                            // numbers for this face.
                            for (int l = 0; l < 3; ++l)
                            {
                                faceNodes[l] = mapPrism[indir[minId][prismTet[j+offset][tetFaceNodes[k][l]]]];
                            }
                            
                            tmp = faceNodes;
                            sort(faceNodes.begin(), faceNodes.end());
                            
                            // If this face matches a triple denoting a
                            // split quad face, add the face to the
                            // expansion list.
                            if ((fid == 0 && (
                                     (faceNodes[0] == 0 && faceNodes[1] == 1 && faceNodes[2] == 2)   ||
                                     (faceNodes[0] == 0 && faceNodes[1] == 2 && faceNodes[2] == 3)   ||
                                     (faceNodes[0] == 0 && faceNodes[1] == 1 && faceNodes[2] == 3)   ||
                                     (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 3))) ||
                                (fid == 2 && (
                                    (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 5)   ||
                                    (faceNodes[0] == 1 && faceNodes[1] == 4 && faceNodes[2] == 5)   ||
                                    (faceNodes[0] == 1 && faceNodes[1] == 2 && faceNodes[2] == 4)   ||
                                    (faceNodes[0] == 2 && faceNodes[1] == 4 && faceNodes[2] == 5))) ||
                                (fid == 4 && (
                                    (faceNodes[0] == 0 && faceNodes[1] == 3 && faceNodes[2] == 5) ||
                                    (faceNodes[0] == 0 && faceNodes[1] == 4 && faceNodes[2] == 5) ||
                                    (faceNodes[0] == 0 && faceNodes[1] == 3 && faceNodes[2] == 4) ||
                                    (faceNodes[0] == 3 && faceNodes[1] == 4 && faceNodes[2] == 5))))
                            {
                                triNodeList[0] = nodeList[mapPrism[tmp[0]]];
                                triNodeList[1] = nodeList[mapPrism[tmp[1]]];
                                triNodeList[2] = nodeList[mapPrism[tmp[2]]];
                                elmt           = GetElementFactory().
                                    CreateInstance(eTriangle,bconf,triNodeList,tagBE);
                                m->element[m->expDim-1].push_back(elmt);
                            }
                        }
                    }
                }
            }
            
            // Remove 2D elements.
            vector<ElementSharedPtr> tmp;
            for (int i = 0; i < m->element[m->expDim-1].size(); ++i)
            {
                set<int>::iterator it = toRemove.find(i);
                if (it == toRemove.end())
                {
                    tmp.push_back(m->element[m->expDim-1][i]);
                }
            }
            
            m->element[m->expDim-1] = tmp;
            
            // Re-process mesh to eliminate duplicate vertices and edges.
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
