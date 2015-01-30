////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBL.cpp
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
//  Description: Refine prismatic boundary layer elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
using namespace std;

#include "MeshElements.h"
#include "ProcessBL.h"

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/HexExp.h>

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey ProcessBL::className =
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eProcessModule, "bl"), ProcessBL::create,
                "Refines a prismatic boundary layer.");

        int **helper2d(int lda, int arr[][2])
        {
            int **ret = new int*[lda];
            for (int i = 0; i < lda; ++i)
            {
                ret[i] = new int[2];
                ret[i][0] = arr[i][0];
                ret[i][1] = arr[i][1];
            }
            return ret;
        }

        int **helper2d(int lda, int arr[][4])
        {
            int **ret = new int*[lda];
            for (int i = 0; i < lda; ++i)
            {
                ret[i] = new int[4];
                ret[i][0] = arr[i][0];
                ret[i][1] = arr[i][1];
                ret[i][2] = arr[i][2];
                ret[i][3] = arr[i][3];
            }
            return ret;
        }

        struct SplitMapHelper
        {
            int size;
            int layerOff;
            int *edge;
            int *offset;
            int *inc;
            int **conn;
            int bfacesSize;
            int *bfaces;
        };

        struct SplitEdgeHelper
        {
            int size;
            int *edge;
            int **edgeVert;
            int *offset;
            int *inc;
        };

        ProcessBL::ProcessBL(MeshSharedPtr m) : ProcessModule(m)
        {
            // BL mesh configuration.
            m_config["layers"]     = ConfigOption(false, "2",
                "Number of layers to refine.");
            m_config["nq"]         = ConfigOption(false, "5",
                "Number of points in high order elements.");
            m_config["surf"]       = ConfigOption(false, "",
                "Tag identifying surface connected to prism.");
            m_config["r"]          = ConfigOption(false, "2.0",
                "Ratio to use in geometry progression.");
        }

        ProcessBL::~ProcessBL()
        {

        }

        void ProcessBL::Process()
        {
            if (m_mesh->m_verbose)
            {
                cout << "ProcessBL: Refining prismatic boundary layer..."
                     << endl;
            }

            // A set containing all element types which are valid.
            set<LibUtilities::ShapeType> validElTypes;
            validElTypes.insert(LibUtilities::ePrism);
            validElTypes.insert(LibUtilities::eHexahedron);

            int nodeId  = m_mesh->m_vertexSet.size();
            int nl      = m_config["layers"].as<int>();
            int nq      = m_config["nq"].    as<int>();

            // determine if geometric ratio is string or a constant.
            LibUtilities::AnalyticExpressionEvaluator rEval;
            NekDouble r             =  1;
            int       rExprId       = -1;
            bool      ratioIsString = false;

            if (m_config["r"].isType<NekDouble>())
            {
                r = m_config["r"].as<NekDouble>();
            }
            else
            {
                std::string rstr = m_config["r"].as<string>();
                rExprId = rEval.DefineFunction("x y z", rstr);
                ratioIsString = true;
            }

            // Prismatic node -> face map.
            int prismFaceNodes[5][4] = {
                {0,1,2,3},{0,1,4,-1},{1,2,5,4},{3,2,5,-1},{0,3,5,4}};
            int hexFaceNodes  [6][4] = {
                {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7}};
            map<LibUtilities::ShapeType, int **> faceNodeMap;
            faceNodeMap[LibUtilities::ePrism]      = helper2d(5, prismFaceNodes);
            faceNodeMap[LibUtilities::eHexahedron] = helper2d(6, hexFaceNodes);

            // Default PointsType.
            LibUtilities::PointsType pt = LibUtilities::eGaussLobattoLegendre;

            // Map which takes element ID to face on surface. This enables
            // splitting to occur in either y-direction of the prism.
            boost::unordered_map<int, int> splitEls;
            boost::unordered_map<int, int>::iterator sIt;

            // Set up maps which takes an edge (in nektar++ ordering) and return
            // their offset and stride in the 3d array of collapsed quadrature
            // points. Note that this map includes only the edges that are on
            // the triangular faces as the edges in the normal direction are
            // linear.
            map<LibUtilities::ShapeType, map<int, SplitMapHelper> > splitMap;
            int po = nq*(nl+1);

            SplitMapHelper splitPrism;
            int splitMapEdgePrism  [6]    = {0, 2,  4,  5,    6,       7};
            int splitMapOffsetPrism[6]    = {0, nq, 0,  nq-1, nq+nq-1, nq};
            int splitMapIncPrism   [6]    = {1, 1,  po, po,   po,      po};
            int splitMapBFacesPrism[3]    = {0, 2, 4};
            int splitMapConnPrism  [6][2] = {{0,0}, {1,0}, {1,1},
                                             {0,1}, {2,0}, {2,1}};
            splitPrism.size       = 6;
            splitPrism.layerOff   = nq;
            splitPrism.edge       = splitMapEdgePrism;
            splitPrism.offset     = splitMapOffsetPrism;
            splitPrism.inc        = splitMapIncPrism;
            splitPrism.conn       = helper2d(6, splitMapConnPrism);
            splitPrism.bfacesSize = 3;
            splitPrism.bfaces     = splitMapBFacesPrism;
            splitMap[LibUtilities::ePrism][1] = splitPrism;
            splitMap[LibUtilities::ePrism][3] = splitPrism;

            int ho = nq*(nq-1);
            int tl = nq*nq;
            SplitMapHelper splitHex0;
            int splitMapEdgeHex0  [8]    = {0, 1,    2,     3,   8,  9,       10,     11};
            int splitMapOffsetHex0[8]    = {0, nq-1, tl-1,  ho,  tl, tl+nq-1, 2*tl-1, tl+ho};
            int splitMapIncHex0   [8]    = {1, nq,   -1,   -nq,  1,  nq,      -1,     -nq};
            int splitMapBFacesHex0[4]    = {1, 2, 3, 4};
            int splitMapConnHex0  [8][2] = {{0,0}, {1,0}, {2,0}, {3,0},
                                            {0,1}, {1,1}, {2,1}, {3,1}};
            splitHex0.size       = 8;
            splitHex0.layerOff   = nq*nq;
            splitHex0.edge       = splitMapEdgeHex0;
            splitHex0.offset     = splitMapOffsetHex0;
            splitHex0.inc        = splitMapIncHex0;
            splitHex0.conn       = helper2d(8, splitMapConnHex0);
            splitHex0.bfacesSize = 4;
            splitHex0.bfaces     = splitMapBFacesHex0;
            splitMap[LibUtilities::eHexahedron][0] = splitHex0;
            splitMap[LibUtilities::eHexahedron][5] = splitHex0;

            // splitEdge enumerates the edges in the standard prism along which
            // new nodes should be generated. These edges are the three between
            // the two triangular faces.
            //
            // edgeVertMap specifies the vertices which comprise those edges in
            // splitEdge; for example splitEdge[0] = 3 which connects vertices 0
            // and 3.
            //
            // edgeOffset holds the offset of each of edges 3, 1 and 8
            // respectively inside the collapsed coordinate system.
            map<LibUtilities::ShapeType, map<int, SplitEdgeHelper> > splitEdge;

            int splitPrismEdges   [3]    = {3,     1,     8};
            int splitPrismEdgeVert[3][2] = {{0,3}, {1,2}, {4,5}};
            int splitPrismOffset  [3]    = {0,     nq-1,  nq*(nl+1)*(nq-1)};
            int splitPrismInc     [3]    = {nq,    nq,    nq};
            SplitEdgeHelper splitPrismEdge;
            splitPrismEdge.size     = 3;
            splitPrismEdge.edge     = splitPrismEdges;
            splitPrismEdge.edgeVert = helper2d(3, splitPrismEdgeVert);
            splitPrismEdge.offset   = splitPrismOffset;
            splitPrismEdge.inc      = splitPrismInc;
            splitEdge[LibUtilities::ePrism][1] = splitPrismEdge;
            splitEdge[LibUtilities::ePrism][3] = splitPrismEdge;

            int splitHex0Edges   [4]    = {4,     5,     6,       7};
            int splitHex0EdgeVert[4][2] = {{0,4}, {1,5}, {2,6},   {3,7}};
            int splitHex0Offset  [4]    = {0,     nq-1,  nq*nq-1, nq*(nq-1) };
            int splitHex0Inc     [4]    = {nq*nq, nq*nq, nq*nq,   nq*nq};
            SplitEdgeHelper splitHex0Edge;
            splitHex0Edge.size     = 4;
            splitHex0Edge.edge     = splitHex0Edges;
            splitHex0Edge.edgeVert = helper2d(4, splitHex0EdgeVert);
            splitHex0Edge.offset   = splitHex0Offset;
            splitHex0Edge.inc      = splitHex0Inc;
            splitEdge[LibUtilities::eHexahedron][0] = splitHex0Edge;
            splitEdge[LibUtilities::eHexahedron][5] = splitHex0Edge;

            map<LibUtilities::ShapeType, map<int, bool> > revPoints;
            revPoints[LibUtilities::ePrism][1] = false;
            revPoints[LibUtilities::ePrism][3] = true;

            revPoints[LibUtilities::eHexahedron][0] = true;
            revPoints[LibUtilities::eHexahedron][5] = false;

            // edgeMap associates geometry edge IDs to the (nl+1) vertices which
            // are generated along that edge when a prism is split, and is used
            // to avoid generation of duplicate vertices. It is stored as an
            // unordered map for speed.
            boost::unordered_map<int, vector<NodeSharedPtr> > edgeMap;
            boost::unordered_map<int, vector<NodeSharedPtr> >::iterator eIt;

            string surf = m_config["surf"].as<string>();
            if (surf.size() > 0)
            {
                vector<unsigned int> surfs;
                ParseUtils::GenerateSeqVector(surf.c_str(), surfs);
                sort(surfs.begin(), surfs.end());

                // If surface is defined, process list of elements to find those
                // that are connected to it.
                for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
                {
                    ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
                    int nSurf = el->GetFaceCount();

                    for (int j = 0; j < nSurf; ++j)
                    {
                        int bl = el->GetBoundaryLink(j);
                        if (bl == -1)
                        {
                            continue;
                        }

                        ElementSharedPtr bEl  = m_mesh->m_element[m_mesh->m_expDim-1][bl];
                        vector<int>      tags = bEl->GetTagList();
                        vector<int>      inter;

                        sort(tags.begin(), tags.end());
                        set_intersection(surfs.begin(), surfs.end(),
                                         tags .begin(), tags .end(),
                                         back_inserter(inter));
                        ASSERTL0(inter.size() <= 1,
                                 "Intersection of surfaces wrong");

                        if (inter.size() == 1)
                        {
                            if (el->GetConf().m_e == LibUtilities::ePrism)
                            {
                                if (j % 2 == 0)
                                {
                                    cerr << "WARNING: Found quadrilateral face "
                                         << j << " on surface " << surf
                                         << " connected to prism; ignoring."
                                         << endl;
                                    continue;
                                }

                                if (splitEls.count(el->GetId()) > 0)
                                {
                                    cerr << "WARNING: prism already found; "
                                         << "ignoring" << endl;
                                }

                                splitEls[el->GetId()] = j;
                            }
                            else if (validElTypes.count(el->GetConf().m_e) == 0)
                            {
                                cerr << "WARNING: Unsupported element type "
                                     << "found in surface " << j << "; "
                                     << "ignoring" << endl;
                                continue;
                            }
                        }
                    }
                }
            }
            else
            {
                // Otherwise, add all prismatic elements and assume face 1 of
                // the prism lies on the surface.
                for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
                {
                    ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];

                    if (el->GetConf().m_e == LibUtilities::ePrism)
                    {
                        splitEls[el->GetId()] = 1;
                    }
                    else if (validElTypes.count(el->GetConf().m_e) > 0)
                    {
                        splitEls[el->GetId()] = 0;
                    }
                    else
                    {
                        continue;
                    }
                }
            }

            if (splitEls.size() == 0)
            {
                cerr << "WARNING: No elements detected to split." << endl;
                return;
            }

            // Erase all elements from the element list. Elements will be
            // re-added as they are split.
            vector<ElementSharedPtr> el = m_mesh->m_element[m_mesh->m_expDim];
            m_mesh->m_element[m_mesh->m_expDim].clear();

            // Iterate over list of elements of expansion dimension.
            for (int i = 0; i < el.size(); ++i)
            {
                const int elId = el[i]->GetId();
                sIt = splitEls.find(elId);

                if (sIt == splitEls.end())
                {
                    m_mesh->m_element[m_mesh->m_expDim].push_back(el[i]);
                    continue;
                }

                const int faceNum = sIt->second;
                LibUtilities::ShapeType elType = el[i]->GetConf().m_e;

                SplitMapHelper  &sMap  = splitMap [elType][faceNum];
                SplitEdgeHelper &sEdge = splitEdge[elType][faceNum];

                // Find quadrilateral boundary faces if any
                std::map<int, int> bLink;
                for (int j = 0; j < sMap.bfacesSize; ++j)
                {
                    int bl = el[i]->GetBoundaryLink(sMap.bfaces[j]);
                    if (bl != -1)
                    {
                        bLink[sMap.bfaces[j]] = bl;
                    }
                }

                // Get elemental geometry object.
                SpatialDomains::Geometry3DSharedPtr geom =
                    boost::dynamic_pointer_cast<SpatialDomains::Geometry3D>(
                        el[i]->GetGeom(m_mesh->m_spaceDim));

                // Determine whether to use reverse points.
                LibUtilities::PointsType t =
                    revPoints[elType][faceNum] ?
                    LibUtilities::eBoundaryLayerPoints :
                    LibUtilities::eBoundaryLayerPointsRev;

                // Determine value of r based on geometry.
                if(ratioIsString)
                {
                    NekDouble x,  y,  z;
                    NekDouble x1, y1, z1;
                    int nverts = geom->GetNumVerts();

                    x = y = z = 0.0;

                    for(int i = 0; i < nverts; ++i)
                    {
                        geom->GetVertex(i)->GetCoords(x1,y1,z1);
                        x += x1; y += y1; z += z1;
                    }
                    x /= (NekDouble) nverts;
                    y /= (NekDouble) nverts;
                    z /= (NekDouble) nverts;
                    r = rEval.Evaluate(rExprId,x,y,z,0.0);
                }

                LocalRegions::ExpansionSharedPtr q;

                if (elType == LibUtilities::ePrism)
                {
                    // Create basis.
                    LibUtilities::BasisKey B0(
                        LibUtilities::eModified_A, nq,
                        LibUtilities::PointsKey(nq,pt));
                    LibUtilities::BasisKey B1(
                        LibUtilities::eModified_A, 2,
                        LibUtilities::PointsKey(nl+1, t, r));
                    LibUtilities::BasisKey B2(
                        LibUtilities::eModified_B, nq,
                        LibUtilities::PointsKey(nq,pt));

                    // Create local region.
                    SpatialDomains::PrismGeomSharedPtr g =
                        boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                            geom);
                    q = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                        B0, B1, B2, g);
                }
                else if (elType == LibUtilities::eHexahedron)
                {
                    // Create basis.
                    LibUtilities::BasisKey B0(
                        LibUtilities::eModified_A, nq,
                        LibUtilities::PointsKey(nq,pt));
                    LibUtilities::BasisKey B1(
                        LibUtilities::eModified_A, 2,
                        LibUtilities::PointsKey(nl+1, t, r));

                    // Create local region.
                    SpatialDomains::HexGeomSharedPtr g =
                        boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(
                            geom);
                    q = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(
                        B0, B0, B1, g);
                }

                // Grab co-ordinates.
                Array<OneD, NekDouble> x(nq*nq*(nl+1));
                Array<OneD, NekDouble> y(nq*nq*(nl+1));
                Array<OneD, NekDouble> z(nq*nq*(nl+1));
                q->GetCoords(x,y,z);

                int nSplitEdge = sEdge.size;
                vector<vector<NodeSharedPtr> > edgeNodes(nSplitEdge);

                // Loop over edges to be split.
                for (int j = 0; j < nSplitEdge; ++j)
                {
                    int locEdge = sEdge.edge[j];
                    int edgeId  = el[i]->GetEdge(locEdge)->m_id;

                    // Determine whether we have already generated vertices
                    // along this edge.
                    eIt = edgeMap.find(edgeId);

                    if (eIt == edgeMap.end())
                    {
                        // If not then resize storage to hold new points.
                        edgeNodes[j].resize(nl+1);

                        // Re-use existing vertices at endpoints of edge to
                        // avoid duplicating the existing vertices.
                        edgeNodes[j][0]  = el[i]->GetVertex(sEdge.edgeVert[j][0]);
                        edgeNodes[j][nl] = el[i]->GetVertex(sEdge.edgeVert[j][1]);

                        // Variable geometric ratio
                        if(ratioIsString)
                        {
                            NekDouble x0,y0,z0;
                            NekDouble x1,y1,z1;
                            NekDouble xm,ym,zm;

                            // -> Find edge end and mid points
                            x0 = x[sEdge.offset[j]];
                            y0 = y[sEdge.offset[j]];
                            z0 = z[sEdge.offset[j]];

                            x1 = x[sEdge.offset[j]+nl*nq];
                            y1 = y[sEdge.offset[j]+nl*nq];
                            z1 = z[sEdge.offset[j]+nl*nq];

                            xm = 0.5*(x0+x1);
                            ym = 0.5*(y0+y1);
                            zm = 0.5*(z0+z1);

                            // evaluate r factor based on mid point value
                            NekDouble rnew;
                            rnew = rEval.Evaluate(rExprId,xm,ym,zm,0.0);

                            // Get basis with new r;
                            LibUtilities::PointsKey Pkey(nl+1, t, rnew);
                            LibUtilities::PointsSharedPtr newP
                                = LibUtilities::PointsManager()[Pkey];

                            const Array<OneD, const NekDouble> z = newP->GetZ();

                            // Create new interior nodes based on this new blend
                            for (int k = 1; k < nl; ++k)
                            {
                                xm = 0.5*(1+z[k])*(x1-x0) + x0;
                                ym = 0.5*(1+z[k])*(y1-y0) + y0;
                                zm = 0.5*(1+z[k])*(z1-z0) + z0;
                                edgeNodes[j][k] = NodeSharedPtr(
                                        new Node(nodeId++, xm,ym,zm));
                            }
                        }
                        else
                        {
                            // Create new interior nodes.
                            for (int k = 1; k < nl; ++k)
                            {
                                int pos = sEdge.offset[j] + k*sEdge.inc[j];
                                edgeNodes[j][k] = NodeSharedPtr(
                                    new Node(nodeId++, x[pos], y[pos], z[pos]));
                            }
                        }

                        // Store these edges in edgeMap.
                        edgeMap[edgeId] = edgeNodes[j];
                    }
                    else
                    {
                        edgeNodes[j] = eIt->second;
                    }
                }

                // Create element layers.
                for (int j = 0; j < nl; ++j)
                {
                    // Offset of this layer within the collapsed coordinate
                    // system.
                    int offset = j * sMap.layerOff;

                    // Get corner vertices.
                    vector<NodeSharedPtr> nodeList(sMap.size);
                    for (int k = 0; k < sMap.size; ++k)
                    {
                        nodeList[k] =
                            edgeNodes[sMap.conn[k][0]][j + sMap.conn[k][1]];
                    }

                    // Create the element.
                    ElmtConfig conf(elType, 1, true, true, false);
                    ElementSharedPtr elmt = GetElementFactory().
                        CreateInstance(
                            elType, conf, nodeList, el[i]->GetTagList());

                    // Add high order nodes to split prismatic edges.
                    for (int l = 0; l < sMap.size; ++l)
                    {
                        EdgeSharedPtr HOedge = elmt->GetEdge(
                            sMap.edge[l]);
                        for (int k = 1; k < nq-1; ++k)
                        {
                            int pos = offset + sMap.offset[l] + k*sMap.inc[l];
                            HOedge->m_edgeNodes.push_back(
                                NodeSharedPtr(
                                    new Node(nodeId++,x[pos],y[pos],z[pos])));
                        }
                        HOedge->m_curveType = pt;
                    }

                    // Change the surface elements to match the layers of
                    // elements on the boundary of the domain.
                    map<int,int>::iterator it;
                    for (it = bLink.begin(); it != bLink.end(); ++it)
                    {
                        int fid = it->first;
                        int bl  = it->second;

                        if (j == 0)
                        {
                            // For first layer reuse existing 2D element.
                            ElementSharedPtr e = m_mesh->m_element[m_mesh->m_expDim-1][bl];
                            for (int k = 0; k < 4; ++k)
                            {
                                e->SetVertex(
                                    k, nodeList[faceNodeMap[elType][fid][k]]);
                            }
                        }
                        else
                        {
                            // For all other layers create new element.
                            vector<NodeSharedPtr> qNodeList(4);
                            for (int k = 0; k < 4; ++k)
                            {
                                qNodeList[k] = nodeList[faceNodeMap[elType][fid][k]];
                            }
                            vector<int> tagBE;
                            tagBE = m_mesh->m_element[m_mesh->m_expDim-1][bl]->GetTagList();
                            ElmtConfig bconf(LibUtilities::eQuadrilateral,1,true,true,false);
                            ElementSharedPtr boundaryElmt = GetElementFactory().
                                CreateInstance(LibUtilities::eQuadrilateral,bconf,
                                               qNodeList,tagBE);
                            m_mesh->m_element[m_mesh->m_expDim-1].push_back(boundaryElmt);
                        }
                    }

                    m_mesh->m_element[m_mesh->m_expDim].push_back(elmt);
                }
            }

            // Re-process mesh to eliminate duplicate vertices and edges.
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
