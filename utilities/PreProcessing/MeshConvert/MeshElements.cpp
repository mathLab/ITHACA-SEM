////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshElements.cpp
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <loki/Singleton.h>

#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdNodalTetExp.h>
#include <StdRegions/StdNodalPrismExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>

#include "MeshElements.h"
using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        ElementFactory& GetElementFactory()
        {
            typedef Loki::SingletonHolder<ElementFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        Element::Element(ElmtConfig   pConf,
                         unsigned int pNumNodes,
                         unsigned int pGotNodes) 
            : m_conf(pConf), m_geom()
        {
            if (pNumNodes != pGotNodes)
            {
                cerr << "Number of modes mismatch for type " 
                     << pConf.e << "! Should be " << pNumNodes 
                     << " but got " << pGotNodes << " nodes." << endl;
                abort();
            }
        }

        /**
         * @brief Return the number of elements of the expansion dimension.
         */
        unsigned int Mesh::GetNumElements()
        {
            return element[expDim].size();
        }
        
        /**
         * @brief Return the number of boundary elements (i.e. one below the
         * expansion dimension).
         */
        unsigned int Mesh::GetNumBndryElements()
        {
            unsigned int i, nElmt = 0;
            
            for (i = 0; i < expDim; ++i)
                nElmt += element[i].size();
            
            return nElmt;
        }
        
        /**
         * @brief Return the total number of entities in the mesh (i.e. all
         * elements, regardless of dimension).
         */
        unsigned int Mesh::GetNumEntities()
        {
            unsigned int nEnt = 0;
            
            for (unsigned int d = 0; d <= expDim; ++d)
            {
                nEnt += element[d].size();
            }
            
            return nEnt;
        }
        
        /**
         * @brief Test equality of two conditions - i.e. compare types, fields
         * and values but _not_ composite ids.
         */
        bool operator==(ConditionSharedPtr const &c1, ConditionSharedPtr const &c2)
        {
            int i, n = c1->type.size();
            
            if (n != c2->type.size())
            {
                return false;
            }
            
            for (i = 0; i < n; ++i)
            {
                if (c1->type[i] != c2->type[i])
                {
                    return false;
                }
                
                if (c1->field[i] != c2->field[i] || c1->value[i] != c2->value[i])
                {
                    return false;
                }
            }
            
            return true;
        }

        /**
         * @brief Defines equality between two #NodeSharedPtr objects.
         */
        bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2)
        {
            return *p1 == *p2;
        }

        /**
         * @brief Defines ordering between two #NodeSharedPtr objects.
         */
        bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2)
        {
            return *p1 < *p2;
        }

        /**
         * @brief Defines equality of two edges (equal if IDs of end nodes
         * match in either ordering).
         */
        bool operator==(EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
        {
            return ( ((*(p1->n1) == *(p2->n1)) && (*(p1->n2) == *(p2->n2)))
                  || ((*(p1->n2) == *(p2->n1)) && (*(p1->n1) == *(p2->n2))));
        }

        /**
         * @brief Defines ordering between two edges (based on ID of edges).
         */
        bool operator< (EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
        {
            return p1->id < p2->id;
        }

        /**
         * @brief Defines equality of two faces (equal if IDs of vertices are
         * the same.)
         */
        bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2)
        {
            std::vector<NodeSharedPtr>::iterator it1, it2;
            for (it1 = p1->vertexList.begin(); it1 != p1->vertexList.end(); ++it1)
            {
                if (find(p2->vertexList.begin(), p2->vertexList.end(), *it1)
                    == p2->vertexList.end())
                {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Defines ordering between two faces (depending on ID of
         * faces).
         */
        bool operator< (FaceSharedPtr const &p1, FaceSharedPtr const &p2)
        {
            return p1->id < p2->id;
        }

        /**
         * @brief Replace a vertex in the element.
         * 
         * When a vertex is replaced, the element edges and faces are also
         * searched and the corresponding edge/face nodes are updated to
         * maintain consistency.
         * 
         * @param  p     Index of the vertex to replace.
         * @param  pNew  New vertex.
         */
        void Element::SetVertex(unsigned int p, NodeSharedPtr pNew)
        {
            NodeSharedPtr vOld = vertex[p];
            vertex[p] = pNew;
            for (unsigned int i = 0; i < edge.size(); ++i)
            {
                if (edge[i]->n1 == vOld)
                {
                    edge[i]->n1 = pNew;
                }
                else if (edge[i]->n2 == vOld)
                {
                    edge[i]->n2 = pNew;
                }
            }
            for (unsigned int i = 0; i < face.size(); ++i)
            {
                // Replace vertices in faces
                for (unsigned int j = 0; j < face[i]->vertexList.size(); ++j)
                {
                    if (face[i]->vertexList[j] == vOld)
                    {
                        face[i]->vertexList[j] = pNew;
                    }
                }
                for (unsigned int j = 0; j < face[i]->edgeList.size(); ++j)
                {
                    if (face[i]->edgeList[j]->n1 == vOld)
                    {
                        face[i]->edgeList[j]->n1 = pNew;
                    }
                    else if (face[i]->edgeList[j]->n2 == vOld)
                    {
                        face[i]->edgeList[j]->n2 = pNew;
                    }
                }
            }
        }

        /**
         * @brief Replace an edge in the element.
         * 
         * When an edge is replaced, the element faces are also searched and
         * the corresponding face edges are updated to maintain consistency.
         * 
         * @param  p     Index of the edge to replace.
         * @param  pNew  New edge.
         */
        void Element::SetEdge(unsigned int p, EdgeSharedPtr pNew)
        {
            EdgeSharedPtr vOld = edge[p];
            edge[p] = pNew;
            for (unsigned int i = 0; i < face.size(); ++i)
            {
                for (unsigned int j = 0; j < face[i]->edgeList.size(); ++j)
                {
                    if (face[i]->edgeList[j] == vOld)
                    {
                        face[i]->edgeList[j] = pNew;
                    }
                }
            }
        }
        
        /**
         * @brief Replace a face in the element.
         * 
         * When a face is replaced, no other consistency checks are required.
         * 
         * @param  p     Index of the face to replace.
         * @param  pNew  New face.
         */
        void Element::SetFace(unsigned int p, FaceSharedPtr pNew)
        {
            face[p] = pNew;
        }
        
        /**
         * @brief Obtain the order of an element by looking at edges.
         */
        int Element::GetMaxOrder()
        {
            int i, ret = 1;
            
            for (i = 0; i < edge.size(); ++i)
            {
                int edgeOrder = edge[i]->GetNodeCount()-1;
                if (edgeOrder > ret)
                {
                    ret = edgeOrder;
                }
            }
            
            return ret;
        }
        
        /**
         * @brief Generate a Nektar++ string describing the composite.
         * 
         * The list of composites may include individual element IDs or ranges
         * of element IDs.
         */
        string Composite::GetXmlString(bool doSort)
        {
            if (doSort)
            {
                element_id_less_than sortOperator;
                sort(items.begin(), items.end(), sortOperator);
            }

            stringstream st;
            vector<ElementSharedPtr>::iterator it;
            bool range = false;
            int vId = items[0]->GetId();
            int prevId = vId;

            st << " " << tag << "[" << vId;

            for (it = items.begin()+1; it != items.end(); ++it){
                // store previous element ID and get current one
                prevId = vId;
                vId = (*it)->GetId();

                // continue an already started range
                if (prevId > -1 && vId == prevId + 1)
                {
                    range = true;
                    // if this is the last element, it's the end of a range, so write
                    if (*it == items.back())
                    {
                        st << "-" << vId;
                    }
                    continue;
                }

                // terminate a range, if present
                if (range)
                {
                    st << "-" << prevId;
                    range = false;
                }

                // write what will be either a single entry or start of new range
                st << "," << vId;
            }
            // terminate
            st << "] ";
            return st.str();
        }

        ElementType Point::type = GetElementFactory().
            RegisterCreatorFunction(ePoint, Point::create, "Point");
        
        /**
         * @brief Create a point element.
         */
        Point::Point(ElmtConfig            pConf,
                     vector<NodeSharedPtr> pNodeList, 
                     vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag     = "";
            m_dim     = 0;
            m_taglist = pTagList;
            vertex.push_back(pNodeList[0]);
        }

        /**
         * @brief Return the number of nodes defining a point (i.e. return 1).
         */
        unsigned int Point::GetNumNodes(ElmtConfig pConf)
        {
            return 1;
        }


        ElementType Line::type = GetElementFactory().
            RegisterCreatorFunction(eLine, Line::create, "Line");
        
        /**
         * @brief Create a line element.
         */
        Line::Line(ElmtConfig            pConf,
                   vector<NodeSharedPtr> pNodeList, 
                   vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag     = "S";
            m_dim     = 1;
            m_taglist = pTagList;
            int n     = m_conf.order-1;
            
            // Add vertices
            for (int i = 0; i < 2; ++i) {
                vertex.push_back(pNodeList[i]);
            }
            vector<NodeSharedPtr> edgeNodes;
            if (m_conf.order > 1) {
                for (int j = 0; j<n; ++j) {
                    edgeNodes.push_back(pNodeList[2+j]);
                }
            }
            edge.push_back(boost::shared_ptr<Edge>(
                new Edge(pNodeList[0], pNodeList[1], edgeNodes, m_conf.edgeCurveType)));
        }
        
        SpatialDomains::GeometrySharedPtr Line::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }
                
            // Create edge vertices.
            SpatialDomains::VertexComponentSharedPtr p[2];
            p[0] = vertex[0]->GetGeom(coordDim);
            p[1] = vertex[1]->GetGeom(coordDim);
            
            if (edge[0]->edgeNodes.size() > 0)
            {
                SpatialDomains::CurveSharedPtr c = 
                    MemoryManager<SpatialDomains::Curve>::
                    AllocateSharedPtr(m_id, edge[0]->curveType);
                
                c->m_points.push_back(p[0]);
                for (int i = 0; i < edge[0]->edgeNodes.size(); ++i)
                {
                    c->m_points.push_back(edge[0]->edgeNodes[i]->GetGeom(coordDim));
                }
                c->m_points.push_back(p[1]);
                
                m_geom = MemoryManager<SpatialDomains::SegGeom>::
                    AllocateSharedPtr(m_id, 2, p, c);
            }
            else
            {
                m_geom = MemoryManager<SpatialDomains::SegGeom>::
                    AllocateSharedPtr(m_id, 2, p);
            }
            
            return m_geom;
        }

        /**
         * @brief Return the number of nodes defining a line.
         */
        unsigned int Line::GetNumNodes(ElmtConfig pConf)
        {
            return pConf.order+1;
        }


        ElementType Triangle::type = GetElementFactory().
            RegisterCreatorFunction(eTriangle, Triangle::create, "Triangle");
        
        /**
         * @brief Create a triangle element.
         */
        Triangle::Triangle(ElmtConfig            pConf,
                           vector<NodeSharedPtr> pNodeList, 
                           vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag     = "T";
            m_dim     = 2;
            m_taglist = pTagList;
            int n     = m_conf.order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge. This is based on the ordering produced by
            // gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 4;
            edgeNodeMap[pair<int,int>(2,3)] = 4 + n;
            edgeNodeMap[pair<int,int>(3,1)] = 4 + 2*n;

            // Add vertices
            for (int i = 0; i < 3; ++i) {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                      pNodeList[it->first.second-1],
                                                      edgeNodes,
                                                      m_conf.edgeCurveType)));
            }

            if (m_conf.faceNodes)
            {
                volumeNodes.insert(volumeNodes.begin(),
                                   pNodeList.begin()+3*m_conf.order,
                                   pNodeList.end());
            }

        }

        SpatialDomains::GeometrySharedPtr Triangle::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }

            SpatialDomains::SegGeomSharedPtr         edges[3];
            SpatialDomains::VertexComponentSharedPtr verts[3];
            
            for (int i = 0; i < 3; ++i)
            {
                edges[i] = edge  [i]->GetGeom(coordDim);
                verts[i] = vertex[i]->GetGeom(coordDim);
            }
            
            StdRegions::Orientation edgeorient[3] = {
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
            };
            
            m_geom = MemoryManager<SpatialDomains::TriGeom>::
                AllocateSharedPtr(m_id, verts, edges, edgeorient);
            
            return m_geom;
        }
        
        /**
         * @brief Return the number of nodes defining a triangle.
         */
        unsigned int Triangle::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (!pConf.faceNodes)
                return (n+1)+2*(n-1)+1;
            else
                return (n+1)*(n+2)/2;
        }

        void Triangle::Complete(int order)
        {
            int i, j;
            
            // Create basis key for a nodal tetrahedron.
            LibUtilities::BasisKey B0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B1(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            
            // Create a standard nodal triangle in order to get the
            // Vandermonde matrix to perform interpolation to nodal points.
            StdRegions::StdNodalTriExpSharedPtr nodalTri = 
                MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
                    B0, B1, LibUtilities::eNodalTriEvenlySpaced);
            
            SpatialDomains::TriGeomSharedPtr geom = 
                boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                    this->GetGeom(3));
            
            // Create basis key for a triangle.
            LibUtilities::BasisKey C0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey C1(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            
            // Create a triangle.
            LocalRegions::TriExpSharedPtr tri = 
                MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(
                    C0, C1, geom);
            
            // Get coordinate array for tetrahedron.
            int nqtot = tri->GetTotPoints();
            Array<OneD, NekDouble> alloc(6*nqtot);
            Array<OneD, NekDouble> xi(alloc);
            Array<OneD, NekDouble> yi(alloc+  nqtot);
            Array<OneD, NekDouble> zi(alloc+2*nqtot);
            Array<OneD, NekDouble> xo(alloc+3*nqtot);
            Array<OneD, NekDouble> yo(alloc+4*nqtot);
            Array<OneD, NekDouble> zo(alloc+5*nqtot);
            Array<OneD, NekDouble> tmp;
            
            tri->GetCoords(xi, yi, zi);
            
            for (i = 0; i < 3; ++i)
            {
                tri->FwdTrans(alloc+i*nqtot, nodalTri->UpdateCoeffs());
                // Apply Vandermonde matrix to project onto nodal space.
                nodalTri->ModalToNodal(nodalTri->GetCoeffs(), tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 3; ++i)
            {
                int pos = 3 + i*(order-1);
                edge[i]->edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    edge[i]->edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Extract face-interior nodes.
            int pos = 3 + 3*(order-1);
            for (i = pos; i < (order+1)*(order+2)/2; ++i)
            {
                volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.order       = order;
            m_conf.faceNodes   = true;
            m_conf.volumeNodes = true;
        }


        ElementType Quadrilateral::type = GetElementFactory().
            RegisterCreatorFunction(eQuadrilateral, Quadrilateral::create, 
                                    "Quadrilateral");
        
        /**
         * @brief Create a quadrilateral element.
         */
        Quadrilateral::Quadrilateral(ElmtConfig            pConf,
                                     vector<NodeSharedPtr> pNodeList,
                                     vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag = "Q";
            m_dim = 2;
            m_taglist = pTagList;
            int n = m_conf.order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge. This is based on the ordering produced by
            // gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 5;
            edgeNodeMap[pair<int,int>(2,3)] = 5 + n;
            edgeNodeMap[pair<int,int>(3,4)] = 5 + 2*n;
            edgeNodeMap[pair<int,int>(4,1)] = 5 + 3*n;

            // Add vertices
            for (int i = 0; i < 4; ++i) 
            {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                      pNodeList[it->first.second-1],
                                                      edgeNodes,
                                                      m_conf.edgeCurveType)));
            }

            if (m_conf.faceNodes)
            {
                volumeNodes.insert(volumeNodes.begin(), 
                                   pNodeList.begin()+4*m_conf.order,
                                   pNodeList.end());
            }
        }

        void Quadrilateral::Complete(int order)
        {
            LibUtilities::BasisKey C0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            
            SpatialDomains::QuadGeomSharedPtr geom = 
                boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(
                    this->GetGeom(3));
            
            // Create a quad.
            LocalRegions::QuadExpSharedPtr quad = 
                MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(
                    C0, C0, geom);
            
            // Get coordinate array for quadrilateral.
            int nqtot = quad->GetTotPoints();
            Array<OneD, NekDouble> alloc(3*nqtot);
            Array<OneD, NekDouble> x    (alloc        );
            Array<OneD, NekDouble> y    (alloc+1*nqtot);
            Array<OneD, NekDouble> z    (alloc+2*nqtot);
            
            quad->GetCoords(x, y, z);
            
            // Now extract points from the co-ordinate arrays into the edge
            // and face nodes. First, extract edge-interior nodes.
            int edgeMap[4][2] = {{0,1},{order,order+1},
                                 {nqtot-1,-1},{order*(order+1),-order-1}};
            
            for (int i = 0; i < 4; ++i)
            {
                int pos = edgeMap[i][0] + edgeMap[i][1];
                edge[i]->edgeNodes.clear();
                /*
                cout << "EDGE " << i << " = " 
                     << edge[i]->n1->x << "," << edge[i]->n1->y << " "
                     << edge[i]->n2->x << "," << edge[i]->n2->y << endl;
                */
                for (int j = 1; j < order; ++j, pos += edgeMap[i][1])
                {
                    //cout << "INSERTING: " << x[pos] << " " << y[pos] << endl;
                    edge[i]->edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, x[pos], y[pos], z[pos])));
                }
            }

            // Extract face-interior nodes.
            volumeNodes.clear();
            for (int i = 1; i < order; ++i)
            {
                int pos = i*(order+1);
                for (int j = 1; j < order; ++j)
                {
                    volumeNodes.push_back(
                        NodeSharedPtr(new Node(0, x[pos+j], y[pos+j], z[pos+j])));
                    //cout << "here" << endl;
                }                
            }
            
            m_conf.order       = order;
            m_conf.faceNodes   = true;
            m_conf.volumeNodes = true;
        }

        SpatialDomains::GeometrySharedPtr Quadrilateral::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }

            SpatialDomains::SegGeomSharedPtr         edges[4];
            SpatialDomains::VertexComponentSharedPtr verts[4];
            
            for (int i = 0; i < 4; ++i)
            {
                edges[i] = edge  [i]->GetGeom(coordDim);
                verts[i] = vertex[i]->GetGeom(coordDim);
            }
            
            StdRegions::Orientation edgeorient[4] = {
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
            };
            
            m_geom = MemoryManager<SpatialDomains::QuadGeom>::
                AllocateSharedPtr(m_id, verts, edges, edgeorient);
            
            return m_geom;
        }

        /**
         * @brief Return the number of nodes defining a quadrilateral.
         */
        unsigned int Quadrilateral::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (!pConf.faceNodes)
                return 4*n;
            else
                return (n+1)*(n+1);
        }
        
        
        ElementType Tetrahedron::type = GetElementFactory().
            RegisterCreatorFunction(eTetrahedron, Tetrahedron::create, "Tetrahedron");

        /**
         * @brief Create a tetrahedron element.
         */
        Tetrahedron::Tetrahedron(ElmtConfig            pConf,
                                 vector<NodeSharedPtr> pNodeList,
                                 vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size())
        {
            m_tag = "A";
            m_dim = 3;
            m_taglist = pTagList;
            int n = m_conf.order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 5;
            edgeNodeMap[pair<int,int>(2,3)] = 5 + n;
            edgeNodeMap[pair<int,int>(1,3)] = 5 + 2*n;
            edgeNodeMap[pair<int,int>(1,4)] = 5 + 3*n;
            edgeNodeMap[pair<int,int>(2,4)] = 5 + 4*n;
            edgeNodeMap[pair<int,int>(3,4)] = 5 + 5*n;
            
            // Add vertices
            for (int i = 0; i < 4; ++i) {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            int eid = 0;
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                      pNodeList[it->first.second-1],
                                                      edgeNodes,
                                                      m_conf.edgeCurveType)));
                edge.back()->id = eid++;
            }
            
            // Reorient the tet to ensure collapsed coordinates align between adjacent
            // elements.
            if (m_conf.reorient)
            {
                OrientTet();
            }
            
            // Create faces
            int face_ids[4][3] = {
                {0,1,2},{0,1,3},{1,2,3},{0,2,3}};
            int face_edges[4][3];
            
            for (int j = 0; j < 4; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                for (int k = 0; k < 3; ++k)
                {
                    faceVertices.push_back(vertex[face_ids[j][k]]);
                    NodeSharedPtr a = vertex[face_ids[j][k]];
                    NodeSharedPtr b = vertex[face_ids[j][(k+1)%3]];
                    for (unsigned int i = 0; i < edge.size(); ++i)
                    {
                        if ( ((*(edge[i]->n1)==*a) && (*(edge[i]->n2)==*b))
                                || ((*(edge[i]->n1)==*b) && (*(edge[i]->n2) == *a)) )
                        {
                            face_edges[j][k] = i;
                            faceEdges.push_back(edge[i]);
                            break;
                        }
                    }
                }

                if (m_conf.faceNodes)
                {
                    int N = 4 + 6*n + j*n*(n-1)/2;
                    for (int i = 0; i < n*(n-1)/2; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }

            vector<EdgeSharedPtr> tmp(6);
            tmp[0] = edge[face_edges[0][0]];
            tmp[1] = edge[face_edges[0][1]];
            tmp[2] = edge[face_edges[0][2]];
            tmp[3] = edge[face_edges[1][2]];
            tmp[4] = edge[face_edges[1][1]];
            tmp[5] = edge[face_edges[2][1]];
            edge = tmp;
        }
        
        SpatialDomains::GeometrySharedPtr Tetrahedron::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }
            
            SpatialDomains::TriGeomSharedPtr tfaces[4];
            
            for (int i = 0; i < 4; ++i)
            {
                tfaces[i] = boost::dynamic_pointer_cast
                    <SpatialDomains::TriGeom>(face[i]->GetGeom(coordDim));
            }

            m_geom = MemoryManager<SpatialDomains::TetGeom>::
                AllocateSharedPtr(tfaces);
            
            return m_geom;
        }
        
        /**
         * @brief Return the number of nodes defining a tetrahedron.
         */
        unsigned int Tetrahedron::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (pConf.volumeNodes && pConf.faceNodes)
                return (n+1)*(n+2)*(n+3)/6;
            else if (!pConf.volumeNodes && pConf.faceNodes)
                return 4*(n+1)*(n+2)/2-6*(n+1)+4;
            else
                return 6*(n+1)-8;
        }

        /**
         * @brief .
         */
        void Tetrahedron::Complete(int order)
        {
            int i, j;
            
            // Create basis key for a nodal tetrahedron.
            LibUtilities::BasisKey B0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B1(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            LibUtilities::BasisKey B2(
                LibUtilities::eOrtho_C, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha2Beta0));
            
            // Create a standard nodal tetrahedron in order to get the
            // Vandermonde matrix to perform interpolation to nodal points.
            StdRegions::StdNodalTetExpSharedPtr nodalTet = 
                MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
                    B0, B1, B2, LibUtilities::eNodalTetEvenlySpaced);
            
            Array<OneD, NekDouble> x, y, z;
            
            nodalTet->GetNodalPoints(x,y,z);
            
            SpatialDomains::TetGeomSharedPtr geom = 
                boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(
                    this->GetGeom(3));
            
            // Create basis key for a tetrahedron.
            LibUtilities::BasisKey C0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey C1(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            LibUtilities::BasisKey C2(
                LibUtilities::eOrtho_C, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha2Beta0));
            
            // Create a tet.
            LocalRegions::TetExpSharedPtr tet = 
                MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(
                    C0, C1, C2, geom);
            
            // Get coordinate array for tetrahedron.
            int nqtot = tet->GetTotPoints();
            Array<OneD, NekDouble> alloc(6*nqtot);
            Array<OneD, NekDouble> xi(alloc);
            Array<OneD, NekDouble> yi(alloc+  nqtot);
            Array<OneD, NekDouble> zi(alloc+2*nqtot);
            Array<OneD, NekDouble> xo(alloc+3*nqtot);
            Array<OneD, NekDouble> yo(alloc+4*nqtot);
            Array<OneD, NekDouble> zo(alloc+5*nqtot);
            Array<OneD, NekDouble> tmp;
            
            tet->GetCoords(xi, yi, zi);
            
            for (i = 0; i < 3; ++i)
            {
                tet->FwdTrans(alloc+i*nqtot, nodalTet->UpdateCoeffs());
                // Apply Vandermonde matrix to project onto nodal space.
                nodalTet->ModalToNodal(nodalTet->GetCoeffs(), tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 6; ++i)
            {
                int pos = 4 + i*(order-1);
                edge[i]->edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    edge[i]->edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Now extract face-interior nodes.
            for (i = 0; i < 4; ++i)
            {
                int pos = 4 + 6*(order-1) + i*(order-2)*(order-1)/2;
                face[i]->faceNodes.clear();
                for (j = 0; j < (order-2)*(order-1)/2; ++j)
                {
                    face[i]->faceNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }
            
            // Finally extract volume nodes.
            int pos = 4 + 6*(order-1) + 4*(order-2)*(order-1)/2;
            for (i = pos; i < (order+1)*(order+2)*(order+3)/6; ++i)
            {
                volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.order       = order;
            m_conf.faceNodes   = true;
            m_conf.volumeNodes = true;
        }

        struct TetOrient
        {
            TetOrient(vector<int> nid, int fid) : nid(nid), fid(fid) {}
            vector<int> nid;
            int fid;
        };
        
        struct TetOrientHash : std::unary_function<struct TetOrient, std::size_t>
        {
            std::size_t operator()(struct TetOrient const& p) const
            {
                return boost::hash_range(p.nid.begin(), p.nid.end());
            }
        };
        typedef boost::unordered_set<struct TetOrient, TetOrientHash> TetOrientSet;

        bool operator==(const struct TetOrient &a, const struct TetOrient &b)
        {
            if (a.nid.size() != b.nid.size())
            {
                return false;
            }
            
            for (int i = 0; i < a.nid.size(); ++i)
            {
                if (a.nid[i] != b.nid[i])
                {
                    return false;
                }
            }
            
            return true;
        }
        
        /**
         * @brief Orient tetrahedron to align degenerate vertices.
         * 
         * Orientation of tetrahedral elements is required so that the
         * singular vertices of triangular faces (which occur as a part of the
         * collapsed co-ordinate system) align. The algorithm is based on that
         * used in T. Warburton's thesis and in the original Nektar source.
         * 
         * First the vertices are ordered with the highest global vertex at
         * the top degenerate point, and the base degenerate point has second
         * lowest ID. These vertices are swapped if the element is incorrectly
         * oriented.
         */
        void Tetrahedron::OrientTet()
        {
            static int face_ids[4][3] = {
                {0,1,2},{0,1,3},{1,2,3},{0,2,3}};
            
            TetOrientSet faces;

            // Create a copy of the original vertex ordering. This is used to
            // construct a mapping, #orientationMap, which maps the original
            // face ordering to the new face ordering.
            for (int i = 0; i < 4; ++i)
            {
                vector<int> nodes(3);
                
                nodes[0] = vertex[face_ids[i][0]]->id;
                nodes[1] = vertex[face_ids[i][1]]->id;
                nodes[2] = vertex[face_ids[i][2]]->id;
                
                sort(nodes.begin(), nodes.end());
                struct TetOrient faceNodes(nodes, i);
                faces.insert(faceNodes);
            }
            
            // Order vertices with highest global vertex at top degenerate
            // point. Place second highest global vertex at base degenerate
            // point.
            sort(vertex.begin(), vertex.end());
            
            // Calculate a.(b x c) to determine tet volume; if negative,
            // reverse order of non-degenerate points to correctly orientate
            // the tet.
            double ax  = vertex[1]->x-vertex[0]->x;
            double ay  = vertex[1]->y-vertex[0]->y;
            double az  = vertex[1]->z-vertex[0]->z;
            double bx  = vertex[2]->x-vertex[0]->x;
            double by  = vertex[2]->y-vertex[0]->y;
            double bz  = vertex[2]->z-vertex[0]->z;
            double cx  = vertex[3]->x-vertex[0]->x;
            double cy  = vertex[3]->y-vertex[0]->y;
            double cz  = vertex[3]->z-vertex[0]->z;
            double vol = cx*(ay*bz-az*by)+cy*(az*bx-ax*bz)+cz*(ax*by-ay*bx);
            vol       /= 6.0;
            
            if (fabs(vol) <= 1e-10)
            {
                cerr << "Warning: degenerate tetrahedron, volume = " << vol << endl;
            }

            if (vol < 0)
            {
                swap(vertex[0], vertex[1]);
            }
            
            TetOrientSet::iterator it;
            
            // Search for the face in the original set of face nodes. Then use
            // this to construct the #orientationMap.
            for (int i = 0; i < 4; ++i)
            {
                vector<int> nodes(3);
                
                nodes[0] = vertex[face_ids[i][0]]->id;
                nodes[1] = vertex[face_ids[i][1]]->id;
                nodes[2] = vertex[face_ids[i][2]]->id;
                sort(nodes.begin(), nodes.end());
                
                struct TetOrient faceNodes(nodes, 0);
                
                it = faces.find(faceNodes);
                orientationMap[it->fid] = i;
            }
        }


        ElementType Prism::type = GetElementFactory().
            RegisterCreatorFunction(ePrism, Prism::create, "Prism");
        
        /**
         * @brief Create a prism element.
         */
        Prism::Prism(ElmtConfig            pConf,
                     vector<NodeSharedPtr> pNodeList,
                     vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size())
        {
            m_tag     = "R";
            m_dim     = 3;
            m_taglist = pTagList;
            int n     = m_conf.order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge. This is based on the ordering produced by
            // gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;

            // This edge-node map is based on Nektar++ ordering.
            edgeNodeMap[pair<int,int>(1,2)] = 7;
            edgeNodeMap[pair<int,int>(2,3)] = 7 + n;
            edgeNodeMap[pair<int,int>(4,3)] = 7 + 2*n;
            edgeNodeMap[pair<int,int>(1,4)] = 7 + 3*n;
            edgeNodeMap[pair<int,int>(1,5)] = 7 + 4*n;
            edgeNodeMap[pair<int,int>(2,5)] = 7 + 5*n;
            edgeNodeMap[pair<int,int>(3,6)] = 7 + 6*n;
            edgeNodeMap[pair<int,int>(4,6)] = 7 + 7*n;
            edgeNodeMap[pair<int,int>(5,6)] = 7 + 8*n;

            // Add vertices
            for (int i = 0; i < 6; ++i)
            {
                vertex.push_back(pNodeList[i]);
            }

            int eid = 0;
            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(
                    new Edge(pNodeList[it->first.first-1],
                             pNodeList[it->first.second-1],
                             edgeNodes,
                             m_conf.edgeCurveType)));
                edge.back()->id = eid++;
            }
            
            if (m_conf.reorient)
            {
                OrientPrism();
            }
            
            // Create faces
            int face_ids[5][4] = {
                {0,1,2,3},{0,1,4,-1},{1,2,5,4},{3,2,5,-1},{0,3,5,4}};
            int face_edges[5][4];
            int faceoffset = 0;
            for (int j = 0; j < 5; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                int nEdge = 3 - (j % 2 - 1);
                
                for (int k = 0; k < nEdge; ++k)
                {
                    faceVertices.push_back(vertex[face_ids[j][k]]);
                    NodeSharedPtr a = vertex[face_ids[j][k]];
                    NodeSharedPtr b = vertex[face_ids[j][(k+1) % nEdge]];
                    for (unsigned int i = 0; i < edge.size(); ++i)
                    {
                        if ((edge[i]->n1 == a && edge[i]->n2 == b) || 
                            (edge[i]->n1 == b && edge[i]->n2 == a))
                        {
                            faceEdges.push_back(edge[i]);
                            face_edges[j][k] = i;
                            break;
                        }
                    }
                }
                
                if (m_conf.faceNodes)
                {
                    int facenodes = j%2==0 ? n*n : n*(n-1)/2;
                    faceoffset   += facenodes;
                    int N = 6 + 9*n + faceoffset;
                    for (int i = 0; i < facenodes; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }
            
            // Re-order edge array to be consistent with Nektar++ ordering.
            vector<EdgeSharedPtr> tmp(9);
            tmp[0] = edge[face_edges[0][0]];
            tmp[1] = edge[face_edges[0][1]];
            tmp[2] = edge[face_edges[0][2]];
            tmp[3] = edge[face_edges[0][3]];
            tmp[4] = edge[face_edges[1][2]];
            tmp[5] = edge[face_edges[1][1]];
            tmp[6] = edge[face_edges[2][1]];
            tmp[7] = edge[face_edges[3][2]];
            tmp[8] = edge[face_edges[4][2]];
            edge = tmp;
        }

        /**
         * @brief Return the number of nodes defining a prism.
         */
        unsigned int Prism::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (pConf.faceNodes && pConf.volumeNodes)
                return (n+1)*(n+1)*(n+2)/2;
            else if (pConf.faceNodes && !pConf.volumeNodes)
                return 3*(n+1)*(n+1)+2*(n+1)*(n+2)/2-9*(n+1)+6;
            else
                return 9*(n+1)-12;
        }

        SpatialDomains::GeometrySharedPtr Prism::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }

            SpatialDomains::Geometry2DSharedPtr faces[5];
            
            for (int i = 0; i < 5; ++i)
            {
                faces[i] = face[i]->GetGeom(coordDim);
            }

            m_geom = MemoryManager<SpatialDomains::PrismGeom>::
                AllocateSharedPtr(faces);
            
            return m_geom;
        }

        /**
         * @brief .
         */
        void Prism::Complete(int order)
        {
            int i, j, pos;
            
            // Create basis key for a nodal tetrahedron.
            LibUtilities::BasisKey B0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B1(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B2(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            
            // Create a standard nodal prism in order to get the Vandermonde
            // matrix to perform interpolation to nodal points.
            StdRegions::StdNodalPrismExpSharedPtr nodalPrism = 
                MemoryManager<StdRegions::StdNodalPrismExp>::AllocateSharedPtr(
                    B0, B1, B2, LibUtilities::eNodalPrismEvenlySpaced);
            
            Array<OneD, NekDouble> x, y, z;
            nodalPrism->GetNodalPoints(x,y,z);
            
            SpatialDomains::PrismGeomSharedPtr geom = 
                boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(
                    this->GetGeom(3));
            
            // Create basis key for a prism.
            LibUtilities::BasisKey C0(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey C1(
                LibUtilities::eOrtho_A, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey C2(
                LibUtilities::eOrtho_B, order+1,
                LibUtilities::PointsKey(
                    order+1,LibUtilities::eGaussRadauMAlpha1Beta0));
            
            // Create a prism.
            LocalRegions::PrismExpSharedPtr prism = 
                MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(
                    C0, C1, C2, geom);
            
            // Get coordinate array for tetrahedron.
            int nqtot = prism->GetTotPoints();
            Array<OneD, NekDouble> alloc(6*nqtot);
            Array<OneD, NekDouble> xi(alloc);
            Array<OneD, NekDouble> yi(alloc+  nqtot);
            Array<OneD, NekDouble> zi(alloc+2*nqtot);
            Array<OneD, NekDouble> xo(alloc+3*nqtot);
            Array<OneD, NekDouble> yo(alloc+4*nqtot);
            Array<OneD, NekDouble> zo(alloc+5*nqtot);
            Array<OneD, NekDouble> tmp;
            
            prism->GetCoords(xi, yi, zi);
            
            for (i = 0; i < 3; ++i)
            {
                prism->FwdTrans(alloc+i*nqtot, nodalPrism->UpdateCoeffs());
                // Apply Vandermonde matrix to project onto nodal space.
                nodalPrism->ModalToNodal(nodalPrism->GetCoeffs(), tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 9; ++i)
            {
                pos = 6 + i*(order-1);
                edge[i]->edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    edge[i]->edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Now extract face-interior nodes.
            pos = 6 + 9*(order-1);
            for (i = 0; i < 5; ++i)
            {
                int facesize = i % 2 ? (order-2)*(order-1)/2 : (order-1)*(order-1);
                face[i]->faceNodes.clear();
                for (j = 0; j < facesize; ++j)
                {
                    face[i]->faceNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
                pos += facesize;
            }
            
            // Finally extract volume nodes.
            for (i = pos; i < (order+1)*(order+1)*(order+2)/2; ++i)
            {
                volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.order       = order;
            m_conf.faceNodes   = true;
            m_conf.volumeNodes = true;
        }

        /**
         * @brief Orient prism to align degenerate vertices.
         * 
         * Orientation of prismatric elements is required so that the singular
         * vertices of triangular faces (which occur as a part of the
         * collapsed co-ordinate system) align. The algorithm is based on that
         * used in T. Warburton's thesis and in the original Nektar source.
         * 
         * First the points are re-ordered so that the highest global IDs
         * represent the two singular points of the prism. Then, if necessary,
         * the nodes are rotated either clockwise or counter-clockwise (w.r.t
         * to the p-r plane) to correctly align the prism. The #orientation
         * variable is set to:
         * 
         * - 0 if the prism is not rotated;
         * - 1 if the prism is rotated clockwise;
         * - 2 if the prism is rotated counter-clockwise.
         * 
         * This is necessary for some input modules (e.g. #InputNek) which add
         * high-order information or bounary conditions to faces.
         */
        void Prism::OrientPrism()
        {
            int lid[6], gid[6];
            
            // Re-order vertices.
            for (int i = 0; i < 6; ++i)
            {
                lid[i] = i;
                gid[i] = vertex[i]->id;
            }

            gid[0] = gid[3] = max(gid[0], gid[3]);
            gid[1] = gid[2] = max(gid[1], gid[2]);
            gid[4] = gid[5] = max(gid[4], gid[5]);
            
            for (int i = 1; i < 6; ++i)
            {
                if (gid[0] < gid[i])
                {
                    swap(gid[i], gid[0]);
                    swap(lid[i], lid[0]);
                }
            }

            if (lid[0] == 4 || lid[0] == 5) 
            {
                orientation = 0;
            } 
            else if (lid[0] == 1 || lid[0] == 2) 
            {
                // Rotate prism clockwise in p-r plane
                vector<NodeSharedPtr> vertexmap(6);
                vertexmap[0] = vertex[4];
                vertexmap[1] = vertex[0];
                vertexmap[2] = vertex[3];
                vertexmap[3] = vertex[5];
                vertexmap[4] = vertex[1];
                vertexmap[5] = vertex[2];
                vertex = vertexmap;
                orientation = 1;
            }
            else if (lid[0] == 0 || lid[0] == 3)
            {
                // Rotate prism counter-clockwise in p-r plane
                vector<NodeSharedPtr> vertexmap(6);
                vertexmap[0] = vertex[1];
                vertexmap[1] = vertex[4];
                vertexmap[2] = vertex[5];
                vertexmap[3] = vertex[2];
                vertexmap[4] = vertex[0];
                vertexmap[5] = vertex[3];
                vertex = vertexmap;
                orientation = 2;
            }
            else
            {
                cerr << "Warning: possible prism orientation problem." << endl;
            }
        }


        ElementType Hexahedron::type = GetElementFactory().
            RegisterCreatorFunction(eHexahedron, Hexahedron::create, "Hexahedron");

        /**
         * @brief Create a hexahedral element.
         */
        Hexahedron::Hexahedron(ElmtConfig            pConf,
                               vector<NodeSharedPtr> pNodeList,
                               vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag = "H";
            m_dim = 3;
            m_taglist = pTagList;
            int n = m_conf.order-1;
            
            // Create a map to relate edge nodes to a pair of vertices defining an edge
            // This is based on the ordering produced by gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 9;
            edgeNodeMap[pair<int,int>(1,4)] = 9 + n;
            edgeNodeMap[pair<int,int>(1,5)] = 9 + 2*n;
            edgeNodeMap[pair<int,int>(2,3)] = 9 + 3*n;
            edgeNodeMap[pair<int,int>(2,6)] = 9 + 4*n;
            edgeNodeMap[pair<int,int>(3,4)] = 9 + 5*n;
            edgeNodeMap[pair<int,int>(3,7)] = 9 + 6*n;
            edgeNodeMap[pair<int,int>(4,8)] = 9 + 7*n;
            edgeNodeMap[pair<int,int>(5,6)] = 9 + 8*n;
            edgeNodeMap[pair<int,int>(5,8)] = 9 + 9*n;
            edgeNodeMap[pair<int,int>(6,7)] = 9 + 10*n;
            edgeNodeMap[pair<int,int>(7,8)] = 9 + 11*n;

            // Add vertices
            for (int i = 0; i < 8; ++i) {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                      pNodeList[it->first.second-1],
                                                      edgeNodes,
                                                      m_conf.edgeCurveType)));
            }

            // Create faces
            int face_ids[6][4] = {
                {0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};
            for (int j = 0; j < 6; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                for (int k = 0; k < 4; ++k)
                {
                    faceVertices.push_back(vertex[face_ids[j][k]]);
                    NodeSharedPtr a = vertex[face_ids[j][k]];
                    NodeSharedPtr b = vertex[face_ids[j][(k+1)%4]];
                    for (unsigned int i = 0; i < edge.size(); ++i)
                    {
                        if ( ((*(edge[i]->n1)==*a) && (*(edge[i]->n2)==*b))
                                || ((*(edge[i]->n1)==*b) && (*(edge[i]->n2) == *a)) )
                        {
                            faceEdges.push_back(edge[i]);
                            break;
                        }
                    }
                }

                if (m_conf.faceNodes)
                {
                    int N = 8 + 12*n + j*n*n;
                    for (int i = 0; i < n*n; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }
        }
        
        SpatialDomains::GeometrySharedPtr Hexahedron::GetGeom(int coordDim)
        {
            if (m_geom)
            {
                return m_geom;
            }

            SpatialDomains::QuadGeomSharedPtr faces[6];
            
            for (int i = 0; i < 6; ++i)
            {
                faces[i] = boost::dynamic_pointer_cast
                    <SpatialDomains::QuadGeom>(face[i]->GetGeom(coordDim));
            }

            m_geom = MemoryManager<SpatialDomains::HexGeom>::
                AllocateSharedPtr(faces);
            
            return m_geom;
        }

        /**
         * @brief Return the number of nodes defining a hexahedron.
         */
        unsigned int Hexahedron::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (pConf.faceNodes && pConf.volumeNodes)
                return (n+1)*(n+1)*(n+1);
            else if (pConf.faceNodes && !pConf.volumeNodes)
                return 6*(n+1)*(n+1)-12*(n+1)+8;
            else
                return 12*(n+1)-16;
        }
    }
}
