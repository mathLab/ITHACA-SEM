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
            : m_conf(pConf), 
              m_curveType(LibUtilities::ePolyEvenlySpaced),
              m_geom()
        {
            if (pNumNodes != pGotNodes)
            {
                cerr << "Number of modes mismatch for type " 
                     << pConf.m_e << "! Should be " << pNumNodes 
                     << " but got " << pGotNodes << " nodes." << endl;
                abort();
            }
        }

        /**
         * @brief Return the number of elements of the expansion dimension.
         */
        unsigned int Mesh::GetNumElements()
        {
            return m_element[m_expDim].size();
        }
        
        /**
         * @brief Return the number of boundary elements (i.e. one below the
         * expansion dimension).
         */
        unsigned int Mesh::GetNumBndryElements()
        {
            unsigned int i, nElmt = 0;
            
            for (i = 0; i < m_expDim; ++i)
                nElmt += m_element[i].size();
            
            return nElmt;
        }
        
        /**
         * @brief Return the total number of entities in the mesh (i.e. all
         * elements, regardless of dimension).
         */
        unsigned int Mesh::GetNumEntities()
        {
            unsigned int nEnt = 0;
            
            for (unsigned int d = 0; d <= m_expDim; ++d)
            {
                nEnt += m_element[d].size();
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

        std::ostream &operator<<(std::ostream &os, const NodeSharedPtr &n)
        {
            os << n->m_x << " " << n->m_y << " " << n->m_z;
            return os;
        }

        /**
         * @brief Defines equality of two edges (equal if IDs of end nodes
         * match in either ordering).
         */
        bool operator==(EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
        {
            return ( ((*(p1->m_n1) == *(p2->m_n1)) && (*(p1->m_n2) == *(p2->m_n2)))
                  || ((*(p1->m_n2) == *(p2->m_n1)) && (*(p1->m_n1) == *(p2->m_n2))));
        }

        /**
         * @brief Defines ordering between two edges (based on ID of edges).
         */
        bool operator< (EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
        {
            return p1->m_id < p2->m_id;
        }

        /**
         * @brief Defines equality of two faces (equal if IDs of vertices are
         * the same.)
         */
        bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2)
        {
            std::vector<NodeSharedPtr>::iterator it1, it2;
            for (it1 = p1->m_vertexList.begin(); it1 != p1->m_vertexList.end(); ++it1)
            {
                if (find(p2->m_vertexList.begin(), p2->m_vertexList.end(), *it1)
                    == p2->m_vertexList.end())
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
            return p1->m_id < p2->m_id;
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
            NodeSharedPtr vOld = m_vertex[p];
            m_vertex[p] = pNew;
            for (unsigned int i = 0; i < m_edge.size(); ++i)
            {
                if (m_edge[i]->m_n1 == vOld)
                {
                    m_edge[i]->m_n1 = pNew;
                }
                else if (m_edge[i]->m_n2 == vOld)
                {
                    m_edge[i]->m_n2 = pNew;
                }
            }
            for (unsigned int i = 0; i < m_face.size(); ++i)
            {
                // Replace vertices in faces
                for (unsigned int j = 0; j < m_face[i]->m_vertexList.size(); ++j)
                {
                    if (m_face[i]->m_vertexList[j] == vOld)
                    {
                        m_face[i]->m_vertexList[j] = pNew;
                    }
                }
                for (unsigned int j = 0; j < m_face[i]->m_edgeList.size(); ++j)
                {
                    if (m_face[i]->m_edgeList[j]->m_n1 == vOld)
                    {
                        m_face[i]->m_edgeList[j]->m_n1 = pNew;
                    }
                    else if (m_face[i]->m_edgeList[j]->m_n2 == vOld)
                    {
                        m_face[i]->m_edgeList[j]->m_n2 = pNew;
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
            EdgeSharedPtr vOld = m_edge[p];
            m_edge[p] = pNew;
            for (unsigned int i = 0; i < m_face.size(); ++i)
            {
                for (unsigned int j = 0; j < m_face[i]->m_edgeList.size(); ++j)
                {
                    if (m_face[i]->m_edgeList[j] == vOld)
                    {
                        m_face[i]->m_edgeList[j] = pNew;
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
            m_face[p] = pNew;
        }
        
        /**
         * @brief Obtain the order of an element by looking at edges.
         */
        int Element::GetMaxOrder()
        {
            int i, ret = 1;
            
            for (i = 0; i < m_edge.size(); ++i)
            {
                int edgeOrder = m_edge[i]->GetNodeCount()-1;
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

#if 0 // turn this option off since causes problem with InputNekpp.cpp 
            if (doSort)
            {
                element_id_less_than sortOperator;
                sort(m_items.begin(), m_items.end(), sortOperator);
            }
#endif

            stringstream st;
            vector<ElementSharedPtr>::iterator it;
            bool range = false;
            int vId = m_items[0]->GetId();
            int prevId = vId;

            st << " " << m_tag << "[" << vId;

            for (it = m_items.begin()+1; it != m_items.end(); ++it){
                // store previous element ID and get current one
                prevId = vId;
                vId = (*it)->GetId();

                // continue an already started range
                if (prevId > -1 && vId == prevId + 1)
                {
                    range = true;
                    // if this is the last element, it's the end of a range, so write
                    if (*it == m_items.back())
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

        LibUtilities::ShapeType Point::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::ePoint, Point::create, "Point");
        
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
            m_vertex.push_back(pNodeList[0]);
        }

        /**
         * @brief Return the number of nodes defining a point (i.e. return 1).
         */
        unsigned int Point::GetNumNodes(ElmtConfig pConf)
        {
            return 1;
        }


        LibUtilities::ShapeType Line::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::eSegment, Line::create, "Line");
        
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
            int n     = m_conf.m_order-1;
            
            // Add vertices
            for (int i = 0; i < 2; ++i) {
                m_vertex.push_back(pNodeList[i]);
            }
            vector<NodeSharedPtr> edgeNodes;
            if (m_conf.m_order > 1) {
                for (int j = 0; j<n; ++j) {
                    edgeNodes.push_back(pNodeList[2+j]);
                }
            }
            m_edge.push_back(boost::shared_ptr<Edge>(
                new Edge(pNodeList[0], pNodeList[1], edgeNodes, m_conf.m_edgeCurveType)));
        }
        
        SpatialDomains::GeometrySharedPtr Line::GetGeom(int coordDim)
        {
            // Create edge vertices.
            SpatialDomains::PointGeomSharedPtr p[2];
            SpatialDomains::SegGeomSharedPtr   ret;

            p[0] = m_vertex[0]->GetGeom(coordDim);
            p[1] = m_vertex[1]->GetGeom(coordDim);
            
            if (m_edge[0]->m_edgeNodes.size() > 0)
            {
                SpatialDomains::CurveSharedPtr c = 
                    MemoryManager<SpatialDomains::Curve>::
                    AllocateSharedPtr(m_id, m_edge[0]->m_curveType);
                
                c->m_points.push_back(p[0]);
                for (int i = 0; i < m_edge[0]->m_edgeNodes.size(); ++i)
                {
                    c->m_points.push_back(m_edge[0]->m_edgeNodes[i]->GetGeom(coordDim));
                }
                c->m_points.push_back(p[1]);
                
                ret = MemoryManager<SpatialDomains::SegGeom>::
                    AllocateSharedPtr(m_id, 2, p, c);
            }
            else
            {
                ret = MemoryManager<SpatialDomains::SegGeom>::
                    AllocateSharedPtr(m_id, 2, p);
            }

            return ret;
        }

        /**
         * @brief Return the number of nodes defining a line.
         */
        unsigned int Line::GetNumNodes(ElmtConfig pConf)
        {
            return pConf.m_order+1;
        }


        LibUtilities::ShapeType Triangle::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::eTriangle, Triangle::create, "Triangle");
        
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
            m_curveType = LibUtilities::eNodalTriEvenlySpaced;
            int n     = m_conf.m_order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge. This is based on the ordering produced by
            // gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 4;
            edgeNodeMap[pair<int,int>(2,3)] = 4 + n;
            edgeNodeMap[pair<int,int>(3,1)] = 4 + 2*n;

            // Add vertices. This logic will determine (in 2D) whether the
            // element is clockwise (sum > 0) or counter-clockwise (sum < 0).
            NekDouble sum = 0.0;
            for (int i = 0; i < 3; ++i) {
                int o = (i+1) % 3;
                m_vertex.push_back(pNodeList[i]);
                sum += (pNodeList[o]->m_x - pNodeList[i]->m_x) *
                       (pNodeList[o]->m_y + pNodeList[i]->m_y);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes,
                                                        m_conf.m_edgeCurveType)));
            }

            if (pConf.m_reorient)
            {
                if (sum > 0.0)
                {
                    reverse(m_edge.begin(), m_edge.end());
                }
            }

            if (m_conf.m_faceNodes)
            {
                m_volumeNodes.insert(m_volumeNodes.begin(),
                                   pNodeList.begin()+3*m_conf.m_order,
                                   pNodeList.end());
            }
        }

        SpatialDomains::GeometrySharedPtr Triangle::GetGeom(int coordDim)
        {
            SpatialDomains::SegGeomSharedPtr   edges[3];
            SpatialDomains::PointGeomSharedPtr verts[3];
            SpatialDomains::TriGeomSharedPtr   ret;
            
            for (int i = 0; i < 3; ++i)
            {
                edges[i] = m_edge  [i]->GetGeom(coordDim);
                verts[i] = m_vertex[i]->GetGeom(coordDim);
            }
            
            StdRegions::Orientation edgeorient[3] = {
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[0])
            };
            
            ret = MemoryManager<SpatialDomains::TriGeom>::
                AllocateSharedPtr(m_id, verts, edges, edgeorient);

            return ret;
        }
        
        /**
         * @brief Return the number of nodes defining a triangle.
         */
        unsigned int Triangle::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            if (!pConf.m_faceNodes)
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
                Array<OneD, NekDouble> coeffs(nodalTri->GetNcoeffs());
                tri->FwdTrans(alloc+i*nqtot, coeffs);
                // Apply Vandermonde matrix to project onto nodal space.
                nodalTri->ModalToNodal(coeffs, tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 3; ++i)
            {
                int pos = 3 + i*(order-1);
                m_edge[i]->m_edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    m_edge[i]->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Extract face-interior nodes.
            int pos = 3 + 3*(order-1);
            for (i = pos; i < (order+1)*(order+2)/2; ++i)
            {
                m_volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.m_order       = order;
            m_conf.m_faceNodes   = true;
            m_conf.m_volumeNodes = true;
        }


        LibUtilities::ShapeType Quadrilateral::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::eQuadrilateral, Quadrilateral::create, 
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
            int n = m_conf.m_order-1;

            // Create a map to relate edge nodes to a pair of vertices
            // defining an edge. This is based on the ordering produced by
            // gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 5;
            edgeNodeMap[pair<int,int>(2,3)] = 5 + n;
            edgeNodeMap[pair<int,int>(3,4)] = 5 + 2*n;
            edgeNodeMap[pair<int,int>(4,1)] = 5 + 3*n;

            // Add vertices. This logic will determine (in 2D) whether the
            // element is clockwise (sum > 0) or counter-clockwise (sum < 0).
            NekDouble sum = 0.0;
            for (int i = 0; i < 4; ++i) {
                int o = (i+1) % 4;
                m_vertex.push_back(pNodeList[i]);
                sum += (pNodeList[o]->m_x - pNodeList[i]->m_x) *
                       (pNodeList[o]->m_y + pNodeList[i]->m_y);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes,
                                                        m_conf.m_edgeCurveType)));
            }

            if (pConf.m_reorient)
            {
                if (sum > 0.0)
                {
                    reverse(m_edge.begin(), m_edge.end());
                }
            }

            if (m_conf.m_faceNodes)
            {
                m_volumeNodes.insert(m_volumeNodes.begin(), 
                                   pNodeList.begin()+4*m_conf.m_order,
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
                m_edge[i]->m_edgeNodes.clear();
                for (int j = 1; j < order; ++j, pos += edgeMap[i][1])
                {
                    m_edge[i]->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, x[pos], y[pos], z[pos])));
                }
            }

            // Extract face-interior nodes.
            m_volumeNodes.clear();
            for (int i = 1; i < order; ++i)
            {
                int pos = i*(order+1);
                for (int j = 1; j < order; ++j)
                {
                    m_volumeNodes.push_back(
                        NodeSharedPtr(new Node(0, x[pos+j], y[pos+j], z[pos+j])));
                }                
            }
            
            m_conf.m_order       = order;
            m_conf.m_faceNodes   = true;
            m_conf.m_volumeNodes = true;
        }

        SpatialDomains::GeometrySharedPtr Quadrilateral::GetGeom(int coordDim)
        {
            SpatialDomains::SegGeomSharedPtr   edges[4];
            SpatialDomains::PointGeomSharedPtr verts[4];
            SpatialDomains::QuadGeomSharedPtr  ret;
            
            for (int i = 0; i < 4; ++i)
            {
                edges[i] = m_edge  [i]->GetGeom(coordDim);
                verts[i] = m_vertex[i]->GetGeom(coordDim);
            }
            
            StdRegions::Orientation edgeorient[4] = {
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[0], *edges[1]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[1], *edges[2]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[2], *edges[3]),
                SpatialDomains::SegGeom::GetEdgeOrientation(*edges[3], *edges[0])
            };

            ret = MemoryManager<SpatialDomains::QuadGeom>::
                AllocateSharedPtr(m_id, verts, edges, edgeorient);

            return ret;
        }

        /**
         * @brief Return the number of nodes defining a quadrilateral.
         */
        unsigned int Quadrilateral::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            if (!pConf.m_faceNodes)
                return 4*n;
            else
                return (n+1)*(n+1);
        }
        
        
        LibUtilities::ShapeType Tetrahedron::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::eTetrahedron, Tetrahedron::create, "Tetrahedron");

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
            int n = m_conf.m_order-1;

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
            for (int i = 0; i < 4; ++i)
            {
                m_vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            int eid = 0;
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(
                    EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                           pNodeList[it->first.second-1],
                                           edgeNodes,
                                           m_conf.m_edgeCurveType)));
                m_edge.back()->m_id = eid++;
            }
            
            // Reorient the tet to ensure collapsed coordinates align between
            // adjacent elements.
            if (m_conf.m_reorient)
            {
                OrientTet();
            }
            else
            {
                // If we didn't need to orient the tet then set up the
                // orientation map as the identity mapping.
                for (int i = 0; i < 4; ++i)
                {
                    m_orientationMap[i] = i;
                }
            }
            
            // Face-vertex IDs
            int face_ids[4][3] = {
                {0,1,2}, {0,1,3}, {1,2,3}, {0,2,3}
            };

            // Face-edge IDs
            int face_edges[4][3];

            // Create faces
            for (int j = 0; j < 4; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;

                // Extract the edges for this face. We need to do this because
                // of the reorientation which might have been applied (see the
                // additional note below).
                for (int k = 0; k < 3; ++k)
                {
                    faceVertices.push_back(m_vertex[face_ids[j][k]]);
                    NodeSharedPtr a = m_vertex[face_ids[j][k]];
                    NodeSharedPtr b = m_vertex[face_ids[j][(k+1)%3]];
                    for (unsigned int i = 0; i < m_edge.size(); ++i)
                    {
                        if ( ((*(m_edge[i]->m_n1)==*a) && (*(m_edge[i]->m_n2)==*b))
                                || ((*(m_edge[i]->m_n1)==*b) && (*(m_edge[i]->m_n2) == *a)) )
                        {
                            face_edges[j][k] = i;
                            faceEdges.push_back(m_edge[i]);
                            break;
                        }
                    }
                }

                // When face curvature is supplied, it may have been the case
                // that our tetrahedron was reoriented. In this case, faces have
                // different vertex IDs and so we have to rotate the face
                // curvature so that the two align appropriately.
                if (m_conf.m_faceNodes)
                {
                    const int nFaceNodes = n*(n-1)/2;

                    // Get the vertex IDs of whatever face we are processing.
                    vector<int> faceIds(3);
                    faceIds[0] = faceVertices[0]->m_id;
                    faceIds[1] = faceVertices[1]->m_id;
                    faceIds[2] = faceVertices[2]->m_id;

                    // Find out the original face number as we were given it in
                    // the constructor using the orientation map.
                    int origFace = -1;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (m_orientationMap[i] == j)
                        {
                            origFace = i;
                            break;
                        }
                    }

                    ASSERTL0(origFace >= 0, "Couldn't find face");

                    // Now get the face nodes for the original face.
                    int N = 4 + 6*n + origFace * nFaceNodes;
                    for (int i = 0; i < nFaceNodes; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }

                    // Find the original face vertex IDs.
                    vector<int> origFaceIds(3);
                    origFaceIds[0] = pNodeList[face_ids[origFace][0]]->m_id;
                    origFaceIds[1] = pNodeList[face_ids[origFace][1]]->m_id;
                    origFaceIds[2] = pNodeList[face_ids[origFace][2]]->m_id;

                    // Construct a HOTriangle object which performs the
                    // orientation magically for us.
                    HOTriangle<NodeSharedPtr> hoTri(origFaceIds, faceNodes);
                    hoTri.Align(faceIds);

                    // Copy the face nodes back again.
                    faceNodes = hoTri.surfVerts;
                }

                m_face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
            }

            vector<EdgeSharedPtr> tmp(6);
            tmp[0] = m_edge[face_edges[0][0]];
            tmp[1] = m_edge[face_edges[0][1]];
            tmp[2] = m_edge[face_edges[0][2]];
            tmp[3] = m_edge[face_edges[1][2]];
            tmp[4] = m_edge[face_edges[1][1]];
            tmp[5] = m_edge[face_edges[2][1]];
            m_edge = tmp;
        }
        
        SpatialDomains::GeometrySharedPtr Tetrahedron::GetGeom(int coordDim)
        {
            SpatialDomains::TriGeomSharedPtr tfaces[4];
            SpatialDomains::TetGeomSharedPtr ret;
            
            for (int i = 0; i < 4; ++i)
            {
                tfaces[i] = boost::dynamic_pointer_cast
                    <SpatialDomains::TriGeom>(m_face[i]->GetGeom(coordDim));
            }

            ret = MemoryManager<SpatialDomains::TetGeom>::
                AllocateSharedPtr(tfaces);

            return ret;
        }
        
        /**
         * @brief Return the number of nodes defining a tetrahedron.
         */
        unsigned int Tetrahedron::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            if (pConf.m_volumeNodes && pConf.m_faceNodes)
                return (n+1)*(n+2)*(n+3)/6;
            else if (!pConf.m_volumeNodes && pConf.m_faceNodes)
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
                Array<OneD, NekDouble> coeffs(nodalTet->GetNcoeffs());
                tet->FwdTrans(alloc+i*nqtot, coeffs);
                // Apply Vandermonde matrix to project onto nodal space.
                nodalTet->ModalToNodal(coeffs, tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 6; ++i)
            {
                int pos = 4 + i*(order-1);
                m_edge[i]->m_edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    m_edge[i]->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Now extract face-interior nodes.
            for (i = 0; i < 4; ++i)
            {
                int pos = 4 + 6*(order-1) + i*(order-2)*(order-1)/2;
                m_face[i]->m_faceNodes.clear();
                for (j = 0; j < (order-2)*(order-1)/2; ++j)
                {
                    m_face[i]->m_faceNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }
            
            // Finally extract volume nodes.
            int pos = 4 + 6*(order-1) + 4*(order-2)*(order-1)/2;
            for (i = pos; i < (order+1)*(order+2)*(order+3)/6; ++i)
            {
                m_volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.m_order       = order;
            m_conf.m_faceNodes   = true;
            m_conf.m_volumeNodes = true;
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

                nodes[0] = m_vertex[face_ids[i][0]]->m_id;
                nodes[1] = m_vertex[face_ids[i][1]]->m_id;
                nodes[2] = m_vertex[face_ids[i][2]]->m_id;
                
                sort(nodes.begin(), nodes.end());
                struct TetOrient faceNodes(nodes, i);
                faces.insert(faceNodes);
            }

            // Store a copy of the original vertex ordering so we can create a
            // permutation map later.
            vector<NodeSharedPtr> origVert = m_vertex;

            // Order vertices with highest global vertex at top degenerate
            // point. Place second highest global vertex at base degenerate
            // point.
            sort(m_vertex.begin(), m_vertex.end());
            
            // Calculate a.(b x c) if negative, reverse order of
            // non-degenerate points to correctly orientate the tet.

            NekDouble ax  = m_vertex[1]->m_x-m_vertex[0]->m_x;
            NekDouble ay  = m_vertex[1]->m_y-m_vertex[0]->m_y;
            NekDouble az  = m_vertex[1]->m_z-m_vertex[0]->m_z;
            NekDouble bx  = m_vertex[2]->m_x-m_vertex[0]->m_x;
            NekDouble by  = m_vertex[2]->m_y-m_vertex[0]->m_y;
            NekDouble bz  = m_vertex[2]->m_z-m_vertex[0]->m_z;
            NekDouble cx  = m_vertex[3]->m_x-m_vertex[0]->m_x;
            NekDouble cy  = m_vertex[3]->m_y-m_vertex[0]->m_y;
            NekDouble cz  = m_vertex[3]->m_z-m_vertex[0]->m_z;

            NekDouble nx = (ay*bz-az*by);
            NekDouble ny = (az*bx-ax*bz);
            NekDouble nz = (ax*by-ay*bx);
            NekDouble nmag = sqrt(nx*nx+ny*ny+nz*nz);
            nx /= nmag; ny /= nmag; nz /= nmag; 

            // distance of top vertex from base 
            NekDouble dist = cx*nx+cy*ny+cz*nz;

            if (fabs(dist) <= 1e-4)
            {
                cerr << "Warning: degenerate tetrahedron, 3rd vertex is = " << dist <<" from face" << endl;
            }

            if (dist < 0)
            {
                swap(m_vertex[0], m_vertex[1]);
            }

            nx = (ay*cz-az*cy);
            ny = (az*cx-ax*cz);
            nz = (ax*cy-ay*cx);
            nmag = sqrt(nx*nx+ny*ny+nz*nz);
            nx /= nmag; ny /= nmag; nz /= nmag; 
            
            // distance of top vertex from base 
            dist = bx*nx+by*ny+bz*nz;

            if (fabs(dist) <= 1e-4)
            {
                cerr << "Warning: degenerate tetrahedron, 2nd vertex is = " << dist <<" from face" << endl;
            }

            nx = (by*cz-bz*cy);
            ny = (bz*cx-bx*cz);
            nz = (bx*cy-by*cx);
            nmag = sqrt(nx*nx+ny*ny+nz*nz);
            nx /= nmag; ny /= nmag; nz /= nmag; 
            
            // distance of top vertex from base 
            dist = ax*nx+ay*ny+az*nz;

            if (fabs(dist) <= 1e-4)
            {
                cerr << "Warning: degenerate tetrahedron, 1st vertex is = " << dist <<" from face" << endl;
            }

            
            TetOrientSet::iterator it;
            
            // Search for the face in the original set of face nodes. Then use
            // this to construct the #orientationMap.
            for (int i = 0; i < 4; ++i)
            {
                vector<int> nodes(3);
                
                nodes[0] = m_vertex[face_ids[i][0]]->m_id;
                nodes[1] = m_vertex[face_ids[i][1]]->m_id;
                nodes[2] = m_vertex[face_ids[i][2]]->m_id;
                sort(nodes.begin(), nodes.end());
                
                struct TetOrient faceNodes(nodes, 0);
                
                it = faces.find(faceNodes);
                m_orientationMap[it->fid] = i;

                for (int j = 0; j < 4; ++j)
                {
                    if (m_vertex[i]->m_id == origVert[j]->m_id)
                    {
                        m_origVertMap[j] = i;
                        break;
                    }
                }
            }
        }

        LibUtilities::ShapeType Pyramid::type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::ePyramid, Pyramid::create, "Pyramid");

        /**
         * @brief Create a pyramidic element.
         */
        Pyramid::Pyramid(ElmtConfig            pConf,
                         vector<NodeSharedPtr> pNodeList,
                         vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size())
        {
            m_tag     = "P";
            m_dim     = 3;
            m_taglist = pTagList;
            int n     = m_conf.m_order-1;

            // This edge-node map is based on Nektar++ ordering.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 6;
            edgeNodeMap[pair<int,int>(2,3)] = 6 + n;
            edgeNodeMap[pair<int,int>(4,3)] = 6 + 2*n;
            edgeNodeMap[pair<int,int>(1,4)] = 6 + 3*n;
            edgeNodeMap[pair<int,int>(1,5)] = 6 + 4*n;
            edgeNodeMap[pair<int,int>(2,5)] = 6 + 5*n;
            edgeNodeMap[pair<int,int>(3,5)] = 6 + 6*n;
            edgeNodeMap[pair<int,int>(4,5)] = 6 + 7*n;
            
            // Add vertices
            for (int i = 0; i < 5; ++i)
            {
                m_vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            int eid = 0;
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1)
                {
                    for (int j = it->second; j < it->second + n; ++j)
                    {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(
                    EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                           pNodeList[it->first.second-1],
                                           edgeNodes,
                                           m_conf.m_edgeCurveType)));
                m_edge.back()->m_id = eid++;
            }

            // Create faces
            int face_ids[5][4] = {
                {0,1,2,3}, {0,1,4,-1}, {1,2,4,-1}, {3,2,4,-1}, {0,3,4,-1}
            };
            int face_edges[5][4];
            int faceoffset = 5 + 8*n;
            for (int j = 0; j < 5; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                int nEdge = j > 0 ? 3 : 4;

                for (int k = 0; k < nEdge; ++k)
                {
                    faceVertices.push_back(m_vertex[face_ids[j][k]]);
                    NodeSharedPtr a = m_vertex[face_ids[j][k]];
                    NodeSharedPtr b = m_vertex[face_ids[j][(k+1) % nEdge]];
                    for (unsigned int i = 0; i < m_edge.size(); ++i)
                    {
                        if ((m_edge[i]->m_n1 == a && m_edge[i]->m_n2 == b) ||
                            (m_edge[i]->m_n1 == b && m_edge[i]->m_n2 == a))
                        {
                            faceEdges.push_back(m_edge[i]);
                            face_edges[j][k] = i;
                            break;
                        }
                    }
                }

                if (m_conf.m_faceNodes)
                {
                    int facenodes = j == 0 ? n*n : n*(n-1)/2;
                    for (int i = 0; i < facenodes; ++i)
                    {
                        faceNodes.push_back(pNodeList[faceoffset+i]);
                    }
                    faceoffset   += facenodes;
                }
                m_face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
            }

            // Reorder edges to align with Nektar++ order.
            vector<EdgeSharedPtr> tmp(8);
            tmp[0] = m_edge[face_edges[0][0]];
            tmp[1] = m_edge[face_edges[0][1]];
            tmp[2] = m_edge[face_edges[0][2]];
            tmp[3] = m_edge[face_edges[0][3]];
            tmp[4] = m_edge[face_edges[1][2]];
            tmp[5] = m_edge[face_edges[1][1]];
            tmp[6] = m_edge[face_edges[3][1]];
            tmp[7] = m_edge[face_edges[3][2]];
            m_edge = tmp;
        }
        
        SpatialDomains::GeometrySharedPtr Pyramid::GetGeom(int coordDim)
        {
            SpatialDomains::Geometry2DSharedPtr faces[5];

            for (int i = 0; i < 5; ++i)
            {
                faces[i] = m_face[i]->GetGeom(coordDim);
            }

            m_geom = MemoryManager<SpatialDomains::PyrGeom>::
                AllocateSharedPtr(faces);

            return m_geom;
        }
        
        /**
         * @brief Return the number of nodes defining a pyramid.
         */
        unsigned int Pyramid::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            return 5 + 8*(n-1);
        }

        LibUtilities::ShapeType Prism::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::ePrism, Prism::create, "Prism");
        
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
            int n     = m_conf.m_order-1;

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
                m_vertex.push_back(pNodeList[i]);
            }

            int eid = 0;
            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(EdgeSharedPtr(
                    new Edge(pNodeList[it->first.first-1],
                             pNodeList[it->first.second-1],
                             edgeNodes,
                             m_conf.m_edgeCurveType)));
                m_edge.back()->m_id = eid++;
            }
            
            if (m_conf.m_reorient)
            {
                OrientPrism();
            }
            else
            {
                m_orientation = 0;
            }

            // Create faces
            int face_ids[5][4] = {
                {0,1,2,3},{0,1,4,-1},{1,2,5,4},{3,2,5,-1},{0,3,5,4}};
            int face_edges[5][4];

            int face_offset[5];
            face_offset[0] = 6 + 9*n;
            for (int j = 0; j < 4; ++j)
            {
                int facenodes = j % 2 == 0 ? n*n : n*(n-1)/2;
                face_offset[j+1] = face_offset[j] + facenodes;
            }

            for (int j = 0; j < 5; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                int nEdge = 3 - (j % 2 - 1);

                for (int k = 0; k < nEdge; ++k)
                {
                    faceVertices.push_back(m_vertex[face_ids[j][k]]);
                    NodeSharedPtr a = m_vertex[face_ids[j][k]];
                    NodeSharedPtr b = m_vertex[face_ids[j][(k+1) % nEdge]];
                    unsigned int i;
                    for (i = 0; i < m_edge.size(); ++i)
                    {
                        if ((m_edge[i]->m_n1->m_id == a->m_id 
                             && m_edge[i]->m_n2->m_id == b->m_id) || 
                            (m_edge[i]->m_n1->m_id == b->m_id 
                             && m_edge[i]->m_n2->m_id == a->m_id))
                        {
                            faceEdges.push_back(m_edge[i]);
                            face_edges[j][k] = i;
                            break;
                        }
                    }
                    
                    if(i == m_edge.size())
                    {
                        face_edges[j][k] = -1;
                    }
                }

                if (m_conf.m_faceNodes)
                {
                    int face = j, facenodes;

                    if (j % 2 == 0)
                    {
                        facenodes = n*n;
                        if (m_orientation == 1)
                        {
                            face = (face+4) % 6;
                        }
                        else if (m_orientation == 2)
                        {
                            face = (face+2) % 6;
                        }
                    }
                    else
                    {
                        // TODO: need to rotate these too.
                        facenodes = n*(n-1)/2;
                    }

                    for (int i = 0; i < facenodes; ++i)
                    {
                        faceNodes.push_back(pNodeList[face_offset[face]+i]);
                    }
                }
                m_face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
            }

            // Re-order edge array to be consistent with Nektar++ ordering.
            vector<EdgeSharedPtr> tmp(9);
            ASSERTL1(face_edges[0][0] != -1,"face_edges[0][0] == -1");
            tmp[0] = m_edge[face_edges[0][0]];
            ASSERTL1(face_edges[0][1] != -1,"face_edges[0][1] == -1");
            tmp[1] = m_edge[face_edges[0][1]];
            ASSERTL1(face_edges[0][2] != -1,"face_edges[0][2] == -1");
            tmp[2] = m_edge[face_edges[0][2]];
            ASSERTL1(face_edges[0][3] != -1,"face_edges[0][3] == -1");
            tmp[3] = m_edge[face_edges[0][3]];
            ASSERTL1(face_edges[1][2] != -1,"face_edges[1][2] == -1");
            tmp[4] = m_edge[face_edges[1][2]];
            ASSERTL1(face_edges[1][1] != -1,"face_edges[1][1] == -1");
            tmp[5] = m_edge[face_edges[1][1]];
            ASSERTL1(face_edges[2][1] != -1,"face_edges[2][1] == -1");
            tmp[6] = m_edge[face_edges[2][1]];
            ASSERTL1(face_edges[3][2] != -1,"face_edges[3][2] == -1");
            tmp[7] = m_edge[face_edges[3][2]];
            ASSERTL1(face_edges[4][2] != -1,"face_edges[4][2] == -1");
            tmp[8] = m_edge[face_edges[4][2]];
            m_edge = tmp;
        }

        /**
         * @brief Return the number of nodes defining a prism.
         */
        unsigned int Prism::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            if (pConf.m_faceNodes && pConf.m_volumeNodes)
                return (n+1)*(n+1)*(n+2)/2;
            else if (pConf.m_faceNodes && !pConf.m_volumeNodes)
                return 3*(n+1)*(n+1)+2*(n+1)*(n+2)/2-9*(n+1)+6;
            else
                return 9*(n+1)-12;
        }

        SpatialDomains::GeometrySharedPtr Prism::GetGeom(int coordDim)
        {
            SpatialDomains::Geometry2DSharedPtr faces[5];
            SpatialDomains::PrismGeomSharedPtr  ret;
            
            for (int i = 0; i < 5; ++i)
            {
                faces[i] = m_face[i]->GetGeom(coordDim);
            }

            ret = MemoryManager<SpatialDomains::PrismGeom>::
                AllocateSharedPtr(faces);

            return ret;
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
                Array<OneD, NekDouble> coeffs(nodalPrism->GetNcoeffs());
                prism->FwdTrans(alloc+i*nqtot, coeffs);
                // Apply Vandermonde matrix to project onto nodal space.
                nodalPrism->ModalToNodal(coeffs, tmp=alloc+(i+3)*nqtot);
            }
            
            // Now extract points from the co-ordinate arrays into the
            // edge/face/volume nodes. First, extract edge-interior nodes.
            for (i = 0; i < 9; ++i)
            {
                pos = 6 + i*(order-1);
                m_edge[i]->m_edgeNodes.clear();
                for (j = 0; j < order-1; ++j)
                {
                    m_edge[i]->m_edgeNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
            }

            // Now extract face-interior nodes.
            pos = 6 + 9*(order-1);
            for (i = 0; i < 5; ++i)
            {
                int facesize = i % 2 ? (order-2)*(order-1)/2 : (order-1)*(order-1);
                m_face[i]->m_faceNodes.clear();
                for (j = 0; j < facesize; ++j)
                {
                    m_face[i]->m_faceNodes.push_back(
                        NodeSharedPtr(new Node(0, xo[pos+j], yo[pos+j], zo[pos+j])));
                }
                pos += facesize;
            }
            
            // Finally extract volume nodes.
            for (i = pos; i < (order+1)*(order+1)*(order+2)/2; ++i)
            {
                m_volumeNodes.push_back(
                    NodeSharedPtr(new Node(0, xo[i], yo[i], zo[i])));
            }
            
            m_conf.m_order       = order;
            m_conf.m_faceNodes   = true;
            m_conf.m_volumeNodes = true;
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
                gid[i] = m_vertex[i]->m_id;
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
                m_orientation = 0;
            } 
            else if (lid[0] == 1 || lid[0] == 2) 
            {
                // Rotate prism clockwise in p-r plane
                vector<NodeSharedPtr> vertexmap(6);
                vertexmap[0] = m_vertex[4];
                vertexmap[1] = m_vertex[0];
                vertexmap[2] = m_vertex[3];
                vertexmap[3] = m_vertex[5];
                vertexmap[4] = m_vertex[1];
                vertexmap[5] = m_vertex[2];
                m_vertex = vertexmap;
                m_orientation = 1;
            }
            else if (lid[0] == 0 || lid[0] == 3)
            {
                // Rotate prism counter-clockwise in p-r plane
                vector<NodeSharedPtr> vertexmap(6);
                vertexmap[0] = m_vertex[1];
                vertexmap[1] = m_vertex[4];
                vertexmap[2] = m_vertex[5];
                vertexmap[3] = m_vertex[2];
                vertexmap[4] = m_vertex[0];
                vertexmap[5] = m_vertex[3];
                m_vertex = vertexmap;
                m_orientation = 2;
            }
            else
            {
                cerr << "Warning: possible prism orientation problem." << endl;
            }
        }


        LibUtilities::ShapeType Hexahedron::m_type = GetElementFactory().
            RegisterCreatorFunction(LibUtilities::eHexahedron, Hexahedron::create, "Hexahedron");

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
            int n = m_conf.m_order-1;
            
            // Create a map to relate edge nodes to a pair of vertices defining an edge
            // This is based on the ordering produced by gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 9;
            edgeNodeMap[pair<int,int>(2,3)] = 9 + n;
            edgeNodeMap[pair<int,int>(3,4)] = 9 + 2*n;
            edgeNodeMap[pair<int,int>(4,1)] = 9 + 3*n;
            edgeNodeMap[pair<int,int>(1,5)] = 9 + 4*n;
            edgeNodeMap[pair<int,int>(2,6)] = 9 + 5*n;
            edgeNodeMap[pair<int,int>(3,7)] = 9 + 6*n;
            edgeNodeMap[pair<int,int>(4,8)] = 9 + 7*n;
            edgeNodeMap[pair<int,int>(5,6)] = 9 + 8*n;
            edgeNodeMap[pair<int,int>(6,7)] = 9 + 9*n;
            edgeNodeMap[pair<int,int>(7,8)] = 9 + 10*n;
            edgeNodeMap[pair<int,int>(8,5)] = 9 + 11*n;

            // Add vertices
            for (int i = 0; i < 8; ++i) {
                m_vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (m_conf.m_order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                m_edge.push_back(EdgeSharedPtr(new Edge(pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes,
                                                        m_conf.m_edgeCurveType)));
            }

            // Create faces
            int face_edges[6][4];
            int face_ids[6][4] = {
                {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7}};
            for (int j = 0; j < 6; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                for (int k = 0; k < 4; ++k)
                {
                    faceVertices.push_back(m_vertex[face_ids[j][k]]);
                    NodeSharedPtr a = m_vertex[face_ids[j][k]];
                    NodeSharedPtr b = m_vertex[face_ids[j][(k+1)%4]];
                    for (unsigned int i = 0; i < m_edge.size(); ++i)
                    {
                        if ( ((*(m_edge[i]->m_n1)==*a) && (*(m_edge[i]->m_n2)==*b))
                                || ((*(m_edge[i]->m_n1)==*b) && (*(m_edge[i]->m_n2) == *a)) )
                        {
                            face_edges[j][k] = i;
                            faceEdges.push_back(m_edge[i]);
                            break;
                        }
                    }
                }

                if (m_conf.m_faceNodes)
                {
                    int N = 8 + 12*n + j*n*n;
                    for (int i = 0; i < n*n; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                m_face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.m_faceCurveType)));
            }

            // Reorder edges to be consistent with Nektar++ ordering.
            vector<EdgeSharedPtr> tmp(12);
            tmp[ 0] = m_edge[face_edges[0][0]];
            tmp[ 1] = m_edge[face_edges[0][1]];
            tmp[ 2] = m_edge[face_edges[0][2]];
            tmp[ 3] = m_edge[face_edges[0][3]];
            tmp[ 4] = m_edge[face_edges[1][3]];
            tmp[ 5] = m_edge[face_edges[1][1]];
            tmp[ 6] = m_edge[face_edges[2][1]];
            tmp[ 7] = m_edge[face_edges[4][1]];
            tmp[ 8] = m_edge[face_edges[5][0]];
            tmp[ 9] = m_edge[face_edges[5][1]];
            tmp[10] = m_edge[face_edges[5][2]];
            tmp[11] = m_edge[face_edges[5][3]];
            m_edge = tmp;
        }
        
        SpatialDomains::GeometrySharedPtr Hexahedron::GetGeom(int coordDim)
        {
            SpatialDomains::QuadGeomSharedPtr faces[6];
            SpatialDomains::HexGeomSharedPtr  ret;
            
            for (int i = 0; i < 6; ++i)
            {
                faces[i] = boost::dynamic_pointer_cast
                    <SpatialDomains::QuadGeom>(m_face[i]->GetGeom(coordDim));
            }

            ret = MemoryManager<SpatialDomains::HexGeom>::
                AllocateSharedPtr(faces);

            return ret;
        }

        /**
         * @brief Return the number of nodes defining a hexahedron.
         */
        unsigned int Hexahedron::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.m_order;
            if (pConf.m_faceNodes && pConf.m_volumeNodes)
                return (n+1)*(n+1)*(n+1);
            else if (pConf.m_faceNodes && !pConf.m_volumeNodes)
                return 6*(n+1)*(n+1)-12*(n+1)+8;
            else
                return 12*(n+1)-16;
        }
    }
}
