////////////////////////////////////////////////////////////////////////////////
//
//  File: MeshElements.h
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

#ifndef UTILITY_PREPROCESSING_MESHCONVERT_MESHELEMENTS
#define UTILITY_PREPROCESSING_MESHCONVERT_MESHELEMENTS

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/ShapeType.hpp>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/MeshComponents.h>

namespace Nektar
{
    namespace Utilities
    {
        // Forwards declaration for Element class.
        class Element;
        /// Shared pointer to an element.
        typedef boost::shared_ptr<Element> ElementSharedPtr;

        /**
         * @brief Represents a point in the domain.
         *
         * Such points may either be element vertices, or simply control
         * points on high-order edges/faces, although this information is not
         * contained within this class.
         */
        class Node {
        public:
            /// Create a new node at a specified coordinate.
            Node(int pId, NekDouble pX, NekDouble pY, NekDouble pZ)
                : m_id(pId), m_x(pX), m_y(pY), m_z(pZ), m_geom() {}
            /// Copy an existing node.
            Node(const Node& pSrc)
                : m_id(pSrc.m_id), m_x(pSrc.m_x), m_y(pSrc.m_y),
                  m_z(pSrc.m_z), m_geom() {}
            Node() : m_id(0), m_x(0.0), m_y(0.0), m_z(0.0), m_geom() {}
            ~Node() {}

            /// Reset the local id;
            void SetID(int pId)
            {
                m_id = pId;
            }

            /// Get the local id;
            int GetID(void)
            {
                return m_id;
            }

            /// Define node ordering based on ID.
            bool operator<(const Node& pSrc)
            {
                return (m_id < pSrc.m_id);
            }
            /// Define node equality based on coordinate.
            bool operator==(const Node& pSrc)
            {
                return m_x == pSrc.m_x && m_y == pSrc.m_y && m_z == pSrc.m_z;
            }

            Node operator+(const Node &pSrc) const
            {
                return Node(m_id, m_x+pSrc.m_x, m_y+pSrc.m_y, m_z+pSrc.m_z);
            }

            Node operator-(const Node &pSrc) const
            {
                return Node(m_id, m_x-pSrc.m_x, m_y-pSrc.m_y, m_z-pSrc.m_z);
            }

            Node operator*(const Node &pSrc) const
            {
                return Node(m_id, m_x*pSrc.m_x, m_y*pSrc.m_y, m_z*pSrc.m_z);
            }

            Node operator*(const NekDouble &alpha) const
            {
                return Node(m_id, alpha*m_x, alpha*m_y, alpha*m_z);
            }

            Node operator/(const NekDouble &alpha) const
            {
                return Node(m_id, m_x/alpha, m_y/alpha, m_z/alpha);
            }

            void operator+=(const Node &pSrc)
            {
                m_x += pSrc.m_x;
                m_y += pSrc.m_y;
                m_z += pSrc.m_z;
            }

            void operator*=(const NekDouble &alpha)
            {
                m_x *= alpha;
                m_y *= alpha;
                m_z *= alpha;
            }

            void operator/=(const NekDouble &alpha)
            {
                m_x /= alpha;
                m_y /= alpha;
                m_z /= alpha;
            }

            NekDouble abs2() const
            {
                return m_x*m_x+m_y*m_y+m_z*m_z;
            }

            NekDouble dot(const Node &pSrc) const
            {
                return m_x*pSrc.m_x + m_y*pSrc.m_y + m_z*pSrc.m_z;
            }


            Node curl(const Node &pSrc) const
            {
                return Node(m_id, m_y*pSrc.m_z - m_z*pSrc.m_y,
                            m_z*pSrc.m_x-m_x*pSrc.m_z, m_x*pSrc.m_y-m_y*pSrc.m_x);
            }

            /// Generate a %SpatialDomains::PointGeom for this node.
            SpatialDomains::PointGeomSharedPtr GetGeom(int coordDim)
            {
                SpatialDomains::PointGeomSharedPtr ret =
                    MemoryManager<SpatialDomains::PointGeom>
                        ::AllocateSharedPtr(coordDim,m_id,m_x,m_y,m_z);

                return ret;
            }

            /// ID of node.
            int m_id;
            /// X-coordinate.
            NekDouble m_x;
            /// Y-coordinate.
            NekDouble m_y;
            /// Z-coordinate.
            NekDouble m_z;

        private:
            SpatialDomains::PointGeomSharedPtr m_geom;
        };
        /// Shared pointer to a Node.
        typedef boost::shared_ptr<Node> NodeSharedPtr;

        bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2);
        bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2);
        std::ostream &operator<<(std::ostream &os, const NodeSharedPtr &n);

        /**
         * @brief Defines a hash function for nodes.
         *
         * The hash of a node is straight-forward; a combination of the x, y,
         * and z co-ordinates in this order.
         */
        struct NodeHash : std::unary_function<NodeSharedPtr, std::size_t>
        {
            std::size_t operator()(NodeSharedPtr const& p) const
            {
                std::size_t seed = 0;
                boost::hash_combine(seed, p->m_x);
                boost::hash_combine(seed, p->m_y);
                boost::hash_combine(seed, p->m_z);
                return seed;
            }
        };
        typedef boost::unordered_set<NodeSharedPtr, NodeHash> NodeSet;


        /**
         * @brief Represents an edge which joins two points.
         *
         * An edge is defined by two nodes (vertices) and, for high-order edges,
         * a set of control nodes defining the shape of the edge.
         */
        class Edge {
        public:
            /// Creates a new edge.
            Edge(NodeSharedPtr pVertex1, NodeSharedPtr pVertex2,
                 std::vector<NodeSharedPtr> pEdgeNodes,
                 LibUtilities::PointsType pCurveType)
                : m_n1(pVertex1), m_n2(pVertex2), m_edgeNodes(pEdgeNodes),
                  m_curveType(pCurveType), m_geom() {}
            /// Copies an existing edge.
            Edge(const Edge& pSrc)
                : m_n1(pSrc.m_n1), m_n2(pSrc.m_n2), m_edgeNodes(pSrc.m_edgeNodes),
                  m_curveType(pSrc.m_curveType), m_geom(pSrc.m_geom) {}
            ~Edge() {}

            /// Returns the total number of nodes defining the edge.
            unsigned int GetNodeCount() const
            {
                return m_edgeNodes.size() + 2;
            }
            void GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
            {
                nodeList.push_back(m_n1);
                for (int k = 0; k < m_edgeNodes.size(); ++k) 
                {
                    nodeList.push_back(m_edgeNodes[k]);
                }
                nodeList.push_back(m_n2);
            }

            /// Creates a Nektar++ string listing the coordinates of all the
            /// nodes.
            std::string GetXmlCurveString() const
            {
                std::vector<NodeSharedPtr> nodeList;

                GetCurvedNodes(nodeList);

                std::stringstream s;
                std::string str;

                // put them into a string and return
                for (int k = 0; k < nodeList.size(); ++k)
                {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  nodeList[k]->m_x << "  " << nodeList[k]->m_y
                      << "  " << nodeList[k]->m_z << "    ";

                }

                return s.str();
            }

            /// Generate a SpatialDomains::SegGeom object for this edge.
            SpatialDomains::SegGeomSharedPtr GetGeom(int coordDim)
            {
                // Create edge vertices.
                SpatialDomains::PointGeomSharedPtr p[2];
                SpatialDomains::SegGeomSharedPtr ret;

                p[0] = m_n1->GetGeom(coordDim);
                p[1] = m_n2->GetGeom(coordDim);

                // Create a curve if high-order information exists.
                if (m_edgeNodes.size() > 0)
                {
                    SpatialDomains::CurveSharedPtr c =
                        MemoryManager<SpatialDomains::Curve>::
                        AllocateSharedPtr(m_id, m_curveType);

                    c->m_points.push_back(p[0]);
                    for (int i = 0; i < m_edgeNodes.size(); ++i)
                    {
                        c->m_points.push_back(m_edgeNodes[i]->GetGeom(coordDim));
                    }
                    c->m_points.push_back(p[1]);

                    ret = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(m_id, coordDim, p, c);
                }
                else
                {
                    ret = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(m_id, coordDim, p);
                }

                return ret;
            }

            /// ID of edge.
            unsigned int m_id;
            /// First vertex node.
            NodeSharedPtr m_n1;
            /// Second vertex node.
            NodeSharedPtr m_n2;
            /// List of control nodes between the first and second vertices.
            std::vector<NodeSharedPtr> m_edgeNodes;
            /// Distributions of points along edge.
            LibUtilities::PointsType m_curveType;
            /// Element(s) which are linked to this edge.
            vector<pair<ElementSharedPtr, int> > m_elLink;

        private:
            SpatialDomains::SegGeomSharedPtr m_geom;
        };
        /// Shared pointer to an edge.
        typedef boost::shared_ptr<Edge> EdgeSharedPtr;

        bool operator==(EdgeSharedPtr const &p1, EdgeSharedPtr const &p2);
        bool operator< (EdgeSharedPtr const &p1, EdgeSharedPtr const &p2);

        /**
         * @brief Defines a hash function for edges.
         *
         * The hash of an edge is defined using the IDs of the two nodes which
         * define it. First the minimum ID is hashed, then the maximum
         * ID, which takes the two possible orientations into account.
         */
        struct EdgeHash : std::unary_function<EdgeSharedPtr, std::size_t>
        {
            std::size_t operator()(EdgeSharedPtr const& p) const
            {
                std::size_t seed = 0;
                unsigned int id1 = p->m_n1->m_id;
                unsigned int id2 = p->m_n2->m_id;
                boost::hash_combine(seed, id1 < id2 ? id1 : id2);
                boost::hash_combine(seed, id2 < id1 ? id1 : id2);
                return seed;
            }
        };
        typedef boost::unordered_set<EdgeSharedPtr, EdgeHash> EdgeSet;


        /**
         * @brief Represents a face comprised of three or more edges.
         *
         * A face is defined by a list of vertices, a list of edges joining
         * these vertices, and a list of control nodes within the interior of
         * the face, defining the shape of the face.
         */
        class Face {
        public:
            /// Create a new face.
            Face(std::vector<NodeSharedPtr> pVertexList,
                 std::vector<NodeSharedPtr> pFaceNodes,
                 std::vector<EdgeSharedPtr> pEdgeList,
                 LibUtilities::PointsType   pCurveType)
                : m_vertexList(pVertexList),
                m_edgeList  (pEdgeList),
                m_faceNodes (pFaceNodes),
                m_curveType (pCurveType),
                 m_geom    () {}

            /// Copy an existing face.
            Face(const Face& pSrc)
                : m_vertexList(pSrc.m_vertexList), m_edgeList  (pSrc.m_edgeList),
                m_faceNodes (pSrc.m_faceNodes),  m_curveType (pSrc.m_curveType),
                  m_geom    (pSrc.m_geom) {}
            ~Face() {}

            /// Equality is defined by matching all vertices.
            bool operator==(Face& pSrc)
            {
                std::vector<NodeSharedPtr>::iterator it1, it2;
                for (it1 = m_vertexList.begin(); it1 != m_vertexList.end(); ++it1)
                {
                    if (find(pSrc.m_vertexList.begin(), pSrc.m_vertexList.end(), *it1)
                            == pSrc.m_vertexList.end())
                    {
                        return false;
                    }
                }
                return true;
            }

            /// Returns the total number of nodes (vertices, edge nodes and
            /// face nodes).
            unsigned int GetNodeCount() const
            {
                unsigned int n = m_faceNodes.size();
                for (int i = 0; i < m_edgeList.size(); ++i)
                {
                    n += m_edgeList[i]->GetNodeCount();
                }
                n -= m_vertexList.size();
                return n;
            }

            /// Assemble a list of nodes on curved face
            void GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
            {
                // Treat 2D point distributions differently to 3D.
                if (m_curveType == LibUtilities::eNodalTriFekete       ||
                    m_curveType == LibUtilities::eNodalTriEvenlySpaced ||
                    m_curveType == LibUtilities::eNodalTriElec)
                {
                    int n = m_edgeList[0]->GetNodeCount();
                    
                    nodeList.insert(nodeList.end(), m_vertexList.begin(), m_vertexList.end());
                    for (int k = 0; k < m_edgeList.size(); ++k) 
                    {
                        nodeList.insert(nodeList.end(), m_edgeList[k]->m_edgeNodes.begin(),
                                   m_edgeList[k]->m_edgeNodes.end());
                        if (m_edgeList[k]->m_n1 != m_vertexList[k])
                        {
                            // If edge orientation is reversed relative to node
                            // ordering, we need to reverse order of nodes.
                            std::reverse(nodeList.begin() + 3 + k*(n-2),
                                         nodeList.begin() + 3 + (k+1)*(n-2));
                        }
                    }
                    nodeList.insert(nodeList.end(), m_faceNodes.begin(), m_faceNodes.end());
                }
                else
                {
                    // Write out in 2D tensor product order.
                    ASSERTL0(m_vertexList.size() == 4,
                             "Face nodes of tensor product only supported "
                             "for quadrilaterals.");

                    int n = (int)sqrt((NekDouble)GetNodeCount());
                    nodeList.resize(n*n);

                    ASSERTL0(n*n == GetNodeCount(), "Wrong number of modes?");

                    // Write vertices
                    nodeList[0]       = m_vertexList[0];
                    nodeList[n-1]     = m_vertexList[1];
                    nodeList[n*n-1]   = m_vertexList[2];
                    nodeList[n*(n-1)] = m_vertexList[3];
                    
                    // Write edge-interior
                    int skips[4][2] = {{0,1}, {n-1,n}, {n*n-1,-1}, {n*(n-1),-n}};
                    for (int i = 0; i < 4; ++i)
                    {
                        bool reverseEdge = m_edgeList[i]->m_n1 == m_vertexList[i];

                        if (!reverseEdge)
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                nodeList[skips[i][0] + j*skips[i][1]] = 
                                    m_edgeList[i]->m_edgeNodes[n-2-j];
                            }
                        }
                        else
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                nodeList[skips[i][0] + j*skips[i][1]] = 
                                    m_edgeList[i]->m_edgeNodes[j-1];
                            }
                        }
                    }
                    
                    // Write interior
                    for (int i = 1; i < n-1; ++i)
                    {
                        for (int j = 1; j < n-1; ++j)
                        {
                            nodeList[i*n+j] = m_faceNodes[(i-1)*(n-2)+(j-1)];
                        }
                    }
                }
            }

            /// Generates a string listing the coordinates of all nodes
            /// associated with this face.
            std::string GetXmlCurveString() const
            {
                std::stringstream s;
                std::string str;
                vector<NodeSharedPtr> nodeList;
                
                // assemble listof nodes
                GetCurvedNodes(nodeList);
                
                // put them into a string 
                for (int k = 0; k < nodeList.size(); ++k) 
                {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  nodeList[k]->m_x << "  " << nodeList[k]->m_y
                      << "  " << nodeList[k]->m_z << "    ";
                }
                
                return s.str();
            }

            /// Generate either SpatialDomains::TriGeom or
            /// SpatialDomains::QuadGeom for this element.
            SpatialDomains::Geometry2DSharedPtr GetGeom(int coordDim)
            {
                int nEdge = m_edgeList.size();

                SpatialDomains::SegGeomSharedPtr    edges[4];
                SpatialDomains::Geometry2DSharedPtr ret;
                StdRegions::Orientation             edgeo[4];

                for (int i = 0; i < nEdge; ++i)
                {
                    edges[i] = m_edgeList[i]->GetGeom(coordDim);
                }

                for (int i = 0; i < nEdge; ++i)
                {
                    edgeo[i] = SpatialDomains::SegGeom::GetEdgeOrientation(
                                            *edges[i], *edges[(i+1) % nEdge]);

                }

                if(m_faceNodes.size() > 0)
                {
                    if (nEdge == 3)
                    {
                        SpatialDomains::CurveSharedPtr c =
                            MemoryManager<SpatialDomains::Curve>::
                            AllocateSharedPtr(m_id, m_curveType);

                        for(int j = 0; j < m_vertexList.size(); j++)
                        {
                            c->m_points.push_back(m_vertexList[j]->GetGeom(coordDim));
                        }
                        for(int j = 0; j < nEdge; j++)
                        {
                            vector<NodeSharedPtr> ed = m_edgeList[j]->m_edgeNodes;
                            if(edgeo[j] == StdRegions::eBackwards)
                            {
                                for(int k = ed.size()-1; k >=0; k--)
                                {
                                    c->m_points.push_back(ed[k]->GetGeom(coordDim));
                                }
                            }
                            else
                            {
                                for(int k = 0; k < ed.size(); k++)
                                {
                                    c->m_points.push_back(ed[k]->GetGeom(coordDim));
                                }
                            }
                        }
                        for(int j = 0; j < m_faceNodes.size(); j++)
                        {
                            c->m_points.push_back(m_faceNodes[j]->GetGeom(coordDim));
                        }

                        ret = MemoryManager<SpatialDomains::TriGeom>::
                            AllocateSharedPtr(m_id, edges, edgeo, c);
                    }
                    else
                    {
                        SpatialDomains::CurveSharedPtr c =
                            MemoryManager<SpatialDomains::Curve>::
                            AllocateSharedPtr(m_id, m_curveType);

                        ASSERTL0(m_vertexList.size() == 4,
                                 "Face nodes of tensor product only supported "
                                 "for quadrilaterals.");

                        int n = (int)sqrt((NekDouble)GetNodeCount());
                        vector<NodeSharedPtr> tmp(n*n);

                        ASSERTL0(n*n == GetNodeCount(), "Wrong number of modes?");

                        // Write vertices
                        tmp[0]       = m_vertexList[0];
                        tmp[n-1]     = m_vertexList[1];
                        tmp[n*n-1]   = m_vertexList[2];
                        tmp[n*(n-1)] = m_vertexList[3];

                        // Write edge-interior
                        int skips[4][2] = {{0,1}, {n-1,n}, {n*n-1,-1}, {n*(n-1),-n}};
                        for (int i = 0; i < 4; ++i)
                        {
                            bool reverseEdge = edgeo[i] == StdRegions::eBackwards;

                            if (reverseEdge)
                            {
                                for (int j = 1; j < n-1; ++j)
                                {
                                    tmp[skips[i][0] + j*skips[i][1]] =
                                        m_edgeList[i]->m_edgeNodes[n-2-j];
                                }
                            }
                            else
                            {
                                for (int j = 1; j < n-1; ++j)
                                {
                                    tmp[skips[i][0] + j*skips[i][1]] =
                                        m_edgeList[i]->m_edgeNodes[j-1];
                                }
                            }
                        }

                        // Write interior
                        for (int i = 1; i < n-1; ++i)
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                tmp[i*n+j] = m_faceNodes[(i-1)*(n-2)+(j-1)];
                            }
                        }

                        for (int k = 0; k < tmp.size(); ++k)
                        {
                            c->m_points.push_back(tmp[k]->GetGeom(coordDim));
                        }

                        ret = MemoryManager<SpatialDomains::QuadGeom>::
                            AllocateSharedPtr(m_id, edges, edgeo, c);
                    }
                }
                else
                {
                    if (nEdge == 3)
                    {
                        ret = MemoryManager<SpatialDomains::TriGeom>::
                            AllocateSharedPtr(m_id, edges, edgeo);
                    }
                    else
                    {
                        ret = MemoryManager<SpatialDomains::QuadGeom>::
                            AllocateSharedPtr(m_id, edges, edgeo);
                    }
                }

                return ret;
            }

            /// ID of the face.
            unsigned int m_id;
            /// List of vertex nodes.
            std::vector<NodeSharedPtr> m_vertexList;
            /// List of corresponding edges.
            std::vector<EdgeSharedPtr> m_edgeList;
            /// List of face-interior nodes defining the shape of the face.
            std::vector<NodeSharedPtr> m_faceNodes;
            /// Distribution of points in this face.
            LibUtilities::PointsType   m_curveType;
            /// Element(s) which are linked to this face.
            vector<pair<ElementSharedPtr, int> > m_elLink;

            SpatialDomains::Geometry2DSharedPtr m_geom;
        };
        /// Shared pointer to a face.
        typedef boost::shared_ptr<Face> FaceSharedPtr;

        bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2);
        bool operator< (FaceSharedPtr const &p1, FaceSharedPtr const &p2);

        struct FaceHash : std::unary_function<FaceSharedPtr, std::size_t>
        {
            std::size_t operator()(FaceSharedPtr const& p) const
            {
                unsigned int              nVert = p->m_vertexList.size();
                std::size_t               seed  = 0;
                std::vector<unsigned int> ids(nVert);

                for (int i = 0; i < nVert; ++i)
                {
                    ids[i] = p->m_vertexList[i]->m_id;
                }

                std::sort(ids.begin(), ids.end());
                boost::hash_range(seed, ids.begin(), ids.end());

                return seed;
            }
        };
        typedef boost::unordered_set<FaceSharedPtr, FaceHash> FaceSet;


        /**
         * @brief Basic information about an element.
         *
         * ElmtConfig contains four member variables which denote the
         * properties of an element when it is created.
         */
        struct ElmtConfig
        {
        ElmtConfig(LibUtilities::ShapeType  pE,
                   unsigned int             pOrder,
                   bool                     pFn,
                   bool                     pVn,
                   bool                     pReorient = true,
                   LibUtilities::PointsType pECt
                   = LibUtilities::ePolyEvenlySpaced,
                   LibUtilities::PointsType pFCt
                   = LibUtilities::ePolyEvenlySpaced)
        : m_e            (pE),
                m_faceNodes    (pFn),
                m_volumeNodes  (pVn),
                m_order        (pOrder),
                m_reorient     (pReorient),
                m_edgeCurveType(pECt),
                m_faceCurveType(pFCt)
            {
            }

        ElmtConfig(ElmtConfig const &p) :
            m_e            (p.m_e),
                m_faceNodes    (p.m_faceNodes),
                m_volumeNodes  (p.m_volumeNodes),
                m_order        (p.m_order),
                m_reorient     (p.m_reorient),
                m_edgeCurveType(p.m_edgeCurveType),
                m_faceCurveType(p.m_faceCurveType)
            {
            }

            ElmtConfig() {}

            /// Element type (e.g. triangle, quad, etc).
            LibUtilities::ShapeType   m_e;
            /// Denotes whether the element contains face nodes. For 2D
            /// elements, if this is true then the element contains interior
            /// nodes.
            bool                      m_faceNodes;
            /// Denotes whether the element contains volume (i.e. interior)
            /// nodes. These are not supported by either the mesh converter or
            /// Nektar++ but are included for completeness and are required
            /// for some output modules (e.g. Gmsh).
            bool                     m_volumeNodes;
            /// Order of the element.
            unsigned int             m_order;
            /// Denotes whether the element needs to be re-orientated for a
            /// spectral element framework.
            bool                     m_reorient;
            /// Distribution of points in edges.
            LibUtilities::PointsType m_edgeCurveType;
            /// Distribution of points in faces.
            LibUtilities::PointsType m_faceCurveType;
        };


        /**
         * @brief Base class for element definitions.
         *
         * An element is defined by a list of vertices, edges and faces
         * (depending on the dimension of the problem). This base class
         * provides the underlying structure.
         */
        class Element {
        public:
            Element(ElmtConfig   pConf,
                    unsigned int pNumNodes,
                    unsigned int pGotNodes);

            /// Returns the ID of the element (or associated edge or face for
            /// boundary elements).
            unsigned int GetId() const {
                if (m_faceLink.get() != 0) return m_faceLink->m_id;
                if (m_edgeLink.get() != 0) return m_edgeLink->m_id;
                return m_id;
            }
            /// Returns the expansion dimension of the element.
            unsigned int GetDim() const {
                return m_dim;
            }
            /// Returns the configuration of the element.
            ElmtConfig GetConf() const {
                return m_conf;
            }
            /// Returns the tag which defines the element shape.
            std::string GetTag() const {
                if (m_faceLink.get() != 0) return "F";
                if (m_edgeLink.get() != 0) return "E";
                return m_tag;
            }
            /// Access a vertex node.
            NodeSharedPtr GetVertex(unsigned int i) const {
                return m_vertex[i];
            }
            /// Access an edge.
            EdgeSharedPtr GetEdge(unsigned int i) const {
                return m_edge[i];
            }
            /// Access a face.
            FaceSharedPtr GetFace(unsigned int i) const {
                return m_face[i];
            }
            /// Access the list of vertex nodes.
            std::vector<NodeSharedPtr> GetVertexList() const {
                return m_vertex;
            }
            /// Access the list of edges.
            std::vector<EdgeSharedPtr> GetEdgeList() const {
                return m_edge;
            }
            /// Access the list of faces.
            std::vector<FaceSharedPtr> GetFaceList() const {
                return m_face;
            }
            /// Access the list of volume nodes.
            std::vector<NodeSharedPtr> GetVolumeNodes() const {
                return m_volumeNodes;
            }
            void SetVolumeNodes(std::vector<NodeSharedPtr> &nodes) {
                m_volumeNodes = nodes;
            }
            LibUtilities::PointsType GetCurveType() const {
                return m_curveType;
            }
            void SetCurveType(LibUtilities::PointsType cT) {
                m_curveType = cT;
            }
            /// Returns the total number of nodes (vertices, edge nodes and
            /// face nodes and volume nodes).
            unsigned int GetNodeCount() const
            {
                unsigned int n = m_volumeNodes.size();
                if (m_dim == 1)
                {
                    n += 2;
                }
                else if (m_dim == 2)
                {
                    for (int i = 0; i < m_edge.size(); ++i)
                    {
                        n += m_edge[i]->GetNodeCount();
                    }
                    n -= m_vertex.size();
                }
                else
                {
                    cerr << "Not supported." << endl;
                    exit(1);
                }
                return n;
            }
            /// Access the list of tags associated with this element.
            std::vector<int> GetTagList() const {
                return m_taglist;
            }
            /// Returns the number of vertices.
            unsigned int GetVertexCount() const {
                return m_vertex.size();
            }
            /// Returns the number of edges.
            unsigned int GetEdgeCount() const {
                return m_edge.size();
            }
            /// Returns the number of faces.
            unsigned int GetFaceCount() const {
                return m_face.size();
            }
            /// Change the ID of the element.
            void SetId(unsigned int p) {
                m_id = p;
            }
            /// Replace a vertex with another vertex object.
            void SetVertex(unsigned int p, NodeSharedPtr pNew);
            /// Replace an edge with another edge object.
            void SetEdge(unsigned int p, EdgeSharedPtr pNew);
            /// Replace a face with another face object.
            void SetFace(unsigned int p, FaceSharedPtr pNew);
            /// Set a correspondence between this element and an edge
            /// (2D boundary element).
            void SetEdgeLink(EdgeSharedPtr pLink) {
                m_edgeLink = pLink;
            }
            /// Get correspondence between this element and an edge.
            EdgeSharedPtr GetEdgeLink() {
                return m_edgeLink;
            }
            /// Set a correspondence between this element and a face
            /// (3D boundary element).
            void SetFaceLink(FaceSharedPtr pLink) {
                m_faceLink = pLink;
            }
            /// Get correspondence between this element and a face.
            FaceSharedPtr GetFaceLink() {
                return m_faceLink;
            }
            /// Set a correspondence between edge or face i and its
            /// representative boundary element m->element[expDim-1][j].
            void SetBoundaryLink(int i, int j) {
                m_boundaryLinks[i] = j;
            }
            /// Get the location of the boundary face/edge i for this element.
            int GetBoundaryLink(int i) {
                std::map<int,int>::iterator it = m_boundaryLinks.find(i);
                if (it == m_boundaryLinks.end())
                {
                    return -1;
                }
                else
                {
                    return it->second;
                }
            }
            /// Set the list of tags associated with this element.
            void SetTagList(const std::vector<int> &tags)
            {
                m_taglist = tags;
            }
            /// Generate a list of vertices (1D), edges (2D), or faces (3D).
            virtual std::string GetXmlString() const
            {
                std::stringstream s;
                switch (m_dim)
                {
                case 1:
                    for(int j=0; j< m_vertex.size(); ++j)
                    {
                        s << std::setw(5) << m_vertex[j]->m_id << " ";
                    }
                    break;
                case 2:
                    for(int j=0; j< m_edge.size(); ++j)
                    {
                        s << std::setw(5) << m_edge[j]->m_id << " ";
                    }
                    break;
                case 3:
                    for(int j=0; j< m_face.size(); ++j)
                    {
                        s << std::setw(5) << m_face[j]->m_id << " ";
                    }
                    break;
                }
                return s.str();
            }

            void GetCurvedNodes(std::vector<NodeSharedPtr> &nodeList) const
            {
                // Node orderings are different for different elements.
                // Triangle
                if (m_vertex.size() == 2)
                {
                    nodeList.push_back(m_vertex[0]);
                    for (int i = 0; i < m_volumeNodes.size(); ++i)
                    {
                        nodeList.push_back(m_volumeNodes[i]);
                    }
                    nodeList.push_back(m_vertex[1]);
                }
                else if (m_vertex.size() == 3)
                {
                    int n = m_edge[0]->GetNodeCount();
                    nodeList.resize(n*(n+1)/2);

                    // Populate nodelist
                    std::copy(m_vertex.begin(), m_vertex.end(), nodeList.begin());
                    for (int i = 0; i < 3; ++i)
                    {
                        std::copy(m_edge[i]->m_edgeNodes.begin(),
                                  m_edge[i]->m_edgeNodes.end(),
                                  nodeList.begin() + 3 + i*(n-2));
                        if (m_edge[i]->m_n1 != m_vertex[i])
                        {
                            // If edge orientation is reversed relative to node
                            // ordering, we need to reverse order of nodes.
                            std::reverse(nodeList.begin() + 3 + i*(n-2),
                                         nodeList.begin() + 3 + (i+1)*(n-2));
                        }
                    }

                    // Copy volume nodes.
                    std::copy(m_volumeNodes.begin(), m_volumeNodes.end(),
                              nodeList.begin() + 3*(n-1));
                }
                // Quad
                else if (m_dim == 2 && m_vertex.size() == 4)
                {
                    int n = m_edge[0]->GetNodeCount();
                    nodeList.resize(n*n);

                    // Write vertices
                    nodeList[0]       = m_vertex[0];
                    nodeList[n-1]     = m_vertex[1];
                    nodeList[n*n-1]   = m_vertex[2];
                    nodeList[n*(n-1)] = m_vertex[3];

                    // Write edge-interior
                    int skips[4][2] = {{0,1}, {n-1,n}, {n*n-1,-1}, {n*(n-1),-n}};
                    for (int i = 0; i < 4; ++i)
                    {
                        bool reverseEdge = m_edge[i]->m_n1 == m_vertex[i];

                        if (!reverseEdge)
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                nodeList[skips[i][0] + j*skips[i][1]] =
                                    m_edge[i]->m_edgeNodes[n-2-j];
                            }
                        }
                        else
                        {
                            for (int j = 1; j < n-1; ++j)
                            {
                                nodeList[skips[i][0] + j*skips[i][1]] =
                                    m_edge[i]->m_edgeNodes[j-1];
                            }
                        }
                    }

                    // Write interior
                    for (int i = 1; i < n-1; ++i)
                    {
                        for (int j = 1; j < n-1; ++j)
                        {
                            nodeList[i*n+j] = m_volumeNodes[(i-1)*(n-2)+(j-1)];
                        }
                    }
                }
                else
                {
                    cerr << "GetXmlCurveString for a " << m_vertex.size()
                         << "-vertex element is not yet implemented." << endl;
                }
            }

            /// Generates a string listing the coordinates of all nodes
            /// associated with this element.
            std::string GetXmlCurveString() const
            {
                // Temporary node list for reordering
                std::vector<NodeSharedPtr> nodeList;

                GetCurvedNodes(nodeList);

                // Finally generate the XML string corresponding to our new
                // node reordering.
                std::stringstream s;
                std::string str;
                for (int k = 0; k < nodeList.size(); ++k)
                {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  nodeList[k]->m_x << "  " << nodeList[k]->m_y
                      << "  " << nodeList[k]->m_z << "    ";

                }
                return s.str();
            }

            /// Generate a Nektar++ geometry object for this element.
            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim)
            {
                ASSERTL0(false, "This function should be implemented at a shape level.");
                return boost::shared_ptr<SpatialDomains::Geometry>();
            }
            int GetMaxOrder();
            /// Complete this object.
            virtual void Complete(int order)
            {
                ASSERTL0(false, "This function should be implemented at a shape level.");
            }
            void Print()
            {
                int i, j;
                for (i = 0; i < m_vertex.size(); ++i)
                {
                    cout << m_vertex[i]->m_x << " " << m_vertex[i]->m_y << " " << m_vertex[i]->m_z << endl;
                }
                for (i = 0; i < m_edge.size(); ++i)
                {
                    for (j = 0; j < m_edge[i]->m_edgeNodes.size(); ++j)
                    {
                        NodeSharedPtr n = m_edge[i]->m_edgeNodes[j];
                        cout << n->m_x << " " << n->m_y << " " << n->m_z << endl;
                    }
                }
                for (i = 0; i < m_face.size(); ++i)
                {
                    for (j = 0; j < m_face[i]->m_faceNodes.size(); ++j)
                    {
                        NodeSharedPtr n = m_face[i]->m_faceNodes[j];
                        cout << n->m_x << " " << n->m_y << " " << n->m_z << endl;
                    }
                }
            }

        protected:
            /// ID of the element.
            unsigned int m_id;
            /// Dimension of the element.
            unsigned int m_dim;
            /// Contains configuration of the element.
            ElmtConfig   m_conf;
            /// Tag character describing the element.
            std::string m_tag;
            /// List of integers specifying properties of the element.
            std::vector<int> m_taglist;
            /// List of element vertex nodes.
            std::vector<NodeSharedPtr> m_vertex;
            /// List of element edges.
            std::vector<EdgeSharedPtr> m_edge;
            /// List of element faces.
            std::vector<FaceSharedPtr> m_face;
            /// List of element volume nodes.
            std::vector<NodeSharedPtr> m_volumeNodes;
            /// Volume curve type
            LibUtilities::PointsType m_curveType;
            /// Pointer to the corresponding edge if element is a 2D boundary.
            EdgeSharedPtr m_edgeLink;
            /// Pointer to the corresponding face if element is a 3D boundary.
            FaceSharedPtr m_faceLink;
            /// Array mapping faces/edges to the location of the appropriate
            /// boundary elements in m->element.
            std::map<int,int> m_boundaryLinks;
            /// Nektar++ geometry object for this element.
            SpatialDomains::GeometrySharedPtr m_geom;
        };
        /// Container for elements; key is expansion dimension, value is
        /// vector of elements of that dimension.
        typedef std::map<unsigned int, std::vector<ElementSharedPtr> > ElementMap;
        /// Element factory definition.
        typedef Nektar::LibUtilities::NekFactory<LibUtilities::ShapeType, Element,
            ElmtConfig, std::vector<NodeSharedPtr>, std::vector<int> > ElementFactory;
        ElementFactory& GetElementFactory();

        /// Define element ordering based on ID.
        struct element_id_less_than
        {
            typedef boost::shared_ptr<Element> pT;
            bool operator()(const pT a, const pT b) const
            {
                // check for 0
                if (a.get() == 0)
                {
                    // if b is also 0, then they are equal, hence a is not
                    // less than b
                    return b.get() != 0;
                }
                else if (b.get() == 0)
                {
                    return false;
                }
                else
                {
                    return a->GetId() < b->GetId();
                }
            }
        };


        /**
         * @brief A composite is a collection of elements.
         *
         * All elements should be of the same type, i.e. have the same tag.
         */
        class Composite {
        public:
            Composite() : m_reorder(true) {}

            /// Generate the list of IDs of elements within this composite.
            std::string GetXmlString(bool doSort=true);

            /// ID of composite.
            unsigned int m_id;
            /// Element type tag.
            std::string m_tag;
            /// boundary label
            std::string m_label;
            /// Determines whether items can be reordered.
            bool m_reorder;
            /// List of elements in this composite.
            std::vector<ElementSharedPtr> m_items;
        };
        /// Shared pointer to a composite.
        typedef boost::shared_ptr<Composite> CompositeSharedPtr;
        /// Container of composites; key is the composite id, value is the
        /// composite.
        typedef std::map<unsigned int, CompositeSharedPtr> CompositeMap;

        /**
         * Enumeration of condition types (Dirichlet, Neumann, etc).
         */
        enum ConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic,
            eHOPCondition,
            SIZE_ConditionType
        };

        /**
         * @brief Defines a boundary condition.
         *
         * A boundary condition is defined by its type (e.g. Dirichlet), the
         * field it applies to, the value imposed on this field and the
         * composite which the boundary condition is applied to.
         */
        struct Condition
        {
        Condition() : type(), field(), value(), m_composite() {}
            std::vector<ConditionType> type;
            std::vector<std::string>   field;
            std::vector<std::string>   value;
            std::vector<int>           m_composite;
        };

        typedef boost::shared_ptr<Condition> ConditionSharedPtr;
        typedef std::map<int,ConditionSharedPtr> ConditionMap;

        bool operator==(ConditionSharedPtr const &c1, ConditionSharedPtr const &c2);

        class Mesh
        {
        public:
            Mesh() : m_verbose(false) {}

            /// Verbose flag
            bool                            m_verbose;
            /// Dimension of the expansion.
            unsigned int                    m_expDim;
            /// Dimension of the space in which the mesh is defined.
            unsigned int                    m_spaceDim;
            /// List of mesh nodes.
            std::vector<NodeSharedPtr>      m_node;
            /// Set of element vertices.
            NodeSet                         m_vertexSet;
            /// Set of element edges.
            EdgeSet                         m_edgeSet;
            /// Set of element faces.
            FaceSet                         m_faceSet;
            /// Map for elements.
            ElementMap                      m_element;
            /// Map for composites.
            CompositeMap                    m_composite;
            /// Boundary conditions maps tag to condition.
            ConditionMap                    m_condition;
            /// List of fields names.
            std::vector<std::string>        m_fields;
            /// Map of vertex normals.
            boost::unordered_map<int, Node> m_vertexNormals;
            /// Set of all pairs of element ID and edge/face number on which to
            /// apply spherigon surface smoothing.
            set<pair<int,int> >             m_spherigonSurfs;
            /// List of face labels for composite annotation
            map<int,string>                 m_faceLabels;

            /// Returns the total number of elements in the mesh with
            /// dimension expDim.
            unsigned int                    GetNumElements();
            /// Returns the total number of elements in the mesh with
            /// dimension < expDim.
            unsigned int                    GetNumBndryElements();
            /// Returns the total number of entities in the mesh.
            unsigned int                    GetNumEntities();

        };
        /// Shared pointer to a mesh.
        typedef boost::shared_ptr<Mesh> MeshSharedPtr;


        /**
         * @brief A 0-dimensional vertex.
         */
        class Point : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                return boost::shared_ptr<Element>(
                    new Point(pConf, pNodeList, pTagList));
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Point(ElmtConfig                 pConf,
                  std::vector<NodeSharedPtr> pNodeList,
                  std::vector<int>           pTagList);
            Point(const Point& pSrc);
            virtual ~Point() {}

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };


        /**
         * @brief A 1-dimensional line between two vertex nodes.
         */
        class Line : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                return boost::shared_ptr<Element>(
                    new Line(pConf, pNodeList, pTagList));
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Line(ElmtConfig                 pConf,
                 std::vector<NodeSharedPtr> pNodeList,
                 std::vector<int>           pTagList);
            Line(const Point& pSrc);
            virtual ~Line() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };

        /**
         * @brief A lightweight struct for dealing with high-order triangle
         * alignment.
         *
         * The logic underlying these routines is taken from the original Nektar
         * code.
         */
        template<typename T>
        struct HOTriangle
        {
            HOTriangle(vector<int> pVertId,
                       vector<T>   pSurfVerts) :
                vertId(pVertId), surfVerts(pSurfVerts) {}
            HOTriangle(vector<int> pVertId) : vertId(pVertId) {}

            /// The triangle vertex IDs
            vector<int> vertId;

            /// The triangle surface vertices -- templated so that this can
            /// either be nodes or IDs.
            vector<T> surfVerts;

            /**
             * @brief Rotates the triangle of data points inside #surfVerts
             * counter-clockwise nrot times.
             *
             * @param nrot Number of times to rotate triangle.
             */
            void Rotate(int nrot)
            {
                int n, i, j, cnt;
                int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
                vector<T> tmp(np*np);

                for (n = 0; n < nrot; ++n)
                {
                    for (cnt = i = 0; i < np; ++i)
                    {
                        for (j = 0; j < np-i; ++j, cnt++)
                        {
                            tmp[i*np+j] = surfVerts[cnt];
                        }
                    }
                    for (cnt = i = 0; i < np; ++i)
                    {
                        for (j = 0; j < np-i; ++j,cnt++)
                        {
                            surfVerts[cnt] = tmp[(np-1-i-j)*np+i];
                        }
                    }
                }
            }

            /**
             * @brief Reflect data points inside #surfVerts.
             *
             * This applies a mapping essentially doing the following
             * reordering:
             *
             * 9          9
             * 7 8    ->  8 7
             * 4 5 6      6 5 4
             * 0 1 2 3    3 2 1 0
             */
            void Reflect()
            {
                int i, j, cnt;
                int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
                vector<T> tmp(np*np);

                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j,cnt++)
                    {
                        tmp[i*np+np-i-1-j] = surfVerts[cnt];
                    }
                }

                for (cnt = i = 0; i < np; ++i)
                {
                    for(j = 0; j < np-i; ++j,cnt++)
                    {
                        surfVerts[cnt] = tmp[i*np+j];
                    }
                }
            }

            /**
             * @brief Align this surface to a given vertex ID.
             */
            void Align(vector<int> vertId)
            {
                if (vertId[0] == this->vertId[0])
                {
                    if (vertId[1] == this->vertId[1] ||
                        vertId[1] == this->vertId[2])
                    {
                        if (vertId[1] == this->vertId[2])
                        {
                            Rotate(1);
                            Reflect();
                        }
                    }
                }
                else if (vertId[0] == this->vertId[1])
                {
                    if (vertId[1] == this->vertId[0] ||
                        vertId[1] == this->vertId[2])
                    {
                        if (vertId[1] == this->vertId[0])
                        {
                            Reflect();
                        }
                        else
                        {
                            Rotate(2);
                        }
                    }
                }
                else if (vertId[0] == this->vertId[2])
                {
                    if (vertId[1] == this->vertId[0] ||
                        vertId[1] == this->vertId[1])
                    {
                        if (vertId[1] == this->vertId[1])
                        {
                            Rotate(2);
                            Reflect();
                        }
                        else
                        {
                            Rotate(1);
                        }
                    }
                }
            }
        };

        /**
         * @brief A 2-dimensional three-sided element.
         */
        class Triangle : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Triangle(pConf, pNodeList, pTagList));
                vector<EdgeSharedPtr> m_edges = e->GetEdgeList();
                for (int i = 0; i < m_edges.size(); ++i)
                {
                    m_edges[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Triangle(ElmtConfig                 pConf,
                     std::vector<NodeSharedPtr> pNodeList,
                     std::vector<int>           pTagList);
            Triangle(const Triangle& pSrc);
            virtual ~Triangle() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            virtual void Complete(int order);

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };


        /**
         * @brief A 2-dimensional four-sided element.
         */
        class Quadrilateral : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Quadrilateral(pConf, pNodeList, pTagList));
                vector<EdgeSharedPtr> m_edges = e->GetEdgeList();
                for (int i = 0; i < m_edges.size(); ++i)
                {
                    m_edges[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Quadrilateral(ElmtConfig                 pConf,
                          std::vector<NodeSharedPtr> pNodeList,
                          std::vector<int>           pTagList);
            Quadrilateral(const Quadrilateral& pSrc);
            virtual ~Quadrilateral() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            virtual void Complete(int order);

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };


        /**
         * @brief A 3-dimensional four-faced element.
         */
        class Tetrahedron : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Tetrahedron(pConf, pNodeList, pTagList));
                vector<FaceSharedPtr> faces = e->GetFaceList();
                for (int i = 0; i < faces.size(); ++i)
                {
                    faces[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Tetrahedron(ElmtConfig                 pConf,
                        std::vector<NodeSharedPtr> pNodeList,
                        std::vector<int>           pTagList);
            Tetrahedron(const Tetrahedron& pSrc);
            virtual ~Tetrahedron() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            virtual void Complete(int order);

            static unsigned int GetNumNodes(ElmtConfig pConf);

            int m_orientationMap[4];
            int m_origVertMap[4];

        protected:
            void OrientTet();
        };


        /**
         * @brief A 3-dimensional square-based pyramidic element
         */
        class Pyramid : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Pyramid(pConf, pNodeList, pTagList));
                vector<FaceSharedPtr> faces = e->GetFaceList();
                for (int i = 0; i < faces.size(); ++i)
                {
                    faces[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType type;

            Pyramid(ElmtConfig                 pConf,
                    std::vector<NodeSharedPtr> pNodeList,
                    std::vector<int>           pTagList);
            Pyramid(const Pyramid& pSrc);
            virtual ~Pyramid() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            static unsigned int GetNumNodes(ElmtConfig pConf);

            /**
             * Orientation of pyramid.
             */
            int orientationMap[5];
        };

        /**
         * @brief A 3-dimensional five-faced element (2 triangles, 3
         * quadrilaterals).
         */
        class Prism : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Prism(pConf, pNodeList, pTagList));
                vector<FaceSharedPtr> faces = e->GetFaceList();
                for (int i = 0; i < faces.size(); ++i)
                {
                    faces[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Prism(ElmtConfig                 pConf,
                  std::vector<NodeSharedPtr> pNodeList,
                  std::vector<int>           pTagList);
            Prism(const Prism& pSrc);
            virtual ~Prism() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            virtual void Complete(int order);

            static unsigned int GetNumNodes(ElmtConfig pConf);

            /**
             * Orientation of prism; unchanged = 0; clockwise = 1;
             * counter-clockwise = 2. This is set by OrientPrism.
             */
            unsigned int m_orientation;

        protected:
            void OrientPrism();
        };


        /**
         * @brief A 3-dimensional six-faced element.
         */
        class Hexahedron : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(
                ElmtConfig                 pConf,
                std::vector<NodeSharedPtr> pNodeList,
                std::vector<int>           pTagList)
            {
                ElementSharedPtr e = boost::shared_ptr<Element>(
                    new Hexahedron(pConf, pNodeList, pTagList));
                vector<FaceSharedPtr> faces = e->GetFaceList();
                for (int i = 0; i < faces.size(); ++i)
                {
                    faces[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static LibUtilities::ShapeType m_type;

            Hexahedron(ElmtConfig                 pConf,
                       std::vector<NodeSharedPtr> pNodeList,
                       std::vector<int>           pTagList);
            Hexahedron(const Hexahedron& pSrc);
            virtual ~Hexahedron() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };
    }
}

#endif
