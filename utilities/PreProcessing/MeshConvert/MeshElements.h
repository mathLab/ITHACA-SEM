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

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/Curve.hpp>
#include <SpatialDomains/MeshComponents.h>

namespace Nektar
{
    namespace Utilities
    {
        enum ElementType {
            ePoint,
            eLine,
            eTriangle,
            eQuadrilateral,
            eTetrahedron,
            ePyramid,
            ePrism,
            eHexahedron,
            SIZE_ElementType
        };

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
            Node(unsigned int pId, double pX, double pY, double pZ)
                : id(pId), x(pX), y(pY), z(pZ), m_geom() {}
            /// Copy an existing node.
            Node(const Node& pSrc)
                : id(pSrc.id), x(pSrc.x), y(pSrc.y), 
                  z(pSrc.z), m_geom() {}
            Node() : id(0), x(0.0), y(0.0), z(0.0), m_geom() {}
            ~Node() {}

            /// Define node ordering based on ID.
            bool operator<(const Node& pSrc)
            {
                return (id < pSrc.id);
            }
            /// Define node equality based on coordinate.
            bool operator==(const Node& pSrc)
            {
                return x == pSrc.x && y == pSrc.y && z == pSrc.z;
            }
            
            Node operator+(const Node &pSrc) const
            {
                return Node(id, x+pSrc.x, y+pSrc.y, z+pSrc.z);
            }
            
            Node operator-(const Node &pSrc) const
            {
                return Node(id, x-pSrc.x, y-pSrc.y, z-pSrc.z);
            }

            Node operator*(const Node &pSrc) const
            {
                return Node(id, x*pSrc.x, y*pSrc.y, z*pSrc.z);
            }
            
            Node operator*(const double &alpha) const
            {
                return Node(id, alpha*x, alpha*y, alpha*z);
            }
            
            Node operator/(const double &alpha) const
            {
                return Node(id, x/alpha, y/alpha, z/alpha);
            }
            
            void operator+=(const Node &pSrc)
            {
                x += pSrc.x;
                y += pSrc.y;
                z += pSrc.z;
            }
            
            void operator*=(const double &alpha)
            {
                x *= alpha;
                y *= alpha;
                z *= alpha;
            }
            
            void operator/=(const double &alpha)
            {
                x /= alpha;
                y /= alpha;
                z /= alpha;
            }
            
            double abs2() const
            {
                return x*x+y*y+z*z;
            }

            double dot(const Node &pSrc) const
            {
                return x*pSrc.x + y*pSrc.y + z*pSrc.z;
            }

            /// Generate a %SpatialDomains::VertexComponent for this node.
            SpatialDomains::VertexComponentSharedPtr GetGeom(int coordDim)
            {
                if (m_geom)
                {
                    return m_geom;
                }
                
                m_geom = MemoryManager<SpatialDomains::VertexComponent>::
                    AllocateSharedPtr(coordDim,id,x,y,z);
                return m_geom;
            }
            
            /// ID of node.
            unsigned int id;
            /// X-coordinate.
            double x;
            /// Y-coordinate.
            double y;
            /// Z-coordinate.
            double z;
            
        private:
            SpatialDomains::VertexComponentSharedPtr m_geom;
        };
        /// Shared pointer to a Node.
        typedef boost::shared_ptr<Node> NodeSharedPtr;

        bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2);
        bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2);

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
                boost::hash_combine(seed, p->x);
                boost::hash_combine(seed, p->y);
                boost::hash_combine(seed, p->z);
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
                : n1(pVertex1), n2(pVertex2), edgeNodes(pEdgeNodes),
                  curveType(pCurveType), m_geom() {}
            /// Copies an existing edge.
            Edge(const Edge& pSrc)
                : n1(pSrc.n1), n2(pSrc.n2), edgeNodes(pSrc.edgeNodes),
                  curveType(pSrc.curveType), m_geom(pSrc.m_geom) {}
            ~Edge() {}

            /// Returns the total number of nodes defining the edge.
            unsigned int GetNodeCount() const
            {
                return edgeNodes.size() + 2;
            }

            /// Creates a Nektar++ string listing the coordinates of all the
            /// nodes.
            std::string GetXmlCurveString() const
            {
                std::stringstream s;
                std::string str;
                s << std::scientific << std::setprecision(8) << "     "
                  <<  n1->x << "  " << n1->y << "  " << n1->z << "     ";
                for (int k = 0; k < edgeNodes.size(); ++k) {
                    s << std::scientific << std::setprecision(8) << "     "
                      <<  edgeNodes[k]->x << "  " << edgeNodes[k]->y
                      << "  " << edgeNodes[k]->z << "     ";
                }
                s << std::scientific << std::setprecision(8) << "     "
                  <<  n2->x << "  " << n2->y << "  " << n2->z;
                return s.str();
            }

            /// Generate a SpatialDomains::SegGeom object for this edge.
            SpatialDomains::SegGeomSharedPtr GetGeom(int coordDim)
            {
                if (m_geom)
                {
                    return m_geom;
                }
                
                // Create edge vertices.
                SpatialDomains::VertexComponentSharedPtr p[2];
                p[0] = n1->GetGeom(coordDim);
                p[1] = n2->GetGeom(coordDim);
                
                // Create a curve if high-order information exists.
                if (edgeNodes.size() > 0)
                {
                    SpatialDomains::CurveSharedPtr c = 
                        MemoryManager<SpatialDomains::Curve>::
                        AllocateSharedPtr(id, curveType);
                    
                    c->m_points.push_back(p[0]);
                    for (int i = 0; i < edgeNodes.size(); ++i)
                    {
                        c->m_points.push_back(edgeNodes[i]->GetGeom(coordDim));
                    }
                    c->m_points.push_back(p[1]);
                    
                    m_geom = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(id, coordDim, p, c);
                }
                else
                {
                    m_geom = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(id, coordDim, p);
                }
                
                return m_geom;
            }

            /// ID of edge.
            unsigned int id;
            /// First vertex node.
            NodeSharedPtr n1;
            /// Second vertex node.
            NodeSharedPtr n2;
            /// List of control nodes between the first and second vertices.
            std::vector<NodeSharedPtr> edgeNodes;
            /// Distributions of points along edge.
            LibUtilities::PointsType curveType;

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
                unsigned int id1 = p->n1->id;
                unsigned int id2 = p->n2->id;
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
                : vertexList(pVertexList), 
                  edgeList  (pEdgeList),
                  faceNodes (pFaceNodes), 
                  curveType (pCurveType),
                  m_geom    () {}
            
            /// Copy an existing face.
            Face(const Face& pSrc)
                : vertexList(pSrc.vertexList), edgeList  (pSrc.edgeList),
                  faceNodes (pSrc.faceNodes),  curveType (pSrc.curveType),
                  m_geom    (pSrc.m_geom) {}
            ~Face() {}

            /// Equality is defined by matching all vertices.
            bool operator==(Face& pSrc)
            {
                std::vector<NodeSharedPtr>::iterator it1, it2;
                for (it1 = vertexList.begin(); it1 != vertexList.end(); ++it1)
                {
                    if (find(pSrc.vertexList.begin(), pSrc.vertexList.end(), *it1)
                            == pSrc.vertexList.end())
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
                unsigned int n = faceNodes.size();
                for (int i = 0; i < edgeList.size(); ++i)
                {
                    n += edgeList[i]->GetNodeCount();
                }
                n -= vertexList.size();
                return n;
            }

            /// Generates a string listing the coordinates of all nodes
            /// associated with this face.
            std::string GetXmlCurveString() const
            {
                std::stringstream s;
                std::string str;
                for (int k = 0; k < vertexList.size(); ++k) {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  vertexList[k]->x << "  " << vertexList[k]->y
                      << "  " << vertexList[k]->z << "    ";
                }
                for (int k = 0; k < edgeList.size(); ++k) {
                    for (int i = 0; i < edgeList[k]->edgeNodes.size(); ++i)
                    s << std::scientific << std::setprecision(8) << "    "
                      << edgeList[k]->edgeNodes[i]->x << "  "
                      << edgeList[k]->edgeNodes[i]->y << "  "
                      << edgeList[k]->edgeNodes[i]->z << "    ";
                }
                for (int k = 0; k < faceNodes.size(); ++k) {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  faceNodes[k]->x << "  " << faceNodes[k]->y
                      << "  " << faceNodes[k]->z << "    ";
                }
                return s.str();
            }

            /// Generate either SpatialDomains::TriGeom or
            /// SpatialDomains::QuadGeom for this element.
            SpatialDomains::Geometry2DSharedPtr GetGeom(int coordDim)
            {
                if (m_geom)
                {
                    return m_geom;
                }
                
                int nEdge = edgeList.size();
                
                SpatialDomains::SegGeomSharedPtr edges[4];
                StdRegions::Orientation      edgeo[4];
                
                for (int i = 0; i < nEdge; ++i)
                {
                    edges[i] = edgeList[i]->GetGeom(coordDim);
                }
                
                for (int i = 0; i < nEdge; ++i)
                {
                    edgeo[i] = SpatialDomains::SegGeom::GetEdgeOrientation(
                        *edges[i], *edges[(i+1) % nEdge]);
                }
                
                if (nEdge == 3)
                {
                    m_geom = MemoryManager<SpatialDomains::TriGeom>::
                        AllocateSharedPtr(id, edges, edgeo);
                }
                else
                {
                    m_geom = MemoryManager<SpatialDomains::QuadGeom>::
                        AllocateSharedPtr(id, edges, edgeo);
                }

                return m_geom;
            }
            
            /// ID of the face.
            unsigned int id;
            /// List of vertex nodes.
            std::vector<NodeSharedPtr> vertexList;
            /// List of corresponding edges.
            std::vector<EdgeSharedPtr> edgeList;
            /// List of face-interior nodes defining the shape of the face.
            std::vector<NodeSharedPtr> faceNodes;
            /// Distribution of points in this face.
            LibUtilities::PointsType   curveType;
            /// Element(s) which are linked to this face.
            vector<pair<ElementSharedPtr, int> > elLink; 
            
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
                unsigned int              nVert = p->vertexList.size();
                std::size_t               seed  = 0;
                std::vector<unsigned int> ids(nVert);
                
                for (int i = 0; i < nVert; ++i)
                {
                    ids[i] = p->vertexList[i]->id;
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
            ElmtConfig(ElementType pE, unsigned int pOrder, 
                       bool pFn, bool pVn, bool pReorient = true,
                       LibUtilities::PointsType pECt=LibUtilities::ePolyEvenlySpaced,
                       LibUtilities::PointsType pFCt=LibUtilities::ePolyEvenlySpaced):
                e(pE), faceNodes(pFn), volumeNodes(pVn), order(pOrder),
                reorient(pReorient), edgeCurveType(pECt), faceCurveType(pFCt) {}
            ElmtConfig(ElmtConfig const &p) :
                e(p.e), faceNodes(p.faceNodes), volumeNodes(p.volumeNodes), 
                order(p.order), reorient(p.reorient), 
                edgeCurveType(p.edgeCurveType), faceCurveType(p.faceCurveType) {}

            ElmtConfig() {}
            
            /// Element type (e.g. triangle, quad, etc).
            ElementType              e;
            /// Denotes whether the element contains face nodes. For 2D
            /// elements, if this is true then the element contains interior
            /// nodes.
            bool                     faceNodes;
            /// Denotes whether the element contains volume (i.e. interior)
            /// nodes. These are not supported by either the mesh converter or
            /// Nektar++ but are included for completeness and are required
            /// for some output modules (e.g. Gmsh).
            bool                     volumeNodes;
            /// Order of the element.
            unsigned int             order;
            /// Denotes whether the element needs to be re-orientated for a
            /// spectral element framework.
            bool                     reorient;
            /// Distribution of points in edges.
            LibUtilities::PointsType edgeCurveType;
            /// Distribution of points in faces.
            LibUtilities::PointsType faceCurveType;
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
                if (m_faceLink.get() != 0) return m_faceLink->id;
                if (m_edgeLink.get() != 0) return m_edgeLink->id;
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
                return vertex[i];
            }
            /// Access an edge.
            EdgeSharedPtr GetEdge(unsigned int i) const {
                return edge[i];
            }
            /// Access a face.
            FaceSharedPtr GetFace(unsigned int i) const {
                return face[i];
            }
            /// Access the list of vertex nodes.
            std::vector<NodeSharedPtr> GetVertexList() const {
                return vertex;
            }
            /// Access the list of edges.
            std::vector<EdgeSharedPtr> GetEdgeList() const {
                return edge;
            }
            /// Access the list of faces.
            std::vector<FaceSharedPtr> GetFaceList() const {
                return face;
            }
            /// Access the list of volume nodes.
            std::vector<NodeSharedPtr> GetVolumeNodes() const {
                return volumeNodes;
            }
            /// Returns the total number of nodes (vertices, edge nodes and
            /// face nodes and volume nodes).
            unsigned int GetNodeCount() const
            {
                unsigned int n = volumeNodes.size();
                if (m_dim == 2)
                {
                    for (int i = 0; i < edge.size(); ++i)
                    {
                        n += edge[i]->GetNodeCount();
                    }
                    n -= vertex.size();
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
                return vertex.size();
            }
            /// Returns the number of edges.
            unsigned int GetEdgeCount() const {
                return edge.size();
            }
            /// Returns the number of faces.
            unsigned int GetFaceCount() const {
                return face.size();
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
            /// Set a correspondence between this element and a face
            /// (3D boundary element).
            void SetFaceLink(FaceSharedPtr pLink) {
                m_faceLink = pLink;
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
            void SetTagList(const std::vector<int> &tags) {
                m_taglist = tags;
            }
            /// Generate a list of vertices (1D), edges (2D), or faces (3D).
            virtual std::string GetXmlString() const 
            {
                std::stringstream s;
                switch (m_dim)
                {
                case 1:
                    for(int j=0; j< vertex.size(); ++j)
                    {
                        s << std::setw(5) << vertex[j]->id << " ";
                    }
                    break;
                case 2:
                    {
                        NekDouble cross;
                        
                        // caclulate sign based on cross product of vertices
                        if(edge[0]->n1 == edge[1]->n1)
                        {
                            cross  = (edge[0]->n2->x - edge[0]->n1->x)*
                                (edge[1]->n2->y - edge[1]->n1->y) - 
                                (edge[0]->n2->y - edge[0]->n1->y)*
                                (edge[1]->n2->x - edge[1]->n1->x); 
                        }
                        else if(edge[0]->n1 == edge[1]->n2)
                        {
                            cross  = (edge[0]->n2->x - edge[0]->n1->x)*
                                (edge[1]->n1->y - edge[1]->n2->y) - 
                                (edge[0]->n2->y - edge[0]->n1->y)*
                                (edge[1]->n1->x - edge[1]->n2->x); 
                        }
                        else if(edge[0]->n2 == edge[1]->n1)
                        {
                            cross  = (edge[0]->n1->x - edge[0]->n2->x)*
                                (edge[1]->n2->y - edge[1]->n1->y) - 
                                (edge[0]->n1->y - edge[0]->n2->y)*
                                (edge[1]->n2->x - edge[1]->n1->x); 

                        }
                        else if(edge[0]->n2 == edge[1]->n2)
                        {
                            cross  = (edge[0]->n1->x - edge[0]->n2->x)*
                                (edge[1]->n1->y - edge[1]->n2->y) - 
                                (edge[0]->n1->y - edge[0]->n2->y)*
                                (edge[1]->n1->x - edge[1]->n2->x); 
                        }
                        
                        // provide edges in anticlockwise sense
                        if(cross  < 0.0)
                        {
                            for(int j=0; j< edge.size(); ++j)
                            {
                                s << std::setw(5) << edge[j]->id << " ";
                            }
                        }
                        else
                        {
                            for(int j=edge.size()-1; j>=0; --j)
                            {
                                s << std::setw(5) << edge[j]->id << " ";
                            }
                        }
                    }
                    break;
                case 3:
                    for(int j=0; j< face.size(); ++j)
                    {
                        s << std::setw(5) << face[j]->id << " ";
                    }
                    break;
                }
                return s.str();
            }

            /**
             * @brief Reorders Gmsh ordered nodes into a row-by-row ordering
             * required for Nektar++ curve tags.
             *
             * The interior nodes of elements in the Gmsh format are ordered
             * as for a lower-order element of the same type. This promotes
             * the recursive approach to the reordering algorithm.
             */
            std::vector<NodeSharedPtr> tensorNodeOrdering(
                    const std::vector<NodeSharedPtr> &nodes,
                    int n) const
            {
                std::vector<NodeSharedPtr> nodeList;
                int cnt2;

                // Triangle
                if (vertex.size() == 3)
                {
                    nodeList.resize(nodes.size());

                    // Vertices
                    nodeList[0] = nodes[0];
                    if (n > 1)
                    {
                        nodeList[n-1] = nodes[1];
                        nodeList[n*(n+1)/2 - 1] = nodes[2];
                    }

                    // Edges
                    int cnt = n;
                    for (int i = 1; i < n-1; ++i)
                    {
                        nodeList[i] = nodes[3+i-1];
                        nodeList[cnt] = nodes[3+3*(n-2)-i];
                        nodeList[cnt+n-i-1] = nodes[3+(n-2)+i-1];
                        cnt += n-i;
                    }

                    // Interior (recursion)
                    if (n > 3)
                    {
                        // Reorder interior nodes
                        std::vector<NodeSharedPtr> interior((n-3)*(n-2)/2);
                        std::copy(nodes.begin() + 3+3*(n-2), nodes.end(), interior.begin());
                        interior = tensorNodeOrdering(interior, n-3);

                        // Copy into full node list
                        cnt = n;
                        cnt2 = 0;
                        for (int j = 1; j < n-2; ++j)
                        {
                            for (int i = 0; i < n-j-2; ++i)
                            {
                                nodeList[cnt+i+1] = interior[cnt2+i];
                            }
                            cnt += n-j;
                            cnt2 += n-2-j;
                        }
                    }
                }
                // Quad
                else if (m_dim == 2 && vertex.size() == 4)
                {
                    nodeList.resize(nodes.size());

                    // Vertices and edges
                    nodeList[0] = nodes[0];
                    if (n > 1)
                    {
                        nodeList[n-1] = nodes[1];
                        nodeList[n*n-1] = nodes[2];
                        nodeList[n*(n-1)] = nodes[3];
                    }
                    for (int i = 1; i < n-1; ++i)
                    {
                        nodeList[i] = nodes[4+i-1];
                    }
                    for (int i = 1; i < n-1; ++i)
                    {
                        nodeList[n*n-1-i] = nodes[4+2*(n-2)+i-1];
                    }

                    // Interior (recursion)
                    if (n > 2)
                    {
                        // Reorder interior nodes
                        std::vector<NodeSharedPtr> interior((n-2)*(n-2));
                        std::copy(nodes.begin() + 4+4*(n-2), nodes.end(), interior.begin());
                        interior = tensorNodeOrdering(interior, n-2);

                        // Copy into full node list
                        for (int j = 1; j < n-1; ++j)
                        {
                            nodeList[j*n] = nodes[4+3*(n-2)+n-2-j];
                            for (int i = 1; i < n-1; ++i)
                            {
                                nodeList[j*n+i] = interior[(j-1)*(n-2)+(i-1)];
                            }
                            nodeList[(j+1)*n-1] = nodes[4+(n-2)+j-1];
                        }
                    }
                }
                else
                {
                    cerr << "TensorNodeOrdering for a " << vertex.size()
                         << "-vertex element is not yet implemented." << endl;
                }
                return nodeList;
            }

            /// Generates a string listing the coordinates of all nodes
            /// associated with this element.
            std::string GetXmlCurveString() const
            {
                // Temporary node list for reordering
                std::vector<NodeSharedPtr> nodeList;

                // Node orderings are different for different elements.
                // Triangle
                if (vertex.size() == 3)
                {
                    int n = edge[0]->GetNodeCount();
                    nodeList.resize(n*(n+1)/2);

                    // Populate nodelist
                    std::copy(vertex.begin(), vertex.end(), nodeList.begin());
                    for (int i = 0; i < 3; ++i)
                    {
                        std::copy(edge[i]->edgeNodes.begin(), edge[i]->edgeNodes.end(), nodeList.begin() + 3 + i*(n-2));
                        if (edge[i]->n1 != vertex[i])
                        {
                            // If edge orientation is reversed relative to node
                            // ordering, we need to reverse order of nodes.
                            std::reverse(nodeList.begin() + 3 + i*(n-2),
                                         nodeList.begin() + 3 + (i+1)*(n-2));
                        }
                    }

                    // Triangle ordering lists vertices, edges then interior.
                    // Interior nodes are row by row from edge 0 up to vertex 2
                    // so need to reorder interior nodes only.
                    std::vector<NodeSharedPtr> interior(volumeNodes.size());
                    std::copy(volumeNodes.begin(), volumeNodes.end(), interior.begin());
                    interior = tensorNodeOrdering(interior, n-3);
                    std::copy(interior.begin(), interior.end(), nodeList.begin() + 3*(n-1));
                }
                // Quad
                else if (m_dim == 2 && vertex.size() == 4)
                {
                    int n = edge[0]->GetNodeCount();
                    nodeList.resize(n*n);

                    // Populate nodelist
                    std::copy(vertex.begin(), vertex.end(), nodeList.begin());
                    for (int i = 0; i < 4; ++i)
                    {
                        std::copy(edge[i]->edgeNodes.begin(),
                                  edge[i]->edgeNodes.end(),
                                  nodeList.begin() + 4 + i*(n-2));

                        if (edge[i]->n1 != vertex[i])
                        {
                            // If edge orientation is reversed relative to node
                            // ordering, we need to reverse order of nodes.
                            std::reverse(nodeList.begin() + 4 + i*(n-2),
                                         nodeList.begin() + 4 + (i+1)*(n-2));
                        }
                    }
                    std::copy(volumeNodes.begin(), volumeNodes.end(), nodeList.begin() + 4*(n-1));

                    // Quadrilateral ordering lists all nodes row by row
                    // starting from edge 0 up to edge 2, so need to reorder
                    // all nodes.
                    nodeList = tensorNodeOrdering(nodeList, n);
                }
                else
                {
                    cerr << "GetXmlCurveString for a " << vertex.size()
                         << "-vertex element is not yet implemented." << endl;
                }

                // Finally generate the XML string corresponding to our new
                // node reordering.
                std::stringstream s;
                std::string str;
                for (int k = 0; k < nodeList.size(); ++k)
                {
                    s << std::scientific << std::setprecision(8) << "    "
                      <<  nodeList[k]->x << "  " << nodeList[k]->y
                      << "  " << nodeList[k]->z << "    ";

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
                for (i = 0; i < vertex.size(); ++i)
                {
                    cout << vertex[i]->x << " " << vertex[i]->y << " " << vertex[i]->z << endl;
                }
                for (i = 0; i < edge.size(); ++i)
                {
                    for (j = 0; j < edge[i]->edgeNodes.size(); ++j)
                    {
                        NodeSharedPtr n = edge[i]->edgeNodes[j];
                        cout << n->x << " " << n->y << " " << n->z << endl;
                    }
                }
                for (i = 0; i < face.size(); ++i)
                {
                    for (j = 0; j < face[i]->faceNodes.size(); ++j)
                    {
                        NodeSharedPtr n = face[i]->faceNodes[j];
                        cout << n->x << " " << n->y << " " << n->z << endl;
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
            std::vector<NodeSharedPtr> vertex;
            /// List of element edges.
            std::vector<EdgeSharedPtr> edge;
            /// List of element faces.
            std::vector<FaceSharedPtr> face;
            /// List of element volume nodes.
            std::vector<NodeSharedPtr> volumeNodes;
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
        typedef Nektar::LibUtilities::NekFactory<ElementType, Element,
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
            Composite() {}

            /// Generate the list of IDs of elements within this composite.
            std::string GetXmlString(bool doSort=true);

            /// ID of composite.
            unsigned int id;
            /// Element type tag.
            std::string tag;
            /// List of elements in this composite.
            std::vector<ElementSharedPtr> items;
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
            Condition() : type(), field(), value(), composite() {}
            std::vector<ConditionType> type;
            std::vector<std::string>   field;
            std::vector<std::string>   value;
            std::vector<int>           composite;
        };

        typedef boost::shared_ptr<Condition> ConditionSharedPtr;
        typedef std::map<int,ConditionSharedPtr> ConditionMap;

        bool operator==(ConditionSharedPtr const &c1, ConditionSharedPtr const &c2);
        
        class Mesh
        {
        public:
            Mesh() : verbose(false) {}
            
            /// Verbose flag
            bool                       verbose;
            /// Dimension of the expansion.
            unsigned int               expDim;
            /// Dimension of the space in which the mesh is defined.
            unsigned int               spaceDim;
            /// List of mesh nodes.
            std::vector<NodeSharedPtr> node;
            /// Set of element vertices.
            NodeSet                    vertexSet;
            /// Set of element edges.
            EdgeSet                    edgeSet;
            /// Set of element faces.
            FaceSet                    faceSet;
            /// Map for elements.
            ElementMap                 element;
            /// Map for composites.
            CompositeMap               composite;
            /// Boundary conditions maps tag to condition.
            ConditionMap               condition;
            /// List of fields names.
            std::vector<std::string>   fields;
            /// Map of vertex normals.
            boost::unordered_map<int, Node> vertexNormals;
            /// Set of all pairs of element ID and face number on which to apply
            /// spherigon surface smoothing.
            set<pair<int,int> > spherigonFaces;
            /// Returns the total number of elements in the mesh with
            /// dimension expDim.
            unsigned int               GetNumElements();
            /// Returns the total number of elements in the mesh with
            /// dimension < expDim.
            unsigned int               GetNumBndryElements();
            /// Returns the total number of entities in the mesh.
            unsigned int               GetNumEntities();
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
            static ElementType type;

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
            static ElementType type;

            Line(ElmtConfig                 pConf,
                 std::vector<NodeSharedPtr> pNodeList, 
                 std::vector<int>           pTagList);
            Line(const Point& pSrc);
            virtual ~Line() {}
            
            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            
            static unsigned int GetNumNodes(ElmtConfig pConf);
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
                return boost::shared_ptr<Element>(
                    new Triangle(pConf, pNodeList, pTagList));
            }
            /// Element type
            static ElementType type;

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
                return boost::shared_ptr<Element>(
                    new Quadrilateral(pConf, pNodeList, pTagList));
            }
            /// Element type
            static ElementType type;

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
                    faces[i]->elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static ElementType type;

            Tetrahedron(ElmtConfig                 pConf,
                        std::vector<NodeSharedPtr> pNodeList,
                        std::vector<int>           pTagList);
            Tetrahedron(const Tetrahedron& pSrc);
            virtual ~Tetrahedron() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
            virtual void Complete(int order);
            
            static unsigned int GetNumNodes(ElmtConfig pConf);

            /**
             * Orientation of tet; unchanged = 0; base vertex swapped = 1.
             */
            int orientationMap[4];

        protected:
            void OrientTet();
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
                    faces[i]->elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static ElementType type;

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
            unsigned int orientation;

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
                    faces[i]->elLink.push_back(pair<ElementSharedPtr, int>(e,i));
                }
                return e;
            }
            /// Element type
            static ElementType type;

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
