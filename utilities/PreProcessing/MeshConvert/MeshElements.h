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
#include <sstream>

#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>

#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include <SpatialDomains/InterfaceComponent.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/TetGeom.h>
#include <SpatialDomains/PyrGeom.h>
#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/HexGeom.h>
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
            ~Node() {}

            /// Define node ordering based on ID.
            bool operator<(const Node& pSrc)
            {
                return (id < pSrc.id);
            }
            /// Define node equality based on coordinate.
            bool operator==(const Node& pSrc)
            {
                return ((x==pSrc.x) && (y==pSrc.y) && (z==pSrc.z));
            }

            /// Generate a %SpatialDomains::VertexComponent for this node.
            SpatialDomains::VertexComponentSharedPtr GetGeom()
            {
                if (!m_geom)
                {
                    m_geom = MemoryManager<SpatialDomains::VertexComponent>::
                        AllocateSharedPtr(3,id,x,y,z);
                }
                
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
        
        struct NodeHash : std::unary_function<NodeSharedPtr, std::size_t>
        {
            std::size_t operator()(NodeSharedPtr const& p) const
            {
                std::size_t seed = 0;
                boost::hash_combine(seed, p -> x);
                boost::hash_combine(seed, p -> y);
                boost::hash_combine(seed, p -> z);
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

            /// Creates a string listing the coordinates of all the nodes.
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

            /// Generate a %SpatialDomains::SegGeom object for this edge.
            SpatialDomains::SegGeomSharedPtr GetGeom()
            {
                if (m_geom)
                {
                    return m_geom;
                }
                
                // Create edge vertices.
                SpatialDomains::VertexComponentSharedPtr p[2];
                p[0] = n1->GetGeom();
                p[1] = n2->GetGeom();
                
                // Create a curve if high-order information exists.
                if (edgeNodes.size() > 0)
                {
                    SpatialDomains::CurveSharedPtr c = 
                        MemoryManager<SpatialDomains::Curve>::
                        AllocateSharedPtr(id, curveType);
                    
                    c->m_points.push_back(p[0]);
                    for (int i = 0; i < edgeNodes.size(); ++i)
                    {
                        c->m_points.push_back(edgeNodes[i]->GetGeom());
                    }
                    c->m_points.push_back(p[1]);
                    
                    m_geom = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(id, 3, p, c);
                }
                else
                {
                    m_geom = MemoryManager<SpatialDomains::SegGeom>::
                        AllocateSharedPtr(id, 3, p);
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
                  faceNodes (pFaceNodes), 
                  edgeList  (pEdgeList),
                  curveType (pCurveType),
                  m_geom    () {}
            
            /// Copy an existing face.
            Face(const Face& pSrc)
                : vertexList(pSrc.vertexList), faceNodes(pSrc.faceNodes), 
                  edgeList  (pSrc.edgeList),   curveType(pSrc.curveType),
                  m_geom    (pSrc.m_geom) {}
            ~Face() {}

            /// Equality is defined by matching all vertices.
            bool operator==(Face& pSrc)
            {
                bool e = true;
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

            /// Generate either %SpatialDomains::TriGeom or
            /// %SpatialDomains::QuadGeom for this element.
            SpatialDomains::Geometry2DSharedPtr GetGeom()
            {
                if (m_geom)
                {
                    return m_geom;
                }
                
                int nEdge = edgeList.size();
                
                SpatialDomains::SegGeomSharedPtr edges[4];
                StdRegions::EdgeOrientation      edgeo[4];
                
                for (int i = 0; i < nEdge; ++i)
                {
                    edges[i] = edgeList[i]->GetGeom();
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
            /// Set the list of tags associated with this element.
            void SetTagList(const std::vector<int> &tags) {
                m_taglist = tags;
            }
            /// Generate a list of vertices (1D), edges (2D), or faces (3D).
            virtual std::string GetXmlString() const {
                std::stringstream s;
                switch (m_dim)
                {
                case 1:
                    for(int j=0; j< vertex.size(); ++j){
                        s << std::setw(5) << vertex[j]->id << " ";
                    }
                    break;
                case 2:
                    for(int j=0; j< edge.size(); ++j){
                        s << std::setw(5) << edge[j]->id << " ";
                    }
                    break;
                case 3:
                    for(int j=0; j< face.size(); ++j){
                        s << std::setw(5) << face[j]->id << " ";
                    }
                    break;
                }
                return s.str();
            }
            /// Generate a Nektar++ geometry object for this element.
            virtual SpatialDomains::GeometrySharedPtr GetGeom()
            {
                ASSERTL0(false, "This function should be implemented on a shape level.");
                return boost::shared_ptr<SpatialDomains::Geometry>();
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
            /// Pointer to the corresponding edge if element is a 2D boundary.
            EdgeSharedPtr m_edgeLink;
            /// Pointer to the corresponding face if element is a 3D boundary.
            FaceSharedPtr m_faceLink;
            /// Nektar++ geometry object for this element.
            SpatialDomains::GeometrySharedPtr m_geom;
        };
        /// Shared pointer to an element.
        typedef boost::shared_ptr<Element> ElementSharedPtr;
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
            const bool operator()(const pT a, const pT b) const
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
            /// Generate the list of IDs of elements within this composite.
            std::string GetXmlString();

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

        enum ConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic,
            eHOPCondition,
            SIZE_ConditionType
        };

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
            Mesh(const std::string inFilename, 
                 const std::string outFilename);

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
            /// Original filename.
            std::string                inFilename;
            /// Intended target.
            std::string                outFilename;
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
            
            virtual SpatialDomains::GeometrySharedPtr GetGeom();
            
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
            
            virtual SpatialDomains::GeometrySharedPtr GetGeom();

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

            virtual SpatialDomains::GeometrySharedPtr GetGeom();

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
                return boost::shared_ptr<Element>(
                    new Tetrahedron(pConf, pNodeList, pTagList));
            }
            /// Element type
            static ElementType type;

            Tetrahedron(ElmtConfig                 pConf,
                        std::vector<NodeSharedPtr> pNodeList,
                        std::vector<int>           pTagList);
            Tetrahedron(const Tetrahedron& pSrc);
            virtual ~Tetrahedron() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom();

            static unsigned int GetNumNodes(ElmtConfig pConf);

            /**
             * Orientation of tet; unchanged = 0; base vertex swapped = 1.
             */
            unsigned int orientation;

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
                return boost::shared_ptr<Element>(
                    new Prism(pConf, pNodeList, pTagList));
            }
            /// Element type
            static ElementType type;

            Prism(ElmtConfig                 pConf,
                  std::vector<NodeSharedPtr> pNodeList,
                  std::vector<int>           pTagList);
            Prism(const Prism& pSrc);
            virtual ~Prism() {}

            virtual SpatialDomains::GeometrySharedPtr GetGeom();

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
                return boost::shared_ptr<Element>(
                    new Hexahedron(pConf, pNodeList, pTagList));
            }
            /// Element type
            static ElementType type;

            Hexahedron(ElmtConfig                 pConf,
                       std::vector<NodeSharedPtr> pNodeList,
                       std::vector<int>           pTagList);
            Hexahedron(const Hexahedron& pSrc);
            virtual ~Hexahedron() {}
            
            virtual SpatialDomains::GeometrySharedPtr GetGeom();

            static unsigned int GetNumNodes(ElmtConfig pConf);
        };
    }
}

#endif
