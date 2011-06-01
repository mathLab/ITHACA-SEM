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

#include <LibUtilities/BasicUtils/NekFactory.hpp>

namespace Nektar
{
    namespace Utilities
    {
        /// Defines a less-than operator between objects referred to using
        /// shared pointers.
        template <typename T>
        struct shared_ptr_less_than
        {
            typedef boost::shared_ptr<T> pT;
            const bool operator()(const pT a, const pT b) const {
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
                    return *a < *b;
                }
            }
        };

        /// Defines an equality operator between objects referred to using
        /// shared pointers.
        template <typename T>
        struct shared_ptr_equality
        {
            typedef boost::shared_ptr<T> pT;
            const bool operator()(const pT a, const pT b) const {
                return *a == *b;
            }
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
                : id(pId), x(pX), y(pY), z(pZ) {}
            /// Copy an existing node.
            Node(const Node& pSrc)
                : id(pSrc.id), x(pSrc.x), y(pSrc.y), z(pSrc.z) {}
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

            /// ID of node.
            unsigned int id;
            /// X-coordinate.
            double x;
            /// Y-coordinate.
            double y;
            /// Z-coordinate.
            double z;
        };
        /// Shared pointer to a Node.
        typedef boost::shared_ptr<Node> NodeSharedPtr;


        /**
         * @brief Represents an edge which joins two points.
         *
         * An edge is defined by two nodes (vertices) and, for high-order edges,
         * a set of control nodes defining the shape of the edge.
         */
        class Edge {
        public:
            /// Creates a new edge.
            Edge(NodeSharedPtr pVertex1, NodeSharedPtr pVertex2, std::vector<NodeSharedPtr> pEdgeNodes)
                : n1(pVertex1), n2(pVertex2), edgeNodes(pEdgeNodes) {}
            /// Copies an existing edge.
            Edge(const Edge& pSrc)
                : n1(pSrc.n1), n2(pSrc.n2), edgeNodes(pSrc.edgeNodes) {}
            ~Edge() {}

            /// Equality is defined based on the vertices.
            bool operator==(const Edge& pSrc)
            {
                return ( ((*n1 == *(pSrc.n1)) && (*n2 == *(pSrc.n2)))
                        || ((*n2 == *(pSrc.n1)) && (*n1 == *(pSrc.n2))));
            }

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
                s << std::scientific << std::setprecision(3) << "     "
                  <<  n1->x << "  " << n1->y << "  " << n1->z << "     ";
                for (int k = 0; k < edgeNodes.size(); ++k) {
                    s << std::scientific << std::setprecision(3) << "     "
                      <<  edgeNodes[k]->x << "  " << edgeNodes[k]->y
                      << "  " << edgeNodes[k]->z << "     ";
                }
                s << std::scientific << std::setprecision(3) << "     "
                  <<  n2->x << "  " << n2->y << "  " << n2->z;
                return s.str();
            }

            /// ID of edge.
            unsigned int id;
            /// First vertex node.
            NodeSharedPtr n1;
            /// Second vertex node.
            NodeSharedPtr n2;
            /// List of control nodes between the first and second vertices.
            std::vector<NodeSharedPtr> edgeNodes;
        };
        /// Shared pointer to an edge.
        typedef boost::shared_ptr<Edge> EdgeSharedPtr;


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
            Face(std::vector<NodeSharedPtr> pVertexList, std::vector<NodeSharedPtr> pFaceNodes,
                 std::vector<EdgeSharedPtr> pEdgeList)
                : vertexList(pVertexList), faceNodes(pFaceNodes), edgeList(pEdgeList) {}
            /// Copy an existing face.
            Face(const Face& pSrc)
                : vertexList (pSrc.vertexList), faceNodes(pSrc.faceNodes), edgeList(pSrc.edgeList) {}
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
                    s << std::scientific << std::setprecision(3) << "     "
                      <<  vertexList[k]->x << "  " << vertexList[k]->y
                      << "  " << vertexList[k]->z << "     ";
                }
                for (int k = 0; k < edgeList.size(); ++k) {
                    for (int i = 0; i < edgeList[k]->edgeNodes.size(); ++i)
                    s << std::scientific << std::setprecision(3) << "     "
                      << edgeList[k]->edgeNodes[i]->x << "  "
                      << edgeList[k]->edgeNodes[i]->y << "  "
                      << edgeList[k]->edgeNodes[i]->z << "     ";
                }
                for (int k = 0; k < faceNodes.size(); ++k) {
                    s << std::scientific << std::setprecision(3) << "     "
                      <<  faceNodes[k]->x << "  " << faceNodes[k]->y
                      << "  " << faceNodes[k]->z << "     ";
                }
                return s.str();
            }

            /// ID of the face.
            unsigned int id;
            /// List of vertex nodes.
            std::vector<NodeSharedPtr> vertexList;
            /// List of corresponding edges.
            std::vector<EdgeSharedPtr> edgeList;
            /// List of face-interior nodes defining the shape of the face.
            std::vector<NodeSharedPtr> faceNodes;
        };
        /// Shared pointer to a face.
        typedef boost::shared_ptr<Face> FaceSharedPtr;


        /**
         * @brief Base class for element definitions.
         *
         * An element is defined by a list of vertices, edges and faces
         * (depending on the dimension of the problem). This base class
         * provides the underlying structure.
         */
        class Element {
        public:
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

        protected:
            /// ID of the element.
            unsigned int m_id;
            /// Dimension of the element.
            unsigned int m_dim;
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
        };
        /// Shared pointer to an element.
        typedef boost::shared_ptr<Element> ElementSharedPtr;
        /// Element factory definition.
        typedef Nektar::LibUtilities::NekFactory< unsigned int, Element,
                std::vector<NodeSharedPtr>, std::vector<int> > ElementFactory;
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


        /**
         * @brief A 0-dimensional vertex.
         */
        class Point : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Point(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Point(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Point(const Point& pSrc);
            virtual ~Point() {}
        };


        /**
         * @brief A 1-dimensional line between two vertex nodes.
         */
        class Line : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Line(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Line(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Line(const Point& pSrc);
            virtual ~Line() {}
        };


        /**
         * @brief A 2-dimensional three-sided element.
         */
        class Triangle : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Triangle(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Triangle(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Triangle(const Triangle& pSrc);
            virtual ~Triangle() {}
        };


        /**
         * @brief A 2-dimensional four-sided element.
         */
        class Quadrilateral : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Quadrilateral(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Quadrilateral(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Quadrilateral(const Quadrilateral& pSrc);
            virtual ~Quadrilateral() {}
        };


        /**
         * @brief A 3-dimensional four-faced element.
         */
        class Tetrahedron : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Tetrahedron(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Tetrahedron(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Tetrahedron(const Tetrahedron& pSrc);
            virtual ~Tetrahedron() {}

        protected:
            void OrientTet();
        };


        /**
         * @brief A 3-dimensional six-faced element.
         */
        class Hexahedron : public Element {
        public:
            /// Creates an instance of this class
            static ElementSharedPtr create(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList) {
                return boost::shared_ptr<Element>(new Hexahedron(pNodeList, pTagList));
            }
            /// Gmsh IDs
            static unsigned int typeIds[];

            Hexahedron(std::vector<NodeSharedPtr> pNodeList, std::vector<int> pTagList);
            Hexahedron(const Hexahedron& pSrc);
            virtual ~Hexahedron() {}
        };
    }
}

#endif
