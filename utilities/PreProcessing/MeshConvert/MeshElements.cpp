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

#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <loki/Singleton.h>

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

        /**
         * When an edge is replaced, the element faces are also searched and
         * the corresponding face edges are updated to maintain consistency.
         */
        void Element::SetEdge(unsigned int p, EdgeSharedPtr pNew)
        {
            EdgeSharedPtr vOld = edge[p];
            edge[p] = pNew;
            for (int i = 0; i < face.size(); ++i)
            {
                for (int j = 0; j < face[i]->edgeList.size(); ++j)
                {
                    if (*face[i]->edgeList[j] == *vOld)
                    {
                        face[i]->edgeList[j] = pNew;
                    }
                }
            }
        }


        /**
         * When a face is replaced, no other consistency checks are required.
         */
        void Element::SetFace(unsigned int p, FaceSharedPtr pNew)
        {
            face[p] = pNew;
        }


        /**
         * The list of composites may include individual element IDs or ranges
         * of element IDs.
         */
        string Composite::GetXmlString()
        {
            element_id_less_than sortOperator;
            sort(items.begin(), items.end(), sortOperator);

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

        unsigned int Point::typeIds[1] = {
                GetElementFactory().RegisterCreatorFunction(15, Point::create, "Point")
            };
        Point::Point(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element() {
            m_tag = "";
            m_dim = 0;
            m_taglist = pTagList;
            vertex.push_back(pNodeList[0]);
        }


        unsigned int Line::typeIds[5] = {
                GetElementFactory().RegisterCreatorFunction(1,  Line::create, "Order 1 Line"),
                GetElementFactory().RegisterCreatorFunction(8,  Line::create, "Order 2 Line"),
                GetElementFactory().RegisterCreatorFunction(26, Line::create, "Order 3 Line"),
                GetElementFactory().RegisterCreatorFunction(27, Line::create, "Order 4 Line"),
                GetElementFactory().RegisterCreatorFunction(28, Line::create, "Order 5 Line")
            };
        Line::Line(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element() {
            m_tag = "S";
            m_dim = 1;
            m_taglist = pTagList;
            int n = 0;
            switch (m_taglist.back())
            {
            case 1:  n = 0; break;
            case 8:  n = 1; break;
            case 26: n = 2; break;
            case 27: n = 3; break;
            case 28: n = 4; break;
            }
            int order = n+1;

            // Add vertices
            for (int i = 0; i < 2; ++i) {
                vertex.push_back(pNodeList[i]);
            }
            vector<NodeSharedPtr> edgeNodes;
            if (order > 1) {
                for (int j = 0; j<n; ++j) {
                    edgeNodes.push_back(pNodeList[2+j]);
                }
            }
            edge.push_back(boost::shared_ptr<Edge>(new Edge(pNodeList[0], pNodeList[1], edgeNodes)));

        }


        unsigned int Triangle::typeIds[5] = {
                GetElementFactory().RegisterCreatorFunction(2,  Triangle::create, "Order 1 Triangle"),
                GetElementFactory().RegisterCreatorFunction(9,  Triangle::create, "Order 2 Triangle"),
                GetElementFactory().RegisterCreatorFunction(21, Triangle::create, "Order 3 Triangle"),
                GetElementFactory().RegisterCreatorFunction(23, Triangle::create, "Order 4 Triangle"),
                GetElementFactory().RegisterCreatorFunction(25, Triangle::create, "Order 5 Triangle")
            };
        Triangle::Triangle(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element() {
            m_tag = "T";
            m_dim = 2;
            m_taglist = pTagList;
            int n = 0;
            switch (m_taglist.back())
            {
            case 2:  n = 0; break;
            case 9:  n = 1; break;
            case 21: n = 2; break;
            case 23: n = 3; break;
            case 25: n = 4; break;
            }
            int order = n+1;

            // Create a map to relate edge nodes to a pair of vertices defining an edge
            // This is based on the ordering produced by gmsh.
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
                if (order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(  pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes)));
            }
        }


        unsigned int Quadrilateral::typeIds[3] = {
                GetElementFactory().RegisterCreatorFunction(3,  Quadrilateral::create, "Order 1 Quadrilateral"),
                GetElementFactory().RegisterCreatorFunction(10, Quadrilateral::create, "Order 2 Quadrilateral"),
                GetElementFactory().RegisterCreatorFunction(16, Quadrilateral::create, "Order 3 Quadrilateral"),
            };
        Quadrilateral::Quadrilateral(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element() {
            m_tag = "Q";
            m_dim = 2;
            m_taglist = pTagList;
            int n = 0;
            switch (m_taglist.back())
            {
            case 3:  n = 0; break;
            case 10: n = 1; break;
            case 16: n = 2; break;
            }
            int order = n+1;

            // Create a map to relate edge nodes to a pair of vertices defining an edge
            // This is based on the ordering produced by gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 5;
            edgeNodeMap[pair<int,int>(2,3)] = 5 + n;
            edgeNodeMap[pair<int,int>(3,4)] = 5 + 2*n;
            edgeNodeMap[pair<int,int>(4,1)] = 5 + 3*n;

            // Add vertices
            for (int i = 0; i < 4; ++i) {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(  pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes)));
            }
        }


        unsigned int Tetrahedron::typeIds[5] = {
                GetElementFactory().RegisterCreatorFunction(4,  Tetrahedron::create, "Order 1 Tetrahedron"),
                GetElementFactory().RegisterCreatorFunction(11, Tetrahedron::create, "Order 2 Tetrahedron"),
                GetElementFactory().RegisterCreatorFunction(29, Tetrahedron::create, "Order 3 Tetrahedron"),
                GetElementFactory().RegisterCreatorFunction(30, Tetrahedron::create, "Order 4 Tetrahedron"),
                GetElementFactory().RegisterCreatorFunction(31, Tetrahedron::create, "Order 5 Tetrahedron")
            };

        Tetrahedron::Tetrahedron(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element()
        {
            m_tag = "A";
            m_dim = 3;
            m_taglist = pTagList;
            int n = 0;
            switch (m_taglist.back())
            {
            case 4:  n = 0; break;
            case 11: n = 1; break;
            case 29: n = 2; break;
            case 30: n = 3; break;
            case 31: n = 4; break;
            }
            int order = n+1;

            // Create a map to relate edge nodes to a pair of vertices defining an edge
            // This is based on the ordering produced by gmsh.
            map<pair<int,int>, int> edgeNodeMap;
            map<pair<int,int>, int>::iterator it;
            edgeNodeMap[pair<int,int>(1,2)] = 5;
            edgeNodeMap[pair<int,int>(2,3)] = 5 + n;
            edgeNodeMap[pair<int,int>(3,1)] = 5 + 2*n;
            edgeNodeMap[pair<int,int>(4,1)] = 5 + 3*n;
            edgeNodeMap[pair<int,int>(4,3)] = 5 + 4*n;
            edgeNodeMap[pair<int,int>(4,2)] = 5 + 5*n;

            // Add vertices
            for (int i = 0; i < 4; ++i) {
                vertex.push_back(pNodeList[i]);
            }

            // Create edges (with corresponding set of edge points)
            for (it = edgeNodeMap.begin(); it != edgeNodeMap.end(); ++it)
            {
                vector<NodeSharedPtr> edgeNodes;
                if (order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(  pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes)));
            }

            // Reorient the tet to ensure collapsed coordinates align between adjacent
            // elements.
            OrientTet();

            // Create faces
            int face_ids[4][3] = {
                {0,1,2},{0,1,3},{1,2,3},{2,0,3}};
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
                            faceEdges.push_back(edge[i]);
                            break;
                        }
                    }

                }
                int N = 4 + 6*n + j*n*(n-1)/2;
                for (unsigned int i = 0; i < n*(n-1)/2; ++i)
                {
                    faceNodes.push_back(pNodeList[N+i]);
                }
                face.push_back(FaceSharedPtr(new Face(faceVertices, faceNodes, faceEdges)));
            }
        }

        void Tetrahedron::OrientTet() {
            // Order vertices with lowest global vertex at top degenerate point
            // Place second lowest global vertex at base degenerate point
            sort(vertex.begin(), vertex.end());
            reverse(vertex.begin(), vertex.end());

            // Check orientation of tet and order remaining two points
            double ax, ay, az, vol, vax, vay, vaz, vbx, vby, vbz, vcx, vcy, vcz;
            vector<NodeSharedPtr> v;
            v.push_back(vertex[0]);
            v.push_back(vertex[1]);
            v.push_back(vertex[2]);
            v.push_back(vertex[3]);

            // Compute cross produc (b x c)
            vax = v[0]->x-v[2]->x;
            vay = v[0]->y-v[2]->y;
            vaz = v[0]->z-v[2]->z;
            vbx = v[1]->x-v[2]->x;
            vby = v[1]->y-v[2]->y;
            vbz = v[1]->z-v[2]->z;
            vcx = v[3]->x-v[2]->x;
            vcy = v[3]->y-v[2]->y;
            vcz = v[3]->z-v[2]->z;
            ax = vay*vbz - vaz*vby;
            ay = vaz*vbx - vbz*vax;
            az = vax*vby - vbx*vay;
            double dot = (ax*vcx + ay*vcy + az*vcz);

            // If negative, reverse order of non-degenerate points to correctly
            // orientate tet.
            if (dot < 0)
            {
                swap(vertex[0], vertex[1]);
                swap(v[0], v[1]);
            }
        }


        unsigned int Hexahedron::typeIds[2] = {
                GetElementFactory().RegisterCreatorFunction(5, Hexahedron::create, "Order 1 Hexahedron"),
                GetElementFactory().RegisterCreatorFunction(12, Hexahedron::create, "Order 2 Hexahedron")
            };

        Hexahedron::Hexahedron(vector<NodeSharedPtr> pNodeList, vector<int> pTagList)
                : Element() {
            m_tag = "H";
            m_dim = 3;
            m_taglist = pTagList;
            int n = 0;
            switch (m_taglist.back())
            {
            case 5:  n = 0; break;
            case 12: n = 1; break;
            }
            int order = n+1;

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
                if (order > 1) {
                    for (int j = it->second; j < it->second + n; ++j) {
                        edgeNodes.push_back(pNodeList[j-1]);
                    }
                }
                edge.push_back(EdgeSharedPtr(new Edge(  pNodeList[it->first.first-1],
                                                        pNodeList[it->first.second-1],
                                                        edgeNodes)));
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
                int N = 8 + 12*n + j*n*n;
                for (unsigned int i = 0; i < n*n; ++i)
                {
                    faceNodes.push_back(pNodeList[N+i]);
                }
                face.push_back(FaceSharedPtr(new Face(faceVertices, faceNodes, faceEdges)));
            }
        }
    }
}
