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

        Element::Element(ElmtConfig   pConf,
                         unsigned int pNumNodes,
                         unsigned int pGotNodes) 
            : m_conf(pConf)
        {
            if (pNumNodes != pGotNodes)
            {
                cerr << "Number of modes mismatch for type " << pConf.e << "! Should be " << pNumNodes 
                     << " but got " << pGotNodes << " nodes." << endl;
                abort();
            }
        }

        Mesh::Mesh(const string pInFilename, const string pOutFilename) :
            inFilename(pInFilename), outFilename(pOutFilename)
        {
            
        }

        unsigned int Mesh::GetNumElements()
        {
            return element[expDim].size();
        }

        unsigned int Mesh::GetNumBndryElements()
        {
            int i, nElmt = 0;
            
            for (i = 0; i < expDim; ++i)
                nElmt += element[i].size();
            
            return nElmt;
        }
        
        unsigned int Mesh::GetNumEntities()
        {
            int nEnt = 0;
            
            for (int d = 0; d < expDim; ++d)
            {
                nEnt += element[d].size();
            }
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
                    if (face[i]->edgeList[j] == vOld)
                    {
                        face[i]->edgeList[j] = pNew;
                    }
                }
            }
        }

        bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2)
        {
            return p1->x == p2->x && p1->y == p2->y && p1->z == p2->z;
        }

        bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2)
        {
            return p1->id < p2->id;
        }

        bool operator==(EdgeSharedPtr const &p1, EdgeSharedPtr const &p2)
        {
            return ( ((*(p1->n1) == *(p2->n1)) && (*(p1->n2) == *(p2->n2)))
                  || ((*(p1->n2) == *(p2->n1)) && (*(p1->n1) == *(p2->n2))));
        }

        bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2)
        {
            bool e = true;
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

        ElementType Point::type = GetElementFactory().
            RegisterCreatorFunction(ePoint, Point::create, "Point");
        
        Point::Point(ElmtConfig            pConf,
                     vector<NodeSharedPtr> pNodeList, 
                     vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag = "";
            m_dim = 0;
            m_taglist = pTagList;
            vertex.push_back(pNodeList[0]);
        }

        unsigned int Point::GetNumNodes(ElmtConfig pConf)
        {
            return 1;
        }


        ElementType Line::type = GetElementFactory().
            RegisterCreatorFunction(eLine, Line::create, "Line");
        
        Line::Line(ElmtConfig            pConf,
                   vector<NodeSharedPtr> pNodeList, 
                   vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag = "S";
            m_dim = 1;
            m_taglist = pTagList;
            int n = m_conf.order-1;
            
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
        
        unsigned int Line::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            return pConf.order+1;
        }


        ElementType Triangle::type = GetElementFactory().
            RegisterCreatorFunction(eTriangle, Triangle::create, "Triangle");
        
        Triangle::Triangle(ElmtConfig            pConf,
                           vector<NodeSharedPtr> pNodeList, 
                           vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size()) 
        {
            m_tag = "T";
            m_dim = 2;
            m_taglist = pTagList;
            int n = m_conf.order-1;

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
        }

        unsigned int Triangle::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (!pConf.faceNodes)
                return (n+1)+2*(n-1)+1;
            else
                return (n+1)*(n+2)/2;
        }


        ElementType Quadrilateral::type = GetElementFactory().
            RegisterCreatorFunction(eQuadrilateral, Quadrilateral::create, 
                                    "Quadrilateral");
        
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
            for (int i = 0; i < 4; ++i) {
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
        }

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
            // defining an edge. This is based on the ordering produced by
            // gmsh.
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

                if (m_conf.faceNodes)
                {
                    int N = 4 + 6*n + j*n*(n-1)/2;
                    for (unsigned int i = 0; i < n*(n-1)/2; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }
        }
        
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

        void Tetrahedron::OrientTet() {
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
        }


        ElementType Prism::type = GetElementFactory().
            RegisterCreatorFunction(ePrism, Prism::create, "Prism");
        
        Prism::Prism(ElmtConfig            pConf,
                     vector<NodeSharedPtr> pNodeList,
                     vector<int>           pTagList)
            : Element(pConf, GetNumNodes(pConf), pNodeList.size())
        {
            m_tag = "R";
            m_dim = 3;
            m_taglist = pTagList;
            int n = m_conf.order-1;

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
            for (int i = 0; i < 6; ++i) {
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
            
            OrientPrism();
            
            // Create faces
            int face_ids[5][4] = {
                {0,1,2,3},{0,1,4,-1},{1,2,5,4},{3,2,5,-1},{0,3,5,4}};
            int faceoffset = 0;
            for (int j = 0; j < 5; ++j)
            {
                vector<NodeSharedPtr> faceVertices;
                vector<EdgeSharedPtr> faceEdges;
                vector<NodeSharedPtr> faceNodes;
                for (int k = 0; k < 4; ++k)
                {
                    if (face_ids[j][k] == -1)
                        break;
                    faceVertices.push_back(vertex[face_ids[j][k]]);
                    NodeSharedPtr a = vertex[face_ids[j][k]];
                    NodeSharedPtr b = vertex[face_ids[j][(k+1)%(j%2==0?4:3)]];
                    for (unsigned int i = 0; i < edge.size(); ++i)
                    {
                        if ((edge[i]->n1 == a && edge[i]->n2 == b) || 
                            (edge[i]->n1 == b && edge[i]->n2 == a))
                        {
                            faceEdges.push_back(edge[i]);
                            break;
                        }
                    }
                }
                
                if (m_conf.faceNodes)
                {
                    int facenodes = j%2==0 ? n*n : n*(n-1)/2;
                    faceoffset   += facenodes;
                    int N = 6 + 9*n + faceoffset;
                    for (unsigned int i = 0; i < facenodes; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }
        }

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

        void Prism::OrientPrism()
        {
            int lid[6], gid[6];
            
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
                    for (unsigned int i = 0; i < n*n; ++i)
                    {
                        faceNodes.push_back(pNodeList[N+i]);
                    }
                }
                face.push_back(FaceSharedPtr(
                    new Face(faceVertices, faceNodes, faceEdges, m_conf.faceCurveType)));
            }
        }
        
        unsigned int Hexahedron::GetNumNodes(ElmtConfig pConf)
        {
            int n = pConf.order;
            if (pConf.faceNodes && pConf.volumeNodes)
                return (n+1)*(n+1)*(n+1);
            else if (pConf.faceNodes && !pConf.volumeNodes)
                return 6*(n+1)*(n+1)-12*(n+1)+8;
            else
                return 12*(n+1)-8;
        }
    }
}
