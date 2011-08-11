////////////////////////////////////////////////////////////////////////////////
//
//  File: ConvertSem.cpp
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
//  Description: Semtex session converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "MeshElements.h"
#include "ConvertSem.h"

namespace Nektar
{
    namespace Utilities
    {
        string ConvertSem::className = GetConvertFactory().RegisterCreatorFunction("sem", ConvertSem::create);

        ConvertSem::ConvertSem()
            : Convert()
        {

        }


        ConvertSem::ConvertSem(const ConvertSem& pSrc)
            : Convert(pSrc)
        {

        }


        ConvertSem::~ConvertSem()
        {

        }


        /**
         * Gmsh file contains a list of nodes and their coordinates, along with
         * a list of elements and those nodes which define them. We read in and
         * store the list of nodes in #m_node and store the list of elements in
         * #m_element. Each new element is supplied with a list of entries from
         * #m_node which defines the element. Finally some mesh statistics are
         * printed.
         *
         * @param   pFilename           Filename of Gmsh file to read.
         */
        void ConvertSem::ReadFile(const string pFilename)
        {
            m_expDim = 0;
            string line;
            int nVertices = 0;
            int nEntities = 0;
            int nElements = 0;
            int nBoundaryElements = 0;
            int elm_type = 0;
            vector<vector<NodeSharedPtr> > elementList;
            fstream mshFile(pFilename.c_str());

            if (!mshFile.is_open())
            {
                cout << "Unable to find semtex file" << endl;
                mshFile.close();
                return;
            }

            cout << "Start reading ConvertSem..." << endl;
            while (!mshFile.eof())
            {
                getline(mshFile, line);
                stringstream s(line);
                string word;
                s >> word;
                if (word == "<NODES")
                {
                    string tag = s.str();
                    int start = tag.find_first_of('=');
                    int end = tag.find_first_of('>');
                    nVertices = atoi(tag.substr(start+1,end).c_str());

                    // Read nodes
                    int id = 0;
                    int i = 0;
                    while (i < nVertices)
                    {
                        getline(mshFile, line);
                        if (line.length() < 7) continue;
                        stringstream st(line);
                        double x = 0, y = 0, z = 0;
                        st >> id >> x >> y >> z;

                        if ((y * y) > 0.000001 && m_spaceDim != 3)
                        {
                            m_spaceDim = 2;
                        }
                        if ((z * z) > 0.000001)
                        {
                            m_spaceDim = 3;
                        }
                        id -= 1; // counter starts at 0
                        m_node.push_back(boost::shared_ptr<Node>(new Node(id, x, y, z)));
                        ++i;
                    }
                }
                if (word == "<ELEMENTS")
                {
                    string tag = s.str();
                    int start = tag.find_first_of('=');
                    int end = tag.find_first_of('>');
                    nEntities = atoi(tag.substr(start+1,end).c_str());

                    // Read elements
                    int zeroDid = 0, oneDid = 0, twoDid = 0, threeDid = 0;
                    int i = 0;
                    while (i < nEntities)
                    {
                        getline(mshFile, line);
                        if (line.length() < 18)
                        {
                            continue;
                        }
                        stringstream st(line);
                        string tmp;
                        int id = 0, num_tag = 0, num_nodes = 0;
                        elm_type = 3;

                        // Create element tags
                        vector<int> tags;
                        tags.push_back(0); // composite
                        tags.push_back(3); // element type

                        // Read element node list
                        st >> id >> tmp;
                        vector<NodeSharedPtr> nodeList;
                        for (int k = 0; k < 4; ++k)
                        {
                            int node = 0;
                            st >> node;
                            node -= 1;
                            nodeList.push_back(m_node[node]);
                        }

                        elementList.push_back(nodeList);
                        // Create element
                        ElementSharedPtr E = GetElementFactory().CreateInstance(elm_type,nodeList,tags);

                        // Determine mesh expansion dimension
                        if (E->GetDim() > m_expDim) {
                            m_expDim = E->GetDim();
                        }
                        m_element.push_back(E);
                        ++i;
                    }

                    // Compute the number of full-dimensional elements and
                    // boundary elements.
                    for (int i = 0; i < m_element.size(); ++i) {
                        if (m_element[i]->GetDim() == m_expDim) {
                            nElements++;
                        }
                        if (m_element[i]->GetDim() == m_expDim - 1) {
                            nBoundaryElements++;
                        }
                    }
                }
                if (word == "<CURVES")
                {
                    string tag = s.str();
                    int start = tag.find_first_of('=');
                    int end = tag.find_first_of('>');
                    int nCurves = atoi(tag.substr(start+1,end).c_str());

                    int i = 0;
                    while (i < nCurves)
                    {
                        getline(mshFile, line);
                        if (line.length() < 18)
                        {
                            continue;
                        }
                        stringstream st(line);
                        string tmp;
                        int id = 0, elmt = 0, side = 0, radius = 0, num_tag = 0, num_nodes = 0;
                        int nodeId = nVertices;
                        st >> id >> elmt >> side >> tmp >> radius;
                        id--;
                        elmt--;

                        // First, let's make the element second order
                        if (elementList[elmt].size() == 4)
                        {
                            for (int j = 0; j < 4; ++j)
                            {
                                NodeSharedPtr n1 = elementList[elmt][j];
                                NodeSharedPtr n2 = elementList[elmt][(j+1)%4];
                                //double x = elementList
                                //m_node.push_back(boost::shared_ptr<Node>(new Node(nodeId++, x, y, z)));
                                //elementList[elmt].push_back();
                            }
                        }
                        ++i;
                    }

                    cout << "Expansion dimension is " << m_expDim << endl;
                    cout << "Space dimension is " << m_spaceDim << endl;
                    cout << "Read " << m_node.size() << " nodes" << endl;
                    cout << "Read " << m_element.size() << " geometric entities" << endl;
                    cout << "Read " << nElements << " " << m_expDim << "-D elements" << endl;
                    cout << "Read " << nBoundaryElements << " boundary entities" << endl;
                }
            }
            mshFile.close();
        }


        /**
         *
         */
        void ConvertSem::Process()
        {
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }


        /**
         * Each element is processed in turn and the vertices extracted. A
         * unique list of vertices is then compiled in #m_vertex. Finally, the
         * vertices are enumerated.
         */
        void ConvertSem::ProcessVertices()
        {
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetVertexCount(); ++j)
                    {
                        NodeSharedPtr x = m_element[i]->GetVertex(j);
                        if ( FindVertexIndex(x) < 0 )
                        {
                            m_vertex.push_back(x);
                        }
                    }
                }
            }
            // Enumerate vertices
            for (int i = 0; i < m_vertex.size(); ++i)
            {
                m_vertex[i]->id = i;
            }
        }


        /**
         * This routine only proceeds if the expansion dimension is 2 or 3.
         *
         * All elements are first scanned and a list of unique edges produced
         * in #m_edge which is then enumerated. Since each element generated
         * its edges independently, we must now ensure that each element only
         * uses edge objects from the #m_edge list. This ensures there are no
         * duplicate edge objects. Finally, we scan the list of elements for
         * 1-D boundary elements which correspond to an edge in #m_edge. For
         * such elements, we set its edgeLink to reference the corresponding
         * edge in #m_edge.
         */
        void ConvertSem::ProcessEdges()
        {
            if (m_expDim < 2) return;

            // Scan all elements and generate list of unique edges
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetEdgeCount(); ++j)
                    {
                        EdgeSharedPtr x = m_element[i]->GetEdge(j);
                        if ( FindEdgeIndex(x) < 0)
                        {
                            m_edge.push_back(x);
                        }
                    }
                }
            }
            // Enumerate edges
            for (int i = 0; i < m_edge.size(); ++i)
            {
                m_edge[i]->id = i;
            }
            // Remove duplicate edges
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetEdgeCount(); ++j)
                    {
                        int q = FindEdgeIndex(m_element[i]->GetEdge(j));
                        if (q < 0)
                        {
                            cout << "ERROR: should be able to find edge!" << endl;
                            abort();
                        }
                        m_element[i]->SetEdge(j,m_edge[q]);
                    }
                }
            }
            // Create links for 1D elements
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == 1)
                {
                    NodeSharedPtr v0 = m_element[i]->GetVertex(0);
                    NodeSharedPtr v1 = m_element[i]->GetVertex(1);
                    vector<NodeSharedPtr> edgeNodes;
                    EdgeSharedPtr E = boost::shared_ptr<Edge>(new Edge(v0, v1, edgeNodes));
                    int q = FindEdgeIndex(E);
                    if (q < 0)
                    {
                        cout << "Cannot find corresponding element face for 1D element " << i << endl;
                        abort();
                    }
                    m_element[i]->SetEdgeLink(m_edge[q]);
                }
            }
        }


        /**
         * This routine only proceeds if the expansion dimension is 3.
         *
         * All elements are scanned and a unique list of faces is produced in
         * #m_face, which is then enumerated. Since elements created their own
         * faces independently, we examine each element only uses face objects
         * from #m_face. Duplicate faces of those in #m_face are replaced with
         * the corresponding entry in #m_face. Finally, we scan the list of
         * elements for 2-D boundary faces which correspond to faces in
         * #m_face. For such elements, we set its faceLink to reference the
         * corresponding face in #m_face.
         */
        void ConvertSem::ProcessFaces()
        {
            if (m_expDim < 3) return;

            // Scan all elements and generate list of unique faces
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetFaceCount(); ++j)
                    {
                        FaceSharedPtr x = m_element[i]->GetFace(j);
                        if ( FindFaceIndex(x) < 0)
                        {
                            m_face.push_back(x);
                        }
                    }
                }
            }
            // Enumerate list of faces
            for (int i = 0; i < m_face.size(); ++i)
            {
                m_face[i]->id = i;
            }
            // Remove duplicate faces
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    for (int j = 0; j < m_element[i]->GetFaceCount(); ++j)
                    {
                        int q = FindFaceIndex(m_element[i]->GetFace(j));
                        if (q < 0)
                        {
                            cout << "ERROR: should be able to find face!" << endl;
                            abort();
                        }
                        m_element[i]->SetFace(j,m_face[q]);
                    }
                }
            }
            // Create links for 2D elements
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == 2)
                {
                    vector<NodeSharedPtr> vertices = m_element[i]->GetVertexList();
                    vector<NodeSharedPtr> faceNodes;
                    vector<EdgeSharedPtr> edgeList = m_element[i]->GetEdgeList();
                    FaceSharedPtr F = boost::shared_ptr<Face>(new Face(vertices, faceNodes, edgeList));
                    int q = FindFaceIndex(F);
                    if (q < 0)
                    {
                        cout << "Cannot find corresponding element face for 2D element " << i << endl;
                        abort();
                    }
                    m_element[i]->SetFaceLink(m_face[q]);
                }
            }
        }


        /**
         * For all elements of equal dimension to the mesh dimension, we
         * enumerate sequentially. All other elements in the list should be of
         * lower dimension and have ID set by a corresponding edgeLink or
         * faceLink (as set in #ProcessEdges or #ProcessFaces).
         */
        void ConvertSem::ProcessElements()
        {
            int cnt = 0;
            for (int i = 0; i < m_element.size(); ++i)
            {
                if (m_element[i]->GetDim() == m_expDim)
                {
                    m_element[i]->SetId(cnt++);
                }
            }
        }


        /**
         * Each element is assigned to a composite ID by Gmsh. First we scan
         * the element list and generate a list of composite IDs. We then
         * generate the composite objects and populate them with a second scan
         * through the element list.
         */
        void ConvertSem::ProcessComposites()
        {
            int p = 0;
            vector<ElementSharedPtr>::iterator it;
            list<int> compIdList;
            list<int>::iterator it2;

            // Generate sorted list of unique composite IDs to create.
            for (it = m_element.begin(); it != m_element.end(); ++it)
            {
                compIdList.push_back((*it)->GetTagList()[0]);
            }
            compIdList.sort();
            compIdList.unique();

            // Create a composite object for each unique composite ID.
            m_composite.resize(compIdList.size());
            it2 = compIdList.begin();
            for (int i = 0; i < m_composite.size(); ++i, ++it2)
            {
                m_composite[i] = boost::shared_ptr<Composite>(new Composite);
                m_composite[i]->id = *it2;
                m_composite[i]->tag = "";
            }

            // Populate composites with elements.
            for (int i = 0; i < m_element.size(); ++i)
            {
                int p = m_element[i]->GetTagList()[0];
                for (int j = 0; j < m_composite.size(); ++j)
                {
                    if (m_composite[j]->id == p)
                    {
                        if (m_composite[j]->tag == "")
                        {
                            // This is the first element in this composite
                            // so assign the composite tag.
                            m_composite[j]->tag = m_element[i]->GetTag();
                        }
                        else
                        {
                            // Otherwise, check element tag matches composite.
                            if (m_element[i]->GetTag() != m_composite[j]->tag)
                            {
                                cout << "Different types of elements in same composite!" << endl;
                                cout << " -> Composite uses " << m_composite[j]->tag << endl;
                                cout << " -> Element uses   " << m_element[i]->GetTag() << endl;
                                cout << "Have you specified physical volumes and surfaces?" << endl;
                            }
                        }
                        m_composite[j]->items.push_back(m_element[i]);
                        break;
                    }
                }
            }
        }

    }
}
