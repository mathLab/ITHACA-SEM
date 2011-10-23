////////////////////////////////////////////////////////////////////////////////
//
//  File: ConvertPly.cpp
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
//  Description: PLY converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "MeshElements.h"
#include "ConvertPly.h"

namespace Nektar
{
    namespace Utilities
    {
        string ConvertPly::className = GetConvertFactory().RegisterCreatorFunction("ply", ConvertPly::create);

        ConvertPly::ConvertPly()
            : Convert()
        {

        }


        ConvertPly::ConvertPly(const ConvertPly& pSrc)
            : Convert(pSrc)
        {

        }


        ConvertPly::~ConvertPly()
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
        void ConvertPly::ReadFile(const string pFilename)
        {
            m_expDim = 0;
            string line;
            int nVertices = 0;
            int nEntities = 0;
            int nElements = 0;
            int nBoundaryElements = 0;
            int elm_type = 0;

            fstream mshFile(pFilename.c_str());

            if (!mshFile.is_open())
            {
                cout << "Unable to find ply file" << endl;
                mshFile.close();
                return;
            }

            cout << "Start reading ConvertPly..." << endl;
            while (!mshFile.eof())
            {
                getline(mshFile, line);
                stringstream s(line);
                string word;
                s >> word;
                if (word == "element")
                {
                    s >> word;
                    if (word == "vertex")
                    {
                        s >> nVertices;
                    }
                    else if (word == "face")
                    {
                        s >> nEntities;
                    }
                    continue;
                }
                else if (word == "end_header")
                {
                    // Read nodes
                    int id = 0;
                    for (int i = 0; i < nVertices; ++i)
                    {
                        getline(mshFile, line);
                        stringstream st(line);
                        double x = 0, y = 0, z = 0;
                        st >> x >> y >> z;

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
                        id++;
                    }

                    // Read elements
                    int zeroDid = 0, oneDid = 0, twoDid = 0, threeDid = 0;
                    for (int i = 0; i < nEntities; ++i)
                    {
                        getline(mshFile, line);
                        stringstream st(line);
                        int id = 0, num_tag = 0, num_nodes = 0;
                        elm_type = 2;

                        // Create element tags
                        vector<int> tags;
                        tags.push_back(0); // composite
                        tags.push_back(2); // element type

                        // Read element node list
                        st >> id;
                        vector<NodeSharedPtr> nodeList;
                        for (int k = 0; k < 3; ++k)
                        {
                            int node = 0;
                            st >> node;
                            nodeList.push_back(m_node[node]);
                        }

                        // Create element
                        ElementSharedPtr E = GetElementFactory().CreateInstance(elm_type,nodeList,tags);

                        // Determine mesh expansion dimension
                        if (E->GetDim() > m_expDim) {
                            m_expDim = E->GetDim();
                        }
                        m_element.push_back(E);
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
        void ConvertPly::Process()
        {
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }
    }
}
