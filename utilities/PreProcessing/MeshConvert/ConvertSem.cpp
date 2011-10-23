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
    }
}
