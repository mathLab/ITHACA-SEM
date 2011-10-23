////////////////////////////////////////////////////////////////////////////////
//
//  File: ConvertGmsh.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "MeshElements.h"
#include "ConvertGmsh.h"

namespace Nektar
{
    namespace Utilities
    {
        string ConvertGmsh::className = GetConvertFactory().RegisterCreatorFunction("msh", ConvertGmsh::create);

        ConvertGmsh::ConvertGmsh()
            : Convert()
        {

        }


        ConvertGmsh::ConvertGmsh(const ConvertGmsh& pSrc)
            : Convert(pSrc)
        {

        }


        ConvertGmsh::~ConvertGmsh()
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
        void ConvertGmsh::ReadFile(const string pFilename)
        {
            m_expDim = 0;
            m_spaceDim = 0;
            string line;
            int nVertices = 0;
            int nEntities = 0;
            int nElements = 0;
            int nBoundaryElements = 0;
            int elm_type = 0;

            fstream mshFile(pFilename.c_str());

            if (!mshFile.is_open())
            {
                cout << "Unable to find msh file" << endl;
                mshFile.close();
                return;
            }

            cout << "Start reading ConvertGmsh..." << endl;
            while (!mshFile.eof())
            {
                getline(mshFile, line);
                stringstream s(line);
                string word;
                s >> word;

                // Process nodes.
                if (word == "$Nodes")
                {
                    getline(mshFile, line);
                    stringstream s(line);
                    s >> nVertices;
                    int id = 0;
                    for (int i = 0; i < nVertices; ++i)
                    {
                        getline(mshFile, line);
                        stringstream st(line);
                        double x = 0, y = 0, z = 0;
                        st >> id >> x >> y >> z;

                        if ((x * x) > 0.000001 && m_spaceDim < 1)
                        {
                            m_spaceDim = 1;
                        }
                        if ((y * y) > 0.000001 && m_spaceDim < 2)
                        {
                            m_spaceDim = 2;
                        }
                        if ((z * z) > 0.000001 && m_spaceDim < 3)
                        {
                            m_spaceDim = 3;
                        }
                        id -= 1; // counter starts at 0
                        m_node.push_back(boost::shared_ptr<Node>(new Node(id, x, y, z)));
                    }
                }
                // Process elements
                else if (word == "$Elements")
                {
                    int zeroDid = 0, oneDid = 0, twoDid = 0, threeDid = 0;
                    getline(mshFile, line);
                    stringstream s(line);
                    s >> nEntities;
                    for (int i = 0; i < nEntities; ++i)
                    {
                        getline(mshFile, line);
                        stringstream st(line);
                        int id = 0, num_tag = 0, num_nodes = 0;

                        st >> id >> elm_type >> num_tag;
                        id -= 1; // counter starts at 0

                        // Read element tags
                        vector<int> tags;
                        for (int j = 0; j < num_tag; ++j)
                        {
                            int tag = 0;
                            st >> tag;
                            tags.push_back(tag);
                        }
                        tags.push_back(elm_type);

                        // Read element node list
                        vector<NodeSharedPtr> nodeList;
                        num_nodes = GetNnodes(elm_type);
                        for (int k = 0; k < num_nodes; ++k)
                        {
                            int node = 0;
                            st >> node;
                            node -= 1; // counter starts at 0
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
        void ConvertGmsh::Process()
        {
            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
        }


        /**
         * This routine aids the reading of Gmsh files only.
         */
        int ConvertGmsh::GetNnodes(int ConvertGmshEntity)
        {
            int nNodes;

            switch(ConvertGmshEntity)
            {
            case 1:  nNodes = 2;  break;
            case 2:  nNodes = 3;  break;
            case 3:  nNodes = 4;  break;
            case 4:  nNodes = 4;  break;
            case 5:  nNodes = 8;  break;
            case 6:  nNodes = 6;  break;
            case 7:  nNodes = 5;  break;
            case 8:  nNodes = 3;  break;
            case 9:  nNodes = 6;  break;
            case 10: nNodes = 9;  break;
            case 11: nNodes = 10; break;
            case 12: nNodes = 27; break;
            case 13: nNodes = 18; break;
            case 14: nNodes = 14; break;
            case 15: nNodes = 1;  break;
            case 16: nNodes = 8;  break;
            case 17: nNodes = 20; break;
            case 18: nNodes = 15; break;
            case 19: nNodes = 13; break;
            case 20: nNodes = 9;  break;
            case 21: nNodes = 10; break;
            case 22: nNodes = 12; break;
            case 23: nNodes = 15; break;
            case 24: nNodes = 15; break;
            case 25: nNodes = 21; break;
            case 26: nNodes = 4;  break;
            case 27: nNodes = 5;  break;
            case 28: nNodes = 6;  break;
            case 29: nNodes = 20; break;
            case 30: nNodes = 35; break;
            case 31: nNodes = 56; break;
            default:
              cout << "unknown Gmsh element type" << endl;
            }

            return nNodes;
        }
    }
}
