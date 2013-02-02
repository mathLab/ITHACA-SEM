////////////////////////////////////////////////////////////////////////////////
//
//  File: InputGmsh.cpp
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
#include <iostream>
using namespace std;

#include "MeshElements.h"
#include "InputGmsh.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputGmsh::className = 
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "msh"), InputGmsh::create,
                "Reads Gmsh msh file.");

        std::map<unsigned int, ElmtConfig> InputGmsh::elmMap = 
            InputGmsh::GenElmMap();


        /**
         * @brief Set up InputGmsh object.
         *
         */
        InputGmsh::InputGmsh(MeshSharedPtr m) : InputModule(m)
        {
            
        }

        InputGmsh::~InputGmsh()
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
        void InputGmsh::Process()
        {
            // Open the file stream.
            OpenStream();
            
            m->expDim = 0;
            m->spaceDim = 0;
            string line;
            int nVertices = 0;
            int nEntities = 0;
            int nElements = 0;
            int nBoundaryElements = 0;
            int elm_type = 0;
            int prevId = -1;
            map<unsigned int, ElmtConfig>::iterator it;

            if (m->verbose)
            {
                cout << "InputGmsh: Start reading file..." << endl;
            }

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

                        if ((x * x) > 0.000001 && m->spaceDim < 1)
                        {
                            m->spaceDim = 1;
                        }
                        if ((y * y) > 0.000001 && m->spaceDim < 2)
                        {
                            m->spaceDim = 2;
                        }
                        if ((z * z) > 0.000001 && m->spaceDim < 3)
                        {
                            m->spaceDim = 3;
                        }
                        
                        id -= 1; // counter starts at 0
                        
                        if (id-prevId != 1)
                        {
                            cerr << "Gmsh vertex ids should be contiguous" << endl;
                            abort();
                        }
                        prevId = id;
                        m->node.push_back(boost::shared_ptr<Node>(new Node(id, x, y, z)));
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

                        it = elmMap.find(elm_type);
                        if (it == elmMap.end())
                        {
                            cerr << "Error: element type " << elm_type
                                 << " not supported" << endl;
                            abort();
                        }

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
                            nodeList.push_back(m->node[node]);
                        }

                        // Prism nodes need re-ordering for Nektar++.
                        if (it->second.e == ePrism)
                        {
                            // Mirror first in uv plane to swap around
                            // triangular faces
                            swap(nodeList[0], nodeList[3]);
                            swap(nodeList[1], nodeList[4]);
                            swap(nodeList[2], nodeList[5]);
                            // Reorder base points so that face/vertices map
                            // correctly.
                            swap(nodeList[4], nodeList[2]);
                            
                            if (it->second.order == 2)
                            {
                                vector<NodeSharedPtr> nodemap(18);
                                
                                // Vertices remain unchanged.
                                nodemap[ 0] = nodeList[ 0];
                                nodemap[ 1] = nodeList[ 1];
                                nodemap[ 2] = nodeList[ 2];
                                nodemap[ 3] = nodeList[ 3];
                                nodemap[ 4] = nodeList[ 4];
                                nodemap[ 5] = nodeList[ 5];
                                // Reorder edge nodes: first mirror in uv
                                // plane and then place in Nektar++ ordering.
                                nodemap[ 6] = nodeList[12];
                                nodemap[ 7] = nodeList[10];
                                nodemap[ 8] = nodeList[ 6];
                                nodemap[ 9] = nodeList[ 8];
                                nodemap[10] = nodeList[13];
                                nodemap[11] = nodeList[14];
                                nodemap[12] = nodeList[ 9];
                                nodemap[13] = nodeList[ 7];
                                nodemap[14] = nodeList[11];
                                // Face vertices remain unchanged.
                                nodemap[15] = nodeList[15];
                                nodemap[16] = nodeList[16];
                                nodemap[17] = nodeList[17];
                                
                                nodeList = nodemap;
                            }
                            else if (it->second.order > 2)
                            {
                                cerr << "Error: gmsh prisms only supported up "
                                     << "to second order." << endl;
                                abort();
                            }
                        }
                        
                        // Create element
                        ElementSharedPtr E = GetElementFactory().
                            CreateInstance(it->second.e,it->second,nodeList,tags);

                        // Determine mesh expansion dimension
                        if (E->GetDim() > m->expDim) {
                            m->expDim = E->GetDim();
                        }
                        m->element[E->GetDim()].push_back(E);
                    }
                }
            }
            mshFile.close();

            // Process rest of mesh.
            ProcessVertices  ();
            ProcessEdges     ();
            ProcessFaces     ();
            ProcessElements  ();
            ProcessComposites();
        }

        /**
         * For a given msh ID, return the corresponding number of nodes.
         */
        int InputGmsh::GetNnodes(unsigned int InputGmshEntity)
        {
            int nNodes;
            map<unsigned int, ElmtConfig>::iterator it;
            
            it = elmMap.find(InputGmshEntity);
            
            if (it == elmMap.end())
            {
                cerr << "Unknown element type " << InputGmshEntity << endl;
                abort();
            }
            
            switch(it->second.e)
            {
                case ePoint: 
                    nNodes = Point::        GetNumNodes(it->second);
                    break;
                case eLine: 
                    nNodes = Line::         GetNumNodes(it->second);
                    break;
                case eTriangle: 
                    nNodes = Triangle::     GetNumNodes(it->second);
                    break;
                case eQuadrilateral: 
                    nNodes = Quadrilateral::GetNumNodes(it->second);
                    break;
                case eTetrahedron: 
                    nNodes = Tetrahedron::  GetNumNodes(it->second);
                    break;
                case ePrism: 
                    nNodes = Prism::        GetNumNodes(it->second);
                    break;
                case eHexahedron:
                    nNodes = Hexahedron::   GetNumNodes(it->second);
                    break;
                default:
                    cerr << "Unknown element type!" << endl;
                    abort();
                    break;
            }
            
            return nNodes;
        }
        
        /*
         * @brief Populate the element map #elmMap.
         * 
         * This function primarily populates the element mapping #elmMap,
         * which takes a msh ID used by Gmsh and translates to element type,
         * element order and whether the element is incomplete (i.e. whether
         * it contains solely boundary nodes, or just face nodes). Note that
         * some of these elements, such as prisms of order >= 3, cannot yet be
         * generated by Gmsh.
         */
        std::map<unsigned int, ElmtConfig> InputGmsh::GenElmMap()
        {
            std::map<unsigned int, ElmtConfig> tmp;
            
            //                    Elmt type,   order,  face, volume
            tmp[  1] = ElmtConfig(eLine,           1,  true,  true);
            tmp[  2] = ElmtConfig(eTriangle,       1,  true,  true);
            tmp[  3] = ElmtConfig(eQuadrilateral,  1,  true,  true);
            tmp[  4] = ElmtConfig(eTetrahedron,    1,  true,  true);
            tmp[  5] = ElmtConfig(eHexahedron,     1,  true,  true);
            tmp[  6] = ElmtConfig(ePrism,          1,  true,  true);
            tmp[  8] = ElmtConfig(eLine,           2,  true,  true);
            tmp[  9] = ElmtConfig(eTriangle,       2,  true,  true);
            tmp[ 10] = ElmtConfig(eQuadrilateral,  2,  true,  true);
            tmp[ 11] = ElmtConfig(eTetrahedron,    2,  true,  true);
            tmp[ 12] = ElmtConfig(eHexahedron,     2,  true,  true);
            tmp[ 13] = ElmtConfig(ePrism,          2,  true,  true);
            tmp[ 15] = ElmtConfig(ePoint,          1,  true, false);
            tmp[ 16] = ElmtConfig(eQuadrilateral,  2, false, false);
            tmp[ 17] = ElmtConfig(eHexahedron,     2, false, false);
            tmp[ 18] = ElmtConfig(ePrism,          2, false, false);
            tmp[ 20] = ElmtConfig(eTriangle,       3, false, false);
            tmp[ 21] = ElmtConfig(eTriangle,       3,  true, false);
            tmp[ 22] = ElmtConfig(eTriangle,       4, false, false);
            tmp[ 23] = ElmtConfig(eTriangle,       4,  true, false);
            tmp[ 24] = ElmtConfig(eTriangle,       5, false, false);
            tmp[ 25] = ElmtConfig(eTriangle,       5,  true, false);
            tmp[ 26] = ElmtConfig(eLine,           3,  true, false);
            tmp[ 27] = ElmtConfig(eLine,           4,  true, false);
            tmp[ 28] = ElmtConfig(eLine,           5,  true, false);
            tmp[ 29] = ElmtConfig(eTetrahedron,    3,  true,  true);
            tmp[ 30] = ElmtConfig(eTetrahedron,    4,  true,  true);
            tmp[ 31] = ElmtConfig(eTetrahedron,    5,  true,  true);
            tmp[ 32] = ElmtConfig(eTetrahedron,    4,  true, false);
            tmp[ 33] = ElmtConfig(eTetrahedron,    5,  true, false);
            tmp[ 36] = ElmtConfig(eQuadrilateral,  3,  true, false);
            tmp[ 37] = ElmtConfig(eQuadrilateral,  4,  true, false);
            tmp[ 38] = ElmtConfig(eQuadrilateral,  5,  true, false);
            tmp[ 39] = ElmtConfig(eQuadrilateral,  3, false, false);
            tmp[ 40] = ElmtConfig(eQuadrilateral,  4, false, false);
            tmp[ 41] = ElmtConfig(eQuadrilateral,  5, false, false);
            tmp[ 42] = ElmtConfig(eTriangle,       6,  true, false);
            tmp[ 43] = ElmtConfig(eTriangle,       7,  true, false);
            tmp[ 44] = ElmtConfig(eTriangle,       8,  true, false);
            tmp[ 45] = ElmtConfig(eTriangle,       9,  true, false);
            tmp[ 46] = ElmtConfig(eTriangle,      10,  true, false);
            tmp[ 47] = ElmtConfig(eQuadrilateral,  6,  true, false);
            tmp[ 48] = ElmtConfig(eQuadrilateral,  7,  true, false);
            tmp[ 49] = ElmtConfig(eQuadrilateral,  8,  true, false);
            tmp[ 50] = ElmtConfig(eQuadrilateral,  9,  true, false);
            tmp[ 51] = ElmtConfig(eQuadrilateral, 10,  true, false);
            tmp[ 52] = ElmtConfig(eTriangle,       6, false, false);
            tmp[ 53] = ElmtConfig(eTriangle,       7, false, false);
            tmp[ 54] = ElmtConfig(eTriangle,       8, false, false);
            tmp[ 55] = ElmtConfig(eTriangle,       9, false, false);
            tmp[ 56] = ElmtConfig(eTriangle,      10, false, false);
            tmp[ 57] = ElmtConfig(eQuadrilateral,  6, false, false);
            tmp[ 58] = ElmtConfig(eQuadrilateral,  7, false, false);
            tmp[ 59] = ElmtConfig(eQuadrilateral,  8, false, false);
            tmp[ 60] = ElmtConfig(eQuadrilateral,  9, false, false);
            tmp[ 61] = ElmtConfig(eQuadrilateral, 10, false, false);
            tmp[ 62] = ElmtConfig(eLine,           6,  true, false);
            tmp[ 63] = ElmtConfig(eLine,           7,  true, false);
            tmp[ 64] = ElmtConfig(eLine,           8,  true, false);
            tmp[ 65] = ElmtConfig(eLine,           9,  true, false);
            tmp[ 66] = ElmtConfig(eLine,          10,  true, false);
            tmp[ 71] = ElmtConfig(eTetrahedron,    6,  true,  true);
            tmp[ 72] = ElmtConfig(eTetrahedron,    7,  true,  true);
            tmp[ 73] = ElmtConfig(eTetrahedron,    8,  true,  true);
            tmp[ 74] = ElmtConfig(eTetrahedron,    9,  true,  true);
            tmp[ 75] = ElmtConfig(eTetrahedron,   10,  true,  true);
            tmp[ 79] = ElmtConfig(eTetrahedron,    6,  true, false);
            tmp[ 80] = ElmtConfig(eTetrahedron,    7,  true, false);
            tmp[ 81] = ElmtConfig(eTetrahedron,    8,  true, false);
            tmp[ 82] = ElmtConfig(eTetrahedron,    9,  true, false);
            tmp[ 83] = ElmtConfig(eTetrahedron,   10,  true, false);
            tmp[ 90] = ElmtConfig(ePrism,          3,  true,  true);
            tmp[ 91] = ElmtConfig(ePrism,          4,  true,  true);
            tmp[ 92] = ElmtConfig(eHexahedron,     3,  true,  true);
            tmp[ 93] = ElmtConfig(eHexahedron,     4,  true,  true);
            tmp[ 94] = ElmtConfig(eHexahedron,     5,  true,  true);
            tmp[ 95] = ElmtConfig(eHexahedron,     6,  true,  true);
            tmp[ 96] = ElmtConfig(eHexahedron,     7,  true,  true);
            tmp[ 97] = ElmtConfig(eHexahedron,     8,  true,  true);
            tmp[ 98] = ElmtConfig(eHexahedron,     9,  true,  true);
            tmp[ 99] = ElmtConfig(eHexahedron,     3,  true, false);
            tmp[100] = ElmtConfig(eHexahedron,     4,  true, false);
            tmp[101] = ElmtConfig(eHexahedron,     5,  true, false);
            tmp[102] = ElmtConfig(eHexahedron,     6,  true, false);
            tmp[103] = ElmtConfig(eHexahedron,     7,  true, false);
            tmp[104] = ElmtConfig(eHexahedron,     8,  true, false);
            tmp[105] = ElmtConfig(eHexahedron,     9,  true, false);
            tmp[106] = ElmtConfig(ePrism,          5,  true,  true);
            tmp[107] = ElmtConfig(ePrism,          6,  true,  true);
            tmp[108] = ElmtConfig(ePrism,          7,  true,  true);
            tmp[109] = ElmtConfig(ePrism,          8,  true,  true);
            tmp[110] = ElmtConfig(ePrism,          9,  true,  true);
            tmp[111] = ElmtConfig(ePrism,          3,  true, false);
            tmp[112] = ElmtConfig(ePrism,          4,  true, false);
            tmp[113] = ElmtConfig(ePrism,          5,  true, false);
            tmp[114] = ElmtConfig(ePrism,          6,  true, false);
            tmp[115] = ElmtConfig(ePrism,          7,  true, false);
            tmp[116] = ElmtConfig(ePrism,          8,  true, false);
            tmp[117] = ElmtConfig(ePrism,          9,  true, false);
            
            return tmp;
        }
    }
}
