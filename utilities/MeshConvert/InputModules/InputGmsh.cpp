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

#include "../MeshElements.h"
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
         * @brief Reorder a quadrilateral to appear in Nektar++ ordering from
         * Gmsh.
         */
        std::vector<int> quadTensorNodeOrdering(
            const std::vector<int> &nodes, int n)
        {
            std::vector<int> nodeList;

            // Triangle
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
                std::vector<int> interior((n-2)*(n-2));
                std::copy(nodes.begin() + 4+4*(n-2), nodes.end(), interior.begin());
                interior = quadTensorNodeOrdering(interior, n-2);

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

            return nodeList;
        }

        std::vector<int> triTensorNodeOrdering(
            const std::vector<int> &nodes, int n)
        {
            std::vector<int> nodeList;
            int cnt2;

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
                nodeList[i]         = nodes[3+i-1];
                nodeList[cnt]       = nodes[3+3*(n-2)-i];
                nodeList[cnt+n-i-1] = nodes[3+(n-2)+i-1];
                cnt += n-i;
            }

            // Interior (recursion)
            if (n > 3)
            {
                // Reorder interior nodes
                std::vector<int> interior((n-3)*(n-2)/2);
                std::copy(nodes.begin() + 3+3*(n-2), nodes.end(), interior.begin());
                interior = triTensorNodeOrdering(interior, n-3);

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

            return nodeList;
        }

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

            m_mesh->m_expDim = 0;
            m_mesh->m_spaceDim = 0;
            string line;
            int nVertices = 0;
            int nEntities = 0;
            int elm_type = 0;
            int prevId = -1;
            int maxTagId = -1;

            map<unsigned int, ElmtConfig>::iterator it;

            // This map takes each element ID and maps it to a permutation map
            // that is required to take Gmsh element node orderings and map them
            // to Nektar++ orderings.
            boost::unordered_map<int, vector<int> > orderingMap;
            boost::unordered_map<int, vector<int> >::iterator oIt;

            if (m_mesh->m_verbose)
            {
                cout << "InputGmsh: Start reading file..." << endl;
            }

            while (!m_mshFile.eof())
            {
                getline(m_mshFile, line);
                stringstream s(line);
                string word;
                s >> word;

                // Process nodes.
                if (word == "$Nodes")
                {
                    getline(m_mshFile, line);
                    stringstream s(line);
                    s >> nVertices;
                    int id = 0;
                    for (int i = 0; i < nVertices; ++i)
                    {
                        getline(m_mshFile, line);
                        stringstream st(line);
                        double x = 0, y = 0, z = 0;
                        st >> id >> x >> y >> z;

                        if ((x * x) > 0.000001 && m_mesh->m_spaceDim < 1)
                        {
                            m_mesh->m_spaceDim = 1;
                        }
                        if ((y * y) > 0.000001 && m_mesh->m_spaceDim < 2)
                        {
                            m_mesh->m_spaceDim = 2;
                        }
                        if ((z * z) > 0.000001 && m_mesh->m_spaceDim < 3)
                        {
                            m_mesh->m_spaceDim = 3;
                        }

                        id -= 1; // counter starts at 0

                        if (id-prevId != 1)
                        {
                            cerr << "Gmsh vertex ids should be contiguous" << endl;
                            abort();
                        }
                        prevId = id;
                        m_mesh->m_node.push_back(boost::shared_ptr<Node>(new Node(id, x, y, z)));
                    }
                }
                // Process elements
                else if (word == "$Elements")
                {
                    getline(m_mshFile, line);
                    stringstream s(line);
                    s >> nEntities;
                    for (int i = 0; i < nEntities; ++i)
                    {
                        getline(m_mshFile, line);
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
                        tags.resize(1);

                        maxTagId = max(maxTagId, tags[0]);

                        // Read element node list
                        vector<NodeSharedPtr> nodeList;
                        num_nodes = GetNnodes(elm_type);
                        for (int k = 0; k < num_nodes; ++k)
                        {
                            int node = 0;
                            st >> node;
                            node -= 1; // counter starts at 0
                            nodeList.push_back(m_mesh->m_node[node]);
                        }

                        // Look up reordering.
                        oIt = orderingMap.find(elm_type);

                        // If it's not created, then create it.
                        if (oIt == orderingMap.end())
                        {
                            oIt = orderingMap.insert(
                                make_pair(elm_type, CreateReordering(elm_type)))
                                .first;
                        }

                        // Apply reordering map where necessary.
                        if (oIt->second.size() > 0)
                        {
                            vector<int> &mapping = oIt->second;
                            vector<NodeSharedPtr> tmp = nodeList;
                            for (int i = 0; i < mapping.size(); ++i)
                            {
                                nodeList[i] = tmp[mapping[i]];
                            }
                        }

                        // Create element
                        ElementSharedPtr E = GetElementFactory().
                            CreateInstance(it->second.m_e,it->second,nodeList,tags);

                        // Determine mesh expansion dimension
                        if (E->GetDim() > m_mesh->m_expDim) {
                            m_mesh->m_expDim = E->GetDim();
                        }
                        m_mesh->m_element[E->GetDim()].push_back(E);
                    }
                }
            }
            m_mshFile.close();

            // Go through element and remap tags if necessary.
            map<int, map<LibUtilities::ShapeType, int> > compMap;
            map<int, map<LibUtilities::ShapeType, int> >::iterator cIt;
            map<LibUtilities::ShapeType, int>::iterator sIt;

            for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
            {
                ElementSharedPtr el = m_mesh->m_element[m_mesh->m_expDim][i];
                LibUtilities::ShapeType type = el->GetConf().m_e;

                vector<int> tags = el->GetTagList();
                int tag = tags[0];

                cIt = compMap.find(tag);

                if (cIt == compMap.end())
                {
                    compMap[tag][type] = tag;
                    continue;
                }

                // Reset tag for this element.
                sIt = cIt->second.find(type);
                if (sIt == cIt->second.end())
                {
                    maxTagId++;
                    cIt->second[type] = maxTagId;
                    tags[0] = maxTagId;
                    el->SetTagList(tags);
                }
                else if (sIt->second != tag)
                {
                    tags[0] = sIt->second;
                    el->SetTagList(tags);
                }
            }

            bool printInfo = false;
            for (cIt = compMap.begin(); cIt != compMap.end(); ++cIt)
            {
                if (cIt->second.size() > 1)
                {
                    printInfo = true;
                    break;
                }
            }

            if (printInfo)
            {
                cout << "Multiple elements in composite detected; remapped:"
                     << endl;
                for (cIt = compMap.begin(); cIt != compMap.end(); ++cIt)
                {
                    if (cIt->second.size() > 1)
                    {
                        sIt = cIt->second.begin();
                        cout << "- Tag " << cIt->first << " => " << sIt->second
                             << " ("
                             << LibUtilities::ShapeTypeMap[sIt->first] << ")";
                        sIt++;

                        for (; sIt != cIt->second.end(); ++sIt)
                        {
                            cout << ", " << sIt->second << " ("
                                 << LibUtilities::ShapeTypeMap[sIt->first]
                                 << ")";
                        }

                        cout << endl;
                    }
                }
            }

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

            switch(it->second.m_e)
            {
                case LibUtilities::ePoint:
                    nNodes = Point::        GetNumNodes(it->second);
                    break;
                case LibUtilities::eSegment:
                    nNodes = Line::         GetNumNodes(it->second);
                    break;
                case LibUtilities::eTriangle:
                    nNodes = Triangle::     GetNumNodes(it->second);
                    break;
                case LibUtilities::eQuadrilateral:
                    nNodes = Quadrilateral::GetNumNodes(it->second);
                    break;
                case LibUtilities::eTetrahedron:
                    nNodes = Tetrahedron::  GetNumNodes(it->second);
                    it->second.m_faceCurveType =
                        LibUtilities::eNodalTriEvenlySpaced;
                    break;
                case LibUtilities::ePyramid:
                    nNodes = Pyramid::      GetNumNodes(it->second);
                    break;
                case LibUtilities::ePrism:
                    nNodes = Prism::        GetNumNodes(it->second);
                    break;
                case LibUtilities::eHexahedron:
                    nNodes = Hexahedron::   GetNumNodes(it->second);
                    break;
                default:
                    cerr << "Unknown element type!" << endl;
                    abort();
                    break;
            }

            return nNodes;
        }

        /**
         * @brief Create a reordering map for a given element.
         *
         * Since Gmsh and Nektar++ have different vertex, edge and face
         * orientations, we need to reorder the nodes in a Gmsh MSH file so that
         * they work with the Nektar++ orderings, since this is what is used in
         * the elements defined in the converter.
         */
        vector<int> InputGmsh::CreateReordering(unsigned int InputGmshEntity)
        {
            map<unsigned int, ElmtConfig>::iterator it;

            it = elmMap.find(InputGmshEntity);

            if (it == elmMap.end())
            {
                cerr << "Unknown element type " << InputGmshEntity << endl;
                abort();
            }

            // For specific elements, call the appropriate function to perform
            // the renumbering.
            switch(it->second.m_e)
            {
                case LibUtilities::eTriangle:
                    return TriReordering(it->second);
                    break;
                case LibUtilities::eQuadrilateral:
                    return QuadReordering(it->second);
                    break;
                case LibUtilities::eTetrahedron:
                    return TetReordering(it->second);
                    break;
                case LibUtilities::ePrism:
                    return PrismReordering(it->second);
                    break;
                case LibUtilities::eHexahedron:
                    return HexReordering(it->second);
                    break;
                default:
                    break;
            }

            // Default: no reordering.
            vector<int> returnVal;
            return returnVal;
        }

        /**
         * @brief Create a reordering for triangles.
         */
        vector<int> InputGmsh::TriReordering(ElmtConfig conf)
        {
            const int order = conf.m_order;
            const int n     = order-1;

            // Copy vertices.
            vector<int> mapping(3);
            for (int i = 0; i < 3; ++i)
            {
                mapping[i] = i;
            }

            if (order == 1)
            {
                return mapping;
            }

            // Curvilinear edges.
            mapping.resize(3 + 3*n);

            for (int i = 3; i < 3+3*n; ++i)
            {
                mapping[i] = i;
            }

            if (!conf.m_faceNodes)
            {
                return mapping;
            }

            // Interior nodes.
            vector<int> interior(n*(n-1)/2);
            for (int i = 0; i < interior.size(); ++i)
            {
                interior[i] = i + 3+3*n;
            }

            if (interior.size() > 0)
            {
                interior = triTensorNodeOrdering(interior, n-1);
            }

            mapping.insert(mapping.end(), interior.begin(), interior.end());
            return mapping;
        }

        /**
         * @brief Create a reordering for quadrilaterals.
         */
        vector<int> InputGmsh::QuadReordering(ElmtConfig conf)
        {
            const int order = conf.m_order;
            const int n     = order-1;

            // Copy vertices.
            vector<int> mapping(4);
            for (int i = 0; i < 4; ++i)
            {
                mapping[i] = i;
            }

            if (order == 1)
            {
                return mapping;
            }

            // Curvilinear edges.
            mapping.resize(4 + 4*n);

            for (int i = 4; i < 4+4*n; ++i)
            {
                mapping[i] = i;
            }

            if (!conf.m_faceNodes)
            {
                return mapping;
            }

            // Interior nodes.
            vector<int> interior(n*n);
            for (int i = 0; i < interior.size(); ++i)
            {
                interior[i] = i + 4+4*n;
            }

            if (interior.size() > 0)
            {
                interior = quadTensorNodeOrdering(interior, n);
            }
            mapping.insert(mapping.end(), interior.begin(), interior.end());
            return mapping;
        }

        /**
         * @brief Create a reordering for tetrahedra.
         */
        vector<int> InputGmsh::TetReordering(ElmtConfig conf)
        {
            const int order = conf.m_order;
            const int n     = order-1;
            const int n2    = n*(n-1)/2;

            int i, j;
            vector<int> mapping(4);

            // Copy vertices.
            for (i = 0; i < 4; ++i)
            {
                mapping[i] = i;
            }

            if (order == 1)
            {
                return mapping;
            }

            // Curvilinear edges.
            mapping.resize(4 + 6*n);

            // Curvilinear edges.
            static int gmshToNekEdge[6] = {0,1,2,3,5,4};
            static int gmshToNekRev [6] = {0,0,1,1,1,1};

            // Reorder edges.
            int offset, cnt = 4;
            for (i = 0; i < 6; ++i)
            {
                offset = 4 + n * gmshToNekEdge[i];

                if (gmshToNekRev[i])
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+n-j-1] = cnt++;
                    }
                }
                else
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+j] = cnt++;
                    }
                }
            }

            if (conf.m_faceNodes == false || n2 == 0)
            {
                return mapping;
            }

            // Curvilinear faces.
            mapping.resize(4 + 6*n + 4*n2);

            static int gmshToNekFace[4] = {0,1,3,2};

            vector<int> triVertId(3);
            triVertId[0] = 0; triVertId[1] = 1; triVertId[2] = 2;

            // Loop over Gmsh faces
            for (i = 0; i < 4; ++i)
            {
                int face    = gmshToNekFace[i];
                int offset2 = 4 + 6*n + i   *n2;
                offset      = 4 + 6*n + face*n2;

                // Create a list of interior face nodes for this face only.
                vector<int> faceNodes(n2);
                vector<int> toAlign  (3);
                for (j = 0; j < n2; ++j)
                {
                    faceNodes[j] = offset2 + j;
                }

                // Now get the reordering of this face, which puts Gmsh
                // recursive ordering into Nektar++ row-by-row order.
                vector<int> tmp = triTensorNodeOrdering(faceNodes, n-1);
                HOTriangle<int> hoTri(triVertId, tmp);

                // Apply reorientation
                if (i == 0 || i == 2)
                {
                    // Triangle verts {0,2,1} --> {0,1,2}
                    toAlign[0] = 0; toAlign[1] = 2; toAlign[2] = 1;
                    hoTri.Align(toAlign);
                }
                else if (i == 3)
                {
                    // Triangle verts {1,2,0} --> {0,1,2}
                    toAlign[0] = 1; toAlign[1] = 2; toAlign[2] = 0;
                    hoTri.Align(toAlign);
                }

                // Fill in mapping.
                for (j = 0; j < n2; ++j)
                {
                    mapping[offset+j] = hoTri.surfVerts[j];
                }
            }

            return mapping;
        }

        /**
         * @brief Create a reordering for prisms.
         *
         * Note that whilst Gmsh MSH files have the capability to support
         * high-order prisms, presently Gmsh does not seem to be capable of
         * generating higher than second-order prismatic meshes, so most of the
         * following is untested.
         */
        vector<int> InputGmsh::PrismReordering(ElmtConfig conf)
        {
            const int order = conf.m_order;
            const int n     = order-1;

            int i;
            vector<int> mapping(6);

            // To get from Gmsh -> Nektar++ prism, coordinates axes are
            // different; need to mirror in the triangular faces, and then
            // reorder vertices to make ordering anticlockwise on base quad.
            static int gmshToNekVerts[6] = {3,4,1,0,5,2};

            for (i = 0; i < 6; ++i)
            {
                mapping[i] = gmshToNekVerts[i];
            }

            if (order == 1)
            {
                return mapping;
            }

            // Curvilinear edges.
            mapping.resize(6 + 9*n);

            static int gmshToNekEdge[9] = {2,7,3,6,1,8,0,4,5};
            static int gmshToNekRev [9] = {1,0,1,0,1,1,0,0,0};

            // Reorder edges.
            int offset, cnt = 6;
            for (i = 0; i < 9; ++i)
            {
                offset = 6 + n * gmshToNekEdge[i];

                if (gmshToNekRev[i])
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+n-j-1] = cnt++;
                    }
                }
                else
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+j] = cnt++;
                    }
                }
            }

            if (conf.m_faceNodes == false)
            {
                return mapping;
            }

            if (order > 2)
            {
                cerr << "Gmsh prisms of order > 2 with face curvature "
                     << "not supported in MeshConvert (or indeed Gmsh at"
                     << "time of writing)." << endl;
                abort();
            }

            mapping.resize(18);
            mapping[15] = 15;
            mapping[16] = 17;
            mapping[17] = 16;

            return mapping;
        }

        /**
         * @brief Create a reordering for hexahedra.
         */
        vector<int> InputGmsh::HexReordering(ElmtConfig conf)
        {
            const int order = conf.m_order;
            const int n     = order-1;
            const int n2    = n * n;
            int i, j, k;

            vector<int> mapping;

            // Map taking Gmsh edges to Nektar++ edges.
            static int gmshToNekEdge[12] = {
                0, -3, 4, 1, 5, 2, 6, 7, 8, -11, 9, 10
            };

            // Push back vertices.
            mapping.resize(8);
            for (i = 0; i < 8; ++i)
            {
                mapping[i] = i;
            }

            if (order == 1)
            {
                return mapping;
            }

            // Curvilinear edges
            mapping.resize(8 + 12*n);

            // Reorder edges.
            int cnt = 8, offset;
            for (i = 0; i < 12; ++i)
            {
                int edge = gmshToNekEdge[i];
                offset = 8 + n * abs(edge);

                if (edge < 0)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+n-j-1] = cnt++;
                    }
                }
                else
                {
                    for (int j = 0; j < n; ++j)
                    {
                        mapping[offset+j] = cnt++;
                    }
                }
            }

            if (conf.m_faceNodes == false || n2 == 0)
            {
                return mapping;
            }

            // Curvilinear face nodes.
            mapping.resize(8 + 12*n + 6*n2);

            // Map which takes Gmsh -> Nektar++ faces in the local element.
            static int gmsh2NekFace[6] = {0,1,4,2,3,5};

            // Map which defines orientation between Gmsh and Nektar++ faces.
            StdRegions::Orientation faceOrient[6] = {
                StdRegions::eDir1FwdDir2_Dir2FwdDir1,
                StdRegions::eDir1FwdDir1_Dir2FwdDir2,
                StdRegions::eDir1FwdDir2_Dir2FwdDir1,
                StdRegions::eDir1FwdDir1_Dir2FwdDir2,
                StdRegions::eDir1BwdDir1_Dir2FwdDir2,
                StdRegions::eDir1FwdDir1_Dir2FwdDir2
            };

            for (i = 0; i < 6; ++i)
            {
                int face = gmsh2NekFace[i];
                int offset2 = 8 + 12 * n + i * n2;
                offset = 8 + 12 * n + face * n2;

                // Create a list of interior face nodes for this face only.
                vector<int> faceNodes(n2);
                for (j = 0; j < n2; ++j)
                {
                    faceNodes[j] = offset2 + j;
                }

                // Now get the reordering of this face, which puts Gmsh
                // recursive ordering into Nektar++ row-by-row order.
                vector<int> tmp = quadTensorNodeOrdering(faceNodes, n);

                // Finally reorient the face according to the geometry
                // differences.
                if (faceOrient[i] == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
                {
                    // Orientation is the same, just copy.
                    for (j = 0; j < n2; ++j)
                    {
                        mapping[offset+j] = tmp[j];
                    }
                }
                else if (faceOrient[i] == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
                {
                    // Tranposed faces
                    for (j = 0; j < n; ++j)
                    {
                        for (k = 0; k < n; ++k)
                        {
                            mapping[offset + j*n + k] = tmp[k*n + j];
                        }
                    }
                }
                else if (faceOrient[i] == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
                {
                    for (j = 0; j < n; ++j)
                    {
                        for (k = 0; k < n; ++k)
                        {
                            mapping[offset + j*n + k] = tmp[j*n + (n-k-1)];
                        }
                    }
                }
            }

            if (conf.m_volumeNodes == false)
            {
                return mapping;
            }

            const int totPoints = (order+1) * (order+1) * (order+1);
            mapping.resize(totPoints);

            // TODO: Fix ordering of volume nodes.
            for (i = 8 + 12*n + 6*n2; i < totPoints; ++i)
            {
                mapping[i] = i;
            }

            return mapping;
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
            using namespace LibUtilities;
            std::map<unsigned int, ElmtConfig> tmp;

            //                    Elmt type,   order,  face, volume
            tmp[  1] = ElmtConfig(eSegment,        1,  true,  true);
            tmp[  2] = ElmtConfig(eTriangle,       1,  true,  true);
            tmp[  3] = ElmtConfig(eQuadrilateral,  1,  true,  true);
            tmp[  4] = ElmtConfig(eTetrahedron,    1,  true,  true);
            tmp[  5] = ElmtConfig(eHexahedron,     1,  true,  true);
            tmp[  6] = ElmtConfig(ePrism,          1,  true,  true);
            tmp[  7] = ElmtConfig(ePyramid,        1,  true,  true);
            tmp[  8] = ElmtConfig(eSegment,        2,  true,  true);
            tmp[  9] = ElmtConfig(eTriangle,       2,  true,  true);
            tmp[ 10] = ElmtConfig(eQuadrilateral,  2,  true,  true);
            tmp[ 11] = ElmtConfig(eTetrahedron,    2,  true,  true);
            tmp[ 12] = ElmtConfig(eHexahedron,     2,  true,  true);
            tmp[ 13] = ElmtConfig(ePrism,          2,  true,  true);
            tmp[ 14] = ElmtConfig(ePyramid,        2,  true,  true);
            tmp[ 15] = ElmtConfig(ePoint,          1,  true, false);
            tmp[ 16] = ElmtConfig(eQuadrilateral,  2, false, false);
            tmp[ 17] = ElmtConfig(eHexahedron,     2, false, false);
            tmp[ 18] = ElmtConfig(ePrism,          2, false, false);
            tmp[ 19] = ElmtConfig(ePyramid,        2, false, false);
            tmp[ 20] = ElmtConfig(eTriangle,       3, false, false);
            tmp[ 21] = ElmtConfig(eTriangle,       3,  true, false);
            tmp[ 22] = ElmtConfig(eTriangle,       4, false, false);
            tmp[ 23] = ElmtConfig(eTriangle,       4,  true, false);
            tmp[ 24] = ElmtConfig(eTriangle,       5, false, false);
            tmp[ 25] = ElmtConfig(eTriangle,       5,  true, false);
            tmp[ 26] = ElmtConfig(eSegment,        3,  true, false);
            tmp[ 27] = ElmtConfig(eSegment,        4,  true, false);
            tmp[ 28] = ElmtConfig(eSegment,        5,  true, false);
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
            tmp[ 62] = ElmtConfig(eSegment,        6,  true, false);
            tmp[ 63] = ElmtConfig(eSegment,        7,  true, false);
            tmp[ 64] = ElmtConfig(eSegment,        8,  true, false);
            tmp[ 65] = ElmtConfig(eSegment,        9,  true, false);
            tmp[ 66] = ElmtConfig(eSegment,       10,  true, false);
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
            tmp[118] = ElmtConfig(ePyramid,        3,  true,  true);
            tmp[119] = ElmtConfig(ePyramid,        4,  true,  true);
            tmp[120] = ElmtConfig(ePyramid,        5,  true,  true);
            tmp[121] = ElmtConfig(ePyramid,        6,  true,  true);
            tmp[122] = ElmtConfig(ePyramid,        7,  true,  true);
            tmp[123] = ElmtConfig(ePyramid,        8,  true,  true);
            tmp[124] = ElmtConfig(ePyramid,        9,  true,  true);
            tmp[125] = ElmtConfig(ePyramid,        3,  true, false);
            tmp[126] = ElmtConfig(ePyramid,        4,  true, false);
            tmp[127] = ElmtConfig(ePyramid,        5,  true, false);
            tmp[128] = ElmtConfig(ePyramid,        6,  true, false);
            tmp[129] = ElmtConfig(ePyramid,        7,  true, false);
            tmp[130] = ElmtConfig(ePyramid,        7,  true, false);
            tmp[131] = ElmtConfig(ePyramid,        8,  true, false);

            return tmp;
        }
    }
}
