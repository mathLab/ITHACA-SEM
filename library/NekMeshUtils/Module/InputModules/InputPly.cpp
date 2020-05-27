////////////////////////////////////////////////////////////////////////////////
//
//  File: InputPly.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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

#include <NekMeshUtils/MeshElements/Element.h>

#include "InputPly.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputPly::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "ply"),
    InputPly::create,
    "Reads ply triangulation format.");

InputPly::InputPly(MeshSharedPtr m) : InputModule(m)
{
}

InputPly::~InputPly()
{
}

/**
 *
 * @param   pFilename           Filename of Gmsh file to read.
 */
void InputPly::Process()
{

    // Open the file stream.
    OpenStream();

    ReadPly(m_mshFile);

    m_mshFile.reset();

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

void InputPly::ReadPly(io::filtering_istream &mshFile, NekDouble scale)
{
    m_mesh->m_expDim = 0;
    string line;
    int nVertices                  = 0;
    int nEntities                  = 0;
    int nProperties                = 0;
    LibUtilities::ShapeType elType = LibUtilities::eTriangle;
    map<string, int> propMap;

    if (m_mesh->m_verbose)
    {
        cout << "InputPly: Start reading file..." << endl;
    }

    while (!mshFile.eof())
    {
        getline(mshFile, line);
        stringstream s(line);
        string word;
        s >> word;
        if (word == "format")
        {
            s >> word;
            if (word != "ascii")
            {
                ASSERTL0(false, "InputPly file currently only set up to read "
                                "ascii formatted ply files");
            }
        }
        else if (word == "element")
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
        else if (word == "property")
        {
            s >> word >> word;
            propMap[word] = nProperties++;
        }
        else if (word == "end_header")
        {
            // Read nodes
            vector<double> data(nProperties);
            for (int i = 0; i < nVertices; ++i)
            {
                getline(mshFile, line);
                stringstream st(line);

                for (int j = 0; j < nProperties; ++j)
                {
                    st >> data[j];
                }

                double x = data[propMap["x"]];
                double y = data[propMap["y"]];
                double z = data[propMap["z"]];

                if ((y * y) > 0.000001 && m_mesh->m_spaceDim != 3)
                {
                    m_mesh->m_spaceDim = 2;
                }
                if ((z * z) > 0.000001)
                {
                    m_mesh->m_spaceDim = 3;
                }

                x *= scale;
                y *= scale;
                z *= scale;

                m_mesh->m_node.push_back(
                    std::shared_ptr<Node>(new Node(i, x, y, z)));

                // Read vertex normals.
                if (propMap.count("nx") > 0)
                {
                    double nx                  = data[propMap["nx"]];
                    double ny                  = data[propMap["ny"]];
                    double nz                  = data[propMap["nz"]];
                    m_mesh->m_vertexNormals[i] = Node(0, nx, ny, nz);
                }
            }

            // Read elements
            for (int i = 0; i < nEntities; ++i)
            {
                getline(mshFile, line);
                stringstream st(line);
                int id = 0;

                // Create element tags
                vector<int> tags;
                tags.push_back(0); // composite

                // Read element node list
                st >> id;
                vector<NodeSharedPtr> nodeList;
                for (int k = 0; k < 3; ++k)
                {
                    int node = 0;
                    st >> node;
                    nodeList.push_back(m_mesh->m_node[node]);
                }

                // Create element
                ElmtConfig conf(elType, 1, false, false);
                ElementSharedPtr E = GetElementFactory().CreateInstance(
                    elType, conf, nodeList, tags);

                // Determine mesh expansion dimension
                if (E->GetDim() > m_mesh->m_expDim)
                {
                    m_mesh->m_expDim = E->GetDim();
                }
                m_mesh->m_element[E->GetDim()].push_back(E);
            }
        }
    }
}
}
}
