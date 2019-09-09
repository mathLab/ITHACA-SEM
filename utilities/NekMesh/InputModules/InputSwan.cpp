////////////////////////////////////////////////////////////////////////////////
//
//  File: InputSwan.cpp
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
//  Description: Swansea session converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "InputSwan.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputSwan::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "plt"),
    InputSwan::create,
    "Reads Swansea plt format for third-order tetrahedra.");

InputSwan::InputSwan(MeshSharedPtr m) : InputModule(m)
{
}

InputSwan::~InputSwan()
{
}

void InputSwan::Process()
{
    // Open the file stream.
    OpenStream();

    vector<vector<NodeSharedPtr> > elementList;
    vector<int> tmp, tets;
    vector<double> pts;

    if (m_mesh->m_verbose)
    {
        cout << "InputSwan: Start reading file..." << endl;
    }

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;

    // First read in header; 4 integers containing number of tets,
    // number of points, nas2 (unknown) and order of the grid
    // (i.e. GridOrder = 3 => cubic mesh).
    tmp.resize(6);
    m_mshFile.read(reinterpret_cast<char *>(&tmp[0]),
                   static_cast<int>(6 * sizeof(int)));

    if (tmp[0] != tmp[5] || tmp[0] != 4 * sizeof(int))
    {
        cout << "Header data broken" << endl;
        m_mshFile.reset();
        return;
    }

    int NB_Tet    = tmp[1];
    int NB_Points = tmp[2];
    int GridOrder = tmp[4];
    int ND        = (GridOrder + 1) * (GridOrder + 2) * (GridOrder + 3) / 6;

    cout << "NB_Tet    = " << NB_Tet << endl;
    cout << "NB_Points = " << NB_Points << endl;
    cout << "GridOrder = " << GridOrder << endl;
    cout << "ND        = " << ND << endl;

    // Now read in list of tetrahedra. Each tet has ND points, and
    // data are ordered with memory traversed fastest with point
    // number.
    tets.resize(ND * NB_Tet + 2);
    m_mshFile.read(reinterpret_cast<char *>(&tets[0]),
                   static_cast<int>((ND * NB_Tet + 2) * sizeof(int)));

    if (tets[0] != tets[ND * NB_Tet + 1] ||
        tets[0] != ND * NB_Tet * sizeof(int))
    {
        cout << "ERROR [InputSwan]: Tetrahedron data broken." << endl;
        m_mshFile.reset();
        return;
    }

    // Finally, read point data: NB_Points tuples (x,y,z).
    tmp.resize(2);
    pts.resize(3 * NB_Points);
    m_mshFile.read(reinterpret_cast<char *>(&tmp[0]),
                   static_cast<int>(sizeof(int)));
    m_mshFile.read(reinterpret_cast<char *>(&pts[0]),
                   static_cast<int>(3 * NB_Points * sizeof(double)));
    m_mshFile.read(reinterpret_cast<char *>(&tmp[1]),
                   static_cast<int>(sizeof(int)));

    if (tmp[0] != tmp[1] || tmp[0] != 3 * NB_Points * sizeof(double))
    {
        cout << "ERROR [InputSwan]: Point data broken." << endl;
        m_mshFile.reset();
        return;
    }

    int vid                        = 0, i, j;
    LibUtilities::ShapeType elType = LibUtilities::eTetrahedron;

    // Read in list of vertices.
    for (i = 0; i < NB_Points; ++i)
    {
        double x = pts[i];
        double y = pts[1 * NB_Points + i];
        double z = pts[2 * NB_Points + i];
        m_mesh->m_node.push_back(
            std::shared_ptr<Node>(new Node(vid, x, y, z)));
        vid++;
    }

    // Iterate over list of tetrahedra: for each, create nodes. At the
    // moment discard high order data and create linear mesh.
    for (i = 0; i < NB_Tet; ++i)
    {
        vector<NodeSharedPtr> nodeList;
        for (j = 0; j < 20; ++j)
        {
            int vid = tets[j * NB_Tet + i + 1];
            nodeList.push_back(m_mesh->m_node[vid - 1]);
        }

        vector<int> tags;
        tags.push_back(0);      // composite
        tags.push_back(elType); // element type

        ElmtConfig conf(elType, 3, true, true);
        ElementSharedPtr E =
            GetElementFactory().CreateInstance(elType, conf, nodeList, tags);
        m_mesh->m_element[3].push_back(E);
    }

    // Attempt to read in composites. Need to determine number of
    // triangles from data.
    tmp.resize(2);
    m_mshFile.read(reinterpret_cast<char *>(&tmp[0]),
                   static_cast<int>(sizeof(int)));
    int n_tri = tmp[0] / sizeof(int) / 5;
    tets.resize(n_tri * 5);
    m_mshFile.read(reinterpret_cast<char *>(&tets[0]),
                   static_cast<int>(tmp[0]));
    m_mshFile.read(reinterpret_cast<char *>(&tmp[1]),
                   static_cast<int>(sizeof(int)));

    if (tmp[0] != tmp[1])
    {
        cout << "ERROR [InputSwan]: Surface data broken." << endl;
        m_mshFile.reset();
        return;
    }

    elType = LibUtilities::eTriangle;

    // Process list of triangles forming surfaces.
    for (i = 0; i < n_tri; ++i)
    {
        vector<NodeSharedPtr> nodeList;

        for (j = 0; j < 3; ++j)
        {
            nodeList.push_back(m_mesh->m_node[tets[i + j * n_tri] - 1]);
        }

        vector<int> tags;
        tags.push_back(1);      // composite
        tags.push_back(elType); // element type

        ElmtConfig conf(elType, 1, false, false);
        ElementSharedPtr E =
            GetElementFactory().CreateInstance(elType, conf, nodeList, tags);
        m_mesh->m_element[2].push_back(E);
    }

    m_mshFile.reset();

    // Process the rest of the mesh.
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}
}
}
