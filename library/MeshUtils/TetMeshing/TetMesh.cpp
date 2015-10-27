////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.cpp
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
//  Description: tet meshing methods
//
////////////////////////////////////////////////////////////////////////////////

#include <MeshUtils/TetMeshing/TetMesh.h>
#include <MeshUtils/ExtLibInterface/TetGenInterface.h>

using namespace std;
namespace Nektar
{
namespace MeshUtils
{

void TetMesh::Mesh()
{
    if(m_mesh->m_verbose)
        cout << endl << endl << "Tetrahdral mesh generation" << endl;

    TetGenInterfaceSharedPtr tetgen =
        MemoryManager<TetGenInterface>::AllocateSharedPtr();

    //build sequentially ordered maps of nodes that exist and there delta value in the octree
    map<int, NodeSharedPtr> nodesintris;
    map<int, NekDouble> nodedelta;

    NodeSet::iterator nit;
    for(nit = m_mesh->m_vertexSet.begin(); nit != m_mesh->m_vertexSet.end(); nit++)
    {
        nodesintris[(*nit)->m_id] = *(nit);
        nodedelta[(*nit)->m_id] = m_octree->Query((*nit)->GetLoc());
    }

    cout << nodesintris.size() << endl;

    tetgen->InitialMesh(nodesintris, m_mesh->m_element[2]);

    int newpb = 20;

    vector<Array<OneD, NekDouble> > newp;

    while(newpb != newp.size())
    {
        newpb = newp.size();
        newp.clear();
        tetgen->GetNewPoints(nodesintris.size(), newp);
        for(int i = 0; i < newp.size(); i++)
        {
            NekDouble d = m_octree->Query(newp[i]);
            nodedelta[nodesintris.size() + i] = d;
        }
        tetgen->RefineMesh(nodedelta);
    }

    //make new map of all nodes to build tets.
    map<int, NodeSharedPtr> nodes = nodesintris;
    for(int i = 0; i < newp.size(); i++)
    {
        NodeSharedPtr n = boost::shared_ptr<Node>(new Node(nodesintris.size() + i,
                                                    newp[i][0],newp[i][1],newp[i][2]));
        nodes[nodesintris.size() + i] = n;
    }

    tetconnect = tetgen->Extract();

    //tetgen->freetet();

    //create tets
    for(int i = 0; i < tetconnect.size(); i++)
    {
        vector<NodeSharedPtr> n;
        n.push_back(nodes[tetconnect[i][0]]);
        n.push_back(nodes[tetconnect[i][1]]);
        n.push_back(nodes[tetconnect[i][2]]);
        n.push_back(nodes[tetconnect[i][3]]);
        ElmtConfig conf(LibUtilities::eTetrahedron,1,false,false);
        vector<int> tags;
        tags.push_back(0);
        ElementSharedPtr E = GetElementFactory().
                    CreateInstance(LibUtilities::eTetrahedron, conf, n, tags);

        m_mesh->m_element[3].push_back(E);
    }

    if(m_mesh->m_verbose)
        cout << "\tTets :" << tetconnect.size() << endl;
}

}
}
