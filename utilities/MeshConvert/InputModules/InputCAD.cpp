////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.cpp
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
//  Description: create mesh from cad using mesh utils
//
////////////////////////////////////////////////////////////////////////////////

#include <string>



#include <LibUtilities/CADSystem/CADSystem.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/SurfaceMeshing.h>
#include <MeshUtils/TetMeshing/TetMesh.h>
#include <MeshUtils/MeshElem.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "../MeshElements.h"
#include "InputCAD.h"

using namespace std;
namespace Nektar{
namespace Utilities{

ModuleKey InputCAD::className =
GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "CAD"), InputCAD::create,
    "Reads CAD geometry and will generate the mesh file.");

/**
 * @brief Set up InputCAD object.
 */
InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
{
    m_config["min"] = ConfigOption(false,"-1",
            "minimum delta to be in mesh");
    m_config["max"] = ConfigOption(false,"-1",
            "maximum delta to be in mesh");
    m_config["eps"] = ConfigOption(false,"-1",
            "sensitivity to curvature");
    m_config["order"] = ConfigOption(false,"-1",
            "order of the mesh to be produced");
}

InputCAD::~InputCAD()
{

}


void InputCAD::Process()
{
    string CADName = m_config["infile"].as<string>();
    NekDouble m_minDelta = m_config["min"].as<NekDouble>();
    NekDouble m_maxDelta = m_config["max"].as<NekDouble>();
    NekDouble m_eps = m_config["eps"].as<NekDouble>();
    int m_order = m_config["order"].as<int>();

    ASSERTL0(!(m_minDelta == -1 || m_maxDelta == -1 ||
               m_eps == -1 || m_order == -1),
             "User parameters required");


    LibUtilities::CADSystemSharedPtr m_cad =
        MemoryManager<LibUtilities::CADSystem>::AllocateSharedPtr(CADName);

    if(m_mesh->m_verbose)
        cout << "Building mesh for: " << m_cad->GetName() << endl;

    ASSERTL0(m_cad->LoadCAD(),
             "Failed to load CAD");

    if(m_mesh->m_verbose)
    {
        cout << "With parameters:" << endl;
        cout << "\tmin delta: " << m_minDelta << endl
             << "\tmax delta: " << m_maxDelta << endl
             << "\tesp: " << m_eps << endl
             << "\torder: " << m_order << endl;
        m_cad->Report();
    }

    //create octree
    MeshUtils::OctreeSharedPtr m_octree =
        MemoryManager<MeshUtils::Octree>::AllocateSharedPtr(m_cad,
                                m_mesh->m_verbose);

    m_octree->Build(m_minDelta, m_maxDelta, m_eps);

    //create surface mesh
    MeshUtils::SurfaceMeshingSharedPtr m_surfacemeshing =
        MemoryManager<MeshUtils::SurfaceMeshing>::
            AllocateSharedPtr(m_mesh->m_verbose,m_cad,m_octree,m_order);

    m_surfacemeshing->Mesh();

    m_surfacemeshing->HOSurf();

    //create tet mesh
    MeshUtils::TetMeshSharedPtr m_tet =
        MemoryManager<MeshUtils::TetMesh>::AllocateSharedPtr(
            m_mesh->m_verbose, m_cad, m_octree, m_surfacemeshing);

    m_tet->Mesh();

    m_mesh->m_expDim = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_order = m_order+1;

    m_mesh->m_fields.push_back("u");
    m_mesh->m_fields.push_back("v");
    m_mesh->m_fields.push_back("w");
    m_mesh->m_fields.push_back("p");

    map<int, MeshUtils::MeshNodeSharedPtr> Nodes;
    map<int, MeshUtils::MeshEdgeSharedPtr> Edges;
    map<int, MeshUtils::MeshTriSharedPtr> Tris;
    map<int, MeshUtils::MeshTetSharedPtr> Tets;

    m_tet->Get(Nodes,Edges,Tris,Tets);

    map<int, MeshUtils::MeshNodeSharedPtr>::iterator nit;
    map<int, NodeSharedPtr> allnodes;

    //extract all nodes and make mesh convert nodes
    for(nit = Nodes.begin(); nit != Nodes.end(); nit++)
    {
        Array<OneD, NekDouble> loc = nit->second->GetLoc();
        NodeSharedPtr nn =
                boost::shared_ptr<Node>(
                            new Node(nit->second->GetId(),loc[0],
                                     loc[1],loc[2]));
        nn->m_mid = nit->second->GetId();
        allnodes[nit->second->GetId()] = nn;
    }

    map<int, MeshUtils::MeshTetSharedPtr>::iterator tetit;
    //extract all tets and make mesh convert tets
    for(tetit = Tets.begin(); tetit != Tets.end(); tetit++)
    {
        Array<OneD, int> n = tetit->second->GetN();

        vector<NodeSharedPtr> localnode;

        for(int j = 0; j < 4; j++)
        {
            localnode.push_back(allnodes[n[j]]);
        }


        ElmtConfig conf(LibUtilities::eTetrahedron,1,false,false);
        vector<int> tags;
        tags.push_back(0);
        ElementSharedPtr E = GetElementFactory().
                    CreateInstance(LibUtilities::eTetrahedron,
                                   conf,localnode,tags);
        m_mesh->m_element[m_mesh->m_expDim].push_back(E);
    }

    map<int, MeshUtils::MeshTriSharedPtr>::iterator trit;
    //extract all triagnles and make mesh covnert surfaces
    for(trit = Tris.begin(); trit != Tris.end(); trit++)
    {
        Array<OneD, int> n = trit->second->GetN();

        vector<NodeSharedPtr> localnode;
        for(int j = 0; j < 3; j++)
        {
            localnode.push_back(allnodes[n[j]]);
        }

        Array<OneD, int> eg = trit->second->GetE();
        for(int j = 0; j < 3; j++)
        {
            MeshUtils::MeshEdgeSharedPtr e;
            e = Edges[eg[j]];
            vector<int> hon = e->GetHONodes(n[j]);
            for(int k = 0; k < hon.size(); k++)
            {
                localnode.push_back(allnodes[hon[k]]);
            }
        }

        ElmtConfig conf(LibUtilities::eTriangle,m_order,false,false,false);

        vector<int> tags;
        tags.push_back(trit->second->Getcid());
        ElementSharedPtr E = GetElementFactory().
                    CreateInstance(LibUtilities::eTriangle,
                                   conf,localnode,tags);
        E->SetTriID(trit->first);
        m_mesh->m_element[m_mesh->m_expDim-1].push_back(E); //needs to be set to -1
    }

    ProcessVertices  ();
    ProcessEdges     ();
    ProcessFaces     ();
    ProcessElements  ();
    ProcessComposites();

    //look over all surface elements
    //get the face between it and the tet
    //for all the edges in the face insert curved information
    //in the face first orientate and then add the face interior nodes
    for(int i = 0; i < m_mesh->m_element[m_mesh->m_expDim-1].size(); i++)
    {
        FaceSharedPtr f =
                m_mesh->m_element[m_mesh->m_expDim-1][i]->GetFaceLink();
        vector<EdgeSharedPtr> egs = f->m_edgeList;
        for(int j = 0; j < egs.size(); j++)
        {
            if(egs[j]->m_edgeNodes.size()>0)
                continue;

            NodeSharedPtr n1 = egs[j]->m_n1;
            NodeSharedPtr n2 = egs[j]->m_n2;
            int n1ID = n1->m_mid;
            int n2ID = n2->m_mid;

            int edgekey = Nodes[n1ID]->EdgeInCommon(Nodes[n2ID]);

            ASSERTL0(edgekey != -1, "no edge found");

            vector<int> honode = Edges[edgekey]->GetHONodes(n1ID);

            vector<NodeSharedPtr> localhonode;
            for(int k = 0; k < honode.size(); k++)
            {
                localhonode.push_back(allnodes[honode[k]]);
            }

            egs[j]->m_edgeNodes = localhonode;
            egs[j]->m_curveType = LibUtilities::ePolyEvenlySpaced;
        }

        int t = m_mesh->m_element[m_mesh->m_expDim-1][i]->GetTriID();

        vector<int> honode = Tris[t]->GetHONodes();
        vector<NodeSharedPtr> localhonode;

        for(int j = 0; j < honode.size(); j++)
        {
            localhonode.push_back(allnodes[honode[j]]);
        }

        f->m_faceNodes = localhonode;

        Array<OneD, int> n = Tris[t]->GetN();
        vector<int> trivert(3);
        vector<int> facevert(3);
        int aligned = 0;
        for(int j = 0; j < 3; j++)
        {
            trivert[j] = n[j];
            facevert[j] = f->m_vertexList[j]->m_mid;
            if(trivert[j] == facevert[j]) aligned++;
        }

        if(aligned != 3)
        {
            HOTriangle<NodeSharedPtr> hoTri(trivert,localhonode);
            hoTri.Align(facevert);

            f->m_faceNodes = hoTri.surfVerts;
        }

        f->m_curveType = LibUtilities::eNodalTriEvenlySpaced;
    }

    if(m_mesh->m_verbose)
        cout << endl;
}

}
}
