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

#include <MeshUtils/MeshElements/MeshElements.h>

#include <MeshUtils/CADSystem/CADSystem.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <MeshUtils/TetMeshing/TetMesh.h>
//#include <MeshUtils/MeshElem.hpp>

#include <LibUtilities/BasicUtils/SessionReader.h>

#include "InputCAD.h"

using namespace std;
using namespace Nektar::MeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputCAD::className =
GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "mcf"), InputCAD::create,
    "Reads CAD geometry and will generate the mesh file.");

/**
 * @brief Set up InputCAD object.
 */
InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
{

}

InputCAD::~InputCAD()
{

}

void InputCAD::Process()
{
    vector<string> filename;
    filename.push_back(m_config["infile"].as<string>());

    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);

    pSession->LoadParameter("MinDelta", m_minDelta);
    pSession->LoadParameter("MaxDelta", m_maxDelta);
    pSession->LoadParameter("EPS",      m_eps);
    pSession->LoadParameter("Order",    m_order);
    m_CADName = pSession->GetSolverInfo("CADFile");

    CADSystemSharedPtr m_cad = MemoryManager<CADSystem>::
                                            AllocateSharedPtr(m_CADName);

    if(m_mesh->m_verbose)
        cout << "Building mesh for: " << m_CADName << endl;

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
    OctreeSharedPtr m_octree = MemoryManager<Octree>::AllocateSharedPtr(m_cad,
                                    m_mesh->m_verbose, m_minDelta,
                                    m_maxDelta, m_eps);

    m_octree->Build();

    m_mesh->m_expDim = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode = m_order+1;

    //create surface mesh
    m_mesh->m_expDim--; //just to make it easier to surface mesh for now
    SurfaceMeshSharedPtr m_surfacemesh = MemoryManager<SurfaceMesh>::
                                AllocateSharedPtr(m_mesh, m_cad, m_octree);

    m_surfacemesh->Mesh();

    ProcessVertices  ();
    ProcessEdges     ();
    ProcessFaces     ();
    ProcessElements  ();
    ProcessComposites();

    m_surfacemesh->Report();

    m_surfacemesh->Validate();

    m_surfacemesh->Optimise();

    m_surfacemesh->HOAwareness();

    ClearElementLinks(); //mesh needs reprocessing to clean element and edge lists, easiest way to do it
    ProcessVertices  ();
    ProcessEdges     ();
    ProcessFaces     ();
    ProcessElements  ();
    ProcessComposites();

    m_surfacemesh->Validate();

    m_mesh->m_expDim = 3;
    m_mesh->m_fields.push_back("u");
    m_mesh->m_fields.push_back("v");
    m_mesh->m_fields.push_back("w");
    m_mesh->m_fields.push_back("p");

    //create tet mesh
    TetMeshSharedPtr m_tet =
                MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh, m_octree);

    m_tet->Mesh();

    ClearElementLinks();
    ProcessVertices  ();
    ProcessEdges     ();
    ProcessFaces     ();
    ProcessElements  ();
    ProcessComposites();

    m_surfacemesh->HOSurf();

    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        FaceSharedPtr f = m_mesh->m_element[2][i]->GetFaceLink();
        SpatialDomains::GeometrySharedPtr geom = f->GetGeom(m_mesh->m_spaceDim);
        SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

        if(!gfac->IsValid())
        {
            cout << "warning: invalid face curavture in boundary element" << endl;
        }
    }

    /*

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
            egs[j]->m_curveType = LibUtilities::eGaussLobattoLegendre;
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
        //f->m_faceNodes.clear();
        f->m_curveType = LibUtilities::eNodalTriFekete;
    }*/

    if(m_mesh->m_verbose)
        cout << endl;
}

}
}
