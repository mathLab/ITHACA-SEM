////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/VolumeMeshing/BLMeshing/BLMesh.h>
#include <NekMeshUtils/VolumeMeshing/TetMeshing/TetMesh.h>
#include <NekMeshUtils/VolumeMeshing/VolumeMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey VolumeMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "volumemesh"),
    VolumeMesh::create,
    "Generates a volume mesh");

VolumeMesh::VolumeMesh(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["blsurfs"] =
        ConfigOption(true, "0", "Generate prisms on these surfs");
}

VolumeMesh::~VolumeMesh()
{
}

void VolumeMesh::Process()
{
    if (m_mesh->m_verbose)
        cout << endl << "Volume meshing" << endl;

    /*map<int, FaceSharedPtr> surftopriface;
    // map of surface element id to opposite prism
    // face for psudo surface in tetmesh
*/
    TetMeshSharedPtr tet;/*
    if (m_makeBL)
    {
        BLMeshSharedPtr blmesh = MemoryManager<BLMesh>::AllocateSharedPtr(
                                        m_mesh->m_cad, m_mesh, blsurfs, m_blthick);

        blmesh->Mesh();*/

        //m_mesh->m_element[2].clear();
        //m_mesh->m_expDim = 3;
        /*ClearElementLinks();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        return;*/

        //m_surfacemesh->Remesh(m_blmesh);

        /*m_mesh->m_nummode = 2;
        m_mesh->m_expDim = 3;
        m_mesh->m_element[2].clear();
        ClearElementLinks();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        return;*/

        //create tet mesh
        /*m_mesh->m_expDim = 3;

        m_tet = MemoryManager<TetMesh>::AllocateSharedPtr(
            m_mesh, blmesh);
    }
    else
    {*/
    m_mesh->m_expDim = 3;
    tet = MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh);

    //}

    tet->Mesh();

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    if (m_mesh->m_verbose)
        cout << endl;
}
}
}
