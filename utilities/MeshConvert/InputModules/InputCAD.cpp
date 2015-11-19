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
    m_orelax = pSession->GetSolverInfo("OCtreeRelax");
    if(boost::iequals(m_orelax,"TRUE"))
        m_octreeRelax = true;
    else
        m_octreeRelax = false;

    CADSystemSharedPtr m_cad = MemoryManager<CADSystem>::
                                            AllocateSharedPtr(m_CADName);

    if(m_mesh->m_verbose)
    {
        cout << "Building mesh for: " << m_CADName << endl;
        if(m_octreeRelax)
            cout << "With a relaxed octree" << endl;
    }

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

    m_cad->SmallFeatureAnalysis(m_minDelta);

    //create octree
    OctreeSharedPtr m_octree = MemoryManager<Octree>::AllocateSharedPtr(m_cad,
                                    m_mesh->m_verbose, m_minDelta,
                                    m_maxDelta, m_eps, m_octreeRelax);

    m_octree->Build();

    m_mesh->m_expDim = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode = m_order+1;
    //m_mesh->m_nummode = 2;

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

    int nNeg = 0;
    for(int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        bool boundary = false;
        vector<FaceSharedPtr> fs = m_mesh->m_element[3][i]->GetFaceList();
        for(int j = 0; j < fs.size(); j++)
        {
            if(fs[j]->m_faceNodes.size() > 0)
            {
                boundary = true;
                break;
            }
        }
        if(boundary)
        {
            SpatialDomains::GeometrySharedPtr geom = m_mesh->m_element[3][i]->GetGeom(m_mesh->m_spaceDim);
            SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

            if(!gfac->IsValid())
            {
                nNeg++;
                cout << "  - " << m_mesh->m_element[3][i]->GetId() << " ("
                     << LibUtilities::ShapeTypeMap[m_mesh->m_element[3][i]->GetConf().m_e]
                     << ")" << endl;
            }
        }
    }

    if(m_mesh->m_verbose)
    {
        if(nNeg > 0)
            cout << "Warning: " << nNeg << " invalid elements" << endl;
        else
            cout << "0 invalid elements" << endl;
    }

    if(m_mesh->m_verbose)
        cout << endl;
}

}
}
