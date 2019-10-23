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
#include <algorithm>

#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey SurfaceMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "surfacemesh"),
    SurfaceMesh::create,
    "Generates a surface mesh");

SurfaceMesh::SurfaceMesh(MeshSharedPtr m) : ProcessModule(m)
{

}

SurfaceMesh::~SurfaceMesh()
{
}

void SurfaceMesh::Process()
{
    m_mesh->m_expDim--; //just to make it easier to surface mesh for now

    if (m_mesh->m_verbose)
        cout << endl << "Surface meshing" << endl;

    if (m_mesh->m_verbose)
        cout << endl << "\tCurve meshing:" << endl << endl;

    m_mesh->m_numNodes = m_mesh->m_cad->GetNumVerts();

    // linear mesh all curves
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_mesh->m_cad->GetNumCurve(), "Curve progress");
        }

        m_curvemeshes[i] = MemoryManager<CurveMesh>::AllocateSharedPtr(i, m_mesh);

        m_curvemeshes[i]->Mesh();
    }

    if (m_mesh->m_verbose)
        cout << endl << "\tFace meshing:" << endl << endl;

    bool validError = false;
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_mesh->m_cad->GetNumSurf(), "Validating curve meshes");
        }
        FaceMeshSharedPtr face =
            MemoryManager<FaceMesh>::AllocateSharedPtr(i,m_mesh,
                m_curvemeshes, i);

        validError = validError ? true : face->ValidateCurves();

        face->ValidateLoops();
    }

    ASSERTL0(!validError,"valdity error in curve meshes");

    // linear mesh all surfaces
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_mesh->m_cad->GetNumSurf(), "Face progress");
        }
        m_facemeshes[i] =
            MemoryManager<FaceMesh>::AllocateSharedPtr(i,m_mesh,
                m_curvemeshes, i);

        m_facemeshes[i]->Mesh();
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    Report();

    EdgeSet::iterator it;
    for(it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        if((*it)->m_elLink.size() != 2)
        {
            ASSERTL0(false,"surface mesh connectivity error");
        }
    }

    m_mesh->m_expDim++; //revert dim
}
void SurfaceMesh::Report()
{
    if (m_mesh->m_verbose)
    {
        int ns = m_mesh->m_vertexSet.size();
        int es = m_mesh->m_edgeSet.size();
        int ts = m_mesh->m_element[2].size();
        int ep = ns - es + ts;
        cout << endl << "\tSurface mesh statistics" << endl;
        cout << "\t\tNodes: " << ns << endl;
        cout << "\t\tEdges: " << es << endl;
        cout << "\t\tTriangles " << ts << endl;
        cout << "\t\tEuler-Poincaré characteristic: " << ep << endl;
    }
}

}
}
