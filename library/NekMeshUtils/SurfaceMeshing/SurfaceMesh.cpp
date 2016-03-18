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
#include <algorithm>

#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>

#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void SurfaceMesh::Mesh()
{
    if (m_mesh->m_verbose)
        cout << endl << "Surface meshing" << endl;

    if (m_mesh->m_verbose)
        cout << endl << "\tCurve meshing:" << endl << endl;

    m_mesh->m_numNodes = m_cad->GetNumVerts();

    // linear mesh all curves
    for (int i = 1; i <= m_cad->GetNumCurve(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_cad->GetNumCurve(), "Curve progress");
        }

        m_curvemeshes[i] = MemoryManager<CurveMesh>::AllocateSharedPtr(
            i, m_mesh, m_cad->GetCurve(i), m_octree);

        m_curvemeshes[i]->Mesh();
    }

    if (m_mesh->m_verbose)
        cout << endl << "\tFace meshing:" << endl << endl;

    // linear mesh all surfaces
    for (int i = 1; i <= m_cad->GetNumSurf(); i++)
    {
        if (m_mesh->m_verbose)
        {
            LibUtilities::PrintProgressbar(
                i, m_cad->GetNumSurf(), "Face progress");
        }

        vector<unsigned int>::iterator it;

        it = find(m_symsurfs.begin(), m_symsurfs.end(), i);
        if (it == m_symsurfs.end())
        {
            m_facemeshes[i] = MemoryManager<FaceMesh>::AllocateSharedPtr(
                i, m_mesh, m_cad->GetSurf(i), m_octree, m_curvemeshes);
        }
        else
        {
            m_facemeshes[i] = MemoryManager<FaceMesh>::AllocateSharedPtr(
                i, m_mesh, m_cad->GetSurf(i), m_octree, m_curvemeshes, m_bl);
        }

        m_facemeshes[i]->Mesh();
    }
}

// this mesh is valided that each egde is listed twice in the triangles
void SurfaceMesh::Validate()
{
    if (m_mesh->m_verbose)
        cout << endl << "\tVerifying surface mesh" << endl;

    EdgeSet::iterator it;

    for (it = m_mesh->m_edgeSet.begin(); it != m_mesh->m_edgeSet.end(); it++)
    {
        if ((*it)->m_elLink.size() != 2)
        {
            if (m_mesh->m_verbose)
                cout << "\t\tFailed" << endl;
            ASSERTL0(false, "edge not listed twice");
        }
    }

    if (m_mesh->m_verbose)
        cout << "\t\tPassed" << endl;
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
