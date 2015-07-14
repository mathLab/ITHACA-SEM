////////////////////////////////////////////////////////////////////////////////
//
//  File: octree.h
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <sstream>

#include <MeshUtils/SurfaceMeshing.h>

using namespace std;
namespace Nektar{
namespace MeshUtils {

    void SurfaceMeshing::Mesh()
    {
        for(int i = 1; i <= m_cad->GetNumCurve(); i++)
        {
            if(m_verbose)
                cout << endl << "Meshing Curve: " << i << endl;

            m_curvemeshes[i] =
                MemoryManager<CurveMesh>::AllocateSharedPtr(
                    m_verbose, i, m_cad->GetCurve(i), m_octree);

            m_curvemeshes[i]->Mesh(Nodes, Edges);

        }

        for(int i = 1; i <= m_cad->GetNumSurf(); i++)
        {
            if(m_verbose)
                cout << endl << "Surface: " <<  i <<  endl;
            m_surfacemeshes[i] =
                MemoryManager<SurfaceMesh>::AllocateSharedPtr(i,m_verbose,
                    m_cad->GetSurf(i), m_octree,
                    m_curvemeshes,m_order);

            m_surfacemeshes[i]->Mesh(Nodes,Edges,Tris);

        }

        if(m_verbose)
        {
            cout << endl << "Surface mesh statistics" << endl;
            cout << "\tNodes: " << Nodes.size() << endl;
            cout << "\tTriangles " << Tris.size() << endl;
        }

    }


}
}
