////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.h
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
//  Description: class for tet meshing
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_TETMESHING_TETMESH_H
#define NEKTAR_MESHUTILS_TETMESHING_TETMESH_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/MeshElements/Mesh.h>

#include <NekMeshUtils/ExtLibInterface/TetGenInterface.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for taking surface mesh and octree spec and making a 3d tet mesh
 */
class TetMesh
{
public:
    friend class MemoryManager<TetMesh>;

    /**
     * @brief default constructor
     */
    TetMesh(MeshSharedPtr m, int id,
            std::vector<ElementSharedPtr> e = std::vector<ElementSharedPtr>())
                : m_mesh(m), m_surface(e), m_id(id)
    {
    };

    /**
     *@brief execute tet meshing
     */
    void Mesh();

private:

    MeshSharedPtr m_mesh;
    std::vector<ElementSharedPtr> m_surface;
    /// number of tetrahedra
    int m_numtet;
    int m_id;
    /// conncetivity of the tets from the interface
    std::vector<Array<OneD, int> > m_tetconnect;

    TetGenInterfaceSharedPtr tetgen;

};

typedef std::shared_ptr<TetMesh> TetMeshSharedPtr;

}
}

#endif
