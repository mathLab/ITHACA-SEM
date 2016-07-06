////////////////////////////////////////////////////////////////////////////////
//
//  File: BLMesh.h
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
//  Description: class for boundary layer meshing
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_BLMESHING_BLMESH_H
#define NEKTAR_MESHUTILS_BLMESHING_BLMESH_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <NekMeshUtils/MeshElements/Mesh.h>

namespace Nektar
{
namespace NekMeshUtils
{

class BLMesh
{
public:
    friend class MemoryManager<BLMesh>;

    /**
     *@brief default constructor
     */
    BLMesh(MeshSharedPtr             m,
           std::vector<unsigned int> bls,
           std::vector<unsigned int> syms,
           NekDouble                 b)
        : m_mesh(m), m_blsurfs(bls), m_symsurfs(syms), m_bl(b){};

    /**
     * @brief Execute boundary layer meshing
     */
    void Mesh();

    /**
     * @brief Get the map of surface element id to pseudo surface prism face
     */
    std::map<int, FaceSharedPtr> GetSurfToPri()
    {
        return m_surftopriface;
    }

private:
    /// Mesh object containing surface mesh
    MeshSharedPtr m_mesh;
    /// List of surfaces onto which boundary layers are placed
    std::vector<unsigned int> m_blsurfs;
    /// List of symmetry surfaces
    std::vector<unsigned int> m_symsurfs;
    /// Thickness of the boundary layer
    NekDouble m_bl;
    /// Map from surface element ID to opposite face of prism
    std::map<int, FaceSharedPtr> m_surftopriface;
};

typedef boost::shared_ptr<BLMesh> BLMeshSharedPtr;
}
}

#endif
