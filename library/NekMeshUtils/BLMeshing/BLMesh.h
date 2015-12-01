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
//  Description: class for tet meshing
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_BLMESHING_BLMESH_H
#define NEKTAR_MESHUTILS_BLMESHING_BLMESH_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/MeshElements/MeshElements.h>
#include <NekMeshUtils/CADSystem/CADSystem.h>

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
    BLMesh(MeshSharedPtr m, const std::vector<unsigned int> bl,
                            const std::vector<unsigned int> sym,
                            const NekDouble b)
                : m_mesh(m), m_blsurfs(bl), m_symsurfs(sym), m_bl(b)
    {
    };

    /**
     *@brief execute tet meshing
     */
    void Mesh();

private:

    /// print stuff to screen?
    MeshSharedPtr m_mesh;

    std::vector<unsigned int> m_blsurfs;
    std::vector<unsigned int> m_symsurfs;
    NekDouble m_bl;
};

typedef boost::shared_ptr<BLMesh> BLMeshSharedPtr;

}
}

#endif
