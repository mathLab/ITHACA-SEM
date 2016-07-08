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

struct blInfo
{
    NodeSharedPtr pNode;
    NodeSharedPtr oNode;
    int bl;
    Array<OneD, NekDouble> N;
    int symsurf;
    bool onSym;
};

class BLMesh
{
public:
    friend class MemoryManager<BLMesh>;

    /**
     *@brief default constructor
     */
    BLMesh(CADSystemSharedPtr s, MeshSharedPtr m, std::vector<unsigned int> bls,
           NekDouble b) :
                        m_cad(s), m_mesh(m), m_blsurfs(bls), m_bl(b)
    {
    };

    /**
     * @brief Execute boundary layer meshing
     */
    void Mesh();

    /**
     * @brief Get the map of surface element id to pseudo surface prism face
     */
    std::vector<ElementSharedPtr> GetPsuedoSurf()
    {
        return m_psuedoSurface;
    }

    /**
     * @brief Get the list of symetry surfaces
     */
    std::vector<int> GetSymSurfs()
    {
        return m_symSurfs;
    }

    /**
     * @brief Get the map of boundary layer surface nodes to the prism top nodes
     */
    std::map<NodeSharedPtr, NodeSharedPtr> GetNodeMap(int s)
    {
        std::map<int, std::map<NodeSharedPtr, NodeSharedPtr> >::iterator f;
        f = m_symNodes.find(s);
        ASSERTL0(f != m_symNodes.end(), "surf not found");
        return f->second;
    }

    /*
     * @brief Get the list of surfaces to have a boundary layer generated
     */
    std::vector<unsigned int> GetBLSurfs()
    {
        return m_blsurfs;
    }

private:

    NekDouble Visability(std::vector<ElementSharedPtr> tris, Array<OneD, NekDouble> N);
    Array<OneD, NekDouble> GetNormal(std::vector<ElementSharedPtr> tris);

    ///CAD
    CADSystemSharedPtr m_cad;
    /// mesh object containing surface mesh
    MeshSharedPtr m_mesh;
    /// List of surfaces onto which boundary layers are placed
    std::vector<unsigned int> m_blsurfs;
    /// thickness of the boundary layer
    NekDouble m_bl;
    /// list of surfaces to be remeshed due to the boundary layer
    std::vector<int> m_symSurfs;
    /// data structure used to store and develop bl information
    std::map<NodeSharedPtr, blInfo> blData;
    /// list of nodes which will lie of symtetry surfaces
    std::map<int, std::map<NodeSharedPtr, NodeSharedPtr> > m_symNodes;
    /// list of elements which form the psuedo surface from the top of prisms
    std::vector<ElementSharedPtr> m_psuedoSurface;
};

typedef boost::shared_ptr<BLMesh> BLMeshSharedPtr;
}
}

#endif
