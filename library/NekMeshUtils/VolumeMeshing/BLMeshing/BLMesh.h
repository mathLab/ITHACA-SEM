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
    BLMesh(MeshSharedPtr m, std::vector<unsigned int> bls, NekDouble b, int l,
           NekDouble p, int id)
        : m_mesh(m), m_blsurfs(bls), m_bl(b), m_prog(p), m_layer(l), m_id(id)
    {
        
    };

    /**
     * @brief Execute boundary layer meshing
     */
    void Mesh();

    std::vector<unsigned int> GetSymSurfs()
    {
        return m_symSurfs;
    }
    std::vector<unsigned int> GetBLSurfs()
    {
        return m_blsurfs;
    }

    std::map<NodeSharedPtr, NodeSharedPtr> GetSymNodes();

    std::vector<ElementSharedPtr> GetPseudoSurface()
    {
        return m_psuedoSurface;
    }

    struct blInfo
    {
        NodeSharedPtr pNode;
        NodeSharedPtr oNode;
        int bl;
        Array<OneD, NekDouble> N;
        int symsurf;
        bool onSym;
        std::vector<ElementSharedPtr> els;
        std::set<int> surfs;

        bool stopped;

        void AlignNode(NekDouble t)
        {
            pNode->m_x = oNode->m_x + t * N[0];
            pNode->m_y = oNode->m_y + t * N[1];
            pNode->m_z = oNode->m_z + t * N[2];
        }
    };
    typedef std::shared_ptr<blInfo> blInfoSharedPtr;

private:
    void Setup();
    void GrowLayers();
    void Shrink();
    void BuildElements();
    bool TestIntersectionEl(ElementSharedPtr e1, ElementSharedPtr e2);
    bool IsPrismValid(ElementSharedPtr el);
    NekDouble Proximity(NodeSharedPtr n, ElementSharedPtr el);

    NekDouble Visability(std::vector<ElementSharedPtr> tris,
                         Array<OneD, NekDouble> N);
    Array<OneD, NekDouble> GetNormal(std::vector<ElementSharedPtr> tris);

    /// mesh object containing surface mesh
    MeshSharedPtr m_mesh;
    /// List of surfaces onto which boundary layers are placed
    std::vector<unsigned int> m_blsurfs;
    /// thickness of the boundary layer
    NekDouble m_bl;
    NekDouble m_prog;
    int m_layer;
    int m_id;
    Array<OneD, NekDouble> m_layerT;
    /// list of surfaces to be remeshed due to the boundary layer
    std::vector<unsigned int> m_symSurfs;
    /// data structure used to store and develop bl information
    std::map<NodeSharedPtr, blInfoSharedPtr> m_blData;
    std::map<NodeSharedPtr, std::vector<blInfoSharedPtr> >
        m_nToNInfo; // node to neighbouring information
    std::map<ElementSharedPtr, ElementSharedPtr> m_priToTri;
    std::vector<ElementSharedPtr> m_psuedoSurface;
    NekMatrix<NekDouble> m_deriv[3];
};

typedef std::shared_ptr<BLMesh> BLMeshSharedPtr;
}
}

#endif
