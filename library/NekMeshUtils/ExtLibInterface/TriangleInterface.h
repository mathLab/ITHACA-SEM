////////////////////////////////////////////////////////////////////////////////
//
//  File: TriangleInterface.h
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
//  Description: class for interfacing with triangle
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_EXTLIBINTERFACE_TRIANGLEINTERFACE_H
#define NEKTAR_MESHUTILS_EXTLIBINTERFACE_TRIANGLEINTERFACE_H

#include <memory>

#include <NekMeshUtils/MeshElements/Node.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C"
{
    #include <triangle.h>
}

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for interfacing with external library triangle
 */
class TriangleInterface
{
public:
    friend class MemoryManager<TriangleInterface>;

    /**
     * @brief default constructor
     */
    TriangleInterface(){};

    /**
     * @brief assign meshing paramters
     */
    void Assign(std::vector<std::vector<NodeSharedPtr> > &boundingloops,
                std::vector<Array<OneD, NekDouble> >     &centers,
                int                                       i,
                NekDouble                                 str = 1.0)
    {
        m_boundingloops = boundingloops;
        m_centers       = centers;
        m_str           = str;
        sid             = i;
    }

    void AssignStiener(std::vector<NodeSharedPtr> stiner)
    {
        m_stienerpoints = stiner;
    }

    /**
     * @brief Execute meshing
     */
    void Mesh(bool Quality = false);

    /**
     * @brief Extract mesh
     */
    void Extract(std::vector<std::vector<NodeSharedPtr> > &Connec);

private:
    /**
     * @brief Clear memory
     */
    void SetUp();

    struct DelaunayTriangle
    {
    public:
        void Run(char* cmd)
        {
            triangulate(cmd, &in, &out, NULL);
        }
        struct triangulateio in, out;
    };

    /// List of bounding nodes to the surface
    std::vector<std::vector<NodeSharedPtr> > m_boundingloops;
    /// List of additional nodes
    std::vector<NodeSharedPtr>               m_stienerpoints;
    /// Coordinates of the centers of the loops
    std::vector<Array<OneD, NekDouble> >     m_centers;
    /// Map from NekMesh id to triangle id
    std::map<int, NodeSharedPtr>             nodemap;
    /// ID of the surface
    int                                      sid;
    /// Stretching factor of parameter plane
    NekDouble                                m_str;
    /// Triangle data strucutres
    DelaunayTriangle dt;
};

typedef std::shared_ptr<TriangleInterface> TriangleInterfaceSharedPtr;
}
}

#endif
