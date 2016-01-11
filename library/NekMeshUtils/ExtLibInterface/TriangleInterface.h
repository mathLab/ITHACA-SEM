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
//  Description: class for interfacing with triangle
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_MESHUTILS_EXTLIBINTERFACE_TRIANGLEINTERFACE_H
#define NEKTAR_MESHUTILS_EXTLIBINTERFACE_TRIANGLEINTERFACE_H

#include <boost/shared_ptr.hpp>

//horible definitions to get triangle to work
#define REAL double
#define ANSI_DECLARATORS
#define TRILIBRARY
#define VOID int

extern "C"{
#include <triangle.h>
}

#include <NekMeshUtils/MeshElements/MeshElements.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

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
    NEKMESHUTILS_EXPORT TriangleInterface()
    {
    };

    /**
     * @brief assign meshing paramters
     */
    NEKMESHUTILS_EXPORT void Assign(std::vector<std::vector<NodeSharedPtr> > &boundingloops,
                std::vector<Array<OneD, NekDouble> > &centers, int i,
                NekDouble str = 1.0)
    {
        m_boundingloops = boundingloops;
        m_centers = centers;
        m_str = str;
        sid = i;
    }

    NEKMESHUTILS_EXPORT void AssignStiener(std::vector<NodeSharedPtr> stiner)
    {
        m_stienerpoints = stiner;
    }

    /**
     * @brief execute meshing
     */
    NEKMESHUTILS_EXPORT void Mesh(bool Quiet = true, bool Quality = false);

    /**
     * @brief extract mesh
     */
    NEKMESHUTILS_EXPORT void Extract(std::vector<std::vector<NodeSharedPtr> > &Connec);

private:

    /**
     * @brief clear memory
     */
    void SetUp();

    /// list of bounding nodes to the surface
    std::vector<std::vector<NodeSharedPtr> > m_boundingloops;
    /// list of additional nodes
    std::vector<NodeSharedPtr> m_stienerpoints;
    /// coordinates of the centers of the loops
    std::vector<Array<OneD, NekDouble> > m_centers;
    /// map from NekMesh id to triangle id
    std::map<int, NodeSharedPtr> nodemap;
    /// id of the surface
    int sid;
    /// Stretching factor of parameter plane
    NekDouble m_str;
    /// triangle data strucutres
    struct triangulateio in,out;
};

typedef boost::shared_ptr<TriangleInterface> TriangleInterfaceSharedPtr;
}
}

#endif
