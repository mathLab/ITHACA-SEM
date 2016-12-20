////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.h
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
//  Description: CAD object curve.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CADSYSTEM_CADVERT
#define NEKMESHUTILS_CADSYSTEM_CADVERT

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/CADSystem/CADObject.h>
#include <NekMeshUtils/MeshElements/Node.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief base class for CAD verticies.
 *
 */
class CADVert : public CADObject
{
public:
    friend class MemoryManager<CADVert>;

    /**
     * @brief Default constructor.
     */
    CADVert()
    {
        m_type = CADType::eVert;
    }

    ~CADVert(){};

    /**
     * @brief Get x,y,z location of the vertex
     */
    Array<OneD, NekDouble> GetLoc()
    {
        Array<OneD, NekDouble> out(3);
        out[0] = m_node->m_x;
        out[1] = m_node->m_y;
        out[2] = m_node->m_z;
        return out;
    }

    /**
     * @brief returns a node object of the cad vertex
     */
    NodeSharedPtr GetNode()
    {
        return m_node;
    }

    /**
     * @brief if the vertex is degenerate manually set uv for that surface
     */
    void SetDegen(int s, CADSurfSharedPtr su, NekDouble u, NekDouble v)
    {
        degen     = true;
        degensurf = s;
        Array<OneD, NekDouble> uv(2);
        uv[0] = u;
        uv[1] = v;
        m_node->SetCADSurf(s, su, uv);
    }

    /**
     * @brief query is degenerate
     */
    int IsDegen()
    {
        if (degen)
        {
            return degensurf;
        }
        else
        {
            return -1;
        }
    }

protected:
    /// mesh convert object of vert
    NodeSharedPtr m_node;
    /// degen marker
    bool degen;
    /// degen surface
    int degensurf;
};

typedef boost::shared_ptr<CADVert> CADVertSharedPtr;

typedef LibUtilities::NekFactory<std::string, CADVert> CADVertFactory;

CADVertFactory& GetCADVertFactory();

}
}

#endif
