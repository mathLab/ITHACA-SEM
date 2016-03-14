////////////////////////////////////////////////////////////////////////////////
//
//  File: CADObj.h
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

#ifndef NEKMESHUTILS_CADSYSTEM_CADOBJ
#define NEKMESHUTILS_CADSYSTEM_CADOBJ

#include <boost/shared_ptr.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <NekMeshUtils/CADSystem/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

enum cadType
{
    vert,
    curve,
    surf
};

/**
 * @brief class for CAD curves.
 *
 * This class wraps the OpenCascade BRepAdaptor_Curve class for use with
 * Nektar++.
 */
class CADObj
{
public:
    friend class MemoryManager<CADObj>;

    /**
     * @brief Default constructor.
     */
    CADObj()
    {
    }

    virtual ~CADObj(){};

    /**
     * @brief Return ID of the vertex
     */
    int GetId()
    {
        return m_id;
    }

    cadType GetType()
    {
        return m_type;
    }

protected:
    /// ID of the vert.
    int m_id;
    /// type of the cad object
    cadType m_type;
};

typedef boost::shared_ptr<CADObj> CADObjSharedPtr;
}
}

#endif
