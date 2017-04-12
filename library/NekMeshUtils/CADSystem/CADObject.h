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

namespace Nektar
{
namespace NekMeshUtils
{

namespace CADType
{
enum cadType
{
    eVert,
    eCurve,
    eSurf
};
}

namespace CADOrientation
{
enum Orientation
{
    eUnknown,
    eForwards,
    eBackwards
};
}

class CADObject
{
public:
    friend class MemoryManager<CADObject>;

    /**
     * @brief Default constructor.
     */
    CADObject()
    {
    }

    virtual ~CADObject()
    {
    }

    /**
     * @brief Return ID of the vertex
     */
    int GetId()
    {
        return m_id;
    }

    CADType::cadType GetType()
    {
        return m_type;
    }

    virtual CADOrientation::Orientation Orientation()
    {
        ASSERTL0(false,"must be implemented at the cad object level");
        return CADOrientation::eUnknown;
    }

protected:
    /// ID of the vert.
    int m_id;
    /// type of the cad object
    CADType::cadType m_type;
    /// orientation of the CADObject
    CADOrientation::Orientation m_orientation;
};

typedef boost::shared_ptr<CADObject> CADObjectSharedPtr;
}
}

#endif
