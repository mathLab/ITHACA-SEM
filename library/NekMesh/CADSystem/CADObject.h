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

#ifndef NEKMESH_CADSYSTEM_CADOBJ
#define NEKMESH_CADSYSTEM_CADOBJ

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <NekMesh/NekMeshDeclspec.h>
#include <NekMesh/Module/Log.hpp>

namespace Nektar
{
namespace NekMesh
{

namespace CADType
{
enum cadType
{
    eVert,
    eCurve,
    eSurf,
    eOther
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
     * @brief Return ID of the CAD object
     */
    int GetId()
    {
        return m_id;
    }

    /**
     * @brief Get the type of the CAD object
     */
    CADType::cadType GetType()
    {
        return m_type;
    }

    /**
     * @brief Get the Orientation of the CAD object
     */
    virtual CADOrientation::Orientation Orientation()
    {
        ASSERTL0(false, "must be implemented at the CAD object level");
        return CADOrientation::eUnknown;
    }

    /**
     * @brief Give the CAD object a string name
     */
    void SetName(std::string i)
    {
        m_name = i;
    }

    /**
     * @brief Get the name of a CAD object
     */
    std::string GetName()
    {
        return m_name;
    }

    /**
     * @brief Set the logger for this CAD object.
     */
    void SetLogger(Logger &log)
    {
        m_log = log;

        if (m_type == CADType::eVert)
        {
            m_log.SetPrefix("CADVert");
        }
        else if (m_type == CADType::eCurve)
        {
            m_log.SetPrefix("CADCurve");
        }
        else if (m_type == CADType::eSurf)
        {
            m_log.SetPrefix("CADSurf");
        }
    }

protected:
    /// ID of the vert.
    int m_id;
    /// type of the cad object
    CADType::cadType m_type;
    /// orientation of the CADObject
    CADOrientation::Orientation m_orientation;
    /// string name of the cad
    std::string m_name;
    /// Logger
    Logger m_log;
};

typedef std::shared_ptr<CADObject> CADObjectSharedPtr;
}
}

#endif
