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

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/NekFactory.hpp>

#include <NekMeshUtils/CADSystem/CADObject.h>

#include <NekMeshUtils/MeshElements/Node.h>

namespace Nektar
{
namespace NekMeshUtils
{

//forward decleration
class Node;
typedef std::shared_ptr<Node> NodeSharedPtr;
class CADSurf;
typedef std::shared_ptr<CADSurf> CADSurfSharedPtr;
class CADCurve;
typedef std::shared_ptr<CADCurve> CADCurveSharedPtr;

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

    virtual ~CADVert(){};

    /**
     * @brief Get x,y,z location of the vertex
     */
    Array<OneD, NekDouble> GetLoc();

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
    void SetDegen(int s, CADSurfSharedPtr su, NekDouble u, NekDouble v);

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

    /**
     * @brief Calcuate the distance to a vertex from a point l(x,y,z)
     */
    virtual NekDouble DistanceTo(Array<OneD, NekDouble> l) = 0;

    void AddAdjCurve(CADCurveSharedPtr c)
    {
        curves.push_back(c);
    }

    /**
     * @brief Get list of CAD curves which are bound by this vertex
     */
    std::vector<std::weak_ptr<CADCurve> > GetAdjCurves()
    {
        return curves;
    }

protected:
    /// mesh convert object of vert
    NodeSharedPtr m_node;
    /// degen marker
    bool degen;
    /// degen surface
    int degensurf;
    /// adjacent curves
    std::vector<std::weak_ptr<CADCurve> > curves;
};

typedef std::shared_ptr<CADVert> CADVertSharedPtr;

typedef LibUtilities::NekFactory<std::string, CADVert> CADVertFactory;

CADVertFactory& GetCADVertFactory();

}
}

#endif
