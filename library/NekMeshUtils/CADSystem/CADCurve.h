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

#ifndef NEKMESHUTILS_CADSYSTEM_CADCURVE
#define NEKMESHUTILS_CADSYSTEM_CADCURVE

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <NekMeshUtils/CADSystem/OpenCascade.h>

#include <NekMeshUtils/CADSystem/CADObj.h>
#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for CAD curves.
 *
 * This class wraps the OpenCascade BRepAdaptor_Curve class for use with
 * Nektar++.
 */
class CADCurve : public CADObj
{
public:
    friend class MemoryManager<CADCurve>;

    /**
     * @brief Default constructor.
     */
    CADCurve(int i, TopoDS_Shape in);

    ~CADCurve(){};

    /**
     * @brief Returns the minimum and maximum parametric coords t of the curve.
     *
     * @return Array of two entries, min and max parametric coordinate.
     */
    Array<OneD, NekDouble> Bounds();

    /**
     * @brief Calculates the arclength between the two paremetric points \p ti
     * and \p tf. \p ti must be less than \p tf.
     *
     * @param ti First parametric coordinate.
     * @param tf Second parametric coordinate.
     * @return Arc length between \p ti and \p tf.
     */
    NekDouble Length(NekDouble ti, NekDouble tf);

    /**
     * @brief Gets the location (x,y,z) in an array out of the curve at
     * point \p t.
     *
     * @param t Parametric coordinate
     * @return Array of x,y,z
     */
    Array<OneD, NekDouble> P(NekDouble t);

    /**
     * @brief Gets the second derivatives at t
     */
    Array<OneD, NekDouble> D2(NekDouble t);

    /**
     * @brief Calculates the parametric coordinate and arclength location
     * defined by \p s.
     *
     * @param s Arclength location.
     * @return Calculated parametric coordinate.
     *
     * @todo This really needs improving for accuracy.
     */
    NekDouble tAtArcLength(NekDouble s);

    /**
     * @brief Gets the start and end of the curve.
     *
     * @return Array with 6 entries of endpoints x1,y1,z1,x2,y2,z2.
     */
    Array<OneD, NekDouble> GetMinMax();

    /// set the ids of the surfaces either side of the curve
    void SetAdjSurf(std::vector<CADSurfSharedPtr> i)
    {
        m_adjSurfs = i;
    }

    /// returns the ids of neigbouring surfaces
    std::vector<CADSurfSharedPtr> GetAdjSurf()
    {
        return m_adjSurfs;
    }

    /// returns lenght of the curve
    NekDouble GetTotLength()
    {
        return m_length;
    }

    /*
     * @brief assign ids of end vertices in main cad
     */
    void SetVert(std::vector<CADVertSharedPtr> &falVert)
    {
        m_mainVerts = falVert;
    }

    /// get the ids of the vertices that are the ends of the curve,
    /// which are in the main cad list
    std::vector<CADVertSharedPtr> GetVertex()
    {
        return m_mainVerts;
    }

private:
    /// OpenCascade object of the curve.
    BRepAdaptor_Curve m_occCurve;
    /// OpenCascade edge
    TopoDS_Edge m_occEdge;
    /// Length of edge
    NekDouble m_length;
    /// List of surfaces which this curve belongs to.
    std::vector<CADSurfSharedPtr> m_adjSurfs;
    /// list of end vertices
    std::vector<CADVertSharedPtr> m_mainVerts;
};

typedef boost::shared_ptr<CADCurve> CADCurveSharedPtr;
}
}

#endif
