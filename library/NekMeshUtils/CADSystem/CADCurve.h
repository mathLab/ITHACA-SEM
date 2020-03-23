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

#ifndef NEKMESHUTILS_CADSYSTEM_CADCURVE
#define NEKMESHUTILS_CADSYSTEM_CADCURVE

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <NekMeshUtils/CADSystem/CADObject.h>

namespace Nektar
{
namespace NekMeshUtils
{

// forward declarations
class CADVert;
typedef std::shared_ptr<CADVert> CADVertSharedPtr;
class CADSurf;
typedef std::shared_ptr<CADSurf> CADSurfSharedPtr;

/**
 * @brief base class for CAD curves.
 *
 */
class CADCurve : public CADObject
{
public:
    friend class MemoryManager<CADCurve>;

    /**
     * @brief Default constructor.
     */
    CADCurve()
    {
        m_type = CADType::eCurve;
    }

    virtual ~CADCurve()
    {
    }

    /**
     * @brief Returns the minimum and maximum parametric coords t of the curve.
     *
     * @return Array of two entries, min and max parametric coordinate.
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> GetBounds() = 0;

    /**
     * @brief Returns the minimum and maximum parametric coords t of the curve.
     */
    NEKMESHUTILS_EXPORT virtual void GetBounds(
        NekDouble &tmin, NekDouble &tmax) = 0;

    /**
     * @brief Calculates the arclength between the two paremetric points \p ti
     * and \p tf. \p ti must be less than \p tf.
     *
     * @param ti First parametric coordinate.
     * @param tf Second parametric coordinate.
     * @return Arc length between \p ti and \p tf.
     */
    NEKMESHUTILS_EXPORT virtual NekDouble Length(
        NekDouble ti, NekDouble tf) = 0;

    /**
     * @brief Gets the location (x,y,z) in an array out of the curve at
     * point \p t.
     *
     * @param t Parametric coordinate
     * @return Array of x,y,z
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> P(NekDouble t) = 0;

    /**
     * @brief Gets the location (x,y,z) in an array out of the curve at
     * point \p t.
     *
     * @param t Parametric coordinate
     */
    NEKMESHUTILS_EXPORT virtual void P(
        NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z) = 0;

    /**
     * @brief Gets the second derivatives at t
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> D2(NekDouble t) = 0;

    /**
     * @brief Calculates the radius of curvature of the curve at point t
     */
    NEKMESHUTILS_EXPORT virtual NekDouble Curvature(NekDouble t) = 0;

    /**
     * @brief Calculates the parametric coordinate and arclength location
     * defined by \p s.
     *
     * @param s Arclength location.
     * @return Calculated parametric coordinate.
     *
     * @todo This really needs improving for accuracy.
     */
    NEKMESHUTILS_EXPORT virtual NekDouble tAtArcLength(NekDouble s) = 0;

    /**
     * @brief Gets the start and end of the curve.
     *
     * @return Array with 6 entries of endpoints x1,y1,z1,x2,y2,z2.
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> GetMinMax() = 0;

    /**
     * @brief set the ids of the surfaces either side of the curve
     */
    void SetAdjSurf(std::pair<CADSurfSharedPtr, CADOrientation::Orientation> i)
    {
        m_adjSurfs.push_back(i);
    }

    /*
     * @brief returns the ids of surfaces bound by this curve as well as their
     *        Orientation with respect to the loop of curves
     */
    std::vector<std::pair<std::weak_ptr<CADSurf>, CADOrientation::Orientation>>
    GetAdjSurf()
    {
        return m_adjSurfs;
    }

    /*
     * @brief returns lenght of the curve
     */
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

    /*
     * @brief get the vertices that are the ends of the curve,
     * which are in the main cad list
     */
    std::vector<CADVertSharedPtr> GetVertex()
    {
        return m_mainVerts;
    }

    /*
     * @brief locates a point in the parametric space. returns the
     * distance to the point and passes t by reference and updates it
     */
    NEKMESHUTILS_EXPORT virtual NekDouble loct(
        Array<OneD, NekDouble> xyz, NekDouble &t) = 0;

    /**
     * @brief Returns the orientation of the curve with respect to a given
     * surface by id surf
     */
    NEKMESHUTILS_EXPORT CADOrientation::Orientation GetOrienationWRT(int surf);

    /**
     * @brief Returns the normal to the curve which is orientate with respect
     * to the surface surf
     */
    NEKMESHUTILS_EXPORT Array<OneD, NekDouble> NormalWRT(NekDouble t, int surf);

    /**
     * @brief Returns the normal to a curve, it will always point in the concave
     * direction
     */
    NEKMESHUTILS_EXPORT virtual Array<OneD, NekDouble> N(NekDouble t) = 0;

protected:
    /// Length of edge
    NekDouble m_length;
    /// List of surfaces which this curve belongs to.
    std::vector<std::pair<std::weak_ptr<CADSurf>, CADOrientation::Orientation>>
        m_adjSurfs;
    /// list of end vertices
    std::vector<CADVertSharedPtr> m_mainVerts;
};

typedef std::shared_ptr<CADCurve> CADCurveSharedPtr;

typedef LibUtilities::NekFactory<std::string, CADCurve> CADCurveFactory;

CADCurveFactory &GetCADCurveFactory();
}
}

#endif
