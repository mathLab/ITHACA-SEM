///////////////////////////////////////////////////////////////////////////////
//
// File Interpolator.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Interpolator
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_BASICUTILS_INTERPOLATOR_H
#define LIBUTILITIES_BASICUTILS_INTERPOLATOR_H

#include <vector>
#include <iostream>
#include <functional>
#include <memory>

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
namespace LibUtilities
{

enum InterpMethod
{
    eNoMethod,
    eNearestNeighbour,
    eQuadratic,
    eShepard,
    eGauss,
};

/// A class that contains algorithms for interpolation between pts fields,
/// expansions and different meshes
class Interpolator
{
public:
    /**
    * @brief Constructor of the Interpolator class
    *
    * @param method    interpolation method, defaults to a sensible value if not
    * set
    * @param coordId   coordinate id along which the interpolation should be
    * performed
    * @param filtWidth filter width, required by some algorithms such as eGauss
    * @param maxPts    limit number of considered points
    *
    * if method is not specified, the best algorithm is chosen autpomatically.
    *
    * If coordId is not specified, a full 1D/2D/3D interpolation is performed
    * without
    * collapsing any coordinate.
    *
    * filtWidth must be specified for the eGauss algorithm only.
    */
    Interpolator(InterpMethod method = eNoMethod,
                 short int coordId   = -1,
                 NekDouble filtWidth = 0.0,
                 int maxPts          = 1000)
        : m_method(method), m_filtWidth(filtWidth), m_maxPts(maxPts),
          m_coordId(coordId)
    {
    }

    /// Compute interpolation weights without doing any interpolation
    LIB_UTILITIES_EXPORT void CalcWeights(
        const LibUtilities::PtsFieldSharedPtr ptsInField,
        LibUtilities::PtsFieldSharedPtr &ptsOutField,
        bool reuseTree = false);

    /// Interpolate from a pts field to a pts field
    LIB_UTILITIES_EXPORT void Interpolate(
        const LibUtilities::PtsFieldSharedPtr ptsInField,
        LibUtilities::PtsFieldSharedPtr &ptsOutField);

    /// returns the dimension of the Interpolator.
    /// Should be higher than the dimensions of the interpolated fields
    LIB_UTILITIES_EXPORT int GetDim() const;

    /// Returns the filter width
    LIB_UTILITIES_EXPORT NekDouble GetFiltWidth() const;

    /// Returns the coordinate id along which the interpolation should be
    /// performed
    LIB_UTILITIES_EXPORT int GetCoordId() const;

    /// Returns the interpolation method used by this interpolator
    LIB_UTILITIES_EXPORT InterpMethod GetInterpMethod() const;

    /// Returns the input field
    LIB_UTILITIES_EXPORT LibUtilities::PtsFieldSharedPtr GetInField() const;

    /// Returns the output field
    LIB_UTILITIES_EXPORT LibUtilities::PtsFieldSharedPtr GetOutField() const;

    /// Returns if the weights have already been computed
    LIB_UTILITIES_EXPORT void PrintStatistics();

    /// sets a callback funtion which gets called every time the interpolation
    /// progresses
    template <typename FuncPointerT, typename ObjectPointerT>
    void SetProgressCallback(FuncPointerT func, ObjectPointerT obj)
    {
        m_progressCallback = std::bind(
            func, obj, std::placeholders::_1, std::placeholders::_2);
    }

protected:

    /// input field
    LibUtilities::PtsFieldSharedPtr m_ptsInField;
    /// output field
    LibUtilities::PtsFieldSharedPtr m_ptsOutField;

    std::function<void(const int position, const int goal)>
    m_progressCallback;

private:

    class PtsPoint
    {
    public:
        int idx;
        Array<OneD, NekDouble> coords;
        NekDouble dist;

        PtsPoint() : idx(-1), coords(Array<OneD, NekDouble>(3)), dist(1E30){};

        PtsPoint(int idx, Array<OneD, NekDouble> coords, NekDouble dist)
            : idx(idx), coords(coords), dist(dist){};

        bool operator<(const PtsPoint &comp) const
        {
            return (dist < comp.dist);
        };
    };

    /// dimension of this interpolator. Hardcoded to 3
    static const int m_dim = 3;
    typedef boost::geometry::model::point<NekDouble, m_dim, boost::geometry::cs::cartesian> BPoint;
    typedef std::pair<BPoint, unsigned int> PtsPointPair;
    typedef boost::geometry::index::rtree<PtsPointPair, boost::geometry::index::rstar<16> > PtsRtree;



    /// Interpolation Method
    InterpMethod m_method;
    /// A tree structure to speed up the neighbour search.
    /// Note that we fill it with an iterator, so instead of rstar, the
    /// packing algorithm is used.
    std::shared_ptr<PtsRtree> m_rtree;
    /// Interpolation weights for each neighbour.
    /// Structure: m_weights[physPtIdx][neighbourIdx]
    Array<TwoD, NekDouble> m_weights;
    /// Indices of the relevant neighbours for each physical point.
    /// Structure: m_neighInds[ptIdx][neighbourIdx]
    Array<TwoD, unsigned int> m_neighInds;
    /// Filter width used for some interpolation algorithms
    NekDouble m_filtWidth;
    /// Max number of interpolation points
    int m_maxPts;
    /// coordinate id along which the interpolation should be performed
    short int m_coordId;

    LIB_UTILITIES_EXPORT void CalcW_Gauss(const PtsPoint &searchPt,
                                          const NekDouble sigma,
                                          const int maxPts = 250);

    LIB_UTILITIES_EXPORT void CalcW_Linear(const PtsPoint &searchPt, int coordId);

    LIB_UTILITIES_EXPORT void CalcW_NNeighbour(const PtsPoint &searchPt);

    LIB_UTILITIES_EXPORT void CalcW_Shepard(const PtsPoint &searchPt, int numPts);

    LIB_UTILITIES_EXPORT void CalcW_Quadratic(const PtsPoint &searchPt,
                                              int coordId);

    LIB_UTILITIES_EXPORT void SetupTree();

    LIB_UTILITIES_EXPORT void FindNeighbours(const PtsPoint &searchPt,
                                             std::vector<PtsPoint> &neighbourPts,
                                             const NekDouble dist,
                                             const unsigned int maxPts = 1);

    LIB_UTILITIES_EXPORT void FindNNeighbours(const PtsPoint &searchPt,
                                              std::vector<PtsPoint> &neighbourPts,
                                              const unsigned int numPts = 1);
};

typedef std::shared_ptr<Interpolator> InterpolatorSharedPtr;

}
}

#endif
