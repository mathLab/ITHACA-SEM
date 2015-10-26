///////////////////////////////////////////////////////////////////////////////
//
// File Interpolator.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_INTERPOLATOR_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_INTERPOLATOR_H

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using namespace std;

namespace Nektar
{
namespace LibUtilities
{










    template<std::size_t DIM>
        class PtsPoint : public bg::model::point<NekDouble, DIM,  bg::cs::cartesian>{
            public:

                int                                         idx;
                Array<OneD, NekDouble>                      coords;
                NekDouble                                   dist;

                LIB_UTILITIES_EXPORT PtsPoint():
                    idx(-1),
                    coords(Array<OneD, NekDouble>(DIM)),
                    dist(1E30)
                {
                };

                LIB_UTILITIES_EXPORT PtsPoint(int idx, Array<OneD, NekDouble> coords,
                                            NekDouble dist):
                    idx(idx),
                    coords(coords),
                    dist(dist)
                {
                };

                LIB_UTILITIES_EXPORT bool operator < (const PtsPoint &comp) const
                {
                    return (dist < comp.dist);
                };
        };

        typedef bg::model::box<PtsPoint<3> > PtsBox;















enum PtsInterpMethod{
    ePtsNoMethod,
    ePtsNearestNeighbour,
    ePtsQuadratic,
    ePtsShepard,
    ePtsGauss,
};

template<std::size_t DIM>
class Interpolator
{
    public:

        LIB_UTILITIES_EXPORT Interpolator(const PtsField &inField, PtsField &outField) :
            m_inField(inField),
            m_outField(outField),
            m_method(ePtsNoMethod)
        {
        };

        LIB_UTILITIES_EXPORT void CalcWeights(
            PtsInterpMethod method,
            short int coordId = -1,
            NekDouble width = 0.0);

        LIB_UTILITIES_EXPORT void Interpolate(
            PtsInterpMethod method = ePtsNoMethod,
            short int coordId = -1,
            NekDouble width = 0.0);

        LIB_UTILITIES_EXPORT int GetDim() const;

        LIB_UTILITIES_EXPORT void GetInField(PtsField &inField) const;

        LIB_UTILITIES_EXPORT void GetOutField(PtsField &outField) const;

        LIB_UTILITIES_EXPORT void GetWeights(
            Array<OneD, Array<OneD, float> > &weights,
            Array<OneD, Array<OneD, unsigned int> > &neighbourInds) const;

        LIB_UTILITIES_EXPORT void SetWeights
        (const Array<OneD, Array<OneD, float> >
         &weights,
         const Array<OneD, Array<OneD, unsigned int> > &neighbourInds);


        template<typename FuncPointerT, typename ObjectPointerT>
        void setProgressCallback(FuncPointerT func,
                ObjectPointerT obj)
        {
            m_progressCallback = boost::bind(func, obj, _1, _2);
        }

    private:

        PtsField                                    m_inField;
        PtsField                                    m_outField;

        /// Interpolation Method
        PtsInterpMethod                             m_method;
        /// A tree structure to speed up the neighbour search.
        /// Note that we fill it with an iterator, so instead of rstar, the
        /// packing algorithm is used.
        bgi::rtree< PtsPoint<DIM>, bgi::rstar<16> > m_rtree;
        /// Interpolation weights for each neighbour.
        /// Structure: m_weights[physPtIdx][neighbourIdx]
        Array<OneD, Array<OneD, float> >            m_weights;
        /// Indices of the relevant neighbours for each physical point.
        /// Structure: m_neighInds[ptIdx][neighbourIdx]
        Array<OneD, Array<OneD, unsigned int> >     m_neighInds;

        boost::function<void (const int position, const int goal)> m_progressCallback;

        LIB_UTILITIES_EXPORT void CalcW_Gauss(
            const PtsPoint<DIM> &searchPt,
            const NekDouble sigma);

        LIB_UTILITIES_EXPORT void CalcW_Linear(const PtsPoint<DIM> &searchPt, int coordId);

        LIB_UTILITIES_EXPORT void CalcW_NNeighbour(const PtsPoint<DIM> &searchPt, int coordId);

        LIB_UTILITIES_EXPORT void CalcW_Shepard(const PtsPoint<DIM> &searchPt);

        LIB_UTILITIES_EXPORT void CalcW_Quadratic(const PtsPoint<DIM> &searchPt, int coordId);

        LIB_UTILITIES_EXPORT void FindNeighbours(const PtsPoint<DIM> &searchPt,
                vector<PtsPoint<DIM> > &neighbourPts,
                const NekDouble dist);

        LIB_UTILITIES_EXPORT void FindNNeighbours(const PtsPoint<DIM> &searchPt,
                vector<PtsPoint<DIM> > &neighbourPts,
                const unsigned int numPts = 1);
};

// typedef boost::shared_ptr<Interpolator> InterpolatorSharedPtr;
// static InterpolatorSharedPtr NullInterpolator;


}
}

// no idea why boost doesnt have a macro for 1d...
namespace boost { namespace geometry { namespace traits {
    BOOST_GEOMETRY_DETAIL_SPECIALIZE_POINT_TRAITS(Nektar::LibUtilities::PtsPoint<1>, 1, Nektar::NekDouble, cs::cartesian)
    BOOST_GEOMETRY_DETAIL_SPECIALIZE_POINT_ACCESS(Nektar::LibUtilities::PtsPoint<1>, 0, Nektar::NekDouble, coords[0], coords[0])
}}}
BOOST_GEOMETRY_REGISTER_POINT_2D(Nektar::LibUtilities::PtsPoint<2>, Nektar::NekDouble, cs::cartesian, coords[0], coords[1])
BOOST_GEOMETRY_REGISTER_POINT_3D(Nektar::LibUtilities::PtsPoint<3>, Nektar::NekDouble, cs::cartesian, coords[0], coords[1], coords[2])

#endif

