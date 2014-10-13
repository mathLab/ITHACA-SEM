///////////////////////////////////////////////////////////////////////////////
//
// File PtsField.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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
// Description: Pts field
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_PTSFIELD_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_PTSFIELD_H

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{

class PtsPoint
{
    public:

        int                                         m_idx;
        Array<OneD, NekDouble>                      m_coords;
        NekDouble                                   m_distSq;

        PtsPoint() {};

        PtsPoint(int idx, Array<OneD, NekDouble> coords, NekDouble distSq):
            m_idx(idx),
            m_coords(coords),
            m_distSq(distSq)
        {
        };

        bool operator < (const PtsPoint &comp) const
        {
            return (m_distSq < comp.m_distSq);
        };
};


class PtsField
{
    public:

        PtsField(const int dim,
                 const Array<OneD, Array<OneD, NekDouble> > &pts) :
            m_dim(dim),
            m_pts(pts)
        {
        };

        PtsField(const int dim,
                 const vector<std::string> fieldnames,
                 const Array<OneD, Array<OneD, NekDouble> > &pts) :
            m_dim(dim),
            m_fieldNames(fieldnames),
            m_pts(pts)
        {
        };

        PtsField(const int dim,
                 const vector<std::string> fieldnames,
                 const Array<OneD, Array<OneD, NekDouble> > &pts,
                 const Array<OneD, Array<OneD, float> > &weights,
                 const Array<OneD, Array<OneD, unsigned int> > &neighInds) :
            m_dim(dim),
            m_fieldNames(fieldnames),
            m_pts(pts),
            m_weights(weights),
            m_neighInds(neighInds)
        {
        };

        void CalcWeights(
            const Array< OneD, Array< OneD, NekDouble > > &physCoords,
            short int coordId = -1);

        void Interpolate(
            const Array< OneD, Array< OneD, NekDouble > > &physCoords,
            Array<OneD, Array<OneD, NekDouble> > &intFields,
            short int coordId = -1);

        void Interpolate(Array<OneD, Array<OneD, NekDouble> > &intFields);

        void SetWeights(const Array<OneD, Array<OneD, float> > &weights,
                        const Array<OneD, Array<OneD, unsigned int> > &neighbourInds);

        void GetWeights(Array<OneD, Array<OneD, float> > &weights,
                        Array<OneD, Array<OneD, unsigned int> > &neighbourInds) const;

        void SetDim(const int ptsDim);

        int GetDim() const;

        int GetNFields() const;

        vector<std::string> GetFieldNames() const;

        std::string GetFieldName(const int i) const;

        void SetFieldNames(const vector<std::string> fieldName);

        void AddFieldName(const std::string fieldName);

        int GetNpoints() const;

        void GetPts(Array<OneD,  Array<OneD,  NekDouble> > &pts) const;

        void SetPoint(const int fieldIdx, const int pointIdx, NekDouble value);

        void SetPointsArray(Array<OneD,  Array<OneD,  NekDouble> > &pts);

        vector<int> GetPointsPerEdge() const;

        int GetPointsPerEdge(const int i) const;

        void SetPointsPerEdge(const vector<int> nPtsPerEdge);

        template<typename FuncPointerT, typename ObjectPointerT>
        void setProgressCallback(FuncPointerT func, ObjectPointerT obj)
        {
            m_progressCallback = boost::bind(func, obj, _1, _2);
        }

    private:

        /// Dimension of the pts field
        int                                     m_dim;
        /// Names of the field variables
        vector<std::string>                     m_fieldNames;
        /// Point data. For a n-dimensional field,  the first n fields are the
        /// points spatial coordinates. Structure: m_pts[fieldIdx][ptIdx]
        Array<OneD, Array<OneD, NekDouble> >    m_pts;
        /// Number of points per edge. Empty if the point data has no
        /// specific shape, size=1 for a line, 2 for a plane, 3 for a box.
        /// Parallel edges have the same number of points.
        vector<int>                             m_nPtsPerEdge;
        /// Interpolation weights for each neighbour.
        /// Structure: m_weights[physPtIdx][neighbourIdx]
        Array<OneD, Array<OneD, float> >        m_weights;
        /// Indices of the relevant neighbours for each physical point.
        /// Structure: m_neighInds[ptIdx][neighbourIdx]
        Array<OneD, Array<OneD, unsigned int> > m_neighInds;

        boost::function<void (const int position, const int goal)> m_progressCallback;

        void CalcW_Linear(const int physPtIdx, const NekDouble coord);

        void CalcW_Shepard(const int physPtIdx,
                           const Array< OneD, NekDouble > &physPtCoords);

        void CalcW_Quadratic(const int physPtIdx, const NekDouble coord);

        NekDouble DistSq(const Array<OneD, NekDouble > &point1,
                         const Array<OneD, NekDouble > &point2) const;

        void FindNeighbours(const Array<OneD, NekDouble> &physPtCoords,
                            vector<PtsPoint> &neighbourPts,
                            const unsigned int numPts = 1);

};

typedef boost::shared_ptr<PtsField> PtsFieldSharedPtr;
static PtsFieldSharedPtr NullPtsField;


}
}
#endif
