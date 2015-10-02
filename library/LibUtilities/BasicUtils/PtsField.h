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
#include <LibUtilities/BasicUtils/VmathArray.hpp>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{

enum PtsType{
    ePtsFile,
    ePtsLine,
    ePtsPlane,
    ePtsTriBlock,
    ePtsTetBlock
};


class PtsPoint
{
    public:

        int                                         m_idx;
        Array<OneD, NekDouble>                      m_coords;
        NekDouble                                   m_distSq;

        LIB_UTILITIES_EXPORT PtsPoint() {};

        LIB_UTILITIES_EXPORT PtsPoint(int idx, Array<OneD, NekDouble> coords,
                                      NekDouble distSq):
            m_idx(idx),
            m_coords(coords),
            m_distSq(distSq)
        {
        };

        LIB_UTILITIES_EXPORT bool operator < (const PtsPoint &comp) const
        {
            return (m_distSq < comp.m_distSq);
        };
};


class PtsField
{
    public:

        LIB_UTILITIES_EXPORT PtsField(
            const int dim,
            const Array<OneD, Array<OneD, NekDouble> > &pts) :
            m_dim(dim),
            m_pts(pts),
            m_ptsType(ePtsFile)
        {
        };

        LIB_UTILITIES_EXPORT PtsField(
            const int dim,
            const vector<std::string> fieldnames,
            const Array<OneD, Array<OneD, NekDouble> > &pts) :
            m_dim(dim),
            m_fieldNames(fieldnames),
            m_pts(pts),
            m_ptsType(ePtsFile)
        {
        };

        LIB_UTILITIES_EXPORT PtsField(
            const int dim,
            const vector<std::string> fieldnames,
            const Array<OneD, Array<OneD, NekDouble> > &pts,
            const Array<OneD, Array<OneD, float> > &weights,
            const Array<OneD, Array<OneD, unsigned int> > &neighInds) :
            m_dim(dim),
            m_fieldNames(fieldnames),
            m_pts(pts),
            m_ptsType(ePtsFile),
            m_weights(weights),
            m_neighInds(neighInds)
        {
        };

        LIB_UTILITIES_EXPORT void CalcWeights(
            const Array< OneD, Array< OneD, NekDouble > > &physCoords,
            short int coordId = -1);

        LIB_UTILITIES_EXPORT void Interpolate(
            const Array< OneD, Array< OneD, NekDouble > > &physCoords,
            Array<OneD, Array<OneD, NekDouble> > &intFields,
            short int coordId = -1);

        LIB_UTILITIES_EXPORT void Interpolate(
            Array<OneD, Array<OneD, NekDouble> >
            &intFields);

        LIB_UTILITIES_EXPORT void SetWeights
        (const Array<OneD, Array<OneD, float> >
         &weights,
         const Array<OneD, Array<OneD, unsigned int> > &neighbourInds);

        LIB_UTILITIES_EXPORT void GetWeights(
            Array<OneD, Array<OneD, float> > &weights,
            Array<OneD, Array<OneD, unsigned int> > &neighbourInds) const;

        LIB_UTILITIES_EXPORT void GetConnectivity(vector<Array<OneD, int> > &conn) const;

        LIB_UTILITIES_EXPORT void SetConnectivity(const vector<Array<OneD, int> > &conn);

        LIB_UTILITIES_EXPORT void SetDim(const int ptsDim);

        LIB_UTILITIES_EXPORT int GetDim() const;

        LIB_UTILITIES_EXPORT int GetNFields() const;

        LIB_UTILITIES_EXPORT vector<std::string> GetFieldNames() const;

        LIB_UTILITIES_EXPORT std::string GetFieldName(const int i) const;

        LIB_UTILITIES_EXPORT void SetFieldNames(
            const vector<std::string> fieldNames);

        LIB_UTILITIES_EXPORT void AddField(const Array<OneD, NekDouble> &pts,
                                           const std::string fieldName);

        LIB_UTILITIES_EXPORT int GetNpoints() const;

        LIB_UTILITIES_EXPORT NekDouble GetPointVal(const int fieldInd, const int ptInd) const;

        LIB_UTILITIES_EXPORT void GetPts(
            Array<OneD, Array<OneD, NekDouble> > &pts) const;

        LIB_UTILITIES_EXPORT Array<OneD, NekDouble> GetPts(const int fieldInd) const;

        LIB_UTILITIES_EXPORT void SetPts(Array<OneD, Array<OneD, NekDouble> > &pts);

        LIB_UTILITIES_EXPORT vector<int> GetPointsPerEdge() const;

        LIB_UTILITIES_EXPORT int GetPointsPerEdge(const int i) const;

        LIB_UTILITIES_EXPORT void SetPointsPerEdge(const vector<int> nPtsPerEdge);

        LIB_UTILITIES_EXPORT PtsType GetPtsType() const;

        LIB_UTILITIES_EXPORT void SetPtsType(const PtsType type);

        template<typename FuncPointerT, typename ObjectPointerT>
        void setProgressCallback(FuncPointerT func,
                ObjectPointerT obj)
        {
            m_progressCallback = boost::bind(func, obj, _1, _2);
        }

    private:

        /// Dimension of the pts field
        int                                     m_dim;
        /// Names of the field variables
        vector<std::string>                     m_fieldNames;
        /// Point data. For a n-dimensional field, the first n fields are the
        /// points spatial coordinates. Structure: m_pts[fieldIdx][ptIdx]
        Array<OneD, Array<OneD, NekDouble> >    m_pts;
        /// Number of points per edge. Empty if the point data has no
        /// specific shape (ePtsLine) or is a block (ePtsTetBlock,
        /// ePtsTriBlock), size=1 for ePtsLine and 2 for a ePtsPlane
        vector<int>                             m_nPtsPerEdge;
        /// Connectivity data needed for ePtsTetBlock and ePtsTriBlock. For n
        /// Blocks with m elements each, m_ptsConn is a vector of n arrays with
        /// 3*m (ePtsTriBlock) or 4*m (ePtsTetBlock) entries.
        vector<Array<OneD, int> >               m_ptsConn;
        /// Type of the PtsField
        PtsType                                 m_ptsType;
        /// Interpolation weights for each neighbour.
        /// Structure: m_weights[physPtIdx][neighbourIdx]
        Array<OneD, Array<OneD, float> >        m_weights;
        /// Indices of the relevant neighbours for each physical point.
        /// Structure: m_neighInds[ptIdx][neighbourIdx]
        Array<OneD, Array<OneD, unsigned int> > m_neighInds;

        boost::function<void (const int position, const int goal)> m_progressCallback;

        LIB_UTILITIES_EXPORT void CalcW_Linear(const int physPtIdx,
                                               const NekDouble coord);

        LIB_UTILITIES_EXPORT void CalcW_Shepard(
            const int physPtIdx,
            const Array< OneD, NekDouble > &physPtCoords);

        LIB_UTILITIES_EXPORT void CalcW_Quadratic(const int physPtIdx,
                const NekDouble coord);

        LIB_UTILITIES_EXPORT NekDouble DistSq(
            const Array<OneD, NekDouble > &point1,
            const Array<OneD, NekDouble > &point2) const;

        LIB_UTILITIES_EXPORT void FindNeighbours(const Array<OneD, NekDouble>
                &physPtCoords,
                vector<PtsPoint> &neighbourPts,
                const unsigned int numPts = 1);

};

typedef boost::shared_ptr<PtsField> PtsFieldSharedPtr;
static PtsFieldSharedPtr NullPtsField;


}
}
#endif
