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
#include <memory>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

namespace Nektar
{
namespace LibUtilities
{

enum PtsType
{
    ePtsFile,
    ePtsLine,
    ePtsPlane,
    ePtsBox,
    ePtsSegBlock,
    ePtsTriBlock,
    ePtsTetBlock
};

enum PtsInfo
{
    eIsEquiSpacedData,
    ePtsPerElmtEdge
};

static std::map<PtsInfo, int> NullPtsInfoMap;

class PtsField
{
public:
    LIB_UTILITIES_EXPORT PtsField(
        const int dim, const Array<OneD, Array<OneD, NekDouble> > &pts);

    LIB_UTILITIES_EXPORT PtsField(
        const int dim,
        const std::vector<std::string> fieldnames,
        const Array<OneD, Array<OneD, NekDouble> > &pts,
        std::map<PtsInfo, int> ptsInfo = NullPtsInfoMap)
        : m_ptsInfo(ptsInfo), m_dim(dim), m_fieldNames(fieldnames), m_pts(pts),
          m_ptsType(ePtsFile){};

    LIB_UTILITIES_EXPORT PtsField(
        const int dim,
        const std::vector<std::string> fieldnames,
        const Array<OneD, Array<OneD, NekDouble> > &pts,
        const Array<OneD, Array<OneD, float> > & weights,
        const Array<OneD, Array<OneD, unsigned int> > & neighInds)
        : m_ptsInfo(NullPtsInfoMap), m_dim(dim), m_fieldNames(fieldnames),
          m_pts(pts), m_ptsType(ePtsFile) { boost::ignore_unused(weights, neighInds); };

    LIB_UTILITIES_EXPORT void GetConnectivity(
        std::vector<Array<OneD, int> > &conn) const;

    LIB_UTILITIES_EXPORT void SetConnectivity(
        const std::vector<Array<OneD, int> > &conn);

    LIB_UTILITIES_EXPORT void SetDim(const int ptsDim);

    LIB_UTILITIES_EXPORT size_t GetDim() const;

    LIB_UTILITIES_EXPORT size_t GetNFields() const;

    LIB_UTILITIES_EXPORT std::vector<std::string> GetFieldNames() const;

    LIB_UTILITIES_EXPORT std::string GetFieldName(const int i) const;

    LIB_UTILITIES_EXPORT void SetFieldNames(
        const std::vector<std::string> fieldNames);

    LIB_UTILITIES_EXPORT void AddField(const Array<OneD, NekDouble> &pts,
                                       const std::string fieldName);

    LIB_UTILITIES_EXPORT void RemoveField(const std::string fieldName);

    LIB_UTILITIES_EXPORT void AddPoints(const Array< OneD, const Array< OneD, NekDouble > > &pts);

    LIB_UTILITIES_EXPORT size_t GetNpoints() const;

    LIB_UTILITIES_EXPORT NekDouble GetPointVal(const size_t fieldInd,
                                               const size_t ptInd) const;

    LIB_UTILITIES_EXPORT void SetPointVal(const size_t fieldInd,
                                          const size_t ptInd,
                                          const NekDouble val);

    LIB_UTILITIES_EXPORT void GetPts(
        Array<OneD, Array<OneD, NekDouble> > &pts) const;

    LIB_UTILITIES_EXPORT Array<OneD, NekDouble> GetPts(
        const int fieldInd) const;

    LIB_UTILITIES_EXPORT void SetPts(Array<OneD, Array<OneD, NekDouble> > &pts);

    LIB_UTILITIES_EXPORT std::vector<size_t> GetPointsPerEdge() const;

    LIB_UTILITIES_EXPORT size_t GetPointsPerEdge(const size_t i) const;

    LIB_UTILITIES_EXPORT void SetPointsPerEdge(
        const std::vector<size_t> nPtsPerEdge);

    LIB_UTILITIES_EXPORT PtsType GetPtsType() const;

    LIB_UTILITIES_EXPORT void SetPtsType(const PtsType type);

    LIB_UTILITIES_EXPORT std::vector<NekDouble> GetBoxSize() const;

    LIB_UTILITIES_EXPORT void SetBoxSize(const std::vector<NekDouble> boxsize);

    /// map for information about points that can be added through PtsInfo enum
    std::map<PtsInfo, int> m_ptsInfo;

private:
    /// Dimension of the pts field
    size_t m_dim;
    /// Names of the field variables
    std::vector<std::string> m_fieldNames;
    /// Point data. For a n-dimensional field, the first m_dim fields are the
    /// points spatial coordinates. Structure: m_pts[fieldIdx][ptIdx]
    Array<OneD, Array<OneD, NekDouble> > m_pts;
    /// Number of points per edge. Empty if the point data has no
    /// specific shape (ePtsLine) or is a block (ePtsTetBlock,
    /// ePtsTriBlock), size=1 for ePtsLine and 2 for a ePtsPlane
    std::vector<size_t> m_nPtsPerEdge;
    /// Connectivity data needed for ePtsTetBlock and ePtsTriBlock. For n
    /// Blocks with m elements each, m_ptsConn is a vector of n arrays with
    /// 3*m (ePtsTriBlock) or 4*m (ePtsTetBlock) entries.
    std::vector<Array<OneD, int> > m_ptsConn;
    /// Type of the PtsField
    PtsType m_ptsType;

    /// vector of box size xmin,xmax,ymin,ymax,zmin,zmax
    std::vector<NekDouble> m_boxSize;
};

typedef std::shared_ptr<PtsField> PtsFieldSharedPtr;
static PtsFieldSharedPtr NullPtsField;
}
}

#endif
