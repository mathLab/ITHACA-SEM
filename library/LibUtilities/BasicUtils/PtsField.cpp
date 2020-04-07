////////////////////////////////////////////////////////////////////////////////
//
// File: PtsField.cpp
//
// For more information, please see: http://www.nektar.info/
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
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/PtsField.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{

PtsField::PtsField(const int dim,
                   const Array<OneD, Array<OneD, NekDouble> > &pts)
    : m_dim(dim), m_pts(pts), m_ptsType(ePtsFile)
{
    for (int i = 0; i < GetNFields(); ++i)
    {
        m_fieldNames.push_back("NA");
    }
}

/**
 * @brief Set the connectivity data for ePtsTetBlock and ePtsTriBlock
 *
 * @param conn Connectivity data
 * Connectivity data needed for ePtsTetBlock and ePtsTriBlock. For n Blocks with
 * m elements each, m_ptsConn is a vector of n arrays with 3*m (ePtsTriBlock) or
 * 4*m (ePtsTetBlock) entries.
 */
void PtsField::GetConnectivity(vector< Array< OneD, int > > &conn) const
{
    conn = m_ptsConn;
}

/**
 * @brief Get the connectivity data for ePtsTetBlock and ePtsTriBlock
 *
 * @param conn Connectivity data
 * Connectivity data needed for ePtsTetBlock and ePtsTriBlock. For n Blocks with
 * m elements each, m_ptsConn is a vector of n arrays with 3*m (ePtsTriBlock) or
 * 4*m (ePtsTetBlock) entries.
 */
void PtsField::SetConnectivity(const vector< Array< OneD, int > > &conn)
{
    ASSERTL1((m_ptsType == ePtsTetBlock || m_ptsType == ePtsTriBlock ||
              m_ptsType == ePtsSegBlock),
             "ptsType must be set before connectivity");

    m_ptsConn = conn;
}

void PtsField::SetDim(const int ptsDim)
{
    m_dim = ptsDim;
}

size_t PtsField::GetDim() const
{
    return m_dim;
}

size_t PtsField::GetNFields() const
{
    return m_pts.size() - m_dim;
}

vector<std::string> PtsField::GetFieldNames() const
{
    return m_fieldNames;
}

std::string PtsField::GetFieldName(const int i) const
{
    return m_fieldNames[i];
}

void PtsField::SetFieldNames(const vector<std::string> fieldNames)
{
    ASSERTL0(fieldNames.size() == m_pts.size() - m_dim,
             "Number of given fieldNames does not match the number of stored "
             "fields");

    m_fieldNames = fieldNames;
}

void PtsField::AddField(const Array< OneD, NekDouble > &pts,
                        const string fieldName)
{
    size_t nTotvars = m_pts.size();

    ASSERTL1(pts.size() ==  m_pts[0].size(),
            "Field size mismatch");

    // redirect existing pts
    Array<OneD, Array<OneD, NekDouble> > newpts(nTotvars + 1);
    for (size_t i = 0; i < nTotvars; ++i)
    {
        newpts[i] = m_pts[i];
    }
    newpts[nTotvars] = pts;

    m_pts = newpts;

    m_fieldNames.push_back(fieldName);
}


void PtsField::RemoveField(const string fieldName)
{
    size_t nTotvars = m_pts.size();

    // redirect existing pts
    Array<OneD, Array<OneD, NekDouble> > newpts(nTotvars - 1);
    for (size_t i = 0, j = 0; i < nTotvars; ++i)
    {
        if (i < GetDim() || m_fieldNames[i - GetDim()] != fieldName)
        {
            newpts[j++] = m_pts[i];
        }
    }

    m_pts = newpts;

    m_fieldNames.erase(remove(m_fieldNames.begin(), m_fieldNames.end(), fieldName), m_fieldNames.end());
}

void PtsField::AddPoints(const Array<OneD, const Array<OneD, NekDouble> > &pts)
{
    ASSERTL1(pts.size() == m_pts.size(),
             "number of variables mismatch");

    // TODO: dont copy, dont iterate
    for (size_t i = 0; i < m_pts.size(); ++i)
    {
        Array<OneD, NekDouble> tmp(m_pts[i].size() + pts[i].size());
        for (size_t j = 0; j < m_pts[i].size(); ++j)
        {
            tmp[j] = m_pts[i][j];
        }
        for (size_t j = 0; j < pts[i].size(); ++j)
        {
            tmp[m_pts[i].size() + j] = pts[i][j];
        }
        m_pts[i] = tmp;
    }
}

size_t PtsField::GetNpoints() const
{
    return m_pts[0].size();
}

NekDouble PtsField::GetPointVal(const size_t fieldInd, const size_t ptInd) const
{
    return m_pts[fieldInd][ptInd];
}

void PtsField::SetPointVal(const size_t fieldInd,
                           const size_t ptInd,
                           const NekDouble val)
{
    m_pts[fieldInd][ptInd] = val;
}

void PtsField::GetPts(Array< OneD, Array< OneD, NekDouble > > &pts) const
{
    pts = m_pts;
}

Array< OneD, NekDouble > PtsField::GetPts(const int fieldInd) const
{
    return m_pts[fieldInd];
}

void PtsField::SetPts(Array< OneD, Array< OneD, NekDouble > > &pts)
{
    ASSERTL1(pts.size() ==  m_pts.size(),
             "Pts field count mismatch");

    m_pts = pts;
}

vector<size_t> PtsField::GetPointsPerEdge() const
{
    return m_nPtsPerEdge;
}

size_t PtsField::GetPointsPerEdge(const size_t i) const
{
    return m_nPtsPerEdge[i];
}

/**
 * @brief Set the number of points per edge
 *
 * @param nPtsPerEdge Number of points per edge. Empty if the point
 * data has no specific shape (ePtsLine) or is a block (ePtsTetBlock,
 * ePtsTriBlock), size=1 for ePtsLine, 2 for ePtsPlane and 3 for ePtsBox
 */
void PtsField::SetPointsPerEdge(const vector< size_t > nPtsPerEdge)
{
    ASSERTL0(
        m_ptsType == ePtsLine || m_ptsType == ePtsPlane || m_ptsType == ePtsBox,
             "SetPointsPerEdge only supported for ePtsLine, ePtsPlane and ePtsBox.");

    m_nPtsPerEdge = nPtsPerEdge;
}

PtsType PtsField::GetPtsType() const
{
    return m_ptsType;
}

void PtsField::SetPtsType(const PtsType type)
{
    m_ptsType = type;
}

vector<NekDouble> PtsField::GetBoxSize() const
{
    return m_boxSize;
}

void PtsField::SetBoxSize(const vector< NekDouble> boxSize)
{
    m_boxSize = boxSize;
}
}
}
