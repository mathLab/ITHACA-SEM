////////////////////////////////////////////////////////////////////////////////
//
// File: Interpolator.cpp
//
// For more information, please see: http://www.nektar.info/
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
////////////////////////////////////////////////////////////////////////////////

#include <boost/geometry.hpp>
#include <LibUtilities/BasicUtils/Interpolator.h>

using namespace std;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

namespace Nektar
{
namespace LibUtilities
{

/**
 * @brief Compute interpolation weights without doing any interpolation
 *
 * @param ptsInField    input field
 * @param ptsOutField   output field
 * @param reuseTree     if an r-tree has been constructed already, reuse it
 *                      (e.g. for repeated calls over the same input points).
 *
 * In and output fields must have the same dimension.  The most suitable
 * algorithm is chosen automatically if it wasnt set explicitly.
 */
void Interpolator::CalcWeights(const LibUtilities::PtsFieldSharedPtr ptsInField,
                               LibUtilities::PtsFieldSharedPtr &ptsOutField,
                               bool reuseTree)
{
    ASSERTL0(ptsInField->GetDim() <= m_dim, "too many dimesions in inField");
    ASSERTL0(ptsOutField->GetDim() <= m_dim, "too many dimesions in outField");

    m_ptsInField  = ptsInField;
    m_ptsOutField = ptsOutField;

    size_t nOutPts  = m_ptsOutField->GetNpoints();
    int lastProg = 0;

    // set a default method
    if (m_method == eNoMethod)
    {
        if (m_ptsInField->GetDim() == 1 || m_coordId >= 0)
        {
            m_method = eQuadratic;
        }
        else
        {
            m_method = eShepard;
        }
    }

    if ((!m_rtree || !reuseTree) && m_method != eQuadratic)
    {
        SetupTree();
    }

    switch (m_method)
    {
        case eNearestNeighbour:
        {
            m_weights   = Array<TwoD, NekDouble>(nOutPts, 1, 0.0);
            m_neighInds = Array<TwoD, unsigned int>(nOutPts, 1, (unsigned int) 0);

            for (size_t i = 0; i < nOutPts; ++i)
            {
                Array<OneD, NekDouble> tmp(m_dim, 0.0);
                for (size_t j = 0; j < m_ptsOutField->GetDim(); ++j)
                {
                    tmp[j] = m_ptsOutField->GetPointVal(j, i);
                }
                PtsPoint searchPt(i, tmp, 1E30);

                CalcW_NNeighbour(searchPt);

                int progress = int(100 * i / nOutPts);
                if (m_progressCallback && progress > lastProg)
                {
                    m_progressCallback(i, nOutPts);
                    lastProg = progress;
                }
            }

            break;
        }

        case eQuadratic:
        {
            ASSERTL0(m_ptsInField->GetDim() == 1 || m_coordId >= 0,
                     "not implemented");

            m_weights   = Array<TwoD, NekDouble>(nOutPts, 3, 0.0);
            m_neighInds = Array<TwoD, unsigned int>(nOutPts, 3, (unsigned int) 0);

            for (size_t i = 0; i < nOutPts; ++i)
            {
                Array<OneD, NekDouble> tmp(m_dim, 0.0);
                for (size_t j = 0; j < m_ptsOutField->GetDim(); ++j)
                {
                    tmp[j] = m_ptsOutField->GetPointVal(j, i);
                }
                PtsPoint searchPt(i, tmp, 1E30);

                if (m_ptsInField->GetNpoints() <= 2)
                {
                    CalcW_Linear(searchPt, m_coordId);
                }
                else
                {
                    CalcW_Quadratic(searchPt, m_coordId);
                }

                int progress = int(100 * i / nOutPts);
                if (m_progressCallback && progress > lastProg)
                {
                    m_progressCallback(i, nOutPts);
                    lastProg = progress;
                }
            }

            break;
        }

        case eShepard:
        {
            int numPts = m_ptsInField->GetDim();
            numPts     = 2 << numPts; // 2 ^ numPts
            numPts     = min(numPts, int(m_ptsInField->GetNpoints() / 2));

            m_weights   = Array<TwoD, NekDouble>(nOutPts, numPts, 0.0);
            m_neighInds = Array<TwoD, unsigned int>(nOutPts, numPts, (unsigned int) 0);

            for (size_t i = 0; i < nOutPts; ++i)
            {
                Array<OneD, NekDouble> tmp(m_dim, 0.0);
                for (size_t j = 0; j < m_ptsOutField->GetDim(); ++j)
                {
                    tmp[j] = m_ptsOutField->GetPointVal(j, i);
                }
                PtsPoint searchPt(i, tmp, 1E30);

                CalcW_Shepard(searchPt, numPts);

                int progress = int(100 * i / nOutPts);
                if (m_progressCallback && progress > lastProg)
                {
                    m_progressCallback(i, nOutPts);
                    lastProg = progress;
                }
            }

            break;
        }

        case eGauss:
        {
            ASSERTL0(m_filtWidth > NekConstants::kNekZeroTol,
                     "No filter width set");
            // use m_filtWidth as FWHM
            NekDouble sigma = m_filtWidth * 0.4246609001;

            m_maxPts = min(m_maxPts, int(m_ptsInField->GetNpoints() / 2));

            m_weights   = Array<TwoD, NekDouble>(nOutPts, m_maxPts, 0.0);
            m_neighInds = Array<TwoD, unsigned int>(nOutPts, m_maxPts, (unsigned int) 0);

            for (size_t i = 0; i < nOutPts; ++i)
            {
                Array<OneD, NekDouble> tmp(m_dim, 0.0);
                for (size_t j = 0; j < m_ptsOutField->GetDim(); ++j)
                {
                    tmp[j] = m_ptsOutField->GetPointVal(j, i);
                }
                PtsPoint searchPt(i, tmp, 1E30);

                CalcW_Gauss(searchPt, sigma, m_maxPts);

                int progress = int(100 * i / nOutPts);
                if (m_progressCallback && progress > lastProg)
                {
                    m_progressCallback(i, nOutPts);
                    lastProg = progress;
                }
            }

            break;
        }

        default:
            NEKERROR(ErrorUtil::efatal, "Invalid interpolation m_method");
            break;
    }
}

/**
 * @brief Interpolate from a pts field to a pts field
 *
 * @param ptsInField    input field
 * @param ptsOutField   output field
 *
 * In and output fields must have the same dimension and number of fields.
 * The most suitable algorithm is chosen automatically if it wasnt set
 * explicitly.
 */
void Interpolator::Interpolate(const LibUtilities::PtsFieldSharedPtr ptsInField,
                               LibUtilities::PtsFieldSharedPtr &ptsOutField)
{
    ASSERTL0(ptsInField->GetNFields() == ptsOutField->GetNFields(),
             "number of fields does not match");
    ASSERTL0(ptsInField->GetDim() <= m_dim, "too many dimesions in inField");
    ASSERTL0(ptsOutField->GetDim() <= m_dim, "too many dimesions in outField");

    m_ptsInField  = ptsInField;
    m_ptsOutField = ptsOutField;

    if (m_weights.GetRows() == 0)
    {
        CalcWeights(m_ptsInField, m_ptsOutField);
    }

    ASSERTL0(m_weights.GetRows() == m_ptsOutField->GetNpoints(),
             "weights dimension mismatch");

    size_t nFields = m_ptsOutField->GetNFields();
    size_t nOutPts = m_ptsOutField->GetNpoints();
    size_t inDim   = m_ptsInField->GetDim();

    // interpolate points and transform
    for (size_t i = 0; i < nFields; ++i)
    {
        for (size_t j = 0; j < nOutPts; ++j)
        {
            size_t nPts = m_weights.GetColumns();

            // skip if there were no neighbours found for this point
            if (nPts == 0)
            {
                continue;
            }

            NekDouble val = 0.0;
            for (size_t k = 0; k < nPts; ++k)
            {
                size_t nIdx = m_neighInds[j][k];
                val += m_weights[j][k] *
                       m_ptsInField->GetPointVal(inDim + i, nIdx);
            }
            m_ptsOutField->SetPointVal(m_ptsOutField->GetDim() + i, j, val);
        }
    }
}

int Interpolator::GetDim() const
{
    return m_dim;
}

int Interpolator::GetCoordId() const
{
    return m_coordId;
}

NekDouble Interpolator::GetFiltWidth() const
{
    return m_filtWidth;
}

InterpMethod Interpolator::GetInterpMethod() const
{
    return m_method;
}

LibUtilities::PtsFieldSharedPtr Interpolator::GetInField() const
{
    return m_ptsInField;
}

LibUtilities::PtsFieldSharedPtr Interpolator::GetOutField() const
{
    return m_ptsOutField;
}

void Interpolator::PrintStatistics()
{
    int meanN = 0;
    for (int i = 0; i < m_neighInds.GetRows(); ++i)
    {
        for (int j = 0; j < m_neighInds.GetColumns(); ++j)
        {
            if (m_neighInds[i][j] > 0)
            {
                meanN +=1;
            }
        }
    }

    cout << "Number of points: " << m_neighInds.GetRows() << endl;
    if (m_neighInds.GetRows() > 0)
    {
        cout << "mean Number of Neighbours per point: "
             << meanN / m_neighInds.GetRows() << endl;
    }
}

/**
 * @brief Computes interpolation weights using gaussian interpolation
 *
 * @param searchPt    point for which the weights are computed
 * @param sigma       standard deviation of the gauss function
 *
 * Performs an interpolation using gauss weighting. Ideal for filtering fields.
 * The filter width should be half the FWHM (= 1.1774 sigma) and must be set in
 * the constructor of the Interpolator class.
 */
void Interpolator::CalcW_Gauss(const PtsPoint &searchPt,
                               const NekDouble sigma,
                               const int maxPts)
{
    // find nearest neighbours
    vector<PtsPoint> neighbourPts;
    FindNeighbours(searchPt, neighbourPts, 4 * sigma, maxPts);
    size_t numPts = neighbourPts.size();

    // handle the cases that there was no or just one point within 4 * sigma
    if (numPts == 0)
    {
        return;
    }
    if (numPts == 1)
    {
        m_neighInds[searchPt.idx][0] = neighbourPts.front().idx;
        m_weights[searchPt.idx][0] = 1.0;

        return;
    }

    NekDouble sigmaNew = 0.25 * neighbourPts.back().dist;

    for (size_t i = 0; i < numPts; i++)
    {
        m_neighInds[searchPt.idx][i] = neighbourPts.at(i).idx;
    }

    NekDouble wSum = 0.0;
    NekDouble ts2  = 2.0 * sigmaNew * sigmaNew;
    for (size_t i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] =
            exp(-1.0 * neighbourPts[i].dist * neighbourPts[i].dist / ts2);
        wSum += m_weights[searchPt.idx][i];
    }

    for (size_t i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] = m_weights[searchPt.idx][i] / wSum;
    }
}

/**
 * @brief Computes interpolation weights using linear interpolation
 *
 * @param searchPt    point for which the weights are computed
 * @param m_coordId   coordinate id along which the interpolation should be
 * performed
 *
 * Currently, only implemented for 1D
 */
void Interpolator::CalcW_Linear(const PtsPoint &searchPt, int m_coordId)
{
    int npts = m_ptsInField->GetNpoints();
    int i;

    NekDouble coord = searchPt.coords[m_coordId];

    for (i = 0; i < npts - 1; ++i)
    {
        if ((m_ptsInField->GetPointVal(0, i) <=
             (coord + NekConstants::kNekZeroTol)) &&
            (coord <=
             (m_ptsInField->GetPointVal(0, i + 1) + NekConstants::kNekZeroTol)))
        {
            NekDouble pdiff = m_ptsInField->GetPointVal(0, i + 1) -
                              m_ptsInField->GetPointVal(0, i);

            m_neighInds[searchPt.idx][0] = i;
            m_neighInds[searchPt.idx][1] = i + 1;

            m_weights[searchPt.idx][0] =
                (m_ptsInField->GetPointVal(0, i + 1) - coord) / pdiff;
            m_weights[searchPt.idx][1] =
                (coord - m_ptsInField->GetPointVal(0, i)) / pdiff;

            break;
        }
    }
    ASSERTL0(i != npts - 1, "Failed to find coordinate " +
                                boost::lexical_cast<string>(coord) +
                                " within provided input points");
}

/**
 * @brief Computes interpolation weights using nearest neighbour interpolation
 *
 * @param searchPt    point for which the weights are computed
 * @param m_coordId   coordinate id along which the interpolation should be
 * performed
 *
 */
void Interpolator::CalcW_NNeighbour(const PtsPoint &searchPt)
{
    // find nearest neighbours
    vector<PtsPoint> neighbourPts;
    // TODO: we currently dont handle the case when there are more than one
    // most distant points (of same distance)
    FindNNeighbours(searchPt, neighbourPts, 1);

    m_neighInds[searchPt.idx][0] = neighbourPts[0].idx;
    m_weights[searchPt.idx][0] = 1.0;
}

/**
* @brief Computes interpolation weights using linear interpolation
*
* @param searchPt    point for which the weights are computed
* @param m_coordId   coordinate id along which the interpolation should be
* performed
*
* The algorithm is based on Shepard, D. (1968). A two-dimensional interpolation
* function for irregularly-spaced data. Proceedings of the 1968 23rd ACM
* National Conference. pp. 517–524.
*
* In order to save memory, for n dimesnions, only 2^n points are considered.
* Contrary to Shepard, we use a fixed number of points with fixed weighting
* factors 1/d^n.
*/
void Interpolator::CalcW_Shepard(const PtsPoint &searchPt, int numPts)
{
    // find nearest neighbours
    vector<PtsPoint> neighbourPts;
    FindNNeighbours(searchPt, neighbourPts, numPts);

    for (int i = 0; i < neighbourPts.size(); i++)
    {
        m_neighInds[searchPt.idx][i] = neighbourPts[i].idx;
    }

    // In case d < kVertexTheSameDouble ( d^2 < kNekSqrtTol), use the exact
    // point and return
    for (int i = 0; i < neighbourPts.size(); ++i)
    {
        if (neighbourPts[i].dist <= NekConstants::kNekZeroTol)
        {
            m_weights[searchPt.idx][i] = 1.0;
            return;
        }
    }

    NekDouble wSum = 0.0;
    for (int i = 0; i < neighbourPts.size(); ++i)
    {
        m_weights[searchPt.idx][i] = 1 / pow(double(neighbourPts[i].dist),
                                             double(m_ptsInField->GetDim()));
        wSum += m_weights[searchPt.idx][i];
    }

    for (int i = 0; i < neighbourPts.size(); ++i)
    {
        m_weights[searchPt.idx][i] = m_weights[searchPt.idx][i] / wSum;
    }
}

/**
 * @brief Computes interpolation weights using quadratic interpolation
 *
 * @param searchPt    point for which the weights are computed
 * @param m_coordId   coordinate id along which the interpolation should be
 * performed
 *
 * Currently, only implemented for 1D. Falls back to linear interpolation if
 * only 2
 * values are available.
 */
void Interpolator::CalcW_Quadratic(const PtsPoint &searchPt, int m_coordId)
{
    int npts = m_ptsInField->GetNpoints();
    int i;

    NekDouble coord = searchPt.coords[m_coordId];

    for (i = 0; i < npts - 1; ++i)
    {
        if ((m_ptsInField->GetPointVal(0, i) <=
             (coord + NekConstants::kNekZeroTol)) &&
            (coord <=
             (m_ptsInField->GetPointVal(0, i + 1) + NekConstants::kNekZeroTol)))
        {
            NekDouble pdiff = m_ptsInField->GetPointVal(0, i + 1) -
                              m_ptsInField->GetPointVal(0, i);
            NekDouble h1, h2, h3;

            if (i < npts - 2)
            {
                // forwards stencil
                NekDouble pdiff2 = m_ptsInField->GetPointVal(0, i + 2) -
                                   m_ptsInField->GetPointVal(0, i + 1);

                h1 = (m_ptsInField->GetPointVal(0, i + 1) - coord) *
                     (m_ptsInField->GetPointVal(0, i + 2) - coord) /
                     (pdiff * (pdiff + pdiff2));
                h2 = (coord - m_ptsInField->GetPointVal(0, i)) *
                     (m_ptsInField->GetPointVal(0, i + 2) - coord) /
                     (pdiff * pdiff2);
                h3 = (coord - m_ptsInField->GetPointVal(0, i)) *
                     (coord - m_ptsInField->GetPointVal(0, i + 1)) /
                     ((pdiff + pdiff2) * pdiff2);

                m_neighInds[searchPt.idx][0] = i;
                m_neighInds[searchPt.idx][1] = i + 1;
                m_neighInds[searchPt.idx][2] = i + 2;
            }
            else
            {
                // backwards stencil
                NekDouble pdiff2 = m_ptsInField->GetPointVal(0, i) -
                                   m_ptsInField->GetPointVal(0, i - 1);

                h1 = (m_ptsInField->GetPointVal(0, i + 1) - coord) *
                     (coord - m_ptsInField->GetPointVal(0, i - 1)) /
                     (pdiff * pdiff2);
                h2 = (coord - m_ptsInField->GetPointVal(0, i)) *
                     (coord - m_ptsInField->GetPointVal(0, i - 1)) /
                     (pdiff * (pdiff + pdiff2));
                h3 = (m_ptsInField->GetPointVal(0, i) - coord) *
                     (m_ptsInField->GetPointVal(0, i + 1) - coord) /
                     ((pdiff + pdiff2) * pdiff);

                m_neighInds[searchPt.idx][0] = i;
                m_neighInds[searchPt.idx][1] = i + 1;
                m_neighInds[searchPt.idx][2] = i - 1;
            }

            m_weights[searchPt.idx][0] = h1;
            m_weights[searchPt.idx][1] = h2;
            m_weights[searchPt.idx][2] = h3;

            break;
        }
    }
    ASSERTL0(i != npts - 1, "Failed to find coordinate " +
                                boost::lexical_cast<string>(coord) +
                                " within provided input points");
}

void Interpolator::SetupTree()
{
    std::vector<PtsPointPair> inPoints;
    for (int i = 0; i < m_ptsInField->GetNpoints(); ++i)
    {
        Array<OneD, NekDouble> coords(3, 0.0);
        for (int j = 0; j < m_ptsInField->GetDim(); ++j)
        {
            coords[j] = m_ptsInField->GetPointVal(j, i);
        }
        inPoints.push_back(
            PtsPointPair(BPoint(coords[0], coords[1], coords[2]), i));
    }
    m_rtree = MemoryManager<PtsRtree>::AllocateSharedPtr();
    m_rtree->insert(inPoints.begin(), inPoints.end());

    // remove duplicates from tree
    for (auto &it : inPoints)
    {
        std::vector<PtsPointPair> result;

        // find nearest 2 points (2 because one of these might be the one we
        // are
        // checking)
        m_rtree->query(bgi::nearest(it.first, 2),
                       std::back_inserter(result));

        // in case any of these 2 points is too close, remove the current
        // point
        // from the tree
        for (auto &it2 : result)
        {
            if (it.second != it2.second &&
                bg::distance(it.first, it2.first) <= NekConstants::kNekZeroTol)
            {
                m_rtree->remove(it);
                break;
            }
        }
    }
}

/**
 * @brief Finds the neares neighbours of a point
 *
 * @param searchPt     point for which the neighbours are searched
 * @param neighbourPts possible neighbour points
 * @param numPts       limits the number of neighbours found to the numPts
 * nearest ones
 *
 */
void Interpolator::FindNNeighbours(const PtsPoint &searchPt,
                                   vector<PtsPoint> &neighbourPts,
                                   const unsigned int numPts)
{
    std::vector<PtsPointPair> result;
    BPoint searchBPoint(searchPt.coords[0], searchPt.coords[1],
                        searchPt.coords[2]);
    m_rtree->query(bgi::nearest(searchBPoint, numPts),
                   std::back_inserter(result));

    // massage into or own format
    for (int i = 0; i < result.size(); ++i)
    {
        int idx = result[i].second;
        Array<OneD, NekDouble> coords(m_dim, 0.0);
        for (int j = 0; j < m_ptsInField->GetDim(); ++j)
        {
            coords[j] = m_ptsInField->GetPointVal(j, idx);
        }
        NekDouble d = bg::distance(searchBPoint, result[i].first);
        neighbourPts.push_back(PtsPoint(idx, coords, d));
    }

    sort(neighbourPts.begin(), neighbourPts.end());
}

/**
 * @brief Finds the neares neighbours of a point
 *
 * @param searchPt     point for which the neighbours are searched
 * @param neighbourPts possible neighbour points
 * @param dist         limits the distance of the neighbours
 *
 */
void Interpolator::FindNeighbours(const PtsPoint &searchPt,
                                  vector<PtsPoint> &neighbourPts,
                                  const NekDouble dist,
                                  const unsigned int maxPts)
{
    BPoint searchBPoint(searchPt.coords[0], searchPt.coords[1],
                        searchPt.coords[2]);
    BPoint bbMin(searchPt.coords[0] - dist, searchPt.coords[1] - dist,
                 searchPt.coords[2] - dist);
    BPoint bbMax(searchPt.coords[0] + dist, searchPt.coords[1] + dist,
                 searchPt.coords[2] + dist);
    bg::model::box<BPoint> bbox(bbMin, bbMax);

    // find points within the distance box
    std::vector<PtsPointPair> result;
    if (maxPts >= 1)
    {
        m_rtree->query(bgi::within(bbox) && bgi::nearest(searchBPoint, maxPts),
                       std::back_inserter(result));
    }
    else
    {
        m_rtree->query(bgi::within(bbox), std::back_inserter(result));
    }

    // massage into or own format
    for (int i = 0; i < result.size(); ++i)
    {
        int idx = result[i].second;
        Array<OneD, NekDouble> coords(m_dim, 0.0);
        for (int j = 0; j < m_ptsInField->GetDim(); ++j)
        {
            coords[j] = m_ptsInField->GetPointVal(j, idx);
        }
        NekDouble d = bg::distance(searchBPoint, result[i].first);

        // discard points beyonf dist
        if (d <= dist)
        {
            neighbourPts.push_back(PtsPoint(idx, coords, d));
        }
    }

    sort(neighbourPts.begin(), neighbourPts.end());
}
}
}
