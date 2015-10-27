////////////////////////////////////////////////////////////////////////////////
//
// File: Interpolator.cpp
//
// For more information, please see: http://www.nektar.info/
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
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/Interpolator.h>

namespace Nektar
{
namespace LibUtilities
{

Interpolator::Interpolator(const PtsFieldSharedPtr inField, PtsFieldSharedPtr &outField) :
            m_inField(inField),
            m_outField(outField),
            m_method(ePtsNoMethod),
            m_dim(3)
{
    ASSERTL0(inField->GetNFields() == outField->GetNFields(), "number of fields does not match");
}



/**
 * @brief Compute the weights for an interpolation of the field values to physical points
 *
 * @param physCoords    coordinates of the physical points
 * @param coord_id      id of the coordinate to use for interpolation.
 *
 * Set coord_id to -1 to use n-D interpolation for an n-dimensional field.
 * The most suitable algorithm is chosen automatically.
 */
void Interpolator::CalcWeights(PtsInterpMethod method, short int coordId, NekDouble width)
{
    int nOutPts = m_outField->GetNpoints();
    int lastProg = 0;

    NekDouble sigma = width * 0.4246609001;;
    if (method == ePtsGauss)
    {
        ASSERTL0(width > NekConstants::kNekZeroTol, "No filter width set");
    }

    m_weights = Array<OneD, Array<OneD, float> >(nOutPts);
    m_neighInds = Array<OneD, Array<OneD, unsigned int> >(nOutPts);

    std::vector<PtsPoint > inPoints;
    for (int i = 0; i < m_inField->GetNpoints(); ++i)
    {
        Array<OneD, NekDouble> coords(m_dim);
        for (int j = 0; j < m_dim; ++j)
        {
            coords[j] = m_inField->GetPointVal(j,i);
        }
        inPoints.push_back(PtsPoint(i, coords, 1E30));
    }
    m_rtree.insert(inPoints.begin(), inPoints.end());

    // interpolate points and transform
    for (int i = 0; i < nOutPts; ++i)
    {
        Array<OneD, NekDouble> tmp(m_dim);
        for (int j = 0; j < m_dim; ++j)
        {
            tmp[j] = m_outField->GetPointVal(j,i);
        }
        PtsPoint searchPt(i, tmp, 1E30);

        if ((m_dim == 1 || coordId >= 0) && ! (method == ePtsGauss || method == ePtsShepard))
        {
            if (m_dim == 1)
            {
                coordId = 0;
            }

            if (m_outField->GetNpoints() <= 2)
            {
                CalcW_Linear(searchPt, coordId);
            }
            else
            {
                CalcW_Quadratic(searchPt, coordId);
            }
        }
        else
        {
            ASSERTL0(method == ePtsGauss || method == ePtsShepard, "interpolation method not implemented for this dimension");
            switch ( method )
            {
                case ePtsGauss:
                    CalcW_Gauss(searchPt, sigma);
                    break;
                case ePtsShepard:
                    CalcW_Shepard(searchPt);
                    break;
                default:
                    ASSERTL0(false, "unknown interpolation method");
            }
        }

        int progress = int(100 * i / nOutPts);
        if (m_progressCallback && progress > lastProg)
        {
            m_progressCallback(i, nOutPts);
            lastProg = progress;
        }
    }
}


/**
 * @brief Perform the interpolate of field values to physical points
 *
 * @param intFields     interpolated field at the physical points
 *
 * The weights must have already been computed by @CalcWeights or set by
 * @SetWeights.
 */
void Interpolator::Interpolate(
            Nektar::LibUtilities::PtsInterpMethod method, short int coordId, Nektar::NekDouble width)
{
    if ( m_weights.num_elements() == 0)
    {
        CalcWeights(method, coordId, width);
    }

    ASSERTL1(m_weights[0].num_elements() == m_neighInds[0].num_elements(),
             "weights / neighInds mismatch")
    int nFields = m_outField->GetNFields();
    int nOutPts = m_outField->GetNpoints();
    int inDim = m_inField->GetDim();

    // interpolate points and transform
    for (int i = 0; i < nFields; ++i)
    {
        for (int j = 0; j < nOutPts; ++j)
        {
            NekDouble val = 0.0;
            int nPts = m_weights[j].num_elements();
            for (int k = 0; k < nPts; ++k)
            {
                unsigned int nIdx = m_neighInds[j][k];
                val += m_weights[j][k] * m_inField->GetPointVal(inDim + i, nIdx);
            }
            m_outField->SetPointVal(i,j, val);
        }
    }
}


int Interpolator::GetDim() const
{
    return m_dim;
}


PtsFieldSharedPtr Interpolator::GetInField() const
{
    return m_inField;
}


PtsFieldSharedPtr Interpolator::GetOutField() const
{
    return m_outField;
}


/**
 * @brief Get the interpolation weights and corresponding neighbour indices
 *
 * @param weights       Interpolation weights for each neighbour.
 * Structure: m_weights[physPtIdx][neighbourIdx]
 * @param neighbourInds Indices of the relevant neighbours for each physical point.
 * Structure: m_neighInds[ptIdx][neighbourIdx]
 */
void Interpolator::GetWeights(
            Array<OneD, Array<OneD, float> > &weights,
            Array<OneD, Array<OneD, unsigned int> > &neighbourInds) const
{
    weights = m_weights;
    neighbourInds = m_neighInds;
}


/**
 * @brief Set the interpolation weights for an interpolation
 *
 * @param weights       Interpolation weights for each neighbour.
 * Structure: m_weights[physPtIdx][neighbourIdx]
 * @param neighbourInds Indices of the relevant neighbours for each physical point.
 * Structure: m_neighInds[ptIdx][neighbourIdx]
 */
void Interpolator::SetWeights(const Array<OneD, Array<OneD, float> >
         &weights,
         const Array<OneD, Array<OneD, unsigned int> > &neighbourInds)
{
    ASSERTL0(weights.num_elements() ==  neighbourInds.num_elements(),
             "weights and neighbourInds not of same number of physical points")

    m_weights = weights;
    m_neighInds = neighbourInds;

}


/**
 * @brief Compute interpolation weights for a physical point using Gauss filtering
 *
 * @param physPtIdx         The index of the physical point in its storage array
 * @param physPt            The coordinates of the physical point
 *
 */
void Interpolator::CalcW_Gauss(const PtsPoint &searchPt, const NekDouble sigma)
{
    NekDouble ts2 = 2 * sigma * sigma;
    NekDouble fac = 1.0 / (sigma * sqrt(2 * M_PI));
    fac = pow(fac, m_dim);

    // find nearest neighbours
    int maxPts = 500;
    vector<PtsPoint > neighbourPts;
    FindNeighbours(searchPt, neighbourPts, 1.96 * sigma);
    int numPts = min( (int) neighbourPts.size(), maxPts);

    // handle the case that there was no point within 1.96 * sigma
    if (numPts == 0)
    {
        m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (0);
        m_weights[searchPt.idx] = Array<OneD, float> (0);

        return;
    }

    m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (numPts);
    for (int i = 0; i < numPts; i++)
    {
        m_neighInds[searchPt.idx][i] = neighbourPts.at(i).idx;
    }

    m_weights[searchPt.idx] = Array<OneD, float> (numPts, 0.0);

    NekDouble wSum = 0.0;
    for (int i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] = fac * exp(-1 * pow(neighbourPts[i].dist, 2.0f) / ts2);
        wSum += m_weights[searchPt.idx][i];
    }

    for (int i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] = m_weights[searchPt.idx][i] / wSum;
    }

    ASSERTL0(Vmath::Nnan(numPts, m_weights[searchPt.idx], 1) == 0, "NaN found in weights");
}


/**
 * @brief Compute interpolation weights for a 1D physical point using linear
 * interpolation.
 *
 * @param physPtIdx         The index of the physical point in its storage array
 * @param coord             The coordinate of the physical point
 */

void Interpolator::CalcW_Linear(const PtsPoint &searchPt, int coordId)
{
    int npts = m_inField->GetNpoints();
    int i;

    NekDouble coord = searchPt.coords[coordId];

    int numPts = 2;
    m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (numPts);
    m_weights[searchPt.idx] = Array<OneD, float> (numPts, 0.0);

    for (i = 0; i < npts - 1; ++i)
    {
        if ((m_inField->GetPointVal(0,i) <= coord) && (coord <= m_inField->GetPointVal(0,i + 1)))
        {
            NekDouble pdiff = m_inField->GetPointVal(0,i+1) - m_inField->GetPointVal(0,i);

            m_neighInds[searchPt.idx][0] = i;
            m_neighInds[searchPt.idx][1] = i + 1;

            m_weights[searchPt.idx][0] = (m_inField->GetPointVal(0,i+1) - coord) / pdiff;
            m_weights[searchPt.idx][1] = (coord - m_inField->GetPointVal(0,i)) / pdiff;

            break;
        }
    }
    ASSERTL0(i != npts - 1, "Failed to find coordinate " +
             boost::lexical_cast<string>(coord) +
             " within provided input points");
};


void Interpolator::CalcW_NNeighbour(const PtsPoint &searchPt, int coordId)
{
    // find nearest neighbours
    vector<PtsPoint > neighbourPts;
    FindNNeighbours(searchPt, neighbourPts, 1);

    m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (1, neighbourPts.at(0).idx);
    m_weights[searchPt.idx] = Array<OneD, float> (1, 1.0);

    ASSERTL0(Vmath::Nnan(1, m_weights[searchPt.idx], 1) == 0, "NaN found in weights");
}


/**
 * @brief Compute interpolation weights for a physical point using a modified
 * Shepard algorithm.
 *
 * @param physPtIdx         The index of the physical point in its storage array
 * @param physPt            The coordinates of the physical point
 *
 * The algorithm is based on Shepard, D. (1968). A two-dimensional interpolation
 * function for irregularly-spaced data. Proceedings of the 1968 23rd ACM
 * National Conference. pp. 517–524.
 *
 * In order to save memory, for n dimesnions, only 2^n points are considered.
 * Contrary to Shepard, we use a fixed number of points with fixed weighting
 * factors 1/d^n.
 */

void Interpolator::CalcW_Shepard(const PtsPoint &searchPt)
{
    // find nearest neighbours
    vector<PtsPoint > neighbourPts;
    int numPts = pow(float(2), m_dim);
    numPts = min(numPts, int(m_inField->GetNpoints() / 2));
    FindNNeighbours(searchPt, neighbourPts, numPts);

    m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (numPts);
    for (int i = 0; i < numPts; i++)
    {
        m_neighInds[searchPt.idx][i] = neighbourPts.at(i).idx;
    }

    m_weights[searchPt.idx] = Array<OneD, float> (numPts, 0.0);

    // In case d < kVertexTheSameDouble ( d^2 < kNekSqrtTol), use the exact
    // point and return
    for (int i = 0; i < numPts; ++i)
    {
        if (neighbourPts[i].dist <= NekConstants::kNekZeroTol)
        {
            m_weights[searchPt.idx][i] = 1.0;
            return;
        }
    }

    NekDouble wSum = 0.0;
    for (int i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] = 1 / pow(double(neighbourPts[i].dist),
                                          double(3)); //TODO: 3 or m_dim?
        wSum += m_weights[searchPt.idx][i];
    }

    for (int i = 0; i < numPts; ++i)
    {
        m_weights[searchPt.idx][i] = m_weights[searchPt.idx][i] / wSum;
    }

    ASSERTL0(Vmath::Nnan(numPts, m_weights[searchPt.idx], 1) == 0, "NaN found in weights");

}

/**
* @brief Compute interpolation weights for a 1D physical point using quadratic
* interpolation.
*
* @param physPtIdx         The index of the physical point in its storage array
* @param coord             The coordinate of the physical point
*/

void Interpolator::CalcW_Quadratic(const PtsPoint &searchPt, int coordId)
{
    int npts = m_inField->GetNpoints();
    int i;

    NekDouble coord = searchPt.coords[coordId];

    int numPts = 3;
    m_neighInds[searchPt.idx] = Array<OneD, unsigned int> (numPts);
    m_weights[searchPt.idx] = Array<OneD, float> (numPts, 0.0);

    for (i = 0; i < npts - 1; ++i)
    {
        if ((m_inField->GetPointVal(0,i) <= coord) && (coord <= m_inField->GetPointVal(0,i+1)))
        {
            NekDouble pdiff = m_inField->GetPointVal(0,i+1) - m_inField->GetPointVal(0,i);
            NekDouble h1, h2, h3;

            if (i < npts - 2)
            {
                // forwards stencil
                NekDouble pdiff2 = m_inField->GetPointVal(0,i+2) - m_inField->GetPointVal(0,i+1);

                h1 = (m_inField->GetPointVal(0,i+1) - coord)
                     * (m_inField->GetPointVal(0,i+2) - coord)
                     / (pdiff * (pdiff + pdiff2));
                h2 = (coord - m_inField->GetPointVal(0,i))
                     * (m_inField->GetPointVal(0,i+2) - coord)
                     / (pdiff * pdiff2);
                h3 = (coord - m_inField->GetPointVal(0,i))
                     * (coord - m_inField->GetPointVal(0,i+1))
                     / ((pdiff + pdiff2) * pdiff2);

                m_neighInds[searchPt.idx][0] = i;
                m_neighInds[searchPt.idx][1] = i + 1;
                m_neighInds[searchPt.idx][2] = i + 2;
            }
            else
            {
                // backwards stencil
                NekDouble pdiff2 = m_inField->GetPointVal(0,i) - m_inField->GetPointVal(0,i-1);

                h1 = (m_inField->GetPointVal(0,i+1) - coord)
                     * (coord - m_inField->GetPointVal(0,i-1))
                     / (pdiff * pdiff2);
                h2 = (coord - m_inField->GetPointVal(0,i))
                     * (coord - m_inField->GetPointVal(0,i-1))
                     / (pdiff * (pdiff + pdiff2));
                h3 = (m_inField->GetPointVal(0,i) - coord)
                     * (m_inField->GetPointVal(0,i+1) - coord)
                     / ((pdiff + pdiff2) * pdiff);

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
};


/**
 * @brief Find nearest neighbours using a brute-force "algorithm".
 *
 * @param physPt              Coordinates of the physical point its neighbours
 * we are looking for
 * @param neighbourPts        The points we found
 * @param numPts              The number of points to find
 *
 * This iterates over all points, computes the (squared) euclidean distance
 * and chooses the numPts closest points. Thus, its very expensive and
 * inefficient.
 */

void Interpolator::FindNNeighbours(const PtsPoint &searchPt,
                               vector< PtsPoint > &neighbourPts,
                               const unsigned int numPts)
{
    // find points within the distance box
    m_rtree.query(bgi::nearest(searchPt, numPts), std::back_inserter(neighbourPts));

    for(typename vector<PtsPoint >::iterator it = neighbourPts.begin(); it != neighbourPts.end(); ++it)
    {
        it->dist = bg::distance(searchPt, *it);
    }

    sort(neighbourPts.begin(), neighbourPts.end());
}


/**
 * @brief Find nearest neighbours using a brute-force "algorithm".
 *
 * @param physPt              Coordinates of the physical point its neighbours
 * we are looking for
 * @param neighbourPts        The points we found
 * @param dist                The max distance of a neighbour point
 *
 * This iterates over all points, computes the (squared) euclidean distance
 * and chooses the points within the defined distance. Thus, its very expensive
 * and inefficient.
 */

void Interpolator::FindNeighbours(const PtsPoint &searchPt,
                              vector<PtsPoint > &neighbourPts,
                              const NekDouble dist)
{
    PtsPoint bbMin, bbMax;
    bg::strategy::transform::translate_transformer<PtsPoint, PtsPoint > t1(- dist, - dist, - dist);
    bg::strategy::transform::translate_transformer<PtsPoint, PtsPoint > t2(dist, dist, dist);
    bg::transform(searchPt, bbMin, t1);
    bg::transform(searchPt, bbMax, t2);
    PtsBox bbox(bbMin, bbMax);

    // find points within the distance box
    m_rtree.query(bgi::within(bbox), std::back_inserter(neighbourPts));

    for(typename vector<PtsPoint >::iterator it = neighbourPts.begin(); it != neighbourPts.end(); ++it)
    {
        it->dist = bg::distance(searchPt, *it);
    }

    sort(neighbourPts.begin(), neighbourPts.end());

    // remove everything beyond the distance
    for(typename vector<PtsPoint >::iterator it = neighbourPts.begin(); it != neighbourPts.end(); ++it)
    {
        if (it->dist > dist)
        {
            neighbourPts.erase(it, neighbourPts.end());
            break;
        }
    }
}


}
}
