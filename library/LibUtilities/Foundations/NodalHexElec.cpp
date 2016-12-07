///////////////////////////////////////////////////////////////////////////////
//
// File NodalHexElec.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: NodalHexElec
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/NodalHexElec.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/Points.h>

namespace Nektar
{
namespace LibUtilities
{
void NodalHexElec::CalculatePoints()
{
    // Allocate the storage for points
    unsigned int numPoints = GetNumPoints();

    PointsKey e(numPoints, eGaussLobattoLegendre);
    PointsManager()[e]->GetPoints(m_e0);
    m_ew = PointsManager()[e]->GetW();

    for (int i = 0; i < 3; i++)
    {
        m_points[i] = Array<OneD, DataType>(numPoints*numPoints*numPoints);
    }

    for (int k = 0, ct = 0; k < numPoints; k++)
    {
        for (int j = 0; j < numPoints; j++)
        {
            for (int i = 0; i < numPoints; i++, ct++)
            {
                m_points[0][ct] = m_e0[i];
                m_points[1][ct] = m_e0[j];
                m_points[2][ct] = m_e0[k];
            }
        }
    }
}

void NodalHexElec::CalculateWeights()
{
    unsigned int numPoints = GetNumPoints();

    m_weights = Array<OneD, DataType>(numPoints*numPoints*numPoints);

    for (int k = 0, ct = 0; k < numPoints; k++)
    {
        for (int j = 0; j < numPoints; j++)
        {
            for (int i = 0; i < numPoints; i++, ct++)
            {
                m_weights[ct] = m_ew[i] * m_ew[j] * m_ew[k];
            }
        }
    }
}

void NodalHexElec::CalculateDerivMatrix()
{
}

boost::shared_ptr<PointsBaseType> NodalHexElec::Create(const PointsKey &key)
{
    boost::shared_ptr<PointsBaseType> returnval(
        MemoryManager<NodalHexElec>::AllocateSharedPtr(key));
    returnval->Initialize();
    return returnval;
}

}
}
