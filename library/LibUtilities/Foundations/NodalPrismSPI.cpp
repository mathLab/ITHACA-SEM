///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriElec.cpp
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
// Description: 2D Nodal Triangle Fekete Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/NodalPrismSPI.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Foundations/NodalUtil.h>

namespace Nektar
{
    namespace LibUtilities
    {
        void NodalPrismSPI::CalculatePoints()
        {
            // Allocate the storage for points
            unsigned int numPoints = GetNumPoints();

            PointsKey e(numPoints,eGaussLobattoLegendre);
            PointsManager()[e]->GetPoints(e0);
            ew = PointsManager()[e]->GetW();

            PointsKey t(numPoints,eNodalTriSPI);
            PointsManager()[t]->GetPoints(t0,t1);
            tw = PointsManager()[t]->GetW();
            numtri = tw.num_elements();

            for(int i = 0; i < 3; i++)
            {
                m_points[i] = Array<OneD, DataType>(numtri*numPoints);
            }

            for(int j = 0, ct = 0; j < numPoints; j++)
            {
                for(int i = 0; i < numtri; i++, ct++)
                {
                    //need to flip y and z because of quad orientation
                    m_points[0][ct] = t0[i];
                    m_points[1][ct] = e0[j];
                    m_points[2][ct] = t1[i];
                }
            }
// exit(0);
        }

        void NodalPrismSPI::CalculateWeights()
        {
            unsigned int numPoints = GetNumPoints();

            m_weights = Array<OneD, DataType>(numtri*numPoints);

            for(int j = 0, ct = 0; j < numPoints; j++)
            {
                for(int i = 0; i < numtri; i++, ct++)
                {
                    m_weights[ct] = tw[i] * ew[j];
                }
            }
        }

        void NodalPrismSPI::CalculateDerivMatrix()
        {

        }

        boost::shared_ptr<PointsBaseType> NodalPrismSPI::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager<NodalPrismSPI>::AllocateSharedPtr(key));
            returnval->Initialize();
            return returnval;
        }
    } // end of namespace stdregion
} // end of namespace stdregion
