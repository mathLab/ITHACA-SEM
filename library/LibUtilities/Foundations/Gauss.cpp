///////////////////////////////////////////////////////////////////////////////
//
// File Points1D.cpp
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
// Description: 1D Points definitions 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/Gauss.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void GaussPolyPoints::CalculatePoints()
        {
            // Allocate the points.
            Points<double>::CalculatePoints();

            // Now set the points.
        }

        void GaussPolyPoints::CalculateWeights()
        {
            // Allocate the weights.
            Points<double>::CalculateWeights();
        }

        void GaussPolyPoints::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix.
            Points<double>::CalculateDerivMatrix();
        }

        GaussPolyPoints::GaussPolyPoints(const PointsKey &key, 
            const double alpha, const double beta):
        Points<double>(key),
            m_alpha(alpha),
            m_beta(beta)
        {
        }

        boost::shared_ptr< Points<double> > GaussPolyPoints::Create(const PointsKey &key)
        {
            boost::shared_ptr< Points<double> > returnval(new GaussPolyPoints(key, 1.0, 1.0));

            returnval->Initialize();

            return returnval;
        }

    } // end of namespace stdregion
} // end of namespace stdregion




