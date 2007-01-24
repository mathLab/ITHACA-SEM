///////////////////////////////////////////////////////////////////////////////
//
// File GaussPoints.cpp
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
#include <LibUtilities/Foundations/GaussPoints.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void GaussPoints::CalculatePoints()
        {
            // Allocate the storage for points and weights
            PointsBaseType::CalculatePoints();
            PointsBaseType::CalculateWeights();

            int numpoints = m_pointsKey.GetNumPoints();

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::zwgj(m_points[0],m_weights,numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::zwgrjm(m_points[0],m_weights,numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::zwgrjp(m_points[0],m_weights,numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre: 
                Polylib::zwglj(m_points[0],m_weights,numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::zwgj(m_points[0],m_weights,numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::zwgrjm(m_points[0],m_weights,numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::zwgrjp(m_points[0],m_weights,numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev: 
                Polylib::zwglj(m_points[0],m_weights,numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::zwgrjm(m_points[0],m_weights,numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::zwgrjm(m_points[0],m_weights,numpoints,0.0,2.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        void GaussPoints::CalculateWeights()
        {
            //For Gauss Quadrature, this is done as part of the points computation
        }

        void GaussPoints::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix.
            Points<double>::CalculateDerivMatrix();
        }

        boost::shared_ptr< PointsBaseType > GaussPoints::Create(const PointsKey &key)
        {
            boost::shared_ptr< PointsBaseType > returnval(new GaussPoints(key));

            returnval->Initialize();

            return returnval;
        }

    } // end of namespace LibUtilities
} // end of namespace Nektar




