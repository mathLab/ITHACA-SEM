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
// Description: GaussPoints Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

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
                Polylib::zwgj(m_points[0].get(),m_weights.get(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::zwgrjp(m_points[0].get(),m_weights.get(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre: 
                Polylib::zwglj(m_points[0].get(),m_weights.get(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::zwgj(m_points[0].get(),m_weights.get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::zwgrjp(m_points[0].get(),m_weights.get(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev: 
                Polylib::zwglj(m_points[0].get(),m_weights.get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::zwgrjm(m_points[0].get(),m_weights.get(),numpoints,2.0,0.0);
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
            Points<NekDouble>::CalculateDerivMatrix();

            int numpoints = m_pointsKey.GetNumPoints();
            int totpoints = m_pointsKey.GetTotNumPoints();
            double * dmtemp = new NekDouble[totpoints*totpoints];

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Dgj(dmtemp,m_points[0].get(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Dgrjp(dmtemp,m_points[0].get(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Dglj(dmtemp,m_points[0].get(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Dgj(dmtemp,m_points[0].get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Dgrjp(dmtemp,m_points[0].get(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Dglj(dmtemp,m_points[0].get(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0].get(),numpoints,2.0,0.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }

            std::copy(dmtemp,dmtemp+totpoints*totpoints,m_derivmatrix->begin());

            delete[] dmtemp;
        }

        void GaussPoints::CalculateInterpMatrix(unsigned int npts, SharedArray<const NekDouble> xpoints, SharedArray<NekDouble> interp)
        {

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Imgj(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Imgrjp(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Imglj(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Imgj(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Imgrjp(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Imglj(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,0.0,2.0);
		        break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::Imgrjm(interp.get(),m_points[0].get(),xpoints.get(),GetNumPoints(),npts,2.0,0.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        boost::shared_ptr< PointsBaseType > GaussPoints::Create(const PointsKey &pkey)
        {
            boost::shared_ptr< PointsBaseType > returnval(new GaussPoints(pkey));

            returnval->Initialize();

            return returnval;
        }

        boost::shared_ptr< NekMatrix<double> > GaussPoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            SharedArray<const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(const PointsKey &pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(SharedArray<const NekDouble> x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(unsigned int numpoints, SharedArray<const NekDouble> x)
        {
            SharedArray<NekDouble> interp = MemoryManager::AllocateSharedArray<double>(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            boost::shared_ptr< NekMatrix<DataType> > returnval(new NekMatrix<DataType>(numpoints,GetNumPoints(),interp.get()));

            return returnval;
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar

/**
* %Log%
*/

