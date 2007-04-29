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
        template<typename T>
        void GaussPoints<T>::CalculatePoints()
        {
            // Allocate the storage for points and weights
            PointsBaseType::CalculatePoints();
            PointsBaseType::CalculateWeights();

            int numpoints = m_pointsKey.GetNumPoints();

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::zwgj(m_points[0]->data(),m_weights->data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::zwgrjp(m_points[0]->data(),m_weights->data(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre: 
                Polylib::zwglj(m_points[0]->data(),m_weights->data(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::zwgj(m_points[0]->data(),m_weights->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::zwgrjp(m_points[0]->data(),m_weights->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev: 
                Polylib::zwglj(m_points[0]->data(),m_weights->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::zwgrjm(m_points[0]->data(),m_weights->data(),numpoints,2.0,0.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        template<typename T>
        void GaussPoints<T>::CalculateWeights()
        {
            //For Gauss Quadrature, this is done as part of the points computation
        }

        template<typename T>
        void GaussPoints<T>::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix.
            Points<T>::CalculateDerivMatrix();

            int numpoints = m_pointsKey.GetNumPoints();
            int totpoints = m_pointsKey.GetTotNumPoints();
            T * dmtemp = new T[totpoints*totpoints];

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Dgj(dmtemp,m_points[0]->data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Dgrjp(dmtemp,m_points[0]->data(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Dglj(dmtemp,m_points[0]->data(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Dgj(dmtemp,m_points[0]->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Dgrjp(dmtemp,m_points[0]->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Dglj(dmtemp,m_points[0]->data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0]->data(),numpoints,2.0,0.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }

            std::copy(dmtemp,dmtemp+totpoints*totpoints,m_derivmatrix->begin());

            delete[] dmtemp;
        }

        template<typename T>
        void GaussPoints<T>::CalculateInterpMatrix(unsigned int npts, typename Nek1DConstSharedArray<T>::type xpoints, typename Nek1DSharedArray<T>::type interp)
        {
            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Imgj(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Imgrjp(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Imglj(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Imgj(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Imgrjp(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Imglj(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,0.0,2.0);
		        break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::Imgrjm(interp.get(),m_points[0]->data(),xpoints.get(),GetNumPoints(),npts,2.0,0.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        template<typename T>
        boost::shared_ptr<typename GaussPoints<T>::PointsBaseType> GaussPoints<T>::Create(const PointsKey &pkey)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager::AllocateSharedPtr<GaussPoints<T> >(pkey));

            returnval->Initialize();

            return returnval;
        }

        template<typename T>
        boost::shared_ptr< NekMatrix<T> > GaussPoints<T>::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            typename Nek1DConstSharedArray<T>::type xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > GaussPoints<T>::GetI(const PointsKey &pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > GaussPoints<T>::GetI(typename Nek1DConstSharedArray<T>::type x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > GaussPoints<T>::GetI(unsigned int numpoints, typename Nek1DConstSharedArray<T>::type x)
        {
            typename Nek1DSharedArray<T>::type interp = MemoryManager::AllocateSharedPtr<typename Nek1DArray<T>::type>(boost::extents[GetNumPoints()*numpoints]);

            CalculateInterpMatrix(numpoints, x, interp);

            boost::shared_ptr< NekMatrix<T> > returnval(MemoryManager::AllocateSharedPtr<NekMatrix<T> >(numpoints,GetNumPoints(),interp->data()));

            return returnval;
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar

/**
* $Log$
*/

