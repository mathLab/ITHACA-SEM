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
        namespace
        {
            const bool gaussInited1 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
            const bool gaussInited2 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMLegendre), GaussPoints::Create);
            const bool gaussInited3 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPLegendre), GaussPoints::Create);
            const bool gaussInited4 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoLegendre), GaussPoints::Create);
            const bool gaussInited5 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussChebyshev), GaussPoints::Create);
            const bool gaussInited6 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMChebyshev), GaussPoints::Create);
            const bool gaussInited7 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPChebyshev), GaussPoints::Create);
            const bool gaussInited8 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1), GaussPoints::Create);
            const bool gaussInited9 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
            const bool gaussInited10 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2), GaussPoints::Create);
        }

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

            int numpoints = m_pointsKey.GetNumPoints();
            int totpoints = m_pointsKey.GetTotNumPoints();
            double * dmtemp = new double[totpoints*totpoints];

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Dgj(dmtemp,m_points[0],numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Dgrjm(dmtemp,m_points[0],numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Dgrjp(dmtemp,m_points[0],numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Dglj(dmtemp,m_points[0],numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Dgj(dmtemp,m_points[0],numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Dgrjm(dmtemp,m_points[0],numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Dgrjp(dmtemp,m_points[0],numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Dglj(dmtemp,m_points[0],numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Dgrjm(dmtemp,m_points[0],numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Dgrjm(dmtemp,m_points[0],numpoints,0.0,2.0);
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }

            std::copy(dmtemp,dmtemp+totpoints*totpoints,m_derivmatrix->begin());

            delete[] dmtemp;
        }

        void GaussPoints::CalculateInterpMatrix(unsigned int npts, const double * xpoints, double * interp)
        {

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Imgj(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Imgrjm(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Imgrjp(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Imglj(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Imgj(interp,m_points[0],xpoints,GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Imgrjm(interp,m_points[0],xpoints,GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Imgrjp(interp,m_points[0],xpoints,GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Imglj(interp,m_points[0],xpoints,GetNumPoints(),npts,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Imgrjm(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Imgrjm(interp,m_points[0],xpoints,GetNumPoints(),npts,0.0,2.0);
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

        boost::shared_ptr< NekMatrix<double> > GaussPoints::CreateMatrix(const PointsKey &ownpkey, const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            const double * xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return PointsManager()[ownpkey]->GetI(numpoints, xpoints);
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(const PointsKey &pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(const double* x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<double> > GaussPoints::GetI(unsigned int numpoints, const double *x)
        {
            double * interp = new double[GetNumPoints()*numpoints];

            CalculateInterpMatrix(numpoints, x, interp);

            boost::shared_ptr< NekMatrix<DataType> > returnval(new NekMatrix<DataType>(numpoints,GetNumPoints(),interp));

            delete[] interp;

            return returnval;
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar

/**
* %Log%
*/

