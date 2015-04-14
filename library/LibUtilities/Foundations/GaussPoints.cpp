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

#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>


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
                Polylib::zwgj(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::zwgrjp(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre: 
                Polylib::zwglj(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::zwgj(m_points[0].data(),m_weights.data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::zwgrjp(m_points[0].data(),m_weights.data(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev: 
                Polylib::zwglj(m_points[0].data(),m_weights.data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::zwgrjm(m_points[0].data(),m_weights.data(),numpoints,2.0,0.0);
                break;

            case eGaussKronrodLegendre:
                Polylib::zwgk(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;
	      
            case eGaussRadauKronrodMLegendre:
                Polylib::zwrk(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;
	      
            case eGaussRadauKronrodMAlpha1Beta0:
                Polylib::zwrk(m_points[0].data(),m_weights.data(),numpoints,1.0,0.0);
                break;
	      
            case eGaussLobattoKronrodLegendre:
                Polylib::zwlk(m_points[0].data(),m_weights.data(),numpoints,0.0,0.0);
                break;
		
            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        void GaussPoints::CalculateWeights()
        {
            // For Gauss Quadrature, this is done as part of the points computation
        }

        void GaussPoints::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix
            PointsBaseType::CalculateDerivMatrix();

            int numpoints = m_pointsKey.GetNumPoints();
            int totpoints = m_pointsKey.GetTotNumPoints();
            Array<OneD, NekDouble> dmtempSharedArray = Array<OneD, NekDouble>(totpoints*totpoints);
            NekDouble *dmtemp = dmtempSharedArray.data();

            switch(m_pointsKey.GetPointsType())
            {
            case eGaussGaussLegendre:
                Polylib::Dgj(dmtemp,m_points[0].data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauMLegendre:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,0.0,0.0);
                break;

            case eGaussRadauPLegendre:
                Polylib::Dgrjp(dmtemp,m_points[0].data(),numpoints,0.0,0.0);
                break;

            case eGaussLobattoLegendre:
                Polylib::Dglj(dmtemp,m_points[0].data(),numpoints,0.0,0.0);
                break;

            case eGaussGaussChebyshev:
                Polylib::Dgj(dmtemp,m_points[0].data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMChebyshev:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauPChebyshev:
                Polylib::Dgrjp(dmtemp,m_points[0].data(),numpoints,-0.5,-0.5);
                break;

            case eGaussLobattoChebyshev:
                Polylib::Dglj(dmtemp,m_points[0].data(),numpoints,-0.5,-0.5);
                break;

            case eGaussRadauMAlpha0Beta1:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,0.0,1.0);
                break;

            case eGaussRadauMAlpha0Beta2:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,0.0,2.0);
                break;

            case eGaussRadauMAlpha1Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,1.0,0.0);
                break;

            case eGaussRadauMAlpha2Beta0:
                Polylib::Dgrjm(dmtemp,m_points[0].data(),numpoints,2.0,0.0);
                break;
		
            case eGaussKronrodLegendre:
            case eGaussRadauKronrodMLegendre:
            case eGaussRadauKronrodMAlpha1Beta0:
            case eGaussLobattoKronrodLegendre:
                {
                for(unsigned int i=0;i<m_pointsKey.GetNumPoints();++i)
                {
                    for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                    {
                        (*m_derivmatrix[0])(i,j) = LagrangePolyDeriv(m_points[0][i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                    }
                }
                return;
                }
                break;
	      
            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
	    
            std::copy(dmtemp,dmtemp+totpoints*totpoints,m_derivmatrix[0]->begin());
        }

        void GaussPoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {
            switch(m_pointsKey.GetPointsType())
            {
                case eGaussGaussLegendre:
                    Polylib::Imgj(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,0.0);
                    break;

                case eGaussRadauMLegendre:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,0.0);
                    break;

                case eGaussRadauPLegendre:
                    Polylib::Imgrjp(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,0.0);
                    break;

                case eGaussLobattoLegendre:
                    Polylib::Imglj(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,0.0);
                    break;

                case eGaussGaussChebyshev:
                    Polylib::Imgj(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,-0.5,-0.5);
                    break;

                case eGaussRadauMChebyshev:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,-0.5,-0.5);
                    break;

                case eGaussRadauPChebyshev:
                    Polylib::Imgrjp(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,-0.5,-0.5);
                    break;

                case eGaussLobattoChebyshev:
                    Polylib::Imglj(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,-0.5,-0.5);
                    break;

                case eGaussRadauMAlpha0Beta1:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,1.0);
                    break;

                case eGaussRadauMAlpha0Beta2:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,0.0,2.0);
                    break;

                case eGaussRadauMAlpha1Beta0:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,1.0,0.0);
                    break;

                case eGaussRadauMAlpha2Beta0:
                    Polylib::Imgrjm(interp.data(),m_points[0].data(),xpoints.data(),GetNumPoints(),npts,2.0,0.0);
                    break;
		
                case eGaussKronrodLegendre:
                case eGaussRadauKronrodMLegendre:
                case eGaussRadauKronrodMAlpha1Beta0:
                case eGaussLobattoKronrodLegendre:
                {
                    for(unsigned int i=0;i<npts;++i)
                    {
                        for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                        {
                            interp[i + j*npts] = LagrangePoly(xpoints[i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                        }
                    }
                }
                break;

            default:
                ASSERTL0(false, "Unknown Gauss quadrature point distribution requested");
            }
        }

        boost::shared_ptr<Points<NekDouble> > GaussPoints::Create(const PointsKey &pkey)
        {
            boost::shared_ptr< Points<NekDouble> > returnval(MemoryManager< GaussPoints >::AllocateSharedPtr(pkey));

            returnval->Initialize();

            return returnval;
        }

        boost::shared_ptr< NekMatrix<NekDouble> > GaussPoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            // Delegate to function below
            return GetI(numpoints, xpoints);
        }


        const boost::shared_ptr<NekMatrix<NekDouble> > GaussPoints::GetI(const PointsKey &pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const boost::shared_ptr<NekMatrix<NekDouble> > GaussPoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            // Delegate to function below
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<NekDouble> > GaussPoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            NekDouble* t = interp.data();
            unsigned int np = GetNumPoints();
            boost::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np,t));

            return returnval;
        }

         NekDouble GaussPoints::LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts,
            const Array<OneD, const NekDouble>& funcvals)
        {
            NekDouble sum = 0.0;

            for(int i=0;i<npts;++i)
            {
                sum += funcvals[i]*LagrangePoly(x,i,npts,xpts);
            }
            return sum;
        }


        NekDouble GaussPoints::LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts)
        {
            NekDouble h=1.0;

            for(int i=0;i<pt; ++i)
            {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            for(int i=pt+1;i<npts;++i)
            {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            return h;
        }

        NekDouble GaussPoints::LagrangePolyDeriv(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts)
        {
            NekDouble h;
            NekDouble y=0.0;

            for(int j=0;j<npts;++j)
            {
                if(j!=pt)
                {
                    h=1.0;
                    for(int i=0;i<npts;++i)
                    {
                        if(i!=pt)
                        {
                            if(i!=j)
                            {
                                h = h*(x-xpts[i]);
                            }
                            h = h/(xpts[pt]-xpts[i]);
                        }
                    }
                    y = y + h;
                }
            }
            return y;
        }
        
        const boost::shared_ptr<NekMatrix<NekDouble> > GaussPoints::GetGalerkinProjection(const PointsKey &pkey)
        {
            return m_GalerkinProjectionManager[pkey];
        }

        boost::shared_ptr< NekMatrix<NekDouble> > GaussPoints::CreateGPMatrix(const PointsKey &pkey)
        {
            boost::shared_ptr< NekMatrix<NekDouble> > returnval = CalculateGalerkinProjectionMatrix(pkey);

            // Delegate to function below
            return  returnval;
        }

        boost::shared_ptr<NekMatrix<NekDouble> > GaussPoints::CalculateGalerkinProjectionMatrix(const PointsKey &pkey)
        {
            int numpointsfrom = pkey.GetNumPoints();
            int numpointsto   = GetNumPoints();
            
            Array<OneD, const NekDouble> weightsfrom;
            
            weightsfrom = PointsManager()[pkey]->GetW();
                                    
            boost::shared_ptr< NekMatrix<NekDouble> > Interp = GetI(pkey);

            Array<OneD, NekDouble> GalProj(numpointsfrom*numpointsto);
            
            // set up inner product matrix and multiply by inverse of
            // diagaonal mass matrix
            for(int i = 0; i < numpointsto; ++i)
            {
                Vmath::Vmul(numpointsfrom,Interp->GetPtr().get() +i*numpointsfrom,1,
                            &weightsfrom[0],1,&GalProj[0] +i,numpointsto);
                Vmath::Smul(numpointsfrom,1.0/m_weights[i],&GalProj[0]+i,numpointsto,
                            &GalProj[0]+i,numpointsto);
            }

            
            NekDouble* t = GalProj.data();
            boost::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpointsto,numpointsfrom,t));

            return returnval;
        }

    } // end of namespace LibUtilities
} // end of namespace Nektar

