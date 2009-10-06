///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.cpp
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
// Description: SegExp routines
//
///////////////////////////////////////////////////////////////////////////////
#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/LocalRegions.hpp>
#include <LocalRegions/SegExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {

        // constructor
        SegExp::SegExp(const LibUtilities::BasisKey &Ba, 
                       const SpatialDomains::Geometry1DSharedPtr &geom):
            StdRegions::StdSegExp(Ba),
            m_geom(geom),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {         
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,  StdRegions::eNoExpansionType,*this),boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i, StdRegions::eNoExpansionType,*this), boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }

            GenMetricInfo();
        }

        SegExp::SegExp(const LibUtilities::BasisKey &Ba):
            StdRegions::StdSegExp(Ba),
            m_geom(),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                          boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors. ; 
            Array<OneD,NekDouble> ndata(1,1.0); 
            m_metricinfo->ResetGmat(ndata,1,1,1);
            m_metricinfo->ResetJac(1,ndata);
        }

        // copy constructor
        SegExp::SegExp(const SegExp &S):
            StdRegions::StdSegExp(S),
            m_geom(S.m_geom),
            m_metricinfo(S.m_metricinfo),
            m_matrixManager(std::string("SegExpMatrix")),
            m_staticCondMatrixManager(std::string("SegExpStaticCondMatrix"))
        {      
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                boost::bind(&SegExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this),  
                                                          boost::bind(&SegExp::CreateStaticCondMatrix, this, _1));
            }       
        }

        // by default the StdExpansion1D destructor will be called

        SegExp::~SegExp()
        {
        }

        // interpolate and possibly generate geometric factors. 
        void SegExp::GenMetricInfo()
        {
            m_metricinfo = m_geom->GetGeomFactors(m_base);
//             SpatialDomains::GeomFactorsSharedPtr Xgfac = m_geom->GetGeomFactors();

//             if( Xgfac->GetGtype() != SpatialDomains::eDeformed)
//             {
//                 m_metricinfo = Xgfac;
//             }
//             else
//             {                
//                 int coordim = m_geom->GetCoordim();
//                 int expdim  = 1;
//                 int  nq = m_base[0]->GetNumPoints();
//                 SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;
//                 Array<OneD,NekDouble> ndata(coordim*nq); 
//                 Array<OneD, const NekDouble> ojac  = Xgfac->GetJac();
//                 Array<TwoD, const NekDouble> gmat  = Xgfac->GetGmat();;

//                 // assume all directiosn of geombasis are same
//                 LibUtilities::BasisSharedPtr CBasis0 = m_geom->GetBasis(0,0); 

//                 m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
//                     AllocateSharedPtr(gtype,expdim,coordim); 

//                 // check to see if basis are different distributions
//                 if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey())))
//                 {
//                     int i,j;
//                     int cnq = CBasis0->GetNumPoints();
//                     Array<OneD, NekDouble> ctmp(cnq);
//                     Array<OneD, NekDouble> tmp(nq);
//                     Array<OneD, NekDouble> jac(nq,0.0);

//                     // interpolate Geometric data   
//                     for(i = 0; i < coordim; ++i)
//                     {
//                         // Calculate  lcoal derivatives for interpolation
//                         // Need check in case points are aligned to x,y axis
//                         for(j = 0; j < cnq; ++j)
//                         {
//                             ctmp[j] = (fabs(gmat[i][j]) > NekConstants::kNekZeroTol)? 1.0/gmat[i][j]: 0.0;
//                         }

//                         LibUtilities::Interp1D(CBasis0->GetBasisKey(),  ctmp, 
//                                  m_base[0]->GetBasisKey(), tmp);
     
//                         for(j = 0; j < nq; ++j)
//                         {
//                             ndata[i*nq+j] = (fabs(tmp[j]) > NekConstants::kNekZeroTol)? 1.0/tmp[j]: 0.0;
//                         }

//                         Vmath::Vvtvp(nq,tmp,1,tmp,1,jac,1,jac,1);

//                     }
//                     Vmath::Vsqrt(nq,jac,1,jac,1);

//                     m_metricinfo->ResetGmat(ndata,nq,1,coordim);
//                     m_metricinfo->ResetJac(nq,jac);                    
//                 }
//                 else  // Same data can be used 
//                 {
//                     // Copy Geometric data
//                     Blas::Dcopy(nq*coordim,&gmat[0][0],1,&ndata[0],1);

//                     m_metricinfo->ResetGmat(ndata,nq,1,coordim);

//                     // Copy Jacobian
//                     ndata = Array<OneD,NekDouble>(nq);    
//                     Blas::Dcopy(nq,ojac,1,ndata,1);

//                     m_metricinfo->ResetJac(nq,ndata);                    
//                 }   
//             }
        }


        //----------------------------
        // Integration Methods
        //----------------------------

        /** \brief Integrate the physical point list \a inarray over region
            and return the value

            Inputs:\n

            - \a inarray: definition of function to be returned at
            quadrature point of expansion.

            Outputs:\n

            - returns \f$\int^1_{-1} u(\xi_1)d \xi_1 \f$ where \f$inarray[i]
            = u(\xi_{1i}) \f$
        */

        NekDouble SegExp::Integral(const Array<OneD, const NekDouble>&  inarray)
        {

            int    nquad0 = m_base[0]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            NekDouble  ival;
            Array<OneD,NekDouble> tmp(nquad0);

            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp,1);
            }
            else
            {
                Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
            }

            // call StdSegExp version;
            ival = StdSegExp::Integral(tmp);
            return ival;
        }


        void SegExp::IProductWRTBase(const Array<OneD, const NekDouble>& base, 
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &outarray, 
                                     int coll_check)
        {
            int   nquad0 = m_base[0]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp(nquad0);


            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0, jac[0], inarray, 1, tmp, 1);
            }
            StdSegExp::IProductWRTBase(base,tmp,outarray,coll_check);
        }

        void SegExp::IProductWRTDerivBase(const int dir, 
                                          const Array<OneD, const NekDouble>& inarray, 
                                          Array<OneD, NekDouble> & outarray)
        {
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD, NekDouble> tmp1(nquad);
            
            switch(dir)
            {
            case 0:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[0],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[0][0], inarray, 1, tmp1, 1);
                    }                 
                }
                break;
            case 1:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[1],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[1][0], inarray, 1, tmp1, 1);
                    }
                }
                break;
            case 2:
                {
                    ASSERTL1(m_geom->GetCoordim() == 3,"input dir is out of range");
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,gmat[2],1,inarray,1,tmp1,1);
                    }
                    else
                    {
                        Vmath::Smul(nquad, gmat[2][0], inarray, 1, tmp1, 1);
                    }
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }   
            IProductWRTBase(m_base[0]->GetDbdata(),tmp1,outarray,1);       
        }

        void SegExp::LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,NekDouble> &outarray)
        {
            int    i;
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);

            BwdTrans(inarray,physValues);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {            
            case 1:
                {
                    PhysDeriv(physValues,dPhysValuesdx);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                    }          
                }
                break;
            case 2:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                    }         
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);
                    Array<OneD,NekDouble> dPhysValuesdz(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz); 

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[2][0],1,dPhysValuesdz.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.get(), 1, dPhysValuesdx.get(), 1);
                    }  
                }
                break;
            default:
                ASSERTL0(false,"Wrong number of dimensions");
                break;
            }
            
            IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);      
        }


        void SegExp::HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,NekDouble> &outarray,
                                       const double lambda)
        {
            int    i;
            int    nquad = m_base[0]->GetNumPoints();
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> physValues(nquad);
            Array<OneD,NekDouble> dPhysValuesdx(nquad);
            Array<OneD,NekDouble> wsp(m_ncoeffs);

            BwdTrans(inarray,physValues);

            // mass matrix operation
            IProductWRTBase((m_base[0]->GetBdata()),physValues,wsp,1);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {            
            case 1:
                {
                    PhysDeriv(physValues,dPhysValuesdx);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                    }          
                }
                break;
            case 2:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                    }         
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdy(nquad);
                    Array<OneD,NekDouble> dPhysValuesdz(nquad);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz); 

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nquad,&gmat[0][0],1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[1][0],1,dPhysValuesdy.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                        Vmath::Vvtvp(nquad,&gmat[2][0],1,dPhysValuesdz.get(),1,dPhysValuesdx.get(),1,dPhysValuesdx.get(),1);
                    }
                    else
                    {
                        Blas::Dscal(nquad, gmat[0][0], dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[1][0], dPhysValuesdy.get(), 1, dPhysValuesdx.get(), 1);
                        Blas::Daxpy(nquad, gmat[2][0], dPhysValuesdz.get(), 1, dPhysValuesdx.get(), 1);
                    }  
                }
                break;
            default:
                ASSERTL0(false,"Wrong number of dimensions");
                break;
            }
            
            IProductWRTBase(m_base[0]->GetDbdata(),dPhysValuesdx,outarray,1);
            Blas::Daxpy(m_ncoeffs, lambda, wsp.get(), 1, outarray.get(), 1);            
        }

        //----------------------------
        // Differentiation Methods
        //-----------------------------

        /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
            physical quadrature points given by \a inarray and return in \a
            outarray.

            This is a wrapper around StdExpansion1D::Tensor_Deriv

            Input:\n

            - \a n: number of derivatives to be evaluated where \f$ n \leq  dim\f$

            - \a inarray: array of function evaluated at the quadrature points

            Output: \n

            - \a outarray: array of the derivatives \f$
            du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dx, 
            du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dy, 
            du/d_{\xi_1}|_{\xi_{1i}} d\xi_1/dz, 
            \f$ depending on value of \a dim
        */


        void SegExp::PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                               Array<OneD,NekDouble> &out_d0,
                               Array<OneD,NekDouble> &out_d1,
                               Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            Array<TwoD, const NekDouble>  gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> diff(nquad0);

            StdExpansion1D::PhysTensorDeriv(inarray,diff);

            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[0][0],1,&diff[0],1,
                                &out_d0[0],1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[1][0],1,&diff[0],1,
                                &out_d1[0],1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Vmul(nquad0,&gmat[2][0],1,&diff[0],1,
                                &out_d2[0],1);
                }
            }
            else 
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[0][0], diff, 1,
                                out_d0, 1);
                }

                if(out_d1.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[1][0], diff, 1,
                                out_d1, 1);
                }

                if(out_d2.num_elements())
                {
                    Vmath::Smul(nquad0, gmat[2][0], diff, 1,
                                out_d2, 1);
                }
            } 
        } 
        
        void SegExp::PhysDeriv(const int dir, 
                               const Array<OneD, const NekDouble>& inarray,
                               Array<OneD, NekDouble> &outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysDeriv(inarray, outarray, NullNekDouble1DArray, NullNekDouble1DArray);   
                }
                break;
            case 1:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, outarray, NullNekDouble1DArray);   
                }
                break;
            case 2:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, NullNekDouble1DArray, outarray);   
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }             
        }


        //----------------------------
        // Evaluation Methods
        //----------------------------


        void SegExp::SetCoeffsToOrientation(StdRegions::EdgeOrientation dir)
        {
            SetCoeffsToOrientation(dir,m_coeffs,m_coeffs);
        }

        void SegExp::SetCoeffsToOrientation(StdRegions::EdgeOrientation dir, 
                                            Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
        {

            if(dir == StdRegions::eBackwards)
            {
                if(&inarray[0] != &outarray[0])
                {
                    Array<OneD,NekDouble> intmp (inarray);
                    ReverseCoeffsAndSign(intmp,outarray);
                }
                else
                {
                    ReverseCoeffsAndSign(inarray,outarray);
                }
            }
        }

        // Reverse the coefficients in a boundary interior expansion
        // this routine is of use when we need the segment
        // coefficients corresponding to a expansion in the reverse
        // coordinate direction
        void SegExp::ReverseCoeffsAndSign(const Array<OneD,NekDouble> &inarray,
                                          Array<OneD,NekDouble> &outarray)
        {

            int m;
            NekDouble sgn = 1;

            ASSERTL1(&inarray[0] != &outarray[0],"inarray and outarray can not be the same");
            switch(GetBasisType(0))
            {
            case LibUtilities::eModified_A:
                
                //Swap vertices
                outarray[0] = inarray[1]; 
                outarray[1] = inarray[0];
                
                // negate odd modes
                for(m = 2; m < m_ncoeffs; ++m)
                {
                    outarray[m] = sgn*inarray[m];
                    sgn = -sgn;
                }
                break;
            case LibUtilities::eGLL_Lagrange:

                for(m = 0; m < m_ncoeffs; ++m)
                {
                    outarray[m_ncoeffs-1-m] = inarray[m];
                }
                break;                
            default:
                ASSERTL0(false,"This basis is not allowed in this method");
                break;
            }
        }

      /** \brief Mass inversion product from \a inarray to \a outarray
	  
	  Multiply by the inverse of the mass matrix
	  \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$
	  
      */ 
       
      void SegExp::MultiplyByElmtInvMass(const Array<OneD, const NekDouble>& inarray, 
					Array<OneD,NekDouble> &outarray)
        {
	  // get Mass matrix inverse
	  MatrixKey             masskey(StdRegions::eInvMass, DetExpansionType(),*this);
	  DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];
          
	  NekVector<const NekDouble> in(m_ncoeffs,inarray,eCopy);
	  NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
          
	  out = (*matsys)*in;
	}

        /** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a outarray

            Perform a forward transform using a Galerkin projection by
            taking the inner product of the physical points and multiplying
            by the inverse of the mass matrix using the Solve method of the
            standard matrix container holding the local mass matrix, i.e.
            \f$ {\bf \hat{u}} = {\bf M}^{-1} {\bf I} \f$ where \f$ {\bf I}[p] =
            \int^1_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1 \f$

            Inputs:\n

            - \a inarray: array of physical quadrature points to be transformed

            Outputs:\n

            - \a outarray: updated array of expansion coefficients. 

        */ 


        // need to sort out family of matrices 
        void SegExp::FwdTrans(const Array<OneD, const NekDouble>& inarray, 
                              Array<OneD,NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else 
            {
                IProductWRTBase(inarray,outarray);

                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass, DetExpansionType(),*this);
                DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }

        void SegExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> &outarray)
        {
            if(m_base[0]->Collocation())
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                int nInteriorDofs = m_ncoeffs-2;
                int offset;

                switch(m_base[0]->GetBasisType())
                {
                case LibUtilities::eGLL_Lagrange:
                    {
                        offset = 1;
                    }
                    break;
                case LibUtilities::eModified_A:
                case LibUtilities::eModified_B:
                    {
                        offset = 2;
                    }
                    break;
                default:
                    ASSERTL0(false,"This type of FwdTrans is not defined for this expansion type");
                }    

                fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

                outarray[GetVertexMap(0)] = inarray[0];
                outarray[GetVertexMap(1)] = inarray[m_base[0]->GetNumPoints()-1];

                if(m_ncoeffs>2)
                {
                    Array<OneD, NekDouble> tmp0(m_ncoeffs); //ideally, we would like to have tmp0 to be replaced by outarray (currently MassMatrixOp does not allow aliasing)
                    Array<OneD, NekDouble> tmp1(m_ncoeffs);
                    
                    StdRegions::StdMatrixKey  stdmasskey(StdRegions::eMass,DetExpansionType(),*this);
                    MassMatrixOp(outarray,tmp0,stdmasskey);
                    IProductWRTBase(inarray,tmp1);

                    Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);
                    
                    // get Mass matrix inverse (only of interior DOF)
                    MatrixKey             masskey(StdRegions::eMass, DetExpansionType(),*this);
                    DNekScalMatSharedPtr  matsys = (m_staticCondMatrixManager[masskey])->GetBlock(1,1);
                    
                    Blas::Dgemv('N',nInteriorDofs,nInteriorDofs, matsys->Scale(), &((matsys->GetOwnedMatrix())->GetPtr())[0],
                                nInteriorDofs,tmp1.get()+offset,1,0.0,outarray.get()+offset,1);       
                }
            }
        }

        void SegExp::GetCoords(Array<OneD, NekDouble> &coords_0,
                               Array<OneD, NekDouble> &coords_1,
                               Array<OneD, NekDouble> &coords_2)
        { 
            Array<OneD,NekDouble>  x;

            LibUtilities::BasisSharedPtr CBasis; 
            ASSERTL0(m_geom, "m_geom not defined");

            // get physical points defined in Geom
            m_geom->FillGeom();

            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0, "output coords_2 is not defined");
                CBasis = m_geom->GetBasis(2,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_2, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
                             m_base[0]->GetBasisKey(),&coords_2[0]);
                }
            case 2:
                ASSERTL0(coords_1.num_elements() != 0, 
                         "output coords_1 is not defined");
                CBasis = m_geom->GetBasis(1,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_1, 1);
                }
                else // LibUtilities::Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(1))[0],
                             m_base[0]->GetBasisKey(),&coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements() != 0, 
                         "output coords_2 is not defined");
                CBasis = m_geom->GetBasis(0,0);

                if(m_base[0]->GetBasisKey().SamePoints(CBasis->GetBasisKey()))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints(), x, 1, coords_0, 1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp1D(CBasis->GetBasisKey(),&(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetBasisKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
            }
        }


        void SegExp::GetCoord(const Array<OneD, const NekDouble>& Lcoords, 
                              Array<OneD,NekDouble> &coords)
        {
            int  i;

            ASSERTL1(Lcoords[0] >= -1.0&& Lcoords[0] <= 1.0,
                     "Local coordinates are not in region [-1,1]");

            m_geom->FillGeom();
            for(i = 0; i < m_geom->GetCoordim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }     
        }

        void SegExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
        {
            if(format==eTecplot)
            {
                int i,j;
                int     nquad = m_base[0]->GetNumPoints();
                
                Array<OneD,NekDouble> coords[3];
                
                ASSERTL0(m_geom,"m_geom not defined");
                
                int     coordim  = m_geom->GetCoordim();
                
                coords[0] = Array<OneD,NekDouble>(nquad);
                coords[1] = Array<OneD,NekDouble>(nquad);
                coords[2] = Array<OneD,NekDouble>(nquad);
                
                GetCoords(coords[0],coords[1],coords[2]);
                
                if(dumpVar)
                {
                    outfile << "Variables = x";
                    
                    if(coordim == 2)
                    {
                        
                        outfile << ", y";
                    }
                    else if (coordim == 3)
                    {
                        outfile << ", y, z";
                    }
                    outfile << ", "<< var << std::endl << std::endl;
                }
                
                outfile << "Zone, I=" << nquad <<", F=Point" << std::endl;
                
                for(i = 0; i < nquad; ++i)
                {
                    for(j = 0; j < coordim; ++j)
                    {
                        outfile << coords[j][i] << " ";
                    }
                    outfile << m_phys[i] << std::endl;
                }
            }
            else if(format==eGmsh)
            {  
                int i,j;
                int     nquad = m_base[0]->GetNumPoints();
                
                Array<OneD,NekDouble> coords[3];
                
                ASSERTL0(m_geom,"m_geom not defined");             
                coords[0] = Array<OneD,NekDouble>(nquad,0.0);
                coords[1] = Array<OneD,NekDouble>(nquad,0.0);
                coords[2] = Array<OneD,NekDouble>(nquad,0.0);
                
                GetCoords(coords[0],coords[1],coords[2]);

                if(dumpVar)
                {
                    outfile<<"View.Type = 2;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
 
                for(i = 0; i < nquad; ++i)
                {
                    outfile << "SP(" << coords[0][i] << ", " << coords[1][i] << ", " << coords[2][i] << ")";
                    outfile << "{" << m_phys[i] << "};" << endl;
                }

                if(dumpVar)
                { 
                    outfile << "};" << endl;
                }
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }

        DNekMatSharedPtr SegExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey = m_base[0]->GetBasisKey();
            StdRegions::StdSegExpSharedPtr tmp = MemoryManager<StdSegExp>::AllocateSharedPtr(bkey);
            
            return tmp->GetStdMatrix(mkey); 
        }

        NekDouble SegExp::PhysEvaluate(const Array<OneD, const NekDouble>& coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(1);

            ASSERTL0(m_geom,"_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);

            return StdSegExp::PhysEvaluate(Lcoord);
        }

        DNekScalMatSharedPtr SegExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            NekDouble fac;

            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0;
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        fac = (m_metricinfo->GetJac())[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eInvMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetExpansionType(),
                                                         *this);
                        DNekMatSharedPtr mat = GenMatrix(masskey);
                        mat->Invert();

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        fac = 1.0/(m_metricinfo->GetJac())[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0; 
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        int dir;
                        switch(mkey.GetMatrixType())
                        {
                        case StdRegions::eWeakDeriv0:
                            dir = 0;
                            break;
                        case StdRegions::eWeakDeriv1:
                            ASSERTL1(m_geom->GetCoordim() >= 2,"Cannot call eWeakDeriv2 in a coordinate system which is not at least two-dimensional");
                            dir = 1;
                            break;
                        case StdRegions::eWeakDeriv2:
                            ASSERTL1(m_geom->GetCoordim() == 3,"Cannot call eWeakDeriv2 in a coordinate system which is not three-dimensional");
                            dir = 2;
                            break;
                        }
                        fac = m_metricinfo->GetGmat()[dir][0]*
                            m_metricinfo->GetJac()[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        fac = 1.0;
                        goto UseLocRegionsMatrix;
                    }
                    else
                    {
                        int coordim = m_geom->GetCoordim();
                        fac = 0.0;
                        for(int i = 0; i < coordim; ++i)
                        {
                            fac += m_metricinfo->GetGmat()[i][0]*
                                m_metricinfo->GetGmat()[i][0];
                        }
                        fac *= m_metricinfo->GetJac()[0];
                        goto UseStdRegionsMatrix;
                    }
                }
                break;
            case StdRegions::eHelmholtz:
                {
                    NekDouble factor = mkey.GetConstant(0);
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetExpansionType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetExpansionType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                    
                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                }
                break;
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    NekDouble one    = 1.0;
                    
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    StdRegions::StdMatrixKey hkey(StdRegions::eHybridDGHelmholtz,
                                                  DetExpansionType(),*this,
                                                  mkey.GetConstant(0),
                                                  mkey.GetConstant(1));
                    DNekMatSharedPtr mat = GenMatrix(hkey);

                    mat->Invert();
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
                break;
            UseLocRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            UseStdRegionsMatrix:
                {
                    DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }
        
        DNekMatSharedPtr SegExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHybridDGHelmholtz:
            case StdRegions::eHybridDGLamToU:
            case StdRegions::eHybridDGLamToQ0:
            case StdRegions::eHybridDGLamToQ1:
            case StdRegions::eHybridDGLamToQ2:
            case StdRegions::eHybridDGHelmBndLam:
                returnval = Expansion1D::GenMatrix(mkey);
                break;
            default:
                returnval = StdSegExp::GenMatrix(mkey);
                break;
            }
            
            return returnval;
        }


        DNekScalBlkMatSharedPtr SegExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;
            
            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            Array<OneD, unsigned int> exp_size(2);
            exp_size[0] = nbdry;
            exp_size[1] = nint;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size); //Really need a constructor which takes Arrays
            NekDouble factor = 1.0;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eHelmholtz: // special case since Helmholtz not defined in StdRegions

                // use Deformed case for both regular and deformed geometries 
                factor = 1.0;
                goto UseLocRegionsMatrix;
                break;
            default:
                if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                {
                    factor = 1.0;
                    goto UseLocRegionsMatrix;
                }
                else
                {
                    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr& mat = GetStdStaticCondMatrix(*(mkey.GetStdMatKey()));                    
                    DNekScalMatSharedPtr Atmp;
                    DNekMatSharedPtr     Asubmat;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(0,0)));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,Asubmat = mat->GetBlock(0,1)));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,Asubmat = mat->GetBlock(1,0)));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,Asubmat = mat->GetBlock(1,1)));
                }
                break;
            UseLocRegionsMatrix:
                {
                    int i,j;
                    int cnt = 0;
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);

                    Array<OneD,unsigned int> bmap(nbdry);
                    Array<OneD,unsigned int> imap(nint);
                    GetBoundaryMap(bmap);
                    GetInteriorMap(imap);

                    for(i = 0; i < nbdry; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*A)(i,j) = mat(bmap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nint; ++j)
                        {
                            (*B)(i,j) = mat(bmap[i],imap[j]);
                        }
                    }
                    
                    for(i = 0; i < nint; ++i)
                    {
                        for(j = 0; j < nbdry; ++j)
                        {
                            (*C)(i,j) = mat(imap[i],bmap[j]);
                        }
                        
                        for(j = 0; j < nint; ++j)
                        {
                            (*D)(i,j) = mat(imap[i],imap[j]);
                        }
                    }
                    
                    // Calculate static condensed system 
                    if(nint)
                    {
                        D->Invert();
                        (*B) = (*B)*(*D);
                        (*A) = (*A) - (*B)*(*C);
                    }
                    
                    DNekScalMatSharedPtr     Atmp;

                    returnval->SetBlock(0,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,A));
                    returnval->SetBlock(0,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,B));
                    returnval->SetBlock(1,0,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(factor,C));
                    returnval->SetBlock(1,1,Atmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(invfactor,D));
                }
            }


            
            return returnval;
        }

    } // end of namespace    
}//end of namespace

// $Log: SegExp.cpp,v $
// Revision 1.61  2009/04/29 11:18:09  pvos
// Made demo codes working with nodal expansions
//
// Revision 1.60  2009/04/27 21:34:07  sherwin
// Updated WriteToField
//
// Revision 1.59  2009/03/04 14:17:38  pvos
// Removed all methods that take and Expansion as argument
//
// Revision 1.58  2009/01/21 16:59:57  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.57  2008/12/18 14:08:24  pvos
// NekConstants update
//
// Revision 1.56  2008/11/05 16:08:15  pvos
// Added elemental optimisation functionality
//
// Revision 1.55  2008/09/23 18:20:25  pvos
// Updates for working ProjectContField3D demo
//
// Revision 1.54  2008/09/09 15:05:09  sherwin
// Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
//
// Revision 1.53  2008/08/14 22:12:56  sherwin
// Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
//
// Revision 1.52  2008/07/29 22:25:35  sherwin
// general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
//
// Revision 1.51  2008/07/19 21:15:38  sherwin
// Removed MapTo function, made orientation anticlockwise, changed enum from BndSys to BndLam
//
// Revision 1.50  2008/07/09 11:44:49  sherwin
// Replaced GetScaleFactor call with GetConstant(0)
//
// Revision 1.49  2008/07/04 10:19:05  pvos
// Some updates
//
// Revision 1.48  2008/07/02 14:09:18  pvos
// Implementation of HelmholtzMatOp and LapMatOp on shape level
//
// Revision 1.47  2008/06/06 14:57:51  pvos
// Minor Updates
//
// Revision 1.46  2008/06/05 20:18:38  ehan
// Fixed undefined function GetGtype() in the ASSERTL2().
//
// Revision 1.45  2008/06/02 23:35:18  ehan
// Fixed warning : no new line at end of file
//
// Revision 1.44  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.43  2008/05/29 21:33:37  pvos
// Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
//
// Revision 1.42  2008/05/29 01:02:13  bnelson
// Added precompiled header support.
//
// Revision 1.41  2008/05/14 18:06:50  sherwin
// mods to fix Seggeom to Geometry1D casting
//
// Revision 1.40  2008/05/10 18:27:33  sherwin
// Modifications necessary for QuadExp Unified DG Solver
//
// Revision 1.39  2008/04/06 05:59:05  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.38  2008/04/02 22:19:26  pvos
// Update for 2D local to global mapping
//
// Revision 1.37  2008/03/18 14:12:53  pvos
// Update for nodal triangular helmholtz solver
//
// Revision 1.36  2008/03/12 15:24:29  pvos
// Clean up of the code
//
// Revision 1.35  2008/02/28 10:04:11  sherwin
// Modes for UDG codes
//
// Revision 1.34  2008/01/21 19:59:32  sherwin
// Updated to take SegGeoms instead of EdgeComponents
//
// Revision 1.33  2007/12/17 13:04:30  sherwin
// Modified GenMatrix to take a StdMatrixKey and removed m_constant from MatrixKey
//
// Revision 1.32  2007/12/06 22:49:09  pvos
// 2D Helmholtz solver updates
//
// Revision 1.31  2007/11/20 16:28:46  sherwin
// Added terms for UDG Helmholtz solver
//
// Revision 1.30  2007/10/04 12:10:03  sherwin
// Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
//
// Revision 1.29  2007/10/03 11:37:50  sherwin
// Updates relating to static condensation implementation
//
// Revision 1.28  2007/08/10 03:38:08  jfrazier
// Updated with new rev of NekManager.
//
// Revision 1.27  2007/07/30 21:00:06  sherwin
// Temporary fix in CreateMatrix
//
// Revision 1.26  2007/07/28 05:09:33  sherwin
// Fixed version with updated MemoryManager
//
// Revision 1.25  2007/07/20 00:45:51  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.24  2007/07/13 15:22:11  sherwin
// Update for Helmholtz (working without bcs )
//
// Revision 1.23  2007/07/13 09:02:22  sherwin
// Mods for Helmholtz solver
//
// Revision 1.22  2007/07/12 12:53:01  sherwin
// Updated to have a helmholtz matrix
//
// Revision 1.21  2007/07/11 19:26:04  sherwin
// update for new Manager structure
//
// Revision 1.20  2007/07/11 06:36:23  sherwin
// Updates with MatrixManager update
//
// Revision 1.19  2007/07/10 17:17:26  sherwin
// Introduced Scaled Matrices into the MatrixManager
//
// Revision 1.18  2007/06/07 15:54:19  pvos
// Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
// Also made corrections to various ASSERTL2 calls
//
// Revision 1.17  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.16  2007/05/30 15:59:05  sherwin
// Fixed bug with new definition of GetLocCoords
//
// Revision 1.15  2007/05/28 08:35:25  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.14  2007/05/27 16:10:28  bnelson
// Update to new Array type.
//
// Revision 1.13  2007/04/26 15:00:16  sherwin
// SJS compiling working version using SHaredArrays
//
// Revision 1.12  2007/04/08 03:33:30  jfrazier
// Minor reformatting and fixing SharedArray usage.
//
// Revision 1.11  2007/04/04 21:49:24  sherwin
// Update for SharedArray
//
// Revision 1.10  2007/03/25 15:48:21  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.9  2007/03/20 09:13:37  kirby
// new geomfactor routines; update for metricinfo; update style
//
// Revision 1.8  2007/03/14 21:24:07  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.7  2007/03/02 12:01:54  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.6  2006/06/01 15:27:27  sherwin
// Modifications to account for LibUtilities reshuffle
//
// Revision 1.5  2006/06/01 14:15:58  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.4  2006/05/30 14:00:03  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.3  2006/05/29 17:05:49  sherwin
// Modified to put shared_ptr around geom definitions
//
// Revision 1.2  2006/05/06 20:36:16  sherwin
// Modifications to get LocalRegions/Project1D working
//
// Revision 1.1  2006/05/04 18:58:46  kirby
// *** empty log message ***
//
// Revision 1.38  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.37  2006/03/13 18:20:33  sherwin
//
// Fixed error in ResetGmat
//
// Revision 1.36  2006/03/12 21:59:48  sherwin
//
// compiling version of LocalRegions
//
