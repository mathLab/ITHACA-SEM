///////////////////////////////////////////////////////////////////////////////
//
// File QuadExp.cpp 
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////
#include <LocalRegions/LocalRegions.h>
#include <LocalRegions/LocalRegions.hpp>
#include <stdio.h>
#include <LocalRegions/QuadExp.h>


namespace Nektar
{
    namespace LocalRegions 
    {    
        QuadExp::QuadExp(const LibUtilities::BasisKey &Ba, 
                         const LibUtilities::BasisKey &Bb, 
                         const SpatialDomains::QuadGeomSharedPtr &geom):
            StdRegions::StdQuadExp(Ba,Bb),
            m_geom(geom),
            m_metricinfo(),
            m_matrixManager(std::string("QuadExpMatrix")),
            m_staticCondMatrixManager(std::string("QuadExpStaticCondMatrix"))
        {      
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i, StdRegions::eNoExpansionType,*this), boost::bind(&QuadExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoExpansionType,*this),
                                                          boost::bind(&QuadExp::CreateStaticCondMatrix, this, _1));
            }            
            GenMetricInfo();
        }
    
    
        QuadExp::QuadExp(const LibUtilities::BasisKey &Ba,
                         const LibUtilities::BasisKey &Bb):
            StdRegions::StdQuadExp(Ba,Bb),
            m_geom(),
            m_metricinfo(MemoryManager<SpatialDomains::GeomFactors>::AllocateSharedPtr()),
            m_matrixManager(std::string("QuadExpMatrix")),
            m_staticCondMatrixManager(std::string("QuadExpStaticCondMatrix"))
        {            
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this), boost::bind(&QuadExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,StdRegions::eNoExpansionType,*this), boost::bind(&QuadExp::CreateStaticCondMatrix, this, _1));
            }
            
            // Set up unit geometric factors. 
            Array<OneD,NekDouble> ndata(4,0.0); 
            ndata[0] = ndata[3] = 1.0;
            m_metricinfo->ResetGmat(ndata,1,2,2);
            m_metricinfo->ResetJac(1,ndata);
        }

        QuadExp::QuadExp(const QuadExp &T):
            StdRegions::StdQuadExp(T),
            m_geom(T.m_geom),
            m_metricinfo(T.m_metricinfo),
            m_matrixManager(std::string("QuadExpMatrix")),
            m_staticCondMatrixManager(std::string("QuadExpStaticCondMatrix"))
        {
        }
    
        // by default the StdQuadExp destructor will be called    
        QuadExp::~QuadExp()
        {
        }
    
        void QuadExp::GenMetricInfo()
        {
            SpatialDomains::GeomFactorsSharedPtr Xgfac;
            
            Xgfac = m_geom->GetGeomFactors();
            
            if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {
                //basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(m_geom->GetBasis(0,0)->GetBasisKey()))||!(m_base[1]->GetBasisKey().SamePoints(m_geom->GetBasis(0,1)->GetBasisKey())))
                {
                    StdRegions::ExpansionType shape = StdRegions::eQuadrilateral;
                    
                    m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                        AllocateSharedPtr(shape,*Xgfac,m_base); 
                }
                else // Same data can be used 
                {                   
                    m_metricinfo = Xgfac;
                }
                
            }
        }
        
        //----------------------------        
        // Integration Methods
        //----------------------------
        
        /** \brief Integrate the physical point list \a inarray over region
            and return the value
            
            Inputs:\n
            
            - \a inarray: definition of function to be returned at quadrature point 
            of expansion. 
            
            Outputs:\n
            
            - returns \f$\int^1_{-1}\int^1_{-1} u(\xi_1, \xi_2) J[i,j] d
            \xi_1 d \xi_2 \f$ where \f$inarray[i,j] =
            u(\xi_{1i},\xi_{2j}) \f$ and \f$ J[i,j] \f$ is the
            Jacobian evaluated at the quadrature point.
        */
        NekDouble QuadExp::Integral(const Array<OneD, const NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            NekDouble ival;
            Array<OneD,NekDouble> tmp(nquad0*nquad1);
            
            // multiply inarray with Jacobian
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1, jac[0], inarray, 1, tmp, 1);
            }
            
            // call StdQuadExp version;
            ival = StdQuadExp::Integral(tmp);            
            return  ival;
        }

        void QuadExp::IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
                                      const Array<OneD, const NekDouble>& base1,
                                      const Array<OneD, const NekDouble>& inarray, 
                                      Array<OneD, NekDouble> &outarray,
                                      int coll_check)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp(nquad0*nquad1);
            
            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1, jac, 1, inarray, 1, tmp, 1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1, jac[0], inarray, 1, tmp, 1);
            }
            
            StdQuadExp::IProductWRTBase(base0,base1,tmp,outarray,coll_check);
            
        }

        void QuadExp::IProductWRTDerivBase(const int dir, 
                                           const Array<OneD, const NekDouble>& inarray, 
                                           Array<OneD, NekDouble> & outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD, NekDouble> tmp1(nqtot);
            Array<OneD, NekDouble> tmp2(nqtot);
            Array<OneD, NekDouble> tmp3(m_ncoeffs);
            
            switch(dir)
            {
            case 0:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nqtot,gmat[0],1,inarray,1,tmp1,1);
                        Vmath::Vmul(nqtot,gmat[1],1,inarray,1,tmp2,1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[0][0], inarray, 1, tmp1, 1);
                        Vmath::Smul(nqtot, gmat[1][0], inarray, 1, tmp2, 1);
                    }                 
                }
                break;
            case 1:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nqtot,gmat[2],1,inarray,1,tmp1,1);
                        Vmath::Vmul(nqtot,gmat[3],1,inarray,1,tmp2,1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[2][0], inarray, 1, tmp1, 1);
                        Vmath::Smul(nqtot, gmat[3][0], inarray, 1, tmp2, 1);
                    }
                }
                break;
            case 2:
                {
                    ASSERTL1(m_geom->GetCoordim() == 3,"input dir is out of range");
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nqtot,gmat[4],1,inarray,1,tmp1,1);
                        Vmath::Vmul(nqtot,gmat[5],1,inarray,1,tmp2,1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[4][0], inarray, 1, tmp1, 1);
                        Vmath::Smul(nqtot, gmat[5][0], inarray, 1, tmp2, 1);
                    }
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }   
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp1,tmp3,1);
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),tmp2,outarray,1);
            Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);                 
        }

        void QuadExp::LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> physValues(nqtot);
            Array<OneD,NekDouble> dPhysValuesdx(nqtot);
            Array<OneD,NekDouble> dPhysValuesdy(nqtot);

            Array<OneD,NekDouble> wsp(m_ncoeffs);
            Array<OneD,NekDouble> tmp(nqtot);

            BwdTrans(inarray,physValues);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {
            case 2:
                {
                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nqtot,&gmat[0][0],1,dPhysValuesdx.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[2][0],1,dPhysValuesdy.get(),1,tmp.get(),1,tmp.get(),1);

                        Vmath::Vmul (nqtot,&gmat[3][0],1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[1][0],1,dPhysValuesdx.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[0][0], dPhysValuesdx.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[2][0], dPhysValuesdy.get(), 1, tmp.get(), 1);

                        Blas::Dscal(nqtot, gmat[3][0], dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[1][0], dPhysValuesdx.get(), 1, dPhysValuesdy.get(), 1);
                    }          
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdz(nqtot);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nqtot,&gmat[0][0],1,dPhysValuesdx.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[2][0],1,dPhysValuesdy.get(),1,tmp.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[4][0],1,dPhysValuesdz.get(),1,tmp.get(),1,tmp.get(),1);

                        Vmath::Vmul (nqtot,&gmat[3][0],1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[1][0],1,dPhysValuesdx.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[5][0],1,dPhysValuesdz.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[0][0], dPhysValuesdx.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[2][0], dPhysValuesdy.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[4][0], dPhysValuesdz.get(), 1, tmp.get(), 1);

                        Blas::Dscal(nqtot, gmat[3][0], dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[1][0], dPhysValuesdx.get(), 1, dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[5][0], dPhysValuesdz.get(), 1, dPhysValuesdy.get(), 1);
                    }        

                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions should be greater than 2");
                break;
            }
            
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp,outarray,1);
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,wsp,1);  
            Vmath::Vadd(m_ncoeffs,wsp.get(),1,outarray.get(),1,outarray.get(),1);      
            
        }

        void QuadExp::HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray,
                                        const double lambda)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            const Array<TwoD, const NekDouble>& gmat = m_metricinfo->GetGmat();

            Array<OneD,NekDouble> physValues(nqtot);
            Array<OneD,NekDouble> dPhysValuesdx(nqtot);
            Array<OneD,NekDouble> dPhysValuesdy(nqtot);

            Array<OneD,NekDouble> wsp(m_ncoeffs);
            Array<OneD,NekDouble> tmp(nqtot);

            BwdTrans(inarray,physValues);

            // mass matrix operation
            IProductWRTBase((m_base[0]->GetBdata()),(m_base[1]->GetBdata()),
                            physValues,wsp,1);

            // Laplacian matrix operation
            switch(m_geom->GetCoordim())
            {
            case 2:
                {
                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nqtot,&gmat[0][0],1,dPhysValuesdx.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[2][0],1,dPhysValuesdy.get(),1,tmp.get(),1,tmp.get(),1);

                        Vmath::Vmul (nqtot,&gmat[3][0],1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[1][0],1,dPhysValuesdx.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[0][0], dPhysValuesdx.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[2][0], dPhysValuesdy.get(), 1, tmp.get(), 1);

                        Blas::Dscal(nqtot, gmat[3][0], dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[1][0], dPhysValuesdx.get(), 1, dPhysValuesdy.get(), 1);
                    }          
                }
                break;
            case 3:
                {
                    Array<OneD,NekDouble> dPhysValuesdz(nqtot);

                    PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy,dPhysValuesdz);

                    // multiply with the proper geometric factors
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul (nqtot,&gmat[0][0],1,dPhysValuesdx.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[2][0],1,dPhysValuesdy.get(),1,tmp.get(),1,tmp.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[4][0],1,dPhysValuesdz.get(),1,tmp.get(),1,tmp.get(),1);

                        Vmath::Vmul (nqtot,&gmat[3][0],1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[1][0],1,dPhysValuesdx.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                        Vmath::Vvtvp(nqtot,&gmat[5][0],1,dPhysValuesdz.get(),1,dPhysValuesdy.get(),1,dPhysValuesdy.get(),1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot, gmat[0][0], dPhysValuesdx.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[2][0], dPhysValuesdy.get(), 1, tmp.get(), 1);
                        Blas::Daxpy(nqtot, gmat[4][0], dPhysValuesdz.get(), 1, tmp.get(), 1);

                        Blas::Dscal(nqtot, gmat[3][0], dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[1][0], dPhysValuesdx.get(), 1, dPhysValuesdy.get(), 1);
                        Blas::Daxpy(nqtot, gmat[5][0], dPhysValuesdz.get(), 1, dPhysValuesdy.get(), 1);
                    }        

                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions should be greater than 2");
                break;
            }
            
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp,outarray,1);
            Blas::Daxpy(m_ncoeffs, lambda, wsp.get(), 1, outarray.get(), 1);

            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,wsp,1);  
            Vmath::Vadd(m_ncoeffs,wsp.get(),1,outarray.get(),1,outarray.get(),1);      
            
        }       
        
        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////
        
        /** 
            \brief Calculate the derivative of the physical points 
            
            For quadrilateral region can use the Tensor_Deriv function
            defined under StdExpansion.
            
        **/
        void QuadExp::PhysDeriv(const Array<OneD, const NekDouble> & inarray,
                                Array<OneD,NekDouble> &out_d0,
                                Array<OneD,NekDouble> &out_d1,
                                Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> diff0(nquad0*nquad1);
            Array<OneD,NekDouble> diff1(nquad0*nquad1);
            
            StdQuadExp::PhysDeriv(inarray, diff0, diff1);
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,gmat[0],1,diff0,1, out_d0, 1);
                    Vmath::Vvtvp (nquad0*nquad1,gmat[1],1,diff1,1, out_d0, 1,
                                  out_d0,1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,gmat[2],1,diff0,1, out_d1, 1);
                    Vmath::Vvtvp (nquad0*nquad1,gmat[3],1,diff1,1, out_d1, 1,
                                  out_d1,1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,gmat[4],1,diff0,1, out_d2, 1);
                    Vmath::Vvtvp (nquad0*nquad1,gmat[5],1,diff1,1, out_d2, 1,
                                  out_d2,1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())                
                {
                    Vmath::Smul (nquad0*nquad1, gmat[0][0], diff0, 1, out_d0, 1);
                    Blas::Daxpy (nquad0*nquad1, gmat[1][0], diff1, 1, out_d0, 1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Smul (nquad0*nquad1, gmat[2][0], diff0, 1, out_d1, 1);
                    Blas::Daxpy (nquad0*nquad1, gmat[3][0], diff1, 1, out_d1, 1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Smul (nquad0*nquad1, gmat[4][0], diff0,1, out_d2, 1);
                    Blas::Daxpy (nquad0*nquad1, gmat[5][0], diff1,1, out_d2, 1);
                }
            }
        }
        
        void QuadExp::PhysDeriv(const int dir, 
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
        
        /** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a (this)->m_coeffs  
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - (this)->_coeffs: updated array of expansion coefficients. 
            
        */    
        void QuadExp::FwdTrans(const Array<OneD, const NekDouble> & inarray, 
                               Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, m_coeffs, 1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);
                
                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetExpansionType(),*this);
                DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }

        void QuadExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                              Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
            {
                int i,j;
                int npoints[2] = {m_base[0]->GetNumPoints(),
                                  m_base[1]->GetNumPoints()};
                int nmodes[2]  = {m_base[0]->GetNumModes(),
                                  m_base[1]->GetNumModes()};

                fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

                Array<OneD, NekDouble> physEdge[4];
                Array<OneD, NekDouble> coeffEdge[4];
                StdRegions::EdgeOrientation orient[4];
                for(i = 0; i < 4; i++)
                {
                    physEdge[i]  = Array<OneD, NekDouble>(npoints[i%2]);
                    coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i%2]);
                    orient[i]    = GetEorient(i);
                }

                for(i = 0; i < npoints[0]; i++)
                {
                    physEdge[0][i] = inarray[i];
                    physEdge[2][i] = inarray[npoints[0]*npoints[1]-1-i];
                }

                for(i = 0; i < npoints[1]; i++)
                {
                    physEdge[1][i] = inarray[npoints[0]-1+i*npoints[0]];
                    physEdge[3][i] = inarray[(npoints[1]-1)*npoints[0]-i*npoints[0]];
                }

                for(i = 0; i < 4; i++)
                {
                    if( orient[i] == StdRegions::eBackwards )
                    {
                        reverse( (physEdge[i]).get() , (physEdge[i]).get() + npoints[i%2] );
                    }
                }

                SegExpSharedPtr segexp[4];
                for(i = 0; i < 4; i++)
                {
                    segexp[i] = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(m_base[i%2]->GetBasisKey(),GetGeom2D()->GetEdge(i));
                }

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;
                NekDouble sign;

                for(i = 0; i < 4; i++)
                {
                    segexp[i%2]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                    GetEdgeToElementMap(i,orient[i],mapArray,signArray);
                    for(j=0; j < nmodes[i%2]; j++)
                    {
                        sign = (NekDouble) signArray[j];
                        outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                    }
                }

                Array<OneD, NekDouble> tmp0(m_ncoeffs);
                Array<OneD, NekDouble> tmp1(m_ncoeffs);
                
                MassMatrixOp(outarray,tmp0);
                IProductWRTBase(inarray,tmp1);
                
                Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);
                
                // get Mass matrix inverse (only of interior DOF)
                // use block (1,1) of the static condensed system
                // note: this block alreay contains the inverse matrix
                MatrixKey             masskey(StdRegions::eMass,DetExpansionType(),*this);
                DNekScalMatSharedPtr  matsys = (m_staticCondMatrixManager[masskey])->GetBlock(1,1);

                int nBoundaryDofs = NumBndryCoeffs();
                int nInteriorDofs = m_ncoeffs - nBoundaryDofs; 

                Array<OneD, NekDouble> rhs(nInteriorDofs);
                Array<OneD, NekDouble> result(nInteriorDofs);

                GetInteriorMap(mapArray);

                for(i = 0; i < nInteriorDofs; i++)
                {
                    rhs[i] = tmp1[ mapArray[i] ];
                }

                Blas::Dgemv('N', nInteriorDofs, nInteriorDofs, matsys->Scale(), &((matsys->GetOwnedMatrix())->GetPtr())[0],
                            nInteriorDofs,rhs.get(),1,0.0,result.get(),1);   

                for(i = 0; i < nInteriorDofs; i++)
                {
                    outarray[ mapArray[i] ] = result[i];
                }
            }

        }        

        void QuadExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                                Array<OneD,NekDouble> &coords_1,
                                Array<OneD,NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            Array<OneD,NekDouble>  x;
            
            ASSERTL0(m_geom, "m_geom not defined");
            
            // get physical points defined in Geom
            m_geom->FillGeom();
            
            switch(m_geom->GetCoordim())
            {
            case 3:
                ASSERTL0(coords_2.num_elements() != 0, 
                         "output coords_2 is not defined");
                
                CBasis0 = m_geom->GetBasis(2,0); 
                CBasis1 = m_geom->GetBasis(2,1);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(2);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_2,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(),&(m_geom->UpdatePhys(2))[0], m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_2[0]);
                }
            case 2:
                ASSERTL0(coords_1.num_elements(), 
                         "output coords_1 is not defined");
                
                CBasis0 = m_geom->GetBasis(1,0); 
                CBasis1 = m_geom->GetBasis(1,1);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(1);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_1,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), &(m_geom->UpdatePhys(1))[0], 
                                           m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_1[0]);
                }
            case 1:
                ASSERTL0(coords_0.num_elements(), 
                         "output coords_0 is not defined");
                
                CBasis0 = m_geom->GetBasis(0,0); 
                CBasis1 = m_geom->GetBasis(0,1);
                
                if((m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))&&
                   (m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    x = m_geom->UpdatePhys(0);
                    Blas::Dcopy(m_base[0]->GetNumPoints()*
                                m_base[1]->GetNumPoints(),
                                x,1,coords_0,1);
                }
                else // Interpolate to Expansion point distribution
                {
                    LibUtilities::Interp2D(CBasis0->GetPointsKey(), CBasis1->GetPointsKey(), &(m_geom->UpdatePhys(0))[0],
                                           m_base[0]->GetPointsKey(),m_base[1]->GetPointsKey(),&coords_0[0]);
                }
                break;
            default:
                ASSERTL0(false,"Number of dimensions are greater than 2");
                break;
            }
        }
        
        // get the coordinates "coords" at the local coordinates "Lcoords"
        
        void QuadExp::GetCoord(const Array<OneD, const NekDouble> &Lcoords, 
                               Array<OneD,NekDouble> &coords)
        {
            int  i;
            
            ASSERTL1(Lcoords[0] >= -1.0 && Lcoords[1] <= 1.0 && 
                     Lcoords[1] >= -1.0 && Lcoords[1]  <=1.0,
                     "Local coordinates are not in region [-1,1]");
            
            m_geom->FillGeom();            
            for(i = 0; i < m_geom->GetCoordDim(); ++i)
            {
                coords[i] = m_geom->GetCoord(i,Lcoords);
            }
        }
               
        void QuadExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar)
        {
            if(format==eTecplot)
            {
                int i,j;
                int nquad0 = m_base[0]->GetNumPoints();
                int nquad1 = m_base[1]->GetNumPoints();
                Array<OneD,NekDouble> coords[3];
                
                ASSERTL0(m_geom,"m_geom not defined");
                
                int     coordim  = m_geom->GetCoordim();
                
                coords[0] = Array<OneD,NekDouble>(nquad0*nquad1);
                coords[1] = Array<OneD,NekDouble>(nquad0*nquad1);
                coords[2] = Array<OneD,NekDouble>(nquad0*nquad1);
                
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
                    outfile << ", v\n" << std::endl;
                }
                
                outfile << "Zone, I=" << nquad0 << ", J=" << 
                    nquad1 <<", F=Point" << std::endl;
                
                for(i = 0; i < nquad0*nquad1; ++i)
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
                if(dumpVar)
                {
                    outfile<<"View.MaxRecursionLevel = 8;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
                outfile<<"SQ("<<endl;                
                // write the coordinates of the vertices of the quadrilateral
                Array<OneD,NekDouble> coordVert1(2);
                Array<OneD,NekDouble> coordVert2(2);
                Array<OneD,NekDouble> coordVert3(2);
                Array<OneD,NekDouble> coordVert4(2);
                coordVert1[0]=-1.0;
                coordVert1[1]=-1.0;
                coordVert2[0]=1.0;
                coordVert2[1]=-1.0;
                coordVert3[0]=1.0;
                coordVert3[1]=1.0;
                coordVert4[0]=-1.0;
                coordVert4[1]=1.0;
                outfile<<m_geom->GetCoord(0,coordVert1)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert1)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert2)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert2)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert3)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert3)<<", 0.0,"<<endl;
                outfile<<m_geom->GetCoord(0,coordVert4)<<", ";
                outfile<<m_geom->GetCoord(1,coordVert4)<<", 0.0"<<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j,k;

                int nModes0 = m_base[0]->GetNumModes();
                int nModes1 = m_base[1]->GetNumModes();

                const LibUtilities::PointsKey Pkey1Gmsh(nModes0,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::PointsKey Pkey2Gmsh(nModes1,LibUtilities::eGaussLobattoLegendre);
                const LibUtilities::BasisKey  Bkey1Gmsh(m_base[0]->GetBasisType(),nModes0,Pkey1Gmsh);
                const LibUtilities::BasisKey  Bkey2Gmsh(m_base[1]->GetBasisType(),nModes1,Pkey2Gmsh);

                StdRegions::StdQuadExpSharedPtr EGmsh;
                EGmsh = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(Bkey1Gmsh,Bkey2Gmsh);

                int nMonomialPolynomials = EGmsh->GetNcoeffs();
              
                Array<OneD,NekDouble> xi1(nMonomialPolynomials);
                Array<OneD,NekDouble> xi2(nMonomialPolynomials);
              
                Array<OneD,NekDouble> x(nMonomialPolynomials);
                Array<OneD,NekDouble> y(nMonomialPolynomials);

                EGmsh->GetCoords(xi1,xi2);
              
                for(i=0;i<nMonomialPolynomials;i++)
                {
                    x[i] = xi1[i];//0.5*(1.0+xi1[i]);
                    y[i] = xi2[i];//0.5*(1.0+xi2[i]);
                }

                int cnt  = 0;
                Array<TwoD, int> exponentMap(nMonomialPolynomials,3,0);
                for(i = 0; i < nModes1; i++)
                {
                    for(j = 0; j < nModes0; j++)
                    {
                        exponentMap[cnt][0] = j;
                        exponentMap[cnt++][1] = i;
                    }         
                }
              
                NekMatrix<NekDouble> vdm(nMonomialPolynomials,nMonomialPolynomials);
                for(i = 0 ; i < nMonomialPolynomials; i++)
                {
                    for(j = 0 ; j < nMonomialPolynomials; j++)
                    {
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*pow(y[i],exponentMap[j][1]);
                    }
                }

                vdm.Invert();

                Array<OneD,NekDouble> rhs(nMonomialPolynomials);
                EGmsh->BwdTrans(m_coeffs,rhs);

                NekVector<const NekDouble> in(nMonomialPolynomials,rhs,eWrapper);
                NekVector<NekDouble> out(nMonomialPolynomials);
                out = vdm*in;

                //write the coefficients
                outfile<<"{";
                for(i = 0; i < nMonomialPolynomials; i++)
                {
                    outfile<<out[i];
                    if(i < nMonomialPolynomials - 1)
                    {
                        outfile<<", ";
                    }
                }
                outfile<<"};"<<endl;
              
                if(dumpVar)
                {   
                    outfile<<"INTERPOLATION_SCHEME"<<endl;
                    outfile<<"{"<<endl;
                    for(i=0; i < nMonomialPolynomials; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < nMonomialPolynomials; j++)
                        {
                            if(i==j)
                            {
                                outfile<<"1.00";
                            }
                            else
                            {
                                outfile<<"0.00";
                            }
                            if(j < nMonomialPolynomials - 1)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nMonomialPolynomials - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"}"<<endl;
                        }
                    }
                    
                    outfile<<"{"<<endl;
                    for(i=0; i < nMonomialPolynomials; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < 3; j++)
                        {
                            outfile<<exponentMap[i][j];
                            if(j < 2)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nMonomialPolynomials - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"};"<<endl;
                        }
                    }
                    outfile<<"};"<<endl;
                }                    
            } 
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }
        
        DNekMatSharedPtr QuadExp::CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
            LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
            StdRegions::StdQuadExpSharedPtr tmp = MemoryManager<StdQuadExp>::AllocateSharedPtr(bkey0,bkey1);            
            return tmp->GetStdMatrix(mkey); 
        }
        
        NekDouble QuadExp::PhysEvaluate(const Array<OneD, const NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);
            
            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);
        
            return StdQuadExp::PhysEvaluate(Lcoord);
        }
        
        DNekMatSharedPtr QuadExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
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
                returnval = Expansion2D::GenMatrix(mkey);
                break;
            default:
                returnval = StdQuadExp::GenMatrix(mkey);
            }
               
            return returnval;            
        }
        
        
        DNekScalMatSharedPtr QuadExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");
            
            switch(mkey.GetMatrixType())
            {
            case StdRegions::eMass:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {   
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,mat);
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
                        NekDouble fac = 1.0/(m_metricinfo->GetJac())[0];
                        DNekMatSharedPtr mat = GetStdMatrix(*mkey.GetStdMatKey());
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(fac,mat);                        
                    }
                }
                break;
            case StdRegions::eWeakDeriv0:
            case StdRegions::eWeakDeriv1:
            case StdRegions::eWeakDeriv2:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                        
                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();
                        int dir;

                        switch(mkey.GetMatrixType())
                        {
                        case StdRegions::eWeakDeriv0:
                            dir = 0;
                            break;
                        case StdRegions::eWeakDeriv1:
                            dir = 1;
                            break;
                        case StdRegions::eWeakDeriv2:
                            dir = 2;
                            break;
                        }                            

                        MatrixKey deriv0key(StdRegions::eWeakDeriv0,
                                            mkey.GetExpansionType(), *this);  
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetExpansionType(), *this);

                        DNekMat &deriv0 = *GetStdMatrix(*deriv0key.GetStdMatKey());
                        DNekMat &deriv1 = *GetStdMatrix(*deriv1key.GetStdMatKey());
                        
                        int rows = deriv0.GetRows();
                        int cols = deriv1.GetColumns();

                        DNekMatSharedPtr WeakDeriv = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        (*WeakDeriv) = gmat[2*dir][0]*deriv0 + gmat[2*dir+1][0]*deriv1;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,WeakDeriv);
                    }
                }
                break;
            case StdRegions::eLaplacian:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        NekDouble one = 1.0;
                        DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                    }
                    else
                    {
                        ASSERTL1(m_geom->GetCoordDim() == 2,"Standard Region Laplacian is only set up for Quads in two-dimensional");
                        MatrixKey lap00key(StdRegions::eLaplacian00,
                                           mkey.GetExpansionType(), *this);  
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetExpansionType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetExpansionType(), *this);  

                        DNekMat &lap00 = *GetStdMatrix(*lap00key.GetStdMatKey());
                        DNekMat &lap01 = *GetStdMatrix(*lap01key.GetStdMatKey());
                        DNekMat &lap11 = *GetStdMatrix(*lap11key.GetStdMatKey());

                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

                        int rows = lap00.GetRows();
                        int cols = lap00.GetColumns();

                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                        (*lap) = (gmat[0][0]*gmat[0][0]+gmat[2][0]*gmat[2][0])*lap00 + 
                            (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0])*(lap01 + Transpose(lap01)) +
                            (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0])*lap11;

                        returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(jac,lap);
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
            case StdRegions::eHybridDGLamToQ1:
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
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }
                
            return returnval;
        }

        StdRegions::StdExpansion1DSharedPtr QuadExp::GetEdgeExp(int edge, bool SetUpNormals)
        {
            GenSegExpSharedPtr returnval; 
            SpatialDomains::Geometry1DSharedPtr edg = m_geom->GetEdge(edge);
            
            returnval = MemoryManager<GenSegExp>::AllocateSharedPtr(DetEdgeBasisKey(edge),edg);
            
            
            if(SetUpNormals)
            {
                int i;
                int coordim = GetCoordim();
                int npoints = returnval->GetNumPoints(0);
                StdRegions::EdgeOrientation edgedir = GetEorient(edge);

                Array<OneD,NekDouble> phys_normals = m_metricinfo->GenNormals2D(StdRegions::eQuadrilateral,edge,returnval->GetBasis(0)->GetPointsKey());

                if(edgedir == StdRegions::eBackwards)
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        for(i = 0; i < coordim; ++i)
                        {
                            Vmath::Reverse(npoints,&phys_normals[i*npoints],1,
                                           &phys_normals[i*npoints],1);
                        }
                    }
                    
                    Vmath::Neg(coordim*npoints,phys_normals,1);
                }
                
                returnval->SetPhysNormals(phys_normals);
            }
            
            return returnval;
        }



        // Get edge values from the 2D Phys space along an edge
        // following a counter clockwise edge convention for definition
        // of edgedir, Note that point distribution is given by QuadExp. 
        void QuadExp::GetEdgePhysVals(const int edge,
                                      const Array<OneD, const NekDouble> &inarray, 
                                      Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();
            
            Array<OneD,const NekDouble> e_tmp;

            StdRegions::EdgeOrientation edgedir = GetEorient(edge);
            switch(edge)
            {
            case 0:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,inarray,1,outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+(nquad0-1),-1,
                                 outarray,1);
                }
                    
                break;
            case 1:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0-1),nquad0,
                                 outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0*nquad1-1),
                                 -nquad0, outarray,1);
                }
                break; 
            case 2:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+(nquad0*nquad1-1),-1,
                                 outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+nquad0*(nquad1-1),1,
                                 outarray,1);
                }
                break; 
            case 3:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray + nquad0*(nquad1-1),
                                 -nquad0,outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,inarray,nquad0,outarray,1);
                }
                break; 
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }
        }



        // Get edge values  following counter clockwise edge
        // convention for definition of edgedir at points defined by EdgeExp.

        void QuadExp::GetEdgePhysVals(const int edge, const StdRegions::StdExpansion1DSharedPtr &EdgeExp, 
                                      const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            Array<OneD,const NekDouble> e_tmp;
            Array<OneD,NekDouble>       outtmp(max(nquad0,nquad1));


            // get points in Cartesian orientation 
            switch(edge)
            {
            case 0:
                Vmath::Vcopy(nquad0,inarray,1,outtmp,1);
                break;
            case 1:
                Vmath::Vcopy(nquad1,e_tmp = inarray+(nquad0-1),nquad0,outtmp,1);
                break; 
            case 2:
                Vmath::Vcopy(nquad0,e_tmp = inarray+nquad0*(nquad1-1),1,
                             outtmp,1);
                break; 
            case 3:
                Vmath::Vcopy(nquad1,inarray,nquad0,outtmp,1);
                break; 
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }

            // Interpolate if required 
            LibUtilities::Interp1D(m_base[edge%2]->GetPointsKey(),outtmp,
                     EdgeExp->GetBasis(0)->GetPointsKey(),outarray);
                
            //Reverse data if necessary
            if(GetCartesianEorient(edge) == StdRegions::eBackwards)
            {
                Vmath::Reverse(EdgeExp->GetNumPoints(0),&outarray[0],1,
                               &outarray[0],1);
            }

        }


        DNekScalBlkMatSharedPtr QuadExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype() == SpatialDomains::eNoGeomType,"Geometric information is not set up");

            // set up block matrix system
            int nbdry = NumBndryCoeffs();
            int nint = m_ncoeffs - nbdry;
            unsigned int exp_size[] = {nbdry,nint};
            int nblks = 2;
            returnval = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Arrays
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
                    DNekScalMatSharedPtr& mat = GetLocMatrix(mkey);
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
                    int cnt2 = 0;
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
        

    }//end of namespace
}//end of namespace

/** 
 *    $Log: QuadExp.cpp,v $
 *    Revision 1.51  2008/09/09 15:05:09  sherwin
 *    Updates related to cuved geometries. Normals have been removed from m_metricinfo and replaced with a direct evaluation call. Interp methods have been moved to LibUtilities
 *
 *    Revision 1.50  2008/08/14 22:12:56  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.49  2008/08/02 18:58:50  sherwin
 *    Version with matrix formulation of AddUDGHelmholtzEdgeTerms formulated instead of edge integrals. Will probably push routine back to StdExpansion2D at next iteration
 *
 *    Revision 1.48  2008/07/31 11:13:22  sherwin
 *    Depracated GetEdgeBasis and replaced with DetEdgeBasisKey
 *
 *    Revision 1.47  2008/07/29 22:25:34  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.46  2008/07/19 21:15:38  sherwin
 *    Removed MapTo function, made orientation anticlockwise, changed enum from BndSys to BndLam
 *
 *    Revision 1.45  2008/07/16 22:20:54  sherwin
 *    Added AddEdgeNormBoundaryInt
 *
 *    Revision 1.44  2008/07/12 19:08:29  sherwin
 *    Modifications for DG advection routines
 *
 *    Revision 1.43  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.42  2008/07/09 11:45:48  sherwin
 *    Modification to make workinvg UDG solver, added BndSysForce matrix, replaced GetScaleFactor with GetConstant(0)
 *
 *    Revision 1.41  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.40  2008/07/02 14:09:18  pvos
 *    Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 *    Revision 1.39  2008/06/06 14:57:51  pvos
 *    Minor Updates
 *
 *    Revision 1.38  2008/06/05 20:18:21  ehan
 *    Fixed undefined function GetGtype() in the ASSERTL2().
 *
 *    Revision 1.37  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.36  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.35  2008/05/29 01:02:13  bnelson
 *    Added precompiled header support.
 *
 *    Revision 1.34  2008/05/14 09:02:35  pvos
 *    fixed bug
 *
 *    Revision 1.33  2008/05/10 18:27:32  sherwin
 *    Modifications necessary for QuadExp Unified DG Solver
 *
 *    Revision 1.32  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.31  2008/04/02 22:19:26  pvos
 *    Update for 2D local to global mapping
 *
 *    Revision 1.30  2008/03/18 14:12:53  pvos
 *    Update for nodal triangular helmholtz solver
 *
 *    Revision 1.29  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.28  2008/02/28 10:04:11  sherwin
 *    Modes for UDG codes
 *
 *    Revision 1.27  2008/01/25 16:46:54  sherwin
 *    Added UDG work
 *
 *    Revision 1.26  2008/01/21 19:59:32  sherwin
 *    Updated to take SegGeoms instead of EdgeComponents
 *
 *    Revision 1.25  2007/12/17 13:04:30  sherwin
 *    Modified GenMatrix to take a StdMatrixKey and removed m_constant from MatrixKey
 *
 *    Revision 1.24  2007/12/06 22:49:08  pvos
 *    2D Helmholtz solver updates
 *
 *    Revision 1.23  2007/11/08 16:54:26  pvos
 *    Updates towards 2D helmholtz solver
 *
 *    Revision 1.22  2007/08/11 23:41:22  sherwin
 *    Various updates
 *
 *    Revision 1.21  2007/07/28 05:09:32  sherwin
 *    Fixed version with updated MemoryManager
 *
 *    Revision 1.20  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.19  2007/07/12 12:53:00  sherwin
 *    Updated to have a helmholtz matrix
 *
 *    Revision 1.18  2007/07/11 19:26:04  sherwin
 *    update for new Manager structure
 *
 *    Revision 1.17  2007/07/11 06:36:23  sherwin
 *    Updates with MatrixManager update
 *
 *    Revision 1.16  2007/07/10 17:17:25  sherwin
 *    Introduced Scaled Matrices into the MatrixManager
 *
 *    Revision 1.15  2007/06/17 19:00:45  bnelson
 *    Removed unused variables.
 *
 *    Revision 1.14  2007/06/07 15:54:19  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.13  2007/06/06 11:29:31  pvos
 *    Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
 *
 *    Revision 1.12  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.11  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.10  2007/05/31 11:38:16  pvos
 *    Updated QuadExp and TriExp
 *
 *    Revision 1.9  2007/04/26 15:00:16  sherwin
 *    SJS compiling working version using SHaredArrays
 *
 *    Revision 1.8  2006/08/05 19:03:47  sherwin
 *    Update to make the multiregions 2D expansion in connected regions work
 *
 *    Revision 1.7  2006/06/02 18:48:39  sherwin
 *    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 *    Revision 1.6  2006/06/01 15:27:27  sherwin
 *    Modifications to account for LibUtilities reshuffle
 *
 *    Revision 1.5  2006/06/01 14:15:57  sherwin
 *    Added typdef of boost wrappers and made GeoFac a boost shared pointer.
 *
 *    Revision 1.4  2006/05/30 14:00:03  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.3  2006/05/29 17:05:49  sherwin
 *    Modified to put shared_ptr around geom definitions
 *
 *    Revision 1.2  2006/05/06 20:36:16  sherwin
 *    Modifications to get LocalRegions/Project1D working
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.23  2006/03/13 19:47:54  sherwin
 *
 *    Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
 *
 *    Revision 1.22  2006/03/13 18:20:33  sherwin
 *
 *    Fixed error in ResetGmat
 *
 *    Revision 1.21  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.20  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
