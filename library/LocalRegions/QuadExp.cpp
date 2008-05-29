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
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i, StdRegions::eNoShapeType,*this), boost::bind(&QuadExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
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
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&QuadExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&QuadExp::CreateStaticCondMatrix, this, _1));
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
                int nq0 = m_base[0]->GetNumPoints();
                int nq1 = m_base[1]->GetNumPoints();
                int nq = nq0*nq1;
                int coordim = m_geom->GetCoordim();
                int expdim = 2;
                SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;

                LibUtilities::BasisSharedPtr CBasis0;
                LibUtilities::BasisSharedPtr CBasis1;
     
                Array<OneD, const NekDouble> ojac = Xgfac->GetJac();   
                Array<TwoD, const NekDouble> ogmat = Xgfac->GetGmat();
                Array<OneD,NekDouble> njac(nq);
                Array<OneD,NekDouble> ngmat(2*coordim*nq);
                
                CBasis0 = m_geom->GetBasis(0,0); // this assumes all goembasis are same
                CBasis1 = m_geom->GetBasis(0,1);
                int Cnq0 = CBasis0->GetNumPoints();
                int Cnq1 = CBasis1->GetNumPoints();

                
                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                    AllocateSharedPtr(gtype,expdim,coordim); 
                
                //basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
                   !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    int i;   

                    // interpolate Jacobian        
                    Interp2D(CBasis0->GetBasisKey(),
                             CBasis1->GetBasisKey(),
                             &ojac[0],
                             m_base[0]->GetBasisKey(),
                             m_base[1]->GetBasisKey(), 
                             &njac[0]);
                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    Array<OneD,NekDouble> dxdxi(nq);
                    for(i = 0; i < 2*coordim; ++i)
                    {
                        Vmath::Vmul(nq,&ojac[0],1,&ogmat[i][0],1,&dxdxi[0],1);
                        //Vmath::Vmul(nq,ogmat[i],1,ojac,1,dxdxi,1);
                        Interp2D(CBasis0->GetBasisKey(),
                                 CBasis1->GetBasisKey(), 
                                 &dxdxi[0], 
                                 m_base[0]->GetBasisKey(),
                                 m_base[1]->GetBasisKey(),
                                 &ngmat[0] + i*nq);
                        Vmath::Vdiv(nq,&ngmat[0]+i*nq,1,&njac[0],1,&ngmat[0]+i*nq,1);
                    }
                    
                    m_metricinfo->ResetGmat(ngmat,nq,2,coordim); 

                    // interpolate normals
                    Array<TwoD,NekDouble> newnorm(4,coordim*max(nq0,nq1));    
                    Array<TwoD, const NekDouble> normals = Xgfac->GetNormals();
                    
                    for(i = 0; i < coordim; ++i)
                    {
                        Interp1D(CBasis0->GetBasisKey(),&(normals[0])[i*Cnq0],
                                 m_base[0]->GetBasisKey(),&(newnorm[0])[i*nq0]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[1])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[1])[i*nq1]);
                        
                        Interp1D(CBasis0->GetBasisKey(),&(normals[2])[i*Cnq0],
                                 m_base[0]->GetBasisKey(),&(newnorm[2])[i*nq0]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[3])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[3])[i*nq1]);
                    }
                    
                    m_metricinfo->ResetNormals(newnorm);

                    //                    NEKERROR(ErrorUtil::ewarning,
                    //"Need to check/debug routine for deformed elements");
                }
                else // Same data can be used 
                {                   
                    // Copy Jacobian
                    Blas::Dcopy(nq, ojac, 1, njac ,1);                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    Blas::Dcopy(2*coordim*nq,&ogmat[0][0],1,ngmat.data(),1);                    
                    m_metricinfo->ResetGmat(ngmat,nq,2,coordim);                 

                    m_metricinfo->ResetNormals(Xgfac->GetNormals());

                    NEKERROR(ErrorUtil::ewarning,
                             "Need to check/debug routine for deformed elements");
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
        
        /** 
            \brief Calculate the inner product of inarray with respect to
            the basis B=base0*base1 and put into outarray:
            
            \f$ \begin{array}{rcl} I_{pq} = (\phi_q \phi_q, u) & = &
            \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1} \phi_p(\xi_{0,i})
            \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j})
            J_{i,j}\\ & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
            \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}
            J_{i,j} \end{array} \f$
            
            where
            
            \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
            
            which can be implemented as
            
            \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} = 
            {\bf B_1 U}  \f$
            \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
            {\bf B_0 F}  \f$
        **/
        void QuadExp::IProductWRTBase(const Array<OneD, const NekDouble> &base0, 
                                      const Array<OneD, const NekDouble> &base1, 
                                      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray, 
                                      const int coll_check)
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
            Array<TwoD, const NekDouble> gmat = m_metricinfo->GetGmat();

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
                    ASSERTL1(dir >= 0 &&dir < 3,"input dir is out of range");
                }
                break;
            }   
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp1,tmp3,1);
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),tmp2,outarray,1);
            Vmath::Vadd(m_ncoeffs, tmp3, 1, outarray, 1, outarray, 1);                 
        }
        
    
        /** \brief  Inner product of \a inarray over region with respect to the 
            expansion basis (this)->_Base[0] and return in \a outarray 
            
            Wrapper call to QuadExp::IProduct_WRT_B
        
            Input:\n
            
            - \a inarray: array of function evaluated at the physical
            collocation points
            
            Output:\n
            
            - \a outarray: array of inner product with respect to each
            basis over region
            
        */    
        void QuadExp::IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
                                      Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                            inarray,outarray,1);
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
                                              DetShapeType(),*this);
                DNekScalMatSharedPtr& matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
                NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
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
                    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(),&(m_geom->UpdatePhys(2))[0],
                             m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_2[0]);
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
                    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(1))[0],
                             m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_1[0]);
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
                    Interp2D(CBasis0->GetBasisKey(), CBasis1->GetBasisKey(), &(m_geom->UpdatePhys(0))[0],
                             m_base[0]->GetBasisKey(),m_base[1]->GetBasisKey(),&coords_0[0]);
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
        
        DNekScalMatSharedPtr QuadExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");
            
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
                        StdRegions::StdMatrixKey masskey(StdRegions::eMass,DetShapeType(),
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
                                           mkey.GetShapeType(), *this);  
                        MatrixKey deriv1key(StdRegions::eWeakDeriv1,
                                            mkey.GetShapeType(), *this);

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
                                           mkey.GetShapeType(), *this);  
                        MatrixKey lap01key(StdRegions::eLaplacian01,
                                           mkey.GetShapeType(), *this);
                        MatrixKey lap11key(StdRegions::eLaplacian11,
                                           mkey.GetShapeType(), *this);  

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
                    NekDouble factor = mkey.GetScaleFactor();
                    MatrixKey masskey(StdRegions::eMass,
                                      mkey.GetShapeType(), *this);    
                    DNekScalMat &MassMat = *(this->m_matrixManager[masskey]);
                    MatrixKey lapkey(StdRegions::eLaplacian,
                                     mkey.GetShapeType(), *this);
                    DNekScalMat &LapMat = *(this->m_matrixManager[lapkey]);

                    int rows = LapMat.GetRows();
                    int cols = LapMat.GetColumns();

                    DNekMatSharedPtr helm = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);

                    NekDouble one = 1.0;
                    (*helm) = LapMat + factor*MassMat;
                    
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,helm);            
                }
                break;
            case StdRegions::eUnifiedDGHelmholtz:
            case StdRegions::eUnifiedDGLamToU:
            case StdRegions::eUnifiedDGLamToQ0:
            case StdRegions::eUnifiedDGLamToQ1:
            case StdRegions::eUnifiedDGHelmBndSys:
                {
                    NekDouble one    = 1.0;
                    
                    DNekMatSharedPtr mat = GenMatrix(*mkey.GetStdMatKey());
                    returnval = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,mat);
                }
            break;
            case StdRegions::eInvUnifiedDGHelmholtz:
                {
                    NekDouble one = 1.0;

                    StdRegions::StdMatrixKey hkey(StdRegions::eUnifiedDGHelmholtz,
                                                  DetShapeType(),*this,
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


        SegExpSharedPtr QuadExp::GetEdgeExp(int edge)
        {
            SegExpSharedPtr returnval; 
            int dir = (edge == 0|| edge == 2)? 0:1;
            
            SpatialDomains::SegGeomSharedPtr edg = m_geom->GetEdge(edge);
            
            returnval = MemoryManager<SegExp>::AllocateSharedPtr(m_base[dir]->GetBasisKey(),edg);

            return returnval; 
        }


        void QuadExp::GetEdgePhysVals(const int edge, const Array<OneD, const NekDouble> &inarray, Array<OneD,NekDouble> &outarray)
        {
            int nquad0 = m_base[0]->GetNumPoints();
            int nquad1 = m_base[1]->GetNumPoints();

            Array<OneD,const NekDouble> e_tmp;

            StdRegions::EdgeOrientation edgedir = GetCartesianEorient(edge);
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
                    Vmath::Vcopy(nquad0,e_tmp = inarray+nquad0*(nquad1-1),1,
                                 outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad0,e_tmp = inarray+(nquad0*nquad1-1),-1,
                                 outarray,1);

                }
                break; 
            case 3:
                if(edgedir == StdRegions::eForwards)
                {
                    Vmath::Vcopy(nquad1,inarray,nquad0,outarray,1);
                }
                else
                {
                    Vmath::Vcopy(nquad1,e_tmp = inarray + nquad0*(nquad1-1),
                                 -nquad0,outarray,1);
                }
                break; 
            default:
                ASSERTL0(false,"edge value (< 3) is out of range");
                break;
            }

        }

        void QuadExp::AddNormBoundaryInt(const int dir,
                                         Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray,
                                         bool InArrayIsTrace) 
       {

            int i,e,cnt;
            int order_e,nquad_e;
            SpatialDomains::GeomType Gtype = m_metricinfo->GetGtype();
            Array<TwoD, const NekDouble> normals = m_metricinfo->GetNormals();
            
            StdRegions::StdExpMap vmap;

            Array<OneD,SegExpSharedPtr>   EdgeExp(4);
            
            // Set up edge segment expansions
            for(i = 0; i < 4; ++i)
            {
                EdgeExp[i] = GetEdgeExp(i);
            }

            
            cnt = 0;
            for(e = 0; e < 4; ++e)
            {
                order_e = EdgeExp[e]->GetNcoeffs();
                nquad_e = EdgeExp[e]->GetNumPoints(0);
                

                // Fill modes into expansion
                if(InArrayIsTrace == true)
                {
                    for(i = 0; i < order_e; ++i)
                    {
                        EdgeExp[e]->SetCoeff(i,inarray[i+cnt]);
                    }
                    cnt += order_e;
                }
                else
                {
                    MapTo(order_e,EdgeExp[e]->GetBasisType(0),e,
                          GetCartesianEorient(e),vmap);
                    for(i = 0; i < order_e; ++i)
                    {
                        EdgeExp[e]->SetCoeff(i,inarray[vmap[i]]);
                    }
                }

                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                     EdgeExp[e]->UpdatePhys());

                if(Gtype == SpatialDomains::eDeformed)
                {
                    if(GetCartesianEorient(e) == StdRegions::eForwards)
                    {
                        Vmath::Vmul(nquad_e,&(normals[e][dir*nquad_e]),1,
                                    &(EdgeExp[e]->GetPhys())[0],1,
                                    &(EdgeExp[e]->UpdatePhys())[0],1);
                    }
                    else
                    {
                        //use reverse ordering of normals to be
                        //consistent with edge orientation
                        Vmath::Vmul(nquad_e,&(normals[e][dir*nquad_e])+nquad_e-1,-1,
                                    &(EdgeExp[e]->GetPhys())[0],1,
                                    &(EdgeExp[e]->UpdatePhys())[0],1);

                    }
                }
                else
                {
                    Vmath::Smul(nquad_e,normals[e][dir],
                                EdgeExp[e]->GetPhys(),1,
                                EdgeExp[e]->UpdatePhys(),1);

                }
                AddEdgeBoundaryInt(e,EdgeExp[e],outarray);
            }

       }

        void QuadExp:: AddEdgeBoundaryInt(const int edge, 
                                          SegExpSharedPtr &EdgeExp,
                                          Array <OneD,NekDouble > &outarray)
            
        {
            
            int i;
            int order_e = EdgeExp->GetNcoeffs();                    
            StdRegions::StdExpMap vmap;

            MapTo(order_e,EdgeExp->GetBasisType(0),edge, 
                  GetCartesianEorient(edge),vmap);
            
            EdgeExp->IProductWRTBase(EdgeExp->GetPhys(),EdgeExp->UpdateCoeffs());
            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[vmap[i]] += EdgeExp->GetCoeff(i);
            }
        }



        // Boundary terms associated with elemental Helmholtz matrix operations
        void QuadExp::AddUDGHelmholtzBoundaryTerms(const NekDouble tau, 
                                                   const Array<OneD,
                                                   const NekDouble> &inarray,
                                                   Array<OneD,NekDouble> &outarray)
        {
            int i,e;
            int nbndry  = NumBndryCoeffs();
            int nquad0  = GetNumPoints(0);
            int nquad1  = GetNumPoints(1);
            int coordim = m_geom->GetCoordim();

            Array<OneD,NekDouble> inval(max(nquad0,nquad1));
            SpatialDomains::GeomType Gtype = m_metricinfo->GetGtype();
            Array<TwoD,const NekDouble> normals = m_metricinfo->GetNormals();

            StdRegions::StdExpMap emap;
            StdRegions::EdgeOrientation edgedir;

            Array<OneD,NekDouble> in_phys(nquad0*nquad1);

            Array<OneD,SegExpSharedPtr>        EdgeExp(4);
            Array<OneD,Array<OneD,NekDouble> > deriv(3);


            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");

            
            // Set up edge segment expansions
            for(i = 0; i < 4; ++i)
            {
                EdgeExp[i] = GetEdgeExp(i);
            }

            //  Get physical solution. 
            BwdTrans(inarray,in_phys);

            // Calculate derivative if needed for matrix terms.
            deriv[0] = Array<OneD,NekDouble>(nquad0*nquad1);
            deriv[1] = Array<OneD,NekDouble>(nquad0*nquad1);
            if(m_geom->GetCoordDim() == 2)
            {
                PhysDeriv(in_phys,deriv[0],deriv[1]);
            }
            else
            {
                deriv[2] = Array<OneD,NekDouble>(nquad0*nquad1);
                PhysDeriv(in_phys,deriv[0],deriv[1],deriv[2]);
            }

            // Loop over edges
            for(e = 0; e < 4; ++e)
            {

                GetEdgePhysVals(e,in_phys,EdgeExp[e]->UpdatePhys());
                AddUDGHelmholtzEdgeTerms(tau,e,EdgeExp[e],outarray);

#if 1
                //=============================================================
                // Add -D^T M^{-1}G operation =-<n phi_i, n.d(in_phys)/dx]>
                //term which arise in matrix formulations but not rhs
                int nquad_e = EdgeExp[e]->GetNumPoints(0);                    
                int order_e = EdgeExp[e]->GetNcoeffs();                    

                edgedir = GetCartesianEorient(e);
                MapTo(order_e,EdgeExp[e]->GetBasisType(0),e, 
                      edgedir,emap);
            
                
                Vmath::Zero(nquad_e,&(EdgeExp[e]->UpdatePhys())[0],1);
                
                if(Gtype == SpatialDomains::eDeformed)
                {
                    for(i = 0; i < coordim; ++i)
                    {
                        GetEdgePhysVals(e,deriv[i],inval);
                        if(edgedir == StdRegions::eForwards)
                        {
                            Vmath::Vvtvp(nquad_e,&inval[0],1,
                                         &(normals[e][i*nquad_e]),1,
                                         &(EdgeExp[e]->UpdatePhys())[0],1,
                                         &(EdgeExp[e]->UpdatePhys())[0],1);
                        }
                        else
                        {
                            Vmath::Vvtvp(nquad_e,&inval[0],1,
                                         &(normals[e][i*nquad_e])+nquad_e-1,-1,
                                         &(EdgeExp[e]->UpdatePhys())[0],1,
                                         &(EdgeExp[e]->UpdatePhys())[0],1);
                        }
                        
                    }
                }
                else
                {
                    for(i = 0; i < coordim; ++i)
                    {
                        GetEdgePhysVals(e,deriv[i],inval);
                        Vmath::Svtvp(nquad_e,normals[e][i],&inval[0],1,
                                     &(EdgeExp[e]->UpdatePhys())[0],1,
                                     &(EdgeExp[e]->UpdatePhys())[0],1);
                    }
                }
                
                // Fill edge and take inner product
                EdgeExp[e]->IProductWRTBase(EdgeExp[e]->GetPhys(),EdgeExp[e]->UpdateCoeffs());
                
                // Put data in out array
                for(i = 0; i < order_e; ++i)
                {
                    outarray[emap[i]] -= emap.GetSign(i)*EdgeExp[e]->GetCoeff(i);
                }                    
#endif
            }
            //================================================================
        }


        // Boundary terms associated with elemental Helmholtz matrix
        // operations from the trace space
        void QuadExp::AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                             const Array<OneD,const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             bool InputDataIsCartesianOrient)
        {
            int i;
            Array<OneD,SegExpSharedPtr>  EdgeExp(4);
            Array<OneD,NekDouble> tmp(inarray);
           
            // Set up edge segment expansions
            for(i = 0; i < 4; ++i)
            {
                EdgeExp[i] = GetEdgeExp(i);
            }

            if(InputDataIsCartesianOrient == true)
            {
                int nedge, cnt = 0;
                for(i = 0; i < 4; ++i)
                {
                    nedge = EdgeExp[i]->GetNcoeffs(); 
                    Vmath::Vcopy(nedge, &inarray[0]+cnt, 1,
                                 &(EdgeExp[i]->UpdateCoeffs())[0],1);
                    EdgeExp[i]->SetCoeffsToOrientation(GetCartesianEorient(i));
                    Vmath::Vcopy(nedge, &(EdgeExp[i]->UpdateCoeffs())[0], 1,
                                 &tmp[0]+cnt,1);
                    cnt += nedge;
                }

            }
            AddUDGHelmholtzTraceTerms(tau,tmp,EdgeExp,outarray);

        }
                                                

        // This method assumes that data in EdgeExp is ordered
        // according to elemental cartesian direction.
        void QuadExp::AddUDGHelmholtzTraceTerms(const NekDouble tau, 
                                                const Array<OneD, const NekDouble> &inarray,
                                                Array<OneD,SegExpSharedPtr> &EdgeExp,
                                                Array<OneD,NekDouble> &outarray)
        {

            ASSERTL0(&inarray[0] != &outarray[0],"Input and output arrays use the same memory");


            int e,cnt;
            int order_e;
            Array<OneD, const NekDouble> tmp;

            cnt = 0;
            for(e = 0; e < 4; ++e)
            {
                
                order_e = EdgeExp[e]->GetNcoeffs();                    
                Vmath::Vcopy(order_e,tmp =inarray+cnt,1,EdgeExp[e]->UpdateCoeffs(),1);
                EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),EdgeExp[e]->UpdatePhys());
                AddUDGHelmholtzEdgeTerms(tau,e,EdgeExp[e],outarray);
                cnt += order_e;
            }
        }
                                              


        // Calculate edge contribution assuming edgeExp values are 
        // defined according to the local orientation of the input file. 

        void QuadExp::AddUDGHelmholtzEdgeTerms(const NekDouble tau, 
                                               const int edge,
                                               SegExpSharedPtr &EdgeExp, 
                                               Array <OneD,NekDouble > &outarray)
        {            
            int i,j,k,n;
            int nquad_e = EdgeExp->GetNumPoints(0);                    
            int order_e = EdgeExp->GetNcoeffs();                    
            int nquad0  = GetNumPoints(0);
            int nquad1  = GetNumPoints(1);
            int order0  = m_base[0]->GetNumModes();
            int order1  = m_base[1]->GetNumModes();
            Array<OneD,NekDouble> inval (max(nquad0,nquad1));
            Array<OneD,NekDouble> inval1(max(nquad0,nquad1));
            Array<OneD,NekDouble> outcoeff(max(order0,order1));

            Array<OneD,NekDouble> tmp;

            SpatialDomains::GeomType Gtype = m_metricinfo->GetGtype();
            int coordim = m_geom->GetCoordim();

            Array<TwoD, const NekDouble> normals = m_metricinfo->GetNormals();
            Array<TwoD, const NekDouble> gmat    = m_metricinfo->GetGmat();
            Array<OneD, const NekDouble> jac     = EdgeExp->GetMetricInfo()->GetJac();

            MatrixKey    imasskey(StdRegions::eInvMass, DetShapeType(),*this);
            DNekScalMat  &invMass = *m_matrixManager[imasskey];
            
            StdRegions::StdExpMap emap;
            StdRegions::EdgeOrientation edgedir = GetCartesianEorient(edge);

            MapTo(order_e,EdgeExp->GetBasisType(0),edge, edgedir,emap);

#if 1
            //================================================================
            // Add F = \tau <phi_i,in_phys>
            // Fill edge and take inner product
            EdgeExp->IProductWRTBase(EdgeExp->GetPhys(),
                                     EdgeExp->UpdateCoeffs());
            // add data to out array
            for(i = 0; i < order_e; ++i)
            {
                outarray[emap[i]] += emap.GetSign(i)*tau*EdgeExp->GetCoeff(i);
            }
            //================================================================
#endif
            
#if 1 
            //================================================================
            //Add -E^T M^{-1}D_i^e = -< d\phi_i/dx.n, in_phys[j]> 
            //
            // since \phi_{ij}(r,s) = \phi^a_i(r)\phi^a_j(s) (will drop ^a in what follows):
            // <d\phi_{ij}/dn, u>|_r = <d\phi_i/dr, (rx.nx+ry.ny+rz.nz).u> \phi_j(s)
            //                       + <\phi_i, (sx.nx+sy.sy+sz.nz).u> d\phi_j(s)/ds
            //
            // So edge inner product can be expressed as two inner
            // products multiplied by factors, Note \phi_j(s) will
            // be one for a interior-boundary expansion. An
            // analogous expression exists for the |_s edges
            
            if(edge == 0|| edge == 2) // edges aligned with r
            {
                int st = (edge == 0)? 0: nquad0*(nquad1-1);
                
                // assemble jac*(rx.nx + ry.ny + rz.nz)*in_phys
                if(Gtype == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nquad0,&gmat[0][st],1,
                                &normals[edge][0],1,&inval[0],1);
                    Vmath::Vvtvp(nquad0,&gmat[2][st],1,
                                 &normals[edge][nquad0],1,
                                 &inval[0],1,&inval[0],1);
                    if(coordim == 3)
                    {
                        Vmath::Vvtvp(nquad0,&gmat[4][st],1,
                                     &normals[edge][2*nquad0],1,
                                     &inval[0],1,&inval[0],1);
                        
                    }

                    Vmath::Vmul(nquad0,jac,1,inval,1,inval1,1);
                    
                    if(edgedir == StdRegions::eForwards)
                    {
                        Vmath::Vmul(nquad0,EdgeExp->GetPhys(),1,inval1,1,inval,1);
                    }
                    else
                    {
                        Vmath::Vmul(nquad0,EdgeExp->GetPhys(),1, 
                                    tmp = inval1 + (nquad0-1),-1,inval,1);
                    }
                }
                else
                {
                    NekDouble fac = 0.0;
                    for(i = 0; i < coordim; ++i)
                    {
                        fac += jac[0]*gmat[2*i][0]*normals[edge][i];
                    }
                    Vmath::Smul(nquad0,fac,EdgeExp->GetPhys(),1,inval,1);
                }


                // Innerproduct wrt deriv base (which is only defined
                // in standard region [-1,1])
                EdgeExp->StdSegExp::IProductWRTDerivBase(0,inval,
                                                         EdgeExp->UpdateCoeffs());
                // add data to out array
                for(i = 0; i < order_e; ++i)
                {
                    outarray[emap[i]] -= emap.GetSign(i)*EdgeExp->GetCoeff(i);
                }
                    
                // assemble (sx.nx + sy.ny + sz.nz)*in_phys
                if(Gtype == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nquad0,&gmat[1][st],1,
                                &normals[edge][0],1,&inval[0],1);

                    Vmath::Vvtvp(nquad0,&gmat[3][st],1,
                                 &normals[edge][0]+nquad0,1,
                                 &inval[0],1,&inval[0],1);

                    if(coordim == 3)
                    {
                        Vmath::Vvtvp(nquad0,&gmat[5][st],1,
                                     &normals[edge][0]+2*nquad0,1,
                                     &inval[0],1,&inval[0],1);
                        
                    }
                    
                    // reverse direction of normal information
                    if(edgedir == StdRegions::eForwards)
                    {
                        Vmath::Vcopy(nquad0,inval,1,inval1,1);
                    }
                    else
                    {
                        Vmath::Vcopy(nquad0,tmp = inval+(nquad0-1),-1,inval1,1);
                    }

                    Vmath::Vmul(nquad0,EdgeExp->GetPhys(),1,inval1,1,inval,1);
                }
                else
                {
                    NekDouble fac = 0.0;
                    for(i = 0; i < coordim; ++i)
                    {
                        fac += gmat[2*i+1][0]*normals[edge][i];
                    }
                    Vmath::Smul(nquad0,fac,EdgeExp->GetPhys(),1,inval,1);
                }
                
                // Innerproduct wrt deriv base 
                EdgeExp->IProductWRTBase(inval,EdgeExp->UpdateCoeffs());
                
                // add data to out array
                Array<OneD, const NekDouble> Dbase = m_base[1]->GetDbdata();
                st = (edge == 0)? 0: nquad0-1;
                
                // unpack edge coefficients 
                if(edgedir == StdRegions::eBackwards)
                {
                    EdgeExp->ReverseCoeffsAndSign(EdgeExp->GetCoeffs(),outcoeff);
                }
                else
                {
                    Vmath::Vcopy(order0,EdgeExp->GetCoeffs(),1,outcoeff,1);
                }

                for(i = 0; i < order0; ++i)
                {
                    for(j = 0; j < order1; ++j)
                    {
                        outarray[j*order0+i] -= Dbase[j*nquad0+st]*outcoeff[i];
                    }
                }
            }
            else // edges aligned with s
            {
                int st = (edge == 1)? nquad0-1: 0 ;
                    
                // assemble jac*(sx.nx + sy.ny + sz.nz)*in_phys
                if(Gtype == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nquad1,&gmat[1][st],nquad0,
                                &normals[edge][0],1,&inval[0],1);
                    
                    Vmath::Vvtvp(nquad1,&gmat[3][st],nquad0,
                                 &normals[edge][nquad1],1,
                                 &inval[0],1,&inval[0],1);

                    if(coordim == 3)
                    {
                        Vmath::Vvtvp(nquad1,&gmat[5][st],nquad0,
                                     &normals[edge][2*nquad1],1,
                                     &inval[0],1,&inval[0],1);
                        
                    }
                    Vmath::Vmul(nquad1,jac,nquad0,inval,1,inval1,1);

                    if(edgedir == StdRegions::eForwards)
                    {
                        Vmath::Vmul(nquad1,EdgeExp->GetPhys(),1,inval1,1,inval,1);
                    }
                    else
                    {
                        Vmath::Vmul(nquad1,EdgeExp->GetPhys(),1,
                                    tmp = inval1+(nquad1-1),-1,inval,1);

                    }

                }
                else
                {
                    NekDouble fac = 0.0;
                    for(i = 0; i < coordim; ++i)
                    {
                        fac += jac[0]*gmat[2*i+1][0]*normals[edge][i];
                    }
                    Vmath::Smul(nquad1,fac,EdgeExp->GetPhys(),1,inval,1);
                }
                
                // Innerproduct wrt deriv base 
                EdgeExp->StdSegExp::IProductWRTDerivBase(0,inval,
                                           EdgeExp->UpdateCoeffs());
                
                // add data to out array
                for(i = 0; i < EdgeExp->GetNcoeffs(); ++i)
                {
                    outarray[emap[i]] -= emap.GetSign(i)*EdgeExp->GetCoeff(i);
                }
                
                // assemble (rx.nx + ry.ny + rz.nz)*in_phys
                if(Gtype == SpatialDomains::eDeformed)
                {
                    Vmath::Vmul(nquad1,&gmat[0][st],nquad0,
                                &normals[edge][0],1,&inval[0],1);
                    Vmath::Vvtvp(nquad1,&gmat[2][st],nquad0,
                                 &normals[edge][nquad1],1,
                                 &inval[0],1,&inval[0],1);
                    if(coordim == 3)
                    {
                        Vmath::Vvtvp(nquad1,&gmat[5][st],nquad0,
                                     &normals[edge][2*nquad1],1,
                                     &inval[0],1,&inval[0],1);
                        
                    }
                    Vmath::Vmul(nquad1,EdgeExp->GetPhys(),1,inval,1,inval,1);
                }
                else
                {
                    NekDouble fac = 0.0;
                    for(i = 0; i < coordim; ++i)
                    {
                        fac += gmat[2*i][0]*normals[edge][i];
                    }
                    Vmath::Smul(nquad1,fac,EdgeExp->GetPhys(),1,inval,1);
                }
                
                // Innerproduct wrt deriv base 
                EdgeExp->IProductWRTBase(inval,EdgeExp->UpdateCoeffs());
                
                // add data to out array
                Array<OneD, const NekDouble> Dbase = m_base[0]->GetDbdata();
                st = (edge == 1)? nquad0-1:0;
                

                // unpack edge coefficients 
                if(edgedir == StdRegions::eBackwards)
                {
                    EdgeExp->ReverseCoeffsAndSign(EdgeExp->GetCoeffs(),outcoeff);
                }
                else
                {
                    Vmath::Vcopy(order1,EdgeExp->GetCoeffs(),1,outcoeff,1);
                }


                for(i = 0; i < order0; ++i)
                {
                    for(j = 0; j < order1; ++j)
                    {
                        outarray[j*order0+i] -= Dbase[i*nquad1+st]*outcoeff[j];
                    }
                }
            }
            //================================================================
#endif            
            
#if 1
            //===============================================================
            // Add E M^{-1} G term 
            for(n = 0; n < coordim; ++n)
            {
                //G;
                if(Gtype == SpatialDomains::eDeformed)
                {
                    if(edgedir == StdRegions::eForwards)
                    {
                        Vmath::Vmul(nquad_e,&normals[edge][n*nquad_e],1,
                                    &(EdgeExp->GetPhys())[0],1, &inval[0],1);
                    }
                    else
                    {
                        Vmath::Vmul(nquad_e,&normals[edge][n*nquad_e]+nquad_e-1,-1,
                                    &(EdgeExp->GetPhys())[0],1, &inval[0],1);

                    }
                }
                else
                {
                    Vmath::Smul(nquad_e,normals[edge][n],
                                EdgeExp->GetPhys(),1,inval,1);
                }
                
                EdgeExp->IProductWRTBase(inval,EdgeExp->UpdateCoeffs());
                
                // M^{-1} G
                for(i = 0; i < order_e; ++i)
                {
                    outcoeff[i] = 0;
                    for(j = 0; j < order_e; ++j)
                    {
                        outcoeff[i] += emap.GetSign(i)*invMass(emap[i],emap[j])*emap.GetSign(j)*EdgeExp->GetCoeff(j);
                    }
                }
                
                EdgeExp->BwdTrans(EdgeExp->GetCoeffs(),inval); 

                if(Gtype == SpatialDomains::eDeformed)
                {
                    if(edgedir == StdRegions::eForwards)
                    {
                        Vmath::Vmul(nquad_e,&normals[edge][n*nquad_e],1,
                                    &(EdgeExp->GetPhys())[0],1,&inval[0],1);
                    }
                    else
                    {
                        Vmath::Vmul(nquad_e,&normals[edge][n*nquad_e] + nquad_e-1,
                                    -1, &(EdgeExp->GetPhys())[0],1,&inval[0],1);
                    }
                }
                else
                {
                    Vmath::Smul(nquad_e,normals[edge][n],
                                EdgeExp->GetPhys(),1,inval,1);
                }
                
                EdgeExp->IProductWRTBase(inval,EdgeExp->UpdateCoeffs());
                
                // Put data in out array
                for(i = 0; i < order_e; ++i)
                {
                    outarray[emap[i]] += emap.GetSign(i)*EdgeExp->GetCoeff(i);
                }        
            }            
            //===============================================================
#endif            
        }
        

        DNekMatSharedPtr QuadExp::GenMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            
            DNekMatSharedPtr returnval;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eUnifiedDGHelmBndSys:
                {
                    int order_e, nquad_e;
                    int i,j,e,cnt;
                    int nbndry = NumDGBndryCoeffs();
                    int nquad0 = GetNumPoints(0);
                    int nquad1 = GetNumPoints(1);
                    Array<OneD,NekDouble>       work(max(nquad0,nquad1));
                    SpatialDomains::GeomType    Gtype   = m_metricinfo->GetGtype();
                    Array<TwoD,const NekDouble> normals = m_metricinfo->GetNormals();
                    int coordim = m_geom->GetCoordim();
                    Array<OneD,SegExpSharedPtr>        EdgeExp(4);

                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);
                    NekDouble lam; 

                    StdRegions::StdExpMap emap;
                    StdRegions::EdgeOrientation edgedir;
                    
                    ASSERTL0(coordim < 3,"Needs to be  set up for expansion in 3 space");

                    // Set up edge segment expansions from local geom info
                    for(i = 0; i < 4; ++i)
                    {
                        EdgeExp[i] = GetEdgeExp(i);
                    }

                    // declare matrix space
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(nbndry, nbndry);
                    DNekMat &BndMat = *returnval;
                    
                    // Matrix to map Lambda to U
                    LocalRegions::MatrixKey Umatkey(StdRegions::eUnifiedDGLamToU,DetShapeType(),*this, lambdaval,tau);
                    DNekScalMat &LamToU = *GetLocMatrix(Umatkey); 
                
                    // Matrix to map Lambda to Q0
                    LocalRegions::MatrixKey Q0matkey(StdRegions::eUnifiedDGLamToQ0, DetShapeType(),*this, lambdaval,tau);
                    DNekScalMat &LamToQ0 = *GetLocMatrix(Q0matkey); 
                    
                    // Matrix to map Lambda to Q1
                    LocalRegions::MatrixKey Q1matkey(StdRegions::eUnifiedDGLamToQ1, DetShapeType(),*this, lambdaval,tau);
                    DNekScalMat &LamToQ1 = *GetLocMatrix(Q1matkey); 

                    // Set up matrix derived from <mu, Q_lam.n - \tau (
                    // U_lam - Lam) > Not we do not need the lambda term in
                    // continuous flux since this is equal and
                    // opposite on assembly
                    for(i = 0; i < nbndry; ++i)
                    {
                        cnt = 0;
                        
                        for(e = 0; e < 4; ++e)
                        {
                            order_e = EdgeExp[e]->GetNcoeffs();  
                            nquad_e = EdgeExp[e]->GetNumPoints(0);    

                            edgedir = GetCartesianEorient(e);
                            
                             MapTo(order_e,EdgeExp[e]->GetBasisType(0),e,
                                   edgedir,emap);

                             // Q0 * n0
                             for(j = 0; j < order_e; ++j)
                             {
                                 EdgeExp[e]->SetCoeff(j,emap.GetSign(j)*LamToQ0(emap[j],i));
                             }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());
                            if(Gtype == SpatialDomains::eDeformed)
                            {
                                if(edgedir == StdRegions::eForwards)
                                {
                                    Vmath::Vmul(nquad_e,&normals[e][0],1,
                                                &EdgeExp[e]->GetPhys()[0],1,
                                                &work[0],1);
                                }
                                else
                                {
                                    Vmath::Vmul(nquad_e,&normals[e][0] +nquad_e-1,-1,
                                                &EdgeExp[e]->GetPhys()[0],1,
                                                &work[0],1);
                                }
                            }
                            else
                            {
                                Vmath::Smul(nquad_e,normals[e][0],
                                            EdgeExp[e]->GetPhys(),1,work,1);
                            }

                            // Q1 * n1
                            for(j = 0; j < order_e; ++j)
                            {
                                EdgeExp[e]->SetCoeff(j,emap.GetSign(j)*LamToQ1(emap[j],i));
                            }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());
                            if(Gtype == SpatialDomains::eDeformed)
                            {
                                
                                if(edgedir == StdRegions::eForwards)
                                {
                                    Vmath::Vvtvp(nquad_e,&normals[e][nquad_e],1,
                                                 &EdgeExp[e]->GetPhys()[0],1,
                                                 &work[0],1,&work[0],1);
                                }
                                else
                                {
                                    Vmath::Vvtvp(nquad_e,&normals[e][nquad_e]+nquad_e-1,-1,
                                                 &EdgeExp[e]->GetPhys()[0],1,
                                                 &work[0],1,&work[0],1);

                                }
                            }
                            else
                            {
                                Vmath::Svtvp(nquad_e,normals[e][1],
                                             EdgeExp[e]->GetPhys(),1,
                                             work,1,work,1);
                            }

                            
                            // - tau (ulam - lam)
                            for(j = 0; j < order_e; ++j)
                            {
                                lam = ((i >= cnt) &&( i < cnt + order_e) &&(j == i - cnt))? 1:0;
                                EdgeExp[e]->SetCoeff(j,emap.GetSign(j)*LamToU(emap[j],i) - lam);
                            }
                            
                            EdgeExp[e]->BwdTrans(EdgeExp[e]->GetCoeffs(),
                                                 EdgeExp[e]->UpdatePhys());

                            Vmath::Svtvp(nquad_e,-tau,EdgeExp[e]->GetPhys(),1,
                                         work,1,work,1);
                      
                            EdgeExp[e]->IProductWRTBase(work,EdgeExp[e]->UpdateCoeffs());
                            
                            EdgeExp[e]->SetCoeffsToOrientation(edgedir);
                            
                            for(j = 0; j < order_e; ++j)
                            {
                                BndMat(cnt+j,i) = EdgeExp[e]->GetCoeff(j);
                            }
                            
                            cnt += order_e;

                        }
                    }

                    // Could symmetrise system here if necessary

                }
                break;
            default:
                returnval = StdQuadExp::GenMatrix(mkey);
                break;
            }

            return returnval;
        }



        DNekScalBlkMatSharedPtr QuadExp::CreateStaticCondMatrix(const MatrixKey &mkey)
        {
            DNekScalBlkMatSharedPtr returnval;

            ASSERTL2(m_metricinfo->GetGtype == SpatialDomains::eNoGeomType,"Geometric information is not set up");

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
