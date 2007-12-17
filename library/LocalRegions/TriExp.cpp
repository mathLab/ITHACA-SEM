///////////////////////////////////////////////////////////////////////////////
//
// File TriExp.cpp  
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

#include <LocalRegions/LocalRegions.hpp>
#include <stdio.h>

#include <LocalRegions/TriExp.h>

namespace Nektar
{
    namespace LocalRegions 
    {
        TriExp::TriExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb,
                       const SpatialDomains::TriGeomSharedPtr &geom):    
            StdRegions::StdTriExp(Ba,Bb),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom = geom;
            
            for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&TriExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&TriExp::CreateStaticCondMatrix, this, _1));
            }
            
            GenMetricInfo();
        }
        
        TriExp::TriExp(const LibUtilities::BasisKey &Ba,
                       const LibUtilities::BasisKey &Bb):
            StdRegions::StdTriExp(Ba,Bb),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {

           for(int i = 0; i < StdRegions::SIZE_MatrixType; ++i)
            {
                m_matrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                          StdRegions::eNoShapeType,*this),
                                                boost::bind(&TriExp::CreateMatrix, this, _1));
                m_staticCondMatrixManager.RegisterCreator(MatrixKey((StdRegions::MatrixType) i,
                                                                    StdRegions::eNoShapeType,*this),
                                                          boost::bind(&TriExp::CreateStaticCondMatrix, this, _1));
            }

            // Set up unit geometric factors. 
            m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                AllocateSharedPtr(); 
            int coordim = 2;
            Array<OneD,NekDouble> ndata = Array<OneD,NekDouble>(coordim*coordim,0.0); 
            ndata[0] = ndata[3] = 1.0;
            m_metricinfo->ResetGmat(ndata,1,2,coordim);
            m_metricinfo->ResetJac(1,ndata);
        }
        
        TriExp::TriExp(const TriExp &T):
            StdRegions::StdTriExp(T),
            m_matrixManager(std::string("StdExp")),
            m_staticCondMatrixManager(std::string("StdExpStdCondMat"))
        {
            m_geom = T.m_geom;
            m_metricinfo = T.m_metricinfo;
        }        
        
        TriExp::~TriExp()
        {
        }
    
    
        void TriExp::GenMetricInfo()
        {
            SpatialDomains::GeomFactorsSharedPtr Xgfac;
            
            Xgfac = m_geom->GetGeomFactors();
            
            if(Xgfac->GetGtype() != SpatialDomains::eDeformed)
            {
                m_metricinfo = Xgfac;
            }
            else
            {
                int nq = GetTotPoints();   
                int coordim = m_geom->GetCoordim();
                int expdim = 2;
                SpatialDomains::GeomType gtype = SpatialDomains::eDeformed;

                LibUtilities::BasisSharedPtr CBasis0;
                LibUtilities::BasisSharedPtr CBasis1;
                CBasis0 = m_geom->GetBasis(0,0); // assumes all goembasis are same
                CBasis1 = m_geom->GetBasis(0,1);
                int Cnq0 = CBasis0->GetNumPoints();
                int Cnq1 = CBasis1->GetNumPoints();
     
                ConstArray<OneD,NekDouble> ojac = Xgfac->GetJac();   
                ConstArray<TwoD,NekDouble> ogmat = Xgfac->GetGmat();
                Array<OneD,NekDouble> njac(nq);
                Array<OneD,NekDouble> ngmat(2*coordim*nq);
                
                m_metricinfo = MemoryManager<SpatialDomains::GeomFactors>::
                    AllocateSharedPtr(gtype,expdim,coordim); 
                
                //basis are different distributions
                if(!(m_base[0]->GetBasisKey().SamePoints(CBasis0->GetBasisKey()))||
                   !(m_base[1]->GetBasisKey().SamePoints(CBasis1->GetBasisKey())))
                {
                    int i;   
                    int nq0 = m_base[0]->GetNumPoints();
                    int nq1 = m_base[1]->GetNumPoints();

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
                    Array<TwoD,NekDouble> newnorm = Array<TwoD,NekDouble>(4,coordim*max(nq0,nq1));    
                    ConstArray<TwoD,NekDouble> normals = Xgfac->GetNormals();
                    
                    for(i = 0; i < coordim; ++i)
                    {
                        Interp1D(CBasis0->GetBasisKey(),&(normals[0])[i*Cnq0],
                                 m_base[0]->GetBasisKey(),&(newnorm[0])[i*nq0]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[1])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[1])[i*nq1]);
                        
                        Interp1D(CBasis1->GetBasisKey(),&(normals[2])[i*Cnq1],
                                 m_base[1]->GetBasisKey(),&(newnorm[2])[i*nq1]);
                    }
                    
                    m_metricinfo->ResetNormals(newnorm);
                    
                    NEKERROR(ErrorUtil::ewarning,
                             "Need to check/debug routine for deformed elements");
                }
                else // Same data can be used 
                {                   
                    // Copy Jacobian
                    Blas::Dcopy(nq,&ojac[0],1,&njac[0],1);                    
                    m_metricinfo->ResetJac(nq,njac);

                    // interpolate Geometric data
                    ngmat = Array<OneD,NekDouble>(2*nq*coordim); 
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
            \xi_1 d \xi_2 \f$ where \f$inarray[i,j] = u(\xi_{1i},\xi_{2j})
            \f$ and \f$ J[i,j] \f$ is the Jacobian evaluated at the
            quadrature point.
        */
        NekDouble TriExp::Integral(const ConstArray<OneD,NekDouble> &inarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            NekDouble ival;
            Array<OneD,NekDouble> tmp   = Array<OneD,NekDouble>(nquad0*nquad1);
            
            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1,&jac[0],1,(NekDouble*)&inarray[0],1,
                            &tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1,(NekDouble) jac[0],
                            (NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            
            // call StdQuadExp version;
            ival = StdTriExp::Integral(tmp);
            
            return ival; 
        }
        
        
        /** 
            \brief Calculate the inner product of inarray with respect to
            the basis B=base0*base1 and put into outarray:
            
            \f$ 
            \begin{array}{rcl}
            I_{pq} = (\phi_q \phi_q, u) & = & \sum_{i=0}^{nq_0} \sum_{j=0}^{nq_1}
            \phi_p(\xi_{0,i}) \phi_q(\xi_{1,j}) w^0_i w^1_j u(\xi_{0,i} \xi_{1,j}) 
            J_{i,j}\                                    \
            & = & \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i})
            \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j}  J_{i,j}
            \end{array}
            \f$ 
            
            where
            
            \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
            
            which can be implemented as
            
            \f$  f_{qi} = \sum_{j=0}^{nq_1} \phi_q(\xi_{1,j}) \tilde{u}_{i,j} = 
            {\bf B_1 U}  \f$
            \f$  I_{pq} = \sum_{i=0}^{nq_0} \phi_p(\xi_{0,i}) f_{qi} = 
            {\bf B_0 F}  \f$
        **/
        
        void TriExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
                                     const ConstArray<OneD,NekDouble> &base1, 
                                     const ConstArray<OneD,NekDouble> &inarray,
                                     Array<OneD,NekDouble> &outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            ConstArray<OneD,NekDouble> jac = m_metricinfo->GetJac();
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(nquad0*nquad1);
            
            // multiply inarray with Jacobian
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                Vmath::Vmul(nquad0*nquad1,&jac[0],1,(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            else
            {
                Vmath::Smul(nquad0*nquad1,jac[0],(NekDouble*)&inarray[0],1,&tmp[0],1);
            }
            
            StdTriExp::IProductWRTBase(base0,base1,tmp,outarray);
        }


        void TriExp::IProductWRTDerivBase(const int dir, 
                                              const ConstArray<OneD, NekDouble>& inarray, 
                                              Array<OneD, NekDouble> & outarray)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            int    numModes = m_ncoeffs;
            ConstArray<TwoD,NekDouble> gmat = m_metricinfo->GetGmat();
            
            Array<OneD, NekDouble> tmp0(nqtot);
            Array<OneD, NekDouble> tmp1(nqtot);
            Array<OneD, NekDouble> tmp2(nqtot);
            Array<OneD, NekDouble> tmp3(numModes);

            ConstArray<OneD, NekDouble> z0,z1;
            Array<OneD, NekDouble> gfac0(nqtot);
            Array<OneD, NekDouble> gfac1(nqtot);

            z0 = ExpPointsProperties(0)->GetZ();
            z1 = ExpPointsProperties(1)->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }
            for(i = 0; i < nquad0; ++i)
            {
                gfac1[i] = 0.5*(1+z0[i]);
            }

            for(i = 0; i < nquad1; ++i)  
            {
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
            }
            
            for(i = 0; i < nquad1; ++i) 
            {
                Vmath::Vmul(nquad0,&gfac1[0],1,&tmp0[0]+i*nquad0,1,&tmp1[0]+i*nquad0,1);
            }
            
            switch(dir)
            {
            case 0:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nqtot,&gmat[0][0],1,&tmp0[0],1,&tmp0[0],1);
                        Vmath::Vmul(nqtot,&gmat[1][0],1,&tmp1[0],1,&tmp1[0],1);
                        Vmath::Vmul(nqtot,&gmat[1][0],1,&inarray[0],1,&tmp2[0],1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot,gmat[0][0],&tmp0[0],1,&tmp0[0],1);
                        Vmath::Smul(nqtot,gmat[1][0],&tmp1[0],1,&tmp1[0],1);
                        Vmath::Smul(nqtot,gmat[1][0],&inarray[0],1,&tmp2[0],1);
                    }
                }
                break;
            case 1:
                {
                    if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vmul(nqtot,&gmat[2][0],1,&tmp0[0],1,&tmp0[0],1);
                        Vmath::Vmul(nqtot,&gmat[3][0],1,&tmp1[0],1,&tmp1[0],1);
                        Vmath::Vmul(nqtot,&gmat[3][0],1,&inarray[0],1,&tmp2[0],1);
                    }
                    else
                    {
                        Vmath::Smul(nqtot,gmat[2][0],&tmp0[0],1,&tmp0[0],1);
                        Vmath::Smul(nqtot,gmat[3][0],&tmp1[0],1,&tmp1[0],1);
                        Vmath::Smul(nqtot,gmat[3][0],&inarray[0],1,&tmp2[0],1);
                    }
                }
                break;
            default:
                {
                    ASSERTL1(dir >= 0 &&dir < 2,"input dir is out of range");
                }
                break;
            }       
            Vmath::Vadd(nqtot,&tmp0[0],1,&tmp1[0],1,&tmp1[0],1); 
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),tmp1,tmp3);
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),tmp2,outarray);
            Vmath::Vadd(numModes,&tmp3[0],1,&outarray[0],1,&outarray[0],1);                 
        }
        
        /** \brief  Inner product of \a inarray over region with respect to the 
            expansion basis (this)->_Base[0] and return in \a outarray 
            
            Wrapper call to TriExp::IProductWRTBase
            
            Input:\n
            
            - \a inarray: array of function evaluated at the physical
            collocation points
            
            Output:\n
            
            - \a outarray: array of inner product with respect to each
            basis over region
            
        */
        void TriExp::IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, 
                                     Array<OneD,NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                            inarray,outarray);
        }    
        
        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////
        
        /** 
            \brief Calculate the deritive of the physical points 
        **/
        void TriExp::PhysDeriv(const ConstArray<OneD,NekDouble> & inarray,
                               Array<OneD,NekDouble> &out_d0,
                               Array<OneD,NekDouble> &out_d1,
                               Array<OneD,NekDouble> &out_d2)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            ConstArray<TwoD,NekDouble> gmat = m_metricinfo->GetGmat();
            Array<OneD,NekDouble> Diff0 = Array<OneD,NekDouble>(nquad0*nquad1);
            Array<OneD,NekDouble> Diff1 = Array<OneD,NekDouble>(nquad0*nquad1);
            
            //         if(m_geom)
            //         {
            //         ASSERTL2(n <= m_geom->GetCoordDim(),
            //              "value of n is larger than the number of coordinates");
            //         }
            
            StdTriExp::PhysDeriv(inarray, Diff0, Diff1);
            
            if(m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
            {
                if(out_d0.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,&gmat[0][0],1,&Diff0[0],1, &out_d0[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1,&gmat[1][0],1,&Diff1[0],1, &out_d0[0], 1,
                                  &out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,&gmat[2][0],1,&Diff0[0],1, &out_d1[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1,&gmat[3][0],1,&Diff1[0],1, &out_d1[0], 1,
                                  &out_d1[0],1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Vmul  (nquad0*nquad1,&gmat[4][0],1,&Diff0[0],1, &out_d2[0], 1);
                    Vmath::Vvtvp (nquad0*nquad1,&gmat[5][0],1,&Diff1[0],1, &out_d2[0], 1,
                                  &out_d2[0],1);
                }
            }
            else // regular geometry
            {
                if(out_d0.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1,gmat[0][0],&Diff0[0],1, &out_d0[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[1][0],&Diff1[0],1, &out_d0[0], 1);
                    //Vmath::Svtvp (nquad0*nquad1,gmat[1][0],&Diff1[0],1, &out_d0[0], 1,
                    //            &out_d0[0],1);
                }
                
                if(out_d1.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1,gmat[2][0],&Diff0[0],1, &out_d1[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[3][0],&Diff1[0],1, &out_d1[0], 1);
                }
                
                if(out_d2.num_elements())
                {
                    Vmath::Smul  (nquad0*nquad1,gmat[4][0],&Diff0[0],1, &out_d2[0], 1);
                    Blas::Daxpy (nquad0*nquad1,gmat[5][0],&Diff1[0],1, &out_d2[0], 1);
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
        
        
        void TriExp::FwdTrans(const ConstArray<OneD,NekDouble> & inarray, 
                              Array<OneD,NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(GetNcoeffs(),&inarray[0],1,&m_coeffs[0],1);
            }
            else
            {
                IProductWRTBase(inarray,outarray);
                
                // get Mass matrix inverse
                MatrixKey             masskey(StdRegions::eInvMass,
                                              DetShapeType(),*this);
                DNekScalMatSharedPtr  matsys = m_matrixManager[masskey];
                
                // copy inarray in case inarray == outarray
                DNekVec in (m_ncoeffs,outarray);
                DNekVec out(m_ncoeffs,outarray,eWrapper);
                
                out = (*matsys)*in;
            }
        }
        
        void TriExp::GetCoords(Array<OneD,NekDouble> &coords_0,
                               Array<OneD,NekDouble> &coords_1,
                               Array<OneD,NekDouble> &coords_2)
        {
            LibUtilities::BasisSharedPtr CBasis0;
            LibUtilities::BasisSharedPtr CBasis1;
            Array<OneD,NekDouble>  x;
            
            ASSERTL0(m_geom, "m_geom not define");
            
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
                                &x[0],1,&coords_2[0],1);
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
                                &x[0],1,&coords_1[0],1);
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
                                &x[0],1,&coords_0[0],1);
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
        void TriExp::GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, 
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
        
        
        void TriExp::WriteToFile(FILE *outfile)
        {
            int i,j;
            Array<OneD,NekDouble> coords[3];
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            
            ASSERTL0(m_geom,"_geom not defined");
            
            int  coordim   = m_geom->GetCoordDim();
            
            coords[0] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[1] = Array<OneD,NekDouble>(nquad0*nquad1);
            coords[2] = Array<OneD,NekDouble>(nquad0*nquad1);
            
            std::fprintf(outfile,"Variables = x");
            if(coordim == 2)
            {
                GetCoords(coords[0],coords[1]);
                fprintf(outfile,", y");
            }
            else if (coordim == 3)
            {
                GetCoords(coords[0],coords[1],coords[2]);
                fprintf(outfile,", y, z");
            }
            
            fprintf(outfile,", v\n");
            
            fprintf(outfile,"Zone, I=%d, J=%d, F=Point\n",nquad0,nquad1);
            for(i = 0; i < nquad0*nquad1; ++i)
            {
                for(j = 0; j < coordim; ++j)
                {
                    fprintf(outfile,"%lf ",coords[j][i]);
                }
                fprintf(outfile,"%lf \n",m_phys[i]);
            }
        }
        
        void TriExp::WriteToFile(std::ofstream &outfile, const int dumpVar)
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
        
        DNekMatSharedPtr TriExp::GetStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
                TriExpSharedPtr tmp = MemoryManager<TriExp>::
                    AllocateSharedPtr(bkey0,bkey1);
                
                return tmp->StdTriExp::GetStdMatrix(mkey);                
            }
            else
            {
                return StdTriExp::GetStdMatrix(mkey);
            }
        }

        DNekBlkMatSharedPtr TriExp::GetStdStaticCondMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            // Need to check if matrix exists in stdMatrixManager.
            // If not then make a local expansion with standard metric info
            // and generate matrix. Otherwise direct call is OK. 
            if(!StdMatManagerAlreadyCreated(mkey))
            {
                LibUtilities::BasisKey bkey0 = m_base[0]->GetBasisKey();
                LibUtilities::BasisKey bkey1 = m_base[1]->GetBasisKey();
                TriExpSharedPtr tmp = MemoryManager<TriExp>::
                    AllocateSharedPtr(bkey0,bkey1);
                
                return tmp->StdTriExp::GetStdStaticCondMatrix(mkey);                
            }
            else
            {
                return StdTriExp::GetStdStaticCondMatrix(mkey);
            }
        }
        
        NekDouble TriExp::PhysEvaluate(const ConstArray<OneD,NekDouble> &coord)
        {
            Array<OneD,NekDouble> Lcoord = Array<OneD,NekDouble>(2);      
            
            ASSERTL0(m_geom,"m_geom not defined");
            m_geom->GetLocCoords(coord,Lcoord);
            
            return StdTriExp::PhysEvaluate(Lcoord);
        }
        
        DNekScalMatSharedPtr TriExp::CreateMatrix(const MatrixKey &mkey)
        {
            DNekScalMatSharedPtr returnval;
            
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
                        
                        DNekMatSharedPtr lap00 = GetStdMatrix(*lap00key.GetStdMatKey());
                        DNekMatSharedPtr lap01 = GetStdMatrix(*lap01key.GetStdMatKey());
                        DNekMatSharedPtr lap11 = GetStdMatrix(*lap11key.GetStdMatKey());
                        
                        NekDouble jac = (m_metricinfo->GetJac())[0];
                        ConstArray<TwoD,NekDouble> gmat = m_metricinfo->GetGmat();
                        
                        int rows = lap00->GetRows();
                        int cols = lap00->GetColumns();
                        
                        DNekMatSharedPtr lap = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols);
                        
                        (*lap) = (gmat[0][0]*gmat[0][0] + gmat[2][0]*gmat[2][0]) * (*lap00) + 
                            (gmat[0][0]*gmat[1][0] + gmat[2][0]*gmat[3][0]) * (*lap01 + Transpose(*lap01)) +
                            (gmat[1][0]*gmat[1][0] + gmat[3][0]*gmat[3][0]) * (*lap11);
                        
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
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

        DNekScalBlkMatSharedPtr TriExp::CreateStaticCondMatrix(const MatrixKey &mkey)
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
                    DNekScalMatSharedPtr mat = GetLocMatrix(mkey);
                    factor = mat->Scale();
                    goto UseStdRegionsMatrix;
                }
                break;
            UseStdRegionsMatrix:
                {
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekBlkMatSharedPtr  mat = GetStdStaticCondMatrix(*(mkey.GetStdMatKey()));                    
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
                    NekDouble            invfactor = 1.0/factor;
                    NekDouble            one = 1.0;
                    DNekScalMat &mat = *GetLocMatrix(mkey);
                    DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
                    DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
                    DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
                    DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);
                    
                    const ConstArray<OneD,int> bmap = GetBoundaryMap();
                    const ConstArray<OneD,int> imap = GetInteriorMap();
                    
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
 *    $Log: TriExp.cpp,v $
 *    Revision 1.23  2007/12/06 22:49:09  pvos
 *    2D Helmholtz solver updates
 *
 *    Revision 1.22  2007/11/08 16:54:27  pvos
 *    Updates towards 2D helmholtz solver
 *
 *    Revision 1.21  2007/08/11 23:41:22  sherwin
 *    Various updates
 *
 *    Revision 1.20  2007/07/28 05:09:33  sherwin
 *    Fixed version with updated MemoryManager
 *
 *    Revision 1.19  2007/07/20 00:45:52  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.18  2007/07/12 12:53:01  sherwin
 *    Updated to have a helmholtz matrix
 *
 *    Revision 1.17  2007/07/11 19:26:04  sherwin
 *    update for new Manager structure
 *
 *    Revision 1.16  2007/07/11 06:36:23  sherwin
 *    Updates with MatrixManager update
 *
 *    Revision 1.15  2007/07/10 17:17:26  sherwin
 *    Introduced Scaled Matrices into the MatrixManager
 *
 *    Revision 1.14  2007/06/17 19:00:45  bnelson
 *    Removed unused variables.
 *
 *    Revision 1.13  2007/06/07 15:54:19  pvos
 *    Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 *    Also made corrections to various ASSERTL2 calls
 *
 *    Revision 1.12  2007/06/06 11:29:31  pvos
 *    Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
 *
 *    Revision 1.11  2007/06/01 17:08:07  pvos
 *    Modification to make LocalRegions/Project2D run correctly (PART1)
 *
 *    Revision 1.10  2007/05/31 19:13:12  pvos
 *    Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 *    Revision 1.9  2007/05/31 11:38:17  pvos
 *    Updated QuadExp and TriExp
 *
 *    Revision 1.8  2006/12/10 18:59:46  sherwin
 *    Updates for Nodal points
 *
 *    Revision 1.7  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.6  2006/06/02 18:48:39  sherwin
 *    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 *    Revision 1.5  2006/06/01 14:15:58  sherwin
 *    Added typdef of boost wrappers and made GeoFac a boost shared pointer.
 *
 *    Revision 1.4  2006/05/30 14:00:04  sherwin
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
 *    Revision 1.17  2006/03/13 19:47:54  sherwin
 *
 *    Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
 *
 *    Revision 1.16  2006/03/13 18:20:33  sherwin
 *
 *    Fixed error in ResetGmat
 *
 *    Revision 1.15  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.14  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/

