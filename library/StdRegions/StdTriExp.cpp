///////////////////////////////////////////////////////////////////////////////
//
// File StdTriExp.cpp
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
// Description: Triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdSegExp.h>       // for StdSegExp, etc

namespace Nektar
{
    namespace StdRegions
    {
        
        StdTriExp::StdTriExp()
        {
        }


        StdTriExp::StdTriExp(
            const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb) :
            StdExpansion (LibUtilities::StdTriData::getNumberOfCoefficients(
                              Ba.GetNumModes(),
                              Bb.GetNumModes()),
                          2,Ba,Bb),
            StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes()),
                           Ba,Bb)
        {    
            ASSERTL0(Ba.GetNumModes() <= Bb.GetNumModes(), 
                     "order in 'a' direction is higher than order "
                     "in 'b' direction");
        }

        StdTriExp::StdTriExp(const StdTriExp &T):
            StdExpansion(T),
            StdExpansion2D(T)
        {
        }

        StdTriExp::~StdTriExp()
        {
        }

        //-------------------------------
        // Integration Methods
        //-------------------------------
        NekDouble StdTriExp::v_Integral(
            const Array<OneD, const NekDouble>& inarray)
        {
            int    i;
            int nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> w1_tmp(nquad1);

            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();
            Array<OneD, const NekDouble> w1 = m_base[1]->GetW();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();

            switch(m_base[1]->GetPointsType())
            {
                case LibUtilities::eGaussLobattoLegendre: // Legendre inner product 
                {
                    for(i = 0; i < nquad1; ++i)
                    {
                        w1_tmp[i] = 0.5*(1-z1[i])*w1[i];
                    }
                    break;
                }
                case LibUtilities::eGaussRadauMAlpha1Beta0: // (0,1) Jacobi Inner product 
                {
                    Vmath::Smul(nquad1, 0.5, w1, 1, w1_tmp,1);      
                    break;
                }
                default:
                {
                    ASSERTL0(false, "populate swith for this point type");
                    break;
                }
            }

            return StdExpansion2D::Integral(inarray,w0,w1_tmp);
        }     

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        /**
         * \brief Calculate the derivative of the physical points.
         *
         * \f$ \frac{\partial u}{\partial  x_1} =  \left . 
         * \frac{2.0}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
         * \right |_{\eta_2}\f$
         *
         * \f$ \frac{\partial u}{\partial  x_2} =  \left . 
         * \frac{1+\eta_1}{1-\eta_2} \frac{\partial u}{\partial d\eta_1}
         * \right |_{\eta_2}  + \left . \frac{\partial u}{\partial d\eta_2}
         * \right |_{\eta_1}  \f$
         */
        void StdTriExp::v_PhysDeriv(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& out_d0, 
                  Array<OneD,       NekDouble>& out_d1,
                  Array<OneD,       NekDouble>& out_d2)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp(nquad0*nquad1);

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // set up geometric factor: 2/(1-z1)
            for (i = 0; i < nquad1; ++i)
            {
                wsp[i] = 2.0/(1-z1[i]);
            }

            if (out_d0.num_elements() > 0)
            {
                PhysTensorDeriv(inarray, out_d0, out_d1);
                
                for (i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,wsp[i],&out_d0[0]+i*nquad0,1);
                }

                // if no d1 required do not need to calculate both deriv
                if (out_d1.num_elements() > 0)
                {
                    // set up geometric factor: (1_z0)/(1-z1)
                    for (i = 0; i < nquad0; ++i)
                    {
                        wsp[i] = 0.5*(1+z0[i]);
                    }
                    
                    for (i = 0; i < nquad1; ++i) 
                    {
                        Vmath::Vvtvp(nquad0,&wsp[0],1,&out_d0[0]+i*nquad0,
                                     1,&out_d1[0]+i*nquad0,
                                     1,&out_d1[0]+i*nquad0,1);
                    }    
                }
            }
            else if (out_d1.num_elements() > 0)
            {
                Array<OneD, NekDouble> diff0(nquad0*nquad1);
                PhysTensorDeriv(inarray, diff0, out_d1);
                
                for (i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,wsp[i],&diff0[0]+i*nquad0,1);
                }

                for (i = 0; i < nquad0; ++i)
                {
                    wsp[i] = 0.5*(1+z0[i]);
                }
                
                for (i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,&wsp[0],1,&diff0[0]+i*nquad0,
                                 1,&out_d1[0]+i*nquad0,
                                 1,&out_d1[0]+i*nquad0,1);
                } 
            }
        }
        
        void StdTriExp::v_PhysDeriv(
            const int                           dir, 
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
                case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray);   
                    break;
                }
                case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray);   
                    break;
                }
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                    break;
                }
            }             
        }

        void StdTriExp::v_StdPhysDeriv(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& out_d0,
                  Array<OneD,       NekDouble>& out_d1,
                  Array<OneD,       NekDouble>& out_d2)
        {
            StdTriExp::v_PhysDeriv(inarray, out_d0, out_d1);
        }

        void StdTriExp::v_StdPhysDeriv(
            const int dir, 
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_PhysDeriv(dir,inarray,outarray);
        }
        

        //---------------------------------------
        // Transforms
        //---------------------------------------

        /** 
         * \brief Backward tranform for triangular elements
         *
         * @note 'q' (base[1]) runs fastest in this element.
         */
        void StdTriExp::v_BwdTrans(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_BwdTrans_SumFac(inarray,outarray);
        }


        void StdTriExp::v_BwdTrans_SumFac(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()* 
                                       m_base[1]->GetNumModes());

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  inarray,outarray,wsp);
        }

        void StdTriExp::BwdTrans_SumFacKernel(
            const Array<OneD, const NekDouble>& base0, 
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray,
                  Array<OneD,       NekDouble>& wsp)
        {
            int  i;
            int  mode;
            int  nquad0  = m_base[0]->GetNumPoints();
            int  nquad1  = m_base[1]->GetNumPoints();
            int  nmodes0 = m_base[0]->GetNumModes();
            int  nmodes1 = m_base[1]->GetNumModes();
            
            ASSERTL1(wsp.num_elements() >= nquad0*nmodes1,
                     "Workspace size is not sufficient");
            ASSERTL2((m_base[1]->GetBasisType() != LibUtilities::eOrtho_B)||
                     (m_base[1]->GetBasisType() != LibUtilities::eModified_B),
                     "Basis[1] is not of general tensor type");

            for (i = mode = 0; i < nmodes0; ++i)
            {
                Blas::Dgemv('N', nquad1,nmodes1-i,1.0,base1.get()+mode*nquad1,
                            nquad1,&inarray[0]+mode,1,0.0,&wsp[0]+i*nquad1,1);
                mode += nmodes1-i;
            }
            
            // fix for modified basis by splitting top vertex mode
            if(m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                Blas::Daxpy(nquad1,inarray[1],base1.get()+nquad1,1,
                            &wsp[0]+nquad1,1);
            }
            
            Blas::Dgemm('N','T', nquad0,nquad1,nmodes0,1.0, base0.get(),nquad0,
                        &wsp[0], nquad1,0.0, &outarray[0], nquad0);
        }
        
        void StdTriExp::v_FwdTrans(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase(inarray,outarray);
            
            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);
            
            // copy inarray in case inarray == outarray
            NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
            
            out = (*matsys)*in;
        }


        void StdTriExp::v_FwdTrans_BndConstrained(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            int i,j;
            int npoints[2] = {m_base[0]->GetNumPoints(),
                              m_base[1]->GetNumPoints()};
            int nmodes[2]  = {m_base[0]->GetNumModes(),
                              m_base[1]->GetNumModes()};

            fill(outarray.get(), outarray.get()+m_ncoeffs, 0.0 );

            Array<OneD, NekDouble> physEdge[3];
            Array<OneD, NekDouble> coeffEdge[3];
            for(i = 0; i < 3; i++)
            {
                physEdge[i]  = Array<OneD, NekDouble>(npoints[i!=0]);
                coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i!=0]);
            }

            for(i = 0; i < npoints[0]; i++)
            {
                physEdge[0][i] = inarray[i];
            }

            for(i = 0; i < npoints[1]; i++)
            {
                physEdge[1][i] = inarray[npoints[0]-1+i*npoints[0]];
                physEdge[2][i] = inarray[(npoints[1]-1)*npoints[0]-i*npoints[0]];
            }

            StdSegExpSharedPtr segexp[2] = {
                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
                    m_base[0]->GetBasisKey()),
                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(
                    m_base[1]->GetBasisKey())
            };

            Array<OneD, unsigned int> mapArray;
            Array<OneD, int>          signArray;
            NekDouble sign;

            for (i = 0; i < 3; i++)
            {
                //segexp[i!=0]->v_FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);
                segexp[i!=0]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                v_GetEdgeToElementMap(i,eForwards,mapArray,signArray);
                for (j = 0; j < nmodes[i != 0]; j++)
                {
                    sign = (NekDouble) signArray[j];
                    outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                }
            }

            Array<OneD, NekDouble> tmp0(m_ncoeffs);
            Array<OneD, NekDouble> tmp1(m_ncoeffs);
                
            StdMatrixKey      masskey(eMass,DetShapeType(),*this);
            MassMatrixOp(outarray,tmp0,masskey);
            v_IProductWRTBase(inarray,tmp1);
                
            Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);
                
            // get Mass matrix inverse (only of interior DOF)
            // use block (1,1) of the static condensed system
            // note: this block alreay contains the inverse matrix
            DNekMatSharedPtr matsys = 
                (m_stdStaticCondMatrixManager[masskey])->GetBlock(1,1);

            int nBoundaryDofs = v_NumBndryCoeffs();
            int nInteriorDofs = m_ncoeffs - nBoundaryDofs; 

            Array<OneD, NekDouble> rhs   (nInteriorDofs);
            Array<OneD, NekDouble> result(nInteriorDofs);

            v_GetInteriorMap(mapArray);

            for (i = 0; i < nInteriorDofs; i++)
            {
                rhs[i] = tmp1[ mapArray[i] ];
            }

            Blas::Dgemv('N',nInteriorDofs,nInteriorDofs,
                        1.0,&(matsys->GetPtr())[0],nInteriorDofs,
                        rhs.get(),1,
                        0.0,result.get(),1);   

            for (i = 0; i < nInteriorDofs; i++)
            {
                outarray[ mapArray[i] ] = result[i];
            }
        }

        //---------------------------------------
        // Inner product functions
        //---------------------------------------
        
        /** 
         * \brief Calculate the inner product of inarray with respect to the
         * basis B=base0[p]*base1[pq] and put into outarray.
         *
         * \f$ 
         * \begin{array}{rcl}
         * I_{pq} = (\phi^A_q \phi^B_{pq}, u) &=& 
         * \sum_{i=0}^{nq_0}\sum_{j=0}^{nq_1}
         * \phi^A_p(\eta_{0,i})\phi^B_{pq}(\eta_{1,j}) w^0_i w^1_j 
         * u(\xi_{0,i} \xi_{1,j}) \\
         * & = & \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i})
         * \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) \tilde{u}_{i,j} 
         * \end{array}
         * \f$ 
         *
         * where
         *
         * \f$  \tilde{u}_{i,j} = w^0_i w^1_j u(\xi_{0,i},\xi_{1,j}) \f$
         *
         * which can be implemented as
         *
         * \f$  f_{pj} = \sum_{i=0}^{nq_0} \phi^A_p(\eta_{0,i}) 
         * \tilde{u}_{i,j} 
         * \rightarrow {\bf B_1 U}  \f$
         * \f$  I_{pq} = \sum_{j=0}^{nq_1} \phi^B_{pq}(\eta_{1,j}) f_{pj} 
         * \rightarrow {\bf B_2[p*skip] f[skip]}  \f$
         *
         * \b Recall: \f$ \eta_{1} = \frac{2(1+\xi_1)}{(1-\xi_2)}-1, \, 
         * \eta_2 = \xi_2\f$
         *
         * \b Note: For the orthgonality of this expansion to be realised the
         * 'q' ordering must run fastest in contrast to the Quad and Hex
         * ordering where 'p' index runs fastest to be consistent with the
         * quadrature ordering.
         *
         * In the triangular space the i (i.e. \f$\eta_1\f$ direction)
         * ordering still runs fastest by convention.
         */        
        void StdTriExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTBase_SumFac(inarray,outarray);
        }
        
        void StdTriExp::v_IProductWRTBase_MatOp(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);
            
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdTriExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
                            
            Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);         
            
            // multiply by integration constants 
            MultiplyByQuadratureMetric(inarray,tmp);
            
            IProductWRTBase_SumFacKernel(
                m_base[0]->GetBdata(),
                m_base[1]->GetBdata(),
                tmp,outarray,wsp);
        }
        
        void StdTriExp::IProductWRTBase_SumFacKernel(
            const Array<OneD, const NekDouble>& base0, 
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray,
                  Array<OneD,       NekDouble>& wsp)
        {
            int    i;
            int    mode;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();   
            int    nmodes0 = m_base[0]->GetNumModes();
            int    nmodes1 = m_base[1]->GetNumModes();
            
            ASSERTL1(wsp.num_elements() >= nquad1*nmodes0,
                     "Workspace size is not sufficient");
                
            Blas::Dgemm('T','N',nquad1,nmodes0,nquad0,1.0,inarray.get(),nquad0,
                        base0.get(),nquad0,0.0,wsp.get(),nquad1);
                
            // Inner product with respect to 'b' direction 
            for (mode=i=0; i < nmodes0; ++i)
            {
                Blas::Dgemv('T',nquad1,nmodes1-i,1.0, base1.get()+mode*nquad1,
                            nquad1,wsp.get()+i*nquad1,1, 0.0, 
                            outarray.get() + mode,1);
                mode += nmodes1 - i;
            }
            
            // fix for modified basis by splitting top vertex mode
            if (m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                outarray[1] += Blas::Ddot(nquad1,base1.get()+nquad1,1,
                                          wsp.get()+nquad1,1);
            }
        }   
        
        void StdTriExp::v_IProductWRTDerivBase(
            const int                           dir, 
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }
        
        void StdTriExp::v_IProductWRTDerivBase_MatOp(
            const int                           dir, 
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            int nq = GetTotPoints();
            MatrixType mtype;

            switch(dir)
            {
                case 0:
                {
                    mtype = eIProductWRTDerivBase0;
                    break;
                }
                case 1:
                {
                    mtype = eIProductWRTDerivBase1;
                    break;
                }
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                    break;
                }
            }  

            StdMatrixKey      iprodmatkey(mtype,DetShapeType(),*this);
            DNekMatSharedPtr  iprodmat = GetStdMatrix(iprodmatkey);
 
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }

        void StdTriExp::v_IProductWRTDerivBase_SumFac(
            const int                           dir, 
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            int    i;
            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nqtot   = nquad0*nquad1; 
            int    nmodes0 = m_base[0]->GetNumModes();
            int    wspsize = max(max(nqtot,m_ncoeffs),nquad1*nmodes0);

            Array<OneD, NekDouble> gfac0(2*wspsize);
            Array<OneD, NekDouble> tmp0 (gfac0+wspsize);

            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();
            
            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }

            for(i = 0; i < nquad1; ++i)  
            {
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,
                            &tmp0[0]+i*nquad0,1);
            }
                 
            MultiplyByQuadratureMetric(tmp0,tmp0);

            switch(dir)
            {
                case 0:
                {                    
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,outarray,gfac0);
                    break;
                }
                case 1:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);    
                    const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();

                    for (i = 0; i < nquad0; ++i)
                    {
                        gfac0[i] = 0.5*(1+z0[i]);
                    }        

                    for (i = 0; i < nquad1; ++i) 
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,
                                    &tmp0[0]+i*nquad0,1);
                    }       
          
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,tmp3,gfac0); 

                    MultiplyByQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 tmp0,outarray,gfac0);  
                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,
                                &outarray[0],1);      
                    break;
                }
                default:
                {
                    ASSERTL1(false, "input dir is out of range");
                    break;
                }
            }             
        }

        //---------------------------------------
        // Evaluation functions
        //---------------------------------------

        NekDouble StdTriExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& coords)
        {
            return PhysEvaluate(coords,m_phys);
        }

        NekDouble StdTriExp::v_PhysEvaluate(
            const Array<OneD, const NekDouble>& coords,
            const Array<OneD, const NekDouble>& physvals)
        {
            Array<OneD, NekDouble> coll(2);

            // set up local coordinate system 
            if (
            	//fabs(coords[0]+1.0) < NekConstants::kNekZeroTol &&
                fabs(coords[1]-1.0) < NekConstants::kNekZeroTol)
            {
                coll[0] = -1.0;
                coll[1] =  1.0;
            }
            else
            {
                coll[0] = 2*(1+coords[0])/(1-coords[1])-1.0; 
                coll[1] = coords[1]; 
            }

            return StdExpansion2D::v_PhysEvaluate(coll,physvals);
        }

        void StdTriExp::v_FillMode(
            const int mode, Array<OneD, NekDouble> &outarray)
        {
            int   i,m;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   order0 = m_base[0]->GetNumModes();
            int   order1 = m_base[1]->GetNumModes();
            int   mode0;
            Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();

            ASSERTL2(mode <= m_ncoeffs,
                     "calling argument mode is larger than "
                     "total expansion order");

            m = order1;
            for (i = 0; i < order0; ++i, m+=order1-i)
            {
                if (m > mode)
                {
                    mode0 = i;
                    break;
                }
            }

            // deal with top vertex mode in modified basis
            if (mode == 1 && 
                m_base[0]->GetBasisType() == LibUtilities::eModified_A)
            {
                Vmath::Fill(nquad0*nquad1 , 1.0, outarray, 1);
            }
            else
            {
                for (i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad0,(NekDouble *)(base0.get()+mode0*nquad0),
                                 1,&outarray[0]+i*nquad0,1);
                }
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode*nquad1),
                            1,&outarray[0]+i,nquad0,&outarray[0]+i,nquad0);
            }
        }       


        int StdTriExp::v_GetNverts() const
        {
            return 3;
        }
        
        int StdTriExp::v_GetNedges() const
        {
            return 3;
        }

        LibUtilities::ShapeType StdTriExp::v_DetShapeType() const
        {
            return LibUtilities::eTriangle;
        }
        
        int StdTriExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                     "BasisType is not a boundary interior form");
            
            return 3 + (GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
        } 
        
        int StdTriExp::v_NumDGBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_B,
                     "BasisType is not a boundary interior form");
            
            return GetBasisNumModes(0) + 2*GetBasisNumModes(1);
        } 
        
        int StdTriExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");

            if (i == 0)
            {
                return GetBasisNumModes(0);
            }
            else
            {
                return GetBasisNumModes(1);
            }
        }

        int StdTriExp::v_GetEdgeNumPoints(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 2),"edge id is out of range");
            
            if (i == 0)
            {
                return GetNumPoints(0);
            }
            else
            {
                return GetNumPoints(1); 
            }
        }

        int StdTriExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int> &nummodes, 
            int                             &modes_offset)
        {
            int nmodes = LibUtilities::StdTriData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1]);
            modes_offset += 2;
                
            return nmodes;
        }

        LibUtilities::BasisType StdTriExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");
            
            if (i == 0)
            {
                return GetBasisType(0);
            }
            else
            {
                return GetBasisType(1);
            }
        }

        void StdTriExp::v_ReadFromFile(
            std::ifstream &infile, 
            OutputFormat   format, 
            const bool     dumpVar)
        {
            if (format == eTecplot)
            {
                int  i,j;
                int  nq0,nq1;
                int  nquad0 = m_base[0]->GetNumPoints();
                int  nquad1 = m_base[1]->GetNumPoints();
                char str[256];

                if(dumpVar)
                {
                    infile.getline(str,sizeof(str));
                    infile.getline(str,sizeof(str));
                }
                infile.getline(str,sizeof(str));
                sscanf(str,"Zone, I=%d, J=%d",&nq0,&nq1);
                ASSERTL1(nq0 == nquad0,"nquad0 does not match");
                ASSERTL1(nq1 == nquad1,"nquad0 does not match");
                
                for(j = 0; j < nquad1; ++j)
                {
                    for(i = 0; i < nquad0; ++i)
                    {
                        infile.getline(str,sizeof(str));
                        sscanf(str,"%*f %*f %lf",&m_phys[0]+j*nquad0+i);
                    }
                }
            } 
            else
            {
                ASSERTL0(false, "Input routine not implemented for "
                         "requested type of output");
            }
        }

        void StdTriExp::v_WriteToFile(
            std::ofstream &outfile,
            OutputFormat   format,
            const bool     dumpVar,
            std::string    var)
        {
            if(format==eTecplot)
            {
                int  i,j;
                int  nquad0 = m_base[0]->GetNumPoints();
                int  nquad1 = m_base[1]->GetNumPoints();
                Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
                Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
                
                if(dumpVar)
                { 
                    outfile << "Variables = z1,  z2"; 
                    outfile << ", "<< var << std::endl << std::endl;
                }
                outfile << "Zone, I=" << nquad0
                        << ", J=" << nquad1 <<", F=Point" << std::endl;
                
                for(j = 0; j < nquad1; ++j)
                {
                    for(i = 0; i < nquad0; ++i)
                    {
                        outfile << 0.5*(1+z0[i])*(1.0-z1[j])-1 <<  " " << 
                            z1[j] << " " << m_phys[j*nquad0+i] << std::endl;
                    }
                }
            }
            else if(format==eGmsh)
            {   
                if(dumpVar)
                {
                    outfile<<"View.MaxRecursionLevel = 4;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
                    outfile<<"View.AdaptVisualizationGrid = 1;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
                
                outfile<<"ST("<<endl;                
                // write the coordinates of the vertices of the triangle
                outfile<<"-1.0, -1.0, 0.0,"<<endl;
                outfile<<" 1.0, -1.0, 0.0,"<<endl;
                outfile<<"-1.0,  1.0, 0.0" <<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j;
                int maxnummodes = max(m_base[0]->GetNumModes(),
                                      m_base[1]->GetNumModes());
                   
                const LibUtilities::PointsKey Pkey1Gmsh(
                    maxnummodes,LibUtilities::eGaussGaussLegendre);
                const LibUtilities::PointsKey Pkey2Gmsh(
                    maxnummodes,LibUtilities::eGaussGaussLegendre);
                const LibUtilities::BasisKey  Bkey1Gmsh(
                    m_base[0]->GetBasisType(),maxnummodes,Pkey1Gmsh);
                const LibUtilities::BasisKey  Bkey2Gmsh(
                    m_base[1]->GetBasisType(),maxnummodes,Pkey2Gmsh);
                LibUtilities::PointsType ptype = LibUtilities::eNodalTriElec;

                StdRegions::StdNodalTriExpSharedPtr EGmsh;
                EGmsh = MemoryManager<StdRegions::StdNodalTriExp>::
                    AllocateSharedPtr(Bkey1Gmsh,Bkey2Gmsh,ptype);
                
                Array<OneD,NekDouble> xi1(EGmsh->GetNcoeffs());
                Array<OneD,NekDouble> xi2(EGmsh->GetNcoeffs());
                EGmsh->GetNodalPoints(xi1,xi2);
                
                Array<OneD,NekDouble> x(EGmsh->GetNcoeffs());
                Array<OneD,NekDouble> y(EGmsh->GetNcoeffs());
                
                for(i=0;i<EGmsh->GetNcoeffs();i++)
                {
                    x[i] = 0.5*(1.0+xi1[i]);
                    y[i] = 0.5*(1.0+xi2[i]);
                }

                int cnt  = 0;
                int cnt2 = 0;
                int nDumpCoeffs = maxnummodes*maxnummodes;
                Array<TwoD, int> dumpExponentMap(nDumpCoeffs,3,0);
                Array<OneD, int> indexMap(EGmsh->GetNcoeffs(),0);
                Array<TwoD, int> exponentMap(EGmsh->GetNcoeffs(),3,0);
                for(i = 0; i < maxnummodes; i++)
                {
                    for(j = 0; j < maxnummodes; j++)
                    {
                        if(j<maxnummodes-i)
                        {
                            exponentMap[cnt][0] = j;
                            exponentMap[cnt][1] = i;
                            indexMap[cnt++]  = cnt2;
                        }

                        dumpExponentMap[cnt2][0]   = j;
                        dumpExponentMap[cnt2++][1] = i;
                    }            
                }

                NekMatrix<NekDouble> vdm(EGmsh->GetNcoeffs(),
                                         EGmsh->GetNcoeffs());
                for(i = 0 ; i < EGmsh->GetNcoeffs(); i++)
                {
                    for(j = 0 ; j < EGmsh->GetNcoeffs(); j++)
                    {
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*
                            pow(y[i],exponentMap[j][1]);
                    }
                } 

                vdm.Invert();  

                Array<OneD, NekDouble> tmp2(EGmsh->GetNcoeffs());
                EGmsh->ModalToNodal(m_coeffs,tmp2);       

                NekVector<NekDouble> in(EGmsh->GetNcoeffs(),tmp2,eWrapper);
                NekVector<NekDouble> out(EGmsh->GetNcoeffs());
                out = vdm*in;

                Array<OneD,NekDouble> dumpOut(nDumpCoeffs,0.0);
                for(i = 0 ; i < EGmsh->GetNcoeffs(); i++)
                {
                    dumpOut[ indexMap[i]  ] = out[i];
                }

                //write the coefficients
                outfile<<"{";
                for(i = 0; i < nDumpCoeffs; i++)
                {
                    outfile<<dumpOut[i];
                    if(i < nDumpCoeffs - 1)
                    {
                        outfile<<", ";
                    }
                }
                outfile<<"};"<<endl;
              
                if(dumpVar)
                {   
                    outfile<<"INTERPOLATION_SCHEME"<<endl;
                    outfile<<"{"<<endl;
                    for(i=0; i < nDumpCoeffs; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < nDumpCoeffs; j++)
                        {
                            if(i==j)
                            {
                                outfile<<"1.00";
                            }
                            else
                            {
                                outfile<<"0.00";
                            }
                            if(j < nDumpCoeffs - 1)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nDumpCoeffs - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"}"<<endl;
                        }
                    }
                    
                    outfile<<"{"<<endl;
                    for(i=0; i < nDumpCoeffs; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < 3; j++)
                        {
                            outfile<<dumpExponentMap[i][j];
                            if(j < 2)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nDumpCoeffs  - 1)
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
                ASSERTL0(false, "Output routine not implemented for "
                         "requested type of output");
            }
        }

        void StdTriExp::v_WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  cnt = 0;
            Array<OneD, NekDouble> wsp(order0*order1,0.0);

            // Put coeffs into matrix and reverse order so that p index is
            // fastest; recall q is fastest for triangles.
            for(i = 0; i < order0; ++i)
            {
                for(j = 0; j < order1-i; ++j,cnt++)
                {
                    wsp[i+j*order1] = m_coeffs[cnt];
                }
            }

            outfile <<"Coeffs = [" << " "; 

            for(j = 0; j < order1; ++j)
            {
                for(i = 0; i < order0; ++i)
                {
                    outfile << wsp[j*order0+i] <<" ";
                }
                outfile << std::endl; 
            }
            outfile << "]" ; 
        }
        
        void StdTriExp::v_GetCoords(Array<OneD, NekDouble> &coords_0, 
                                    Array<OneD, NekDouble> &coords_1,
                                    Array<OneD, NekDouble> &coords_2)
        {
            Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
            int nq0 = GetNumPoints(0);
            int nq1 = GetNumPoints(1);
            int i,j;

            for(i = 0; i < nq1; ++i)
            {
                for(j = 0; j < nq0; ++j)
                {
                    coords_0[i*nq0+j] = (1+z0[j])*(1-z1[i])/2.0 - 1.0;
                }
                Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
            }
        }

        bool StdTriExp::v_IsBoundaryInteriorExpansion()
        {
            return m_base[0]->GetBasisType() == LibUtilities::eModified_A &&
                   m_base[1]->GetBasisType() == LibUtilities::eModified_B;
        }

        int StdTriExp::v_DetCartesianDirOfEdge(const int edge)
        {
            ASSERTL2(edge >= 0 && edge <= 2, "edge id is out of range");
            
            return edge == 0 ? 0 : 1;
        }

        const LibUtilities::BasisKey StdTriExp::v_DetEdgeBasisKey(
            const int i) const
        {
            ASSERTL2(i >= 0 && i <= 2, "edge id is out of range");
            
            if (i == 0)
            {
                return GetBasis(0)->GetBasisKey();
            }
            else
            {
                switch(m_base[1]->GetBasisType())
                {
                case LibUtilities::eModified_B:
                {
                    switch(m_base[1]->GetPointsType())
                    {
                    case LibUtilities::eGaussRadauMAlpha1Beta0:
                    {
                        LibUtilities::PointsKey pkey(
                                m_base[1]->GetBasisKey().GetPointsKey().
                                GetNumPoints()+1,
                                LibUtilities::eGaussLobattoLegendre);
                        return LibUtilities::BasisKey(
                                LibUtilities::eModified_A,
                                m_base[1]->GetNumModes(),pkey);
                        break;
                    }

                    default:
                        ASSERTL0(false,"unexpected points distribution");
                        break;
                    }
                }
                default:
                    ASSERTL0(false,"Information not available to set edge key");
                    break;
                }
            }
            return LibUtilities::NullBasisKey;
        }



        //--------------------------
        // Mappings
        //--------------------------
        
        void StdTriExp::v_GetEdgeToElementMap(
            const int                  eid, 
            const Orientation      edgeOrient,
            Array<OneD, unsigned int>& maparray,
            Array<OneD,          int>& signarray)
        {
            ASSERTL0(GetEdgeBasisType(eid) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(eid) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            
            int i;
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeCoeffs = GetEdgeNcoeffs(eid);

            if(maparray.num_elements() != nEdgeCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
            }

            if(signarray.num_elements() != nEdgeCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeCoeffs,1);
            }
            else
            {
                fill(signarray.get() , signarray.get()+nEdgeCoeffs, 1);
            }
            
            switch(eid)
            {
                case 0:
                {         
                    int cnt = 0;
                    for(i = 0; i < nEdgeCoeffs; cnt+=nummodes1-i, ++i)
                    {
                        maparray[i] = cnt; 
                    }  

                    if(edgeOrient==eBackwards)
                    {
                        swap( maparray[0] , maparray[1] );

                        for(i = 3; i < nEdgeCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        } 
                    }        
                    break;
                }
                case 1:
                {
                    maparray[0] = nummodes1;
                    maparray[1] = 1;
                    for(i = 2; i < nEdgeCoeffs; i++)
                    {
                        maparray[i] = nummodes1-1+i; 
                    } 

                    if(edgeOrient==eBackwards)
                    {
                        swap( maparray[0] , maparray[1] );

                        for(i = 3; i < nEdgeCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }                         
                    break;
                }
                case 2:
                {
                    for(i = 0; i < nEdgeCoeffs; i++)
                    {
                        maparray[i] = i;
                    }   

                    if(edgeOrient==eForwards)
                    {
                        swap( maparray[0] , maparray[1] );
                        
                        for(i = 3; i < nEdgeCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }                         
                    }                      
                    break;
                }
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
                break;
            }  
        }
        
        int StdTriExp::v_GetVertexMap(const int localVertexId)
        {
            ASSERTL0(
                GetEdgeBasisType(localVertexId) == LibUtilities::eModified_A ||
                GetEdgeBasisType(localVertexId) == LibUtilities::eModified_B,
                "Mapping not defined for this type of basis");
            
            int localDOF;
            switch(localVertexId)
            {
                case 0:
                { 
                    localDOF = 0;    
                    break;
                }
                case 1:
                {   
                    localDOF = m_base[1]->GetNumModes();                 
                    break;
                }
                case 2:
                { 
                    localDOF = 1;    
                    break;
                }
                default:
                {
                    ASSERTL0(false,"eid must be between 0 and 2");
                    break;
                }
            }
            
            return localDOF;
        }

        void StdTriExp::v_GetEdgeInteriorMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int>& maparray,
            Array<OneD,          int>& signarray)
        {
            ASSERTL0(GetEdgeBasisType(eid) == LibUtilities::eModified_A||
                     GetEdgeBasisType(eid) == LibUtilities::eModified_B,
                     "Mapping not defined for this type of basis");
            int i;
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;

            if(maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }

            if(signarray.num_elements() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nEdgeIntCoeffs, 1 );
            }

            switch(eid)
            {
                case 0:
                {         
                    int cnt = 2*nummodes1 - 1;
                    for(i = 0; i < nEdgeIntCoeffs; cnt+=nummodes1-2-i, ++i)
                    {
                        maparray[i] = cnt; 
                    }   

                    if(edgeOrient==eBackwards)
                    {                            
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }       
                    break;
                }
                case 1:
                {
                    for(i = 0; i < nEdgeIntCoeffs; i++)
                    {
                        maparray[i] = nummodes1+1+i; 
                    }            

                    if(edgeOrient==eBackwards)
                    {                            
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }            
                    break;
                }
                case 2:
                {
                    for(i = 0; i < nEdgeIntCoeffs; i++)
                    {
                        maparray[i] = 2+i;
                    }  

                    if(edgeOrient==eForwards)
                    {                            
                        for(i = 1; i < nEdgeIntCoeffs; i+=2)
                        {
                            signarray[i] = -1;
                        }
                    }                      
                    break;
                }
                default:
                {
                    ASSERTL0(false,"eid must be between 0 and 2");
                    break;
                }
            }
        }

        void StdTriExp::v_GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                     GetBasisType(1) == LibUtilities::eModified_B,
                     "Expansion not of a proper type");
            
            int i,j;
            int cnt = 0;
            int nummodes0, nummodes1;
            int startvalue;
            if(outarray.num_elements()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            startvalue = 2*nummodes1;

            for(i = 0; i < nummodes0-2; i++)
            {
                for(j = 0; j < nummodes1-3-i; j++)
                {
                    outarray[cnt++]=startvalue+j;
                }
                startvalue+=nummodes1-2-i;
            }
        }

        void StdTriExp::v_GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A &&
                     GetBasisType(1) == LibUtilities::eModified_B,
                     "Expansion not of a proper type");
            int i;
            int cnt;
            int nummodes0, nummodes1;
            int value;

            if (outarray.num_elements()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            value = 2*nummodes1-1;
            for(i = 0; i < value; i++)
            {
                outarray[i]=i;
            }
            cnt = value;

            for(i = 0; i < nummodes0-2; i++)
            {
                outarray[cnt++]=value;
                value += nummodes1-2-i;
            }
        }


        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdTriExp::v_GenMatrix(const StdMatrixKey &mkey)
        {

            MatrixType mtype   = mkey.GetMatrixType();
            
            DNekMatSharedPtr Mat; 
            
            switch(mtype)
            {
            default:
                {
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
                }
                break;
            }
            
            return Mat;
        }
        
        DNekMatSharedPtr StdTriExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }

        
        //---------------------------------------
        // Operator evaluation functions
        //---------------------------------------
        
        void StdTriExp::v_MassMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }
        
        void StdTriExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdTriExp::v_LaplacianMatrixOp_MatFree(inarray,outarray,mkey);
        }
        
        void StdTriExp::v_LaplacianMatrixOp(
            const int                           k1,
            const int                           k2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::LaplacianMatrixOp_MatFree(
                k1,k2,inarray,outarray,mkey);
        }
        
        void StdTriExp::v_WeakDerivMatrixOp(
            const int                           i,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
        }
        
        void StdTriExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdTriExp::v_HelmholtzMatrixOp_MatFree(inarray,outarray,mkey);
        }
        
        void StdTriExp::v_LaplacianMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            if(mkey.GetNVarCoeff() == 0)
            {
                // This implementation is only valid when there are no
                // coefficients associated to the Laplacian operator
                int    i;
                int    nquad0 = m_base[0]->GetNumPoints();
                int    nquad1 = m_base[1]->GetNumPoints();
                int    nquadmax = max(nquad0,nquad1);
                int    nqtot = nquad0*nquad1; 

                Array<OneD,NekDouble> physValues(3*nqtot+m_ncoeffs+nquadmax);
                Array<OneD,NekDouble> dPhysValuesdx(physValues+nqtot);
                Array<OneD,NekDouble> dPhysValuesdy(physValues+2*nqtot);
                Array<OneD,NekDouble> tmp(physValues+3*nqtot);
                Array<OneD,NekDouble> gfac0(physValues+3*nqtot+m_ncoeffs);

                BwdTrans_SumFac(inarray,physValues);

                // Laplacian matrix operation
                PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);
                // multiply with metric terms of collapsed coordinate system
                const Array<OneD,const NekDouble>& z0 = m_base[0]->GetZ();
                const Array<OneD,const NekDouble>& z1 = m_base[1]->GetZ();

                for(i = 0; i < nquad0; ++i)
                {
                    gfac0[i] = 0.5*(1+z0[i]);
                }        
            
                for(i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,&gfac0[0],1,
                                 dPhysValuesdy.get()+i*nquad0,1,
                                 dPhysValuesdx.get()+i*nquad0,1,
                                 dPhysValuesdx.get()+i*nquad0,1);
                } 

                for(i = 0; i < nquad1; ++i)
                {
                    gfac0[i] = 2.0/(1-z1[i]);
                }

                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,gfac0[i],dPhysValuesdx.get()+i*nquad0,1);
                }
             
                MultiplyByQuadratureMetric(dPhysValuesdx,dPhysValuesdx);
                MultiplyByQuadratureMetric(dPhysValuesdy,dPhysValuesdy);
                
                IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                             m_base[1]->GetBdata(),
                                             dPhysValuesdx,outarray,physValues);
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetDbdata(),
                                             dPhysValuesdy,tmp,physValues);  
                Vmath::Vadd(m_ncoeffs,tmp.get(),1,outarray.get(),1,
                            outarray.get(),1);          
            }    
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
                    inarray,outarray,mkey);
            }    
        }       


        void StdTriExp::v_SVVLaplacianFilter(Array<OneD, NekDouble> &array,
                                             const StdMatrixKey &mkey)
        {
            int qa = m_base[0]->GetNumPoints();
            int qb = m_base[1]->GetNumPoints();
            int nmodes_a = m_base[0]->GetNumModes();
            int nmodes_b = m_base[1]->GetNumModes();

            // Declare orthogonal basis. 
            LibUtilities::PointsKey pa(qa,m_base[0]->GetPointsType());
            LibUtilities::PointsKey pb(qb,m_base[1]->GetPointsType());
            
            LibUtilities::BasisKey Ba(LibUtilities::eOrtho_A,nmodes_a,pa);
            LibUtilities::BasisKey Bb(LibUtilities::eOrtho_B,nmodes_b,pb);
            StdTriExp OrthoExp(Ba,Bb);
            
            Array<OneD, NekDouble> orthocoeffs(OrthoExp.GetNcoeffs());
            int j,k;
            
            int cnt;
            int cuttoff = (int) (mkey.GetConstFactor(StdRegions::eFactorSVVCutoffRatio)*nmodes_a);
            NekDouble  SvvDiffCoeff = mkey.GetConstFactor(StdRegions::eFactorSVVDiffCoeff);
            
            // project onto physical space.
            OrthoExp.FwdTrans(array,orthocoeffs);
            
            // apply SVV filter. 
            for(cnt = j = 0; j < nmodes_a; ++j)
            {
                for(k = 0; k < nmodes_b-j; ++k)
                {
                    if(j + k >= cuttoff)
                    {
                        orthocoeffs[cnt] *= (1.0+SvvDiffCoeff*exp(-(j+k-nmodes_a)*(j+k-nmodes_a)/((NekDouble)((j+k-cuttoff+1)*(j+k-cuttoff+1)))));
                    }
                    cnt++;
                }
            }
            
            // backward transform to physical space
            OrthoExp.BwdTrans(orthocoeffs,array);
        }
        
        void StdTriExp::v_HelmholtzMatrixOp_MatFree(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            int    nquadmax = max(nquad0,nquad1);
            NekDouble lambda = mkey.GetConstFactor(eFactorLambda);

            Array<OneD,NekDouble> physValues(3*nqtot+m_ncoeffs+nquadmax);
            Array<OneD,NekDouble> dPhysValuesdx(physValues+nqtot);
            Array<OneD,NekDouble> dPhysValuesdy(physValues+2*nqtot);
            Array<OneD,NekDouble> tmp(physValues+3*nqtot);
            Array<OneD,NekDouble> gfac0(physValues+3*nqtot+m_ncoeffs);

            BwdTrans_SumFac(inarray,physValues);

            // mass matrix operation
            IProductWRTBase_SumFac(physValues,tmp);

            // Laplacian matrix operation
            PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);
            // multiply with metric terms of collapsed coordinate system
            const Array<OneD,const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD,const NekDouble>& z1 = m_base[1]->GetZ();

            for(i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }        
            
            for(i = 0; i < nquad1; ++i) 
            {
                Vmath::Vvtvp(nquad0,&gfac0[0],1,
                             dPhysValuesdy.get()+i*nquad0,1,
                             dPhysValuesdx.get()+i*nquad0,1,
                             dPhysValuesdx.get()+i*nquad0,1);
            } 

            for(i = 0; i < nquad1; ++i)
            {
                gfac0[i] = 2.0/(1-z1[i]);
            }

            for(i = 0; i < nquad1; ++i)  
            {
                Blas::Dscal(nquad0,gfac0[i],dPhysValuesdx.get()+i*nquad0,1);
            }
             
            MultiplyByQuadratureMetric(dPhysValuesdx,dPhysValuesdx);
            MultiplyByQuadratureMetric(dPhysValuesdy,dPhysValuesdy);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                         m_base[1]->GetBdata(),
                                         dPhysValuesdx,outarray,physValues);
            Blas::Daxpy(m_ncoeffs, lambda, tmp.get(), 1, outarray.get(), 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetDbdata(),
                                         dPhysValuesdy,tmp,physValues);  
            Vmath::Vadd(m_ncoeffs,tmp.get(),1,outarray.get(),1,
                        outarray.get(),1);                  
        }

        void StdTriExp::v_GeneralMatrixOp_MatOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            DNekMatSharedPtr mat = m_stdMatrixManager[mkey];
            
            if(inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);
                
                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, tmp.get(), 1, 0.0, outarray.get(), 1); 
            }
            else
            {
                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1); 
            }
        }

        //---------------------------------------
        // Private helper functions
        //---------------------------------------
        
        void StdTriExp::MultiplyByQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            int    i; 
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
                
            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,inarray.get()+i*nquad0,1,
                            w0.get(),1, outarray.get()+i*nquad0,1);
            }
                
            switch(m_base[1]->GetPointsType())
            {
                // Legendre inner product 
                case LibUtilities::eGaussLobattoLegendre: 
                    for(i = 0; i < nquad1; ++i)
                    {
                        Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i],
                                    outarray.get()+i*nquad0,1);
                    }
                    break;

                // (1,0) Jacobi Inner product 
                case LibUtilities::eGaussRadauMAlpha1Beta0:
                    for(i = 0; i < nquad1; ++i)
                    {
                        Blas::Dscal(nquad0,0.5*w1[i],outarray.get()+i*nquad0,1);
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Unsupported quadrature points type.");
                    break;
            }
        }
    }//end namespace
}//end namespace
