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

namespace Nektar
{
    namespace StdRegions
    {
        
        StdTriExp::StdTriExp() // default constructor of StdExpansion is directly called. 
        {
        } //default constructor


        StdTriExp::StdTriExp(const LibUtilities::BasisKey &Ba, const LibUtilities::BasisKey &Bb):
            StdExpansion2D(StdTriData::getNumberOfCoefficients(Ba.GetNumModes(),(Ba.GetNumModes())),Ba,Bb)
        {    
            ASSERTL0(Ba.GetNumModes() <=  Bb.GetNumModes(), "order in 'a' direction is higher than order in 'b' direction");
        }

        StdTriExp::StdTriExp(const StdTriExp &T):
            StdExpansion2D(T)
        {
        }

        StdTriExp::~StdTriExp()
        {
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////
        NekDouble StdTriExp::Integral(const Array<OneD, const NekDouble>& inarray)
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
                for(i = 0; i < nquad1; ++i)
                {
                    w1_tmp[i] = 0.5*(1-z1[i])*w1[i];
                }
                break;
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (0,1) Jacobi Inner product 
                Vmath::Smul(nquad1, 0.5, w1, 1, w1_tmp,1);      
                break;
            default:
                ASSERTL0(false, "populate swith for this point type");
            }

            return StdExpansion2D::Integral(inarray,w0,w1_tmp);
        }     

        void StdTriExp::IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray, 
                                              Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetExpansionType(),*this);
            DNekMatSharedPtr& iprodmat = GetStdMatrix(iprodmatkey);            
            
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1.0, 0.0, outarray.get(), 1.0);
        }

        void StdTriExp::IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                                Array<OneD, NekDouble> &outarray)
        {
            int i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
                            
            Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);         
            
            // multiply by integration constants 
            MultiplyByQuadratureMetric(inarray,tmp);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),m_base[1]->GetBdata(),tmp,outarray,wsp);
        }
     
        void StdTriExp::IProductWRTDerivBase_SumFac(const int dir, 
                                                    const Array<OneD, const NekDouble>& inarray, 
                                                    Array<OneD, NekDouble> & outarray)
        {
            int    i;
            int    nquad0  = m_base[0]->GetNumPoints();
            int    nquad1  = m_base[1]->GetNumPoints();
            int    nqtot   = nquad0*nquad1; 
            int    wspsize = max(nqtot,m_ncoeffs);

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
                Vmath::Smul(nquad0,gfac0[i],&inarray[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
            }
                 
            MultiplyByQuadratureMetric(tmp0,tmp0);

            switch(dir)
            {
            case 0:
                {                    
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,outarray,gfac0);
                }
                break;
            case 1:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);    
                    const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();

                    for(i = 0; i < nquad0; ++i)
                    {
                        gfac0[i] = 0.5*(1+z0[i]);
                    }        

                    for(i = 0; i < nquad1; ++i) 
                    {
                        Vmath::Vmul(nquad0,&gfac0[0],1,&tmp0[0]+i*nquad0,1,&tmp0[0]+i*nquad0,1);
                    }       
          
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 tmp0,tmp3,gfac0); 

                    MultiplyByQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 tmp0,outarray,gfac0);  
                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,&outarray[0],1);      
                }
                break;
            default:
                {
                    ASSERTL1(dir >= 0 &&dir < 2,"input dir is out of range");
                }
                break;
            }             
        }

        void StdTriExp::IProductWRTDerivBase_MatOp(const int dir, 
                                                   const Array<OneD, const NekDouble>& inarray, 
                                                   Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            MatrixType mtype;

            switch(dir)
            {
            case 0:
                {
                    mtype = eIProductWRTDerivBase0;
                }
                break;
            case 1:
                {
                    mtype = eIProductWRTDerivBase1;
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }  

            StdMatrixKey      iprodmatkey(mtype,DetExpansionType(),*this);
            DNekMatSharedPtr& iprodmat = GetStdMatrix(iprodmatkey);            
 
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1.0, 0.0, outarray.get(), 1.0);            
        }

        void StdTriExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int    i,m;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            int   order0 = m_base[0]->GetNumModes();
            int   order1 = m_base[1]->GetNumModes();
            Array<OneD, const NekDouble> base0 = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1 = m_base[1]->GetBdata();
            int   mode0;

            ASSERTL2(mode >= m_ncoeffs, 
                     "calling argument mode is larger than total expansion order");

            m= order1;
            for(i = 0; i < order0; ++i, m+=order1-i)
            {
                if(m > mode)
                {
                    mode0 = i;
                    break;
                }
            }

            // deal with top vertex mode in modified basis
            if((mode == 1)&&(m_base[0]->GetBasisType() == LibUtilities::eModified_A))
            {
                Vmath::Fill(nquad0*nquad1 , 1.0, outarray, 1);
            }
            else
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Vmath::Vcopy(nquad0,(NekDouble *)(base0.get()+mode0*nquad0),
                                 1,&outarray[0]+i*nquad0,1);
                }
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode*nquad1),1,&outarray[0]+i,
                            nquad0,&outarray[0]+i,nquad0);
            }
        }       

        void StdTriExp::GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray,
                                               const StdMatrixKey &mkey)
        {
            DNekMatSharedPtr& mat = m_stdMatrixManager[mkey];
            
            if(inarray.get() == outarray.get())
            {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                Vmath::Vcopy(m_ncoeffs,inarray.get(),1,tmp.get(),1);
                
                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, tmp.get(), 1.0, 0.0, outarray.get(), 1.0); 
            }
            else
            {
                Blas::Dgemv('N', m_ncoeffs, m_ncoeffs, 1.0, mat->GetPtr().get(),
                            m_ncoeffs, inarray.get(), 1.0, 0.0, outarray.get(), 1.0); 
            }
        }
        
        void StdTriExp::MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
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
            case LibUtilities::eGaussLobattoLegendre: // Legendre inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*(1-z1[i])*w1[i], outarray.get()+i*nquad0,1);
                }
                break;
            case LibUtilities::eGaussRadauMAlpha1Beta0: // (1,0) Jacobi Inner product 
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,0.5*w1[i], outarray.get()+i*nquad0,1);      
                }
                break;
            }          
        }

        void StdTriExp::LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdMatrixKey &mkey)
        {
            if(mkey.GetNvariableLaplacianCoefficients() == 0)
            {
                // This implementation is only valid when there are no coefficients
                // associated to the Laplacian operator
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
                    Vmath::Vvtvp(nquad0,&gfac0[0],1,dPhysValuesdy.get()+i*nquad0,1,dPhysValuesdx.get()+i*nquad0,1,dPhysValuesdx.get()+i*nquad0,1);
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
                
                IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),dPhysValuesdx,outarray,physValues);
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,tmp,physValues);  
                Vmath::Vadd(m_ncoeffs,tmp.get(),1,outarray.get(),1,outarray.get(),1);          
            }    
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }    
        }       
        
        void StdTriExp::HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdMatrixKey &mkey)
        {
            cout << "hie" << endl;
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 
            int    nquadmax = max(nquad0,nquad1);
            NekDouble lambda = mkey.GetConstant(0);

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
                Vmath::Vvtvp(nquad0,&gfac0[0],1,dPhysValuesdy.get()+i*nquad0,1,dPhysValuesdx.get()+i*nquad0,1,dPhysValuesdx.get()+i*nquad0,1);
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
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),dPhysValuesdx,outarray,physValues);
            Blas::Daxpy(m_ncoeffs, lambda, tmp.get(), 1, outarray.get(), 1);

            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,tmp,physValues);  
            Vmath::Vadd(m_ncoeffs,tmp.get(),1,outarray.get(),1,outarray.get(),1);                  
        }

        //-----------------------------
        // Differentiation Methods
        //-----------------------------

        void StdTriExp::PhysDeriv(const Array<OneD, const NekDouble>& inarray, 
                                  Array<OneD, NekDouble> &out_d0, 
                                  Array<OneD, NekDouble> &out_d1,
                                  Array<OneD, NekDouble> &out_d2)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> wsp(nquad0*nquad1);

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();

            // set up geometric factor: 2/(1-z1)
            for(i = 0; i < nquad1; ++i)
            {
                wsp[i] = 2.0/(1-z1[i]);
            }

            if(out_d0.num_elements() > 0)
            {
                PhysTensorDeriv(inarray, out_d0, out_d1);
                
                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,wsp[i],&out_d0[0]+i*nquad0,1);
                }

                if(out_d1.num_elements() > 0)// if no d1 required do not need to calculate both deriv
                {
                    // set up geometric factor: (1_z0)/(1-z1)
                    for(i = 0; i < nquad0; ++i)
                    {
                        wsp[i] = 0.5*(1+z0[i]);
                    }
                    
                    for(i = 0; i < nquad1; ++i) 
                    {
                        Vmath::Vvtvp(nquad0,&wsp[0],1,&out_d0[0]+i*nquad0,1,&out_d1[0]+i*nquad0,1,
                                     &out_d1[0]+i*nquad0,1);
                    }    
                }
            }
            else if(out_d1.num_elements() > 0)
            {
                Array<OneD, NekDouble> diff0(nquad0*nquad1);
                PhysTensorDeriv(inarray, diff0, out_d1);
                
                for(i = 0; i < nquad1; ++i)  
                {
                    Blas::Dscal(nquad0,wsp[i],&diff0[0]+i*nquad0,1);
                }

                for(i = 0; i < nquad0; ++i)
                {
                    wsp[i] = 0.5*(1+z0[i]);
                }
                
                for(i = 0; i < nquad1; ++i) 
                {
                    Vmath::Vvtvp(nquad0,&wsp[0],1,&diff0[0]+i*nquad0,1,&out_d1[0]+i*nquad0,1,
                                 &out_d1[0]+i*nquad0,1);
                } 
            }
        }
        
        void StdTriExp::PhysDeriv(const int dir, 
                                  const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD, NekDouble> &outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    PhysDeriv(inarray, outarray, NullNekDouble1DArray);   
                }
                break;
            case 1:
                {
                    PhysDeriv(inarray, NullNekDouble1DArray, outarray);   
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }             
        }


        ///////////////////////////////
        // Evaluation Methods
        ///////////////////////////////

        void StdTriExp::BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                        Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()* 
                                       m_base[1]->GetNumModes());

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  inarray,outarray,wsp);
        }

        void StdTriExp::FwdTrans(const Array<OneD, const NekDouble>& inarray, 
                                 Array<OneD, NekDouble> &outarray)
        {

            IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetExpansionType(),*this);
            DNekMatSharedPtr& matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);

            out = (*matsys)*in;
        }

        void StdTriExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                                Array<OneD, NekDouble> &outarray)
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

            StdSegExpSharedPtr segexp[2] = {MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[0]->GetBasisKey()),
                                            MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[1]->GetBasisKey())};

            Array<OneD, unsigned int> mapArray;
            Array<OneD, int>          signArray;
            NekDouble sign;

            for(i = 0; i < 3; i++)
            {
                segexp[i!=0]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                GetEdgeToElementMap(i,eForwards,mapArray,signArray);
                for(j=0; j < nmodes[i!=0]; j++)
                {
                    sign = (NekDouble) signArray[j];
                    outarray[ mapArray[j] ] = sign * coeffEdge[i][j];
                }
            }

            Array<OneD, NekDouble> tmp0(m_ncoeffs);
            Array<OneD, NekDouble> tmp1(m_ncoeffs);
                
            StdMatrixKey      masskey(eMass,DetExpansionType(),*this);
            MassMatrixOp(outarray,tmp0,masskey);
            IProductWRTBase(inarray,tmp1);
                
            Vmath::Vsub(m_ncoeffs, tmp1, 1, tmp0, 1, tmp1, 1);
                
            // get Mass matrix inverse (only of interior DOF)
            // use block (1,1) of the static condensed system
            // note: this block alreay contains the inverse matrix
            DNekMatSharedPtr  matsys = (m_stdStaticCondMatrixManager[masskey])->GetBlock(1,1);

            int nBoundaryDofs = NumBndryCoeffs();
            int nInteriorDofs = m_ncoeffs - nBoundaryDofs; 

            Array<OneD, NekDouble> rhs(nInteriorDofs);
            Array<OneD, NekDouble> result(nInteriorDofs);

            GetInteriorMap(mapArray);

            for(i = 0; i < nInteriorDofs; i++)
            {
                rhs[i] = tmp1[ mapArray[i] ];
            }

            Blas::Dgemv('N',nInteriorDofs,nInteriorDofs,1.0, &(matsys->GetPtr())[0],
                        nInteriorDofs,rhs.get(),1,0.0,result.get(),1);   

            for(i = 0; i < nInteriorDofs; i++)
            {
                outarray[ mapArray[i] ] = result[i];
            }
        }

        NekDouble StdTriExp::PhysEvaluate(const Array<OneD, const NekDouble>& coords)
        {
            Array<OneD, NekDouble> coll(2);

            // set up local coordinate system 
            if((fabs(coords[0]+1.0) < NekConstants::kNekZeroTol)
               &&(fabs(coords[1]-1.0) < NekConstants::kNekZeroTol))
            {
                coll[0] = 0.0;
                coll[1] = 1.0;
            }
            else
            {
                coll[0] = 2*(1+coords[0])/(1-coords[1])-1.0; 
                coll[1] = coords[1]; 
            }

            return  StdExpansion2D::PhysEvaluate(coll); 
        }

        void StdTriExp::GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            ASSERTL1((Btype0 == LibUtilities::eModified_A)&&
                     (Btype1 == LibUtilities::eModified_B),
                     "Expansion not of a proper type");
            int i;
            int cnt;
            int nummodes0, nummodes1;
            int value;
            if(outarray.num_elements()!=NumBndryCoeffs())
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

        void StdTriExp::GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            ASSERTL1((Btype0 == LibUtilities::eModified_A)&&
                     (Btype1 == LibUtilities::eModified_B),
                     "Expansion not of a proper type");

            int i,j;
            int cnt=0;
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

        int StdTriExp::GetVertexMap(const int localVertexId)
        {
            ASSERTL0((GetEdgeBasisType(localVertexId)==LibUtilities::eModified_A)||
                     (GetEdgeBasisType(localVertexId)==LibUtilities::eModified_B),
                     "Mapping not defined for this type of basis");
            
            int localDOF;
            switch(localVertexId)
            {
            case 0:
                { 
                    localDOF = 0;    
                }
                break;
            case 1:
                {   
                    localDOF = m_base[1]->GetNumModes();                 
                }
                break;
            case 2:
                { 
                    localDOF = 1;    
                }
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
                break;
            }

            return localDOF;
        }
 
        void StdTriExp::GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                           Array<OneD, unsigned int> &maparray,
                                           Array<OneD, int> &signarray)
        {
            ASSERTL0((GetEdgeBasisType(eid)==LibUtilities::eModified_A)||
                     (GetEdgeBasisType(eid)==LibUtilities::eModified_B),
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
                }
                break;
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
                }
                break;
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
                }
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
                break;
            }  
        }

        void StdTriExp::GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int> &signarray)
        {           
            ASSERTL0((GetEdgeBasisType(eid)==LibUtilities::eModified_A)||
                     (GetEdgeBasisType(eid)==LibUtilities::eModified_B),
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
                fill( signarray.get() , signarray.get()+nEdgeCoeffs, 1 );
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
                }
                break;
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
                }
                break;
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
                }
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 2");
                break;
            }  
        }

        void StdTriExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
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
                outfile << "Zone, I=" << nquad0 <<", J=" << nquad1 <<", F=Point" << std::endl;
                
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
                    outfile<<"View.MaxRecursionLevel = 8;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }
                
                outfile<<"ST("<<endl;                
                // write the coordinates of the vertices of the triangle
                outfile<<"-1.0, -1.0, 0.0,"<<endl;
                outfile<<" 1.0, -1.0, 0.0,"<<endl;
                outfile<<"-1.0,  1.0, 0.0" <<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j,k;
                int maxnummodes = max(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
                   
                const LibUtilities::PointsKey Pkey1Gmsh(maxnummodes,LibUtilities::eGaussGaussLegendre);
                const LibUtilities::PointsKey Pkey2Gmsh(maxnummodes,LibUtilities::eGaussGaussLegendre);
                const LibUtilities::BasisKey  Bkey1Gmsh(m_base[0]->GetBasisType(),maxnummodes,Pkey1Gmsh);
                const LibUtilities::BasisKey  Bkey2Gmsh(m_base[1]->GetBasisType(),maxnummodes,Pkey2Gmsh);
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

                NekMatrix<NekDouble> vdm(EGmsh->GetNcoeffs(),EGmsh->GetNcoeffs());
                for(i = 0 ; i < EGmsh->GetNcoeffs(); i++)
                {
                    for(j = 0 ; j < EGmsh->GetNcoeffs(); j++)
                    {
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*pow(y[i],exponentMap[j][1]);
                    }
                } 

                vdm.Invert();  

                Array<OneD, NekDouble> tmp2(EGmsh->GetNcoeffs());
                EGmsh->ModalToNodal(m_coeffs,tmp2);       

                NekVector<const NekDouble> in(EGmsh->GetNcoeffs(),tmp2,eWrapper);
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
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }

        
        void StdTriExp::ReadFromFile(std::ifstream &infile, OutputFormat format, const bool dumpVar)
        {
            if(format==eTecplot)
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
                        sscanf(str,"%*lf %*lf %lf",&m_phys[0]+j*nquad0+i);
                    }
                }
            } 
            else
            {
                ASSERTL0(false, "Input routine not implemented for requested type of output");
            }
        }

        //   I/O routine
        void StdTriExp::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int  i,j;
            int  order0 = m_base[0]->GetNumModes();
            int  order1 = m_base[1]->GetNumModes();
            int  cnt = 0;
            Array<OneD, NekDouble> wsp(order0*order1,0.0);

            // put coeffs into matrix and reverse order so that p index is fastest
            // recall q is fastest for tri's

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
        
        const LibUtilities::BasisKey StdTriExp::DetEdgeBasisKey(const int i) const
        {
            ASSERTL2((i >= 0)&&(i <= 2),"edge id is out of range");
            
            if(i == 0)
            {
                return  GetBasis(0)->GetBasisKey();
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
                                LibUtilities::PointsKey pkey(m_base[1]->GetBasisKey().GetPointsKey().GetNumPoints()+1,LibUtilities::eGaussLobattoLegendre);
                                
                                return LibUtilities::BasisKey(LibUtilities::eModified_A,m_base[1]->GetNumModes(),pkey);
                            }
                            break;
                            
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
        }


        void StdTriExp::GetCoords(Array<OneD, NekDouble> &coords_0, 
                                  Array<OneD, NekDouble> &coords_1)
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

    }//end namespace
}//end namespace


/** 
 * $Log: StdTriExp.cpp,v $
 * Revision 1.54  2009/04/22 22:30:48  sherwin
 * Added ReadFromFile method to read back in .dat file
 *
 * Revision 1.53  2009/04/20 16:11:47  sherwin
 * Mods to handle output and optimise DG work
 *
 * Revision 1.52  2009/01/21 16:58:39  pvos
 * Added additional geometric factors to improve efficiency
 *
 * Revision 1.51  2008/12/18 14:11:35  pvos
 * NekConstants Update
 *
 * Revision 1.50  2008/12/16 11:31:52  pvos
 * Performance updates
 *
 * Revision 1.49  2008/11/24 10:31:14  pvos
 * Changed name from _PartitionedOp to _MatFree
 *
 * Revision 1.48  2008/11/19 16:02:47  pvos
 * Added functionality for variable Laplacian coeffcients
 *
 * Revision 1.47  2008/11/05 16:08:15  pvos
 * Added elemental optimisation functionality
 *
 * Revision 1.46  2008/09/23 18:19:26  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.45  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.44  2008/09/08 08:04:18  pvos
 * removed NEKTAR_USING_DIRECT_BLAS_CALLS
 *
 * Revision 1.43  2008/08/28 15:03:54  pvos
 * small efficiency updates
 *
 * Revision 1.42  2008/07/31 21:20:56  sherwin
 * Updates to make DG advection run with tris
 *
 * Revision 1.41  2008/07/31 11:10:15  sherwin
 * Updates for handling EdgeBasisKey for use with DG advection. Depracated GetEdgeBasis and added DetEdgeBasisKey
 *
 * Revision 1.40  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.39  2008/07/04 10:18:41  pvos
 * Some updates
 *
 * Revision 1.38  2008/07/02 14:08:56  pvos
 * Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 * Revision 1.37  2008/06/05 20:13:12  ehan
 * Removed undefined function GetAlpha() in the ASSERTL2().
 *
 * Revision 1.36  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.35  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.34  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.33  2008/04/22 05:22:15  bnelson
 * Speed enhancements.
 *
 * Revision 1.32  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.31  2008/04/03 16:12:11  pvos
 * updates for NEKTAR_USING_DIRECT_BLAS_CALLS
 *
 * Revision 1.30  2008/04/02 22:18:10  pvos
 * Update for 2D local to global mapping
 *
 * Revision 1.29  2008/03/12 15:25:09  pvos
 * Clean up of the code
 *
 * Revision 1.28  2008/01/23 09:09:46  sherwin
 * Updates for Hybrized DG
 *
 * Revision 1.27  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.26  2007/12/06 22:44:47  pvos
 * 2D Helmholtz solver updates
 *
 * Revision 1.25  2007/11/20 10:16:23  pvos
 * IProductWRTDerivBase routine fix
 *
 * Revision 1.24  2007/11/08 16:55:14  pvos
 * Updates towards 2D helmholtz solver
 *
 * Revision 1.23  2007/10/15 20:40:24  ehan
 * Make changes of column major matrix
 *
 * Revision 1.22  2007/07/20 02:16:55  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.21  2007/07/11 13:35:18  kirby
 * *** empty log message ***
 *
 * Revision 1.20  2007/07/10 21:05:18  kirby
 * even more fixes
 *
 * Revision 1.19  2007/07/09 15:19:15  sherwin
 * Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
 *
 * Revision 1.18  2007/06/07 15:54:19  pvos
 * Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 * Also made corrections to various ASSERTL2 calls
 *
 * Revision 1.17  2007/05/31 19:13:12  pvos
 * Updated NodalTriExp + LocalRegions/Project2D + some other modifications
 *
 * Revision 1.16  2007/05/22 02:01:50  bnelson
 * Changed Array::size to Array::num_elements.
 *
 * Fixed some compiler errors in assertions.
 *
 * Revision 1.15  2007/05/15 05:18:24  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.14  2007/04/10 14:00:46  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.13  2007/04/06 08:44:43  sherwin
 * Update to make 2D regions work at StdRegions level
 *
 * Revision 1.12  2007/04/05 15:20:11  sherwin
 * Updated 2D stuff to comply with SharedArray philosophy
 *
 * Revision 1.11  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.10  2007/01/17 16:36:58  pvos
 * updating doxygen documentation
 *
 * Revision 1.9  2007/01/17 16:05:41  pvos
 * updated doxygen documentation
 *
 * Revision 1.8  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.7  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.6  2006/07/02 17:16:19  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.5  2006/06/13 18:05:02  sherwin
 * Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 * Revision 1.4  2006/06/06 15:25:21  jfrazier
 * Removed unreferenced variables and replaced ASSERTL0(false, ....) with
 * NEKERROR.
 *
 * Revision 1.3  2006/06/02 18:48:40  sherwin
 * Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
 *
 * Revision 1.2  2006/06/01 14:13:37  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.50  2006/04/25 20:23:34  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.49  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.48  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.47  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.46  2006/03/06 12:39:59  sherwin
 *
 * Added NekConstants class for all constants in this library
 *
 * Revision 1.45  2006/03/03 23:04:54  sherwin
 *
 * Corrected Mistake in StdBasis.cpp to do with eModified_B
 *
 * Revision 1.44  2006/03/01 08:25:05  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.43  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/ 


