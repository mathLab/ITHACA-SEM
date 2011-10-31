///////////////////////////////////////////////////////////////////////////////
//
// File StdQuadExp.cpp
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
// Description: Quadrilateral routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
    namespace StdRegions
    {


        StdQuadExp::StdQuadExp()
        {
        }

        StdQuadExp::StdQuadExp(const LibUtilities::BasisKey &Ba, 
                               const LibUtilities::BasisKey &Bb):
            StdExpansion  (Ba.GetNumModes()*Bb.GetNumModes(),2,Ba,Bb),
            StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb)
        { 
        }

        StdQuadExp::StdQuadExp(const StdQuadExp &T):
            StdExpansion(T),
            StdExpansion2D(T)
        {
        }

        StdQuadExp::~StdQuadExp()
        {
        }

        //////////////////////////////
        // Integration Methods
        //////////////////////////////

        NekDouble StdQuadExp::Integral(const Array<OneD, const NekDouble>& inarray)
        {
            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();
            Array<OneD, const NekDouble> w1 = m_base[1]->GetW();
            
            return StdExpansion2D::Integral(inarray,w0,w1);
        }

        void StdQuadExp::IProductWRTBase_SumFac(const Array<OneD, const NekDouble>& inarray, 
                                                Array<OneD, NekDouble> &outarray)
        {
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    order0 = m_base[0]->GetNumModes();
                            
            Array<OneD,NekDouble> tmp(nquad0*nquad1+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nquad0*nquad1);         
            
            // multiply by integration constants 
            MultiplyByQuadratureMetric(inarray,tmp);
            
            IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                         m_base[1]->GetBdata(),
                                         tmp,outarray,wsp,true,true);
        }
        
        void StdQuadExp::IProductWRTBase_MatOp(const Array<OneD, const NekDouble>& inarray, 
                                               Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      iprodmatkey(eIProductWRTBase,DetExpansionType(),*this);
            DNekMatSharedPtr& iprodmat = GetStdMatrix(iprodmatkey);            
            
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);
        }
        
        void StdQuadExp::IProductWRTDerivBase_SumFac(const int dir, 
                                                     const Array<OneD, const NekDouble>& inarray, 
                                                     Array<OneD, NekDouble> &outarray)
        {   
            ASSERTL0((dir==0)||(dir==1),"input dir is out of range");
 
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int     nqtot = nquad0*nquad1;
            int    order0 = m_base[0]->GetNumModes();
                            
            Array<OneD,NekDouble> tmp(nqtot+nquad1*order0);
            Array<OneD,NekDouble> wsp(tmp+nqtot);         
            
            // multiply by integration constants 
            MultiplyByQuadratureMetric(inarray,tmp);
            
            if(dir) // dir == 1
            {
                IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                             m_base[1]->GetDbdata(),
                                             tmp,outarray,wsp,true,false);
            }
            else    // dir == 0
            {
                IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                             m_base[1]->GetBdata(),
                                             tmp,outarray,wsp,false,true);
            } 
        }  
        
        void StdQuadExp::IProductWRTDerivBase_MatOp(const int dir, 
                                                    const Array<OneD, const NekDouble>& inarray, 
                                                    Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0((dir==0)||(dir==1),"input dir is out of range");

            int nq = GetTotPoints();
            MatrixType mtype;

            if(dir) // dir == 1
            {
                mtype = eIProductWRTDerivBase1;
            }
            else    // dir == 0
            {
                mtype = eIProductWRTDerivBase0;
            } 
            
            StdMatrixKey      iprodmatkey(mtype,DetExpansionType(),*this);
            DNekMatSharedPtr& iprodmat = GetStdMatrix(iprodmatkey);            
 
            Blas::Dgemv('N',m_ncoeffs,nq,1.0,iprodmat->GetPtr().get(),
                        m_ncoeffs, inarray.get(), 1, 0.0, outarray.get(), 1);            
        }
        
        void StdQuadExp::FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int   nquad0 = m_base[0]->GetNumPoints();
            int   nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, const NekDouble> base0  = m_base[0]->GetBdata();
            Array<OneD, const NekDouble> base1  = m_base[1]->GetBdata();
            int   btmp0 = m_base[0]->GetNumModes();
            int   mode0 = mode%btmp0;
            int   mode1 = mode/btmp0;


            ASSERTL2(mode1 == (int)floor((1.0*mode)/btmp0),
                     "Integer Truncation not Equiv to Floor");

            ASSERTL2(m_ncoeffs <= mode, 
                     "calling argument mode is larger than total expansion order");

            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vcopy(nquad0,(NekDouble *)(base0.get() + mode0*nquad0),
                             1, &outarray[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,(NekDouble *)(base1.get() + mode1*nquad1),1,
                            &outarray[0]+i,nquad0,&outarray[0]+i,nquad0);
            }
        }

        DNekMatSharedPtr StdQuadExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            int      i;
            int      order0    = GetBasisNumModes(0);
            int      order1    = GetBasisNumModes(1);
            MatrixType mtype   = mkey.GetMatrixType();
            
            DNekMatSharedPtr Mat; 

            switch(mtype)
            {
            case eMass:
				{
					Mat = StdExpansion::CreateGeneralMatrix(mkey);
					// For Fourier basis set the imaginary component of mean mode
					// to have a unit diagonal component in mass matrix 
					if(m_base[0]->GetBasisType() == LibUtilities::eFourier)
					{
						for(i = 0; i < order1; ++i)
						{
							(*Mat)(order0*i+1,i*order0+1) = 1.0;
						}
					}
                
					if(m_base[1]->GetBasisType() == LibUtilities::eFourier)
					{
						for(i = 0; i < order0; ++i)
						{
							(*Mat)(order0+i ,order0+i) = 1.0;
						}
					}
				}
				break;
			case eFwdTrans:
				{
					Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);
                    StdMatrixKey iprodkey(eIProductWRTBase,DetExpansionType(),*this);
                    DNekMat &Iprod = *GetStdMatrix(iprodkey);
                    StdMatrixKey imasskey(eInvMass,DetExpansionType(),*this);
                    DNekMat &Imass = *GetStdMatrix(imasskey);
                    
                    (*Mat) = Imass*Iprod;
				}
				break;
			default:
                {                
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
				}
				break;
            }
	
            return Mat;
        }
        
        void StdQuadExp::GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD,NekDouble> &outarray,
                                               const StdMatrixKey &mkey)
        {
            DNekMatSharedPtr& mat = m_stdMatrixManager[mkey];
            
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
        
        void StdQuadExp::MultiplyByQuadratureMetric(const Array<OneD, const NekDouble>& inarray,
                                                    Array<OneD, NekDouble> &outarray)
        {         
            int i; 
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
                
            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0,1,
                            w0.get(),1,outarray.get()+i*nquad0,1);
            }
                
            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,outarray.get()+i,nquad0,w1.get(),1,
                            outarray.get()+i,nquad0);
            }                
        }

        void StdQuadExp::LaplacianMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdMatrixKey &mkey)
        {
            if(mkey.GetNvariableLaplacianCoefficients() == 0)
            {
                // This implementation is only valid when there are no coefficients
                // associated to the Laplacian operator
                int       nquad0  = m_base[0]->GetNumPoints();
                int       nquad1  = m_base[1]->GetNumPoints();
                int       nqtot   = nquad0*nquad1; 
                int       nmodes0 = m_base[0]->GetNumModes();
                int       nmodes1 = m_base[1]->GetNumModes();
                int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);
                
                const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
                const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
                const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
                const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();
                
                // Allocate temporary storage
                Array<OneD,NekDouble> wsp0(3*wspsize);
                Array<OneD,NekDouble> wsp1(wsp0+wspsize);
                Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
                
                if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
                {  
                    // LAPLACIAN MATRIX OPERATION
                    // wsp0 = u       = B   * u_hat 
                    // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                    // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                    BwdTrans_SumFacKernel(base0,base1,inarray,wsp0,wsp1,true,true);
                    StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
                }
                else
                {
                    StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);                    
                }
                
                // wsp1 = k = wsp1 * w0 * w1
                // wsp2 = l = wsp2 * w0 * w1
                MultiplyByQuadratureMetric(wsp1,wsp1);
                MultiplyByQuadratureMetric(wsp2,wsp2);
                
                // outarray = m = (D_xi1 * B)^T * k 
                // wsp1     = n = (D_xi2 * B)^T * l 
                IProductWRTBase_SumFacKernel(dbase0,base1,wsp1,outarray,wsp0,false,true);
                IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,    wsp0,true,false);
                
                // outarray = outarray + wsp1
                //          = L * u_hat
                Vmath::Vadd(m_ncoeffs,wsp1.get(),1,outarray.get(),1,outarray.get(),1);     
            }
            else
            {
                StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(inarray,outarray,mkey);
            }
        }

        void StdQuadExp::HelmholtzMatrixOp_MatFree(const Array<OneD, const NekDouble> &inarray,
                                                         Array<OneD,NekDouble> &outarray,
                                                         const StdMatrixKey &mkey)
        { 
            int       nquad0  = m_base[0]->GetNumPoints();
            int       nquad1  = m_base[1]->GetNumPoints();
            int       nqtot   = nquad0*nquad1; 
            int       nmodes0 = m_base[0]->GetNumModes();
            int       nmodes1 = m_base[1]->GetNumModes();
            int       wspsize = max(max(max(nqtot,m_ncoeffs),nquad1*nmodes0),nquad0*nmodes1);
            NekDouble lambda  = mkey.GetConstant(0);
                                
            const Array<OneD, const NekDouble>& base0  = m_base[0]->GetBdata();
            const Array<OneD, const NekDouble>& base1  = m_base[1]->GetBdata();
            const Array<OneD, const NekDouble>& dbase0 = m_base[0]->GetDbdata();
            const Array<OneD, const NekDouble>& dbase1 = m_base[1]->GetDbdata();

            // Allocate temporary storage
            Array<OneD,NekDouble> wsp0(4*wspsize);
            Array<OneD,NekDouble> wsp1(wsp0+wspsize);
            Array<OneD,NekDouble> wsp2(wsp0+2*wspsize);
            Array<OneD,NekDouble> wsp3(wsp0+3*wspsize);

            if(!(m_base[0]->Collocation() && m_base[1]->Collocation()))
            {  
                // MASS MATRIX OPERATION
                // The following is being calculated:
                // wsp0     = B   * u_hat = u
                // wsp1     = W   * wsp0
                // outarray = B^T * wsp1  = B^T * W * B * u_hat = M * u_hat 
                BwdTrans_SumFacKernel       (base0,base1,inarray,wsp0,    wsp1,true,true);
                MultiplyByQuadratureMetric  (wsp0,wsp2);
                IProductWRTBase_SumFacKernel(base0,base1,wsp2,   outarray,wsp1,true,true);
                
                // LAPLACIAN MATRIX OPERATION
                // wsp1 = du_dxi1 = D_xi1 * wsp0 = D_xi1 * u
                // wsp2 = du_dxi2 = D_xi2 * wsp0 = D_xi2 * u
                StdExpansion2D::PhysTensorDeriv(wsp0,wsp1,wsp2);
            }
            else
            {
                // specialised implementation for the classical spectral element method
                StdExpansion2D::PhysTensorDeriv(inarray,wsp1,wsp2);
                MultiplyByQuadratureMetric(inarray,outarray);
            }
            
            // wsp1 = k = wsp1 * w0 * w1
            // wsp2 = l = wsp2 * w0 * w1
            MultiplyByQuadratureMetric(wsp1,wsp1);
            MultiplyByQuadratureMetric(wsp2,wsp2);
            
            // wsp1 = m = (D_xi1 * B)^T * k 
            // wsp0 = n = (D_xi2 * B)^T * l 
            IProductWRTBase_SumFacKernel(dbase0,base1,wsp1,wsp0,wsp3,false,true);
            IProductWRTBase_SumFacKernel(base0,dbase1,wsp2,wsp1,wsp3,true,false);

            // outarray = lambda * outarray + (wsp0 + wsp1)
            //          = (lambda * M + L ) * u_hat
            Vmath::Vstvpp(m_ncoeffs,lambda,&outarray[0],1,&wsp1[0],1,&wsp0[0],1,&outarray[0],1);
        }
        
        ///////////////////////////////
        /// Differentiation Methods
        ///////////////////////////////

        void StdQuadExp::PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                   Array<OneD, NekDouble> &out_d0,
                                   Array<OneD, NekDouble> &out_d1,
                                   Array<OneD, NekDouble> &out_d2)
        {
            PhysTensorDeriv(inarray, out_d0, out_d1);
        }
        
        void StdQuadExp::PhysDeriv(const int dir, 
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

        //------------------------------
        // Evaluation Methods
        //-----------------------------
        
        void StdQuadExp::BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> wsp(m_base[0]->GetNumPoints()* 
                                       m_base[1]->GetNumModes());

            BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                  m_base[1]->GetBdata(),
                                  inarray,outarray,wsp,true,true);
        }

        void StdQuadExp::v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD, NekDouble> &outarray)
        {
            if((m_base[0]->Collocation())&&(m_base[1]->Collocation()))
            {
                Vmath::Vcopy(m_ncoeffs, inarray, 1, outarray, 1);
            }
            else
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
        }

        void StdQuadExp::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
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
                for(i = 0; i < 4; i++)
                {
                    physEdge[i]  = Array<OneD, NekDouble>(npoints[i%2]);
                    coeffEdge[i] = Array<OneD, NekDouble>(nmodes[i%2]);
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

                StdSegExpSharedPtr segexp[2] = {MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[0]->GetBasisKey()),
                                                MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(m_base[1]->GetBasisKey())};

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;
                NekDouble sign;

                for(i = 0; i < 4; i++)
                {
                    //segexp[i%2]->v_FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);
                    segexp[i%2]->FwdTrans_BndConstrained(physEdge[i],coeffEdge[i]);

                    GetEdgeToElementMap(i,eForwards,mapArray,signArray);
                    for(j=0; j < nmodes[i%2]; j++)
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

        }

        void StdQuadExp::GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            int i;
            int cnt=0;
            int nummodes0, nummodes1;
            int value1, value2;
            if(outarray.num_elements()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }

            nummodes0 = m_base[0]->GetNumModes();
            nummodes1 = m_base[1]->GetNumModes();

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
                value1 = nummodes0;
                break;
            case LibUtilities::eModified_A:
                value1 = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < value1; i++)
            {
                outarray[i]=i;
            }
            cnt=value1;

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
                value2 = value1+nummodes0-1;
                break;
            case LibUtilities::eModified_A:
                value2 = value1+1;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                outarray[cnt++]=value1+i*nummodes0;
                outarray[cnt++]=value2+i*nummodes0;
            }


            if(Btype1 == LibUtilities::eGLL_Lagrange)
            {
                for(i = nummodes0*(nummodes1-1);i < GetNcoeffs(); i++)
                { 
                    outarray[cnt++] = i;
                }
            }
        }

        void StdQuadExp::GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
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

            const LibUtilities::BasisType Btype0 = GetBasisType(0);
            const LibUtilities::BasisType Btype1 = GetBasisType(1);

            switch(Btype1)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue = nummodes0;
                break;
            case LibUtilities::eModified_A:
                startvalue = 2*nummodes0;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            switch(Btype0)
            {
            case LibUtilities::eGLL_Lagrange:
                startvalue++;
                break;
            case LibUtilities::eModified_A:
                startvalue+=2;
                break;
            default:
                ASSERTL0(0,"Mapping array is not defined for this expansion");
                break;
            }

            for(i = 0; i < nummodes1-2; i++)
            {
                for(j = 0; j < nummodes0-2; j++)
                {
                    outarray[cnt++]=startvalue+j;
                }
                startvalue+=nummodes0;
            }
        }

        int StdQuadExp::GetVertexMap(const int localVertexId)
        {
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
                    if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                    {
                        localDOF = m_base[0]->GetNumModes()-1;
                    }
                    else
                    {
                        localDOF = 1;
                    }
                }
                break;
            case 2:
                {   
                    if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                    {
                        localDOF = m_base[0]->GetNumModes()*m_base[1]->GetNumModes()-1;
                    }
                    else
                    {
                        localDOF = m_base[0]->GetNumModes()+1;
                    }                    
                }
                break;
            case 3:
                { 
                    if(m_base[0]->GetBasisType()==LibUtilities::eGLL_Lagrange)
                    {
                        localDOF = m_base[0]->GetNumModes() * (m_base[1]->GetNumModes()-1);
                    }
                    else
                    {
                        localDOF = m_base[0]->GetNumModes();
                    }
                }
                break;
            default:
                ASSERTL0(false,"eid must be between 0 and 3");
                break;
            }

            return localDOF;
        }
 
        void StdQuadExp::GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                            Array<OneD, unsigned int> &maparray,
                                            Array<OneD, int> &signarray)
        {
            int i;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;
            const LibUtilities::BasisType bType = GetEdgeBasisType(eid);

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

            if(bType == LibUtilities::eModified_A)
            {
                switch(eid)
                {
                case 0:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = i+2;
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
                            maparray[i] = (i+2)*nummodes0 + 1;
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
                            maparray[i] = nummodes0+i+2;
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
                case 3:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0;
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
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }  
            }
            else if(bType == LibUtilities::eGLL_Lagrange)
            {
                switch(eid)
                {
                case 0:
                    {                  
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = i+1;
                        }
                    }
                    break;
                case 1:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0 - 1;
                        }                     
                    }
                    break;
                case 2:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0*nummodes1 - 2 - i;
                        }                                               
                    }
                    break;
                case 3:
                    {   
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0*(nummodes1-2-i);
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                } 
                if(edgeOrient == eBackwards)
                {
                    reverse( maparray.get() , maparray.get()+nEdgeIntCoeffs );
                }
            } 
            else
            {
                ASSERTL0(false,"Mapping not defined for this type of basis");
            }

        }

        void StdQuadExp::GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                             Array<OneD, unsigned int> &maparray,
                                             Array<OneD, int> &signarray)
        {
            int i;
            const int nummodes0 = m_base[0]->GetNumModes();
            const int nummodes1 = m_base[1]->GetNumModes();
            const int nEdgeCoeffs = GetEdgeNcoeffs(eid);
            const LibUtilities::BasisType bType = GetEdgeBasisType(eid);
            
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

            if(bType == LibUtilities::eModified_A)
            {
                switch(eid)
                {
                case 0:
                    {
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = i;
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
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = i*nummodes0 + 1;
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
                            maparray[i] = nummodes0+i;
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
                case 3:
                    {
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = i*nummodes0;
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
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                }
            }
            else if(bType == LibUtilities::eGLL_Lagrange)
            {
                switch(eid)
                {
                case 0:
                    {                  
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = i;
                        }
                    }
                    break;
                case 1:
                    {
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = (i+1)*nummodes0 - 1;
                        }                     
                    }
                    break;
                case 2:
                    {
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = nummodes0*nummodes1 - 1 - i;
                        }                                               
                    }
                    break;
                case 3:
                    {   
                        for(i = 0; i < nEdgeCoeffs; i++)
                        {
                            maparray[i] = nummodes0*(nummodes1-1-i);
                        }
                    }
                    break;
                default:
                    ASSERTL0(false,"eid must be between 0 and 3");
                    break;
                } 
                if(edgeOrient == eBackwards)
                {
                    reverse( maparray.get() , maparray.get()+nEdgeCoeffs );
                }
            } 
            else
            {
                ASSERTL0(false,"Mapping not defined for this type of basis");
            }
        }       

        void StdQuadExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
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
                    outfile << "Variables = z1,  z2" ;      
                    outfile << ", "<< var << std::endl << std::endl;
                }
                outfile << "Zone, I=" << nquad0 <<", J=" << nquad1 <<", F=Point" << std::endl;
                
                for(j = 0; j < nquad1; ++j)
                {
                    for(i = 0; i < nquad0; ++i)
                    {
                        outfile << z0[i] <<  " " << 
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
                outfile<<"SQ("<<endl;                
                // write the coordinates of the vertices of the quadrilateral
                outfile<<"-1.0, -1.0, 0.0,"<<endl;
                outfile<<" 1.0, -1.0, 0.0,"<<endl;
                outfile<<" 1.0, 1.0, 0.0,"<<endl;
                outfile<<"-1.0,  1.0, 0.0" <<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j;

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
                    x[i] = xi1[i];
                    y[i] = xi2[i];
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


        void StdQuadExp::ReadFromFile(std::ifstream &infile, OutputFormat format, const bool dumpVar)
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
                        sscanf(str,"%*f %*f %lf",&m_phys[0]+j*nquad0+i);
                    }
                }
            } 
            else
            {
                ASSERTL0(false, "Input routine not implemented for requested type of output");
            }
        }

        //   I/O routine
        void StdQuadExp::WriteCoeffsToFile(std::ofstream &outfile)
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

        void StdQuadExp::GetCoords(Array<OneD, NekDouble> &coords_0, 
                                   Array<OneD, NekDouble> &coords_1)
        {
            Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
            Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
            int nq0 = GetNumPoints(0);
            int nq1 = GetNumPoints(1);
            int i;

            for(i = 0; i < nq1; ++i)
            {
                Blas::Dcopy(nq0,z0.get(), 1,&coords_0[0] + i*nq0,1);
                Vmath::Fill(nq0,z1[i],&coords_1[0] + i*nq0,1);
            }
        }

        void StdQuadExp::IProductWRTBase_SumFacKernel(const Array<OneD, const NekDouble>& base0, 
                                                     const Array<OneD, const NekDouble>& base1,
                                                     const Array<OneD, const NekDouble>& inarray,
                                                     Array<OneD, NekDouble> &outarray,
                                                     Array<OneD, NekDouble> &wsp,
                                                     bool doCheckCollDir0,
                                                     bool doCheckCollDir1)
            {
                int    nquad0 = m_base[0]->GetNumPoints();
                int    nquad1 = m_base[1]->GetNumPoints();   
                int    nmodes0 = m_base[0]->GetNumModes();
                int    nmodes1 = m_base[1]->GetNumModes();

                bool colldir0 = doCheckCollDir0?(m_base[0]->Collocation()):false;
                bool colldir1 = doCheckCollDir1?(m_base[1]->Collocation()):false;
                                
                if(colldir0 && colldir1)
                {
                    Vmath::Vcopy(m_ncoeffs,inarray.get(),1,outarray.get(),1);
                }
                else if(colldir0)
                {
                    Blas::Dgemm('N', 'N',nmodes0,nmodes1, nquad1,1.0, inarray.get(),                            
                                nmodes0, base1.get(), nquad1, 0.0,outarray.get(),nmodes0); 
                }
                else if(colldir1)
                {
                    Blas::Dgemm('T','N',nmodes0,nquad1,nquad0,1.0,base0.get(),               
                                nquad0,inarray.get(),nquad0,0.0,outarray.get(),nmodes0);
                }
                else
                { 
                    ASSERTL1(wsp.num_elements()>=nquad1*nmodes0,"Workspace size is not sufficient");

#if 0 
                    Blas::Dgemm('T','N',nmodes0,nquad1,nquad0,1.0,base0.get(),               
                                nquad0,inarray.get(),nquad0,0.0,wsp.get(),nmodes0);
                    
#else
                    for(int i = 0; i < nmodes0; ++i)
                    {
                        for(int j = 0; j < nquad1; ++j)
                        {
                            wsp[j*nmodes0+i] = Blas::Ddot(nquad0,
                                                          base0.get()+i*nquad0,1,
                                                          inarray.get()+j*nquad0,1);
                        }
                    }
#endif              
                    Blas::Dgemm('N', 'N',nmodes0,nmodes1, nquad1,1.0, wsp.get(),                            
                                nmodes0, base1.get(), nquad1, 0.0,outarray.get(),nmodes0); 
                }
            }        



    } //end namespace            
}//end namespace

/** 
 * $Log: StdQuadExp.cpp,v $
 * Revision 1.54  2009/12/17 01:37:54  bnelson
 * Fixed visual studio compiler warning.
 *
 * Revision 1.53  2009/12/15 18:09:02  cantwell
 * Split GeomFactors into 1D, 2D and 3D
 * Added generation of tangential basis into GeomFactors
 * Updated ADR2DManifold solver to use GeomFactors for tangents
 * Added <GEOMINFO> XML session section support in MeshGraph
 * Fixed const-correctness in VmathArray
 * Cleaned up LocalRegions code to generate GeomFactors
 * Removed GenSegExp
 * Temporary fix to SubStructuredGraph
 * Documentation for GlobalLinSys and GlobalMatrix classes
 *
 * Revision 1.52  2009/09/23 12:42:09  pvos
 * Updates for variable order expansions
 *
 * Revision 1.51  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.50  2009/04/22 22:30:48  sherwin
 * Added ReadFromFile method to read back in .dat file
 *
 * Revision 1.49  2009/01/21 16:58:39  pvos
 * Added additional geometric factors to improve efficiency
 *
 * Revision 1.48  2008/12/16 11:31:52  pvos
 * Performance updates
 *
 * Revision 1.47  2008/11/24 10:31:14  pvos
 * Changed name from _PartitionedOp to _MatFree
 *
 * Revision 1.46  2008/11/19 16:02:47  pvos
 * Added functionality for variable Laplacian coeffcients
 *
 * Revision 1.45  2008/11/05 16:08:15  pvos
 * Added elemental optimisation functionality
 *
 * Revision 1.44  2008/09/23 18:19:26  pvos
 * Updates for working ProjectContField3D demo
 *
 * Revision 1.43  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.42  2008/08/14 22:09:51  sherwin
 * Modifications to remove HDG routines from StdRegions and removed StdExpMap
 *
 * Revision 1.41  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.40  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.39  2008/07/02 14:08:56  pvos
 * Implementation of HelmholtzMatOp and LapMatOp on shape level
 *
 * Revision 1.38  2008/06/05 20:12:51  ehan
 * Removed undefined function GetAlpha() in the ASSERTL2().
 *
 * Revision 1.37  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.36  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.35  2008/05/10 18:27:33  sherwin
 * Modifications necessary for QuadExp Unified DG Solver
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
 * Revision 1.29  2008/03/18 14:15:45  pvos
 * Update for nodal triangular helmholtz solver
 *
 * Revision 1.28  2008/03/12 15:25:09  pvos
 * Clean up of the code
 *
 * Revision 1.27  2008/02/29 19:15:19  sherwin
 * Update for UDG stuff
 *
 * Revision 1.26  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.25  2007/12/06 22:44:47  pvos
 * 2D Helmholtz solver updates
 *
 * Revision 1.24  2007/11/08 16:55:14  pvos
 * Updates towards 2D helmholtz solver
 *
 * Revision 1.23  2007/10/15 20:38:54  ehan
 * Make changes of column major matrix
 *
 * Revision 1.22  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.21  2007/07/12 12:55:16  sherwin
 * Simplified Matrix Generation
 *
 * Revision 1.20  2007/07/11 13:35:18  kirby
 * *** empty log message ***
 *
 * Revision 1.19  2007/07/10 21:05:17  kirby
 * even more fixes
 *
 * Revision 1.17  2007/07/09 15:19:15  sherwin
 * Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
 *
 * Revision 1.16  2007/06/07 15:54:19  pvos
 * Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
 * Also made corrections to various ASSERTL2 calls
 *
 * Revision 1.15  2007/05/15 05:18:24  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.14  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.13  2007/04/06 08:44:43  sherwin
 * Update to make 2D regions work at StdRegions level
 *
 * Revision 1.12  2007/04/05 15:20:11  sherwin
 * Updated 2D stuff to comply with SharedArray philosophy
 *
 * Revision 1.11  2007/04/05 11:39:47  pvincent
 * quad_edited
 *
 * Revision 1.10  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.9  2007/01/20 22:35:21  sherwin
 * Version with StdExpansion compiling
 *
 * Revision 1.8  2007/01/18 18:44:45  bnelson
 * Updates to compile on Visual Studio 2005.
 *
 * Revision 1.7  2007/01/17 16:36:58  pvos
 * updating doxygen documentation
 *
 * Revision 1.6  2007/01/17 16:05:41  pvos
 * updated doxygen documentation
 *
 * Revision 1.5  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.4  2006/08/05 19:03:48  sherwin
 * Update to make the multiregions 2D expansion in connected regions work
 *
 * Revision 1.3  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.2  2006/06/01 14:13:36  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:32  kirby
 * *** empty log message ***
 *
 * Revision 1.38  2006/04/25 20:23:34  jfrazier
 * Various fixes to correct bugs, calls to ASSERT, etc.
 *
 * Revision 1.37  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.36  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.35  2006/03/05 22:11:03  sherwin
 *
 * Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
 *
 * Revision 1.34  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.33  2006/02/26 23:37:30  sherwin
 *
 * Updates and compiling checks upto StdExpansions1D
 *
 **/ 





