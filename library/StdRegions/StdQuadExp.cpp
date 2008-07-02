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
        StdExpansion2D(Ba.GetNumModes()*Bb.GetNumModes(),Ba,Bb)
        { 
        }

        StdQuadExp::StdQuadExp(const StdQuadExp &T):
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


        void StdQuadExp::IProductWRTBase(const Array<OneD, const NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetBdata(),
                inarray,outarray,1);
        }

        void StdQuadExp::IProductWRTDerivBase(const int dir, 
                                             const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> & outarray)
        {
            switch(dir)
            {
            case 0:
                {
                    IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),
                                    inarray,outarray,1);
                }
                break;
            case 1:
                {
                    IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),
                                    inarray,outarray,1);  
                }
                break;
            default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }             
        }
        

        void StdQuadExp:: IProductWRTBase(const Array<OneD, const NekDouble>& base0, 
            const Array<OneD, const NekDouble>& base1,
            const Array<OneD, const NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray,
            int coll_check)
        {
            int i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            Array<OneD, NekDouble> tmp(nquad0*nquad1);

            Array<OneD, const NekDouble> w0 = m_base[0]->GetW();
            Array<OneD, const NekDouble> w1 = m_base[1]->GetW();

//              Does Not Compile
//             ASSERTL2((m_base[0]->GetAlpha() == 0.0)&&(m_base[1]->GetAlpha() == 0.0),
//                      "Basis has non-zero alpha weight");
//             ASSERTL2((m_base[0]->GetBeta() == 0.0)&&(m_base[1]->GetBeta() == 0.0),
//                      "Basis has non-zero beta weight");

            // multiply by integration constants 
            for(i = 0; i < nquad1; ++i)
            {
                Vmath::Vmul(nquad0,(NekDouble*)&inarray[0]+i*nquad0,1,
                    w0.get(),1,&tmp[0]+i*nquad0,1);
            }

            for(i = 0; i < nquad0; ++i)
            {
                Vmath::Vmul(nquad1,&tmp[0]+i,nquad0,w1.get(),1,
                    &tmp[0]+i,nquad0);
            }

            if(coll_check && m_base[0]->Collocation() && m_base[1]->Collocation())
            {               
                Vmath::Vcopy(nquad0*nquad1, tmp, 1, outarray, 1);
            }
            else if( !coll_check || !(m_base[0]->Collocation() || m_base[1]->Collocation()) )
            {
                int    order0 = m_base[0]->GetNumModes();
                int    order1 = m_base[1]->GetNumModes();
                
#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

                Array<OneD, NekDouble> tmp1(nquad0*order1);
                Blas::Dgemm('T','N',order0,nquad1,nquad0,1.0,base0.get(),               
                            nquad0,&tmp[0],nquad0,0.0,&tmp1[0],order0);
                Blas::Dgemm('N', 'N',order0,order1, nquad1,1.0, &tmp1[0],                            
                            order0, base1.get(), nquad1, 0.0,&outarray[0],order0);
                
#else //NEKTAR_USING_DIRECT_BLAS_CALLS

                NekMatrix<const double> in(nquad0,nquad1,tmp,eWrapper);
                NekMatrix<const double> B0(nquad0,order0,base0,eWrapper);
                NekMatrix<const double> B1(nquad1,order1,base1,eWrapper);
                DNekMat out(order0,order1,outarray,eWrapper);
                // out = Transpose(B0)*in*B1; //currently not working with expression templates
                DNekMat tmpMat(nquad0,order1);
                tmpMat = in*B1;
                out = Transpose(B0)*tmpMat;

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS   
 
            }
            else
            {    

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

                if(m_base[0]->Collocation())
                {
                    int order1 = m_base[1]->GetNumModes(); 
                    Blas::Dgemm('N', 'N',nquad0, order1, nquad1, 1.0, &tmp[0],                            
                                nquad0, base1.get(), nquad1, 0.0,&outarray[0],nquad0);
                }
                else
                {
                    int order0 = m_base[0]->GetNumModes();
                    Blas::Dgemm('T','N',order0,nquad1,nquad0,1.0,base0.get(),               
                                nquad0,&tmp[0],nquad0,0.0,&outarray[0],order0);
                }

#else //NEKTAR_USING_DIRECT_BLAS_CALLS

                if(m_base[0]->Collocation())
                {
                    int order1 = m_base[1]->GetNumModes(); 
                    NekMatrix<const double> in(nquad0,nquad1,tmp,eWrapper);
                    NekMatrix<const double> B1(nquad1,order1,base1,eWrapper);
                    DNekMat out(nquad0,order1,outarray,eWrapper);
                    
                    out = in*B1;
                }
                else
                {
                    int order0 = m_base[0]->GetNumModes();
                    NekMatrix<const double> in(nquad0,nquad1,tmp,eWrapper);
                    NekMatrix<const double> B0(nquad0,order0,base0,eWrapper);
                    DNekMat out(order0,nquad1,outarray,eWrapper);
                    
                    out = Transpose(B0)*in;
                }

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS    
            }
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

        DNekMatSharedPtr StdQuadExp::GenMatrix(const StdMatrixKey &mkey)
        {
            int      i;
            int      order0    = GetBasisNumModes(0);
            int      order1    = GetBasisNumModes(1);
            MatrixType mtype   = mkey.GetMatrixType();
            
            DNekMatSharedPtr Mat = StdExpansion::CreateGeneralMatrix(mkey);

            switch(mtype)
            {
            case eMass:
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

            return Mat;
        }

        void StdQuadExp::LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
        {
            cout<<"StdQuadExp::LaplacianMatrixOp"<<endl;
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 

            Array<OneD,NekDouble> physValues(nqtot);
            Array<OneD,NekDouble> dPhysValuesdx(nqtot);
            Array<OneD,NekDouble> dPhysValuesdy(nqtot);

            Array<OneD,NekDouble> wsp(m_ncoeffs);

            BwdTrans(inarray,physValues);

            // Laplacian matrix operation
            PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);
            
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),dPhysValuesdx,outarray,1);
            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,wsp,1);  
            Vmath::Vadd(m_ncoeffs,wsp.get(),1,outarray.get(),1,outarray.get(),1);               
        }

        void StdQuadExp::HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray,
                                           const double lambda)
        {
            int    i;
            int    nquad0 = m_base[0]->GetNumPoints();
            int    nquad1 = m_base[1]->GetNumPoints();
            int    nqtot = nquad0*nquad1; 

            Array<OneD,NekDouble> physValues(nqtot);
            Array<OneD,NekDouble> dPhysValuesdx(nqtot);
            Array<OneD,NekDouble> dPhysValuesdy(nqtot);

            Array<OneD,NekDouble> wsp(m_ncoeffs);

            BwdTrans(inarray,physValues);

            // mass matrix operation
            IProductWRTBase((m_base[0]->GetBdata()),(m_base[1]->GetBdata()),
                            physValues,wsp,1);

            // Laplacian matrix operation
            PhysDeriv(physValues,dPhysValuesdx,dPhysValuesdy);
            
            IProductWRTBase(m_base[0]->GetDbdata(),m_base[1]->GetBdata(),dPhysValuesdx,outarray,1);
            Blas::Daxpy(m_ncoeffs, lambda, wsp.get(), 1, outarray.get(), 1);

            IProductWRTBase(m_base[0]->GetBdata(),m_base[1]->GetDbdata(),dPhysValuesdy,wsp,1);  
            Vmath::Vadd(m_ncoeffs,wsp.get(),1,outarray.get(),1,outarray.get(),1);               
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

        void StdQuadExp::BwdTrans(const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            int           nquad0 = m_base[0]->GetNumPoints();
            int           nquad1 = m_base[1]->GetNumPoints();

            if(m_base[0]->Collocation() && m_base[1]->Collocation())
            {  
                Vmath::Vcopy(nquad0*nquad1, inarray, 1, outarray, 1);
            }
            else if( !(m_base[0]->Collocation() || m_base[1]->Collocation()) )
            { 
                int           order0 = m_base[0]->GetNumModes();
                int           order1 = m_base[1]->GetNumModes();

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

                Array<OneD, NekDouble> tmp(nquad0*order1);
                Blas::Dgemm('N','N', nquad0,order1,order0,1.0, (m_base[0]->GetBdata()).get(),
                            nquad0, &inarray[0], order0,0.0,&tmp[0], nquad0);
                Blas::Dgemm('N','T', nquad0, nquad1,order1, 1.0, &tmp[0],
                            nquad0, (m_base[1]->GetBdata()).get(), nquad1, 0.0, &outarray[0], 
                            nquad0);

#else //NEKTAR_USING_DIRECT_BLAS_CALLS

                NekMatrix<const double> in(order0,order1,inarray,eWrapper);
                NekMatrix<const double> B0(nquad0,order0,m_base[0]->GetBdata(),eWrapper);
                NekMatrix<const double> B1(nquad1,order1,m_base[1]->GetBdata(),eWrapper);
                DNekMat out(nquad0,nquad1,outarray,eWrapper);
                //out = B0*in*Transpose(B1); //(currently not working with expression templates)
                DNekMat tmpM(nquad0,order1);
                tmpM =  B0*in;
                out = tmpM*Transpose(B1);

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS   

            }
            else
            {

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS

                if(m_base[0]->Collocation())
                {
                    int order1 = m_base[1]->GetNumModes();
                    Blas::Dgemm('N','T', nquad0, nquad1,order1, 1.0, &inarray[0],
                                nquad0, (m_base[1]->GetBdata()).get(), nquad1, 0.0, &outarray[0], 
                                nquad0);
                }
                else
                {
                    int order0 = m_base[0]->GetNumModes();
                    Blas::Dgemm('N','N', nquad0,nquad1,order0,1.0, (m_base[0]->GetBdata()).get(),
                            nquad0, &inarray[0], order0,0.0,&outarray[0], nquad0);
                }

#else //NEKTAR_USING_DIRECT_BLAS_CALLS

                if(m_base[0]->Collocation())
                {
                    int order1 = m_base[1]->GetNumModes();
                    NekMatrix<const double> in(nquad0,order1,inarray,eWrapper);
                    NekMatrix<const double> B1(nquad1,order1,m_base[1]->GetBdata(),eWrapper);
                    DNekMat out(nquad0,nquad1,outarray,eWrapper);

                    out = in*Transpose(B1);
                }
                else
                {
                    int order0 = m_base[0]->GetNumModes();
                    NekMatrix<const double> in(order0,nquad1,inarray,eWrapper);
                    NekMatrix<const double> B0(nquad0,order0,m_base[0]->GetBdata(),eWrapper);
                    DNekMat out(nquad0,nquad1,outarray,eWrapper);

                    out = B0*in;
                }

#endif //NEKTAR_USING_DIRECT_BLAS_CALLS 

            }
        }

        void StdQuadExp::FwdTrans(const Array<OneD, const NekDouble>& inarray,
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

        /// Single Point Evaluation
        NekDouble StdQuadExp::PhysEvaluate(Array<OneD, const NekDouble>& coords)
        {
            return  StdExpansion2D::PhysEvaluate2D(coords); 
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
                                            Array<OneD, unsigned int> &maparray)
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
                        
                    }
                    break;
                case 1:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0 + 1;
                        }                        
                    }
                    break;
                case 2:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = nummodes0+i+2;
                        }                        
                    }
                    break;
                case 3:
                    {
                        for(i = 0; i < nEdgeIntCoeffs; i++)
                        {
                            maparray[i] = (i+2)*nummodes0;
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

        // For a specified edge 'eid' this function updates a class
        // StdExpMap which contains the mapping of the edge degrees of
        // freedom back into the elemental domain which is also
        // dependent upon the edge orientation. The vertex and edge
        // ordering of the mapping is dependent upon which basis is
        // being considered, i.e. modal expansions the vertices will
        // be first, nodal expansions the vertices will be the two
        // end points
        void StdQuadExp::MapTo(const int edge_ncoeffs, 
                               const LibUtilities::BasisType Btype, 
                               const int eid, 
                               const EdgeOrientation eorient,
                               StdExpMap &Map)
        {

            int i, start, skip;
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();


            ASSERTL2(eid>=0&&eid <=3,"eid must be between 0 and 3");
            // make sure have correct memory storage
            if(edge_ncoeffs != Map.GetLen())
            {
                Map.SetMapMemory(edge_ncoeffs);
            }

            Array<OneD, int> wsp(edge_ncoeffs);
            Array<OneD, int> sgn = Map.UpdateSign();

            if(eorient == eForwards)
            {
                for(i = 0; i < edge_ncoeffs; ++i)
                {
                    wsp[i] = i;
                    sgn[i] = 1;
                }
            }
            else
            {
                if(Btype == LibUtilities::eGLL_Lagrange)
                {
                    for(i = 0; i < edge_ncoeffs; ++i)
                    {
                        wsp[i] = edge_ncoeffs-i-1;
                        sgn[i] = 1;
                    }
                }
                else{
                    wsp[1] = 0; 
                    wsp[0] = 1;
                    sgn[0] = sgn[1] = 1;
                    
                    int neg = 1;
                    for(i = 2; i < edge_ncoeffs; ++i)
                    {
                        wsp[i] = i;
                        sgn[i] =  neg;
                        neg    = -neg;
                    }
                }
            }
         
            // Set up Mapping details
            if((eid == 0)||(eid == 2))
            { 
                ASSERTL2(Btype == m_base[0]->GetBasisType(),
                    "Expansion type of edge and StdQuadExp are different");

                switch(Btype)
                {
                case LibUtilities::eGLL_Lagrange:
                    ASSERTL2(edge_ncoeffs == order0,
                        "Expansion order of edge and StdQuadExp are different");

                    if(eid == 0)
                    {
                        start = 0;
                        skip  = 1;
                    }
                    else
                    {
                        start = order0*(order1-1);
                        skip = 1;
                    }
                    break;

                case LibUtilities::eModified_A:
                    if(eid == 0)
                    {
                        start = 0;
                        skip  = 1;
                    }
                    else
                    {
                        start = order0;
                        skip = 1;
                    }
                    break;
                default:
                    ASSERTL0(0,"Mapping array is not defined for this expansion");
                    break;
                }
            }
            else
            {
                ASSERTL2(Btype == m_base[1]->GetBasisType(),
                    "Expansion type of edge and StdQuadExp are different");      

                switch(Btype)
                {
                case LibUtilities::eGLL_Lagrange:
                    ASSERTL2(edge_ncoeffs == order1,
                        "Expansion order of edge and StdQuadExp are different");
                    if(eid == 1)
                    {
                        start = order0-1;
                        skip  = order0;
                    }
                    else
                    {
                        start = 0;
                        skip = order0;
                    }
                    break;

                case LibUtilities::eModified_A:    
                    if(eid == 1)
                    {
                        start = 1;
                        skip  = order0;
                    }
                    else
                    {
                        start = 0;
                        skip = order0;
                    }
                    break;
                default:
                    ASSERTL0(0,"Mapping array is not defined for this expansion");
                    break;
                }
            }

            for(i = 0; i < edge_ncoeffs; ++i)
            {
                Map[wsp[i]] = start + i*skip; 
            }

        }

        // same as MapTo but assume that mapping is provided in modal
        // basis type format where vertices are listed first followed
        // by edges degrees of freedom

        void StdQuadExp::MapTo_ModalFormat(const int edge_ncoeffs, 
            const LibUtilities::BasisType Btype, 
            const int eid, 
            const EdgeOrientation eorient,
            StdExpMap &Map)
        {
            MapTo(edge_ncoeffs,Btype,eid,eorient,Map);
            if(Btype == LibUtilities::eGLL_Lagrange)
            {
                int i;
                int vert = Map[edge_ncoeffs-1];
                for(i = edge_ncoeffs-1; i > 1; --i)
                {
                    Map.SetMap(i,Map[i-1]);
                }
                Map.SetMap(1,vert);
            }
        }

        void StdQuadExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar)
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
                    outfile << "Variables = z1,  z2, Coeffs \n" << std::endl;      
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
                    outfile<<"View.MaxRecursionLevel = 8;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
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
    } //end namespace            
}//end namespace

/** 
* $Log: StdQuadExp.cpp,v $
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





