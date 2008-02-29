///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion.cpp
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
// Description: Definition of methods in class StdExpansion which is
// the base class to all expansion shapes
//
///////////////////////////////////////////////////////////////////////////////


#include <StdRegions/StdExpansion.h>

namespace Nektar
{
    namespace StdRegions
    {

        /** define list of number of vertices corresponding to each ShapeType */
        const int g_shapenverts[SIZE_ShapeType] = {0,2,3,4,4,5,6,8};

        /** define list of number of edges corresponding to each ShapeType */
        const int g_shapenedges[SIZE_ShapeType] = {0,1,3,4,6,8,9,12};

        /** define list of number of faces corresponding to each ShapeType */
        const int g_shapenfaces[SIZE_ShapeType] = {0,0,0,0,4,5,5,6};

        StdExpansion::StdExpansion(void): 
        m_ncoeffs(0),
            m_numbases(0)
        {
        }

        StdExpansion::StdExpansion(const int numcoeffs, const int numbases,
            const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb, 
            const LibUtilities::BasisKey &Bc):
            m_ncoeffs(numcoeffs),
            m_numbases(numbases),
            m_stdMatrixManager(std::string("StdExpansion")),
            m_stdStaticCondMatrixManager(std::string("StdExpansionStaticCondMat"))
        {

            m_base = Array<OneD, LibUtilities::BasisSharedPtr>(m_numbases);

            switch(m_numbases)
            {
            case 3:
                ASSERTL2(Bc!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[2] = LibUtilities::BasisManager()[Bc];

            case 2:
                ASSERTL2(Bb!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");

                m_base[1] = LibUtilities::BasisManager()[Bb];
            case 1:
                ASSERTL2(Ba!=LibUtilities::NullBasisKey,
                    "NULL Basis attempting to be used.");
                m_base[0] = LibUtilities::BasisManager()[Ba];
                break;
            default:
                ASSERTL0(false, "numbases incorrectly specified");
            };

            //allocate memory for coeffs
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());

            // Register Creators for  Managers
            for(int i = 0; i < SIZE_MatrixType; ++i)
            {
                m_stdMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoShapeType,*this),
                                  boost::bind(&StdExpansion::CreateStdMatrix, this, _1));
                m_stdStaticCondMatrixManager.RegisterCreator(StdMatrixKey((MatrixType) i,eNoShapeType,*this),  boost::bind(&StdExpansion::CreateStdStaticCondMatrix, this, _1));
            }

        } //end constructor


        StdExpansion::StdExpansion(const StdExpansion &T)
        {
            CopyArray(T.m_base, m_base);

            // NOTE: Copy Constructor produces a deep copy
            // allocate memory for coeffs
            // need to check allocation for variable order. 
            m_ncoeffs = T.m_ncoeffs;
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            Vmath::Vcopy(m_ncoeffs,&T.m_coeffs[0],1,&m_coeffs[0],1);

            //allocate memory for phys
            m_phys = Array<OneD, NekDouble>(GetTotPoints());
            Vmath::Vcopy(GetTotPoints(),&T.m_phys[0],1,&m_phys[0],1);
        }

        StdExpansion::~StdExpansion()
        {
        }


        DNekMatSharedPtr StdExpansion::CreateStdMatrix(const StdMatrixKey &mkey) 
        {
            DNekMatSharedPtr returnval;
            MatrixType mtype = mkey.GetMatrixType();
            
            switch(mtype)
            {
            case eMass: 
            case eLaplacian:
            case eLaplacian00:
            case eLaplacian01:
            case eLaplacian02:
            case eLaplacian11:
            case eLaplacian12:
            case eLaplacian22:
            case eWeakDeriv0:
            case eWeakDeriv1:
            case eWeakDeriv2:
            case eBwdTrans:
            case eNBasisTrans:
                returnval = StdExpansion::GenMatrix(mkey);
                break;

            case eInvMass:
                {
                    StdMatrixKey masskey(eMass,mkey.GetShapeType(),mkey.GetBase(),
                        mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr mmat = m_stdMatrixManager[masskey];
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(*mmat); //Populate standard mass matrix.
                    returnval->Invert();
                }
                break;
                            
            case eInvNBasisTrans:
                {
                    StdMatrixKey tmpkey(eNBasisTrans,mkey.GetShapeType(),mkey.GetBase(),
                        mkey.GetNcoeffs(),mkey.GetNodalPointsType());
                    DNekMatSharedPtr tmpmat = m_stdMatrixManager[tmpkey];
                    DNekLinSys invtmpmat(tmpmat);

                    int dim = tmpmat->GetRows(); //assume square matrix 
                    Array<OneD,NekDouble> invdata = Array<OneD,NekDouble>(dim*dim);
                    Array<OneD,NekDouble> data_offset;

                    Vmath::Zero(dim*dim,&invdata[0],1);

                    for(int i = 0; i < dim; ++i)
                    {
                        // set array to be identity matrix
                        invdata[i*dim + i] = 1.0;

                        //call inverse on symmetric matrix
                        DNekVec   v(dim,data_offset = invdata+i*dim,eWrapper);
                        invtmpmat.Solve(v,v);
                    }

                    NekDouble* t = invdata.data();
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(dim,dim,t);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Matrix creation not defined");
                break;
            }

            return returnval;
        }

        NekDouble StdExpansion::Linf(const ConstArray<OneD, NekDouble>& sol)
        {
            int     ntot;
            NekDouble  val;
            Array<OneD, NekDouble>  wsp;

            ntot = GetTotPoints();
            wsp  = Array<OneD, NekDouble>(ntot);

            Vmath::Vsub(ntot,sol.get(),1,&m_phys[0],1,&wsp[0],1);
            Vmath::Vabs(ntot,&wsp[0],1,&wsp[0],1);
            val = Vmath::Vamax(ntot,&wsp[0],1);

            return  val;
        }

        DNekBlkMatSharedPtr StdExpansion::CreateStdStaticCondMatrix(const StdMatrixKey &mkey) 
        {
            DNekBlkMatSharedPtr returnval;
            MatrixType mtype = mkey.GetMatrixType();
            
            DNekMatSharedPtr  mat = GetStdMatrix(mkey);
            int nbdry = NumBndryCoeffs(); // also checks to see if this is a boundary interior decomposed expansion
            int nint = m_ncoeffs - nbdry;
            DNekMatSharedPtr A = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nbdry);
            DNekMatSharedPtr B = MemoryManager<DNekMat>::AllocateSharedPtr(nbdry,nint);
            DNekMatSharedPtr C = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nbdry);
            DNekMatSharedPtr D = MemoryManager<DNekMat>::AllocateSharedPtr(nint,nint);
            
            int i,j;

            const ConstArray<OneD,int> bmap = GetBoundaryMap();
            const ConstArray<OneD,int> imap = GetInteriorMap();
            
            for(i = 0; i < nbdry; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*A)(i,j) = (*mat)(bmap[i],bmap[j]);
                }
                
                for(j = 0; j < nint; ++j)
                {
                    (*B)(i,j) = (*mat)(bmap[i],imap[j]);
                }
            }
            
            for(i = 0; i < nint; ++i)
            {
                for(j = 0; j < nbdry; ++j)
                {
                    (*C)(i,j) = (*mat)(imap[i],bmap[j]);
                }
                
                for(j = 0; j < nint; ++j)
                {
                    (*D)(i,j) = (*mat)(imap[i],imap[j]);
                }
            }
            
            // Calculate static condensed system 
            if(nint)
            {
                D->Invert();
                (*B) = (*B)*(*D);
                (*A) = (*A) - (*B)*(*C);
            }

            // set up block matrix system
            unsigned int exp_size[] = { nbdry,nint};
            int nblks = 2;
            returnval = MemoryManager<DNekBlkMat>::AllocateSharedPtr(nblks,nblks,exp_size,exp_size); //Really need a constructor which takes Array<OneD,int>
            
            returnval->SetBlock(0,0,A);
            returnval->SetBlock(0,1,B);
            returnval->SetBlock(1,0,C);
            returnval->SetBlock(1,1,D);
            
            return returnval;
        }

        NekDouble StdExpansion::Linf()
        {
            return Vmath::Vamax(GetTotPoints(),&m_phys[0],1);    
        }

        NekDouble StdExpansion::L2(const ConstArray<OneD, NekDouble>& sol)
        {
            int     ntot = GetTotPoints();
            NekDouble  val;
            Array<OneD, NekDouble> wsp;

            wsp = Array<OneD, NekDouble>(ntot);

            Vmath::Vsub(ntot, sol.get(), 1, m_phys.get(), 1, wsp.get(), 1);
            Vmath::Vmul(ntot, wsp.get(), 1, wsp.get(),  1, wsp.get(), 1);

            val = sqrt(v_Integral(wsp));

            return val;
        }

        NekDouble StdExpansion::L2()
        {
            int     ntot = GetTotPoints();
            NekDouble  val;
            Array<OneD, NekDouble> wsp;

            wsp = Array<OneD, NekDouble>(ntot);

            Vmath::Vmul(ntot, &m_phys[0], 1, &m_phys[0], 1, &wsp[0], 1);
            val   = sqrt(v_Integral(wsp));

            return val;
        }

        DNekMatSharedPtr StdExpansion::CreateGeneralMatrix(const StdMatrixKey &mkey)
        {
            int     i;
            DNekMatSharedPtr  returnval;


            switch(mkey.GetMatrixType())
            {
            case StdRegions::eUnifiedDGHelmholtz:
                {

                    ASSERTL1(IsBoundaryInteriorExpansion(),
                             "UnifiedDGHelmholtz matrix not set up "
                             "for non boundary-interior expansions");

                    int j;
                    NekDouble tau = mkey.GetConstant(1);

                    // Get basic Galerkin Helmholtz matrix 
                    StdRegions::StdMatrixKey hkey(eHelmholtz,DetShapeType(),
                                                  *this,mkey.GetConstant(0));
                    returnval  = GenMatrix(hkey);
                    DNekMat &Mat = *returnval;

                    Array<OneD,NekDouble> inarray(m_ncoeffs);
                    Array<OneD,NekDouble> outarray(m_ncoeffs);

                    for(j = 0; j < m_ncoeffs; ++j)
                    {
                        Vmath::Zero(m_ncoeffs,&inarray[0],1);
                        Vmath::Zero(m_ncoeffs,&outarray[0],1);
                        inarray[j] = 1.0;
                        
                        v_AddUDGHelmholtzBoundaryTerms(tau,inarray,outarray);

                        for(i = 0; i < m_ncoeffs; ++i)
                        {
                            Mat(i,j) += outarray[i];
                        }
                    }
                }
                break;
            case StdRegions::eUnifiedDGLamToU:
                {
                    int i,j,k;
                    int nbndry = NumDGBndryCoeffs();
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);

                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    
                    Array<OneD,NekDouble> ulam(m_ncoeffs);
                    DNekVec Ulam(m_ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(m_ncoeffs);
                    DNekVec F(m_ncoeffs,f,eWrapper);


                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,nbndry); 
                    DNekMat &Umat = *returnval;


                    // get mapping for boundary dof
                    ConstArray<OneD,int> bmap = GetBoundaryMap();

                    // Helmholtz matrix
                    DNekScalMat  &invHmat = *GetLocMatrix(eInvUnifiedDGHelmholtz,
                                                          lambdaval, tau);

                    // for each degree of freedom of the lambda space
                    // calculate Umat entry 
                    // Generate Lambda to U_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(m_ncoeffs,&f[0],1);
                        lambda[j] = 1.0;
                        
                        v_AddUDGHelmholtzTraceTerms(tau,lambda,f);
                        
                        Ulam = invHmat*F; // generate Ulam from lambda

                        // fill column of matrix
                        for(k = 0; k < m_ncoeffs; ++k)
                        {
                            Umat(k,j) = Ulam[k]; 
                        }
                    }
                }
                break;
            case StdRegions::eUnifiedDGLamToQ0:
            case StdRegions::eUnifiedDGLamToQ1:
            case StdRegions::eUnifiedDGLamToQ2:
                {
                    int i,j,k,dir;
                    int nbndry = NumDGBndryCoeffs();
                    int nquad  = GetNumPoints(0);
                    const ConstArray<OneD,NekDouble> &Basis  = m_base[0]->GetBdata();
                    Array<OneD,NekDouble> lambda(nbndry);
                    DNekVec Lambda(nbndry,lambda,eWrapper);                    

                    Array<OneD,NekDouble> ulam(m_ncoeffs);
                    DNekVec Ulam(m_ncoeffs,ulam,eWrapper);
                    Array<OneD,NekDouble> f(m_ncoeffs);
                    DNekVec F(m_ncoeffs,f,eWrapper);
                    NekDouble lambdaval = mkey.GetConstant(0);
                    NekDouble tau       = mkey.GetConstant(1);

                    // declare matrix space
                    returnval  = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,nbndry); 
                    DNekMat &Qmat = *returnval;
                    
                    // get mapping for boundary dof
                    ConstArray<OneD,int> bmap = GetBoundaryMap();
                    
                    // Helmholtz matrix
                    DNekScalMat &invHmat = *GetLocMatrix(eInvUnifiedDGHelmholtz, lambdaval,tau);

                    // Inverse mass matrix 
                    DNekScalMat &invMass = *GetLocMatrix(eInvMass);
                    
                    //Weak Derivative matrix 
                    DNekScalMatSharedPtr Dmat;
                    switch(mkey.GetMatrixType())
                    {
                    case eUnifiedDGLamToQ0:
                        dir = 0;
                        Dmat = GetLocMatrix(eWeakDeriv0); 
                        break;
                    case eUnifiedDGLamToQ1:
                        dir = 1;
                        Dmat = GetLocMatrix(eWeakDeriv1); 
                        break;
                    case eUnifiedDGLamToQ2:
                        dir = 2;
                        Dmat = GetLocMatrix(eWeakDeriv2); 
                        break;
                    }

                    // for each degree of freedom of the lambda space
                    // calculate Qmat entry 
                    // Generate Lambda to Q_lambda matrix 
                    for(j = 0; j < nbndry; ++j)
                    {
                        Vmath::Zero(nbndry,&lambda[0],1);
                        Vmath::Zero(m_ncoeffs,&f[0],1);
                        lambda[j] = 1.0;
                        
                        v_AddUDGHelmholtzTraceTerms(tau,lambda,f);
                        
                        Ulam = invHmat*F; // generate Ulam from lambda

                        F = (*Dmat)*Ulam; 
                        

                        v_AddNormTraceInt(dir,lambda,f); // add \tilde{G}\lambda

                        Vmath::Neg(m_ncoeffs,&ulam[0],1);
                        v_AddNormBoundaryInt(dir,ulam,f); // subtrace G Ulam

                        Ulam = invMass*F; // multiply by inverse mass matrix
                        
                        // fill column of matrix (Qmat is in column major format)
                        Vmath::Vcopy(m_ncoeffs,&ulam[0],1,&(Qmat.GetPtr())[0]+j*m_ncoeffs,1);

                    }
                }
                break;
            default:
                {
                    Array<OneD, NekDouble> tmp = Array<OneD, NekDouble>(m_ncoeffs);
                    returnval = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);            
                    DNekMat &Mat = *returnval; 
                    
                    for(i=0; i < m_ncoeffs; ++i)
                    {
                        Vmath::Zero(m_ncoeffs,&tmp[0],1);
                        tmp[i] = 1.0;
                        
                        GeneralMatrixOp(mkey,tmp,tmp);
                        
                        Vmath::Vcopy(m_ncoeffs,&tmp[0],1,&(Mat.GetPtr())[0]+i*m_ncoeffs,1);
                    }
                }
                break;
            }
            
            return returnval;
        }

        void StdExpansion::GeneralMatrixOp(const StdMatrixKey &mkey, 
                                           const ConstArray<OneD,NekDouble> &inarray,
                                           Array<OneD,NekDouble> &outarray)
        {
            switch(mkey.GetMatrixType())
            {
            case eMass:
                MassMatrixOp(inarray,outarray);
                break;
            case eWeakDeriv0:
                WeakDerivMatrixOp(0,inarray,outarray);
                break;
            case eWeakDeriv1:
                WeakDerivMatrixOp(1,inarray,outarray);
                break;
            case eWeakDeriv2:
                WeakDerivMatrixOp(2,inarray,outarray);
                break;
            case eLaplacian:
                LaplacianMatrixOp(inarray,outarray);
                break;
            case eLaplacian00:
                LaplacianMatrixOp(0,0,inarray,outarray);
                break;
            case eLaplacian01:
                LaplacianMatrixOp(0,1,inarray,outarray);
                break;
            case eLaplacian11:
                LaplacianMatrixOp(1,1,inarray,outarray);
                break;
            case eLaplacian22:
                LaplacianMatrixOp(2,2,inarray,outarray);
                break;
            case eBwdTrans:
                BwdTransMatrixOp(inarray,outarray);
                break;
            case eHelmholtz:
                HelmholtzMatrixOp(inarray,outarray,mkey.GetConstant(0));
                break;
            case eNBasisTrans:
                NEKERROR(ErrorUtil::efatal,"This is a specialised matrix for nodal "
                         "expansions only and does not have an operator");
                break;

            default:
                NEKERROR(ErrorUtil::efatal, "Need to populate switch");
                break;
            }
        }
            
        void StdExpansion::MassMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                        Array<OneD,NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp   = Array<OneD, NekDouble>(GetTotPoints());

            v_BwdTrans(inarray,tmp);            
            v_IProductWRTBase(tmp, outarray);
        }

        void StdExpansion::LaplacianMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {

            switch(ShapeTypeDimMap[v_DetShapeType()])
            {
            case 1:
                LaplacianMatrixOp(0,0,inarray,outarray);
                break;

            case 2:
                {
                    Array<OneD, NekDouble> store = Array<OneD, NekDouble>(m_ncoeffs);
                    
                    LaplacianMatrixOp(0,0,inarray,store);
                    LaplacianMatrixOp(1,1,inarray,outarray);
                   
                    Vmath::Vadd(m_ncoeffs,&store[0],1,&outarray[0],1,&outarray[0],1);
                }
                break;
            case 3:
                {
                    Array<OneD, NekDouble> store0 = Array<OneD, NekDouble>(m_ncoeffs);
                    Array<OneD, NekDouble> store1 = Array<OneD, NekDouble>(m_ncoeffs);
                    
                    LaplacianMatrixOp(0,0,inarray,store0);
                    LaplacianMatrixOp(1,1,inarray,store1);
                    LaplacianMatrixOp(2,2,inarray,outarray);
                    
                    Vmath::Vadd(m_ncoeffs,&store0[0],1,&outarray[0],1,&outarray[0],1);
                    Vmath::Vadd(m_ncoeffs,&store1[0],1,&outarray[0],1,&outarray[0],1);
                }
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "Dimenion not recognised.");
                break;
            }
        }


        void StdExpansion::LaplacianMatrixOp(const int k1, const int k2, 
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
                
            ASSERTL1(k1 >= 0 && k1 < ShapeTypeDimMap[v_DetShapeType()],"invalid first  argument");
            ASSERTL1(k2 >= 0 && k2 < ShapeTypeDimMap[v_DetShapeType()],"invalid second argument");
                                  
            Array<OneD, NekDouble> tmp   = Array<OneD, NekDouble>(GetTotPoints());
            Array<OneD, NekDouble> dtmp  = Array<OneD, NekDouble>(GetTotPoints());
            
            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k2,tmp,dtmp);
            v_IProductWRTDerivBase(k1, dtmp, outarray);
        }


        void StdExpansion::WeakDerivMatrixOp(const int k1,
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
            ASSERTL1(k1 >= 0 && k1 < ShapeTypeDimMap[v_DetShapeType()],"invalid first  argument");
                                  
            Array<OneD, NekDouble> tmp   = Array<OneD, NekDouble>(GetTotPoints());
            v_BwdTrans(inarray,tmp);
            v_PhysDeriv(k1,tmp,tmp);
            v_IProductWRTBase(tmp, outarray);
        }


        void StdExpansion::BwdTransMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
        {
            v_BwdTrans(inarray,outarray);
        }


        void StdExpansion::HelmholtzMatrixOp(const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const double lambda)
        {
            Array<OneD,NekDouble> tmp = Array<OneD,NekDouble>(m_ncoeffs);
            MassMatrixOp(inarray,tmp);
            LaplacianMatrixOp(inarray,outarray);

            Blas::Daxpy(m_ncoeffs,lambda,&tmp[0],1,&outarray[0],1);
        }

        void StdExpansion::TensProdBwdTrans(const ConstArray<OneD, NekDouble>& inarray, 
            Array<OneD, NekDouble> &outarray)
        {
            int nq = GetTotPoints();
            StdMatrixKey      bwdtransmatkey(eBwdTrans,DetShapeType(),*this);
            DNekMatSharedPtr  bwdtransmat = m_stdMatrixManager[bwdtransmatkey];

            DNekVec v_in(m_ncoeffs,inarray);
            DNekVec v_out(nq,outarray,eWrapper);

            v_out = (*bwdtransmat) * v_in;
            // This line below should be removed once the eWrapper method of NekVEctor works properly
            Vmath::Vcopy(nq,&((v_out).GetPtr())[0],1,&outarray[0],1);
        }


        // 3D interpolation
	 void StdExpansion::Interp3D(const LibUtilities::BasisKey &fbasis0, 
				     const LibUtilities::BasisKey &fbasis1, 
				     const LibUtilities::BasisKey &fbasis2, 
				     const ConstArray<OneD, NekDouble>& from,  
				     const LibUtilities::BasisKey &tbasis0,
				     const LibUtilities::BasisKey &tbasis1,
	                             const LibUtilities::BasisKey &tbasis2,
				     Array<OneD, NekDouble> &to)
        {
            Interp3D(fbasis0,fbasis1,fbasis2,from.data(),tbasis0,tbasis1,tbasis2,to.data());
        }

	void StdExpansion::Interp3D(const LibUtilities::BasisKey &fbasis0, 
				    const LibUtilities::BasisKey &fbasis1,
				    const LibUtilities::BasisKey &fbasis2,  
				    const NekDouble *from,  
				    const LibUtilities::BasisKey &tbasis0,
				    const LibUtilities::BasisKey &tbasis1,
				    const LibUtilities::BasisKey &tbasis2,
				    NekDouble *to)
        {
            DNekMatSharedPtr I0, I1, I2;

            Array<OneD, NekDouble> wsp1 = Array<OneD, NekDouble>(tbasis1.GetNumPoints()*fbasis0.GetNumPoints());
            Array<OneD, NekDouble> wsp2 = Array<OneD, NekDouble>(tbasis2.GetNumPoints()*fbasis1.GetNumPoints());

            I0 = LibUtilities::PointsManager()[fbasis0.GetPointsKey()]->GetI(tbasis0.GetPointsKey());
            I1 = LibUtilities::PointsManager()[fbasis1.GetPointsKey()]->GetI(tbasis1.GetPointsKey());
            I2 = LibUtilities::PointsManager()[fbasis2.GetPointsKey()]->GetI(tbasis2.GetPointsKey());

            //TODO : following matrix storiges (rows, columns)
            Blas::Dgemm('N', 'T',
                         tbasis2.GetNumPoints(),      //  3
                         fbasis1.GetNumPoints(),      //  4
                         fbasis2.GetNumPoints(),      //  5
                         1.0,                         //  6
                         I2->GetPtr().get(),          //  7
                         tbasis1.GetNumPoints(),      //  8
                         from,                        //  9
                         fbasis0.GetNumPoints(),      // 10 
                         0.0,                         // 11
                         wsp2.get(),                  // 12
                         tbasis1.GetNumPoints());     // 13
              //TODO : following matrix storiges (rows, columns)
            Blas::Dgemm('N', 'T',
                         tbasis1.GetNumPoints(),      //  3
                         fbasis0.GetNumPoints(),      //  4
                         fbasis1.GetNumPoints(),      //  5
                         1.0,                         //  6
                         I1->GetPtr().get(),          //  7
                         tbasis1.GetNumPoints(),      //  8
                         wsp2.get(),                  //  9
                         fbasis0.GetNumPoints(),      // 10 
                         0.0,                         // 11
                         wsp1.get(),                  // 12
                         tbasis1.GetNumPoints());     // 13

              //TODO : following matrix storiges (rows, columns)
            Blas::Dgemm('N', 'T',
                        tbasis0.GetNumPoints(),      //  3
                        tbasis1.GetNumPoints(),      //  4
                        fbasis0.GetNumPoints(),      //  5
                        1.0,                         //  6
                        I0->GetPtr().get(),          //  7
                        tbasis0.GetNumPoints(),      //  8 
                        wsp1.get(),                  //  9
                        tbasis1.GetNumPoints(),      // 10
                        0.0,                         // 11
                        to,                          // 12
                        tbasis0.GetNumPoints());     // 13
        }

        // 2D Interpolation
        void StdExpansion::Interp2D(const LibUtilities::BasisKey &fbasis0, 
            const LibUtilities::BasisKey &fbasis1, 
            const ConstArray<OneD, NekDouble>& from,  
            const LibUtilities::BasisKey &tbasis0,
            const LibUtilities::BasisKey &tbasis1,
            Array<OneD, NekDouble> &to)
        {
            Interp2D(fbasis0,fbasis1,from.data(),tbasis0,tbasis1,to.data());
        }

        void StdExpansion::Interp2D(const LibUtilities::BasisKey &fbasis0, 
            const LibUtilities::BasisKey &fbasis1, 
            const NekDouble *from,  
            const LibUtilities::BasisKey &tbasis0,
            const LibUtilities::BasisKey &tbasis1,
            NekDouble *to)
        {
            DNekMatSharedPtr I0,I1;
            Array<OneD, NekDouble> wsp = Array<OneD, NekDouble>(tbasis1.GetNumPoints()*fbasis0.GetNumPoints());

            I0 = LibUtilities::PointsManager()[fbasis0.GetPointsKey()]->GetI(tbasis0.GetPointsKey());
            I1 = LibUtilities::PointsManager()[fbasis1.GetPointsKey()]->GetI(tbasis1.GetPointsKey());

            Blas::Dgemm('N', 'T', 
                         tbasis1.GetNumPoints(),
                         fbasis0.GetNumPoints(),      
                         fbasis1.GetNumPoints(),     
                         1.0,                        
                         I1->GetPtr().get(),          
                         tbasis1.GetNumPoints(),      
                         from,                       
                         fbasis0.GetNumPoints(),      
                         0.0,                         
                         wsp.get(),                  
                         tbasis1.GetNumPoints());     

            Blas::Dgemm('N', 'T',
                        tbasis0.GetNumPoints(),     
                        tbasis1.GetNumPoints(),      
                        fbasis0.GetNumPoints(),      
                        1.0,                         
                        I0->GetPtr().get(),         
                        tbasis0.GetNumPoints(),      
                        wsp.get(),                   
                        tbasis1.GetNumPoints(),      
                        0.0,                         
                        to,                          
                        tbasis0.GetNumPoints());     
        }

        // 1D Interpolation
        void StdExpansion::Interp1D(const LibUtilities::BasisKey &fbasis0, 
            const ConstArray<OneD, NekDouble>& from,  
            const LibUtilities::BasisKey &tbasis0, 
            Array<OneD, NekDouble> &to)
        {
            Interp1D(fbasis0, from.get(), tbasis0, to.get());
        }

        void StdExpansion::Interp1D(const LibUtilities::BasisKey &fbasis0, 
            const NekDouble *from,  
            const LibUtilities::BasisKey &tbasis0, 
            NekDouble *to)
        {
            DNekMatSharedPtr I0;

            I0 = LibUtilities::PointsManager()[fbasis0.GetPointsKey()]
            ->GetI(tbasis0.GetPointsKey());

            //DNekVec in(fbasis0.GetNumPoints(),from);
            //DNekVec out(tbasis0.GetNumPoints(),to,eWrapper);

            //out  = (*I0)*in;
            // this line should not be needed
            //Vmath::Vcopy(tbasis0.GetNumPoints(),&out[0],1,&to[0],1);

            Blas::Dgemv('N', tbasis0.GetNumPoints(), fbasis0.GetNumPoints(), 
                1.0, I0->GetPtr().get(), tbasis0.GetNumPoints(), 
                from, 1, 0.0, to, 1);
        }

        //   I/O routine
        void StdExpansion::WriteCoeffsToFile(std::ofstream &outfile)
        {
            int i;
            for(i=0; i<m_ncoeffs; ++i)
            {
                outfile << m_coeffs[i] << std::endl;
            }
        }

    }//end namespace
}//end namespace

/**
* $Log: StdExpansion.cpp,v $
* Revision 1.62  2008/02/16 05:59:14  ehan
* Added interpolation 3D.
*
* Revision 1.61  2008/01/23 09:09:46  sherwin
* Updates for Hybrized DG
*
* Revision 1.60  2007/12/17 13:03:45  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.59  2007/12/06 22:44:46  pvos
* 2D Helmholtz solver updates
*
* Revision 1.58  2007/11/29 21:40:20  sherwin
* updates for MultiRegions and DG solver
*
* Revision 1.57  2007/11/08 16:55:12  pvos
* Updates towards 2D helmholtz solver
*
* Revision 1.56  2007/10/15 20:38:32  ehan
* Tested standard mass matrix
*
* Revision 1.55  2007/10/04 12:10:04  sherwin
* Update for working version of static condensation in Helmholtz1D and put lambda coefficient on the mass matrix rather than the Laplacian operator.
*
* Revision 1.54  2007/10/03 11:37:51  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.53  2007/09/27 12:55:57  pvos
* Column major Blas calls corrections
*
* Revision 1.52  2007/09/25 14:25:56  pvos
* Update for helmholtz1D with different expansion orders
*
* Revision 1.51  2007/08/29 23:26:48  jfrazier
* Created non-static manager that shares data across instances.
*
* Revision 1.50  2007/07/27 16:56:50  jfrazier
* Changed manager to static.
*
* Revision 1.49  2007/07/27 00:22:53  bnelson
* Memory manager now accepts non-const parameters to the allocate methods.
*
* Revision 1.48  2007/07/22 23:04:25  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.47  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.46  2007/07/15 19:28:28  bnelson
* *** empty log message ***
*
* Revision 1.45  2007/07/13 15:20:19  kirby
* *** empty log message ***
*
* Revision 1.43  2007/07/13 09:02:25  sherwin
* Mods for Helmholtz solver
*
* Revision 1.42  2007/07/12 12:55:14  sherwin
* Simplified Matrix Generation
*
* Revision 1.41  2007/07/10 20:41:52  kirby
* more fixes
*
* Revision 1.40  2007/07/10 19:27:58  kirby
* Update for new matrix structures
*
* Revision 1.39  2007/07/09 15:19:14  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.38  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.37  2007/05/30 23:56:54  sherwin
* Silly errors
*
* Revision 1.36  2007/05/30 20:49:12  sherwin
* Updates to do with LocalRegions and SpatialDomains
*
* Revision 1.35  2007/05/23 15:12:45  pvos
* removed some obsolete lines
*
* Revision 1.34  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.33  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.32  2007/04/18 16:09:12  pvos
* Added some new Tensor Operations routines
*
* Revision 1.31  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.30  2007/04/08 03:36:57  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.29  2007/03/31 00:40:02  bnelson
* *** empty log message ***
*
* Revision 1.28  2007/03/29 19:35:08  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.27  2007/03/25 15:48:22  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.26  2007/03/21 20:56:42  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.25  2007/03/20 16:58:42  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.24  2007/03/14 21:24:09  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.23  2007/03/02 12:01:51  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.22  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.21  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.20  2007/02/23 19:26:07  jfrazier
* General bug fix and formatting.
*
* Revision 1.19  2007/02/22 22:02:27  sherwin
* Update with executing StdMatManager
*
* Revision 1.18  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.17  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
* Revision 1.16  2007/02/17 04:03:22  jfrazier
* Added NekManager for holding matrices.  Need to finish the create function.
*
* Revision 1.15  2007/02/14 16:35:50  pvos
* Corrected an error in the code
*
* Revision 1.14  2007/02/13 09:52:27  sherwin
* Updates to fix mass matrix inverse issues
*
* Revision 1.13  2007/02/07 12:51:52  sherwin
* Compiling version of Project1D
*
* Revision 1.12  2007/02/06 02:23:28  jfrazier
* Minor cleanup.
*
* Revision 1.11  2007/01/30 20:01:35  sherwin
* Update for first compiling Project1D routine
*
* Revision 1.10  2007/01/29 15:04:53  sherwin
* StdBasis.h moved to LibUtilities. Other minor mods
*
* Revision 1.9  2007/01/28 18:34:18  sherwin
* More modifications to make Demo Project1D compile
*
* Revision 1.8  2007/01/23 23:20:20  sherwin
* New version after Jan 07 update
*
* Revision 1.7  2007/01/20 22:35:20  sherwin
* Version with StdExpansion compiling
*
* Revision 1.6  2007/01/15 11:08:37  pvos
* Updating doxygen documentation
*
* Revision 1.5  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.4  2006/08/05 19:03:48  sherwin
* Update to make the multiregions 2D expansion in connected regions work
*
* Revision 1.3  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.2  2006/05/29 19:03:08  sherwin
* Modifications to wrap geometric information in shared_ptr
*
* Revision 1.1  2006/05/04 18:58:31  kirby
* *** empty log message ***
*
* Revision 1.54  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.53  2006/04/01 21:59:26  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.52  2006/03/21 09:21:31  sherwin
* Introduced NekMemoryManager
*
* Revision 1.51  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.50  2006/03/06 12:39:59  sherwin
*
* Added NekConstants class for all constants in this library
*
* Revision 1.49  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.48  2006/03/03 23:04:54  sherwin
*
* Corrected Mistake in StdBasis.cpp to do with eModified_B
*
* Revision 1.47  2006/03/02 16:20:20  sherwin
*
* Introduced method GetPointsTot
*
* Revision 1.46  2006/03/01 08:25:03  sherwin
*
* First compiling version of StdRegions
*
* Revision 1.45  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.44  2006/02/26 21:23:20  bnelson
* Fixed a variety of compiler errors caused by updates to the coding standard.
*
* Revision 1.43  2006/02/15 08:06:36  sherwin
*
* Put files into coding standard (although they do not compile)
*
* Revision 1.42  2006/02/12 21:51:42  sherwin
*
* Added licence
*
* Revision 1.41  2006/02/10 16:44:10  sherwin
*
* Updated to comply with coding standard
*
**/
