///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.cpp
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
// Description: Expansion list definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

namespace Nektar
{
    namespace MultiRegions
    {
        ExpList::ExpList(void):
            m_ncoeffs(0),
            m_npoints(0),
            m_coeffs(),
            m_phys(),
            m_transState(eNotSet),
            m_physState(false),
            m_exp(MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr()),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr())
        {            
        }
        
        ExpList::ExpList(const ExpList &in):
            m_ncoeffs(in.m_ncoeffs),
            m_npoints(in.m_npoints),
            m_coeffs(m_ncoeffs),
            m_phys(m_npoints),
            m_transState(eNotSet),
            m_physState(false),
            m_exp(in.m_exp),
            m_coeff_offset(in.m_coeff_offset), // Need to check if we need these
            m_phys_offset(in.m_phys_offset)    // or at least use shared pointer
        {
        }


        
        void ExpList::SetCoeffPhys(void)
        {
            int i;
            
            // Set up offset information and array sizes
            m_coeff_offset = Array<OneD,int>(m_exp->size());
            m_phys_offset  = Array<OneD,int>(m_exp->size());
            
            m_ncoeffs = m_npoints = 0;
            
            for(i = 0; i < m_exp->size(); ++i)
            {
                m_coeff_offset[i] = m_ncoeffs;
                m_phys_offset [i] = m_npoints; 
                m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                m_npoints += (*m_exp)[i]->GetNumPoints(0);
            }
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }



        void ExpList::PutCoeffsInToElmtExp(void)
        {
            int i, order_e;
            int cnt = 0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                order_e = (*m_exp)[i]->GetNcoeffs();
                Vmath::Vcopy(order_e,&m_coeffs[cnt], 1, 
                             &((*m_exp)[i]->UpdateCoeffs())[0],1);
                cnt += order_e;
            }
        }

        void ExpList::PutCoeffsInToElmtExp(int eid)
        {
            int order_e;
            int cnt = 0;

            order_e = (*m_exp)[eid]->GetNcoeffs();
            cnt = m_coeff_offset[eid];
            Vmath::Vcopy(order_e,&m_coeffs[cnt], 1, 
                         &((*m_exp)[eid]->UpdateCoeffs())[0],1);
        }


        void ExpList::PutElmtExpInToCoeffs(void)
        {
            int i, order_e;
            int cnt = 0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                order_e = (*m_exp)[i]->GetNcoeffs();
                Vmath::Vcopy(order_e, &((*m_exp)[i]->UpdateCoeffs())[0],1,
                             &m_coeffs[cnt],1);
                cnt += order_e;
            }
        }


        void ExpList::PutElmtExpInToCoeffs(int eid)
        {
            int order_e;
            int cnt = 0;

            order_e = (*m_exp)[eid]->GetNcoeffs();
            cnt = m_coeff_offset[eid];
            
            Vmath::Vcopy(order_e, &((*m_exp)[eid]->UpdateCoeffs())[0],1,
                             &m_coeffs[cnt],1);
        }
        
        void ExpList::PutPhysInToElmtExp(Array<OneD,const NekDouble> &in)
        {
            int i, npoints_e;
            int cnt = 0;
            
            for(i = 0; i < (*m_exp).size(); ++i)
            {
                npoints_e = (*m_exp)[i]->GetTotPoints();
                Vmath::Vcopy(npoints_e, &in[cnt],1, &((*m_exp)[i]->UpdatePhys())[0],1);
                cnt += npoints_e;
            }
        }

        void ExpList::PutElmtExpInToPhys(Array<OneD,NekDouble> &out)
        {
            int i, npoints_e;
            int cnt = 0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                npoints_e = (*m_exp)[i]->GetTotPoints();
                Vmath::Vcopy(npoints_e, &((*m_exp)[i]->GetPhys())[0],1,
                             &out[cnt],1);
                cnt += npoints_e;
            }
        }

        void ExpList::PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out)
        {
            int npoints_e;
            int cnt = m_phys_offset[eid];

            npoints_e = (*m_exp)[eid]->GetTotPoints();
            Vmath::Vcopy(npoints_e, &((*m_exp)[eid]->GetPhys())[0],1,
                         &out[cnt],1);            
        }

        ExpList::~ExpList()
        {
        }
    
        NekDouble ExpList::PhysIntegral(void)
        {
            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            return PhysIntegral(m_phys);
        }
        
        NekDouble ExpList::PhysIntegral(const Array<OneD, const NekDouble> &inarray)
        {
            int       i;
            int       cnt = 0;
            NekDouble sum = 0.0;
            
            for(i = 0; i < GetExpSize(); ++i)
            {
                sum += (*m_exp)[i]->Integral(inarray + cnt);
                cnt += (*m_exp)[i]->GetTotPoints();
            }
            
            return sum; 
        }

        void ExpList::MultiplyByBlockMatrix(const GlobalMatrixKey             &gkey,
                                            const Array<OneD,const NekDouble> &inarray, 
                                                  Array<OneD,      NekDouble> &outarray)
        {
            int nrows;
            int ncols;
            const DNekScalBlkMatSharedPtr& blockmat = GetBlockMatrix(gkey);

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS     
     
            int n,cnt1,cnt2;
            int nblocks = blockmat->GetNumberOfBlockRows();

            for(n = cnt1 = cnt2 = 0; n < nblocks; n++)
            {
                DNekScalMat& elmMat = *(blockmat->GetBlock(n,n));

                nrows = elmMat.GetRows();
                ncols = elmMat.GetColumns();

                Blas::Dgemv('N',nrows,ncols,elmMat.Scale(),elmMat.GetRawPtr(),
                            nrows, inarray.get()+cnt2, 1.0, 0.0, outarray.get()+cnt1, 1.0);

                cnt1 += nrows;
                cnt2 += ncols;
            }

#else

            nrows = blockmat->GetRows();
            ncols = blockmat->GetColumns();

            NekVector<const NekDouble> in (ncols,inarray, eWrapper);
            NekVector<      NekDouble> out(nrows,outarray,eWrapper); 
            out = (*blockmat)*in;

#endif
        }
        
        void ExpList::IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, 
                                      Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int    cnt  = 0;
            int    cnt1 = 0;
            
            Array<OneD,NekDouble> e_outarray;
            
            for(i = 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->IProductWRTBase(inarray+cnt,
                                             e_outarray = outarray+cnt1);
                cnt  += (*m_exp)[i]->GetTotPoints();
                cnt1 += (*m_exp)[i]->GetNcoeffs();
            }
            m_transState = eLocal;
        }
              
        void ExpList::IProductWRTDerivBase(const int dir, 
                                           const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD, NekDouble> &outarray)
        {
            int    i;
            int    cnt  = 0;
            int    cnt1 = 0;
            
            Array<OneD,NekDouble> e_outarray;
            
            for(i = 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->IProductWRTDerivBase(dir,inarray+cnt,
                                             e_outarray = outarray+cnt1);
                cnt  += (*m_exp)[i]->GetTotPoints();
                cnt1 += (*m_exp)[i]->GetNcoeffs();
            }
            m_transState = eLocal;
        }
            
        void ExpList::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &out_d0, 
                                Array<OneD, NekDouble> &out_d1, 
                                Array<OneD, NekDouble> &out_d2)
        {
            int  cnt = 0;
            int  i;
            Array<OneD, NekDouble> e_out_d0;
            Array<OneD, NekDouble> e_out_d1;
            Array<OneD, NekDouble> e_out_d2;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                e_out_d0 = out_d0 + cnt;
                if(out_d1.num_elements())
                {
                    e_out_d1 = out_d1 + cnt;
                }
                
                if(out_d2.num_elements())
                {
                    e_out_d2 = out_d2 + cnt;
                }
                
                (*m_exp)[i]->PhysDeriv(inarray+cnt,e_out_d0,e_out_d1,e_out_d2);
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
        }



        void ExpList::PhysDeriv(const int dir, 
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD, NekDouble> &out_d)
        {
            int  cnt = 0;
            int  i;
            Array<OneD, NekDouble> e_out_d;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                e_out_d = out_d + cnt;
                (*m_exp)[i]->PhysDeriv(dir, inarray+cnt, e_out_d);
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
        }  

        void ExpList::MultiplyByElmtInvMass(const Array<OneD, const NekDouble> &inarray, 
                                            Array<OneD, NekDouble> &outarray)
        {
            GlobalMatrixKey mkey(StdRegions::eInvMass);
            const DNekScalBlkMatSharedPtr& InvMass = GetBlockMatrix(mkey);

            // Inverse mass matrix
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);            
            if(inarray.get() == outarray.get())
            {
                NekVector<const NekDouble> in(m_ncoeffs,inarray); // copy data
                out = (*InvMass)*in;
            }
            else
            {
                NekVector<const NekDouble> in(m_ncoeffs,inarray,eWrapper);
                out = (*InvMass)*in;
            }

        }

        void ExpList::FwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, 
                                          Array<OneD, NekDouble> &outarray)
        {
            Array<OneD,NekDouble> f(m_ncoeffs);

            IProductWRTBase_IterPerExp(inarray,f);
            MultiplyByElmtInvMass(f,outarray);

        }

        void ExpList::FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> &outarray)
        {
            int cnt  = 0;
            int cnt1 = 0;
            int i;
            
            Array<OneD,NekDouble> e_outarray;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->FwdTrans_BndConstrained(inarray+cnt, 
                                      e_outarray = outarray+cnt1);
                cnt  += (*m_exp)[i]->GetTotPoints();
                cnt1 += (*m_exp)[i]->GetNcoeffs();
            }
        }

        const DNekScalBlkMatSharedPtr ExpList::GenBlockMatrix(const GlobalMatrixKey &gkey)
        {
            int i,j,cnt1;
            int n_exp = GetExpSize();
            Array<OneD,unsigned int> nrows(n_exp);
            Array<OneD,unsigned int> ncols(n_exp);
            DNekScalMatSharedPtr    loc_mat;
            DNekScalBlkMatSharedPtr BlkMatrix;            

            switch(gkey.GetMatrixType())
            {      
            case StdRegions::eBwdTrans:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[i]->GetTotPoints();
                        ncols[i] = (*m_exp)[i]->GetNcoeffs();
                    }
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[i]->GetNcoeffs();
                        ncols[i] = (*m_exp)[i]->GetTotPoints();
                    }
                }
                break;                           
            case StdRegions::eMass: 
            case StdRegions::eInvMass:
            case StdRegions::eHelmholtz:
            case StdRegions::eLaplacian:
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[i]->GetNcoeffs();
                        ncols[i] = (*m_exp)[i]->GetNcoeffs();
                    }
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal, "Global Matrix creation not defined for this type of matrix");
                }
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            int nvarcoeffs = gkey.GetNvariableCoefficients();
            Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);
            
            for(i = cnt1 = 0; i < n_exp; ++i)
            {
                if(nvarcoeffs>0)
                {
                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs[j] = gkey.GetVariableCoefficient(j) + cnt1;
                    }
                    cnt1  += (*m_exp)[i]->GetTotPoints();
                }

                LocalRegions::MatrixKey matkey(gkey.GetMatrixType(),
                                               (*m_exp)[i]->DetExpansionType(),
                                               *(*m_exp)[i],
                                               gkey.GetConstants(),
                                               varcoeffs);

                loc_mat = (*m_exp)[i]->GetLocMatrix(matkey);
                BlkMatrix->SetBlock(i,i,loc_mat);
            }
            
            return BlkMatrix;
        }

        const DNekScalBlkMatSharedPtr& ExpList::GetBlockMatrix(const GlobalMatrixKey &gkey)
        {
            BlockMatrixMap::iterator matrixIter = m_blockMat->find(gkey);

            if(matrixIter == m_blockMat->end())
            {
                return ((*m_blockMat)[gkey] = GenBlockMatrix(gkey));
            }
            else
            {
                return matrixIter->second;
            }
        }
       
        void ExpList::GeneralMatrixOp(const GlobalLinSysKey &gkey,
                                      const Array<OneD, const NekDouble> &inarray,                     
                                      Array<OneD, NekDouble>    &outarray)
        {
            int  i,j;
            int  cnt  = 0;
            int  cnt1 = 0;
            Array<OneD,NekDouble>      e_outarray;

            int nvarcoeffs = gkey.GetNvariableCoefficients();
            Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                if(nvarcoeffs>0)
                {
                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs[j] = gkey.GetVariableCoefficient(j) + cnt1;
                    }
                    cnt1  += (*m_exp)[i]->GetTotPoints();
                }
                
                StdRegions::StdMatrixKey mkey(gkey.GetMatrixType(),
                                              (*m_exp)[i]->DetExpansionType(),
                                              *((*m_exp)[i]),
                                              gkey.GetConstants(),varcoeffs);
               
                (*m_exp)[i]->GeneralMatrixOp(inarray + cnt, 
                                             e_outarray = outarray+cnt,
                                             mkey);
                
                cnt   += (*m_exp)[i]->GetNcoeffs();
            }      
        }
        
	
	GlobalLinSysSharedPtr ExpList::GenGlobalLinSysFullDirect(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
	{
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            DNekMatSharedPtr      Gmat;
            
            Gmat = GenGlobalMatrixFull(mkey, locToGloMap);

            if(Gmat->GetRows())
            {
                PointerWrapper w = eWrapper;
                linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,w);
            }
            
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }


	void ExpList::GlobalEigenSystem(const DNekMatSharedPtr &Gmat, 
                                        Array<OneD, NekDouble> &EigValReal, 
                                        Array<OneD, NekDouble> &EigValImag, 
                                        Array<OneD, NekDouble> &EigVecs)
	{
            switch(Gmat->GetStorageType())
            {
            case eFULL:
                {
                    int lda = Gmat->GetRows(),info = 0;
                    NekDouble dum, lwork = 3*lda;
                    Array<OneD,NekDouble> work(3*lda);

                    if(EigVecs == NullNekDouble1DArray) // calculate Right Eigen Vectors
                    {
                        Lapack::Dgeev('N','V',lda, Gmat->GetRawPtr(),lda,&EigValReal[0],
                                      &EigValImag[0],&dum,1,&EigVecs[0],lda,&work[0],lwork,info);
                    }
                    else
                    {
                        Lapack::Dgeev('N','N',lda, Gmat->GetRawPtr(),lda,&EigValReal[0],
                                      &EigValImag[0],&dum,1,&dum,1,&work[0],lwork,info);
                    }
                    ASSERTL0(info == 0,"Info is not zero");
                }
                break;
            default:
                ASSERTL0(false,"Eigenvalue evaluation is not set up for this matrix storage type");
                break;
            }
        }


        DNekMatSharedPtr ExpList::GenGlobalMatrixFull(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
        {
            int i,j,n,gid1,gid2,loc_lda,cnt,cnt1;
            NekDouble sign1,sign2,value;;
            DNekScalMatSharedPtr loc_mat;

            int totDofs     = locToGloMap->GetNumGlobalCoeffs();
            int NumDirBCs   = locToGloMap->GetNumGlobalDirBndCoeffs();

            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0;

            DNekMatSharedPtr Gmat;
            int bwidth = locToGloMap->GetFullSystemBandWidth();

            int nvarcoeffs = mkey.GetNvariableCoefficients();
            Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);
            MatrixStorage matStorage;

            switch(mkey.GetMatrixType())
            {
                // case for all symmetric matices
            case StdRegions::eHelmholtz: 
            case StdRegions::eLaplacian:
                if( (2*(bwidth+1)) < rows)
                {
                    matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage,bwidth,bwidth);
                }
                else
                {
                    matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
                }  
                
                break;
            default: // Assume general matrix - currently only set up for full invert
                {
                    matStorage = eFULL;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
            
                }       
            }             

            
            // fill global symmetric matrix 
            for(n = cnt = cnt1 = 0; n < (*m_exp).size(); ++n)
            {
                if(nvarcoeffs>0)
                {
                        for(j = 0; j < nvarcoeffs; j++)
                        {
                            varcoeffs[j] = mkey.GetVariableCoefficient(j) + cnt1;
                        }
                        cnt1  += (*m_exp)[n]->GetTotPoints();
                }
                
                LocalRegions::MatrixKey matkey(mkey.GetMatrixType(),
                                               (*m_exp)[n]->DetExpansionType(),
                                               *(*m_exp)[n],
                                               mkey.GetConstants(),
                                               varcoeffs);
                
                loc_mat = (*m_exp)[n]->GetLocMatrix(matkey);               
                loc_lda = loc_mat->GetColumns();
                
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetLocalToGlobalMap(cnt + i) - NumDirBCs;
                    sign1 =  locToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetLocalToGlobalMap(cnt + j) - NumDirBCs;
                            sign2 = locToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                // When global matrix is symmetric,
                                // only add the value for the upper
                                // triangular part in order to avoid
                                // entries to be entered twice
                                if((matStorage == eFULL)||(gid2 >= gid1))
                                {
                                    value = Gmat->GetValue(gid1,gid2) + sign1*sign2*(*loc_mat)(i,j);
                                    Gmat->SetValue(gid1,gid2,value);
                                }
                            }		
                        }
                    }
                }
                cnt   += (*m_exp)[n]->GetNcoeffs();
            }        

//             for(i = NumDirBCs; i < NumDirBCs+NumRobinBCs; ++i)
//             {
//                 // Find a way to deal with second parameter of the Robin BC
//                 NekDouble b=1.0;
//                 (*Gmat)((locToGloMap->GetBndCondGlobalID())[i]-NumDirBCs,
//                         (locToGloMap->GetBndCondGlobalID())[i]-NumDirBCs)
//                     -= mkey.GetScaleFactor() * b;
//             }
            

            return Gmat;

            }

        GlobalLinSysSharedPtr ExpList::GenGlobalLinSysStaticCond(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt,cnt1;
            NekDouble sign1,sign2,value;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            
            int nBndDofs = locToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = locToGloMap->GetNumGlobalDirBndCoeffs();

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;

            DNekMatSharedPtr Gmat;
            int bwidth = locToGloMap->GetBndSystemBandWidth();
            if( (2*(bwidth+1)) < rows)
            {
                MatrixStorage matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage,bwidth,bwidth);
            }
            else
            {
                MatrixStorage matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
            }                

            // Setup Block Matrix systems
            int n_exp = GetExpSize();
            Array<OneD,unsigned int> nbdry_size(n_exp);
            Array<OneD,unsigned int> nint_size(n_exp);
            DNekScalBlkMatSharedPtr BinvD;
            DNekScalBlkMatSharedPtr invD;
            DNekScalBlkMatSharedPtr C;

            // set up an array of integers for block matrix construction
            for(i = 0; i < n_exp; ++i)
            {
                nbdry_size[i] = (*m_exp)[i]->NumBndryCoeffs();
                nint_size[i]  = (*m_exp)[i]->GetNcoeffs() - (*m_exp)[i]->NumBndryCoeffs();
            }
            
            MatrixStorage blkmatStorage = eDIAGONAL;
            BinvD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size,nint_size,blkmatStorage);
            invD  = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size,nint_size, blkmatStorage);
            C     = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size,nbdry_size,blkmatStorage);

            DNekScalMatSharedPtr tmp_mat; 

            int nvarcoeffs = mkey.GetNvariableCoefficients();
            Array<OneD, Array<OneD,NekDouble> > varcoeffs(nvarcoeffs);

            // fill global matrix 
            for(n = cnt = cnt1 = 0; n < (*m_exp).size(); ++n)
            {
                if(nvarcoeffs>0)
                {
                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs[j] = mkey.GetVariableCoefficient(j) + cnt1;
                    }
                    cnt1  += (*m_exp)[n]->GetTotPoints();
                }

                LocalRegions::MatrixKey matkey(mkey.GetMatrixType(),
                                               (*m_exp)[n]->DetExpansionType(),
                                               *(*m_exp)[n],
                                               mkey.GetConstants(),
                                               varcoeffs);

                loc_mat = (*m_exp)[n]->GetLocStaticCondMatrix(matkey);                   
                loc_lda = (*m_exp)[n]->NumBndryCoeffs(); 

                BinvD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                invD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                C->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                
                // Set up  Matrix; 
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetLocalToGlobalBndMap(cnt + i) - NumDirBCs;
                    sign1 = locToGloMap->GetLocalToGlobalBndSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetLocalToGlobalBndMap(cnt + j) - NumDirBCs;
                            sign2 = locToGloMap->GetLocalToGlobalBndSign(cnt + j);
                            if(gid2 >= 0)
                            {
                                // As the global matrix should be
                                // symmetric, only add the value for
                                // the upper triangular part in order
                                // to avoid entries to be entered
                                // twice
                                if(gid2 >= gid1)
                                {
                                    value = Gmat->GetValue(gid1,gid2) + sign1*sign2*(*loc_mat)(i,j);
                                    Gmat->SetValue(gid1,gid2,value);
                                }
                            }
                        }		
                    }
                }
                cnt += (*m_exp)[n]->NumBndryCoeffs();
            }
            
            if(rows)
            {
                PointerWrapper w = eWrapper;
                linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,w);
            }
            
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys,BinvD,C,invD);
            return returnlinsys;
        }


	GlobalLinSysSharedPtr ExpList::GenGlobalLinSys(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
	{
            GlobalLinSysSharedPtr returnlinsys; 

            switch(mkey.GetGlobalSysSolnType())
            {
            case eDirectFullMatrix:
                returnlinsys = GenGlobalLinSysFullDirect(mkey, locToGloMap);
                break;
            case eDirectStaticCond:
                returnlinsys = GenGlobalLinSysStaticCond(mkey, locToGloMap);
                break;
            default:
                ASSERTL0(false,"Matrix solution type not defined");
                break;
            }
            
            return returnlinsys;
        }

        GlobalLinSysSharedPtr ExpList::GenGlobalBndLinSys(const GlobalLinSysKey     &mkey, const LocalToGlobalBaseMap &LocToGloBaseMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            
            int totDofs       = LocToGloBaseMap.GetNumGlobalBndCoeffs();
            int NumDirBCs     = LocToGloBaseMap.GetNumLocalDirBndCoeffs();
            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero    = 0.0,sign1,sign2,value; 
            NekDouble factor1 = mkey.GetConstant(0);
            NekDouble factor2 = mkey.GetConstant(1);
            StdRegions::MatrixType linsystype = mkey.GetMatrixType();

            DNekMatSharedPtr Gmat;
            int bwidth = LocToGloBaseMap.GetBndSystemBandWidth();

            if( (2*(bwidth+1)) < rows)
            {
                MatrixStorage matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage,bwidth,bwidth);
            }
            else
            {
                MatrixStorage matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
            }                

            ASSERTL0(linsystype == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");
            
            // fill global matrix 
            for(n = cnt = 0; n < (*m_exp).size(); ++n)
            {
                // Matrix to Bnd Sys
                LocalRegions::MatrixKey Umatkey(linsystype, (*m_exp)[n]->DetExpansionType(),*((*m_exp)[n]), factor1,factor2);
                DNekScalMat &BndSys = *((*m_exp)[n]->GetLocMatrix(Umatkey)); 
                

                loc_lda = BndSys.GetColumns();
                
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = LocToGloBaseMap.GetLocalToGlobalBndMap (cnt + i) 
                        - NumDirBCs;
                    sign1 = LocToGloBaseMap.GetLocalToGlobalBndSign(cnt + i); 

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = LocToGloBaseMap.GetLocalToGlobalBndMap(cnt+j) - NumDirBCs;
                            sign2 = LocToGloBaseMap.GetLocalToGlobalBndSign(cnt+j); 
                            if(gid2 >= 0)
                            {
                                if(gid2 >= gid1)
                                {
                                    value = (*Gmat)(gid1,gid2) 
                                        + sign1*sign2*(BndSys)(i,j);
                                    
                                    Gmat->SetValue(gid1, gid2, value);
                                }
                            }
                        }		
                    }
                }
                cnt += loc_lda;
            }
            
            PointerWrapper w = eWrapper;
            if(rows)
            {
                linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,w);
            }

            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }

        void ExpList::BwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD, NekDouble> &outarray)
        {
            int  i;
            int  cnt  = 0;
            int  cnt1 = 0;
            Array<OneD,NekDouble> e_outarray;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->BwdTrans(inarray + cnt, 
                                      e_outarray = outarray+cnt1);
                cnt   += (*m_exp)[i]->GetNcoeffs();
                cnt1  += (*m_exp)[i]->GetTotPoints();
            }        
        }
        
        void ExpList::GetCoords(Array<OneD, NekDouble> &coord_0,
                                Array<OneD, NekDouble> &coord_1,
                                Array<OneD, NekDouble> &coord_2)
        {
            int    i, cnt = 0;
            Array<OneD, NekDouble> e_coord_0;
            Array<OneD, NekDouble> e_coord_1;
            Array<OneD, NekDouble> e_coord_2;
            
            switch(GetExp(0)->GetCoordim())
            {
            case 1:
                for(i= 0; i < GetExpSize(); ++i)
                {
                    e_coord_0 = coord_0 + cnt;
                    (*m_exp)[i]->GetCoords(e_coord_0);
                    cnt  += (*m_exp)[i]->GetTotPoints();
                }
                break;
            case 2: 
                ASSERTL0(coord_1.num_elements() != 0, 
                         "output coord_1 is not defined");
                
                for(i= 0; i < GetExpSize(); ++i)
                {
                    e_coord_0 = coord_0 + cnt;
                    e_coord_1 = coord_1 + cnt;
                    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1);
                    cnt  += (*m_exp)[i]->GetTotPoints();
                }
                break;
            case 3: 
                ASSERTL0(coord_1.num_elements() != 0,
                         "output coord_1 is not defined");
                ASSERTL0(coord_2.num_elements() != 0,
                         "output coord_2 is not defined");
                
                for(i= 0; i < GetExpSize(); ++i)
                {
                    e_coord_0 = coord_0 + cnt;
                    e_coord_1 = coord_1 + cnt;
                    e_coord_2 = coord_2 + cnt;
                    
                    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1,e_coord_2);
                    cnt  += (*m_exp)[i]->GetTotPoints();
                }
                break;
            }
        }
    
        void ExpList::WriteToFile(std::ofstream &out, OutputFormat format)
        {  
            if(format==eTecplot)
            {
                int i,cnt = 0;

                Array<OneD, const NekDouble> phys = m_phys;
                
                if(m_physState == false)
                {
                    BwdTrans(m_coeffs,m_phys);
                }
                
                (*m_exp)[0]->SetPhys(phys);
                (*m_exp)[0]->WriteToFile(out,eTecplot,1);
                cnt  += (*m_exp)[0]->GetTotPoints();
                
                for(i= 1; i < GetExpSize(); ++i)
                {
                    (*m_exp)[i]->SetPhys(phys+cnt);
                    (*m_exp)[i]->WriteToFile(out,eTecplot,false); 
                    cnt  += (*m_exp)[i]->GetTotPoints();
                }
            }
            else if(format==eGmsh)
            {   
               
                out<<"View.MaxRecursionLevel = 8;"<<endl;
                out<<"View.TargetError = 0.00;"<<endl;
                 
                int i,j,k;
                int nElementalCoeffs =  (*m_exp)[0]->GetBasisNumModes(0);
                StdRegions::ExpansionType locShape = (*m_exp)[0]->DetExpansionType();

                int nDumpCoeffs =  nElementalCoeffs*nElementalCoeffs;
                Array<TwoD, int> exponentMap(nDumpCoeffs,3,0);
                int cnt = 0;
                for(i = 0; i < nElementalCoeffs; i++)
                {
                    for(j = 0; j < nElementalCoeffs; j++)
                    {
                        exponentMap[cnt][0] = j;
                        exponentMap[cnt++][1] = i;
                    }         
                }

                PutCoeffsInToElmtExp();
                bool dumpNewView = true;
                bool closeView = false;
                for(i= 0; i < GetExpSize(); ++i)
                {                 
                    if(nElementalCoeffs != (*m_exp)[i]->GetBasisNumModes(0))
                    {
                        ASSERTL0(false,"Not all elements have the same number of expansions, this will"
                                 "probably lead to a corrupt Gmsh-output file")
                    }

                    if(i>0)
                    {
                        if ( ((*m_exp)[i]->DetExpansionType())!=((*m_exp)[i-1]->DetExpansionType()) )
                        {
                            dumpNewView = true;
                        }
                        else
                        {
                            dumpNewView = false;
                        }
                    }
                    if(i<GetExpSize()-1)
                    {
                        if ( ((*m_exp)[i]->DetExpansionType())!=((*m_exp)[i+1]->DetExpansionType()) )
                        {
                            closeView = true;
                        }
                        else
                        {
                            closeView = false;
                        }
                    }
                    else
                    {
                            closeView = true;
                    }

                    if(dumpNewView)
                    {
                        out<<"View \" \" {"<<endl;
                    }

                    (*m_exp)[i]->WriteToFile(out,eGmsh,false); 

                    if(closeView)
                    {
                        out<<"INTERPOLATION_SCHEME"<<endl;
                        out<<"{"<<endl;
                        for(k=0; k < nDumpCoeffs; k++)
                        {
                            out<<"{";
                            for(j = 0; j < nDumpCoeffs; j++)
                            {
                                if(k==j)
                                {
                                    out<<"1.00";
                                }
                                else
                                {
                                    out<<"0.00";
                                }
                                if(j < nDumpCoeffs - 1)
                                {
                                    out<<", ";
                                }
                            }
                            if(k < nDumpCoeffs - 1)
                            {
                                out<<"},"<<endl;
                            }
                            else
                            {
                                out<<"}"<<endl<<"}"<<endl;
                            }
                        }
                        
                        out<<"{"<<endl;
                        for(k=0; k < nDumpCoeffs; k++)
                        {
                            out<<"{";
                            for(j = 0; j < 3; j++)
                            {
                                out<<exponentMap[k][j];
                                if(j < 2)
                                {
                                    out<<", ";
                                }
                            }
                            if(k < nDumpCoeffs - 1)
                            {
                                out<<"},"<<endl;
                            }
                            else
                            {
                                out<<"}"<<endl<<"};"<<endl;
                            }
                        }
                        out<<"};"<<endl;
                    }                 
                }       
                out<<"Combine ElementsFromAllViews;"<<endl;
                out<<"View.Name = \"\";"<<endl;
            }                    
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }
    
        NekDouble  ExpList::Linf(const Array<OneD, const NekDouble> &soln)
        {            
            NekDouble err = 0.0;
            int       i,cnt = 0;
            Array<OneD, const NekDouble> phys = m_phys;

            for(i= 0; i < GetExpSize(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(phys+cnt);
                err  = std::max(err,(*m_exp)[i]->Linf(soln + cnt));
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
            
            return err;
        }
        
        NekDouble ExpList::L2(const Array<OneD, const NekDouble> &soln)
        {
            
            NekDouble err = 0.0,errl2;
            int    i,cnt = 0;
            Array<OneD, const NekDouble> phys = m_phys;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(phys+cnt);
                errl2 = (*m_exp)[i]->L2(soln+cnt);
                err += errl2*errl2;
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
            
            return sqrt(err);
        }
    
    } //end of namespace
} //end of namespace

