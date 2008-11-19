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
            m_physState(false)
        {            
            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
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
        
        void ExpList::IProductWRTBase(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetPhysState() == true,
                     "Physical space is not set to true ");
            
            IProductWRTBase(Sin.GetPhys(),m_coeffs);
            m_physState = false;
        }
        
        void ExpList::IProductWRTBase(const Array<OneD, const NekDouble> &inarray, 
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
        

        void ExpList::IProductWRTDerivBase(const int dir, const ExpList &Sin)
        {
            ASSERTL2(Sin.GetPhysState() == true,
                     "Physical space is not set to true ");
            
            IProductWRTDerivBase(dir,Sin.GetPhys(),m_coeffs);
            m_physState = false;
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
        
        
        void ExpList::PhysDeriv(ExpList &out_d0, 
                                ExpList &out_d1, 
                                ExpList &out_d2)
        {
            ASSERTL2(m_physState == true,
                     "local physical space is not true ");
            PhysDeriv(m_phys,
                      out_d0.UpdatePhys(), 
                      out_d1.UpdatePhys(), 
                      out_d2.UpdatePhys());
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

        void ExpList::MultiplyByElmtInvMass(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetTransState() == eLocal ||
                     Sin.GetTransState() == eLocalCont, 
                     "Error input state is not in transformed space");
            
            MultiplyByElmtInvMass(Sin.GetPhys(),m_coeffs);
            m_transState = eLocal;
        }


        void ExpList::FwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetPhysState() == true,
                     "Sin physical space is not true ");
            
            FwdTrans(Sin.GetPhys(),m_coeffs);
            m_transState = eLocal;
        }

        void ExpList::FwdTrans_BndConstrained(const ExpList &Sin)
        {
            ASSERTL1(Sin.GetPhysState() == true,
                     "Sin physical space is not true ");
            
            FwdTrans_BndConstrained(Sin.GetPhys(),m_coeffs);
            m_transState = eLocal;
        }
        

        void ExpList::MultiplyByElmtInvMass(const Array<OneD, const NekDouble> &inarray, 
                                            Array<OneD, NekDouble> &outarray)
        {
            static DNekScalBlkMatSharedPtr InvMass;
            if(!InvMass.get())
            {
                InvMass = SetupBlockMatrix(StdRegions::eInvMass);
            }
            
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

        void ExpList::FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                               Array<OneD, NekDouble> &outarray)
        {
            Array<OneD,NekDouble> f(m_ncoeffs);

            IProductWRTBase(inarray,f);
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
        
        DNekScalBlkMatSharedPtr  ExpList::SetupBlockMatrix(StdRegions::MatrixType mtype, NekDouble scalar, NekDouble constant)
        {
            int i;
            int n_exp = GetExpSize();
            Array<OneD,unsigned int> exp_size(n_exp);
            DNekScalMatSharedPtr loc_mat;
            DNekScalBlkMatSharedPtr BlkMatrix;
            
            // set up an array of integers for block matrix construction
            for(i = 0; i < n_exp; ++i)
            {
                exp_size[i] = (*m_exp)[i]->GetNcoeffs();
            }

            BlkMatrix = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size);
            
            for(i = 0; i < n_exp; ++i)
            {
                LocalRegions::MatrixKey mkey(mtype,(*m_exp)[i]->DetExpansionType(),*((*m_exp)[i]),scalar,constant);
                loc_mat = (*m_exp)[i]->GetLocMatrix(mkey);
                BlkMatrix->SetBlock(i,i,loc_mat);
            }
            
            return BlkMatrix;
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
                
                StdRegions::StdMatrixKey mkey(gkey.GetLinSysType(),(*m_exp)[i]->DetExpansionType(),
                                              *((*m_exp)[i]),gkey.GetConstants(),varcoeffs);
               
                (*m_exp)[i]->GeneralMatrixOp(inarray + cnt, 
                                             e_outarray = outarray+cnt,
                                             mkey);
                
                cnt   += (*m_exp)[i]->GetNcoeffs();
            }      
        }
        
	
	GlobalLinSysSharedPtr ExpList::GenGlobalLinSysFullDirect(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt,cnt1;
            NekDouble sign1,sign2;
            DNekScalMatSharedPtr loc_mat;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;

            int totDofs     = locToGloMap->GetNumGlobalCoeffs();
            int NumDirBCs   = locToGloMap->GetNumDirichletBndCoeffs();

            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0;
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);

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

                LocalRegions::MatrixKey matkey(mkey.GetLinSysType(),
                                               (*m_exp)[n]->DetExpansionType(),
                                               *(*m_exp)[n],
                                               mkey.GetConstants(),
                                               varcoeffs);
                
                loc_mat = (*m_exp)[n]->GetLocMatrix(matkey);               
                loc_lda = loc_mat->GetColumns();

                //cout<<"loc mat "<<endl<<*loc_mat<<endl;
		    
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetLocalToGlobalMap(cnt + i);
                    sign1 =  locToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetLocalToGlobalMap(cnt + j);
                            sign2 = locToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += sign1*sign2*(*loc_mat)(i,j);
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

            //linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }


	GlobalLinSysSharedPtr ExpList::GenGlobalLinSysStaticCond(const GlobalLinSysKey &mkey, const LocalToGlobalC0ContMapSharedPtr &locToGloMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt,cnt1;
            NekDouble sign1,sign2;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            
            int nBndDofs = locToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = locToGloMap->GetNumDirichletBndCoeffs();

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);
            
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
            
            BinvD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size,nint_size);
            invD  = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size,nint_size);
            C     = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size,nbdry_size);

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

                LocalRegions::MatrixKey matkey(mkey.GetLinSysType(),
                                               (*m_exp)[n]->DetExpansionType(),
                                               *(*m_exp)[n],
                                               mkey.GetConstants(),
                                               varcoeffs);

                loc_mat = (*m_exp)[n]->GetLocStaticCondMatrix(matkey);                   
                loc_lda = (*m_exp)[n]->NumBndryCoeffs(); 

                //cout<<"loc mat "<<endl<<*loc_mat<<endl;
		    

                BinvD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                invD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                C->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                
                // Set up  Matrix; 
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetLocalToGlobalBndMap(cnt + i);
                    sign1 = locToGloMap->GetLocalToGlobalBndSign(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetLocalToGlobalBndMap(cnt + j);
                            sign2 = locToGloMap->GetLocalToGlobalBndSign(cnt + j);
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += sign1*sign2*(*loc_mat)(i,j);
                            }
                        }		
                    }
                }
                cnt += (*m_exp)[n]->NumBndryCoeffs();
            }
            
            // Believe that we need a call of the type:
            //linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            if(rows)
            {
                linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
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
            int NumDirBCs     = LocToGloBaseMap.GetNumDirichletBndCoeffs();
            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero    = 0.0,sign1,sign2; 
            NekDouble factor1 = mkey.GetConstant(0);
            NekDouble factor2 = mkey.GetConstant(1);
            StdRegions::MatrixType linsystype = mkey.GetLinSysType();
            
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero); 
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
                    gid1  = LocToGloBaseMap.GetLocalToGlobalBndMap(cnt + i);
                    sign1 = LocToGloBaseMap.GetLocalToGlobalBndSign(cnt + i); 

                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = LocToGloBaseMap.GetLocalToGlobalBndMap(cnt + j);
                            sign2 = LocToGloBaseMap.GetLocalToGlobalBndSign(cnt + j); 
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += sign1*sign2*(BndSys)(i,j);
                            }
                        }		
                    }
                }
                cnt += loc_lda;
            }
            
            // Believe that we need a call of the type:
            // linsys=MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            linsys       = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }


        void ExpList::BwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetTransState() == eLocal ||
                     Sin.GetTransState() == eLocalCont, 
                     "Error input state is not in transformed space");
            
            BwdTrans(Sin.GetCoeffs(),m_phys);
            m_physState = true;
        }
        
        void ExpList::BwdTrans(const Array<OneD, const NekDouble> &inarray,
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
                    BwdTrans(*this);
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
    
        NekDouble  ExpList::Linf(const ExpList &Sol)
        {
            ASSERTL2(Sol.GetPhysState() == true,
                     "local physical space is not true ");
            
            NekDouble err = 0.0;
            int       i,cnt = 0;
            Array<OneD, const NekDouble> soln = Sol.GetPhys();
            Array<OneD, const NekDouble> phys = m_phys;
            
            if(m_physState == false)
            {
                BwdTrans(*this);
            }

            for(i= 0; i < GetExpSize(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(phys+cnt);
                err  = std::max(err,(*m_exp)[i]->Linf(soln + cnt));
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
            
            return err;
        }
    
        NekDouble  ExpList::L2(const ExpList &Sol)
        {
            ASSERTL2(Sol.GetPhysState() == true,
                     "local physical space is not true ");

            if(m_physState == false)
            {
                BwdTrans(*this);
            }
            
            return L2(Sol.GetPhys());
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

