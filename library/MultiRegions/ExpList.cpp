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

namespace Nektar
{
    namespace MultiRegions
    {

    ExpList::ExpList(void):
        m_ncoeffs(0),
        m_npoints(0),
        m_transState(eNotSet),
        m_physState(false)
    {
            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
    }
    

        ExpList::ExpList(const ExpList &in):
            m_ncoeffs(in.m_ncoeffs),
            m_npoints(in.m_npoints),
            m_transState(eNotSet),
            m_physState(false),
            m_exp(in.m_exp)
        {
        }

    ExpList::~ExpList()
    {
    }
    
    
    /** \brief Integrate the physical point list \a inarray over region
        and return the value
        
        Inputs:\n
        
        - \a inarray: definition of function to be returned at quadrature point 
        of expansion. 
        
        Outputs:\n
        
        - returns \f$ \sum_{i=1}^{n_{el}} \int_{\Omega_i} u(\xi_1)d \xi_1 \f$ 
    */
        NekDouble ExpList::PhysIntegral(void)
        {
            ASSERTL2(m_physState == true,
                     "local physical space is not true ");

            return PhysIntegral(m_phys);
        }
        
    NekDouble ExpList::PhysIntegral(const ConstArray<OneD, NekDouble> &inarray)
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

    void ExpList::IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray, 
                      Array<OneD, NekDouble> &outarray)
    {
        int    i;
        int    cnt  = 0;
        int    cnt1 = 0;

        ConstArray<OneD,NekDouble> e_inarray;
        Array<OneD,NekDouble> e_outarray;
        
        for(i = 0; i < GetExpSize(); ++i)
        {
                (*m_exp)[i]->IProductWRTBase(e_inarray = inarray+cnt,
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

    
        void ExpList::PhysDeriv(const ConstArray<OneD, NekDouble> &inarray,
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
                Array<OneD, NekDouble> myout1((*m_exp)[i]->GetNcoeffs());
                Array<OneD, NekDouble> myout2((*m_exp)[i]->GetNcoeffs());

                e_out_d0 = out_d0 + cnt;
                if(out_d1.num_elements())
                {
                    e_out_d1 = out_d1 + cnt;
                }
                
                if(out_d2.num_elements())
                {
                    e_out_d2 = out_d2 + cnt;
                }
                
                (*m_exp)[i]->PhysDeriv(inarray,e_out_d0,e_out_d1,e_out_d2);
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
        }

        void ExpList::FwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetPhysState() == true,
                     "Sin physical space is not true ");
            
            FwdTrans(Sin.GetPhys(),m_coeffs);
            m_transState = eLocal;
        }
        
        void ExpList::FwdTrans(const ConstArray<OneD, NekDouble> &inarray, 
                               Array<OneD, NekDouble> &outarray)
        {
#if 0  // elemental matrix inverse
            int cnt  = 0;
            int cnt1 = 0;
            int i;
            ConstArray<OneD,NekDouble> e_inarray;
            Array<OneD,NekDouble> e_outarray;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->FwdTrans(e_inarray  = inarray+cnt, 
                                      e_outarray = outarray+cnt1);
                cnt  += (*m_exp)[i]->GetTotPoints();
                cnt1 += (*m_exp)[i]->GetNcoeffs();
            }
#else  // block matrix inverse 
            static DNekScalBlkMatSharedPtr InvMass;
            Array<OneD,NekDouble> f(m_ncoeffs);
            if(!InvMass.get())
            {
                InvMass = SetupBlockMatrix(StdRegions::eInvMass);
            }
            
            // Inner product
            IProductWRTBase(inarray,f);

            // Inverse mass matrix
            DNekVec in (m_ncoeffs,f);
            DNekVec out(m_ncoeffs,outarray,eWrapper);            
            out = (*InvMass)*in;
#endif
        }
        
        DNekScalBlkMatSharedPtr  ExpList::SetupBlockMatrix(StdRegions::MatrixType mtype, NekDouble scalar, NekDouble constant)
        {
            int i;
            int n_exp = GetExpSize();
            //Array<OneD,unsigned int> exp_size(n_exp);
            Array<OneD,unsigned int> exp_size(n_exp);
            DNekScalMatSharedPtr loc_mat;
            DNekScalBlkMatSharedPtr BlkMatrix;
            
            // set up an array of integers for block matrix construction
            for(i = 0; i < n_exp; ++i)
            {
                exp_size[i] = (*m_exp)[i]->GetNcoeffs();
            }

            unsigned int *exp_size_2 = new unsigned int[n_exp];
            exp_size_2 = &exp_size[0];
            BlkMatrix = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,exp_size_2,exp_size_2);
            
            //BlkMatrix = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,exp_size,exp_size);
            //Cannot get this call to work with array of integers
            
            for(i = 0; i < n_exp; ++i)
            {
                LocalRegions::MatrixKey mkey(mtype,(*m_exp)[i]->DetShapeType(),*((*m_exp)[i]),scalar,constant);
                loc_mat = (*m_exp)[i]->GetLocMatrix(mkey);
                
                // cout << loc_mat->GetOwnedMatrix() << endl;
                BlkMatrix->SetBlock(i,i,loc_mat);
            }
            
            return BlkMatrix;
        }
        
        
        void ExpList::GeneralMatrixOp(const StdRegions::MatrixType mtype,
                                      const ConstArray<OneD,NekDouble> &inarray,                     
                                      Array<OneD, NekDouble>    &outarray,
                                      NekDouble lambda)
        {
            int  i;
            int  cnt  = 0;
            int  cnt1 = 0;
            ConstArray<OneD,NekDouble> e_inarray;
            Array<OneD,NekDouble>      e_outarray;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->GeneralMatrixOp(mtype, e_inarray = inarray + cnt, 
                                             e_outarray = outarray+cnt,lambda);
                cnt   += (*m_exp)[i]->GetNcoeffs();
            }        
        }
        
	
	GlobalLinSysSharedPtr ExpList::GenGlobalLinSysFullDirect(const GlobalLinSysKey &mkey, LocalToGlobalMapSharedPtr &locToGloMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            NekDouble sign1,sign2;
            DNekScalMatSharedPtr loc_mat;
            StdRegions::StdExpansionVectorIter def;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;

            int totDofs     = locToGloMap->GetTotGloDofs();
            int NumDirBCs   = locToGloMap->GetNumDirichletDofs();

            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0;
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);
            
            // fill global matrix 
            for(n = cnt = 0; n < (*m_exp).size(); ++n)
            {
                LocalRegions::MatrixKey matkey(mkey.GetLinSysType(),
                                          (*m_exp)[n]->DetShapeType(),
                                          *(*m_exp)[n],mkey.GetFactor1());
                
                loc_mat = (*m_exp)[n]->GetLocMatrix(matkey);
                loc_lda = loc_mat->GetColumns();

                //cout<<"loc mat "<<endl<<*loc_mat<<endl;
		    
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetMap(cnt + i);
                    sign1 =  locToGloMap->GetSign(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetMap(cnt + j);
                            sign2 = locToGloMap->GetSign(cnt + j);
                            if(gid2 >= NumDirBCs)
                            {
                                (*Gmat)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                    += sign1*sign2*(*loc_mat)(i,j);
                            }
                        }		
                    }
                }
                cnt += (*m_exp)[n]->GetNcoeffs();
            }            
            
//             for(i = NumDirBCs; i < NumDirBCs+NumRobinBCs; ++i)
//             {
//                 // Find a way to deal with second parameter of the Robin BC
//                 NekDouble b=1.0;
//                 (*Gmat)((locToGloMap->GetBndCondGlobalID())[i]-NumDirBCs,
//                         (locToGloMap->GetBndCondGlobalID())[i]-NumDirBCs)
//                     -= mkey.GetScaleFactor() * b;
//             }

            // Believe that we need a call of the type:
            //linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,eWrapper);
            linsys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat);
            
            returnlinsys = MemoryManager<GlobalLinSys>::AllocateSharedPtr(mkey,linsys);
            return returnlinsys;
        }


	GlobalLinSysSharedPtr ExpList::GenGlobalLinSysStaticCond(const GlobalLinSysKey &mkey, LocalToGlobalMapSharedPtr &locToGloMap)
	{
            int i,j,n,gid1,gid2,loc_lda,cnt;
            NekDouble sign1,sign2;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekLinSysSharedPtr   linsys;
            GlobalLinSysSharedPtr returnlinsys;
            
            int nBndDofs = locToGloMap->GetTotGloBndDofs();
            int NumDirBCs = locToGloMap->GetNumDirichletDofs();

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;
            NekDouble zero = 0.0;
            DNekMatSharedPtr Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero);
            
            // Setup Block Matrix systems
            int n_exp = GetExpSize();
            //Array<OneD,unsigned int> nbdry_size(n_exp);
            //Array<OneD,unsigned int> nint_size(n_exp);
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

            unsigned int *nbdry_size_2 = new unsigned int[n_exp];
            unsigned int *nint_size_2  = new unsigned int[n_exp];
            nbdry_size_2 = &nbdry_size[0];
            nint_size_2  = &nint_size[0];
            
            BinvD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,nbdry_size_2,nint_size_2);
            invD  = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,nint_size_2,nint_size_2);
            C     = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,nint_size_2,nbdry_size_2);

            // needs to be formed as: 
            //BinvD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,bndry_size,nint_size);
            //invDC = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,nint_size,nint_size]);
            //invD = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(n_exp,n_exp,nint_size,nint_bndry);

            DNekScalMatSharedPtr tmp_mat; 

            // fill global matrix 
            for(n = cnt = 0; n < (*m_exp).size(); ++n)
            {
                LocalRegions::MatrixKey matkey(mkey.GetLinSysType(),
                                          (*m_exp)[n]->DetShapeType(),
                                          *(*m_exp)[n],mkey.GetFactor1());

                loc_mat = (*m_exp)[n]->GetLocStaticCondMatrix(matkey);
                loc_lda = (*m_exp)[n]->NumBndryCoeffs(); 		    
                
                BinvD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                invD->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
                C->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                
                // Set up interior Matrix; 
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetBndMap(cnt + i);
                    sign1 = locToGloMap->GetBndSign(cnt + i);
                    if(gid1 >= NumDirBCs)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetBndMap(cnt + j);
                            sign2 = locToGloMap->GetBndSign(cnt + j);
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


	GlobalLinSysSharedPtr ExpList::GenGlobalLinSys(const GlobalLinSysKey &mkey, LocalToGlobalMapSharedPtr &locToGloMap)
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

        void ExpList::BwdTrans(const ExpList &Sin)
        {
            ASSERTL2(Sin.GetTransState() == eLocal ||
                     Sin.GetTransState() == eLocalCont, 
                     "Error input state is not in transformed space");
            
            BwdTrans(Sin.GetCoeffs(),m_phys);
            m_physState = true;
        }
        
        void ExpList::BwdTrans(const ConstArray<OneD, NekDouble> &inarray,
                               Array<OneD, NekDouble> &outarray)
        {
            int  i;
            int  cnt  = 0;
            int  cnt1 = 0;
            ConstArray<OneD,NekDouble> e_inarray;
            Array<OneD,NekDouble> e_outarray;
            
            for(i= 0; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->BwdTrans(e_inarray = inarray + cnt, 
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
    
        void ExpList::WriteToFile(std::ofstream &out)
        {
            int i,cnt = 0;
            std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
            StdRegions::StdExpansionVectorIter def;
            ConstArray<OneD, NekDouble> phys = m_phys;
            
            if(m_physState == false)
            {
                BwdTrans(*this);
            }
            
            (*m_exp)[0]->SetPhys(phys);
            (*m_exp)[0]->WriteToFile(out,1);
            cnt  += (*m_exp)[0]->GetTotPoints();
            
            for(i= 1; i < GetExpSize(); ++i)
            {
                (*m_exp)[i]->SetPhys(phys+cnt);
                (*m_exp)[i]->WriteToFile(out,0); 
                cnt  += (*m_exp)[i]->GetTotPoints();
            }
        }
    
        NekDouble  ExpList::Linf(const ExpList &Sol)
        {
            ASSERTL2(Sol.GetPhysState() == true,
                     "local physical space is not true ");
            
            std::vector<StdRegions::StdExpansionVector>::iterator  sdef;
            StdRegions::StdExpansionVectorIter def;
            NekDouble err = 0.0;
            int       i,cnt = 0;
            ConstArray<OneD, NekDouble> soln = Sol.GetPhys();
            ConstArray<OneD, NekDouble> phys = m_phys;
            
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
            
            NekDouble err = 0.0,errl2;
            int    i,cnt = 0;
            ConstArray<OneD, NekDouble> soln = Sol.GetPhys();
            ConstArray<OneD, NekDouble> phys = m_phys;
            
            if(m_physState == false)
            {
                BwdTrans(*this);
            }
            
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

