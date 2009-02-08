///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/LocalToGlobalC0ContMap.h>

namespace Nektar
{
    namespace MultiRegions
    {

	GlobalLinSys::GlobalLinSys(const GlobalLinSysKey &mkey, 
                                   const DNekLinSysSharedPtr linsys):
            m_linSysKey(mkey),
            m_linSys(linsys)
	{
	}

	GlobalLinSys::GlobalLinSys(const GlobalLinSysKey &mkey, 
                                   const DNekLinSysSharedPtr linsys,
                                   const DNekScalBlkMatSharedPtr BinvD,
                                   const DNekScalBlkMatSharedPtr C,
                                   const DNekScalBlkMatSharedPtr invD):
            m_linSysKey(mkey),
            
            m_linSys(linsys)
	{
            m_blkMatrices = Array<OneD,DNekScalBlkMatSharedPtr>(3);

            m_blkMatrices[0] = BinvD;
            m_blkMatrices[1] = C;
            m_blkMatrices[2] = invD;
	}


        void GlobalLinSys::Solve(const Array<OneD,const NekDouble> &in, 
                                 Array<OneD,NekDouble>  &out)
        {
            DNekVec Vin(in.num_elements(),in);
            DNekVec Vout(out.num_elements(),out,eWrapper);
            m_linSys->Solve(Vin,Vout);
        }

        void GlobalLinSys::Solve(const Array<OneD, const NekDouble> &in, 
                                 Array<OneD,NekDouble>  &out,
                                 LocalToGlobalC0ContMap  &locToGloMap)
        {
            switch(m_linSysKey.GetGlobalSysSolnType())
            {
            case eDirectFullMatrix:
                {
                    Solve(in,out); 
                }
                break;
            case eDirectStaticCond:
                {

#ifdef NEKTAR_USING_DIRECT_BLAS_CALLS
                    int nDirDofs  = locToGloMap.GetNumGlobalDirBndCoeffs();
                    int nbndry    = locToGloMap.GetNumGlobalBndCoeffs() - nDirDofs;
                    int nlocbndry = locToGloMap.GetNumLocalBndCoeffs();
                    int nint      = in.num_elements() -nbndry; 
                    

                    int wspsize = (in.get() == out.get())?(nlocbndry+nbndry+nint):(nlocbndry+nint);
                    Array<OneD,NekDouble> wsp(wspsize);
                    Vmath::Vcopy(nbndry+nint,in.get(),1,wsp.get()+nlocbndry,1);

                    DNekVec Vbnd(nbndry,out,eWrapper);

                    int i;
                    int cnt1;
                    int cnt2;
                    int nbndry_el;
                    int nint_el;
                    int nblocks = m_blkMatrices[0]->GetNumberOfBlockRows();

                    if(nbndry)
                    {
                        if(nint)
                        {
                            DNekScalBlkMat &BinvD = *m_blkMatrices[0];

                            // construct boundary forcing                                 
                            for(i = cnt1 = cnt2 = 0; i < nblocks; i++)
                            {
                                nbndry_el = BinvD.GetNumberOfRowsInBlockRow(i);
                                nint_el   = BinvD.GetNumberOfColumnsInBlockColumn(i);

                                if(nint_el)
                                {
                                    DNekScalMat& BinvD_el = *(BinvD.GetBlock(i,i));
                                    Blas::Dgemv('N',nbndry_el,nint_el,BinvD_el.Scale(),BinvD_el.GetRawPtr(),
                                                nbndry_el, in.get()+nbndry+cnt2, 1.0, 0.0, wsp.get()+cnt1, 1.0);
                                }

                                cnt1 += nbndry_el;
                                cnt2 += nint_el;
                            }      
 
                            locToGloMap.AssembleBnd(wsp,Vbnd.GetPtr(),nDirDofs);
                            Vmath::Vsub(nbndry,wsp.get()+nlocbndry,1,out.get(),1,out.get(),1); 
                        }
                        
                        // solve boundary system 
                        m_linSys->Solve(Vbnd,Vbnd);
                    }

                    // solve interior system 
                    if(nint)
                    {
                        DNekScalBlkMat &C     = *m_blkMatrices[1];
                        DNekScalBlkMat &invD  = *m_blkMatrices[2];

                        Vmath::Zero(nlocbndry,wsp.get(),1);
                        locToGloMap.GlobalToLocalBnd(Vbnd.GetPtr(),wsp,nDirDofs);
                             
                        for(i = cnt1 = cnt2 = 0; i < nblocks; i++)
                        {       
                            nint_el   = C.GetNumberOfRowsInBlockRow(i);
                            nbndry_el = C.GetNumberOfColumnsInBlockColumn(i);                     

                            if(nbndry_el)
                            {                                
                                DNekScalMat& C_el = *(C.GetBlock(i,i));
                                
                                Blas::Dgemv('N',nint_el,nbndry_el,(-1.0*C_el.Scale()),C_el.GetRawPtr(),
                                            nint_el, wsp.get()+cnt2, 1.0, 1.0, wsp.get()+nlocbndry+nbndry+cnt1, 1.0);
                            }      
                    
                            if(nint_el)
                            {              
                                DNekScalMat& invD_el = *(invD.GetBlock(i,i));
                                Blas::Dgemv('N',nint_el,nint_el,invD_el.Scale(),invD_el.GetRawPtr(),
                                            nint_el, wsp.get()+nlocbndry+nbndry+cnt1, 1.0, 0.0, out.get()+nbndry+cnt1, 1.0);
                            }
                            
                            cnt1 += nint_el;
                            cnt2 += nbndry_el;
                        }                
                    }                
#else
                    int nDirDofs = locToGloMap.GetNumGlobalDirBndCoeffs();
                    int nbndry  = locToGloMap.GetNumGlobalBndCoeffs() - nDirDofs;
                    int nlocbndry = locToGloMap.GetNumLocalBndCoeffs();
                    int nint    = in.num_elements() -nbndry; 

                    Array<OneD,NekDouble>  offset;  
                    DNekVec Fbnd(nbndry,in);
                    DNekVec Vloc(nlocbndry);
                    DNekVec Vbnd(nbndry,out,eWrapper);

                    DNekVec Fint(nint,in + nbndry);
                    DNekVec Vint(nint,offset = out + nbndry,eWrapper);

                    if(nbndry)
                    {
                        if(nint)
                        {
                            DNekScalBlkMat &BinvD = *m_blkMatrices[0];

                            // construct boundary forcing 
                            Vloc = BinvD*Fint;
                            locToGloMap.AssembleBnd(Vloc,Vbnd,nDirDofs);
                            Fbnd = Fbnd - Vbnd;
                        }
                        
                        // solve boundary system 
                        m_linSys->Solve(Fbnd,Vbnd);
                    }

                    // solve interior system 
                    if(nint)
                    {
                        DNekScalBlkMat &invD  = *m_blkMatrices[2];

                        if(nbndry)
                         {
                             DNekScalBlkMat &C     = *m_blkMatrices[1];

                             Vmath::Zero(Vloc.GetDimension(),Vloc.GetRawPtr(),1);
                             locToGloMap.GlobalToLocalBnd(Vbnd,Vloc,nDirDofs);
                             Fint = Fint - C*Vloc;
                         }

                        Vint = invD*Fint;
                    }
#endif
                }
                break;
            default:
                ASSERTL0(false,"Matrix solution type not defined");
                break;
            }
        }
	
    } //end of namespace
} //end of namespace

