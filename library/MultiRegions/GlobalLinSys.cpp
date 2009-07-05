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
                                   const DNekScalBlkMatSharedPtr SchurCompl,
                                   const DNekScalBlkMatSharedPtr BinvD,
                                   const DNekScalBlkMatSharedPtr C,
                                   const DNekScalBlkMatSharedPtr invD):
            m_linSysKey(mkey),            
            m_linSys(linsys),
            m_blkMatrices(4)
	{
            m_blkMatrices[0] = SchurCompl;
            m_blkMatrices[1] = BinvD;
            m_blkMatrices[2] = C;
            m_blkMatrices[3] = invD;
	}


        void GlobalLinSys::Solve(const Array<OneD,const NekDouble> &in, 
                                 Array<OneD,NekDouble>  &out)
        {
            DNekVec Vin(in.num_elements(),in);
            DNekVec Vout(out.num_elements(),out,eWrapper);
            m_linSys->Solve(Vin,Vout);
        }

        void GlobalLinSys::Solve(const Array<OneD, const NekDouble> &in, 
                                       Array<OneD,       NekDouble> &out,
                                       LocalToGlobalC0ContMap  &locToGloMap,
                                       ExpList* exp,
                                 const Array<OneD, const NekDouble> &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            
            switch(m_linSysKey.GetGlobalSysSolnType())
            {
            case eDirectFullMatrix:
                {
                    int nDirDofs  = locToGloMap.GetNumGlobalDirBndCoeffs();
                    
                    if(nDirDofs)
                    {
                        // calculate the dirichlet forcing
                        int nGlobDofs      = locToGloMap.GetNumGlobalCoeffs();
                        Array<OneD, NekDouble> tmp(nGlobDofs);
                        if(dirForcCalculated)
                        {
                            Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,tmp.get(),1);
                        }
                        else
                        {
                            exp->GeneralMatrixOp(*(m_linSysKey.GetGlobalMatrixKey()),out,tmp,true);
                            Vmath::Vsub(nGlobDofs,in.get(),1,tmp.get(),1,tmp.get(),1);
                        }
                        Array<OneD, NekDouble> offsetarray;
                        Solve(tmp+nDirDofs,offsetarray = out+nDirDofs); 
                    }
                    else
                    {
                        Solve(in,out); 
                    }
                }
                break;
            case eDirectStaticCond:
                {
                    int nGlobDofs        = locToGloMap.GetNumGlobalCoeffs();
                    int nGlobBndDofs     = locToGloMap.GetNumGlobalBndCoeffs();
                    int nDirBndDofs      = locToGloMap.GetNumGlobalDirBndCoeffs();
                    int nGlobHomBndDofs  = nGlobBndDofs - nDirBndDofs;
                    int nLocBndDofs      = locToGloMap.GetNumLocalBndCoeffs();
                    int nIntDofs         = locToGloMap.GetNumGlobalCoeffs() - nGlobBndDofs; 

                    NekVector<NekDouble> F_HomBnd;
                    NekVector<NekDouble> F_Int;
                    Array<OneD, NekDouble> tmp;
                    if(nDirBndDofs && dirForcCalculated)
                    {
                        tmp = Array<OneD, NekDouble>(nGlobDofs);
                        Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,tmp.get(),1);
                        F_HomBnd = NekVector<NekDouble>(nGlobHomBndDofs,tmp+nDirBndDofs,eWrapper);
                        F_Int    = NekVector<NekDouble>(nIntDofs,tmp+nGlobBndDofs,eWrapper);
                    }
                    else
                    {
                        F_HomBnd = NekVector<NekDouble>(nGlobHomBndDofs,in+nDirBndDofs,eCopy);
                        F_Int    = NekVector<NekDouble>(nIntDofs,in+nGlobBndDofs,eCopy);
                    }

                    NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
                    NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,out+nDirBndDofs,eWrapper);
                    NekVector<NekDouble> V_Int(nIntDofs,out+nGlobBndDofs,eWrapper);
                    NekVector<NekDouble> V_LocBnd(nLocBndDofs);

                    if(nGlobHomBndDofs)
                    {
                        if(nIntDofs || ((nDirBndDofs) && (!dirForcCalculated)) )
                        {
                            // construct boundary forcing 
                            DNekScalBlkMat &BinvD = *m_blkMatrices[1];
                            DNekScalBlkMat &SchurCompl = *m_blkMatrices[0];

                            if( nIntDofs  && ((nDirBndDofs) && (!dirForcCalculated)) )
                            {
                                //include dirichlet boundary forcing
                                locToGloMap.GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                                V_LocBnd = SchurCompl*V_LocBnd + BinvD*F_Int;
                            }
                            else if((nDirBndDofs) && (!dirForcCalculated))
                            {
                                //include dirichlet boundary forcing
                                locToGloMap.GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                                V_LocBnd = SchurCompl*V_LocBnd;
                            }
                            else
                            {
                                V_LocBnd = BinvD*F_Int;
                            }
                            locToGloMap.AssembleBnd(V_LocBnd,V_GlobHomBnd,nDirBndDofs);
                            F_HomBnd = F_HomBnd - V_GlobHomBnd;
                        }
                        
                        // solve boundary system 
                        m_linSys->Solve(F_HomBnd,V_GlobHomBnd);
                    }

                    // solve interior system 
                    if(nIntDofs)
                    {
                        DNekScalBlkMat &invD  = *m_blkMatrices[3];

                        if(nGlobHomBndDofs)
                         {
                             DNekScalBlkMat &C     = *m_blkMatrices[2];

                             if(dirForcCalculated && nDirBndDofs)
                             {
                                 locToGloMap.GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,nDirBndDofs);
                             }
                             else
                             {                            
                                 locToGloMap.GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                             }
                             F_Int = F_Int - C*V_LocBnd;
                         }

                        V_Int = invD*F_Int;
                    }
                }
                break;
            default:
                ASSERTL0(false,"Matrix solution type not defined");
                break;
            }
        }
	
    } //end of namespace
} //end of namespace

