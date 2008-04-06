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


        void GlobalLinSys::Solve(const Array<OneD, const NekDouble> &in, 
                                 Array<OneD,NekDouble>  &out,
                                 LocalToGlobalBndryMap  &locToGloMap)
        {
            switch(m_linSysKey.GetGlobalSysSolnType())
            {
            case eDirectFullMatrix:
                {
                    DNekVec Vin(in.num_elements(),in);
                    DNekVec Vout(out.num_elements(),out,eWrapper);
                    m_linSys->Solve(Vin,Vout);
                }
                break;
            case eDirectStaticCond:
                {
                    int i;

                    int nDirDofs = locToGloMap.GetNumDirichletDofs();
                    int nbndry  = locToGloMap.GetTotGloBndDofs() - nDirDofs;
                    int nint    = in.num_elements() -nbndry; 

                    Array<OneD,NekDouble>  offset;  
                    DNekVec Fbnd(nbndry,in);
                    DNekVec Vloc;
                    DNekVec Vbnd(nbndry,out,eWrapper);

                    DNekVec Fint(nint,in + nbndry);
                    DNekVec Vint(nint,offset = out + nbndry,eWrapper);

                    if(nbndry)
                    {
                        DNekScalBlkMat &BinvD = *m_blkMatrices[0];
                    
                        // construct boundary forcing 
                        Vloc = BinvD*Fint;
                        locToGloMap.AssembleBnd(Vloc,Vbnd,nDirDofs);
                        Fbnd = Fbnd - Vbnd;
                        
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
                             locToGloMap.ContToLocalBnd(Vbnd,Vloc,nDirDofs);
                             Fint = Fint - C*Vloc;
                         }
                        Vint = invD*Fint;
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

