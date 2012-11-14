///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysXxtFull.cpp
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
// Description: GlobalLinSysXxtFull definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Xxt.hpp>
#include <MultiRegions/GlobalLinSysXxtFull.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysXxtFull
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysXxtFull::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "XxtFull",
                    GlobalLinSysXxtFull::create,
                    "Xxt Full Matrix.");


        /// Constructor for full direct matrix solve.
        GlobalLinSysXxtFull::GlobalLinSysXxtFull(
                    const GlobalLinSysKey &pLinSysKey,
                    const boost::weak_ptr<ExpList> &pExp,
                    const boost::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
                : GlobalLinSysXxt(pLinSysKey, pExp, pLocToGloMap)
        {

            ASSERTL1(m_linSysKey.GetGlobalSysSolnType()==eXxtFullMatrix,
                     "This routine should only be used when using a Full XXT"
                     " matrix solve");

            CreateMap(pLocToGloMap);
            AssembleMatrixArrays(pLocToGloMap);
        }


        GlobalLinSysXxtFull::~GlobalLinSysXxtFull()
        {

        }


        /**
         * Solve the linear system using a full global matrix system.
         */
        void GlobalLinSysXxtFull::v_Solve(
                    const Array<OneD, const NekDouble>  &pInput,
                          Array<OneD,       NekDouble>  &pOutput,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &pDirForcing)
        {
            bool dirForcCalculated = (bool) pDirForcing.num_elements();
            int nDirDofs  = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            if(nDirDofs)
            {
                // calculate the dirichlet forcing
                int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();
                Array<OneD, NekDouble> tmp(nGlobDofs);
                Array<OneD, NekDouble> tmp2(nGlobDofs);
                if(dirForcCalculated)
                {
                    Vmath::Vsub(nGlobDofs, pInput.get(), 1,
                                pDirForcing.get(), 1,
                                tmp.get(), 1);
                }
                else
                {
                    // Calculate the dirichlet forcing and substract it
                    // from the rhs
                    int nLocDofs = pLocToGloMap->GetNumLocalCoeffs();

                    m_expList.lock()->GeneralMatrixOp(
                            m_linSysKey,
                            pOutput, tmp, true);

                    Vmath::Vsub( nGlobDofs, pInput.get(),1,
                                            tmp.get(),   1,
                                            tmp.get(),   1);
                }

                SolveLinearSystem(nGlobDofs, tmp, tmp2, pLocToGloMap);

                // Put back the Dirichlet boundary conditions
                Vmath::Vadd(nGlobDofs, pOutput.get(), 1, tmp2.get(), 1, pOutput.get(), 1);
            }
            else
            {
                SolveLinearSystem(pLocToGloMap->GetNumGlobalCoeffs(), pInput,pOutput, pLocToGloMap);
            }
        }


        void GlobalLinSysXxtFull::CreateMap(const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            const Array<OneD, const int> &vMap = pLocToGloMap->GetLocalToGlobalMap();
            unsigned int nGloBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int i,j;
            unsigned int nEntries = pLocToGloMap->GetNumLocalCoeffs();

            Array<OneD, int> vCounts(pLocToGloMap->GetNumGlobalCoeffs(), 0);
            for (i = 0; i < nEntries; ++i)
            {
                vCounts[vMap[i]]++;
            }

            m_locToGloSignMult = Array<OneD, NekDouble>(nEntries);
            for (i = 0; i < nEntries; ++i)
            {
                m_locToGloSignMult[i] = 1.0/vCounts[vMap[i]];
            }
        }


        /**
         * Assemble a full matrix from the block matrix stored in
         * #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtFull::AssembleMatrixArrays(
                        const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            ExpListSharedPtr vExp = m_expList.lock();
            unsigned int nElmt = vExp->GetNumElmts();
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount = 0;
            unsigned int rCount = 0;
            unsigned int nRows = 0;
            unsigned int i = 0, j = 0, k = 0, n = 0, cnt = 0, iCnt = 0, jCnt = 0;
            unsigned int nEntries = 0;
            unsigned int numDirBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int nLocal = pLocToGloMap->GetNumLocalCoeffs();
            const Array<OneD, NekDouble> &vMapSign = pLocToGloMap->GetLocalToGlobalSign();
            bool doSign = pLocToGloMap->GetSignChange();
            int gid1, gid2;
            Array<OneD, unsigned int> vSizes(nElmt);

            for (n = 0; n < nElmt; ++n)
            {
                nRows = vExp->GetExp(vExp->GetOffset_Elmt_Id(n))->GetNcoeffs();
                vSizes[n] = nRows;
                nEntries += vSizes[n]*vSizes[n];
            }

            m_Ai = Array<OneD, unsigned int>(nEntries);
            m_Aj = Array<OneD, unsigned int>(nEntries);
            m_Ar = Array<OneD, double>(nEntries, 0.0);
            Array<OneD, unsigned long> vId(nLocal);

            iCount = 0;
            for(n = cnt = 0; n < nElmt; ++n)
            {
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                nRows = loc_mat->GetRows();

                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i);
                    for(j = 0; j < nRows; ++j)
                    {
                            k = rCount + i*vSizes[n] + j;
                            m_Ai[k] = iCount + i;
                            m_Aj[k] = iCount + j;
                            m_Ar[k] = (*loc_mat)(i,j);
                            //m_Ar[k] /= m_locToGloSignMult[cnt+i]*m_locToGloSignMult[cnt+j];
                            if (doSign)
                            {
                                m_Ar[k] *= vMapSign[cnt+i]*vMapSign[cnt+j];
                            }
                    }
                    if (gid1 < numDirBnd)
                        vId[iCount + i] = 0;
                    else
                    {
                        //vId[iCount + i] = pLocToGloMap->GetGlobalToUniversalMap(gid1);
                        vId[iCount + i] = gid1;
                    }
                }
                cnt   += nRows;
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }

            LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm();
            m_crsData = Xxt::Init(nLocal, vId, m_Ai, m_Aj, m_Ar, vComm);
            Xxt::nektar_crs_stats(m_crsData);
        }
    }
}
