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

                //Array<OneD, NekDouble> offsetarray;
                //SolveLinearSystem(nGlobDofs, tmp + nDirDofs,
				//  offsetarray = pOutput + nDirDofs, pLocToGloMap, nDirDofs);
                unsigned int nCoeffs = pLocToGloMap->GetNumLocalCoeffs();
                Array<OneD, NekDouble> vLocalIn(nLocNondir, 0.0);
                Array<OneD, NekDouble> vLocalOut(nLocNondir, 0.0);
                GlobalToLocalNonDir(tmp, vLocalIn, pLocToGloMap);
                Xxt::Solve(vLocalOut, m_crsData, vLocalIn);
                LocalNonDirToGlobal(vLocalOut, pOutput, pLocToGloMap);
            }
            else
            {
                SolveLinearSystem(pLocToGloMap->GetNumGlobalCoeffs(), pInput,pOutput, pLocToGloMap);
            }
        }


        void GlobalLinSysXxtFull::CreateMap(const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            const Array<OneD, const int> &vMap = pLocToGloMap->GetLocalToGlobalMap();
            const Array<OneD, NekDouble> &vMapSign = pLocToGloMap->GetLocalToGlobalSign();
            unsigned int nGloBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int i,j;
            unsigned int nEntries = 0;

            Array<OneD, int> vCounts(pLocToGloMap->GetNumGlobalCoeffs());
            for (i = j = 0; i < pLocToGloMap->GetNumLocalCoeffs(); ++i)
            {
                if (vMap[i] >= nGloBnd)
                {
                    nEntries++;
                }
                vCounts[vMap[i]]++;
            }
            m_locNonDirToGloNonDir = Array<OneD, int>(nEntries);
            m_locNonDirToGloNonDirSign = Array<OneD, NekDouble>(nEntries);
            m_locNonDirToGloNonDirSignMultiplicity = Array<OneD, NekDouble>(nEntries);

            for (i = j = 0; i < pLocToGloMap->GetNumLocalCoeffs(); ++i)
            {
                if (vMap[i] >= nGloBnd)
                {
                    m_locNonDirToGloNonDir[j] = vMap[i];
                    m_locNonDirToGloNonDirSign[j] = vMapSign[i];
                    m_locNonDirToGloNonDirSignMultiplicity[j] = 1.0/vCounts[vMap[i]]*vMapSign[i];
                    j++;
                }
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
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount = 0;
            unsigned int rCount = 0;
            unsigned int nRows = 0;
            unsigned int i = 0, j = 0, k = 0, n = 0, cnt = 0, iCnt = 0, jCnt = 0;
            unsigned int nEntries = 0;
            unsigned int numDirBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int gid1, gid2;
            Array<OneD, unsigned int> vSizes(vExp->GetNumElmts());

            for (n = 0; n < vExp->GetNumElmts(); ++n)
            {
                // todo: no need to actually get the matrix here - just number of coeffs
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                nRows = loc_mat->GetRows();
                for (i = j = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i) - numDirBnd;
                    if (gid1 >= 0)
                    {
                        j++;
                    }
                }
                vSizes[n] = j;
                nEntries += vSizes[n]*vSizes[n];
                nLocNondir += vSizes[n];
                cnt += nRows;
            }

            m_Ai = Array<OneD, unsigned int>(nEntries);
            m_Aj = Array<OneD, unsigned int>(nEntries);
            m_Ar = Array<OneD, double>(nEntries);
            Array<OneD, unsigned long> vId(nLocNondir);

            iCount = 0;
            for(n = cnt = 0; n < vExp->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                nRows = loc_mat->GetRows();

                iCnt = 0;
                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i)-numDirBnd;
                    if(gid1 >= 0)
                    {
                        jCnt = 0;
                        for(j = 0; j < nRows; ++j)
                        {
                            gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                    - numDirBnd;
                            if(gid2 >= 0)
                            {
                                k = rCount + iCnt*vSizes[n] + jCnt;
                                m_Ai[k] = iCount + iCnt;
                                m_Aj[k] = iCount + jCnt;
                                m_Ar[k] = (*loc_mat)(i,j);
                                jCnt++;
                            }

                        }
                        vId[iCount + iCnt] = gid1 + numDirBnd;
                        iCnt++;
                    }
                }
                cnt   += nRows;
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }

            LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm();
            m_crsData = Xxt::Init(nLocNondir, vId, m_Ai, m_Aj, m_Ar, vComm);
            nektar_crs_stats(m_crsData);
/*
            ExpListSharedPtr vExp = m_expList.lock();
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount = 0;
            unsigned int rCount = 0;
            unsigned int nRows = 0;
            unsigned int i = 0, j = 0, k = 0, n = 0, cnt = 0, iCnt = 0, jCnt = 0;
            unsigned int nEntries = 0;
            unsigned int nLocal = pLocToGloMap->GetNumLocalCoeffs();
            unsigned int nDirBnd = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int gid1, gid2;
            Array<OneD, unsigned int> vSizes(vExp->GetNumElmts());

            for (n = 0; n < vExp->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                vSizes[n] = loc_mat->GetRows();
                cout << "vSize " << n << " = " << vSizes[n] << endl;
                nEntries += vSizes[n]*vSizes[n];
            }

            m_Ai = Array<OneD, unsigned int>(nEntries);
            m_Aj = Array<OneD, unsigned int>(nEntries);
            m_Ar = Array<OneD, double>(nEntries);
            Array<OneD, unsigned long> vId(nLocal);

            iCount = 0;
            for(n = cnt = 0; n < vExp->GetNumElmts(); ++n)
            {
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                nRows = loc_mat->GetRows();

                iCnt = 0;
                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(cnt + i)-nDirBnd;
                    jCnt = 0;
                    for(j = 0; j < nRows; ++j)
                    {
                        gid2 = pLocToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                - nDirBnd;
                        k = rCount + iCnt*vSizes[n] + jCnt;
                        m_Ai[k] = iCount + iCnt;
                        m_Aj[k] = iCount + jCnt;
                        if(gid1 >= 0 && gid2 >= 0)
                        {
                            m_Ar[k] = (*loc_mat)(i,j);
                        }
                        else
                        {
                            m_Ar[k] = 0.0;
                        }
                        jCnt++;

                    }
                    vId[iCount + iCnt] = gid1 + numDirBnd;
                    iCnt++;
                }
                cnt   += nRows;
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }
*/
        }


        void GlobalLinSysXxtFull::GlobalToLocalNonDir(const Array<OneD, const NekDouble> &global,
                                       Array<OneD, NekDouble> &local,
                               const boost::shared_ptr<AssemblyMap>
                                                      &pLocToGloMap)
        {
            int n = m_locNonDirToGloNonDir.num_elements();
            Vmath::Gathr(n, m_locNonDirToGloNonDirSignMultiplicity.get(), global.get(), m_locNonDirToGloNonDir.get(), local.get());
        }

        void GlobalLinSysXxtFull::LocalNonDirToGlobal(const Array<OneD, const NekDouble> &local,
                                       Array<OneD, NekDouble> &global,
                               const boost::shared_ptr<AssemblyMap>
                                                      &pLocToGloMap)
        {
            int n = m_locNonDirToGloNonDir.num_elements();
            Vmath::Scatr(n, m_locNonDirToGloNonDirSign.get(), local.get(), m_locNonDirToGloNonDir.get(), global.get());
        }

    }
}
