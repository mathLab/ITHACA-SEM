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
            int nGlobDofs = pLocToGloMap->GetNumGlobalCoeffs();

            Array<OneD, NekDouble> tmp(nGlobDofs);
            Array<OneD, NekDouble> tmp2(nGlobDofs);
            Array<OneD, NekDouble> tmp3 = pOutput + nDirDofs;

            if(nDirDofs)
            {
                // calculate the dirichlet forcing
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
                    //int nLocDofs = pLocToGloMap->GetNumLocalCoeffs();

                    m_expList.lock()->GeneralMatrixOp(
                            m_linSysKey,
                            pOutput, tmp, eGlobal);

                    Vmath::Vsub( nGlobDofs, pInput.get(),1,
                                            tmp.get(),   1,
                                            tmp.get(),   1);
                }

                SolveLinearSystem(pLocToGloMap->GetNumLocalCoeffs(),
                                    tmp, tmp2, pLocToGloMap);

                // Enforce the Dirichlet boundary conditions on the solution
                // array as XXT discards them.
                Vmath::Vcopy(nDirDofs, pOutput, 1,
                                       tmp2,    1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs, pInput, 1, tmp, 1);
                SolveLinearSystem(pLocToGloMap->GetNumLocalCoeffs(),
                                    tmp,tmp2, pLocToGloMap);
            }

            // Perturb the output array (previous solution) by the result of
            // this solve to get full solution.
            Vmath::Vadd(nGlobDofs - nDirDofs,
                        tmp2 + nDirDofs, 1, tmp3, 1, tmp3, 1);

        }


        /**
         * Create the inverse multiplicity map.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtFull::CreateMap(
                    const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            const Array<OneD, const int> &vMap
                                    = pLocToGloMap->GetLocalToGlobalMap();
            unsigned int nGloBnd    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int nGlo       = pLocToGloMap->GetNumGlobalCoeffs();
            unsigned int nEntries   = pLocToGloMap->GetNumLocalCoeffs();
            unsigned int i,j;

            // Count the multiplicity of each global DOF on this process
            Array<OneD, NekDouble> vCounts(nGlo, 0.0);
            for (i = 0; i < nEntries; ++i)
            {
                vCounts[vMap[i]] += 1.0;
            }

            // Get universal multiplicity by globally assembling counts
            pLocToGloMap->UniversalAssemble(vCounts);

            // Construct a map of 1/multiplicity for use in XXT solve
            m_locToGloSignMult = Array<OneD, NekDouble>(nEntries);
            for (i = 0; i < nEntries; ++i)
            {
                m_locToGloSignMult[i] = 1.0/vCounts[vMap[i]];
            }

            m_map = pLocToGloMap->GetLocalToGlobalMap();
        }


        /**
         * Construct the local matrix row index, column index and value index
         * arrays and initialize the XXT data structure with this information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtFull::AssembleMatrixArrays(
                        const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            ExpListSharedPtr vExp = m_expList.lock();
            unsigned int nElmt = vExp->GetNumElmts();
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount     = 0;
            unsigned int rCount     = 0;
            unsigned int nRows      = 0;
            unsigned int nEntries   = 0;
            unsigned int numDirBnd  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int nLocal     = pLocToGloMap->GetNumLocalCoeffs();
            const Array<OneD, NekDouble> &vMapSign
                                    = pLocToGloMap->GetLocalToGlobalSign();
            bool doSign = pLocToGloMap->GetSignChange();
            unsigned int i = 0, j = 0, k = 0, n = 0;
            int gid1;
            Array<OneD, unsigned int> vSizes(nElmt);

            // First construct a map of the number of local DOFs in each block
            // and the number of matrix entries for each block
            for (n = 0; n < nElmt; ++n)
            {
                i = vExp->GetOffset_Elmt_Id(n);
                vSizes[n] = vExp->GetExp(i)->GetNcoeffs();
                nEntries += vSizes[n]*vSizes[n];
            }

            // Set up i-index, j-index and value arrays
            m_Ai = Array<OneD, unsigned int>(nEntries);
            m_Aj = Array<OneD, unsigned int>(nEntries);
            m_Ar = Array<OneD, double>(nEntries, 0.0);

            // Set up the universal ID array for XXT
            Array<OneD, unsigned long> vId(nLocal);

            // Loop over each elemental block, extract matrix indices and value
            // and set the universal ID array
            for(n = iCount = 0; n < nElmt; ++n)
            {
                loc_mat = GetBlock(vExp->GetOffset_Elmt_Id(n));
                nRows = loc_mat->GetRows();

                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalMap(iCount + i);
                    for(j = 0; j < nRows; ++j)
                    {
                            k = rCount + i*vSizes[n] + j;
                            m_Ai[k] = iCount + i;
                            m_Aj[k] = iCount + j;
                            m_Ar[k] = (*loc_mat)(i,j);
                            if (doSign)
                            {
                                m_Ar[k] *= vMapSign[iCount+i]*vMapSign[iCount+j];
                            }
                    }

                    // Dirichlet DOFs are not included in the solve, so we set
                    // these to the special XXT id=0.
                    if (gid1 < numDirBnd)
                    {
                        vId[iCount + i] = 0;
                    }
                    else
                    {
                        vId[iCount + i]
                            = pLocToGloMap->GetGlobalToUniversalMap(gid1);
                    }
                }
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }

            // Set up XXT and output some stats
            LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm();
            m_crsData = Xxt::Init(nLocal, vId, m_Ai, m_Aj, m_Ar, vComm);
            Xxt::nektar_crs_stats(m_crsData);
        }
    }
}
