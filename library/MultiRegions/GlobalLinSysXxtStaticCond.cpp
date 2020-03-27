///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeStaticCond.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: GlobalLinSysIterativeStaticCond definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Xxt.hpp>
#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/GlobalLinSysXxtStaticCond.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysIterativeStaticCond
         *
         * Solves a linear system iteratively using single- or multi-level
         * static condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        string GlobalLinSysXxtStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "XxtStaticCond",
                    GlobalLinSysXxtStaticCond::create,
                    "Iterative static condensation.");

        string GlobalLinSysXxtStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "XxtMultiLevelStaticCond",
                    GlobalLinSysXxtStaticCond::create,
                    "Xxt multi-level static condensation.");

        /**
         * For a matrix system of the form @f[
         * \left[ \begin{array}{cc}
         * \boldsymbol{A} & \boldsymbol{B}\\
         * \boldsymbol{C} & \boldsymbol{D}
         * \end{array} \right]
         * \left[ \begin{array}{c} \boldsymbol{x_1}\\ \boldsymbol{x_2}
         * \end{array}\right]
         * = \left[ \begin{array}{c} \boldsymbol{y_1}\\ \boldsymbol{y_2}
         * \end{array}\right],
         * @f]
         * where @f$\boldsymbol{D}@f$ and
         * @f$(\boldsymbol{A-BD^{-1}C})@f$ are invertible, store and assemble
         * a static condensation system, according to a given local to global
         * mapping. #m_linSys is constructed by AssembleSchurComplement().
         * @param   mKey        Associated matrix key.
         * @param   pLocMatSys  LocalMatrixSystem
         * @param   locToGloMap Local to global mapping.
         */
        GlobalLinSysXxtStaticCond::GlobalLinSysXxtStaticCond(
                     const GlobalLinSysKey &pKey,
                     const std::weak_ptr<ExpList> &pExpList,
                     const std::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysXxt       (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eXxtStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eXxtMultiLevelStaticCond),
                     "This constructor is only valid when using static "
                     "condensation");
            ASSERTL1(pKey.GetGlobalSysSolnType()
                        == pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested "
                     "solution type");
        }


        /**
         *
         */
        GlobalLinSysXxtStaticCond::GlobalLinSysXxtStaticCond(
                     const GlobalLinSysKey &pKey,
                     const std::weak_ptr<ExpList> &pExpList,
                     const DNekScalBlkMatSharedPtr pSchurCompl,
                     const DNekScalBlkMatSharedPtr pBinvD,
                     const DNekScalBlkMatSharedPtr pC,
                     const DNekScalBlkMatSharedPtr pInvD,
                     const std::shared_ptr<AssemblyMap>
                                                            &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysXxt       (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl  = pSchurCompl;
            m_BinvD       = pBinvD;
            m_C           = pC;
            m_invD        = pInvD;
            m_locToGloMap = pLocToGloMap;
        }


        /**
         *
         */
        GlobalLinSysXxtStaticCond::~GlobalLinSysXxtStaticCond()
        {

        }

        /**
         * Create the inverse multiplicity map.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtStaticCond::CreateMap(
                    const std::shared_ptr<AssemblyMap> &pLocToGloMap)
        {
            const Array<OneD, const int> &vMap
                                    = pLocToGloMap->GetLocalToGlobalBndMap();
            unsigned int nGlo       = pLocToGloMap->GetNumGlobalBndCoeffs();
            unsigned int nEntries   = pLocToGloMap->GetNumLocalBndCoeffs();
            unsigned int i;

            // Count the multiplicity of each global DOF on this process
            Array<OneD, NekDouble> vCounts(nGlo, 0.0);
            for (i = 0; i < nEntries; ++i)
            {
                vCounts[vMap[i]] += 1.0;
            }

            // Get universal multiplicity by globally assembling counts
            pLocToGloMap->UniversalAssembleBnd(vCounts);

            // Construct a map of 1/multiplicity for use in XXT solve
            m_locToGloSignMult = Array<OneD, NekDouble>(nEntries);
            for (i = 0; i < nEntries; ++i)
            {
                m_locToGloSignMult[i] = 1.0/vCounts[vMap[i]];
            }

            m_map = pLocToGloMap->GetLocalToGlobalBndMap();
        }

        /**
         * Construct the local matrix row index, column index and value index
         * arrays and initialize the XXT data structure with this information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysXxtStaticCond::v_AssembleSchurComplement(
            std::shared_ptr<AssemblyMap> pLocToGloMap)
        {
            CreateMap(pLocToGloMap);

            ExpListSharedPtr vExp = m_expList.lock();
            unsigned int nElmt = m_schurCompl->GetNumberOfBlockRows();
            DNekScalMatSharedPtr loc_mat;
            unsigned int iCount     = 0;
            unsigned int rCount     = 0;
            unsigned int nRows      = 0;
            unsigned int nEntries   = 0;
            unsigned int numDirBnd  = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int nLocal     = pLocToGloMap->GetNumLocalBndCoeffs();
            const Array<OneD, NekDouble> &vMapSign
                                    = pLocToGloMap->GetLocalToGlobalBndSign();
            bool doSign = pLocToGloMap->GetSignChange();
            unsigned int i = 0, j = 0, k = 0, n = 0;
            int gid1;
            Array<OneD, unsigned int> vSizes(nElmt);

            // First construct a map of the number of local DOFs in each block
            // and the number of matrix entries for each block
            for (n = 0; n < nElmt; ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
                vSizes[n] = loc_mat->GetRows();
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
                loc_mat = m_schurCompl->GetBlock(n,n);
                nRows = loc_mat->GetRows();

                for(i = 0; i < nRows; ++i)
                {
                    gid1 = pLocToGloMap->GetLocalToGlobalBndMap(iCount + i);
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
                            = pLocToGloMap->GetGlobalToUniversalBndMap()[gid1];
                    }
                }
                iCount += vSizes[n];
                rCount += vSizes[n]*vSizes[n];
            }

            // Set up XXT and output some stats
            LibUtilities::CommSharedPtr vComm = pLocToGloMap->GetComm()->GetRowComm();
            m_crsData = Xxt::Init(nLocal, vId, m_Ai, m_Aj, m_Ar, vComm);
            if (m_verbose)
            {
                Xxt::nektar_crs_stats(m_crsData);
            }
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysXxtStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &l2gMap)
        {
            GlobalLinSysXxtStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysXxtStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
