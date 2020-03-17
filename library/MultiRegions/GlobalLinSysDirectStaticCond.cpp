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

#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysDirect
         *
         * Solves a linear system using single- or multi-level static
         * condensation.
         */

        /**
         * Registers the class with the Factory.
         */
        std::string GlobalLinSysDirectStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "DirectStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct static condensation.");

        std::string GlobalLinSysDirectStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "DirectMultiLevelStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct multi-level static condensation.");

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
        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(
                     const GlobalLinSysKey          &pKey,
                     const std::weak_ptr<ExpList> &pExpList,
                     const std::shared_ptr<AssemblyMap>
                     &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysDirect    (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eDirectStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eDirectMultiLevelStaticCond),
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
        GlobalLinSysDirectStaticCond::GlobalLinSysDirectStaticCond(
                     const GlobalLinSysKey                &pKey,
                     const std::weak_ptr<ExpList>         &pExpList,
                     const DNekScalBlkMatSharedPtr         pSchurCompl,
                     const DNekScalBlkMatSharedPtr         pBinvD,
                     const DNekScalBlkMatSharedPtr         pC,
                     const DNekScalBlkMatSharedPtr         pInvD,
                     const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysDirect    (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl = pSchurCompl;
            m_BinvD      = pBinvD;
            m_C          = pC;
            m_invD       = pInvD;
        }


        /**
         *
         */
        GlobalLinSysDirectStaticCond::~GlobalLinSysDirectStaticCond()
        {

        }


        MatrixStorage GlobalLinSysDirectStaticCond::DetermineMatrixStorage(
            const AssemblyMapSharedPtr &pLocToGloMap)
        {
            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int rows = nBndDofs - NumDirBCs;
            int bwidth = pLocToGloMap->GetBndSystemBandWidth();

            MatrixStorage matStorage;
            
            switch(m_linSysKey.GetMatrixType())
            {
                // case for all symmetric matices
            case StdRegions::eMass:
            case StdRegions::eLaplacian:
            case StdRegions::eHelmholtz:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    if( (2*(bwidth+1)) < rows)
                    {
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED; 
                    }
                    else
                    {
                        matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                    }
                }
                break;
            case StdRegions::eLinearAdvectionReaction:
            case StdRegions::eLinearAdvectionDiffusionReaction:
            default:
                {
                    // Current inversion techniques do not seem to
                    // allow banded matrices to be used as a linear
                    // system
                    matStorage = eFULL;
                }
                break;
            }

            return matStorage;
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysDirectStaticCond::v_AssembleSchurComplement(
            const AssemblyMapSharedPtr pLocToGloMap)
        {
            int i, j, n, cnt, gid1, gid2;
            NekDouble sign1, sign2, value;

            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();

            DNekScalBlkMatSharedPtr SchurCompl = m_schurCompl;
            DNekScalBlkMatSharedPtr BinvD      = m_BinvD;
            DNekScalBlkMatSharedPtr C          = m_C;
            DNekScalBlkMatSharedPtr invD       = m_invD;

            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;

            DNekMatSharedPtr Gmat;
            int bwidth = pLocToGloMap->GetBndSystemBandWidth();

            MatrixStorage matStorage = DetermineMatrixStorage(pLocToGloMap);

            switch(matStorage)
            {
                case ePOSITIVE_DEFINITE_SYMMETRIC_BANDED:
                {
                    try {
                        Gmat = MemoryManager<DNekMat>
                            ::AllocateSharedPtr(rows, cols, 0.0,
                                                matStorage,
                                                bwidth, bwidth);
                    }
                    catch (...) {
                        NEKERROR(ErrorUtil::efatal,
                                 "Insufficient memory for GlobalLinSys.");
                    }
                    break;
                }

                case ePOSITIVE_DEFINITE_SYMMETRIC:
                case eFULL:
                {
                    Gmat = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(rows, cols, 0.0, matStorage);
                    break;
                }

                default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Unknown matrix storage type of type not set up");
                }
            }
            
            // fill global matrix
            DNekScalMatSharedPtr loc_mat;
            int loc_lda;
            for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = SchurCompl->GetBlock(n,n);
                loc_lda = loc_mat->GetRows();

                // Set up  Matrix;
                for(i = 0; i < loc_lda; ++i)
                {
                    gid1  = pLocToGloMap->GetLocalToGlobalBndMap (cnt + i)
                                                                    - NumDirBCs;
                    sign1 = pLocToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = pLocToGloMap->GetLocalToGlobalBndMap(cnt+j)
                                                                 - NumDirBCs;
                            sign2 = pLocToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if(gid2 >= 0)
                            {
                                // As the global matrix should be symmetric,
                                // only add the value for the upper triangular
                                // part in order to avoid entries to be entered
                                // twice
                                if((matStorage == eFULL)||(gid2 >= gid1))
                                {
                                    value = Gmat->GetValue(gid1,gid2)
                                                + sign1*sign2*(*loc_mat)(i,j);
                                    Gmat->SetValue(gid1,gid2,value);
                                }
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            if(rows)
            {
                PointerWrapper w = eWrapper;
                m_linSys = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmat,w);
            }
        }

        GlobalLinSysStaticCondSharedPtr GlobalLinSysDirectStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &l2gMap)
        {
            GlobalLinSysDirectStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysDirectStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
