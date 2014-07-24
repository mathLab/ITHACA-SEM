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
        string GlobalLinSysDirectStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "DirectStaticCond",
                    GlobalLinSysDirectStaticCond::create,
                    "Direct static condensation.");

        string GlobalLinSysDirectStaticCond::className2
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
                     const boost::weak_ptr<ExpList> &pExpList,
                     const boost::shared_ptr<AssemblyMap>
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
                     const boost::weak_ptr<ExpList>       &pExpList,
                     const DNekScalBlkMatSharedPtr         pSchurCompl,
                     const DNekScalBlkMatSharedPtr         pBinvD,
                     const DNekScalBlkMatSharedPtr         pC,
                     const DNekScalBlkMatSharedPtr         pInvD,
                     const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
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

#if 0
        /**
         *
         */
        void GlobalLinSysDirectStaticCond::v_Solve(
                    const Array<OneD, const NekDouble>  &in,
                          Array<OneD,       NekDouble>  &out,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble>  &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = pLocToGloMap->AtLastLevel();

            int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = pLocToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = pLocToGloMap->GetNumGlobalCoeffs()
                                                                - nGlobBndDofs;

            Array<OneD, NekDouble> F(nGlobDofs);
            Array<OneD, NekDouble> tmp;

            if(nDirBndDofs && dirForcCalculated)
            {
                Vmath::Vsub(nGlobDofs,in.get(),1,dirForcing.get(),1,F.get(),1);
            }
            else
            {
                Vmath::Vcopy(nGlobDofs,in.get(),1,F.get(),1);
            }

            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,tmp=F+nDirBndDofs,
                                          eWrapper);
            NekVector<NekDouble> F_Int(nIntDofs,tmp=F+nGlobBndDofs,eWrapper);

            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,
                                              tmp=out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,tmp=out+nGlobBndDofs,eWrapper);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,0.0);

            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            //zero GlobHomBnd so that we ensure we are solving for
            //full problem rather than perturbation from initial
            //condition in this case

            Vmath::Zero(nGlobHomBndDofs,tmp = out+nDirBndDofs,1);

            if(nGlobHomBndDofs)
            {
                if(nIntDofs || ((nDirBndDofs) && (!dirForcCalculated)
                                              && (atLastLevel)) )
                {
                    // construct boundary forcing
                    if( nIntDofs  && ((nDirBndDofs) && (!dirForcCalculated)
                                                    && (atLastLevel)) )
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);

                        V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;

                    }
                    else if((nDirBndDofs) && (!dirForcCalculated)
                                          && (atLastLevel))
                    {
                        //include dirichlet boundary forcing
                        DNekScalBlkMat &SchurCompl = *m_schurCompl;
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                        V_LocBnd = SchurCompl*V_LocBnd;
                    }
                    else
                    {
                        DNekScalBlkMat &BinvD      = *m_BinvD;
                        V_LocBnd = BinvD*F_Int;
                    }
                    pLocToGloMap->AssembleBnd(V_LocBnd,V_GlobHomBndTmp,
                                             nDirBndDofs);
                    F_HomBnd = F_HomBnd - V_GlobHomBndTmp;
                }

                // solve boundary system
                if(atLastLevel)
                {
                    m_linSys->Solve(F_HomBnd,V_GlobHomBnd);
                }
                else
                {
                    m_recursiveSchurCompl->Solve(F,
                                V_GlobBnd.GetPtr(),
                                pLocToGloMap->GetNextLevelLocalToGlobalMap());
                }
            }

            // solve interior system
            if(nIntDofs)
            {
                DNekScalBlkMat &invD  = *m_invD;

                if(nGlobHomBndDofs || nDirBndDofs)
                {
                    DNekScalBlkMat &C     = *m_C;

                    if(dirForcCalculated && nDirBndDofs)
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd,
                                                      nDirBndDofs);
                    }
                    else
                    {
                        pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    }
                    F_Int = F_Int - C*V_LocBnd;
                }

                V_Int = invD*F_Int;
            }
        }
#endif

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
            const boost::weak_ptr<ExpList>       &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const boost::shared_ptr<AssemblyMap> &l2gMap)
        {
            GlobalLinSysDirectStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysDirectStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
