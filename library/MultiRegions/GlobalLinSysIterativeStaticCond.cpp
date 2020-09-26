///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterativeStaticCond.cpp
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
// Description: Implementation to linear solver using single-
//              or multi-level static condensation
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

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
        string GlobalLinSysIterativeStaticCond::className
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative static condensation.");

        string GlobalLinSysIterativeStaticCond::className2
                = GetGlobalLinSysFactory().RegisterCreatorFunction(
                    "IterativeMultiLevelStaticCond",
                    GlobalLinSysIterativeStaticCond::create,
                    "Iterative multi-level static condensation.");


        std::string GlobalLinSysIterativeStaticCond::storagedef = 
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "LocalMatrixStorageStrategy",
                "Sparse");
        std::string GlobalLinSysIterativeStaticCond::storagelookupIds[3] = {
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Contiguous",
                MultiRegions::eContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Non-contiguous",
                MultiRegions::eNonContiguous),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LocalMatrixStorageStrategy",
                "Sparse",
                MultiRegions::eSparse),
        };

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
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
            const GlobalLinSysKey                &pKey,
            const std::weak_ptr<ExpList>         &pExpList,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysIterative (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            ASSERTL1((pKey.GetGlobalSysSolnType()==eIterativeStaticCond)||
                     (pKey.GetGlobalSysSolnType()==eIterativeMultiLevelStaticCond),
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
        GlobalLinSysIterativeStaticCond::GlobalLinSysIterativeStaticCond(
            const GlobalLinSysKey                &pKey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &pLocToGloMap,
            const PreconditionerSharedPtr         pPrecon)
            : GlobalLinSys          (pKey, pExpList, pLocToGloMap),
              GlobalLinSysIterative (pKey, pExpList, pLocToGloMap),
              GlobalLinSysStaticCond(pKey, pExpList, pLocToGloMap)
        {
            m_schurCompl  = pSchurCompl;
            m_BinvD       = pBinvD;
            m_C           = pC;
            m_invD        = pInvD;
            m_precon      = pPrecon;
        }


        void GlobalLinSysIterativeStaticCond::v_InitObject()
        {
            auto asmMap = m_locToGloMap.lock();

            m_precon = CreatePrecon(asmMap);

            // Allocate memory for top-level structure
            SetupTopLevel(asmMap);

            // Setup Block Matrix systems
            int n, n_exp = m_expList.lock()->GetNumElmts();

            // Build preconditioner
            m_precon->BuildPreconditioner();

            // Do transform of Schur complement matrix
            int cnt = 0; 
            for (n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() !=
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr mat = m_schurCompl->GetBlock(n, n);
                    DNekScalMatSharedPtr t = m_precon->TransformedSchurCompl(
                                                n, cnt, mat);
                    m_schurCompl->SetBlock(n, n, t);
                    cnt += mat->GetRows();
                }
            }

            // Construct this level
            Initialise(asmMap);
        }
        
        /**
         *
         */
        GlobalLinSysIterativeStaticCond::~GlobalLinSysIterativeStaticCond()
        {
            
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCond::
            v_GetStaticCondBlock(unsigned int n)
        {
            DNekScalBlkMatSharedPtr schurComplBlock;
            DNekScalMatSharedPtr    localMat = m_schurCompl->GetBlock(n,n);
            unsigned int nbdry    = localMat->GetRows();
            unsigned int nblks    = 1;
            unsigned int esize[1] = {nbdry};

            schurComplBlock = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nblks, nblks, esize, esize);
            schurComplBlock->SetBlock(0, 0, localMat);

            return schurComplBlock;
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysIterativeStaticCond::v_AssembleSchurComplement(
            const AssemblyMapSharedPtr pLocToGloMap)
        {
            boost::ignore_unused(pLocToGloMap);
            // Set up unique map
            v_UniqueMap();

            // Build precon again if we in multi-level static condensation (a
            // bit of a hack)
            if (m_linSysKey.GetGlobalSysSolnType() ==
                    eIterativeMultiLevelStaticCond)
            {
                m_precon = CreatePrecon(m_locToGloMap.lock());
                m_precon->BuildPreconditioner();
            }

            PrepareLocalSchurComplement();
        }


        /**
         * Populates sparse block-diagonal schur complement matrix from
         * the block matrices stored in #m_blkMatrices.
         */
        void GlobalLinSysIterativeStaticCond::PrepareLocalSchurComplement()
        {
            LocalMatrixStorageStrategy storageStrategy =
                m_expList.lock()->GetSession()->
                    GetSolverInfoAsEnum<LocalMatrixStorageStrategy>(
                                       "LocalMatrixStorageStrategy");

            switch(storageStrategy)
            {
                case MultiRegions::eContiguous:
                case MultiRegions::eNonContiguous:
                {
                    size_t storageSize = 0;
                    int nBlk           = m_schurCompl->GetNumberOfBlockRows();

                    m_scale = Array<OneD, NekDouble> (nBlk, 1.0);
                    m_rows  = Array<OneD, unsigned int> (nBlk, 0U);

                    // Determine storage requirements for dense blocks.
                    for (int i = 0; i < nBlk; ++i)
                    {
                        m_rows[i]    = m_schurCompl->GetBlock(i,i)->GetRows();
                        m_scale[i]   = m_schurCompl->GetBlock(i,i)->Scale();
                        storageSize += m_rows[i] * m_rows[i];
                    }

                    // Assemble dense storage blocks.
                    DNekScalMatSharedPtr loc_mat;
                    m_denseBlocks.resize(nBlk);
                    double *ptr = 0;

                    if (MultiRegions::eContiguous == storageStrategy)
                    {
                        m_storage.resize    (storageSize);
                        ptr = &m_storage[0];
                    }

                    for (unsigned int n = 0; n < nBlk; ++n)
                    {
                        loc_mat = m_schurCompl->GetBlock(n,n);

                        if (MultiRegions::eContiguous == storageStrategy)
                        {
                            int loc_lda      = loc_mat->GetRows();
                            int blockSize    = loc_lda * loc_lda;
                            m_denseBlocks[n] = ptr;
                            for(int i = 0; i < loc_lda; ++i)
                            {
                                for(int j = 0; j < loc_lda; ++j)
                                {
                                    ptr[j*loc_lda+i] = (*loc_mat)(i,j);
                                }
                            }
                            ptr += blockSize;
                            GlobalLinSys::v_DropStaticCondBlock(n);
                        }
                        else
                        {
                            m_denseBlocks[n] = loc_mat->GetRawPtr();
                        }
                    }
                    break;
                }
                case MultiRegions::eSparse:
                {
                    DNekScalMatSharedPtr loc_mat;
                    int loc_lda;
                    int blockSize = 0;
                    
                    // First run through to split the set of local matrices into
                    // partitions of fixed block size, and count number of local
                    // matrices that belong to each partition.
                    std::vector<std::pair<int,int> > partitions;
                    for(int n = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
                    {
                        loc_mat = m_schurCompl->GetBlock(n,n);
                        loc_lda = loc_mat->GetRows();

                        ASSERTL1(loc_lda >= 0,
                                 boost::lexical_cast<std::string>(n) + "-th "
                                 "matrix block in Schur complement has "
                                 "rank 0!");

                        if (blockSize == loc_lda)
                        {
                            partitions[partitions.size()-1].first++;
                        }
                        else
                        {
                            blockSize = loc_lda;
                            partitions.push_back(make_pair(1,loc_lda));
                        }
                    }

                    MatrixStorage matStorage = eFULL;

                    // Create a vector of sparse storage holders
                    DNekSmvBsrDiagBlkMat::SparseStorageSharedPtrVector
                            sparseStorage (partitions.size());

                    for (int part = 0, n = 0; part < partitions.size(); ++part)
                    {
                        BCOMatType partMat;

                        for(int k = 0; k < partitions[part].first; ++k, ++n)
                        {
                            loc_mat = m_schurCompl->GetBlock(n,n);
                            loc_lda = loc_mat->GetRows();

                            ASSERTL1(loc_lda == partitions[part].second,
                                     boost::lexical_cast<std::string>(n) + "-th"
                                     " matrix block in Schur complement has "
                                     "unexpected rank");

                            NekDouble scale = loc_mat->Scale();
                            if(fabs(scale-1.0) > NekConstants::kNekZeroTol)
                            {
                                Array<OneD, NekDouble>  matarray(loc_lda*loc_lda);
                                Vmath::Smul(loc_lda*loc_lda,scale,
                                            loc_mat->GetRawPtr(),1,&matarray[0],1);
                                partMat[make_pair(k,k)] = BCOEntryType(matarray);
                            }
                            else // scale factor is 1.0
                            {
                                partMat[make_pair(k,k)] = BCOEntryType(
                                loc_lda*loc_lda, loc_mat->GetRawPtr());
                            }

                            GlobalLinSys::v_DropStaticCondBlock(n);
                        }

                        sparseStorage[part] =
                        MemoryManager<DNekSmvBsrDiagBlkMat::StorageType>::
                            AllocateSharedPtr(
                                partitions[part].first, partitions[part].first,
                                partitions[part].second, partMat, matStorage );
                    }

                    // Create block diagonal matrix
                    m_sparseSchurCompl = MemoryManager<DNekSmvBsrDiagBlkMat>::
                                            AllocateSharedPtr(sparseStorage);

                    break;
                }
                default:
                    ErrorUtil::NekError("Solver info property \
                        LocalMatrixStorageStrategy takes values \
                        Contiguous, Non-contiguous and Sparse");
            }
        }

        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nLocal = m_locToGloMap.lock()->GetNumLocalBndCoeffs();
            AssemblyMapSharedPtr asmMap = m_locToGloMap.lock();

            if (m_sparseSchurCompl)
            {
                // Do matrix multiply locally using block-diagonal sparse matrix
                Array<OneD, NekDouble> tmp = m_wsp + nLocal;

                asmMap->GlobalToLocalBnd(pInput, m_wsp);
                m_sparseSchurCompl->Multiply(m_wsp,tmp);
                asmMap->AssembleBnd(tmp, pOutput);
            }
            else
            {
                // Do matrix multiply locally, using direct BLAS calls
                asmMap->GlobalToLocalBnd(pInput, m_wsp);
                int i, cnt;
                Array<OneD, NekDouble> tmpout = m_wsp + nLocal;
                for (i = cnt = 0; i < m_denseBlocks.size(); cnt += m_rows[i], ++i)
                {
                    const int rows = m_rows[i];
                    Blas::Dgemv('N', rows, rows,
                                m_scale[i], m_denseBlocks[i], rows,
                                m_wsp.get()+cnt, 1,
                                0.0, tmpout.get()+cnt, 1);
                }
                asmMap->AssembleBnd(tmpout, pOutput);
            }
        }

        void GlobalLinSysIterativeStaticCond::v_UniqueMap()
        {
            m_map = m_locToGloMap.lock()->GetGlobalToUniversalBndMapUnique();
        }

        void GlobalLinSysIterativeStaticCond::v_PreSolve(
            int                      scLevel,
            Array<OneD, NekDouble>   &F_bnd)
        {
            if (scLevel == 0)
            {
                // When matrices are supplied to the constructor at the top
                // level, the preconditioner is never set up.
                if (!m_precon)
                {
                    m_precon = CreatePrecon(m_locToGloMap.lock());
                    m_precon->BuildPreconditioner();
                }
                
#if 1 // to be consistent with original 

                int nGloBndDofs = m_locToGloMap.lock()->GetNumGlobalBndCoeffs();
                Array<OneD, NekDouble> F_gloBnd(nGloBndDofs);
                NekVector<NekDouble> F_GloBnd(nGloBndDofs,F_gloBnd,eWrapper);
                
                m_locToGloMap.lock()->AssembleBnd(F_bnd,F_gloBnd);
                Set_Rhs_Magnitude(F_GloBnd);
                
#else
                int nLocBndDofs   = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

                //Set_Rhs_Magnitude - version using local array

                Array<OneD, NekDouble> vExchange(1, 0.0);

                vExchange[0] += Blas::Ddot(nLocBndDofs, F_bnd,1,F_bnd,1);
                m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

                // To ensure that very different rhs values are not being
                // used in subsequent solvers such as the velocit solve in
                // INC NS. If this works we then need to work out a better
                // way to control this.
                NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;
                
                if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
                {
                    m_rhs_magnitude = new_rhs_mag;
                }
                else
                {
                    m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) + 
                                       (1.0-m_rhs_mag_sm)*new_rhs_mag); 
                }
#endif
            }
            else
            {
                // for multilevel iterative solver always use rhs
                // vector value with no weighting
                m_rhs_magnitude = NekConstants::kNekUnsetDouble;
            }
        }
        
        void GlobalLinSysIterativeStaticCond::v_BasisFwdTransform(
                                     Array<OneD, NekDouble>& pInOut)
        {
            m_precon->DoTransformBasisToLowEnergy(pInOut);            
        }

        void GlobalLinSysIterativeStaticCond::v_CoeffsBwdTransform(
            Array<OneD, NekDouble>& pInOut)
        {
	    m_precon->DoTransformCoeffsFromLowEnergy(pInOut);
        }

        void GlobalLinSysIterativeStaticCond::v_CoeffsFwdTransform(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
	    m_precon->DoTransformCoeffsToLowEnergy(pInput,pOutput);
        }
        
        GlobalLinSysStaticCondSharedPtr GlobalLinSysIterativeStaticCond::v_Recurse(
            const GlobalLinSysKey                &mkey,
            const std::weak_ptr<ExpList>         &pExpList,
            const DNekScalBlkMatSharedPtr         pSchurCompl,
            const DNekScalBlkMatSharedPtr         pBinvD,
            const DNekScalBlkMatSharedPtr         pC,
            const DNekScalBlkMatSharedPtr         pInvD,
            const std::shared_ptr<AssemblyMap>   &l2gMap)
        {
            GlobalLinSysIterativeStaticCondSharedPtr sys = MemoryManager<
                GlobalLinSysIterativeStaticCond>::AllocateSharedPtr(
                    mkey, pExpList, pSchurCompl, pBinvD, pC, pInvD, l2gMap,
                    m_precon);
            sys->Initialise(l2gMap);
            return sys;
        }
    }
}
