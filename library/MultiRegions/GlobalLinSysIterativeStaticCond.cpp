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
// Description: Implementation to linear solver using single-
//              or multi-level static condensation
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>

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
            const boost::weak_ptr<ExpList>       &pExpList,
            const boost::shared_ptr<AssemblyMap> &pLocToGloMap)
                : GlobalLinSysIterative(pKey, pExpList, pLocToGloMap),
                  m_locToGloMap (pLocToGloMap)
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
            const GlobalLinSysKey &pKey,
            const boost::weak_ptr<ExpList> &pExpList,
            const DNekScalBlkMatSharedPtr pSchurCompl,
            const DNekScalBlkMatSharedPtr pBinvD,
            const DNekScalBlkMatSharedPtr pC,
            const DNekScalBlkMatSharedPtr pInvD,
            const boost::shared_ptr<AssemblyMap>
            &pLocToGloMap,
            const PreconditionerSharedPtr pPrecon)
            : GlobalLinSysIterative(pKey, pExpList, pLocToGloMap),
              m_schurCompl ( pSchurCompl ),
              m_BinvD      ( pBinvD ),
              m_C          ( pC ),
              m_invD       ( pInvD ),
              m_locToGloMap( pLocToGloMap ),
              m_precon     ( pPrecon )
        {
            // Construct this level
            Initialise(pLocToGloMap);
        }


        void GlobalLinSysIterativeStaticCond::v_InitObject()
        {
            MultiRegions::PreconditionerType pType
                = m_locToGloMap->GetPreconType();
            std::string PreconType
                = MultiRegions::PreconditionerTypeMap[pType];            
            v_UniqueMap();
            m_precon = GetPreconFactory().CreateInstance(
                PreconType,GetSharedThisPtr(),m_locToGloMap);
            
            // Allocate memory for top-level structure
            SetupTopLevel(m_locToGloMap);

            // Construct this level
            Initialise(m_locToGloMap);
        }
        
        /**
         *
         */
        GlobalLinSysIterativeStaticCond::~GlobalLinSysIterativeStaticCond()
        {
            
        }
        
        
        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::v_Solve(
            const Array<OneD, const NekDouble> &in,
                  Array<OneD,       NekDouble> &out,
            const AssemblyMapSharedPtr         &pLocToGloMap,
            const Array<OneD, const NekDouble> &dirForcing)
        {
            bool dirForcCalculated = (bool) dirForcing.num_elements();
            bool atLastLevel       = pLocToGloMap->AtLastLevel();
            int  scLevel           = pLocToGloMap->GetStaticCondLevel();
            
            int nGlobDofs          = pLocToGloMap->GetNumGlobalCoeffs();
            int nGlobBndDofs       = pLocToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = pLocToGloMap->GetNumLocalBndCoeffs();
            int nIntDofs           = pLocToGloMap->GetNumGlobalCoeffs()
                - nGlobBndDofs;
            
            Array<OneD, NekDouble> F = m_wsp + 2*nLocBndDofs;
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
            NekVector<NekDouble> F_GlobBnd(nGlobBndDofs,F,eWrapper);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,0.0);
            NekVector<NekDouble> F_Int(nIntDofs,tmp=F+nGlobBndDofs,eWrapper);
            
            NekVector<NekDouble> V_GlobBnd(nGlobBndDofs,out,eWrapper);
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,
                                              tmp=out+nDirBndDofs,
                                              eWrapper);
            NekVector<NekDouble> V_Int(nIntDofs,tmp=out+nGlobBndDofs,eWrapper);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,m_wsp,eWrapper);
            
            NekVector<NekDouble> V_GlobHomBndTmp(nGlobHomBndDofs,0.0);

            // set up normalisation factor for right hand side on first SC level
            if(scLevel == 0)
            {
                Set_Rhs_Magnitude(F_GlobBnd);
            }

            // Select correct matrix to use at difference levels: top level
            // should use the transformed matrix, all other levels should use
            // original matrix since the transformed matrix is recursively
            // passed down.
            DNekScalBlkMatSharedPtr sc = scLevel == 0 ? m_S1Blk : m_schurCompl;

            if(nGlobHomBndDofs)
            {
                // construct boundary forcing
                if( nIntDofs  && ((!dirForcCalculated) && (atLastLevel)) )
                {
                    DNekScalBlkMat &BinvD      = *m_BinvD;
                    DNekScalBlkMat &SchurCompl = *sc;
                    
                    // include dirichlet boundary forcing 
                    pLocToGloMap->GlobalToLocalBnd(V_GlobBnd,V_LocBnd);
                    V_LocBnd = BinvD*F_Int + SchurCompl*V_LocBnd;
                    
                }
                else if((!dirForcCalculated) && (atLastLevel))
                {
                    // include dirichlet boundary forcing
                    DNekScalBlkMat &SchurCompl = *sc;
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

                
                //transform from original basis to low energy
                Array<OneD, NekDouble> tmp;
                m_precon->DoTransformToLowEnergy(F,nDirBndDofs);

                // For parallel multi-level static condensation some
                // processors may have different levels to others. This
                // routine receives contributions to partition vertices from
                // those lower levels, whilst not sending anything to the
                // other partitions, and includes them in the modified right
                // hand side vector.
                int lcLevel = pLocToGloMap->GetLowestStaticCondLevel();
                if(atLastLevel && scLevel < lcLevel)
                {
                    // If this level is not the lowest level across all
                    // processes, we must do dummy communication for the
                    // remaining levels
                    Array<OneD, NekDouble> tmp(nGlobBndDofs);
                    for (int i = scLevel; i < lcLevel; ++i)
                    {
                        Vmath::Fill(nGlobBndDofs, 0.0, tmp, 1);
                        pLocToGloMap->UniversalAssembleBnd(tmp);
                        Vmath::Vcopy(nGlobHomBndDofs,
                                     tmp.get()+nDirBndDofs,          1,
                                     V_GlobHomBndTmp.GetPtr().get(), 1);
                        F_HomBnd = F_HomBnd - V_GlobHomBndTmp;
                    }
                }

                // solve boundary system
                if(atLastLevel)
                {
                    Array<OneD, NekDouble> pert(nGlobBndDofs,0.0);
                    NekVector<NekDouble>   Pert(nGlobBndDofs,pert,eWrapper);

                    // Solve for difference from initial solution given inout;
                    SolveLinearSystem(
                        nGlobBndDofs, F, pert, pLocToGloMap, nDirBndDofs);

                    //transform back to original basis
                    m_precon->DoTransformFromLowEnergy(pert);

                    // Add back initial conditions onto difference
                    Vmath::Vadd(nGlobHomBndDofs,&out[nDirBndDofs],1,
                                &pert[nDirBndDofs],1,&out[nDirBndDofs],1);
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



        /**
         * If at the last level of recursion (or the only level in the case of
         * single-level static condensation), assemble the Schur complement.
         * For other levels, in the case of multi-level static condensation,
         * the next level of the condensed system is computed.
         * @param   pLocToGloMap    Local to global mapping.
         */
        void GlobalLinSysIterativeStaticCond::Initialise(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int nLocalBnd = m_locToGloMap->GetNumLocalBndCoeffs();
            int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
            m_wsp = Array<OneD, NekDouble>(2*nLocalBnd + nGlobal);

            if(pLocToGloMap->AtLastLevel())
            {
                // decide whether to assemble schur complement globally
                // based on global optimisation parameter to the
                // full system matrix (current operator)
                bool doGlobalOp = m_expList.lock()->GetGlobalOptParam()->
                    DoGlobalMatOp(m_linSysKey.GetMatrixType());

                if(doGlobalOp)
                {
                    AssembleSchurComplement(pLocToGloMap);
                }
                else
                {
                    PrepareLocalSchurComplement();
                }
            }
            else
            {
                ConstructNextLevelCondensedSystem(
                        pLocToGloMap->GetNextLevelLocalToGlobalMap());
            }
        }

        int GlobalLinSysIterativeStaticCond::v_GetNumBlocks()
        {
            return m_schurCompl->GetNumberOfBlockRows();
        }

        DNekScalBlkMatSharedPtr GlobalLinSysIterativeStaticCond::
            v_GetStaticCondBlock(unsigned int n)
        {
            DNekScalBlkMatSharedPtr schurComplBlock;
            int  scLevel           = m_locToGloMap->GetStaticCondLevel();
            DNekScalBlkMatSharedPtr sc = scLevel == 0 ? m_S1Blk : m_schurCompl;
            DNekScalMatSharedPtr    localMat = sc->GetBlock(n,n);
            int nbdry    = localMat->GetRows();
            int nblks    = 1;
            unsigned int esize[1] = {nbdry};

            schurComplBlock = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nblks, nblks, esize, esize);
            schurComplBlock->SetBlock(0, 0, localMat);

            return schurComplBlock;
        }

        /**
         * For the first level in multi-level static condensation, or the only
         * level in the case of single-level static condensation, allocate the
         * condensed matrices and populate them with the local matrices
         * retrieved from the expansion list.
         * @param
         */
        void GlobalLinSysIterativeStaticCond::SetupTopLevel(
                const boost::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            int n;
            int n_exp = m_expList.lock()->GetNumElmts();

            const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            // Setup Block Matrix systems
            MatrixStorage blkmatStorage = eDIAGONAL;
            m_schurCompl = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_BinvD      = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_C          = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_invD       = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            //Original schur complement matrix
            m_S1Blk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);

            for(n = 0; n < n_exp; ++n)
            {
                if (m_linSysKey.GetMatrixType() ==
                        StdRegions::eHybridDGHelmBndLam)
                {
                    DNekScalMatSharedPtr loc_mat
                        = GlobalLinSys::v_GetBlock(
                            m_expList.lock()->GetOffset_Elmt_Id(n));
                    m_schurCompl->SetBlock(n,n,loc_mat);
                    m_S1Blk     ->SetBlock(n,n,loc_mat);
                }
                else
                {
                    DNekScalBlkMatSharedPtr loc_S1
                        = GlobalLinSys::v_GetStaticCondBlock(
                            m_expList.lock()->GetOffset_Elmt_Id(n));
                    DNekScalBlkMatSharedPtr loc_schur
                        = m_precon->TransformedSchurCompl(
                            m_expList.lock()->GetOffset_Elmt_Id(n),loc_S1);

                    DNekScalMatSharedPtr t;
                    m_schurCompl->SetBlock(n, n, t = loc_schur->GetBlock(0,0));
                    m_BinvD     ->SetBlock(n, n, t = loc_S1   ->GetBlock(0,1));
                    m_C         ->SetBlock(n, n, t = loc_S1   ->GetBlock(1,0));
                    m_invD      ->SetBlock(n, n, t = loc_S1   ->GetBlock(1,1));
                    m_S1Blk     ->SetBlock(n, n, t = loc_S1   ->GetBlock(0,0));
                }
            }
        }

        /**
         * Assemble the schur complement matrix from the block matrices stored
         * in #m_blkMatrices and the given local to global mapping information.
         * @param   locToGloMap Local to global mapping information.
         */
        void GlobalLinSysIterativeStaticCond::AssembleSchurComplement(
                    const AssemblyMapSharedPtr &pLocToGloMap)
        {
            int i,j,n,cnt,gid1,gid2;
            NekDouble sign1,sign2;

            int nBndDofs  = pLocToGloMap->GetNumGlobalBndCoeffs();
            int NumDirBCs = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            unsigned int rows = nBndDofs - NumDirBCs;
            unsigned int cols = nBndDofs - NumDirBCs;

            // COO sparse storage to assist in assembly
            COOMatType gmat_coo;

            // Get the matrix storage structure
            // (whether to store only one triangular part, if symmetric)
            MatrixStorage matStorage = eFULL;

            // assemble globally
            DNekScalMatSharedPtr loc_mat;
            int loc_lda;
            for(n = cnt = 0; n < m_schurCompl->GetNumberOfBlockRows(); ++n)
            {
                loc_mat = m_schurCompl->GetBlock(n,n);
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

                            if (gid2 >= 0)
                            {
                                gmat_coo[std::make_pair(gid1,gid2)] +=
                                    sign1*sign2*(*loc_mat)(i,j);
                            }
                        }
                    }
                }
                cnt += loc_lda;
            }

            DNekSmvBsrDiagBlkMat::SparseStorageSharedPtrVector
                sparseStorage (1);

            BCOMatType partMat;
            convertCooToBco(rows, cols, 1, gmat_coo, partMat);

            sparseStorage[0] =
                 MemoryManager<DNekSmvBsrDiagBlkMat::StorageType>::
                    AllocateSharedPtr(rows, cols, 1, partMat, matStorage );

            // Create block diagonal matrix
            m_sparseSchurCompl = MemoryManager<DNekSmvBsrDiagBlkMat>::
                                            AllocateSharedPtr(sparseStorage);

            cout << "global SchurCompl: row density = " 
                 << gmat_coo.size()/cols << endl;
            cout << "global SchurCompl: matrix rows = " 
                 << rows << endl;
            cout << "global SchurCompl: matrix nnzs = " 
                 << gmat_coo.size() << endl;
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

            bool verbose = (m_expList.lock()->GetSession()->
                    DefinesCmdLineArgument("verbose"))? true : false;

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
                            GlobalLinSys::v_DropStaticCondBlock(
                                m_expList.lock()->GetOffset_Elmt_Id(n));
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

                    if (verbose)
                    {
                        cout << "sizes of local matrices in order: " << endl;
                        for (int i = 0; i < partitions.size(); i++)
                        {
                            cout << " (" << partitions[i].first << ", "
                                 << partitions[i].second << ")";
                        }
                        cout << endl;
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

                            partMat[make_pair(k,k)] = BCOEntryType(
                                loc_lda*loc_lda, loc_mat->GetRawPtr());

                            GlobalLinSys::v_DropStaticCondBlock(
                                m_expList.lock()->GetOffset_Elmt_Id(n));
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

                    size_t matBytes, bsruBlockBytes;

                    matBytes      = m_sparseSchurCompl->GetMemoryFootprint();
                    bsruBlockBytes = m_sparseSchurCompl->GetMemoryFootprint(0);

                    if (verbose)
                    {
                        cout << "Local matrix memory, bytes = " << matBytes;
                        if (matBytes/(1024*1024) > 0)
                        {
                            std::cout << " ("<< matBytes/(1024*1024) << " MB)"
                                      << std::endl;
                        }
                        else
                        {
                            std::cout << " ("<< matBytes/1024 << " KB)"
                                      << std::endl;
                        }

                        std::cout << "First BSRU submatrix memory, bytes = "
                                  << bsruBlockBytes;

                        if (bsruBlockBytes/(1024*1024) > 0)
                        {
                            std::cout << " ("<< bsruBlockBytes/(1024*1024)
                                      << " MB)" << std::endl;
                        }
                        else
                        {
                            std::cout << " ("<< bsruBlockBytes/1024
                                      << " KB)" << std::endl;
                        }
                    }
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
        void GlobalLinSysIterativeStaticCond::ConstructNextLevelCondensedSystem(
                        const AssemblyMapSharedPtr& pLocToGloMap)
        {
            int i,j,n,cnt;
            DNekScalBlkMatSharedPtr blkMatrices[4];

            // Create temporary matrices within an inner-local scope to ensure
            // any references to the intermediate storage is lost before
            // the recursive step, rather than at the end of the routine.
            // This allows the schur complement matrix from this level to be
            // disposed of in the next level after use without being retained
            // due to lingering shared pointers.
            {

                const Array<OneD,const unsigned int>& nBndDofsPerPatch
                                = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nIntDofsPerPatch
                                = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

                // STEP 1:
                // Based upon the schur complement of the the current level we
                // will substructure this matrix in the form
                //      --     --
                //      | A   B |
                //      | C   D |
                //      --     --
                // All matrices A,B,C and D are (diagonal) blockmatrices.
                // However, as a start we will use an array of DNekMatrices as
                // it is too hard to change the individual entries of a
                // DNekScalBlkMatSharedPtr.

                // In addition, we will also try to ensure that the memory of
                // the blockmatrices will be contiguous. This will probably
                // enhance the efficiency
                // - Calculate the total number of entries in the blockmatrices
                int nPatches  = pLocToGloMap->GetNumPatches();
                int nEntriesA = 0; int nEntriesB = 0;
                int nEntriesC = 0; int nEntriesD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    nEntriesA += nBndDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesB += nBndDofsPerPatch[i]*nIntDofsPerPatch[i];
                    nEntriesC += nIntDofsPerPatch[i]*nBndDofsPerPatch[i];
                    nEntriesD += nIntDofsPerPatch[i]*nIntDofsPerPatch[i];
                }

                // Now create the DNekMatrices and link them to the memory
                // allocated above
                Array<OneD, DNekMatSharedPtr> substructuredMat[4]
                    = {Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix A
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix B
                       Array<OneD, DNekMatSharedPtr>(nPatches),  //Matrix C
                       Array<OneD, DNekMatSharedPtr>(nPatches)}; //Matrix D

                // Initialise storage for the matrices. We do this separately
                // for each matrix so the matrices may be independently
                // deallocated when no longer required.
                Array<OneD, NekDouble> storageA(nEntriesA,0.0);
                Array<OneD, NekDouble> storageB(nEntriesB,0.0);
                Array<OneD, NekDouble> storageC(nEntriesC,0.0);
                Array<OneD, NekDouble> storageD(nEntriesD,0.0);

                Array<OneD, NekDouble> tmparray;
                PointerWrapper wType = eWrapper;
                int cntA = 0;
                int cntB = 0;
                int cntC = 0;
                int cntD = 0;

                for(i = 0; i < nPatches; i++)
                {
                    // Matrix A
                    tmparray = storageA+cntA;
                    substructuredMat[0][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix B
                    tmparray = storageB+cntB;
                    substructuredMat[1][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nBndDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix C
                    tmparray = storageC+cntC;
                    substructuredMat[2][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nBndDofsPerPatch[i],
                                                        tmparray, wType);
                    // Matrix D
                    tmparray = storageD+cntD;
                    substructuredMat[3][i] = MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nIntDofsPerPatch[i],
                                                        nIntDofsPerPatch[i],
                                                        tmparray, wType);

                    cntA += nBndDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntB += nBndDofsPerPatch[i] * nIntDofsPerPatch[i];
                    cntC += nIntDofsPerPatch[i] * nBndDofsPerPatch[i];
                    cntD += nIntDofsPerPatch[i] * nIntDofsPerPatch[i];
                }

                // Then, project SchurComplement onto
                // the substructured matrices of the next level
                DNekScalBlkMatSharedPtr SchurCompl  = m_schurCompl;
                DNekScalMatSharedPtr schurComplSubMat;
                int       schurComplSubMatnRows;
                Array<OneD, const int>       patchId, dofId;
                Array<OneD, const unsigned int>      isBndDof;
                Array<OneD, const NekDouble> sign;
                NekDouble scale;

                for(n = cnt = 0; n < SchurCompl->GetNumberOfBlockRows(); ++n)
                {
                    schurComplSubMat      = SchurCompl->GetBlock(n,n);
                    schurComplSubMatnRows = schurComplSubMat->GetRows();
                    
                    scale = SchurCompl->GetBlock(n,n)->Scale();
                    Array<OneD, NekDouble> schurSubMat
                        = SchurCompl->GetBlock(n,n)->GetOwnedMatrix()->GetPtr();
                    
                    patchId  = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetPatchId() + cnt;
                    dofId    = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetDofId()   + cnt;
                    isBndDof = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->IsBndDof() + cnt;
                    sign     = pLocToGloMap->GetPatchMapFromPrevLevel()
                               ->GetSign() + cnt;

                    // Set up  Matrix;
                    for(i = 0; i < schurComplSubMatnRows; ++i)
                    {
                        int pId = patchId[i];
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][pId]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][patchId[i]]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][patchId[i]]->GetPtr();
                        Array<OneD, NekDouble> subMat3
                            = substructuredMat[3][patchId[i]]->GetPtr();
                        int subMat0rows = substructuredMat[0][pId]->GetRows();
                        int subMat1rows = substructuredMat[1][pId]->GetRows();
                        int subMat2rows = substructuredMat[2][pId]->GetRows();
                        int subMat3rows = substructuredMat[3][pId]->GetRows();
                        
                        if(isBndDof[i])
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat0[dofId[i]+dofId[j]*subMat0rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat1[dofId[i]+dofId[j]*subMat1rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                            }
                        }
                        else
                        {
                            for(j = 0; j < schurComplSubMatnRows; ++j)
                            {
                                ASSERTL0(patchId[i]==patchId[j],
                                         "These values should be equal");
                                
                                if(isBndDof[j])
                                {
                                    subMat2[dofId[i]+dofId[j]*subMat2rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                                else
                                {
                                    subMat3[dofId[i]+dofId[j]*subMat3rows] +=
                                        sign[i]*sign[j]*(
                                            scale*schurSubMat[
                                                i+j*schurComplSubMatnRows]);
                                }
                            }
                        }
                    }
                    cnt += schurComplSubMatnRows;
                }

                // STEP 2: condense the system
                // This can be done elementally (i.e. patch per patch)
                for(i = 0; i < nPatches; i++)
                {
                    if(nIntDofsPerPatch[i])
                    {
                        Array<OneD, NekDouble> subMat0
                            = substructuredMat[0][i]->GetPtr();
                        Array<OneD, NekDouble> subMat1
                            = substructuredMat[1][i]->GetPtr();
                        Array<OneD, NekDouble> subMat2
                            = substructuredMat[2][i]->GetPtr();
                        int subMat0rows = substructuredMat[0][i]->GetRows();
                        int subMat1rows = substructuredMat[1][i]->GetRows();
                        int subMat2rows = substructuredMat[2][i]->GetRows();
                        int subMat2cols = substructuredMat[2][i]->GetColumns();

                        // 1. D -> InvD
                        substructuredMat[3][i]->Invert();
                        // 2. B -> BInvD
                        (*substructuredMat[1][i]) = (*substructuredMat[1][i])*
                            (*substructuredMat[3][i]);
                        // 3. A -> A - BInvD*C (= schurcomplement)
                        // (*substructuredMat[0][i]) =(*substructuredMat[0][i])-
                        // (*substructuredMat[1][i])*(*substructuredMat[2][i]);
                        // Note: faster to use blas directly
                        Blas::Dgemm('N','N', subMat1rows, subMat2cols,
                                    subMat2rows, -1.0, &subMat1[0], subMat1rows,
                                    &subMat2[0], subMat2rows, 1.0, &subMat0[0],
                                    subMat0rows);
                    }
                }

                // STEP 3: fill the block matrices. However, do note that we
                // first have to convert them to a DNekScalMat in order to be
                // compatible with the first level of static condensation
                const Array<OneD,const unsigned int>& nbdry_size
                    = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
                const Array<OneD,const unsigned int>& nint_size
                    = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
                MatrixStorage blkmatStorage = eDIAGONAL;

                blkMatrices[0] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
                blkMatrices[1] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
                blkMatrices[2] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
                blkMatrices[3] = MemoryManager<DNekScalBlkMat>
                    ::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

                DNekScalMatSharedPtr tmpscalmat;
                for(i = 0; i < nPatches; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        tmpscalmat= MemoryManager<DNekScalMat>
                            ::AllocateSharedPtr(1.0,substructuredMat[j][i]);
                        blkMatrices[j]->SetBlock(i,i,tmpscalmat);
                    }
                }
            }

            // We've finished with the Schur complement matrix passed to this
            // level, so return the memory to the system. The Schur complement
            // matrix need only be retained at the last level. Save the other
            // matrices at this level though.
            m_schurCompl.reset();

            m_recursiveSchurCompl =
                MemoryManager<GlobalLinSysIterativeStaticCond>
                ::AllocateSharedPtr(
                    m_linSysKey, m_expList, blkMatrices[0], blkMatrices[1],
                    blkMatrices[2], blkMatrices[3], pLocToGloMap, m_precon);
        }

        /**
         *
         */
        void GlobalLinSysIterativeStaticCond::v_DoMatrixMultiply(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
            int nDir = m_locToGloMap->GetNumGlobalDirBndCoeffs();


            bool doGlobalOp = m_expList.lock()->GetGlobalOptParam()->
                    DoGlobalMatOp(m_linSysKey.GetMatrixType());

            if(doGlobalOp)
            {
                // Do matrix multiply globally
                Array<OneD, NekDouble> in  = pInput + nDir;
                Array<OneD, NekDouble> out = pOutput+ nDir;

                m_sparseSchurCompl->Multiply(in,out);
                m_locToGloMap->UniversalAssembleBnd(pOutput, nDir);
            }
            else if (m_sparseSchurCompl)
            {
                // Do matrix multiply locally using block-diagonal sparse matrix
                Array<OneD, NekDouble> tmp = m_wsp + nLocal;

                m_locToGloMap->GlobalToLocalBnd(pInput, m_wsp);
                m_sparseSchurCompl->Multiply(m_wsp,tmp);
                m_locToGloMap->AssembleBnd(tmp, pOutput);
            }
            else
            {
                // Do matrix multiply locally, using direct BLAS calls
                m_locToGloMap->GlobalToLocalBnd(pInput, m_wsp);
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
                m_locToGloMap->AssembleBnd(tmpout, pOutput);
            }
        }

        void GlobalLinSysIterativeStaticCond::v_UniqueMap()
        {
            m_map = m_locToGloMap->GetGlobalToUniversalBndMapUnique();
        }
    }
}
