/*
 * LocalMatrixSystem.cpp
 *
 *  Created on: 21 Nov 2010
 *      Author: cc
 */

#include "LocalRegions/MatrixKey.h"
#include <MultiRegions/LocalToGlobalC0ContMap.h>
#include "MultiRegions/LocalMatrixSystem.h"

namespace Nektar
{
    namespace MultiRegions
    {
        LocalMatrixSystem::LocalMatrixSystem()
        {

        }

        LocalMatrixSystem::LocalMatrixSystem(const GlobalLinSysKey& mkey,
                    boost::shared_ptr<StdRegions::StdExpansionVector>& exp,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                    const Array<OneD, int>& pOffsets,
                    const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo)
        {

        }

        Array<OneD, DNekScalBlkMatSharedPtr> LocalMatrixSystem::GetLocalSystem()
        {
            return m_sys;
        }

        void LocalMatrixSystem::DeallocateSystem(unsigned int pIndex)
        {
            ASSERTL0(pIndex < m_sys.num_elements(), "Index out of range.");
            m_sys[pIndex].reset();
        }

        // FULL MATRIX SYSTEM
        LocalMatrixSystemFull::LocalMatrixSystemFull(
                    const GlobalLinSysKey& pMatrixKey,
                    boost::shared_ptr<StdRegions::StdExpansionVector>& pExp,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                    const Array<OneD, int>& pOffsets,
                    const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo)
                : LocalMatrixSystem(pMatrixKey, pExp, pLocToGloMap, pOffsets, pRobinBCInfo)
        {
            int n,j;
            int cnt1;
            int n_exp = (*pExp).size();
            int vNPoints = 0;

            Array<OneD, unsigned int> nCoeffsPerElmt(n_exp);

            for(j = 0; j < n_exp; j++)
            {
                nCoeffsPerElmt[j] = (*pExp)[pOffsets[j]]->GetNcoeffs();
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            m_sys = Array<OneD, DNekScalBlkMatSharedPtr>(1);
            m_sys[0] = MemoryManager<DNekScalBlkMat>::
                AllocateSharedPtr(nCoeffsPerElmt,nCoeffsPerElmt,blkmatStorage);

            DNekScalMatSharedPtr loc_mat;

            int nel;
            int nvarcoeffs = pMatrixKey.GetNvariableCoefficients();

            for(n = cnt1 = 0; n < n_exp; ++n)
            {
                nel = pOffsets[n];

                // need to be initialised with zero size for non variable coefficient case
                Array<OneD, Array<OneD,const NekDouble> > varcoeffs;

                if(nvarcoeffs>0)
                {
                    varcoeffs = Array<OneD, Array<OneD,const NekDouble> >(nvarcoeffs);
                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs[j] = pMatrixKey.GetVariableCoefficient(j) + cnt1;
                    }

                    cnt1  +=  (*pExp)[n]->GetTotPoints();
                }

                LocalRegions::MatrixKey matkey(pMatrixKey.GetMatrixType(),
                                               (*pExp)[nel]->DetExpansionType(),
                                               *(*pExp)[nel],
                                               pMatrixKey.GetConstants(),
                                               varcoeffs);

                loc_mat = (*pExp)[nel]->GetLocMatrix(matkey);

                if(pRobinBCInfo.count(nel) != 0) // add robin mass matrix
                {
                    RobinBCInfoSharedPtr rBC;

                    // declare local matrix from scaled matrix.
                    int rows = loc_mat->GetRows();
                    int cols = loc_mat->GetColumns();
                    const NekDouble *dat = loc_mat->GetRawPtr();
                    DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,dat);
                    Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                    // add local matrix contribution
                    for(rBC = pRobinBCInfo.find(nel)->second;rBC; rBC = rBC->next)
                    {
                        (*pExp)[nel]->AddRobinMassMatrix(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs,new_mat);
                    }

                    NekDouble one = 1.0;
                    // redeclare loc_mat to point to new_mat plus the scalar.
                    loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,new_mat);
                }

                m_sys[0]->SetBlock(n,n,loc_mat);
            }
        }


        // FULL MATRIX SYSTEM
        LocalMatrixSystemStaticCond::LocalMatrixSystemStaticCond(
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap)
                : LocalMatrixSystem()
        {
            const Array<OneD,const unsigned int>& nbdry_size = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size  = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();
            MatrixStorage blkmatStorage = eDIAGONAL;

            m_sys    = Array<OneD, DNekScalBlkMatSharedPtr>(4);
            m_sys[0] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_sys[1] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_sys[2] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_sys[3] = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);
        }

        LocalMatrixSystemStaticCond::LocalMatrixSystemStaticCond(
                    const GlobalLinSysKey& pLinSysKey,
                    boost::shared_ptr<StdRegions::StdExpansionVector>& pExp,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                    const Array<OneD, int>& pOffsets,
                    const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo)
                : LocalMatrixSystem(pLinSysKey, pExp, pLocToGloMap, pOffsets, pRobinBCInfo)
        {
            StdRegions::MatrixType linsystype = pLinSysKey.GetMatrixType();
            switch (linsystype)
            {
            case StdRegions::eHybridDGHelmBndLam:
                SetupBnd(pLinSysKey, pExp, pLocToGloMap, pOffsets, pRobinBCInfo);
                break;
            default:
                Setup(pLinSysKey, pExp, pLocToGloMap, pOffsets, pRobinBCInfo);
                break;
            }
        }


        void LocalMatrixSystemStaticCond::Setup(const GlobalLinSysKey &pLinSysKey,
                boost::shared_ptr<StdRegions::StdExpansionVector>& pExp,
                const boost::shared_ptr<LocalToGlobalBaseMap>
                                                        &pLocToGloMap,
                const Array<OneD, int>& pOffsets,
                const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo)
        {
            int n,j;
            int cnt1;

            // Setup Block Matrix systems
            int n_exp = (*pExp).size();
            const Array<OneD,const unsigned int>& nbdry_size = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size  = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            MatrixStorage blkmatStorage = eDIAGONAL;
            m_sys      = Array<OneD, DNekScalBlkMatSharedPtr>(4);
            m_sys[0]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_sys[1]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_sys[2]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_sys[3]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            DNekScalBlkMatSharedPtr SchurCompl  = m_sys[0];
            DNekScalBlkMatSharedPtr BinvD       = m_sys[1];
            DNekScalBlkMatSharedPtr C           = m_sys[2];
            DNekScalBlkMatSharedPtr invD        = m_sys[3];

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    tmp_mat;

            int eid;
            int nvarcoeffs = pLinSysKey.GetNvariableCoefficients();

            //map<int, RobinBCInfoSharedPtr> RobinBCInfo = GetRobinBCInfo();

            for(n = cnt1 = 0; n < n_exp; ++n)
            {
                eid = pOffsets[n];

                // need to be initialised with zero size for non variable coefficient case
                Array<OneD, Array<OneD,const NekDouble> > varcoeffs;

                // set up elemental coefficient if necessary
                if(nvarcoeffs>0)
                {
                    varcoeffs = Array<OneD, Array<OneD,const NekDouble> > (nvarcoeffs);

                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs[j] = pLinSysKey.GetVariableCoefficient(j) + cnt1;
                    }
                    cnt1  += (*pExp)[n]->GetTotPoints();
                }

                LocalRegions::MatrixKey matkey(pLinSysKey.GetMatrixType(),
                                               (*pExp)[eid]->DetExpansionType(),
                                               *(*pExp)[eid],
                                               pLinSysKey.GetConstants(),
                                               varcoeffs);

                loc_mat = (*pExp)[eid]->GetLocStaticCondMatrix(matkey);

                if(pRobinBCInfo.count(eid) != 0) // add robin mass matrix
                {
                    RobinBCInfoSharedPtr rBC;

                    tmp_mat = loc_mat->GetBlock(0,0);

                    // declare local matrix from scaled matrix.
                    int rows = tmp_mat->GetRows();
                    int cols = tmp_mat->GetColumns();
                    const NekDouble *dat = tmp_mat->GetRawPtr();
                    DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,dat);
                    Blas::Dscal(rows*cols,tmp_mat->Scale(),new_mat->GetRawPtr(),1);

                    // add local matrix contribution
                    for(rBC = pRobinBCInfo.find(eid)->second;rBC; rBC = rBC->next)
                    {
                        (*pExp)[eid]->AddRobinMassMatrix(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs,new_mat);
                    }

                    NekDouble one = 1.0;
                    // redeclare loc_mat to point to new_mat plus the scalar.
                    tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,new_mat);
                    loc_mat->SetBlock(0,0,tmp_mat);
                }

                SchurCompl->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,0));
                BinvD     ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(0,1));
                C         ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,0));
                invD      ->SetBlock(n,n, tmp_mat = loc_mat->GetBlock(1,1));
            }
        }


        void LocalMatrixSystemStaticCond::SetupBnd(const GlobalLinSysKey &mkey,
                    boost::shared_ptr<StdRegions::StdExpansionVector>& pExp,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &pLocToGloMap,
                    const Array<OneD, int>& pOffsets,
                    const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo)
        {
            StdRegions::MatrixType linsystype = mkey.GetMatrixType();
            ASSERTL0(linsystype == StdRegions::eHybridDGHelmBndLam,
                     "Routine currently only tested for HybridDGHelmholtz");
            ASSERTL1(mkey.GetGlobalSysSolnType()!=eDirectFullMatrix,
                     "This BndLinSys cannot be constructed in case of a full matrix global solve");
            ASSERTL1(mkey.GetGlobalSysSolnType()==pLocToGloMap->GetGlobalSysSolnType(),
                     "The local to global map is not set up for the requested solution type");

            // We will set up this matrix as a statically condensed system
            // where the interior blocks are zero
            int n,j;
            int cnt1;

            NekDouble factor1, factor2;

            // Setup Block Matrix systems
            int n_exp = (*pExp).size();
            const Array<OneD,const unsigned int>& nbdry_size = pLocToGloMap->GetNumLocalBndCoeffsPerPatch();
            const Array<OneD,const unsigned int>& nint_size  = pLocToGloMap->GetNumLocalIntCoeffsPerPatch();

            MatrixStorage blkmatStorage = eDIAGONAL;
            m_sys      = Array<OneD, DNekScalBlkMatSharedPtr>(4);
            m_sys[0]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nbdry_size, blkmatStorage);
            m_sys[1]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nbdry_size, nint_size , blkmatStorage);
            m_sys[2]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nbdry_size, blkmatStorage);
            m_sys[3]   = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(nint_size , nint_size , blkmatStorage);

            DNekScalBlkMatSharedPtr SchurCompl  = m_sys[0];
            DNekScalBlkMatSharedPtr BinvD       = m_sys[1];
            DNekScalBlkMatSharedPtr C           = m_sys[2];
            DNekScalBlkMatSharedPtr invD        = m_sys[3];

            DNekScalMatSharedPtr loc_mat;

            int eid;
            int totnq, nvarcoeffs = mkey.GetNvariableCoefficients();

            for(n = cnt1 = 0; n < n_exp; ++n)
            {
                eid = pOffsets[n];
                totnq = (*pExp)[eid]->GetCoordim()*( (*pExp)[eid]->GetTotPoints() );

                // need to be initialised with zero size for non variable coefficient case
                Array<OneD, Array<OneD,const NekDouble> > varcoeffs;
                Array<OneD, NekDouble> varcoeffs_wk;

                if(nvarcoeffs>0)
                {
                    varcoeffs = Array<OneD, Array<OneD,const NekDouble> > (nvarcoeffs);

                    // When two varcoeffs in a specific order
                    for(j = 0; j < nvarcoeffs; j++)
                    {
                        varcoeffs_wk = Array<OneD, NekDouble>(totnq,0.0);
                        Vmath::Vcopy(totnq, &(mkey.GetVariableCoefficient(j))[cnt1], 1, &varcoeffs_wk[0],1);
                        varcoeffs[j] = varcoeffs_wk;
                    }

                    cnt1  += totnq;
                }

                int Nconstants = mkey.GetNconstants();

                if(Nconstants>2)
                {
                    factor1 = mkey.GetConstant(eid);
                    factor2 = mkey.GetConstant(Nconstants-1);
                }

                else
                {
                    factor1 = mkey.GetConstant(0);
                    factor2 = mkey.GetConstant(1);
                }

                LocalRegions::MatrixKey matkey(linsystype,
                                               (*pExp)[eid]->DetExpansionType(),
                                               *(*pExp)[eid],factor1,factor2,varcoeffs);

                loc_mat = (*pExp)[eid]->GetLocMatrix(matkey);

                if(pRobinBCInfo.count(eid) != 0) // add robin mass matrix
                {
                    RobinBCInfoSharedPtr rBC;

                    // declare local matrix from scaled matrix.
                    int rows = loc_mat->GetRows();
                    int cols = loc_mat->GetColumns();
                    const NekDouble *dat = loc_mat->GetRawPtr();
                    DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,dat);
                    Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                    // add local matrix contribution
                    for(rBC = pRobinBCInfo.find(eid)->second;rBC; rBC = rBC->next)
                    {
                        (*pExp)[eid]->AddRobinMassMatrix(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs,new_mat);
                    }

                    NekDouble one = 1.0;
                    // redeclare loc_mat to point to new_mat plus the scalar.
                    loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,new_mat);
                }

                SchurCompl->SetBlock(n,n,loc_mat);
            }
        }
    }
}

