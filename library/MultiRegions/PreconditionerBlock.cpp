///////////////////////////////////////////////////////////////////////////////
//
// File PreconditionerBlock.cpp
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
// Description: Block Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerBlock.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/SegExp.h>
#include <cmath>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */
        string PreconditionerBlock::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "Block",
                    PreconditionerBlock::create,
                    "Block Diagonal Preconditioning");
       /**
        * @class Block Preconditioner
        *
        * This class implements Block preconditioning for the conjugate
        * gradient matrix solver.
        */

        PreconditionerBlock::PreconditionerBlock(
            const std::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }

        void PreconditionerBlock::v_InitObject()
        {
            GlobalSysSolnType solvertype =
                m_locToGloMap.lock()->GetGlobalSysSolnType();
            ASSERTL0(solvertype == MultiRegions::eIterativeStaticCond ||
                     solvertype == MultiRegions::ePETScStaticCond,
                     "Solver type not valid");
        }

        void PreconditionerBlock::v_BuildPreconditioner()
        {
            GlobalLinSysKey key = m_linsys.lock()->GetKey();

            // Different setup for HDG and CG.
            if (key.GetMatrixType() == StdRegions::eHybridDGHelmBndLam)
            {
                BlockPreconditionerHDG();
            }
            else
            {
                BlockPreconditionerCG();
            }
        }

        /**
         * \brief Construct a block preconditioner from \f$\mathbf{S}_{1}\f$ for
         * the continuous Galerkin system.
         *
         * The preconditioner is defined as:
         *
         * \f[\mathbf{M}^{-1}=\left[\begin{array}{ccc}
         *  \mathrm{Diag}[(\mathbf{S_{1}})_{vv}] & & \\
         *  & (\mathbf{S}_{1})_{eb} & \\
         *  & & (\mathbf{S}_{1})_{fb} \end{array}\right] \f]
         *
         * where \f$\mathbf{S}_{1}\f$ is the local Schur complement matrix for
         * each element and the subscript denotes the portion of the Schur
         * complement associated with the vertex, edge and face blocks
         * respectively.
         */
        void PreconditionerBlock::BlockPreconditionerCG()
        {
            ExpListSharedPtr expList = m_linsys.lock()->GetLocMat().lock();
            LocalRegions::ExpansionSharedPtr exp;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            int i, j, k, n, cnt, gId;
            int meshVertId, meshEdgeId, meshFaceId;

            auto asmMap = m_locToGloMap.lock();

            const int nExp = expList->GetExpSize();
            const int nDirBnd = asmMap->GetNumGlobalDirBndCoeffs();

            // Grab periodic geometry information.
            PeriodicMap periodicVerts, periodicEdges, periodicFaces;
            expList->GetPeriodicEntities(
                periodicVerts, periodicEdges, periodicFaces);

            // The vectors below are of size 3 to have separate storage for
            // vertices, edges and faces.

            // Maps from geometry ID to the matrix representing the extracted
            // portion of S_1. For example idMats[2] folds the S_1 face blocks.
            vector<map<int, vector<NekDouble> > > idMats(3);

            // Maps from the global ID, as obtained from AssemblyMapCG's
            // localToGlobalMap, to the geometry ID.
            vector<map<int, int> > gidMeshIds(3);

            // Maps from the global ID to the number of degrees of freedom for
            // this geometry object.
            vector<map<int, int> > gidDofs(3);

            // Array containing maximum information needed for the universal
            // numbering later. For i = 0,1,2 for each geometry dimension:
            //
            // maxVertIds[2*i]   = maximum geometry ID at dimension i
            // maxVertIds[2*i+1] = maximum number of degrees of freedom for all
            //                     elements of dimension i.
            Array<OneD, int> maxVertIds(6, -1);

            // Figure out mapping from each elemental contribution to offset in
            // (vert,edge,face) triples.
            for (cnt = n = 0; n < nExp; ++n)
            {
                exp = expList->GetExp(n);

                // Grab reference to local Schur complement matrix.
                DNekScalMatSharedPtr schurMat =
                    m_linsys.lock()->GetStaticCondBlock(n)->GetBlock(0,0);

                // Process vertices to extract relevant portion of the Schur
                // complement matrix.
                for (i = 0; i < exp->GetNverts(); ++i)
                {
                    meshVertId = exp->GetGeom()->GetVid(i);
                    int locId = exp->GetVertexMap(i);

                    // Get the global ID of this vertex.
                    gId = asmMap->GetLocalToGlobalMap(
                        cnt + locId) - nDirBnd;

                    // Ignore all Dirichlet vertices.
                    if (gId < 0)
                    {
                        continue;
                    }

                    gidDofs[0][gId] = 1;

                    // Extract vertex value from Schur complement matrix.
                    NekDouble vertVal = (*schurMat)(locId,locId);

                    // See if we have processed this vertex from another
                    // element.
                    auto gIt = idMats[0].find(gId);

                    if (gIt == idMats[0].end())
                    {
                        // If not then put our 'matrix' inside idMats.
                        idMats[0][gId] = vector<NekDouble>(1, vertVal);
                    }
                    else
                    {
                        // Otherwise combine with the value that is already
                        // there (i.e. do assembly on this degree of freedom).
                        gIt->second[0] += vertVal;
                    }

                    // Now check to see if the vertex is periodic. If it is,
                    // then we change meshVertId to be the minimum of all the
                    // other periodic vertices, so that we don't end up
                    // duplicating the matrix in our final block matrix.
                    auto pIt = periodicVerts.find(meshVertId);
                    if (pIt != periodicVerts.end())
                    {
                        for (j = 0; j < pIt->second.size(); ++j)
                        {
                            meshVertId = min(meshVertId, pIt->second[j].id);
                        }
                    }

                    // Finally record the other information we need into the
                    // other matrices.
                    gidMeshIds[0][gId] = meshVertId;
                    maxVertIds[0] = max(maxVertIds[0], meshVertId);
                    maxVertIds[1] = 1;
                }

                // Process edges. This logic is mostly the same as the previous
                // block.
                for (i = 0; i < exp->GetNedges(); ++i)
                {
                    meshEdgeId = exp->GetGeom()->GetEid(i);

                    Array<OneD, unsigned int> bmap, bmap2;
                    Array<OneD, int> sign;
                    StdRegions::Orientation edgeOrient =
                        exp->GetGeom()->GetEorient(i);

                    // Check if this edge is periodic. We may need to flip
                    // orientation if it is.
                    auto pIt = periodicEdges.find(meshEdgeId);
                    if (pIt != periodicEdges.end())
                    {
                        pair<int, StdRegions::Orientation> idOrient =
                            DeterminePeriodicEdgeOrientId(
                                meshEdgeId, edgeOrient, pIt->second);
                        meshEdgeId = idOrient.first;
                        edgeOrient = idOrient.second;
                    }

                    // Grab edge interior map, and the edge inverse boundary
                    // map, so that we can extract this edge from the Schur
                    // complement matrix.
                    exp->GetEdgeInteriorMap(i, edgeOrient, bmap, sign);
                    bmap2 = exp->GetEdgeInverseBoundaryMap(i);

                    // Allocate temporary storage for the extracted edge matrix.
                    const int nEdgeCoeffs = bmap.size();
                    vector<NekDouble> tmpStore(nEdgeCoeffs*nEdgeCoeffs);

                    gId = asmMap->GetLocalToGlobalMap(cnt + bmap[0]);

                    for (j = 0; j < nEdgeCoeffs; ++j)
                    {
                        // We record the minimum ID from the edge for our
                        // maps. This follows the logic that the assembly map
                        // ordering will always give us a contiguous ordering of
                        // global degrees of freedom for edge interior
                        // coefficients.
                        gId = min(gId,
                                  asmMap->GetLocalToGlobalMap(
                                      cnt + bmap[j])
                                  - nDirBnd);

                        // Ignore Dirichlet edges.
                        if (gId < 0)
                        {
                            continue;
                        }

                        const NekDouble sign1 = sign[j];

                        // Extract this edge, along with sign array for assembly
                        // later.
                        for (k = 0; k < nEdgeCoeffs; ++k)
                        {
                            tmpStore[k+j*nEdgeCoeffs] =
                                sign1*sign[k]*(*schurMat)(bmap2[j], bmap2[k]);
                        }
                    }

                    if (gId < 0)
                    {
                        continue;
                    }

                    gidDofs[1][gId] = nEdgeCoeffs;

                    // Assemble this edge matrix with another one, if it exists.
                    auto gIt = idMats[1].find(gId);
                    if (gIt == idMats[1].end())
                    {
                        idMats[1][gId] = tmpStore;
                    }
                    else
                    {
                        ASSERTL1(tmpStore.size() == gIt->second.size(),
                                 "Number of modes mismatch");
                        Vmath::Vadd(nEdgeCoeffs*nEdgeCoeffs, &gIt->second[0], 1,
                                    &tmpStore[0], 1, &gIt->second[0], 1);
                    }

                    gidMeshIds[1][gId] = meshEdgeId;
                    maxVertIds[2] = max(maxVertIds[2], meshEdgeId);
                    maxVertIds[3] = max(maxVertIds[3], nEdgeCoeffs);
                }

                // Process faces. This logic is mostly the same as the previous
                // block.
                for (i = 0; i < exp->GetNfaces(); ++i)
                {
                    meshFaceId = exp->GetGeom()->GetFid(i);

                    Array<OneD, unsigned int> bmap, bmap2;
                    Array<OneD, int> sign;
                    StdRegions::Orientation faceOrient =
                        exp->GetGeom()->GetForient(i);

                    // Check if this face is periodic. We may need to flip
                    // orientation if it is.
                    auto pIt = periodicFaces.find(meshFaceId);
                    if (pIt != periodicFaces.end())
                    {
                        meshFaceId = min(meshFaceId, pIt->second[0].id);
                        faceOrient = DeterminePeriodicFaceOrient(
                            faceOrient, pIt->second[0].orient);
                    }

                    exp->GetFaceInteriorMap(i, faceOrient, bmap, sign);
                    bmap2 = exp->GetFaceInverseBoundaryMap(i);

                    // Allocate temporary storage for the extracted face matrix.
                    const int nFaceCoeffs = bmap.size();
                    vector<NekDouble> tmpStore(nFaceCoeffs*nFaceCoeffs);

                    gId = asmMap->GetLocalToGlobalMap(cnt + bmap[0]);

                    for (j = 0; j < nFaceCoeffs; ++j)
                    {
                        gId = min(gId,
                                  asmMap->GetLocalToGlobalMap(cnt + bmap[j])
                                  - nDirBnd);

                        // Ignore Dirichlet faces.
                        if (gId < 0)
                        {
                            continue;
                        }

                        const NekDouble sign1 = sign[j];

                        // Extract this face, along with sign array for assembly
                        // later.
                        for (k = 0; k < nFaceCoeffs; ++k)
                        {
                            tmpStore[k+j*nFaceCoeffs] =
                                sign1*sign[k]*(*schurMat)(bmap2[j], bmap2[k]);
                        }
                    }

                    if (gId < 0)
                    {
                        continue;
                    }

                    gidDofs[2][gId] = nFaceCoeffs;

                    // Assemble this face matrix with another one, if it exists.
                    auto gIt = idMats[2].find(gId);
                    if (gIt == idMats[2].end())
                    {
                        idMats[2][gId] = tmpStore;
                    }
                    else
                    {
                        ASSERTL1(tmpStore.size() == gIt->second.size(),
                                 "Number of modes mismatch");
                        Vmath::Vadd(nFaceCoeffs*nFaceCoeffs, &gIt->second[0], 1,
                                    &tmpStore[0], 1, &gIt->second[0], 1);
                    }

                    gidMeshIds[2][gId] = meshFaceId;
                    maxVertIds[4] = max(maxVertIds[4], meshFaceId);
                    maxVertIds[5] = max(maxVertIds[5], nFaceCoeffs);
                }

                cnt += exp->GetNcoeffs();
            }

            // Perform a reduction to find maximum vertex, edge and face
            // geometry IDs.
            m_comm = expList->GetSession()->GetComm()->GetRowComm();
            m_comm->AllReduce(maxVertIds, LibUtilities::ReduceMax);

            // Concatenate all matrices into contiguous storage and figure out
            // universal ID numbering.
            vector<NekDouble> storageBuf;
            vector<long> globalToUniversal;

            for (i = 0, cnt = 1; i < 3; ++i)
            {
                const int maxDofs = maxVertIds[2*i+1];

                // Note that iterating over the map uses the property that keys
                // are accessed in order of ascending order, putting everything
                // in the right order for the global system.
                for (auto &gIt : idMats[i])
                {
                    // Copy matrix into storage.
                    storageBuf.insert(storageBuf.end(),
                                      gIt.second.begin(), gIt.second.end());

                    // Get mesh ID from global ID number.
                    ASSERTL1(gidMeshIds[i].count(gIt.first) > 0,
                             "Unable to find global ID " +
                             boost::lexical_cast<string>(gIt.first) +
                             " inside map");
                    meshVertId = gidMeshIds[i][gIt.first];

                    for (j = 0; j < gIt.second.size(); ++j)
                    {
                        globalToUniversal.push_back(
                            cnt + meshVertId*maxDofs*maxDofs + j);
                    }

                    // Free up the temporary storage.
                    gIt.second.clear();
                }

                cnt += (maxVertIds[2*i]+1)*maxDofs*maxDofs;
            }

            ASSERTL1(storageBuf.size() == globalToUniversal.size(),
                     "Storage buffer and global to universal map size does "
                     "not match");

            Array<OneD, NekDouble> storageData(
                storageBuf.size(), &storageBuf[0]);
            Array<OneD, long> globalToUniversalMap(
                globalToUniversal.size(), &globalToUniversal[0]);

            // Use GS to assemble data between processors.
            Gs::gs_data *tmpGs = Gs::Init(
                globalToUniversalMap, m_comm,
                expList->GetSession()->DefinesCmdLineArgument("verbose"));
            Gs::Gather(storageData, Gs::gs_add, tmpGs);

            // Figure out what storage we need in the block matrix.
            Array<OneD, unsigned int> n_blks(
                1 + idMats[1].size() + idMats[2].size());

            // Vertex block is a diagonal matrix.
            n_blks[0] = idMats[0].size();

            // Now extract number of rows in each edge and face block from the
            // gidDofs map.
            cnt = 1;
            for (i = 1; i < 3; ++i)
            {
                for (auto &gIt : idMats[i])
                {
                    ASSERTL1(gidDofs[i].count(gIt.first) > 0,
                             "Unable to find number of degrees of freedom for "
                             "global ID " + boost::lexical_cast<string>(
                                 gIt.first));

                    n_blks[cnt++] = gidDofs[i][gIt.first];
                }
            }

            // Allocate storage for the block matrix.
            m_blkMat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);

            // We deal with the vertex matrix separately since all vertices can
            // be condensed into a single, block-diagonal matrix.
            DNekMatSharedPtr vertMat = MemoryManager<DNekMat>
                ::AllocateSharedPtr(n_blks[0], n_blks[0], 0.0, eDIAGONAL);

            // Fill the vertex matrix with the inverse of each vertex value.
            cnt = 0;
            for (auto gIt = idMats[0].begin(); gIt != idMats[0].end();
                 ++gIt, ++cnt)
            {
                (*vertMat)(cnt, cnt) = 1.0/storageData[cnt];
            }

            // Put the vertex matrix in the block matrix.
            m_blkMat->SetBlock(0,0,vertMat);

            // Finally, grab the matrices from the block storage, invert them
            // and place them in the correct position inside the block matrix.
            int cnt2 = 1;
            for (i = 1; i < 3; ++i)
            {
                for (auto &gIt : idMats[i])
                {
                    int nDofs = gidDofs[i][gIt.first];

                    DNekMatSharedPtr tmp = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nDofs, nDofs);

                    for (j = 0; j < nDofs; ++j)
                    {
                        for (k = 0; k < nDofs; ++k)
                        {
                            (*tmp)(j,k) = storageData[k+j*nDofs + cnt];
                        }
                    }

                    cnt += nDofs*nDofs;

                    tmp->Invert();
                    m_blkMat->SetBlock(cnt2, cnt2, tmp);
                    ++cnt2;
                }
            }
        }

        /**
         * @brief Construct a block preconditioner for the hybridized
         * discontinuous Galerkin system.
         *
         * This system is built in a similar fashion to its continuous variant
         * found in PreconditionerBlock::BlockPreconditionerCG. In this setting
         * however, the matrix is constructed as:
         *
         * \f[ M^{-1} = \mathrm{Diag}[ (\mathbf{S_{1}})_{f}^{-1} ] \f]
         *
         * where each matrix is the Schur complement system restricted to a
         * single face of the trace system.
         */
        void PreconditionerBlock::BlockPreconditionerHDG()
        {
            std::shared_ptr<MultiRegions::ExpList>
                expList = ((m_linsys.lock())->GetLocMat()).lock();
            std::shared_ptr<MultiRegions::ExpList> trace = expList->GetTrace();
            LocalRegions::ExpansionSharedPtr locExpansion;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            AssemblyMapDGSharedPtr asmMap = std::dynamic_pointer_cast<
                AssemblyMapDG>(m_locToGloMap.lock());

            int i, j, k, n, cnt, cnt2;

            // Figure out number of Dirichlet trace elements
            int nTrace = expList->GetTrace()->GetExpSize();
            int nDir   = asmMap->GetNumGlobalDirBndCoeffs();

            for (cnt = n = 0; n < nTrace; ++n)
            {
                if (cnt >= nDir)
                {
                    break;
                }

                cnt += trace->GetExp(n)->GetNcoeffs();
            }

            nDir = n;

            // Allocate storage for block matrix. Need number of unique faces in
            // trace space.
            int maxTraceSize = -1;
            map<int, int> arrayOffsets;

            for (cnt = 0, n = nDir; n < nTrace; ++n)
            {
                int nCoeffs  = trace->GetExp(n)->GetNcoeffs();
                int nCoeffs2 = nCoeffs * nCoeffs;

                arrayOffsets[n]  = cnt;
                cnt             += nCoeffs2;

                if (maxTraceSize < nCoeffs2)
                {
                    maxTraceSize = nCoeffs2;
                }
            }

            // Find maximum trace size.
            m_comm = expList->GetSession()->GetComm()->GetRowComm();
            m_comm->AllReduce(maxTraceSize, LibUtilities::ReduceMax);

            // Zero matrix storage.
            Array<OneD, NekDouble> tmpStore(cnt, 0.0);

            // Assemble block matrices for each trace element.
            for (cnt = n = 0; n < expList->GetExpSize(); ++n)
            {
                int elmt = n;
                locExpansion = expList->GetExp(elmt);

                Array<OneD, LocalRegions::ExpansionSharedPtr> &elmtToTraceMap =
                    asmMap->GetElmtToTrace()[elmt];

                // Block matrix (lambda)
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);
                bnd_mat = loc_mat->GetBlock(0,0);

                int nFacets = locExpansion->GetNumBases() == 2 ?
                    locExpansion->GetNedges() : locExpansion->GetNfaces();

                for (cnt2 = i = 0; i < nFacets; ++i)
                {
                    int nCoeffs = elmtToTraceMap[i]->GetNcoeffs();
                    int elmtId  = elmtToTraceMap[i]->GetElmtId ();

                    // Ignore Dirichlet edges.
                    if (elmtId < nDir)
                    {
                        cnt  += nCoeffs;
                        cnt2 += nCoeffs;
                        continue;
                    }

                    NekDouble *off = &tmpStore[arrayOffsets[elmtId]];

                    for (j = 0; j < nCoeffs; ++j)
                    {
                        NekDouble sign1 = asmMap->GetLocalToGlobalBndSign(
                            cnt + j);
                        for (k = 0; k < nCoeffs; ++k)
                        {
                            NekDouble sign2 = asmMap->GetLocalToGlobalBndSign(
                                cnt + k);
                            off[j*nCoeffs + k] +=
                                (*bnd_mat)(cnt2+j, cnt2+k) * sign1 * sign2;
                        }
                    }

                    cnt  += nCoeffs;
                    cnt2 += nCoeffs;
                }
            }

            // Set up IDs for universal numbering.
            Array<OneD, long> uniIds(tmpStore.size());
            for (cnt = 0, n = nDir; n < nTrace; ++n)
            {
                LocalRegions::ExpansionSharedPtr traceExp = trace->GetExp(n);
                int nCoeffs = traceExp->GetNcoeffs();
                int geomId  = traceExp->GetGeom()->GetGlobalID();

                for (i = 0; i < nCoeffs*nCoeffs; ++i)
                {
                    uniIds[cnt++] = geomId * maxTraceSize + i + 1;
                }
            }

            // Assemble matrices across partitions.
            Gs::gs_data *gsh = Gs::Init(
                uniIds, m_comm,
                expList->GetSession()->DefinesCmdLineArgument("verbose"));
            Gs::Gather(tmpStore, Gs::gs_add, gsh);

            // Set up diagonal block matrix
            Array<OneD, unsigned int> n_blks(nTrace - nDir);
            for (n = 0; n < nTrace - nDir; ++n)
            {
                n_blks[n] = trace->GetExp(n + nDir)->GetNcoeffs();
            }

            m_blkMat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);

            for (n = 0; n < nTrace - nDir; ++n)
            {
                int nCoeffs = trace->GetExp(n + nDir)->GetNcoeffs();
                DNekMatSharedPtr tmp = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nCoeffs, nCoeffs);
                NekDouble *off = &tmpStore[arrayOffsets[n+nDir]];

                for (i = 0; i < nCoeffs; ++i)
                {
                    for (j = 0; j < nCoeffs; ++j)
                    {
                        (*tmp)(i,j) = off[i*nCoeffs + j];
                    }
                }

                tmp->Invert();
                m_blkMat->SetBlock(n, n, tmp);
            }
        }

        /**
         * @brief Apply preconditioner to \p pInput and store the result in \p
         * pOutput.
         */
        void PreconditionerBlock::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nDir    = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();
            int nGlobal = m_locToGloMap.lock()->GetNumGlobalBndCoeffs();
            int nNonDir = nGlobal-nDir;
            DNekBlkMat &M = (*m_blkMat);
            NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
            NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
            z = M * r;
        }
    }
}
