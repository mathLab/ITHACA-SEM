///////////////////////////////////////////////////////////////////////////////
//
// File Preconditioner.cpp
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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerLowEnergy.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalLinSys.h>
#include <LocalRegions/MatrixKey.h>
#include <cmath>

using namespace std;

namespace Nektar
{
    using namespace LibUtilities;

    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */
        string PreconditionerLowEnergy::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "LowEnergyBlock",
                    PreconditionerLowEnergy::create,
                    "LowEnergy Preconditioning");

       /**
         * @class PreconditionerLowEnergy
         *
         * This class implements low energy preconditioning for the conjugate
	 * gradient matrix solver.
	 */

        PreconditionerLowEnergy::PreconditionerLowEnergy(
            const std::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }

        void PreconditionerLowEnergy::v_InitObject()
        {
            GlobalSysSolnType solvertype =
                m_locToGloMap.lock()->GetGlobalSysSolnType();

            ASSERTL0(solvertype == eIterativeStaticCond ||
                     solvertype == ePETScStaticCond, "Solver type not valid");

            std::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();

            m_comm = expList->GetComm();

            LocalRegions::ExpansionSharedPtr locExpansion;

            locExpansion = expList->GetExp(0);

            int nDim = locExpansion->GetShapeDimension();

            ASSERTL0(nDim==3,
                     "Preconditioner type only valid in 3D");

            //Set up block transformation matrix
            SetupBlockTransformationMatrix();

            //Sets up variable p mask 
            CreateVariablePMask();
	}


        /**
	 * \brief Construct the low energy preconditioner from
	 * \f$\mathbf{S}_{2}\f$
	 *
	 * \f[\mathbf{M}^{-1}=\left[\begin{array}{ccc}
	 *  Diag[(\mathbf{S_{2}})_{vv}] & & \\ & (\mathbf{S}_{2})_{eb} & \\ & &
	 *  (\mathbf{S}_{2})_{fb} \end{array}\right] \f]
	 *
	 * where \f$\mathbf{R}\f$ is the transformation matrix and
	 * \f$\mathbf{S}_{2}\f$ the Schur complement of the modified basis,
	 * given by
	 *
	 * \f[\mathbf{S}_{2}=\mathbf{R}\mathbf{S}_{1}\mathbf{R}^{T}\f]
	 *
	 * where \f$\mathbf{S}_{1}\f$ is the local schur complement matrix for
	 * each element.
	 */
       void PreconditionerLowEnergy::v_BuildPreconditioner()
        {
            std::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            LocalRegions::ExpansionSharedPtr locExpansion;
            GlobalLinSysKey linSysKey=(m_linsys.lock())->GetKey();

            int i, j, k;
            int nVerts, nEdges,nFaces;
            int eid, fid, n, cnt, nmodes, nedgemodes, nfacemodes;
            int nedgemodesloc;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol, nCoeffs;

            // Periodic information
            PeriodicMap periodicVerts;
            PeriodicMap periodicEdges;
            PeriodicMap periodicFaces;
            expList->GetPeriodicEntities(periodicVerts,periodicEdges,periodicFaces);

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr    pRSRT;

            DNekMat RS;
            DNekMat RSRT;

            auto asmMap = m_locToGloMap.lock();

            int nDirBnd      = asmMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts = asmMap->GetNumNonDirVertexModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> VertBlockToUniversalMap(nNonDirVerts,-1);

            //maps for different element types
            int n_exp = expList->GetNumElmts();
            int nNonDirEdgeIDs=asmMap->GetNumNonDirEdges();
            int nNonDirFaceIDs=asmMap->GetNumNonDirFaces();

            set<int> edgeDirMap;
            set<int> faceDirMap;
            map<int,int> uniqueEdgeMap;
            map<int,int> uniqueFaceMap;

            //this should be of size total number of local edges + faces
            Array<OneD, int> modeoffset(1 +  nNonDirEdgeIDs + nNonDirFaceIDs,0);
            Array<OneD, int> globaloffset(1 + nNonDirEdgeIDs + nNonDirFaceIDs,0);

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            LocalRegions::Expansion2DSharedPtr bndCondFaceExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>&
		bndConditions = expList->GetBndConditions();

            int meshVertId;
            int meshEdgeId;
            int meshFaceId;

            const Array<OneD, const int> &extradiredges
                = asmMap->GetExtraDirEdges();
            for(i=0; i<extradiredges.size(); ++i)
            {
                meshEdgeId=extradiredges[i];
                edgeDirMap.insert(meshEdgeId);
            }

            //Determine which boundary edges and faces have dirichlet values
            for(i = 0; i < bndCondExp.size(); i++)
            {
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp = std::dynamic_pointer_cast<
                    LocalRegions::Expansion2D>(bndCondExp[i]->GetExp(j));
                    if (bndConditions[i]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                        {
                            meshEdgeId = bndCondFaceExp->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetEid(k);
                            if(edgeDirMap.count(meshEdgeId) == 0)
                            {
                                edgeDirMap.insert(meshEdgeId);
                            }
                        }
                        meshFaceId = bndCondFaceExp->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetGlobalID();
                        faceDirMap.insert(meshFaceId);
                    }
                }
            }

            int dof=0;
            int maxFaceDof=0;
            int maxEdgeDof=0;
            int nlocalNonDirEdges=0;
            int nlocalNonDirFaces=0;
            int ntotalentries=0;

            map<int,int> EdgeSize;
            map<int,int> FaceSize;
            map<int,pair<int,int> >FaceModes;

            /// -  Count  edges, face and add up min edges and min face sizes
            for(n = 0; n < n_exp; ++n)
            {
                locExpansion = expList->GetExp(n);

                nEdges = locExpansion->GetNedges();
                for(j = 0; j < nEdges; ++j)
                {
                    int nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j) - 2;
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetEid(j);
                    if(EdgeSize.count(meshEdgeId) == 0)
                    {
                        EdgeSize[meshEdgeId] = nEdgeInteriorCoeffs;
                    }
                    else
                    {
                        EdgeSize[meshEdgeId] = min(EdgeSize[meshEdgeId],
                                                   nEdgeInteriorCoeffs);
                    }
                }

                nFaces = locExpansion->GetNfaces();
                for(j = 0; j < nFaces; ++j)
                {
                    int nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetFid(j);
                    if(FaceSize.count(meshFaceId) == 0)
                    {
                        FaceSize[meshFaceId] = nFaceInteriorCoeffs;

                        int m0,m1;
                        locExpansion->GetFaceNumModes(j,locExpansion->GetForient(j),m0,m1);
                        FaceModes[meshFaceId] = pair<int,int>(m0,m1);
                    }
                    else
                    {
                        if(nFaceInteriorCoeffs < FaceSize[meshFaceId])
                        {
                            FaceSize[meshFaceId] =  nFaceInteriorCoeffs;
                            int m0,m1;
                            locExpansion->GetFaceNumModes(j,locExpansion->GetForient(j),m0,m1);
                            FaceModes[meshFaceId] = pair<int,int>(m0,m1);
                        }
                    }
                }
            }

            bool verbose =
                expList->GetSession()->DefinesCmdLineArgument("verbose");

            // For parallel runs need to check have minimum of edges and faces over
            // partition boundaries
            if(m_comm->GetSize() > 1)
            {
                int EdgeSizeLen = EdgeSize.size();
                int FaceSizeLen = FaceSize.size();
                Array<OneD, long>      FacetMap(EdgeSizeLen+FaceSizeLen,-1);
                Array<OneD, NekDouble> FacetLen(EdgeSizeLen+FaceSizeLen,-1);

                map<int,int>::iterator it;

                cnt = 0;
                int maxid = 0;
                for(it = EdgeSize.begin(); it!=EdgeSize.end(); ++it,++cnt)
                {
                    FacetMap[cnt] = it->first;
                    maxid = max(it->first,maxid);
                    FacetLen[cnt] = it->second;
                }
                maxid++;

                m_comm->AllReduce(maxid,ReduceMax);

                for(it = FaceSize.begin(); it!=FaceSize.end(); ++it,++cnt)
                {
                    FacetMap[cnt] = it->first + maxid;
                    FacetLen[cnt] = it->second;
                }

                //Exchange vertex data over different processes
                Gs::gs_data *tmp = Gs::Init(FacetMap, m_comm, verbose);
                Gs::Gather(FacetLen, Gs::gs_min, tmp);

                cnt = 0;
                for(it = EdgeSize.begin(); it!=EdgeSize.end(); ++it,++cnt)
                {
                    it->second = (int) FacetLen[cnt];
                }

                for(it = FaceSize.begin(); it!=FaceSize.end(); ++it,++cnt)
                {
                    it->second = (int)FacetLen[cnt];
                }
            }

            // Loop over all the elements in the domain and compute total edge
            // DOF and set up unique ordering.
            map<int,int> nblks;
            int matrixlocation = 0;

            // First do periodic edges
            for (auto &pIt : periodicEdges)
            {
                meshEdgeId = pIt.first;

                if(edgeDirMap.count(meshEdgeId)==0)
                {
                    dof = EdgeSize[meshEdgeId];
                    if(uniqueEdgeMap.count(meshEdgeId)==0 && dof > 0)
                    {
                        bool SetUpNewEdge = true;


                        for (i = 0; i < pIt.second.size(); ++i)
                        {
                            if (!pIt.second[i].isLocal)
                            {
                                continue;
                            }

                            int meshEdgeId2 = pIt.second[i].id;
                            if(edgeDirMap.count(meshEdgeId2)==0)
                            {
                                if(uniqueEdgeMap.count(meshEdgeId2)!=0)
                                {
                                    // set unique map to same location
                                    uniqueEdgeMap[meshEdgeId] =
                                        uniqueEdgeMap[meshEdgeId2];
                                    SetUpNewEdge = false;
                                }
                            }
                            else
                            {
                                edgeDirMap.insert(meshEdgeId);
                                SetUpNewEdge = false;
                            }
                        }

                        if(SetUpNewEdge)
                        {
                            uniqueEdgeMap[meshEdgeId]=matrixlocation;
                            globaloffset[matrixlocation]+=ntotalentries;
                            modeoffset[matrixlocation]=dof*dof;
                            ntotalentries+=dof*dof;
                            nblks [matrixlocation++]  = dof;
                        }
                    }
                }
            }

            for(cnt=n=0; n < n_exp; ++n)
            {
                locExpansion = expList->GetExp(n);

                for (j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);
                    dof    = EdgeSize[meshEdgeId];
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        if(uniqueEdgeMap.count(meshEdgeId)==0 && dof > 0)

                        {
                            uniqueEdgeMap[meshEdgeId]=matrixlocation;

                            globaloffset[matrixlocation]+=ntotalentries;
                            modeoffset[matrixlocation]=dof*dof;
                            ntotalentries+=dof*dof;
                            nblks[matrixlocation++]   = dof;
                        }
                        nlocalNonDirEdges+=dof*dof;
                    }
                }
            }

            // Loop over all the elements in the domain and compute max face
            // DOF. Reduce across all processes to get universal maximum.
            // - Periodic faces
            for (auto &pIt : periodicFaces)
            {
                meshFaceId = pIt.first;

                if(faceDirMap.count(meshFaceId)==0)
                {
                    dof = FaceSize[meshFaceId];

                    if(uniqueFaceMap.count(meshFaceId) == 0 && dof > 0)
                    {
                        bool SetUpNewFace = true;

                        if(pIt.second[0].isLocal)
                        {
                            int meshFaceId2 = pIt.second[0].id;

                            if(faceDirMap.count(meshFaceId2)==0)
                            {
                                if(uniqueFaceMap.count(meshFaceId2)!=0)
                                {
                                    // set unique map to same location
                                    uniqueFaceMap[meshFaceId] =
                                        uniqueFaceMap[meshFaceId2];
                                    SetUpNewFace = false;
                                }
                            }
                            else // set face to be a Dirichlet face
                            {
                                faceDirMap.insert(meshFaceId);
                                SetUpNewFace = false;
                            }
                        }

                        if(SetUpNewFace)
                        {
                            uniqueFaceMap[meshFaceId]=matrixlocation;

                            modeoffset[matrixlocation]=dof*dof;
                            globaloffset[matrixlocation]+=ntotalentries;
                            ntotalentries+=dof*dof;
                            nblks[matrixlocation++] = dof;
                        }
                    }
                }
            }

            for(cnt=n=0; n < n_exp; ++n)
            {
                locExpansion = expList->GetExp(n);

                for (j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->
                        GetGeom3D()->GetFid(j);

                    dof        = FaceSize[meshFaceId];
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);

                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        if(uniqueFaceMap.count(meshFaceId)==0 && dof > 0)
                        {
                            uniqueFaceMap[meshFaceId]=matrixlocation;
                            modeoffset[matrixlocation]=dof*dof;
                            globaloffset[matrixlocation]+=ntotalentries;
                            ntotalentries+=dof*dof;
                            nblks[matrixlocation++] = dof;
                        }
                        nlocalNonDirFaces+=dof*dof;
                    }
                }
            }

            m_comm->AllReduce(maxEdgeDof, ReduceMax);
            m_comm->AllReduce(maxFaceDof, ReduceMax);

            //Allocate arrays for block to universal map (number of expansions * p^2)
            Array<OneD, long> BlockToUniversalMap(ntotalentries,-1);
            Array<OneD, int> localToGlobalMatrixMap(nlocalNonDirEdges +
                                                    nlocalNonDirFaces,-1);

            //Allocate arrays to store matrices (number of expansions * p^2)
            Array<OneD, NekDouble> BlockArray(nlocalNonDirEdges +
                                              nlocalNonDirFaces,0.0);

            int matrixoffset=0;
            int vGlobal;
            int uniEdgeOffset = 0;

            // Need to obtain a fixed offset for the universal number
            // of the faces which come after the edge numbering
            for(n=0; n < n_exp; ++n)
            {
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetEid(j);

                    uniEdgeOffset = max(meshEdgeId, uniEdgeOffset);
                }
            }
            uniEdgeOffset++;

            m_comm->AllReduce(uniEdgeOffset,ReduceMax);
            uniEdgeOffset = uniEdgeOffset*maxEdgeDof*maxEdgeDof;

            for(n=0; n < n_exp; ++n)
            {
                locExpansion = expList->GetExp(n);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetEid(j);

                    nedgemodes = EdgeSize[meshEdgeId];

                    if(edgeDirMap.count(meshEdgeId)==0 && nedgemodes > 0)
                    {
                        // Determine the Global edge offset
                        int edgeOffset = globaloffset[uniqueEdgeMap[meshEdgeId]];

                        // Determine a universal map offset
                        int uniOffset = meshEdgeId;
                        auto pIt = periodicEdges.find(meshEdgeId);
                        if (pIt != periodicEdges.end())
                        {
                            for (int l = 0; l < pIt->second.size(); ++l)
                            {
                                uniOffset = min(uniOffset, pIt->second[l].id);
                            }
                        }
                        uniOffset = uniOffset*maxEdgeDof*maxEdgeDof;

                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=edgeOffset+k;
                            localToGlobalMatrixMap[matrixoffset+k]=vGlobal;
                            BlockToUniversalMap[vGlobal] = uniOffset + k + 1;
                        }
                        matrixoffset+=nedgemodes*nedgemodes;
                    }
                }

                Array<OneD, unsigned int>           faceInteriorMap;
                Array<OneD, int>                    faceInteriorSign;
                //loop over the faces of the expansion
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    //get mesh face id
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetFid(j);

                    nfacemodes = FaceSize[meshFaceId];

                    //Check if face has dirichlet values
                    if(faceDirMap.count(meshFaceId)==0 && nfacemodes > 0)
                    {
                        // Determine the Global edge offset
                        int faceOffset = globaloffset[uniqueFaceMap[meshFaceId]];
                        // Determine a universal map offset
                        int uniOffset = meshFaceId;
                        // use minimum face edge when periodic
                        auto pIt = periodicFaces.find(meshFaceId);
                        if (pIt != periodicFaces.end())
                        {
                            uniOffset = min(uniOffset, pIt->second[0].id);
                        }
                        uniOffset = uniOffset * maxFaceDof * maxFaceDof;

                        for(k=0; k<nfacemodes*nfacemodes; ++k)
                        {
                            vGlobal=faceOffset+k;

                            localToGlobalMatrixMap[matrixoffset+k]
                                = vGlobal;

                            BlockToUniversalMap[vGlobal] = uniOffset +
                                uniEdgeOffset + k + 1;
                        }
                        matrixoffset+=nfacemodes*nfacemodes;
                    }
                }
            }

            matrixoffset=0;

            map<int,int>::iterator it;
            Array<OneD, unsigned int> n_blks(nblks.size()+1);
            n_blks[0] = nNonDirVerts;
            for(i =1, it = nblks.begin(); it != nblks.end(); ++it)
            {
                n_blks[i++] = it->second;
            }

            m_BlkMat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                locExpansion = expList->GetExp(n);
                nCoeffs=locExpansion->NumBndryCoeffs();

                //Get correct transformation matrix for element type
                DNekMat &R = (*m_RBlk->GetBlock(n,n));

                pRSRT = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nCoeffs, nCoeffs, zero, storage);
                RSRT = (*pRSRT);

                nVerts=locExpansion->GetGeom()->GetNumVerts();
                nEdges=locExpansion->GetGeom()->GetNumEdges();
                nFaces=locExpansion->GetGeom()->GetNumFaces();

                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

                //Extract boundary block (elemental S1)
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                offset = bnd_mat->GetRows();

                DNekScalMat &S=(*bnd_mat);

                DNekMat Sloc(nCoeffs,nCoeffs);

                // For variable p we need to just use the active modes
                NekDouble val;

                for(int i = 0; i < nCoeffs; ++i)
                {
                    for(int j = 0; j < nCoeffs; ++j)
                    {
                        val = S(i,j)*m_variablePmask[cnt+i]*m_variablePmask[cnt+j];
                        Sloc.SetValue(i,j,val);
                    }
                }

                //Calculate R*S*trans(R)
                RSRT = R*Sloc*Transpose(R);

                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    vMap1=locExpansion->GetVertexMap(v);

                    //Get vertex map
                    globalrow = asmMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=locExpansion->GetVertexMap(m);

                            //global matrix location (without offset due to
                            //dirichlet values)
                            globalcol = asmMap->
                                GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol == globalrow)
                            {
                                //modal connectivity between elements
                                sign1 = asmMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = asmMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);

                                vertArray[globalrow]
                                    += sign1*sign2*RSRT(vMap1,vMap2);


                                meshVertId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetVid(v);

                                auto pIt = periodicVerts.find(meshVertId);
                                if (pIt != periodicVerts.end())
                                {
                                    for (k = 0; k < pIt->second.size(); ++k)
                                    {
                                        meshVertId = min(meshVertId, pIt->second[k].id);
                                    }
                                }

                                VertBlockToUniversalMap[globalrow]
                                    = meshVertId + 1;
                            }
                        }
                    }
                }

                //loop over edges of the element and return the edge map
                for (eid=0; eid<nEdges; ++eid)
                {

                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetEid(eid);


                    nedgemodes    = EdgeSize[meshEdgeId];
                    if(nedgemodes)
                    {
                        nedgemodesloc = locExpansion->GetEdgeNcoeffs(eid)-2;
                        DNekMatSharedPtr m_locMat =
                            MemoryManager<DNekMat>::AllocateSharedPtr
                            (nedgemodes,nedgemodes,zero,storage);

                        Array<OneD, unsigned int> edgemodearray = locExpansion->GetEdgeInverseBoundaryMap(eid);

                        if(edgeDirMap.count(meshEdgeId)==0)
                        {
                            for (v=0; v<nedgemodesloc; ++v)
                            {
                                eMap1=edgemodearray[v];
                                sign1 = asmMap->
                                    GetLocalToGlobalBndSign(cnt + eMap1);

                                if(sign1 == 0)
                                {
                                    continue;
                                }

                                for (m=0; m<nedgemodesloc; ++m)
                                {
                                    eMap2=edgemodearray[m];

                                    //modal connectivity between elements
                                    sign2 = asmMap->
                                        GetLocalToGlobalBndSign(cnt + eMap2);

                                    NekDouble globalEdgeValue = sign1*sign2*RSRT(eMap1,eMap2);

                                    if(sign2 != 0)
                                    {
                                        //if(eMap1 == eMap2)
                                        BlockArray[matrixoffset+v*nedgemodes+m]=globalEdgeValue;
                                    }
                                }
                            }
                            matrixoffset+=nedgemodes*nedgemodes;
                        }
                    }
                }

                //loop over faces of the element and return the face map
                for (fid=0; fid<nFaces; ++fid)
                {
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()
                        ->GetGeom3D()->GetFid(fid);

                    nfacemodes   = FaceSize[meshFaceId];
                    if(nfacemodes > 0)
                    {
                        DNekMatSharedPtr m_locMat =
                            MemoryManager<DNekMat>::AllocateSharedPtr
                            (nfacemodes,nfacemodes,zero,storage);

                        if(faceDirMap.count(meshFaceId) == 0)
                        {
                            Array<OneD, unsigned int> facemodearray;
                            StdRegions::Orientation faceOrient =
                                locExpansion->GetForient(fid);

                            auto pIt = periodicFaces.find(meshFaceId);
                            if (pIt != periodicFaces.end())
                            {
                                if(meshFaceId == min(meshFaceId, pIt->second[0].id))
                                {
                                    faceOrient = DeterminePeriodicFaceOrient
                                        (faceOrient,pIt->second[0].orient);
                                }
                            }

                            facemodearray = locExpansion->GetFaceInverseBoundaryMap
                                (fid,faceOrient,FaceModes[meshFaceId].first,
                                 FaceModes[meshFaceId].second);

                            for (v=0; v<nfacemodes; ++v)
                            {
                                fMap1=facemodearray[v];

                                sign1 = asmMap->
                                    GetLocalToGlobalBndSign(cnt + fMap1);

                                ASSERTL1(sign1 != 0,"Something is wrong since we "
                                         "shoudl only be extracting modes for "
                                         "lowest order expansion");

                                for (m=0; m<nfacemodes; ++m)
                                {
                                    fMap2=facemodearray[m];

                                    //modal connectivity between elements
                                    sign2 = asmMap->
                                        GetLocalToGlobalBndSign(cnt + fMap2);

                                    ASSERTL1(sign2 != 0,"Something is wrong since "
                                             "we shoudl only be extracting modes for "
                                             "lowest order expansion");

                                    // Get the face-face value from the
                                    // low energy matrix (S2)
                                    NekDouble globalFaceValue = sign1*sign2*
                                        RSRT(fMap1,fMap2);

                                    //local face value to global face value
                                    //if(fMap1 == fMap2)
                                    BlockArray[matrixoffset+v*nfacemodes+m]=
                                        globalFaceValue;
                                }
                            }
                            matrixoffset+=nfacemodes*nfacemodes;
                        }
                    }
                }

                //offset for the expansion
                cnt+=nCoeffs;
            }

            if(nNonDirVerts!=0)
            {
                //Exchange vertex data over different processes
                Gs::gs_data *tmp = Gs::Init(VertBlockToUniversalMap, m_comm, verbose);
                Gs::Gather(vertArray, Gs::gs_add, tmp);

            }

            Array<OneD, NekDouble> GlobalBlock(ntotalentries,0.0);
            if(ntotalentries)
            {
                //Assemble edge matrices of each process
                Vmath::Assmb(BlockArray.size(),
                             BlockArray,
                             localToGlobalMatrixMap,
                             GlobalBlock);
            }

            //Exchange edge & face data over different processes
            Gs::gs_data *tmp1 = Gs::Init(BlockToUniversalMap, m_comm, verbose);
            Gs::Gather(GlobalBlock, Gs::gs_add, tmp1);

            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }

            //Set the first block to be the diagonal of the vertex space
            m_BlkMat->SetBlock(0,0, VertBlk);

            //Build the edge and face matrices from the vector
            DNekMatSharedPtr gmat;

            offset=0;
            // -1 since we ignore vert block
            for(int loc=0; loc<n_blks.size()-1; ++loc)
            {
                nmodes = n_blks[1+loc];
                gmat = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                for (v=0; v<nmodes; ++v)
                {
                    for (m=0; m<nmodes; ++m)
                    {
                        NekDouble Value = GlobalBlock[offset+v*nmodes+m];
                        gmat->SetValue(v,m,Value);

                    }
                }
                m_BlkMat->SetBlock(1+loc,1+loc, gmat);
                offset+=modeoffset[loc];
            }

            // invert blocks.
            int totblks=m_BlkMat->GetNumberOfBlockRows();
            for (i=1; i< totblks; ++i)
            {
                unsigned int nmodes=m_BlkMat->GetNumberOfRowsInBlockRow(i);
                if(nmodes)
                {
                    DNekMatSharedPtr tmp_mat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                    tmp_mat=m_BlkMat->GetBlock(i,i);
                    tmp_mat->Invert();

                    m_BlkMat->SetBlock(i,i,tmp_mat);
                }
            }
        }


        /**
         * Apply the low energy preconditioner during the conjugate gradient
         * routine
         */
        void PreconditionerLowEnergy::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nDir    = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();
            int nGlobal = m_locToGloMap.lock()->GetNumGlobalBndCoeffs();
            int nNonDir = nGlobal-nDir;
            DNekBlkMat &M = (*m_BlkMat);

            NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
            NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);

            z = M * r;
	}


        /**
         * Set a block transformation matrices for each element type. These are
         * needed in routines that transform the schur complement matrix to and
         * from the low energy basis.
         */
        void PreconditionerLowEnergy::SetupBlockTransformationMatrix(void)
        {
            std::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            StdRegions::StdExpansionSharedPtr locExp;
            StdRegions::StdExpansionSharedPtr locExpSav;
            map<int,int> EdgeSize;

            int n;

            std::map<ShapeType, DNekScalMatSharedPtr>         maxRmat;
            map<ShapeType, LocalRegions::ExpansionSharedPtr > maxElmt;
            map<ShapeType, Array<OneD, unsigned int> >        vertMapMaxR;
            map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > edgeMapMaxR;


            //Sets up reference element and builds transformation matrix for
            // maximum polynomial order meshes
            SetUpReferenceElements(maxRmat,maxElmt,vertMapMaxR,edgeMapMaxR);

            const Array<OneD,const unsigned int>& nbdry_size
                = m_locToGloMap.lock()->GetNumLocalBndCoeffsPerPatch();

            int n_exp=expList->GetNumElmts();

            MatrixStorage blkmatStorage = eDIAGONAL;

            //Variants of R matrices required for low energy preconditioning
            m_RBlk      = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
            m_InvRBlk      = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);

            DNekMatSharedPtr rmat, invrmat;

            int offset = 0;

            // Set up transformation matrices whilst checking to see if
            // consecutive matrices are the same and if so reuse the
            // matrices and store how many consecutive offsets there
            // are
            for(n=0; n < n_exp; ++n)
            {
                locExp = expList->GetExp(n);
                ShapeType eltype = locExp->DetShapeType();

                int nbndcoeffs = locExp->NumBndryCoeffs();

                if(m_sameBlock.size() == 0)
                {
                    rmat = ExtractLocMat(locExp,maxRmat[eltype],
                                         maxElmt[eltype],
                                         vertMapMaxR[eltype],
                                         edgeMapMaxR[eltype]);
                    //Block R matrix
                    m_RBlk->SetBlock(n, n, rmat);

                    invrmat = MemoryManager<DNekMat>::AllocateSharedPtr(*rmat);
                    invrmat->Invert();

                    //Block inverse R matrix
                    m_InvRBlk->SetBlock(n, n, invrmat);

                    m_sameBlock.push_back(pair<int,int>(1,nbndcoeffs));
                    locExpSav = locExp;
                }
                else
                {
                    bool reuse = true;

                    // check to see if same as previous matrix and
                    // reuse if we can
                    for(int i = 0; i < 3; ++i)
                    {
                        if(locExpSav->GetBasis(i) != locExp->GetBasis(i))
                        {
                            reuse = false;
                            break;
                        }
                    }

                    if(reuse)
                    {
                        m_RBlk->SetBlock(n, n, rmat);
                        m_InvRBlk->SetBlock(n, n, invrmat);

                        m_sameBlock[offset] =
                            (pair<int,int>(m_sameBlock[offset].first+1,nbndcoeffs));
                    }
                    else
                    {
                        rmat = ExtractLocMat(locExp,maxRmat[eltype],
                                             maxElmt[eltype],
                                             vertMapMaxR[eltype],
                                             edgeMapMaxR[eltype]);

                        //Block R matrix
                        m_RBlk->SetBlock(n, n, rmat);

                        invrmat = MemoryManager<DNekMat>::AllocateSharedPtr(*rmat);
                        invrmat->Invert();
                        //Block inverse R matrix
                        m_InvRBlk->SetBlock(n, n, invrmat);

                        m_sameBlock.push_back(pair<int,int>(1,nbndcoeffs));
                        offset++;
                        locExpSav = locExp;
                    }
                }
            }
        }

        /**
         * \brief Transform the basis  vector to low energy space
         *
         * As the conjugate gradient system is solved for the low energy basis,
         * the solution vector \f$\mathbf{x}\f$ must be transformed to the low
         * energy basis i.e. \f$\overline{\mathbf{x}}=\mathbf{R}\mathbf{x}\f$.
         */
        void PreconditionerLowEnergy::v_DoTransformBasisToLowEnergy
        (Array<OneD, NekDouble>& pInOut)
        {
            int nLocBndDofs  = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            //Block transformation matrix
            DNekBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs,pInOut.get());

            //Apply mask in case of variable P
	    Vmath::Vmul(nLocBndDofs,pLocalIn, 1, m_variablePmask,1,
			pLocalIn,1);

	    //Multiply by the block transformation matrix
	    int cnt = 0;
	    int cnt1 = 0;
	    for(int i = 0; i < m_sameBlock.size(); ++i)
	    {
		int nexp    = m_sameBlock[i].first;
		int nbndcoeffs = m_sameBlock[i].second;
		Blas::Dgemm('N','N', nbndcoeffs, nexp, nbndcoeffs,
			    1.0, &(R.GetBlock(cnt1,cnt1)->GetPtr()[0]),
			    nbndcoeffs,pLocalIn.get() + cnt,  nbndcoeffs,
			    0.0, pInOut.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }
        }

        /**
         * \brief transform the solution coeffiicents from low energy
         * back to the original coefficient space.
         *
         * After the conjugate gradient routine the output vector is in the low
         * energy basis and must be trasnformed back to the original basis in
         * order to get the correct solution out. the solution vector
         * i.e. \f$\mathbf{x}=\mathbf{R^{T}}\mathbf{\overline{x}}\f$.
         */
        void PreconditionerLowEnergy::v_DoTransformCoeffsFromLowEnergy(
                              Array<OneD, NekDouble>& pInOut)
        {
            int nLocBndDofs   = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            ASSERTL1(pInOut.size() >= nLocBndDofs,
                     "Output array is not greater than the nLocBndDofs");

            //Block transposed transformation matrix
            DNekBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, pInOut.get());

            //Multiply by the transpose of block transformation matrix
	    int cnt = 0;
	    int cnt1 = 0;
	    for(int i = 0; i < m_sameBlock.size(); ++i)
            {
	        int nexp    = m_sameBlock[i].first;
		int nbndcoeffs = m_sameBlock[i].second;
		Blas::Dgemm('T','N', nbndcoeffs, nexp, nbndcoeffs,
			    1.0, &(R.GetBlock(cnt1,cnt1)->GetPtr()[0]),
                            nbndcoeffs,pLocalIn.get() + cnt,  nbndcoeffs, 
			    0.0, pInOut.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;

	    }

            Vmath::Vmul(nLocBndDofs,pInOut, 1, m_variablePmask,1,
                        pInOut,1);
        }

        /**
         * \brief Multiply by the block inverse transformation matrix
         * This transforms the bassi from Low Energy to original basis
         *
         * Note; not typically required
         */ 
        
        void PreconditionerLowEnergy::v_DoTransformBasisFromLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput)
        {
            int nLocBndDofs    = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            ASSERTL1(pInput.size() >= nLocBndDofs,
                     "Input array is smaller than nLocBndDofs");
            ASSERTL1(pOutput.size() >= nLocBndDofs,
                     "Output array is smaller than nLocBndDofs");

            //Block inverse transformation matrix
            DNekBlkMat &invR = *m_InvRBlk;

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, pInput.get());

            //Multiply by the inverse transformation matrix
	    int cnt = 0;
	    int cnt1 = 0;
	    for(int i = 0; i < m_sameBlock.size(); ++i)
            {
	        int nexp    = m_sameBlock[i].first;
		int nbndcoeffs = m_sameBlock[i].second;
		Blas::Dgemm('N','N', nbndcoeffs, nexp, nbndcoeffs,
			    1.0, &(invR.GetBlock(cnt1,cnt1)->GetPtr()[0]),
                            nbndcoeffs,pLocalIn.get() + cnt,  nbndcoeffs, 
			    0.0, pOutput.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }
	}

        /**
         * \brief Multiply by the block tranposed inverse
         * transformation matrix (R^T)^{-1} which is equivlaent to
         * transforming coefficients to LowEnergy space
         *
         * In JCP 2001 paper on low energy this is seen as (C^T)^{-1}
         */ 
        void PreconditionerLowEnergy::v_DoTransformCoeffsToLowEnergy(
                        const Array<OneD, NekDouble>& pInput,
                        Array<OneD, NekDouble>& pOutput)
        {
            int nLocBndDofs     = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            ASSERTL1(pInput.size() >= nLocBndDofs,
                     "Input array is less than nLocBndDofs");
            ASSERTL1(pOutput.size() >= nLocBndDofs,
                     "Output array is less than nLocBndDofs");

            //Block inverse transformation matrix
            DNekBlkMat &invR = *m_InvRBlk;

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, pInput.get());

            //Multiply by the transpose of block transformation matrix
	    int cnt = 0;
	    int cnt1 = 0;
	    for(int i = 0; i < m_sameBlock.size(); ++i)
            {
	        int nexp    = m_sameBlock[i].first;
		int nbndcoeffs = m_sameBlock[i].second;
		Blas::Dgemm('T','N', nbndcoeffs, nexp, nbndcoeffs,
			    1.0, &(invR.GetBlock(cnt1,cnt1)->GetPtr()[0]),
                            nbndcoeffs,pLocalIn.get() + cnt,  nbndcoeffs, 
			    0.0, pOutput.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }
	}

        /**
         * \brief Set up the transformed block  matrix system
         *
         * Sets up a block elemental matrix in which each of the block matrix is
         * the low energy equivalent
         * i.e. \f$\mathbf{S}_{2}=\mathbf{R}\mathbf{S}_{1}\mathbf{R}^{T}\f$
         */
        DNekScalMatSharedPtr PreconditionerLowEnergy::
        v_TransformedSchurCompl( int n,
                                 int bndoffset,
                                 const std::shared_ptr<DNekScalMat > &loc_mat)
	{
            std::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();

            LocalRegions::ExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(n);
            unsigned int nbnd=locExpansion->NumBndryCoeffs();

            MatrixStorage storage = eFULL;
            DNekMatSharedPtr pS2 = MemoryManager<DNekMat>::
                AllocateSharedPtr(nbnd,nbnd,0.0,storage);

            //transformation matrices
            DNekMat &R = (*m_RBlk->GetBlock(n,n));

            // Original Schur Complement
            DNekScalMat &S1 = (*loc_mat);

            DNekMat Sloc(nbnd,nbnd);
            
            // For variable p we need to just use the active modes 
            NekDouble val;

            for(int i = 0; i < nbnd; ++i)
            {
                for(int j = 0; j < nbnd; ++j)
                {
                    val = S1(i,j)*m_variablePmask[bndoffset+i]*
                        m_variablePmask[bndoffset+j];
                    Sloc.SetValue(i,j,val);
                }
            }

            //create low energy matrix
            DNekMat &S2 = (*pS2);

            S2= R*Sloc*Transpose(R);

            DNekScalMatSharedPtr return_val;
            return_val = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, pS2);

	    return return_val;
	}

        /**
         * Create the inverse multiplicity map.
         */
        void PreconditionerLowEnergy::CreateVariablePMask(void)
        {
            unsigned int nLocBnd  = m_locToGloMap.lock()->GetNumLocalBndCoeffs();
            unsigned int i;
            auto asmMap = m_locToGloMap.lock();
            
            const Array< OneD, const NekDouble > &sign 
                = asmMap->GetLocalToGlobalBndSign();

            m_signChange=asmMap->GetSignChange();

            // Construct a map of 1/multiplicity
            m_variablePmask = Array<OneD, NekDouble>(nLocBnd);
            for (i = 0; i < nLocBnd; ++i)
            {
                if(m_signChange)
                {
                    m_variablePmask[i] = fabs(sign[i]);
                }
                else
                {
                    m_variablePmask[i] = 1.0; 
                }
            }
        }

        /**
         *\brief Sets up the reference prismatic element needed to construct
         *a low energy basis
         */
        SpatialDomains::PrismGeomSharedPtr PreconditionerLowEnergy::CreateRefPrismGeom()
        {
            //////////////////////////
            // Set up Prism element //
            //////////////////////////

	    const int three=3;
            const int nVerts = 6;
            const double point[][3] = {
                {-1,-1,0}, {1,-1,0}, {1,1,0},
                {-1,1,0}, {0,-1,sqrt(double(3))}, {0,1,sqrt(double(3))},
            };

            //std::shared_ptr<SpatialDomains::PointGeom> verts[6];
            SpatialDomains::PointGeomSharedPtr verts[6];
            for(int i=0; i < nVerts; ++i)
            {
                verts[i] =  MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr
                    ( three, i, point[i][0], point[i][1], point[i][2] );
            }
            const int nEdges = 9;
            const int vertexConnectivity[][2] = {
                {0,1}, {1,2}, {3,2}, {0,3}, {0,4},
                {1,4}, {2,5}, {3,5}, {4,5}
            };

            // Populate the list of edges
            SpatialDomains::SegGeomSharedPtr edges[nEdges];
            for(int i=0; i < nEdges; ++i){
                SpatialDomains::PointGeomSharedPtr vertsArray[2];
                for(int j=0; j<2; ++j)
                {
                    vertsArray[j] = verts[vertexConnectivity[i][j]];
                }
                edges[i] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(i, three, vertsArray);
            }

            ////////////////////////
            // Set up Prism faces //
            ////////////////////////

            const int nFaces = 5;
            //quad-edge connectivity base-face0, vertical-quadface2, vertical-quadface4
            const int quadEdgeConnectivity[][4] = { {0,1,2,3}, {1,6,8,5}, {3,7,8,4} };
            // QuadId ordered as 0, 1, 2, otherwise return false
            const int                  quadId[] = { 0,-1,1,-1,2 };

            //triangle-edge connectivity side-triface-1, side triface-3
            const int  triEdgeConnectivity[][3] = { {0,5,4}, {2,6,7} };
            // TriId ordered as 0, 1, otherwise return false
            const int                   triId[] = { -1,0,-1,1,-1 };

            // Populate the list of faces
            SpatialDomains::Geometry2DSharedPtr faces[nFaces];
            for(int f = 0; f < nFaces; ++f){
                if(f == 1 || f == 3) {
                    int i = triId[f];
                    SpatialDomains::SegGeomSharedPtr edgeArray[3];
                    for(int j = 0; j < 3; ++j){
                        edgeArray[j] = edges[triEdgeConnectivity[i][j]];
                    }
                    faces[f] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(f, edgeArray);
                }
                else {
                    int i = quadId[f];
                    SpatialDomains::SegGeomSharedPtr edgeArray[4];
                    for(int j=0; j < 4; ++j){
                        edgeArray[j] = edges[quadEdgeConnectivity[i][j]];
                    }
                    faces[f] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(f, edgeArray);
                }
            }

            SpatialDomains::PrismGeomSharedPtr geom = MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(0, faces);

            return geom;
        }

        /**
         *\brief Sets up the reference prismatic element needed to construct
         *a low energy basis mapping arrays
         */
        SpatialDomains::PyrGeomSharedPtr PreconditionerLowEnergy::CreateRefPyrGeom()
        {
            //////////////////////////
            // Set up Pyramid element //
            //////////////////////////

            const int nVerts = 5;
            const double point[][3] = {
                {-1,-1,0}, {1,-1,0}, {1,1,0},
                {-1,1,0}, {0,0,sqrt(double(2))}
            };

            //boost::shared_ptr<SpatialDomains::PointGeom> verts[6];
	    const int three=3;
            SpatialDomains::PointGeomSharedPtr verts[5];
            for(int i=0; i < nVerts; ++i)
            {
                verts[i] =  MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr
                    ( three, i, point[i][0], point[i][1], point[i][2] );
            }
            const int nEdges = 8;
            const int vertexConnectivity[][2] = {
                {0,1}, {1,2}, {2,3}, {3,0},
                {0,4}, {1,4}, {2,4}, {3,4}
            };

            // Populate the list of edges
            SpatialDomains::SegGeomSharedPtr edges[nEdges];
            for(int i=0; i < nEdges; ++i)
            {
                SpatialDomains::PointGeomSharedPtr vertsArray[2];
                for(int j=0; j<2; ++j)
                {
                    vertsArray[j] = verts[vertexConnectivity[i][j]];
                }
                edges[i] = MemoryManager<SpatialDomains::SegGeom>::AllocateSharedPtr(i, three, vertsArray);
            }

            ////////////////////////
            // Set up Pyramid faces //
            ////////////////////////

            const int nFaces = 5;
            //quad-edge connectivity base-face0,
            const int quadEdgeConnectivity[][4] = {{0,1,2,3}};

            //triangle-edge connectivity side-triface-1, side triface-2
            const int  triEdgeConnectivity[][3] = { {0,5,4}, {1,6,5}, {2,7,6}, {3,4,7}};

            // Populate the list of faces
            SpatialDomains::Geometry2DSharedPtr faces[nFaces];
            for(int f = 0; f < nFaces; ++f)
            {
                if(f == 0)
                {
                    SpatialDomains::SegGeomSharedPtr edgeArray[4];
                    for(int j=0; j < 4; ++j)
                    {
                        edgeArray[j] = edges[quadEdgeConnectivity[f][j]];
                    }

                    faces[f] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(f,edgeArray);
                }
                else {
                    int i = f-1;
                    SpatialDomains::SegGeomSharedPtr edgeArray[3];
                    for(int j = 0; j < 3; ++j)
                    {
                        edgeArray[j] = edges[triEdgeConnectivity[i][j]];
                    }
                    faces[f] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(f, edgeArray);
                }
            }

            SpatialDomains::PyrGeomSharedPtr geom =
                MemoryManager<SpatialDomains::PyrGeom>::AllocateSharedPtr(0,faces);

            return geom;
        }

        /**
         *\brief Sets up the reference tretrahedral element needed to construct
         *a low energy basis
         */
        SpatialDomains::TetGeomSharedPtr PreconditionerLowEnergy::CreateRefTetGeom()
        {
            /////////////////////////////////
            // Set up Tetrahedron vertices //
            /////////////////////////////////

	    int i,j;
	    const int three=3;
            const int nVerts = 4;
            const double point[][3] = {
                {-1,-1/sqrt(double(3)),-1/sqrt(double(6))},
                {1,-1/sqrt(double(3)),-1/sqrt(double(6))},
                {0,2/sqrt(double(3)),-1/sqrt(double(6))},
                {0,0,3/sqrt(double(6))}};

            std::shared_ptr<SpatialDomains::PointGeom> verts[4];
	    for(i=0; i < nVerts; ++i)
	    {
	        verts[i] =
                    MemoryManager<SpatialDomains::PointGeom>::
                    AllocateSharedPtr
                    ( three, i, point[i][0], point[i][1], point[i][2] );
	    }

            //////////////////////////////
            // Set up Tetrahedron Edges //
            //////////////////////////////

            // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
            const int nEdges = 6;
            const int vertexConnectivity[][2] = {
                {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
            };

            // Populate the list of edges
            SpatialDomains::SegGeomSharedPtr edges[nEdges];
            for(i=0; i < nEdges; ++i)
            {
                std::shared_ptr<SpatialDomains::PointGeom>
                    vertsArray[2];
                for(j=0; j<2; ++j)
                {
                    vertsArray[j] = verts[vertexConnectivity[i][j]];
                }

                edges[i] = MemoryManager<SpatialDomains::SegGeom>
                    ::AllocateSharedPtr(i, three, vertsArray);
            }

            //////////////////////////////
            // Set up Tetrahedron faces //
            //////////////////////////////

            const int nFaces = 4;
            const int edgeConnectivity[][3] = {
                {0,1,2}, {0,4,3}, {1,5,4}, {2,5,3}
            };

            // Populate the list of faces
            SpatialDomains::TriGeomSharedPtr faces[nFaces];
            for(i=0; i < nFaces; ++i)
            {
                SpatialDomains::SegGeomSharedPtr edgeArray[3];
                for(j=0; j < 3; ++j)
                {
                    edgeArray[j] = edges[edgeConnectivity[i][j]];
                }


                faces[i] = MemoryManager<SpatialDomains::TriGeom>
                    ::AllocateSharedPtr(i, edgeArray);
            }

            SpatialDomains::TetGeomSharedPtr geom =
                MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr
                (0, faces);

            return geom;
        }

        /**
         *\brief Sets up the reference hexahedral element needed to construct
         *a low energy basis
         */
        SpatialDomains::HexGeomSharedPtr PreconditionerLowEnergy::CreateRefHexGeom()
        {
            ////////////////////////////////
            // Set up Hexahedron vertices //
            ////////////////////////////////

	    const int three=3;

            const int nVerts = 8;
            const double point[][3] = {
                {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
                {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}
            };

            // Populate the list of verts
            SpatialDomains::PointGeomSharedPtr verts[8];
            for( int i = 0; i < nVerts; ++i ) {
                verts[i] = MemoryManager<SpatialDomains::PointGeom>
                    ::AllocateSharedPtr(three,  i,   point[i][0],
                                        point[i][1], point[i][2]);
            }

            /////////////////////////////
            // Set up Hexahedron Edges //
            /////////////////////////////

            // SegGeom (int id, const int coordim), EdgeComponent(id, coordim)
            const int nEdges = 12;
            const int vertexConnectivity[][2] = {
                {0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,5},
                {2,6}, {3,7}, {4,5}, {5,6}, {6,7}, {4,7}
            };

            // Populate the list of edges
            SpatialDomains::SegGeomSharedPtr edges[nEdges];
            for( int i = 0; i < nEdges; ++i ) {
                SpatialDomains::PointGeomSharedPtr vertsArray[2];
                for( int j = 0; j < 2; ++j ) {
                    vertsArray[j] = verts[vertexConnectivity[i][j]];
                }
                edges[i] = MemoryManager<SpatialDomains::SegGeom>::
                    AllocateSharedPtr( i, three, vertsArray);
            }

            /////////////////////////////
            // Set up Hexahedron faces //
            /////////////////////////////

            const int nFaces = 6;
            const int edgeConnectivity[][4] = {
                {0,1,2,3}, {0,5,8,4}, {1,6,9,5},
                {2,7,10,6}, {3,7,11,4}, {8,9,10,11}
            };

            // Populate the list of faces
            SpatialDomains::QuadGeomSharedPtr faces[nFaces];
            for( int i = 0; i < nFaces; ++i ) {
                SpatialDomains::SegGeomSharedPtr edgeArray[4];
                for( int j = 0; j < 4; ++j ) {
                    edgeArray[j]    = edges[edgeConnectivity[i][j]];
                }
                faces[i] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(i, edgeArray);
            }

            SpatialDomains::HexGeomSharedPtr geom =
                MemoryManager<SpatialDomains::HexGeom>::AllocateSharedPtr
                (0, faces);

            return geom;
        }


        /**
	 * \brief Loop expansion and determine different variants of the
	 * transformation matrix
	 *
         * Sets up multiple reference elements based on the element expansion.
	 */
        void PreconditionerLowEnergy::SetUpReferenceElements(
                 std::map<ShapeType, DNekScalMatSharedPtr> &maxRmat,
                 map<ShapeType, LocalRegions::ExpansionSharedPtr > &maxElmt,
                 map<ShapeType, Array<OneD, unsigned int> >        &vertMapMaxR,
                 map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &edgeMapMaxR)
        {

            std::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            GlobalLinSysKey linSysKey=(m_linsys.lock())->GetKey();

            LocalRegions::ExpansionSharedPtr locExp;

            // face maps for pyramid and hybrid meshes - not need to return.
            map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > faceMapMaxR;

            /* Determine the maximum expansion order for all elements */
            int nummodesmax = 0;
            Array<OneD, int> Shapes(LibUtilities::SIZE_ShapeType,0);

            for(int n = 0; n < expList->GetNumElmts(); ++n)
            {
                locExp = expList->GetExp(n);

                nummodesmax = max(nummodesmax, locExp->GetBasisNumModes(0));
                nummodesmax = max(nummodesmax, locExp->GetBasisNumModes(1));
                nummodesmax = max(nummodesmax, locExp->GetBasisNumModes(2));

                Shapes[locExp->DetShapeType()] = 1;
            }


            m_comm->AllReduce(nummodesmax, ReduceMax);
            m_comm->AllReduce(Shapes, ReduceMax);

            if(Shapes[ePyramid]) // if Pyramids used then need Tet and Hex expansion
	    {
                Shapes[eTetrahedron] = 1;
                Shapes[eHexahedron]  = 1;
            }

            StdRegions::MatrixType PreconR;
            if(linSysKey.GetMatrixType() == StdRegions::eMass)
            {
                PreconR  = StdRegions::ePreconRMass;
            }
            else
            {
                PreconR  = StdRegions::ePreconR;
            }

            Array<OneD, unsigned int>  vmap;
            Array<OneD, Array<OneD, unsigned int> > emap;
            Array<OneD, Array<OneD, unsigned int> > fmap;

            /*
             * Set up a transformation matrices for equal max order
             * polynomial meshes
             */

            if(Shapes[eHexahedron])
            {
                SpatialDomains::HexGeomSharedPtr   hexgeom   = CreateRefHexGeom();
                //Bases for Hex element
                const BasisKey HexBa(eModified_A, nummodesmax,
                                 PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey HexBb(eModified_A, nummodesmax,
                                 PointsKey(nummodesmax+1,  eGaussLobattoLegendre));
                const BasisKey HexBc(eModified_A, nummodesmax,
                               PointsKey(nummodesmax+1,  eGaussLobattoLegendre));

                //Create reference Hexahdedral expansion
                LocalRegions::HexExpSharedPtr HexExp;

                HexExp = MemoryManager<LocalRegions::HexExp>
                    ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                        hexgeom);

                maxElmt[eHexahedron] = HexExp;

                // Hex:
                HexExp->GetInverseBoundaryMaps(vmap,emap,fmap);
                vertMapMaxR[eHexahedron] = vmap;
                edgeMapMaxR[eHexahedron] = emap;
                faceMapMaxR[eHexahedron] = fmap;

                //Get hexahedral transformation matrix
                LocalRegions::MatrixKey HexR
                    (PreconR, eHexahedron,
                     *HexExp, linSysKey.GetConstFactors());
                maxRmat[eHexahedron] = HexExp->GetLocMatrix(HexR);
            }

            if(Shapes[eTetrahedron])
            {
                SpatialDomains::TetGeomSharedPtr   tetgeom   = CreateRefTetGeom();
                //Bases for Tetrahedral element
                const BasisKey TetBa(eModified_A, nummodesmax,
                                     PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey TetBb(eModified_B, nummodesmax,
                                     PointsKey(nummodesmax,  eGaussRadauMAlpha1Beta0));
                const BasisKey TetBc(eModified_C, nummodesmax,
                                     PointsKey(nummodesmax,  eGaussRadauMAlpha2Beta0));

                //Create reference tetrahedral expansion
                LocalRegions::TetExpSharedPtr TetExp;

                TetExp = MemoryManager<LocalRegions::TetExp>
                    ::AllocateSharedPtr(TetBa,TetBb,TetBc,tetgeom);

                maxElmt[eTetrahedron] = TetExp;

                TetExp->GetInverseBoundaryMaps(vmap,emap,fmap);
                vertMapMaxR[eTetrahedron] = vmap;
                edgeMapMaxR[eTetrahedron] = emap;
                faceMapMaxR[eTetrahedron] = fmap;

                //Get tetrahedral transformation matrix
                LocalRegions::MatrixKey TetR
                    (PreconR, eTetrahedron,
                     *TetExp, linSysKey.GetConstFactors());
                maxRmat[eTetrahedron] = TetExp->GetLocMatrix(TetR);

                if((Shapes[ePyramid])||(Shapes[eHexahedron]))
                {
                    ReSetTetMaxRMat(nummodesmax, TetExp, maxRmat,
                                    vertMapMaxR, edgeMapMaxR, faceMapMaxR);
                }
            }

            if(Shapes[ePyramid])
            {
                SpatialDomains::PyrGeomSharedPtr   pyrgeom   = CreateRefPyrGeom();

                //Bases for Pyramid element
                const BasisKey PyrBa(eModified_A, nummodesmax,
                                     PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey PyrBb(eModified_A, nummodesmax,
                                     PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey PyrBc(eModifiedPyr_C, nummodesmax,
                                     PointsKey(nummodesmax,  eGaussRadauMAlpha2Beta0));

                //Create reference pyramid expansion
                LocalRegions::PyrExpSharedPtr PyrExp;

                PyrExp = MemoryManager<LocalRegions::PyrExp>
                    ::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,pyrgeom);

                maxElmt[ePyramid] = PyrExp;

                // Pyramid:
                PyrExp->GetInverseBoundaryMaps(vmap,emap,fmap);
                vertMapMaxR[ePyramid] = vmap;
                edgeMapMaxR[ePyramid] = emap;
                faceMapMaxR[ePyramid] = fmap;

                // Set up Pyramid Transformation Matrix based on Tet
                // and Hex expansion
                SetUpPyrMaxRMat(nummodesmax,PyrExp,maxRmat,vertMapMaxR,
                                edgeMapMaxR,faceMapMaxR);
            }

            if(Shapes[ePrism])
            {
                SpatialDomains::PrismGeomSharedPtr prismgeom = CreateRefPrismGeom();
                //Bases for Prism element
                const BasisKey PrismBa(eModified_A, nummodesmax,
                                  PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey PrismBb(eModified_A, nummodesmax,
                                  PointsKey(nummodesmax+1, eGaussLobattoLegendre));
                const BasisKey PrismBc(eModified_B, nummodesmax,
                                  PointsKey(nummodesmax,   eGaussRadauMAlpha1Beta0));

                //Create reference prismatic expansion
                LocalRegions::PrismExpSharedPtr PrismExp;

                PrismExp = MemoryManager<LocalRegions::PrismExp>
                    ::AllocateSharedPtr(PrismBa,PrismBb,PrismBc,prismgeom);
                maxElmt[ePrism] = PrismExp;

                // Prism:
                PrismExp->GetInverseBoundaryMaps(vmap,emap,fmap);
                vertMapMaxR[ePrism] = vmap;
                edgeMapMaxR[ePrism] = emap;

                faceMapMaxR[ePrism] = fmap;

                if((Shapes[ePyramid])||(Shapes[eHexahedron]))
                {
		    ReSetPrismMaxRMat(nummodesmax, PrismExp, maxRmat,
                                      vertMapMaxR, edgeMapMaxR,
                                      faceMapMaxR, false);
                }
                else
                {

                    //Get prismatic transformation matrix
                    LocalRegions::MatrixKey PrismR
                        (PreconR, ePrism,
                         *PrismExp, linSysKey.GetConstFactors());
                    maxRmat[ePrism] =
                        PrismExp->GetLocMatrix(PrismR);

                    if(Shapes[eTetrahedron]) // reset triangular faces from Tet
                    {
                        ReSetPrismMaxRMat(nummodesmax, PrismExp, maxRmat,
                                          vertMapMaxR, edgeMapMaxR,
                                          faceMapMaxR, true);
                    }
                }
            }
        }

        void PreconditionerLowEnergy::SetUpPyrMaxRMat(int nummodesmax,
                             LocalRegions::PyrExpSharedPtr &PyrExp,
                             std::map<ShapeType, DNekScalMatSharedPtr> &maxRmat,
                             std::map<ShapeType, Array<OneD, unsigned int> >        &vertMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &edgeMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &faceMapMaxR)
        {
            int nRows = PyrExp->NumBndryCoeffs();
            NekDouble val;
            NekDouble zero = 0.0;
            DNekMatSharedPtr newmat = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);

            // set diagonal to 1
            for(int i = 0; i < nRows; ++i)
            {
                newmat->SetValue(i,i,1.0);
            }

            // The following lists specify the number of adjacent
            // edges to each vertex (nadj) then the Hex vert to
            // use for each pyramid ver in the vert-edge map (VEVert)
            // followed by the hex edge to use for each pyramid edge
            // in the vert-edge map (VEEdge)
            const int nadjedge[]     = {3,3,3,3,4};
            const int VEHexVert[][4] = {{0,0,0,-1},{1,1,1,-1},{2,2,2,-1},{3,3,3,-1},{4,5,6,7}};
            const int VEHexEdge[][4] = {{0,3,4,-1},{0,1,5,-1},{1,2,6,-1},{2,3,7,-1},{4,5,6,7}};
            const int VEPyrEdge[][4] = {{0,3,4,-1},{0,1,5,-1},{1,2,6,-1},{2,3,7,-1},{4,5,6,7}};

            // fill vertex to edge coupling
            for(int v = 0; v < 5; ++v)
            {
                for(int e = 0; e < nadjedge[v]; ++e)
                {
                    for(int i = 0; i < nummodesmax-2; ++i)
                    {
                        // Note this is using wrong shape but gives
                        // answer that seems to have correct error!
                        val = (*maxRmat[eHexahedron])(
                                    vertMapMaxR[eHexahedron][VEHexVert[v][e]],
                                    edgeMapMaxR[eHexahedron][VEHexEdge[v][e]][i]);
                        newmat->SetValue(vertMapMaxR[ePyramid][v],
                                         edgeMapMaxR[ePyramid][VEPyrEdge[v][e]][i],val);
                    }
                }
            }

            int nfacemodes;
            nfacemodes = (nummodesmax-2)*(nummodesmax-2);
            // First four verties are all adjacent to base face
            for(int v = 0; v < 4; ++v)
            {
                for(int i = 0; i < nfacemodes; ++i)
                {
                    val = (*maxRmat[eHexahedron])(vertMapMaxR[eHexahedron][v],
                                                  faceMapMaxR[eHexahedron][0][i]);
                    newmat->SetValue(vertMapMaxR[ePyramid][v],
                                     faceMapMaxR[ePyramid][0][i],val);
                }
            }


            const int nadjface[]     = {2,2,2,2,4};
            const int VFTetVert[][4] = {{0,0,-1,-1},{1,1,-1,-1},{2,2,-1,-1},{0,2,-1,-1},{3,3,3,3}};
            const int VFTetFace[][4] = {{1,3,-1,-1},{1,2,-1,-1},{2,3,-1,-1},{1,3,-1,-1},{1,2,1,2}};
            const int VFPyrFace[][4] = {{1,4,-1,-1},{1,2,-1,-1},{2,3,-1,-1},{3,4,-1,-1},{1,2,3,4}};

            // next handle all triangular faces from tetrahedron
            nfacemodes = (nummodesmax-3)*(nummodesmax-2)/2;
            for(int v = 0; v < 5; ++v)
            {
                for(int f = 0; f < nadjface[v]; ++f)
                {
                    for(int i = 0; i < nfacemodes; ++i)
                    {
                        val = (*maxRmat[eTetrahedron])(vertMapMaxR[eTetrahedron][VFTetVert[v][f]],
                                                       faceMapMaxR[eTetrahedron][VFTetFace[v][f]][i]);
                        newmat->SetValue(vertMapMaxR[ePyramid][v],
                                         faceMapMaxR[ePyramid][VFPyrFace[v][f]][i],val);
                    }

                }
            }

            // Edge to face coupling
            // all base edges are coupled to face 0
            nfacemodes = (nummodesmax-2)*(nummodesmax-2);
            for(int e = 0; e < 4; ++e)
            {
                for(int i = 0; i < nummodesmax-2; ++i)
                {
                    for(int j = 0; j < nfacemodes; ++j)
                    {
                        int edgemapid = edgeMapMaxR[eHexahedron][e][i];
                        int facemapid = faceMapMaxR[eHexahedron][0][j];

                        val = (*maxRmat[eHexahedron])(edgemapid,facemapid);
                        newmat->SetValue(edgeMapMaxR[ePyramid][e][i],
                                         faceMapMaxR[ePyramid][0][j],val);
                    }

                }
            }

            const int nadjface1[]    = {1,1,1,1,2,2,2,2};
            const int EFTetEdge[][2] = {{0,-1},{1,-1},{0,-1},{2,-1},{3,3},{4,4},{5,5},{3,5}};
            const int EFTetFace[][2] = {{1,-1},{2,-1},{1,-1},{3,-1},{1,3},{1,2},{2,3},{1,3}};
            const int EFPyrFace[][2] = {{1,-1},{2,-1},{3,-1},{4,-1},{1,4},{1,2},{2,3},{3,4}};

            // next handle all triangular faces from tetrahedron
            nfacemodes = (nummodesmax-3)*(nummodesmax-2)/2;
            for(int e = 0; e < 8; ++e)
            {
                for(int f = 0; f < nadjface1[e]; ++f)
                {
                    for(int i = 0; i < nummodesmax-2; ++i)
                    {
                        for(int j = 0; j < nfacemodes; ++j)
                        {
                            int edgemapid = edgeMapMaxR[eTetrahedron][EFTetEdge[e][f]][i];
                            int facemapid = faceMapMaxR[eTetrahedron][EFTetFace[e][f]][j];

                            val = (*maxRmat[eTetrahedron])(edgemapid,facemapid);
                            newmat->SetValue(edgeMapMaxR[ePyramid][e][i],
                                             faceMapMaxR[ePyramid][EFPyrFace[e][f]][j],val);
                        }
                    }
                }
            }

            DNekScalMatSharedPtr PyrR;
            PyrR = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, newmat);
            maxRmat[ePyramid] =PyrR;
        }


        void PreconditionerLowEnergy::ReSetTetMaxRMat(int nummodesmax,
                             LocalRegions::TetExpSharedPtr &TetExp,
                             std::map<ShapeType, DNekScalMatSharedPtr> &maxRmat,
                             std::map<ShapeType, Array<OneD, unsigned int> >        &vertMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &edgeMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &faceMapMaxR)
        {
            boost::ignore_unused(faceMapMaxR);

            int nRows = TetExp->NumBndryCoeffs();
            NekDouble val;
            NekDouble zero = 0.0;
            DNekMatSharedPtr newmat = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);

            // copy existing system
            for(int i = 0; i < nRows; ++i)
            {
                for(int j = 0; j < nRows; ++j)
                {
                    val = (*maxRmat[eTetrahedron])(i,j);
                    newmat->SetValue(i,j,val);
                }
            }

            // The following lists specify the number of adjacent
            // edges to each vertex (nadj) then the Hex vert to
            // use for each pyramid ver in the vert-edge map (VEVert)
            // followed by the hex edge to use for each Tet edge
            // in the vert-edge map (VEEdge)
            const int VEHexVert[][4] = {{0,0,0},{1,1,1},{2,2,2},{4,5,6}};
            const int VEHexEdge[][4] = {{0,3,4},{0,1,5},{1,2,6},{4,5,6}};
            const int VETetEdge[][4] = {{0,2,3},{0,1,4},{1,2,5},{3,4,5}};

            // fill vertex to edge coupling
            for(int v = 0; v < 4; ++v)
            {
                for(int e = 0; e < 3; ++e)
                {
                    for(int i = 0; i < nummodesmax-2; ++i)
                    {
                        // Note this is using wrong shape but gives
                        // answer that seems to have correct error!
                        val = (*maxRmat[eHexahedron])(
                                    vertMapMaxR[eHexahedron][VEHexVert[v][e]],
                                    edgeMapMaxR[eHexahedron][VEHexEdge[v][e]][i]);
                        newmat->SetValue(vertMapMaxR[eTetrahedron][v],
                                    edgeMapMaxR[eTetrahedron][VETetEdge[v][e]][i],
                                                        val);
                    }
                }
            }

            DNekScalMatSharedPtr TetR =
                MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, newmat);

            maxRmat[eTetrahedron] = TetR;
        }

        void PreconditionerLowEnergy::ReSetPrismMaxRMat(int nummodesmax,
                             LocalRegions::PrismExpSharedPtr &PrismExp,
                             std::map<ShapeType, DNekScalMatSharedPtr> &maxRmat,
                             std::map<ShapeType, Array<OneD, unsigned int> >        &vertMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &edgeMapMaxR,
                             std::map<ShapeType, Array<OneD, Array<OneD, unsigned int> > > &faceMapMaxR,
                             bool UseTetOnly)
        {
            int nRows = PrismExp->NumBndryCoeffs();
            NekDouble val;
            NekDouble zero = 0.0;
            DNekMatSharedPtr newmat = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);


            int nfacemodes;

            if(UseTetOnly)
            {
                // copy existing system
                for(int i = 0; i < nRows; ++i)
                {
                    for(int j = 0; j < nRows; ++j)
                    {
                        val = (*maxRmat[ePrism])(i,j);
                        newmat->SetValue(i,j,val);
                    }
                }

                // Reset vertex to edge mapping from tet.
                const int VETetVert[][2]   = {{0,0},{1,1},{1,1},{0,0},{3,3},{3,3}};
                const int VETetEdge[][2]   = {{0,3},{0,4},{0,4},{0,3},{3,4},{4,3}};
                const int VEPrismEdge[][2] = {{0,4},{0,5},{2,6},{2,7},{4,5},{6,7}};

                // fill vertex to edge coupling
                for(int v = 0; v < 6; ++v)
                {
                    for(int e = 0; e < 2; ++e)
                    {
                        for(int i = 0; i < nummodesmax-2; ++i)
                        {
                            // Note this is using wrong shape but gives
                            // answer that seems to have correct error!
                            val = (*maxRmat[eTetrahedron])(
                                   vertMapMaxR[eTetrahedron][VETetVert[v][e]],
                                   edgeMapMaxR[eTetrahedron][VETetEdge[v][e]][i]);
                            newmat->
                                SetValue(vertMapMaxR[ePrism][v],
                                         edgeMapMaxR[ePrism][VEPrismEdge[v][e]][i],
                                         val);
                        }
                    }
                }
            }
            else
            {

                // set diagonal to 1
                for(int i = 0; i < nRows; ++i)
                {
                    newmat->SetValue(i,i,1.0);
                }


                // Set vertex to edge mapping from Hex.

                // The following lists specify the number of adjacent
                // edges to each vertex (nadj) then the Hex vert to
                // use for each prism ver in the vert-edge map (VEVert)
                // followed by the hex edge to use for each prism edge
                // in the vert-edge map (VEEdge)
                const int VEHexVert[][3]   = {{0,0,0},{1,1,1},{2,2,2},{3,3,3},
                                              {4,5,5},{6,7,7}};
                const int VEHexEdge[][3]   = {{0,3,4},{0,1,5},{1,2,6},{2,3,7},
                                              {4,5,9},{6,7,11}};
                const int VEPrismEdge[][3] = {{0,3,4},{0,1,5},{1,2,6},{2,3,7},
                                              {4,5,8},{6,7,8}};

                // fill vertex to edge coupling
                for(int v = 0; v < 6; ++v)
                {
                    for(int e = 0; e < 3; ++e)
                    {
                        for(int i = 0; i < nummodesmax-2; ++i)
                        {
                            // Note this is using wrong shape but gives
                            // answer that seems to have correct error!
                            val = (*maxRmat[eHexahedron])(
                                    vertMapMaxR[eHexahedron][VEHexVert[v][e]],
                                    edgeMapMaxR[eHexahedron][VEHexEdge[v][e]][i]);
                            newmat->SetValue(vertMapMaxR[ePrism][v],
                                             edgeMapMaxR[ePrism][VEPrismEdge[v][e]][i],
                                             val);
                        }
                    }
                }


                // Setup vertex to face mapping from Hex
                const int VFHexVert[][2]   = {{0,0},{1,1},{4,5},{2,2},{3,3},{6,7}};
                const int VFHexFace[][2]   = {{0,4},{0,2},{4,2},{0,2},{0,4},{2,4}};

                const int VQFPrismVert[][2] = {{0,0},{1,1},{4,4},{2,2},{3,3},{5,5}};
                const int VQFPrismFace[][2] = {{0,4},{0,2},{4,2},{0,2},{0,4},{2,4}};

                nfacemodes = (nummodesmax-2)*(nummodesmax-2);
                // Replace two Quad faces  on every vertex
                for(int v = 0; v < 6; ++v)
                {
                    for(int f = 0; f < 2; ++f)
                    {
                        for(int i = 0; i < nfacemodes; ++i)
                        {
                            val = (*maxRmat[eHexahedron])(
                                             vertMapMaxR[eHexahedron][VFHexVert[v][f]],
                                             faceMapMaxR[eHexahedron][VFHexFace[v][f]][i]);
                            newmat->SetValue(vertMapMaxR[ePrism][VQFPrismVert[v][f]],
                                             faceMapMaxR[ePrism][VQFPrismFace[v][f]][i],val);
                        }
                    }
                }

                // Mapping of Hex Edge-Face mappings to Prism Edge-Face Mappings
                const int nadjface[] = {1,2,1,2,1,1,1,1,2};
                const int EFHexEdge[][2]    = {{0,-1},{1,1},{2,-1},{3,3},{4,-1},{5,-1},{6,-1},{7,-1},{9,11}};
                const int EFHexFace[][2]    = {{0,-1},{0,2},{0,-1},{0,4},{4,-1},{2,-1},{2,-1},{4,-1},{2,4}};
                const int EQFPrismEdge[][2] = {{0,-1},{1,1},{2,-1},{3,3},{4,-1},{5,-1},{6,-1},{7,-1},{8,8}};
                const int EQFPrismFace[][2] = {{0,-1},{0,2},{0,-1},{0,4},{4,-1},{2,-1},{2,-1},{4,-1},{2,4}};

                // all base edges are coupled to face 0
                nfacemodes = (nummodesmax-2)*(nummodesmax-2);
                for(int e = 0; e < 9; ++e)
                {
                    for(int f = 0; f < nadjface[e]; ++f)
                    {
                        for(int i = 0; i < nummodesmax-2; ++i)
                        {
                            for(int j = 0; j < nfacemodes; ++j)
                            {
                                int edgemapid = edgeMapMaxR[eHexahedron][EFHexEdge[e][f]][i];
                                int facemapid = faceMapMaxR[eHexahedron][EFHexFace[e][f]][j];

                                val = (*maxRmat[eHexahedron])(edgemapid,facemapid);

                                int edgemapid1 = edgeMapMaxR[ePrism][EQFPrismEdge[e][f]][i];
                                int facemapid1 = faceMapMaxR[ePrism][EQFPrismFace[e][f]][j];
                                newmat->SetValue(edgemapid1, facemapid1, val);
                            }
                        }
                    }
                }
            }

            const int VFTetVert[]    = {0,1,3,1,0,3};
            const int VFTetFace[]    = {1,1,1,1,1,1};
            const int VTFPrismVert[] = {0,1,4,2,3,5};
            const int VTFPrismFace[] = {1,1,1,3,3,3};

            //  Handle all triangular faces from tetrahedron
            nfacemodes = (nummodesmax-3)*(nummodesmax-2)/2;
            for(int v = 0; v < 6; ++v)
            {
                for(int i = 0; i < nfacemodes; ++i)
                {
                    val = (*maxRmat[eTetrahedron])
                        (vertMapMaxR[eTetrahedron][VFTetVert[v]],
                         faceMapMaxR[eTetrahedron][VFTetFace[v]][i]);

                    newmat->SetValue(vertMapMaxR[ePrism][VTFPrismVert[v]],
                                     faceMapMaxR[ePrism][VTFPrismFace[v]][i],val);
                }
            }

            // Mapping of Tet Edge-Face mappings to Prism Edge-Face Mappings
            const int EFTetEdge[]    = {0,3,4,0,4,3};
            const int EFTetFace[]    = {1,1,1,1,1,1};
            const int ETFPrismEdge[] = {0,4,5,2,6,7};
            const int ETFPrismFace[] = {1,1,1,3,3,3};

            // handle all edge to triangular faces from tetrahedron
            // (only 6 this time)
            nfacemodes = (nummodesmax-3)*(nummodesmax-2)/2;
            for(int e = 0; e < 6; ++e)
            {
                for(int i = 0; i < nummodesmax-2; ++i)
                {
                    for(int j = 0; j < nfacemodes; ++j)
                    {
                        int edgemapid = edgeMapMaxR[eTetrahedron][EFTetEdge[e]][i];
                        int facemapid = faceMapMaxR[eTetrahedron][EFTetFace[e]][j];
                        val = (*maxRmat[eTetrahedron])(edgemapid,facemapid);

                        newmat->SetValue(edgeMapMaxR[ePrism][ETFPrismEdge[e]][i],
                                         faceMapMaxR[ePrism][ETFPrismFace[e]][j],val);
                    }
                }
            }


            DNekScalMatSharedPtr PrismR
                = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, newmat);
            maxRmat[ePrism] = PrismR;
        }

        DNekMatSharedPtr PreconditionerLowEnergy::
        ExtractLocMat(StdRegions::StdExpansionSharedPtr &locExp,
                      DNekScalMatSharedPtr              &maxRmat,
                      LocalRegions::ExpansionSharedPtr  &maxExp,
                      Array<OneD, unsigned int>         &vmap,
                      Array<OneD, Array<OneD, unsigned int> > &emap)
        {
            NekDouble val;
            NekDouble zero = 0.0;

            int nRows = locExp->NumBndryCoeffs();
            DNekMatSharedPtr newmat = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);

            Array<OneD, unsigned int>  vlocmap;
            Array<OneD, Array<OneD, unsigned int> > elocmap;
            Array<OneD, Array<OneD, unsigned int> > flocmap;

            locExp->GetInverseBoundaryMaps(vlocmap,elocmap,flocmap);

            // fill diagonal
            for(int i = 0; i < nRows; ++i)
            {
                val = 1.0;
                newmat->SetValue(i,i,val);
            }

            int nverts = locExp->GetNverts();
            int nedges = locExp->GetNedges();
            int nfaces = locExp->GetNfaces();

            // fill vertex to edge coupling
            for(int e = 0; e < nedges; ++e)
            {
                int nEdgeInteriorCoeffs = locExp->GetEdgeNcoeffs(e) -2;

                for(int v = 0; v < nverts; ++v)
                {
                    for(int i = 0; i < nEdgeInteriorCoeffs; ++i)
                    {
                        val = (*maxRmat)(vmap[v],emap[e][i]);
                        newmat->SetValue(vlocmap[v],elocmap[e][i],val);
                    }
                }
            }

            for(int f = 0; f < nfaces; ++f)
            {
                // Get details to extrac this face from max reference matrix
                StdRegions::Orientation FwdOrient = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
                int m0,m1; //Local face expansion orders.

                int nFaceInteriorCoeffs = locExp->GetFaceIntNcoeffs(f);

                locExp->GetFaceNumModes(f,FwdOrient,m0,m1);

                Array<OneD, unsigned int> fmapRmat = maxExp->
                    GetFaceInverseBoundaryMap(f,FwdOrient, m0,m1);

                // fill in vertex to face coupling
                for(int v = 0; v < nverts; ++v)
                {
                    for(int i = 0; i < nFaceInteriorCoeffs; ++i)
                    {
                        val = (*maxRmat)(vmap[v],fmapRmat[i]);
                        newmat->SetValue(vlocmap[v],flocmap[f][i],val);
                    }
                }

                // fill in edges to face coupling
                for(int e = 0; e < nedges; ++e)
                {
                    int nEdgeInteriorCoeffs = locExp->GetEdgeNcoeffs(e) -2;

                    for(int j = 0; j < nEdgeInteriorCoeffs; ++j)
                    {

                        for(int i = 0; i < nFaceInteriorCoeffs; ++i)
                        {
                            val = (*maxRmat)(emap[e][j],fmapRmat[i]);
                            newmat->SetValue(elocmap[e][j],flocmap[f][i],val);
                        }
                    }
                }
            }

            return newmat;
        }
    }
}

