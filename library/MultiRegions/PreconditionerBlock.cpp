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
#include <math.h>

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
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
                         const AssemblyMapSharedPtr &pLocToGloMap)
           : Preconditioner(plinsys, pLocToGloMap),
             m_linsys(plinsys),
             m_preconType(pLocToGloMap->GetPreconType()),
             m_locToGloMap(pLocToGloMap)
         {
         }

        void PreconditionerBlock::v_InitObject()
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            ASSERTL0(solvertype == MultiRegions::eIterativeStaticCond,"Solver type not valid");
        }


        void PreconditionerBlock::v_BuildPreconditioner()
        {
            // Different setup for HDG
            GlobalLinSysKey key = m_linsys.lock()->GetKey();
            if (key.GetMatrixType() == StdRegions::eHybridDGHelmBndLam)
            {
                BlockPreconditionerHDG();
                return;
            }

            boost::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            StdRegions::StdExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(0);
            int nDim = locExpansion->GetShapeDimension();

            if (nDim == 1)
            {
                ASSERTL0(0,"Unknown preconditioner");
            }
            else if (nDim == 2)
            {
                BlockPreconditioner2D();
            }
            else if (nDim == 3)
            {
                BlockPreconditioner3D();
            }
        }

       /**
         * \brief Construct a block preconditioner from
         * \f$\mathbf{S}_{1}\f$
         *
         *\f[\mathbf{M}^{-1}=\left[\begin{array}{ccc}
         *  Diag[(\mathbf{S_{1}})_{vv}] & & \\ & (\mathbf{S}_{1})_{eb} & \\ & &
         *  (\mathbf{S}_{1})_{fb} \end{array}\right] \f]
         *
         * where \f$\mathbf{S}_{1}\f$ is the local schur complement matrix for
         * each element.
         */
        void PreconditionerBlock::BlockPreconditioner2D()
        {
            boost::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            LocalRegions::ExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            int i, j, k, nel;
            int nVerts, nEdges;
            int eid, n, cnt, nedgemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2;
            int m, v, eMap1, eMap2;
            int offset, globalrow, globalcol;

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr VertBlk;

            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

            //Vertex and edge preconditioner matrices
            VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> VertBlockToUniversalMap(nNonDirVerts,-1);

            int n_exp = expList->GetNumElmts();
            int nNonDirEdgeIDs=m_locToGloMap->GetNumNonDirEdges();

            //set the number of blocks in the matrix
            Array<OneD,unsigned int> n_blks(1+nNonDirEdgeIDs);
            n_blks[0]=nNonDirVerts;

            map<int,int> edgeDirMap;
            map<int,int> uniqueEdgeMap;

            //this should be of size total number of local edges
            Array<OneD, int> edgemodeoffset(nNonDirEdgeIDs,0);
            Array<OneD, int> edgeglobaloffset(nNonDirEdgeIDs,0);

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            StdRegions::StdExpansion1DSharedPtr bndCondFaceExp;
            LocalRegions::SegExpSharedPtr       bndSegExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>&
                bndConditions = expList->GetBndConditions();

            int meshVertId;
            int meshEdgeId;

            //Determine which boundary edges and faces have dirichlet values
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndSegExp = bndCondExp[i]->GetExp(j)
                                        ->as<LocalRegions::SegExp>();
                    if (bndConditions[i]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        meshEdgeId = (bndSegExp->GetGeom1D())->GetEid();
                        edgeDirMap[meshEdgeId] = 1;
                    }
                }
            }

            int dof=0;
            int maxEdgeDof=0;
            int nlocalNonDirEdges=0;

            int edgematrixlocation=0;
            int ntotaledgeentries=0;

            // Loop over all the elements in the domain and compute max edge
            // DOF. Reduce across all processes to get universal maximum.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                for (j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    dof    = locExpansion->GetEdgeNcoeffs(j)-2;
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion2D>()
                                    ->GetGeom2D()->GetEid(j);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        if(uniqueEdgeMap.count(meshEdgeId)==0)
                        {
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;

                            edgeglobaloffset[edgematrixlocation]+=ntotaledgeentries;

                            edgemodeoffset[edgematrixlocation]=dof*dof;

                            ntotaledgeentries+=dof*dof;

                            n_blks[1+edgematrixlocation++]=dof;

                        }

                        nlocalNonDirEdges+=dof*dof;
                    }
                }
            }

            m_comm = expList->GetComm()->GetRowComm();
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);

            //Allocate arrays for block to universal map (number of expansions * p^2)
            Array<OneD, long> EdgeBlockToUniversalMap(ntotaledgeentries,-1);

            Array<OneD, int> localEdgeToGlobalMatrixMap(nlocalNonDirEdges,-1);

            //Allocate arrays to store matrices (number of expansions * p^2)
            Array<OneD, NekDouble> EdgeBlockArray(nlocalNonDirEdges,-1);

            int edgematrixoffset=0;
            int vGlobal;

            for(n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion2D>()
                                            ->GetGeom2D()->GetEid(j);

                    nedgemodes=locExpansion->GetEdgeNcoeffs(j)-2;

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=edgeglobaloffset[uniqueEdgeMap[meshEdgeId]]+k;


                            localEdgeToGlobalMatrixMap[edgematrixoffset+k]=vGlobal;

                            EdgeBlockToUniversalMap[vGlobal]
                                = meshEdgeId * maxEdgeDof * maxEdgeDof + k + 1;
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }
            }

            edgematrixoffset=0;

            m_blkMat = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                nVerts=locExpansion->GetGeom()->GetNumVerts();
                nEdges=locExpansion->GetGeom()->GetNumEdges();

                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

                //Extract boundary block (elemental S1)
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                offset = bnd_mat->GetRows();

                DNekScalMat &S=(*bnd_mat);

                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    vMap1=locExpansion->GetVertexMap(v);

                    //Get vertex map
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=locExpansion->GetVertexMap(m);

                            //global matrix location (without offset due to
                            //dirichlet values)
                            globalcol = m_locToGloMap->
                                GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol == globalrow)
                            {
                                meshVertId = locExpansion->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetVid(v);

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);

                                vertArray[globalrow]
                                    += sign1*sign2*S(vMap1,vMap2);

                                VertBlockToUniversalMap[globalrow]
                                = meshVertId * maxEdgeDof * maxEdgeDof + 1;
                            }
                        }
                    }
                }

                //loop over edges of the element and return the edge map
                for (eid=0; eid<nEdges; ++eid)
                {
                    nedgemodes=locExpansion->GetEdgeNcoeffs(eid)-2;

                    DNekMatSharedPtr locMat =
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nedgemodes,nedgemodes,zero,storage);

                    meshEdgeId = locExpansion->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetEid(eid);
                    Array<OneD, unsigned int> edgemodearray =
                        locExpansion->GetEdgeInverseBoundaryMap(eid);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        for (v=0; v<nedgemodes; ++v)
                        {
                            eMap1=edgemodearray[v];

                            for (m=0; m<nedgemodes; ++m)
                            {
                                eMap2=edgemodearray[m];

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + eMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + eMap2);

                                NekDouble globalEdgeValue =
                                    sign1*sign2*S(eMap1,eMap2);

                                EdgeBlockArray[edgematrixoffset+v*nedgemodes+m]=
                                    globalEdgeValue;
                            }
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }

                //offset for the expansion
                cnt+=offset;
            }

            //Assemble edge matrices of each process
            Array<OneD, NekDouble> GlobalEdgeBlock(ntotaledgeentries);
            Vmath::Zero(ntotaledgeentries, GlobalEdgeBlock.get(), 1);
            Vmath::Assmb(EdgeBlockArray.num_elements(),
                         EdgeBlockArray.get(),
                         localEdgeToGlobalMatrixMap.get(),
                         GlobalEdgeBlock.get());

            //Exchange vertex data over different processes
            if(nNonDirVerts!=0)
            {
                Gs::gs_data *tmp = Gs::Init(VertBlockToUniversalMap, m_comm);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
            }

            //Exchange edge data over different processes
            Gs::gs_data *tmp1 = Gs::Init(EdgeBlockToUniversalMap, m_comm);
            Gs::Gather(GlobalEdgeBlock, Gs::gs_add, tmp1);

            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                  VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }

            //Set the first block to be the diagonal of the vertex space
            m_blkMat->SetBlock(0,0, VertBlk);

            offset=0;
            //Build the edge matrices from the vector
            for(int loc=0; loc<nNonDirEdgeIDs; ++loc)
            {
                DNekMatSharedPtr gmat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nedgemodes,nedgemodes,zero,storage);

                for (v=0; v<nedgemodes; ++v)
                {
                    for (m=0; m<nedgemodes; ++m)
                    {
                        NekDouble EdgeValue = GlobalEdgeBlock[offset+v*nedgemodes+m];
                        gmat->SetValue(v,m,EdgeValue);
                    }
                }

                m_blkMat->SetBlock(1+loc,1+loc, gmat);

                offset+=edgemodeoffset[loc];
            }

            int totblks=m_blkMat->GetNumberOfBlockRows();
            for (i=1; i< totblks; ++i)
            {
                unsigned int nmodes=m_blkMat->GetNumberOfRowsInBlockRow(i);
                DNekMatSharedPtr tmp_mat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                tmp_mat=m_blkMat->GetBlock(i,i);
                tmp_mat->Invert();
                m_blkMat->SetBlock(i,i,tmp_mat);
            }
        }

        /**
         *
         */
        void PreconditionerBlock::BlockPreconditioner3D()
        {
            boost::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            LocalRegions::ExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            int i, j, k, nel;
            int nVerts, nEdges,nFaces;
            int eid, fid, n, cnt, nedgemodes, nfacemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr VertBlk;

            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

            //Vertex, edge and face preconditioner matrices
            VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> VertBlockToUniversalMap(nNonDirVerts,-1);

            int n_exp = expList->GetNumElmts();
            int nNonDirEdgeIDs=m_locToGloMap->GetNumNonDirEdges();
            int nNonDirFaceIDs=m_locToGloMap->GetNumNonDirFaces();

            //set the number of blocks in the matrix
            Array<OneD,unsigned int> n_blks(1+nNonDirEdgeIDs+nNonDirFaceIDs);
            n_blks[0]=nNonDirVerts;

            map<int,int> edgeDirMap;
            map<int,int> faceDirMap;
            map<int,int> uniqueEdgeMap;
            map<int,int> uniqueFaceMap;

            //this should be of size total number of local edges
            Array<OneD, int> edgemodeoffset(nNonDirEdgeIDs,0);
            Array<OneD, int> facemodeoffset(nNonDirFaceIDs,0);

            Array<OneD, int> edgeglobaloffset(nNonDirEdgeIDs,0);
            Array<OneD, int> faceglobaloffset(nNonDirFaceIDs,0);

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>& bndConditions = expList->GetBndConditions();

            int meshVertId;
            int meshEdgeId;
            int meshFaceId;

            const Array<OneD, const int> &extradiredges
                = m_locToGloMap->GetExtraDirEdges();
            for(i=0; i<extradiredges.num_elements(); ++i)
            {
                meshEdgeId=extradiredges[i];
                edgeDirMap[meshEdgeId] = 1;
            }

            //Determine which boundary edges and faces have dirichlet values
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp = bndCondExp[i]->GetExp(j)->as<StdRegions::StdExpansion2D>();
                    if (bndConditions[i]->GetBoundaryConditionType() ==
                        SpatialDomains::eDirichlet)
                    {
                        for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                        {
                            meshEdgeId = bndCondFaceExp->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetEid(k);
                            if(edgeDirMap.count(meshEdgeId) == 0)
                            {
                                edgeDirMap[meshEdgeId] = 1;
                            }
                        }
                        meshFaceId = bndCondFaceExp->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetFid();
                        faceDirMap[meshFaceId] = 1;
                    }
                }
            }

            int dof=0;
            int maxFaceDof=0;
            int maxEdgeDof=0;
            int nlocalNonDirEdges=0;
            int nlocalNonDirFaces=0;

            int edgematrixlocation=0;
            int ntotaledgeentries=0;

            // Loop over all the elements in the domain and compute max edge
            // DOF. Reduce across all processes to get universal maximum.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                for (j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    dof    = locExpansion->GetEdgeNcoeffs(j)-2;
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        if(uniqueEdgeMap.count(meshEdgeId)==0 && dof > 0)
                        {
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;

                            edgeglobaloffset[edgematrixlocation]+=ntotaledgeentries;

                            edgemodeoffset[edgematrixlocation]=dof*dof;

                            ntotaledgeentries+=dof*dof;

                            n_blks[1+edgematrixlocation++]=dof;
                        }

                        nlocalNonDirEdges+=dof*dof;
                    }
                }
            }

            int facematrixlocation=0;
            int ntotalfaceentries=0;

            // Loop over all the elements in the domain and compute max face
            // DOF. Reduce across all processes to get universal maximum.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                for (j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    dof    = locExpansion->GetFaceIntNcoeffs(j);
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);

                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(j);

                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        if(uniqueFaceMap.count(meshFaceId)==0 && dof > 0)
                        {
                            uniqueFaceMap[meshFaceId]=facematrixlocation;

                            facemodeoffset[facematrixlocation]=dof*dof;

                            faceglobaloffset[facematrixlocation]+=ntotalfaceentries;

                            ntotalfaceentries+=dof*dof;

                            n_blks[1+nNonDirEdgeIDs+facematrixlocation++]=dof;
                        }
                        nlocalNonDirFaces+=dof*dof;
                    }

                }
            }

            m_comm = expList->GetComm();
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);
            m_comm->AllReduce(maxFaceDof, LibUtilities::ReduceMax);

            //Allocate arrays for block to universal map (number of expansions * p^2)
            Array<OneD, long> EdgeBlockToUniversalMap(ntotaledgeentries,-1);
            Array<OneD, long> FaceBlockToUniversalMap(ntotalfaceentries,-1);

            Array<OneD, int> localEdgeToGlobalMatrixMap(nlocalNonDirEdges,-1);
            Array<OneD, int> localFaceToGlobalMatrixMap(nlocalNonDirFaces,-1);

            //Allocate arrays to store matrices (number of expansions * p^2)
            Array<OneD, NekDouble> EdgeBlockArray(nlocalNonDirEdges,-1);
            Array<OneD, NekDouble> FaceBlockArray(nlocalNonDirFaces,-1);

            int edgematrixoffset=0;
            int facematrixoffset=0;
            int vGlobal;

            for(n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);

                    nedgemodes=locExpansion->GetEdgeNcoeffs(j)-2;

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=edgeglobaloffset[uniqueEdgeMap[meshEdgeId]]+k;


                            localEdgeToGlobalMatrixMap[edgematrixoffset+k]=vGlobal;

                            EdgeBlockToUniversalMap[vGlobal]
                                = meshEdgeId * maxEdgeDof * maxEdgeDof + k + 1;
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }

                //loop over the faces of the expansion
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    //get mesh face id
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(j);

                    nfacemodes = locExpansion->GetFaceIntNcoeffs(j);

                    //Check if face is has dirichlet values
                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        for(k=0; k<nfacemodes*nfacemodes; ++k)
                        {
                            vGlobal=faceglobaloffset[uniqueFaceMap[meshFaceId]]+k;

                            localFaceToGlobalMatrixMap[facematrixoffset+k]
                                = vGlobal;

                            FaceBlockToUniversalMap[vGlobal]
                                = meshFaceId * maxFaceDof * maxFaceDof + k + 1;
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                    }
                }
            }

            edgematrixoffset=0;
            facematrixoffset=0;

            m_blkMat = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

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

                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    vMap1=locExpansion->GetVertexMap(v);

                    //Get vertex map
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=locExpansion->GetVertexMap(m);

                            //global matrix location (without offset due to
                            //dirichlet values)
                            globalcol = m_locToGloMap->
                                GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol == globalrow)
                            {
                                meshVertId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetVid(v);

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);

                                vertArray[globalrow]
                                    += sign1*sign2*S(vMap1,vMap2);

                                VertBlockToUniversalMap[globalrow]
                                = meshVertId * maxEdgeDof * maxEdgeDof + 1;
                            }
                        }
                    }
                }

                //loop over edges of the element and return the edge map
                for (eid=0; eid<nEdges; ++eid)
                {
                    nedgemodes=locExpansion->GetEdgeNcoeffs(eid)-2;

                    DNekMatSharedPtr locMat =
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nedgemodes,nedgemodes,zero,storage);

                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(eid);
                    Array<OneD, unsigned int> edgemodearray = locExpansion->GetEdgeInverseBoundaryMap(eid);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        for (v=0; v<nedgemodes; ++v)
                        {
                            eMap1=edgemodearray[v];

                            for (m=0; m<nedgemodes; ++m)
                            {
                                eMap2=edgemodearray[m];

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + eMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + eMap2);

                                NekDouble globalEdgeValue = sign1*sign2*S(eMap1,eMap2);

                                EdgeBlockArray[edgematrixoffset+v*nedgemodes+m]=globalEdgeValue;
                            }
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }

                //loop over faces of the element and return the face map
                for (fid=0; fid<nFaces; ++fid)
                {
                    nfacemodes = locExpansion->GetFaceIntNcoeffs(fid);

                    DNekMatSharedPtr locMat =
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nfacemodes,nfacemodes,zero,storage);

                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(fid);

                    Array<OneD, unsigned int> facemodearray = locExpansion->GetFaceInverseBoundaryMap(fid);

                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        for (v=0; v<nfacemodes; ++v)
                        {
                            fMap1=facemodearray[v];

                            for (m=0; m<nfacemodes; ++m)
                            {
                                fMap2=facemodearray[m];

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + fMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + fMap2);

                                // Get the face-face value from the low energy matrix (S2)
                                NekDouble globalFaceValue = sign1*sign2*S(fMap1,fMap2);

                                //local face value to global face value
                                FaceBlockArray[facematrixoffset+v*nfacemodes+m]=globalFaceValue;
                            }
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                    }
                }

                //offset for the expansion
                cnt+=offset;
            }

            //Assemble edge matrices of each process
            Array<OneD, NekDouble> GlobalEdgeBlock(ntotaledgeentries);
            Vmath::Zero(ntotaledgeentries, GlobalEdgeBlock.get(), 1);
            Vmath::Assmb(EdgeBlockArray.num_elements(),
                         EdgeBlockArray.get(),
                         localEdgeToGlobalMatrixMap.get(),
                         GlobalEdgeBlock.get());

            //Assemble face matrices of each process
            Array<OneD, NekDouble> GlobalFaceBlock(ntotalfaceentries);
            Vmath::Zero(ntotalfaceentries, GlobalFaceBlock.get(), 1);
            Vmath::Assmb(FaceBlockArray.num_elements(),
                         FaceBlockArray.get(),
                         localFaceToGlobalMatrixMap.get(),
                         GlobalFaceBlock.get());

            //Exchange vertex data over different processes
            if(nNonDirVerts!=0)
            {
                Gs::gs_data *tmp = Gs::Init(VertBlockToUniversalMap, m_comm);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
            }

            //Exchange edge data over different processes
            Gs::gs_data *tmp1 = Gs::Init(EdgeBlockToUniversalMap, m_comm);
            Gs::Gather(GlobalEdgeBlock, Gs::gs_add, tmp1);

            //Exchange face data over different processes
            Gs::gs_data *tmp2 = Gs::Init(FaceBlockToUniversalMap, m_comm);
            Gs::Gather(GlobalFaceBlock, Gs::gs_add, tmp2);

            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                  VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }

            //Set the first block to be the diagonal of the vertex space
            m_blkMat->SetBlock(0,0, VertBlk);

            offset=0;
            //Build the edge matrices from the vector
            for(int loc=0; loc<nNonDirEdgeIDs; ++loc)
            {
                DNekMatSharedPtr gmat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nedgemodes,nedgemodes,zero,storage);

                for (v=0; v<nedgemodes; ++v)
                {
                    for (m=0; m<nedgemodes; ++m)
                    {
                        NekDouble EdgeValue = GlobalEdgeBlock[offset+v*nedgemodes+m];
                        gmat->SetValue(v,m,EdgeValue);
                    }
                }

                m_blkMat->SetBlock(1+loc,1+loc, gmat);

                offset+=edgemodeoffset[loc];
            }

            offset=0;
            //Build the face matrices from the vector
            for(int loc=0; loc<nNonDirFaceIDs; ++loc)
            {
                nfacemodes=n_blks[1+nNonDirEdgeIDs+loc];

                DNekMatSharedPtr gmat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nfacemodes,nfacemodes,zero,storage);

                for (v=0; v<nfacemodes; ++v)
                {
                    for (m=0; m<nfacemodes; ++m)
                    {
                        NekDouble FaceValue = GlobalFaceBlock[offset+v*nfacemodes+m];
                        gmat->SetValue(v,m,FaceValue);
                    }
                }

                m_blkMat->SetBlock(1+nNonDirEdgeIDs+loc,1+nNonDirEdgeIDs+loc, gmat);

                offset+=facemodeoffset[loc];
            }


            int totblks=m_blkMat->GetNumberOfBlockRows();
            for (i=1; i< totblks; ++i)
            {
                unsigned int nmodes=m_blkMat->GetNumberOfRowsInBlockRow(i);
                DNekMatSharedPtr tmp_mat =
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                tmp_mat=m_blkMat->GetBlock(i,i);
                tmp_mat->Invert();

                m_blkMat->SetBlock(i,i,tmp_mat);
            }
        }

        void PreconditionerBlock::BlockPreconditionerHDG()
        {
            boost::shared_ptr<MultiRegions::ExpList>
                expList=((m_linsys.lock())->GetLocMat()).lock();
            boost::shared_ptr<MultiRegions::ExpList> trace = expList->GetTrace();
            LocalRegions::ExpansionSharedPtr locExpansion;
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            AssemblyMapDGSharedPtr asmMap = boost::dynamic_pointer_cast<
                AssemblyMapDG>(m_locToGloMap);

            int i, j, k, n, cnt, cnt2;

            // Figure out number of Dirichlet trace elements
            int nTrace = expList->GetTrace()->GetExpSize();
            int nDir   = m_locToGloMap->GetNumGlobalDirBndCoeffs();

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
            LibUtilities::CommSharedPtr comm =
                expList->GetSession()->GetComm()->GetRowComm();
            comm->AllReduce(maxTraceSize, LibUtilities::ReduceMax);

            // Zero matrix storage.
            Array<OneD, NekDouble> tmpStore(cnt, 0.0);

            // Assemble block matrices for each trace element.
            for (cnt = n = 0; n < expList->GetExpSize(); ++n)
            {
                int elmt = expList->GetOffset_Elmt_Id(n);
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
            Array<OneD, long> uniIds(tmpStore.num_elements());
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
            Gs::gs_data *gsh = Gs::Init(uniIds, comm);
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
         *
         */
        void PreconditionerBlock::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nNonDir = nGlobal-nDir;
            DNekBlkMat &M = (*m_blkMat);
            NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
            NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
            z = M * r;
        }
    }
}
