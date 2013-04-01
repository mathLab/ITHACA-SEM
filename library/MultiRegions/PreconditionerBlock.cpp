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
                    "Block Preconditioning");
 
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
             m_locToGloMap(pLocToGloMap)
         {
	 }

        void PreconditionerBlock::v_InitObject()
        {
            BlockPreconditioner3D();
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
            StdRegions::StdExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;

            int i, j, k, nel;
            int nVerts, nEdges; 
            int eid, n, cnt;
            int nedgemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2;
            int m, v, eMap1, eMap2;
            int offset, globalrow, globalcol, nCoeffs;
            NekDouble vertValue;

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;
            
            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;

            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> m_VertBlockToUniversalMap(nNonDirVerts,-1);

            int n_exp = expList->GetNumElmts();

            map<int,int> edgeDirMap;

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
                    bndSegExp = boost::dynamic_pointer_cast<
                        LocalRegions::SegExp>(bndCondExp[i]->GetExp(j));
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
                    meshEdgeId = (locExpansion->GetGeom2D())->GetEid(j);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        nlocalNonDirEdges++;
                    }
                }
            }

            int nNonDirEdgeIDs=m_locToGloMap->GetNumNonDirEdges();

            //communicator
            m_comm = expList->GetComm()->GetRowComm();
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);

            //Allocate arrays for block to universal map (number of expansions * p^2)
            Array<OneD, long> m_EdgeBlockToUniversalMap(nNonDirEdgeIDs*maxEdgeDof*maxEdgeDof,-1);

            Array<OneD, int> m_localEdgeToGlobalMatrixMap(nlocalNonDirEdges*maxEdgeDof*maxEdgeDof,-1);

            //Allocate arrays to store matrices (number of expansions * p^2)
            Array<OneD, NekDouble> m_EdgeBlockArray(nlocalNonDirEdges*maxEdgeDof*maxEdgeDof,-1);

            //Set up mappings
            map<int,int> uniqueEdgeMap;

            //this should be of size total number of local edges
            Array<OneD, unsigned int> m_edgemodeoffset(nNonDirEdgeIDs);

            //set the number of blocks in the matrix
            Array<OneD,unsigned int> n_blks(1+nNonDirEdgeIDs);
            n_blks[0]=nNonDirVerts;

            int edgematrixlocation=0;
            int ecnt;
            int edgematrixoffset=0;
            int vGlobal;
            int ntotaledgeentries=0;

            for(ecnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->GetGeom2D()->GetEid(j);

                    nedgemodes=locExpansion->GetEdgeNcoeffs(j)-2;

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        //if this mesh id has not already been visited then we have
                        //a new block location. Otherwise we point the array
                        //location to the already assigned edge.
                        if(uniqueEdgeMap.count(meshEdgeId)==0)
                        {
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;

                            m_edgemodeoffset[edgematrixlocation]=nedgemodes*nedgemodes;

                            ntotaledgeentries+=nedgemodes*nedgemodes;

                            n_blks[1+edgematrixlocation++]=
                                locExpansion->GetEdgeNcoeffs(j)-2;

                        }

                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=(uniqueEdgeMap[meshEdgeId])*nedgemodes*nedgemodes+k;

                            m_localEdgeToGlobalMatrixMap[edgematrixoffset+k]
                                = (uniqueEdgeMap[meshEdgeId])*nedgemodes*nedgemodes+k;

                            m_EdgeBlockToUniversalMap[vGlobal]
                                = meshEdgeId * maxEdgeDof * maxEdgeDof + k + 1;
                        }
                        edgematrixoffset+=maxEdgeDof*maxEdgeDof;
                    }
                }
            }

            edgematrixoffset=0;

            BlkMat = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);
                nCoeffs=locExpansion->NumBndryCoeffs();

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
                    DNekScalMatSharedPtr tmp_mat;
                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nVerts,nVerts,zero,storage);

                    vMap1=locExpansion->GetVertexMap(v);

                    meshVertId = locExpansion->GetGeom2D()->GetVid(v);
                    
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
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);
                                cout<<globalrow<<endl;
                                vertValue = vertArray[globalrow]
                                      + sign1*sign2*S(vMap1,vMap2);

                                vertArray[globalrow] = vertValue;

                                m_VertBlockToUniversalMap[globalrow]
                                = meshVertId * nVerts * nVerts + 1;
                            }
                        }
                    }
                }
                
                //loop over edges of the element and return the edge map
                for (eid=0; eid<nEdges; ++eid)
                {
                    nedgemodes=locExpansion->GetEdgeNcoeffs(eid)-2;

                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nedgemodes,nedgemodes,zero,storage);
                    
                    meshEdgeId = locExpansion->GetGeom2D()->GetEid(eid);
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

                                NekDouble globalEdgeValue = sign1*sign2*S(eMap1,eMap2);

                                m_EdgeBlockArray[edgematrixoffset+v*nedgemodes+m]=globalEdgeValue;
                            }
                        }
                        edgematrixoffset+=maxEdgeDof*maxEdgeDof;
                    }
                }

                //offset for the expansion
                cnt+=offset;
            }

            //Assemble edge matrices of each process
            Array<OneD, NekDouble> m_GlobalEdgeBlock(ntotaledgeentries);
            Vmath::Zero(ntotaledgeentries, m_GlobalEdgeBlock.get(), 1);
            Vmath::Assmb(m_EdgeBlockArray.num_elements(), 
                         m_EdgeBlockArray.get(), 
                         m_localEdgeToGlobalMatrixMap.get(), 
                         m_GlobalEdgeBlock.get());

            //Exchange edge data over different processes
            if(nNonDirVerts != 0)
            {
                Gs::gs_data *tmp = Gs::Init(m_VertBlockToUniversalMap, m_comm);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
            }

            //Exchange edge data over different processes
            if(nNonDirEdgeIDs != 0)
            {
                Gs::gs_data *tmp1 = Gs::Init(m_EdgeBlockToUniversalMap, m_comm);
                Gs::Gather(m_GlobalEdgeBlock, Gs::gs_add, tmp1);
            }

            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }
            
            //Set the first block to be the diagonal of the vertex space
            BlkMat->SetBlock(0,0, VertBlk);

            offset=0;
            //Build the edge matrices from the vector
            for(int loc=0; loc<nNonDirEdgeIDs; ++loc)
            {
                DNekMatSharedPtr m_gmat = 
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nedgemodes,nedgemodes,zero,storage);

                for (v=0; v<nedgemodes; ++v)
                {
                    for (m=0; m<nedgemodes; ++m)
                    {
                        NekDouble EdgeValue = m_GlobalEdgeBlock[offset+v*nedgemodes+m];
                        m_gmat->SetValue(v,m,EdgeValue);
                    }
                }
    
                BlkMat->SetBlock(1+loc,1+loc, m_gmat);

                offset+=m_edgemodeoffset[loc];
            }

            int totblks=BlkMat->GetNumberOfBlockRows();

            for (i=1; i< totblks; ++i)
            {
                unsigned int nmodes=BlkMat->GetNumberOfRowsInBlockRow(i);
                DNekMatSharedPtr tmp_mat = 
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                tmp_mat=BlkMat->GetBlock(i,i);
                for (j=0; j<tmp_mat->GetRows(); ++j)
                {
                    for (k=0; k<tmp_mat->GetRows(); ++k)
                    {
                        cout<<(*tmp_mat)(j,k)<<" ";
                    }
                    cout<<endl;
                }
                cout<<endl;

                tmp_mat->Invert();
                BlkMat->SetBlock(i,i,tmp_mat);
            }
        }

        void PreconditionerBlock::BlockPreconditioner3D()
        {
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            StdRegions::StdExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;

            int i, j, k, nel;
            int nVerts, nEdges,nFaces; 
            int eid, fid, n, cnt, nedgemodes, nfacemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol, nCoeffs;
            NekDouble vertValue;

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;
            
            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> m_VertBlockToUniversalMap(nNonDirVerts,-1);

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
            Array<OneD, unsigned int> m_edgemodeoffset(nNonDirEdgeIDs);
            Array<OneD, unsigned int> m_facemodeoffset(nNonDirFaceIDs);

            Array<OneD, unsigned int> m_edgeglobaloffset(nNonDirEdgeIDs);
            Array<OneD, unsigned int> m_faceglobaloffset(nNonDirFaceIDs);

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>& 
                bndConditions = expList->GetBndConditions();

            int meshVertId;
            int meshEdgeId;
            int meshFaceId;

            //Determine which boundary edges and faces have dirichlet values
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp = boost::dynamic_pointer_cast<
                        StdRegions::StdExpansion2D>(
                            bndCondExp[i]->GetExp(j));
                    if (bndConditions[i]->GetBoundaryConditionType() == 
                        SpatialDomains::eDirichlet)
                    {
                        for(k = 0; k < bndCondFaceExp->GetNedges(); k++)
                        {
                            meshEdgeId = (bndCondFaceExp->GetGeom2D())->GetEid(k);
                            if(edgeDirMap.count(meshEdgeId) == 0)
                            {
                                edgeDirMap[meshEdgeId] = 1;
                            }
                        }
                        meshFaceId = (bndCondFaceExp->GetGeom2D())->GetFid();
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
                    meshEdgeId = locExpansion->GetGeom3D()->GetEid(j);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        if(uniqueEdgeMap.count(meshEdgeId)==0)
                        {
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;

                            m_edgeglobaloffset[edgematrixlocation]+=ntotaledgeentries;

                            m_edgemodeoffset[edgematrixlocation]=dof*dof;

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

                    meshFaceId = locExpansion->GetGeom3D()->GetFid(j);

                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        if(uniqueFaceMap.count(meshFaceId)==0)
                        {
                            uniqueFaceMap[meshFaceId]=facematrixlocation;

                            m_facemodeoffset[facematrixlocation]=dof*dof;
                            
                            m_faceglobaloffset[facematrixlocation]+=ntotalfaceentries;

                            ntotalfaceentries+=dof*dof;
                            
                            n_blks[1+nNonDirEdgeIDs+facematrixlocation++]=dof;
                            
                        }
                        nlocalNonDirFaces+=dof*dof;
                    }

                }
            }

            m_comm = expList->GetComm();
            //m_comm = expList->GetComm()->GetRowComm(); 2D
            m_comm->AllReduce(maxEdgeDof, LibUtilities::ReduceMax);
            m_comm->AllReduce(maxFaceDof, LibUtilities::ReduceMax);

            //Allocate arrays for block to universal map (number of expansions * p^2)
            Array<OneD, long> m_EdgeBlockToUniversalMap(ntotaledgeentries,-1);
            Array<OneD, long> m_FaceBlockToUniversalMap(ntotalfaceentries,-1);

            Array<OneD, int> m_localEdgeToGlobalMatrixMap(nlocalNonDirEdges,-1);
            Array<OneD, int> m_localFaceToGlobalMatrixMap(nlocalNonDirFaces,-1);

            //Allocate arrays to store matrices (number of expansions * p^2)
            Array<OneD, NekDouble> m_EdgeBlockArray(nlocalNonDirEdges,-1);
            Array<OneD, NekDouble> m_FaceBlockArray(nlocalNonDirFaces,-1);

            int ecnt;
            int fcnt;
            int edgematrixoffset=0;
            int facematrixoffset=0;
            int vGlobal;
            int nbndCoeffs=0;

            for(ecnt=fcnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->GetGeom3D()->GetEid(j);

                    nedgemodes=locExpansion->GetEdgeNcoeffs(j)-2;

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=m_edgeglobaloffset[uniqueEdgeMap[meshEdgeId]]+k;

                            m_localEdgeToGlobalMatrixMap[edgematrixoffset+k]=vGlobal;

                            m_EdgeBlockToUniversalMap[vGlobal]
                                = meshEdgeId * maxEdgeDof * maxEdgeDof + k + 1;
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }

                //loop over the faces of the expansion
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    //get mesh face id
                    meshFaceId = locExpansion->GetGeom3D()->GetFid(j);

                    nfacemodes = locExpansion->GetFaceIntNcoeffs(j);

                    //Check if face is has dirichlet values
                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        for(k=0; k<nfacemodes*nfacemodes; ++k)
                        {
                            vGlobal=m_faceglobaloffset[uniqueFaceMap[meshFaceId]]+k;
                            
                            m_localFaceToGlobalMatrixMap[facematrixoffset+k]
                                = vGlobal;
                            
                            m_FaceBlockToUniversalMap[vGlobal]
                                = meshFaceId * maxFaceDof * maxFaceDof + k + 1;
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                    }
                }
                nbndCoeffs=+locExpansion->NumBndryCoeffs();
            }

            edgematrixoffset=0;
            facematrixoffset=0;

            BlkMat = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);
                nCoeffs=locExpansion->NumBndryCoeffs();

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
                    DNekScalMatSharedPtr tmp_mat;
                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nVerts,nVerts,zero,storage);

                    vMap1=locExpansion->GetVertexMap(v);

                    meshVertId = locExpansion->GetGeom3D()->GetVid(v);
                    
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
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);
                                
                                vertValue = vertArray[globalrow]
                                      + sign1*sign2*S(vMap1,vMap2);

                                vertArray[globalrow] = vertValue;

                                m_VertBlockToUniversalMap[globalrow]
                                = meshVertId * nVerts * nVerts + 1;
                            }
                        }
                    }
                }
                
                //loop over edges of the element and return the edge map
                for (eid=0; eid<nEdges; ++eid)
                {
                    nedgemodes=locExpansion->GetEdgeNcoeffs(eid)-2;

                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nedgemodes,nedgemodes,zero,storage);
                    
                    meshEdgeId = locExpansion->GetGeom3D()->GetEid(eid);
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

                                NekDouble globalEdgeValue = sign1*sign2*S(eMap1,eMap2);

                                m_EdgeBlockArray[edgematrixoffset+v*nedgemodes+m]=globalEdgeValue;
                            }
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }
                
                //loop over faces of the element and return the face map
                for (fid=0; fid<nFaces; ++fid)
                {
                    nfacemodes = locExpansion->GetFaceIntNcoeffs(fid);

                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nfacemodes,nfacemodes,zero,storage);

                    meshFaceId = locExpansion->GetGeom3D()->GetFid(fid);
                    
                    Array<OneD, unsigned int> facemodearray = 
                        locExpansion->GetFaceInverseBoundaryMap(fid);

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

                                //test here with local to global map
                                m_FaceBlockArray[facematrixoffset+v*nfacemodes+m]=globalFaceValue;
                            }
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                    }
                }

                //offset for the expansion
                cnt+=offset;
            }

            //Assemble edge matrices of each process
            Array<OneD, NekDouble> m_GlobalEdgeBlock(ntotaledgeentries);
            Vmath::Zero(ntotaledgeentries, m_GlobalEdgeBlock.get(), 1);
            Vmath::Assmb(m_EdgeBlockArray.num_elements(), 
                         m_EdgeBlockArray.get(), 
                         m_localEdgeToGlobalMatrixMap.get(), 
                         m_GlobalEdgeBlock.get());

            //Assemble face matrices of each process
            Array<OneD, NekDouble> m_GlobalFaceBlock(ntotalfaceentries);
            Vmath::Zero(ntotalfaceentries, m_GlobalFaceBlock.get(), 1);
            Vmath::Assmb(m_FaceBlockArray.num_elements(), 
                         m_FaceBlockArray.get(), 
                         m_localFaceToGlobalMatrixMap.get(), 
                         m_GlobalFaceBlock.get());

            //Exchange edge data over different processes
            if(nNonDirVerts != 0)
            {
                Gs::gs_data *tmp = Gs::Init(m_VertBlockToUniversalMap, m_comm);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
            }

            //Exchange edge data over different processes
            if(nNonDirEdgeIDs != 0)
            {
                Gs::gs_data *tmp1 = Gs::Init(m_EdgeBlockToUniversalMap, m_comm);
                Gs::Gather(m_GlobalEdgeBlock, Gs::gs_add, tmp1);
            }

            //Exchange face data over different processes
            if(nNonDirFaceIDs != 0)
            {
                Gs::gs_data *tmp2 = Gs::Init(m_FaceBlockToUniversalMap, m_comm);
                Gs::Gather(m_GlobalFaceBlock, Gs::gs_add, tmp2);
            }


            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }
            
            //Set the first block to be the diagonal of the vertex space
            BlkMat->SetBlock(0,0, VertBlk);

            offset=0;
            //Build the edge matrices from the vector
            for(int loc=0; loc<nNonDirEdgeIDs; ++loc)
            {
                nedgemodes=n_blks[1+loc];

                DNekMatSharedPtr m_gmat = 
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nedgemodes,nedgemodes,zero,storage);

                for (v=0; v<nedgemodes; ++v)
                {
                    for (m=0; m<nedgemodes; ++m)
                    {
                        NekDouble EdgeValue = m_GlobalEdgeBlock[offset+v*nedgemodes+m];
                        m_gmat->SetValue(v,m,EdgeValue);
                    }
                }
    
                BlkMat->SetBlock(1+loc,1+loc, m_gmat);

                offset+=m_edgemodeoffset[loc];
            }

            offset=0;
            //Build the face matrices from the vector
            for(int loc=0; loc<nNonDirFaceIDs; ++loc)
            {
                nfacemodes=n_blks[1+nNonDirEdgeIDs+loc];

                DNekMatSharedPtr m_gmat = 
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nfacemodes,nfacemodes,zero,storage);

                for (v=0; v<nfacemodes; ++v)
                {
                    for (m=0; m<nfacemodes; ++m)
                    {
                        NekDouble FaceValue = m_GlobalFaceBlock[offset+v*nfacemodes+m];
                        m_gmat->SetValue(v,m,FaceValue);
                    }
                }

                BlkMat->SetBlock(1+nNonDirEdgeIDs+loc,1+nNonDirEdgeIDs+loc, m_gmat);

                offset+=m_facemodeoffset[loc];
            }

            
            int totblks=BlkMat->GetNumberOfBlockRows();

            for (i=1; i< totblks; ++i)
            {
                unsigned int nmodes=BlkMat->GetNumberOfRowsInBlockRow(i);
                DNekMatSharedPtr tmp_mat = 
                    MemoryManager<DNekMat>::AllocateSharedPtr
                    (nmodes,nmodes,zero,storage);

                tmp_mat=BlkMat->GetBlock(i,i);
                for (j=0; j<tmp_mat->GetRows(); ++j)
                {
                    for (k=0; k<tmp_mat->GetRows(); ++k)
                    {
                        //cout<<(*tmp_mat)(j,k)<<" ";
                    }
                    //cout<<endl;
                }

                tmp_mat->Invert();
                BlkMat->SetBlock(i,i,tmp_mat);
            }
            cout<<endl;
        }

  
        /**
         *
         */
        void PreconditionerBlock::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(solvertype)
            {
                case MultiRegions::eIterativeFull:
                {
                    int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                    int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                    int nNonDir = nGlobal-nDir;
                    DNekScalBlkMat &M = (*GloBlkMat);
                    
                    NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                    NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                    z = M * r;
                }
                break;
                case MultiRegions::eIterativeStaticCond:
                {
                    int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                    int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                    int nNonDir = nGlobal-nDir;
                    DNekBlkMat &M = (*BlkMat);
                    
                    NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                    NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                    z = M * r;
                }
                break;
                default:
                    ASSERTL0(0,"Unknown preconditioner");
                    break;
	    }
        }
    }
}






