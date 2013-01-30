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
           m_locToGloMap(pLocToGloMap),
           m_preconType(pLocToGloMap->GetPreconType())
         {
	 }

        void PreconditionerBlock::v_InitObject()
        {
            BlockPreconditioner2D();
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
            const StdRegions::StdExpansionVector &locExpVector = 
                *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;

            int nRow, i, j, k, nmodes, ntotaledgemodes, ntotalfacemodes, nel;
            int nVerts, nEdges,nFaces; 
            int eid, fid, eid2, fid2, n, cnt, nedgemodes, nfacemodes;
            int nEdgeCoeffs, nFaceCoeffs;
            NekDouble zero = 0.0;
            NekDouble MatrixValue;

            int vMap1, vMap2, sign1, sign2, gid1, gid2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol, bnd_rows, nCoeffs;
            NekDouble globalMatrixValue, globalRValue;

            //matrix storage
            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            MatrixStorage blkmatStorage = eDIAGONAL;

            //local element static condensed matrices
            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr    m_RS;
            DNekMatSharedPtr    m_RSRT;

            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDir = nGlobal-nDirBnd;
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();
            int nNonDirEdges  = m_locToGloMap->GetNumNonDirEdgeModes();
            int nNonDirFaces  = m_locToGloMap->GetNumNonDirFaceModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            int n_exp = expList->GetNumElmts();

            map<int,int> edgeDirMap;
            map<int,int> faceDirMap;

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>& bndConditions = expList->GetBndConditions();

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


            int nDirFaceIDs  = m_locToGloMap->GetNumDirFaces();
            int nDirEdgeIDs  = m_locToGloMap->GetNumDirEdges();
            int nNonDirFaceIDs  = m_locToGloMap->GetNumNonDirFaces();
            int nNonDirEdgeIDs  = m_locToGloMap->GetNumNonDirEdges();

            int dof;
            int maxDof;

            // Loop over all the elements in the domain and compute max edge
            // DOF. Reduce across all processes to get universal maximum.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);

                for (j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    dof    = locExpansion->GetFaceNcoeffs(j);
                    maxDof = (dof > maxDof ? dof : maxDof);
                }
            }

            //Set up mappings
            map<int,int> uniqueEdgeMap;
            map<int,int> uniqueFaceMap;

            //this should be of size total number of local edges
            Array<OneD, unsigned int> m_edgeIDToMatrixBlockLocation(nNonDirEdgeIDs);
            Array<OneD, unsigned int> m_faceIDToMatrixBlockLocation(nNonDirFaceIDs);

            //set the number of blocks in the matrix
            Array<OneD,unsigned int> n_blks(1+nNonDirEdgeIDs+nNonDirFaceIDs);
            n_blks[0]=nNonDirVerts;

            int edgematrixlocation=0;
            int facematrixlocation=0;
            int ecnt;
            int fcnt;

            for(ecnt=fcnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);

                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->GetGeom3D()->GetEid(j);

                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        //if this mesh id has not already been visited then we have
                        //a new block location. Otherwise we point the array
                        //location to the already assigned edge.
                        if(uniqueEdgeMap.count(meshEdgeId)==0)
                        {
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;
                            n_blks[1+edgematrixlocation++]=
                                locExpansion->GetEdgeNcoeffs(j)-2;
                        }
                    }
                }

                //loop over the faces of the expansion
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    //get mesh face id
                    meshFaceId = locExpansion->GetGeom3D()->GetFid(j);

                    //Check if face is has dirichlet values
                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        //if this mesh id has not already been visited then we have
                        //a new block location. Otherwise we point the array
                        //location to the already assigned edge.
                        if(uniqueFaceMap.count(meshFaceId)==0)
                        {
                            uniqueFaceMap[meshFaceId]=
                                nNonDirEdgeIDs+facematrixlocation;
                            n_blks[1+nNonDirEdgeIDs+facematrixlocation++]=
                                locExpansion->GetFaceIntNcoeffs(j);
                        }
                    }
                }
            }

            BlkMat = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);

            const Array<OneD,const unsigned int>& nbdry_size
                    = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);
                nCoeffs=locExpansion->NumBndryCoeffs();
                StdRegions::ExpansionType eType=
                    locExpansion->DetExpansionType();

                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);
		
                //Extract boundary block (elemental S1)
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                offset = bnd_mat->GetRows();

                DNekScalMat &S=(*bnd_mat);

                nVerts=locExpansion->GetGeom()->GetNumVerts();
                nEdges=locExpansion->GetGeom()->GetNumEdges();
                nFaces=locExpansion->GetGeom()->GetNumFaces();

                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {

                    DNekScalMatSharedPtr tmp_mat;
                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nVerts,nVerts,zero,storage);

                    //Get vertex map
                    vMap1 = locExpansion->GetVertexMap(v);
                    
                    //Get vertex map
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;
                    
                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            //Get vertex map
                            vMap1 = locExpansion->GetVertexMap(m);

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
                                
                                //Global matrix value
                                globalMatrixValue = VertBlk->
                                    GetValue(globalrow,globalcol)
                                    + sign1*sign2*S(vMap1,vMap2);

                                //build matrix containing the linear finite
                                //element space
                                VertBlk->SetValue
                                    (globalrow,globalcol,globalMatrixValue);
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

                                m_locMat->SetValue(v,m,globalEdgeValue);
                            }
                        }
                        
                        DNekMatSharedPtr tmp_mat = 
                            MemoryManager<DNekMat>::AllocateSharedPtr
                            (nedgemodes,nedgemodes,zero,storage);
                        int loc = uniqueEdgeMap[meshEdgeId]+1;

                        //Get the current matrix in this location
                        tmp_mat=BlkMat->GetBlock(loc,loc);

                        //if is already a matrix then add this matrix with the
                        //new one
                        if(tmp_mat != NullDNekMatSharedPtr)
                        {
                            (*m_locMat)=(*tmp_mat)+(*m_locMat);
                        }

                        //Set the matrix block at the correct location
                        BlkMat->SetBlock(loc,loc, m_locMat); 

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
                                    
                                //Set this value in the local matrix (which
                                //will later form a block of the preconditioner)
                                m_locMat->SetValue(v,m,globalFaceValue);
                            }
                        }

                        DNekMatSharedPtr tmp_mat = 
                            MemoryManager<DNekMat>::AllocateSharedPtr
                            (nfacemodes,nfacemodes,zero,storage);
                        int loc = uniqueFaceMap[meshFaceId]+1;

                        tmp_mat=BlkMat->GetBlock(loc,loc);
                        if(tmp_mat != NullDNekMatSharedPtr)
                        {
                            (*m_locMat)=(*tmp_mat)+(*m_locMat);
                        }

                        BlkMat->SetBlock(loc,loc, m_locMat);
                    }
                }

                //offset for the expansion
                cnt+=offset;

            }

            //Set the first block to be the diagonal of the vertex space
            BlkMat->SetBlock(0,0, VertBlk);
            
            int totblks=BlkMat->GetNumberOfBlockRows();

            for (i=0; i< totblks; ++i)
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

                tmp_mat->Invert();
                BlkMat->SetBlock(i,i,tmp_mat);
            }
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






