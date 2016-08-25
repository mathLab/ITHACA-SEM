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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerLowEnergy.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalLinSys.h>
#include <LocalRegions/MatrixKey.h>
#include <math.h>

using namespace std;

namespace Nektar
{
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
            const boost::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap),
              m_linsys(plinsys),
              m_locToGloMap(pLocToGloMap)
        {
        }
        
        void PreconditionerLowEnergy::v_InitObject()
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            ASSERTL0(solvertype == eIterativeStaticCond ||
                     solvertype == ePETScStaticCond, "Solver type not valid");

            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            
            StdRegions::StdExpansionSharedPtr locExpansion;

            locExpansion = expList->GetExp(0);
            
            int nDim = locExpansion->GetShapeDimension();
            
            ASSERTL0(nDim==3,
                     "Preconditioner type only valid in 3D");
            
            //Sets up reference element and builds transformation matrix for
            // maximum polynomial order meshes
            SetUpReferenceElements();

            //Set up block transformation matrix
            SetupBlockTransformationMatrix();

            //Sets up multiplicity map for transformation from global to local
            CreateMultiplicityMap();
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
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            LocalRegions::ExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            int i, j, k;
            int nVerts, nEdges,nFaces; 
            int eid, fid, n, cnt, nedgemodes, nfacemodes;
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
            
            int nDirBnd      = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts = m_locToGloMap->GetNumNonDirVertexModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);
            
            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> VertBlockToUniversalMap(nNonDirVerts,-1);

            //maps for different element types
            int n_exp = expList->GetNumElmts();
            int nNonDirEdgeIDs=m_locToGloMap->GetNumNonDirEdges();
            int nNonDirFaceIDs=m_locToGloMap->GetNumNonDirFaces();
            
            //set the number of blocks in the matrix
            int numBlks = 1+nNonDirEdgeIDs+nNonDirFaceIDs;
            Array<OneD,unsigned int> n_blks(numBlks);
            for(i = 0; i < numBlks; ++i)
            {
                n_blks[i] = 0;
            }
            n_blks[0]=nNonDirVerts;

            set<int> edgeDirMap;  
            set<int> faceDirMap;  
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
                edgeDirMap.insert(meshEdgeId);
            }
            
            //Determine which boundary edges and faces have dirichlet values
            for(i = 0; i < bndCondExp.num_elements(); i++)
            {
                cnt = 0;
                for(j = 0; j < bndCondExp[i]->GetNumElmts(); j++)
                {
                    bndCondFaceExp = boost::dynamic_pointer_cast<
                    StdRegions::StdExpansion2D>(bndCondExp[i]->GetExp(j));
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
                        meshFaceId = bndCondFaceExp->as<LocalRegions::Expansion2D>()->GetGeom2D()->GetFid();
                        faceDirMap.insert(meshFaceId);
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

            map<int,int> EdgeSize;
            map<int,int> FaceSize;
            
            /// -  Count  edges, face and add up edges and face sizes
            for(n = 0; n < n_exp; ++n)
            {
                eid = expList->GetOffset_Elmt_Id(n);
                locExpansion = expList->GetExp(eid);

                nEdges = locExpansion->GetNedges();
                for(j = 0; j < nEdges; ++j)
                {
                    int nEdgeInteriorCoeffs = locExpansion->GetEdgeNcoeffs(j) - 2;
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);
                    EdgeSize[meshEdgeId] = nEdgeInteriorCoeffs;
                }
                
                nFaces = locExpansion->GetNfaces();
                for(j = 0; j < nFaces; ++j)
                {
                    int nFaceInteriorCoeffs = locExpansion->GetFaceIntNcoeffs(j);
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(j);
                    FaceSize[meshFaceId] = nFaceInteriorCoeffs;
                }
            }

            m_comm = expList->GetComm();

            // Loop over all the elements in the domain and compute max edge
            // DOF and set up unique ordering. 

            // First do periodic edges 
            PeriodicMap::const_iterator pIt;
            for (pIt = periodicEdges.begin(); pIt != periodicEdges.end(); ++pIt)
            {
                meshEdgeId = pIt->first;

                if(edgeDirMap.count(meshEdgeId)==0)
                {
                    dof = EdgeSize[meshEdgeId]; 
                    if(uniqueEdgeMap.count(meshEdgeId)==0 && dof > 0)
                    {
                        bool SetUpNewEdge = true;
                        
                        
                        for (i = 0; i < pIt->second.size(); ++i)
                        {
                            if (!pIt->second[i].isLocal)
                            {
                                continue;
                            }
                            
                            int meshEdgeId2 = pIt->second[i].id;
                            
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
                            uniqueEdgeMap[meshEdgeId]=edgematrixlocation;
                            
                            edgeglobaloffset[edgematrixlocation]+=ntotaledgeentries;
                            
                            edgemodeoffset[edgematrixlocation]=dof*dof;

                            ntotaledgeentries+=dof*dof;
                            
                            n_blks[1+edgematrixlocation++]=dof;   
                            
                        }
                    }
                }
            }

            
            for(cnt=n=0; n < n_exp; ++n)
            {
                eid = expList->GetOffset_Elmt_Id(n);
                locExpansion = expList->GetExp(eid);

                for (j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);
                    dof    = EdgeSize[meshEdgeId];
                    maxEdgeDof = (dof > maxEdgeDof ? dof : maxEdgeDof);

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
            // - Periodic faces
            for (pIt = periodicFaces.begin(); pIt != periodicFaces.end(); ++pIt)
            {
                meshFaceId = pIt->first;
                
                if(faceDirMap.count(meshFaceId)==0)
                {
                    dof = FaceSize[meshFaceId];
                    
                    if(uniqueFaceMap.count(meshFaceId) == 0 && dof > 0)
                    {
                        bool SetUpNewFace = true;
                        
                        if(pIt->second[0].isLocal)
                        {
                            int meshFaceId2 = pIt->second[0].id;
                            
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
                            uniqueFaceMap[meshFaceId]=facematrixlocation;
                            
                            facemodeoffset[facematrixlocation]=dof*dof;
                            
                            faceglobaloffset[facematrixlocation]+=ntotalfaceentries;
                            
                            ntotalfaceentries+=dof*dof;
                            
                            n_blks[1+nNonDirEdgeIDs+facematrixlocation++]=dof;
                        }
                    }
                }
            }


            for(cnt=n=0; n < n_exp; ++n)
            {
                eid = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(eid);

                for (j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(j);

                    dof        = FaceSize[meshFaceId];
                    maxFaceDof = (dof > maxFaceDof ? dof : maxFaceDof);
                 
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
                eid = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(eid);
                
                //loop over the edges of the expansion
                for(j = 0; j < locExpansion->GetNedges(); ++j)
                {
                    //get mesh edge id
                    meshEdgeId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetEid(j);
                    
                    nedgemodes=locExpansion->GetEdgeNcoeffs(j)-2;
                    
                    if(edgeDirMap.count(meshEdgeId)==0)
                    {
                        // Determine the Global edge offset
                        int edgeOffset = edgeglobaloffset[uniqueEdgeMap[meshEdgeId]];
                        
                        // Determine a universal map offset 
                        int uniOffset = meshEdgeId;                            
                        pIt = periodicEdges.find(meshEdgeId);
                        if (pIt != periodicEdges.end())
                        {
                            for (int l = 0; l < pIt->second.size(); ++l)
                            {
                                uniOffset = min(uniOffset, pIt->second[l].id);
                            }
                        }
                        uniOffset = uniOffset *maxEdgeDof*maxEdgeDof; 
                        
                        for(k=0; k<nedgemodes*nedgemodes; ++k)
                        {
                            vGlobal=edgeOffset+k;
                            localEdgeToGlobalMatrixMap[edgematrixoffset+k]=vGlobal;
                            EdgeBlockToUniversalMap[vGlobal] = uniOffset + k + 1;
                        }
                        edgematrixoffset+=nedgemodes*nedgemodes;
                    }
                }
                
                Array<OneD, unsigned int>           faceInteriorMap;
                Array<OneD, int>                    faceInteriorSign;
                //loop over the faces of the expansion
                for(j = 0; j < locExpansion->GetNfaces(); ++j)
                {
                    //get mesh face id
                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(j);

                    nfacemodes = locExpansion->GetFaceIntNcoeffs(j);

                    //Check if face has dirichlet values
                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        // Determine the Global edge offset
                        int faceOffset = faceglobaloffset[uniqueFaceMap[meshFaceId]];                      
                        // Determine a universal map offset 
                        int uniOffset = meshFaceId;                            
                        // use minimum face edge when periodic 
                        pIt = periodicFaces.find(meshFaceId);
                        if (pIt != periodicFaces.end())
                        {
                            uniOffset = min(uniOffset, pIt->second[0].id);
                        }
                        uniOffset = uniOffset * maxFaceDof * maxFaceDof; 
                        
                        for(k=0; k<nfacemodes*nfacemodes; ++k)
                        {
                            vGlobal=faceOffset+k;
                            
                            localFaceToGlobalMatrixMap[facematrixoffset+k]
                                = vGlobal;
                            
                            FaceBlockToUniversalMap[vGlobal] = uniOffset + k + 1;
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                    }
                }
            }
                
            edgematrixoffset=0;
            facematrixoffset=0;

            m_BlkMat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(n_blks, n_blks, blkmatStorage);
            
            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                eid = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(eid);
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

                //Calculate R*S*trans(R)
                RSRT = R*S*Transpose(R);

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
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);
                                
                                vertArray[globalrow]
                                    += sign1*sign2*RSRT(vMap1,vMap2);


                                meshVertId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetVid(v);
                            
                                pIt = periodicVerts.find(meshVertId);
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
                    nedgemodes=locExpansion->GetEdgeNcoeffs(eid)-2;

                    DNekMatSharedPtr m_locMat = 
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

                                NekDouble globalEdgeValue = sign1*sign2*RSRT(eMap1,eMap2);

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

                    DNekMatSharedPtr m_locMat = 
                        MemoryManager<DNekMat>::AllocateSharedPtr
                        (nfacemodes,nfacemodes,zero,storage);

                    meshFaceId = locExpansion->as<LocalRegions::Expansion3D>()->GetGeom3D()->GetFid(fid);
                    
                    if(faceDirMap.count(meshFaceId)==0)
                    {
                        Array<OneD, unsigned int> facemodearray;
                        StdRegions::Orientation faceOrient = locExpansion->GetForient(fid);
                        
                        pIt = periodicFaces.find(meshFaceId);
                        if (pIt != periodicFaces.end())
                        {
                            if(meshFaceId == min(meshFaceId, pIt->second[0].id))
                            {
                                facemodearray = locExpansion->GetFaceInverseBoundaryMap(fid,faceOrient);
                                faceOrient = DeterminePeriodicFaceOrient(faceOrient,pIt->second[0].orient);
                            }
                        }
                        
                        facemodearray = locExpansion->GetFaceInverseBoundaryMap(fid,faceOrient);
                        
                        
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

                                // Get the face-face value from the
                                // low energy matrix (S2)
                                NekDouble globalFaceValue = sign1*sign2*RSRT(fMap1,fMap2);

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

            bool verbose =
                expList->GetSession()->DefinesCmdLineArgument("verbose");

            if(nNonDirVerts!=0)
            {
                //Exchange vertex data over different processes
                Gs::gs_data *tmp = Gs::Init(VertBlockToUniversalMap, m_comm, verbose);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
                
            }
            
            Array<OneD, NekDouble> GlobalEdgeBlock(ntotaledgeentries,0.0);
            if(ntotaledgeentries)
            {
                //Assemble edge matrices of each process
                Vmath::Assmb(EdgeBlockArray.num_elements(),  
                             EdgeBlockArray, 
                             localEdgeToGlobalMatrixMap, 
                             GlobalEdgeBlock);
            }

            //Exchange edge data over different processes
            Gs::gs_data *tmp1 = Gs::Init(EdgeBlockToUniversalMap, m_comm, verbose);
            Gs::Gather(GlobalEdgeBlock, Gs::gs_add, tmp1);

            Array<OneD, NekDouble> GlobalFaceBlock(ntotalfaceentries,0.0);
            if(ntotalfaceentries)
            {
                //Assemble face matrices of each process
                Vmath::Assmb(FaceBlockArray.num_elements(),
                             FaceBlockArray, 
                             localFaceToGlobalMatrixMap, 
                             GlobalFaceBlock);
            }

            //Exchange face data over different processes
            Gs::gs_data *tmp2 = Gs::Init(FaceBlockToUniversalMap, m_comm, verbose);
            Gs::Gather(GlobalFaceBlock, Gs::gs_add, tmp2);
            
            // Populate vertex block
            for (int i = 0; i < nNonDirVerts; ++i)
            {
                VertBlk->SetValue(i,i,1.0/vertArray[i]);
            }

            //Set the first block to be the diagonal of the vertex space
            m_BlkMat->SetBlock(0,0, VertBlk);
            
            offset=0;
            //Build the edge matrices from the vector
            DNekMatSharedPtr gmat;
            for(int loc=0; loc<nNonDirEdgeIDs; ++loc)
            {
                nedgemodes = n_blks[1+loc];
                gmat = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nedgemodes,nedgemodes,zero,storage);
                
                for (v=0; v<nedgemodes; ++v)
                {
                    for (m=0; m<nedgemodes; ++m)
                    {
                        NekDouble EdgeValue = GlobalEdgeBlock[offset+v*nedgemodes+m];
                        gmat->SetValue(v,m,EdgeValue);

                    }
                }

                m_BlkMat->SetBlock(1+loc,1+loc, gmat);
                offset+=edgemodeoffset[loc];
            }
            
            offset=0;
            
            Array<OneD, int> globalToUniversalMap = m_locToGloMap->GetGlobalToUniversalBndMap();
            //Build the face matrices from the vector
            for(int loc=0; loc<nNonDirFaceIDs; ++loc)
            {
                nfacemodes=n_blks[1+nNonDirEdgeIDs+loc];
                gmat = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nfacemodes,nfacemodes,zero,storage);
                
                for (v=0; v<nfacemodes; ++v)
                {
                    for (m=0; m<nfacemodes; ++m)
                    {
                        NekDouble FaceValue = GlobalFaceBlock[offset+v*nfacemodes+m];
                        gmat->SetValue(v,m,FaceValue);
                        
                    }
                }
                m_BlkMat->SetBlock(1+nNonDirEdgeIDs+loc,1+nNonDirEdgeIDs+loc, gmat);
                offset+=facemodeoffset[loc];
            }
               
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
            int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
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
       void PreconditionerLowEnergy::SetupBlockTransformationMatrix()
       {
           boost::shared_ptr<MultiRegions::ExpList> 
               expList=((m_linsys.lock())->GetLocMat()).lock();
           StdRegions::StdExpansionSharedPtr locExpansion;
           StdRegions::StdExpansionSharedPtr locExpansionSav;

           int n;
 
           const Array<OneD,const unsigned int>& nbdry_size
               = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

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
               int eid = expList->GetOffset_Elmt_Id(n);
               locExpansion = expList->GetExp(eid);

               int nbndcoeffs = locExpansion->NumBndryCoeffs();

               if(m_sameBlock.size() == 0)
               {
                   rmat = ExtractLocMat(locExpansion);
                   //Block R matrix
                   m_RBlk->SetBlock(n, n, rmat);

                   invrmat = MemoryManager<DNekMat>::AllocateSharedPtr(*rmat);
                   invrmat->Invert();

                   //Block inverse R matrix
                   m_InvRBlk->SetBlock(n, n, invrmat);
                   
                   m_sameBlock.push_back(pair<int,int>(1,nbndcoeffs));
                   locExpansionSav = locExpansion; 
               }
               else
               {
                   bool reuse = true; 
                   
                   // check to see if same as previous matrix and
                   // reuse if we can
                   for(int i = 0; i < 3; ++i)
                   {
                       if(locExpansionSav->GetBasis(i) != locExpansion->GetBasis(i))
                       {
                           reuse = false; 
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
                       rmat = ExtractLocMat(locExpansion);
                       //Block R matrix
                       m_RBlk->SetBlock(n, n, rmat);

                       invrmat = MemoryManager<DNekMat>::AllocateSharedPtr(*rmat);
                       invrmat->Invert();
                       //Block inverse R matrix
                       m_InvRBlk->SetBlock(n, n, invrmat);
                       
                       m_sameBlock.push_back(pair<int,int>(1,nbndcoeffs));
                       offset++;
                       locExpansionSav = locExpansion; 
                       
                   }
               }
           }
       }
        


        /**
         * \brief Transform the solution vector vector to low energy.
         *
         * As the conjugate gradient system is solved for the low energy basis,
         * the solution vector \f$\mathbf{x}\f$ must be transformed to the low
         * energy basis i.e. \f$\overline{\mathbf{x}}=\mathbf{R}\mathbf{x}\f$.
         */
        void PreconditionerLowEnergy::v_DoTransformToLowEnergy(
            Array<OneD, NekDouble>& pInOut,
            int offset)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            //Non-dirichlet boundary dofs
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,pInOut+offset,
                                          eWrapper);

            //Block transformation matrix
            DNekBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, 0.0);

            //Not actually needed but we should only work with the Global boundary dofs
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobBndDofs, pInOut.get(), 1, tmp.get(), 1);

            //Global boundary (with dirichlet values) to local boundary with multiplicity
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(),
                         tmp.get(), m_map.get(), pLocalIn.get());

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
			    0.0, pLocal.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }

            //Assemble local boundary to global non-dirichlet Dofs
            m_locToGloMap->AssembleBnd(F_LocBnd,F_HomBnd, nDirBndDofs);
        }

        /**
         * \brief Transform the solution vector to low energy form.
         *
         * As the conjugate gradient system is solved for the low energy basis,
         * the solution vector \f$\mathbf{x}\f$ must be transformed to the low
         * energy basis i.e. \f$\overline{\mathbf{x}}=\mathbf{R}\mathbf{x}\f$.
         */
        void PreconditionerLowEnergy::v_DoTransformToLowEnergy(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            //Input/output vectors should be length nGlobHomBndDofs
            ASSERTL1(pInput.num_elements() >= nGlobHomBndDofs,
                     "Input array is greater than the nGlobHomBndDofs");
            ASSERTL1(pOutput.num_elements() >= nGlobHomBndDofs,
                     "Output array is greater than the nGlobHomBndDofs");

            //vectors of length number of non-dirichlet boundary dofs
            NekVector<NekDouble> F_GlobBnd(nGlobHomBndDofs,pInput,eWrapper);
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,pOutput,
                                          eWrapper);
            //Block transformation matrix
            DNekBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, 0.0);

            // Allocated array of size number of global boundary dofs and copy
            // the input array to the tmp array offset by Dirichlet boundary
            // conditions.
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobHomBndDofs, pInput.get(), 1, tmp.get()
                         + nDirBndDofs, 1);
            
            //Global boundary dofs (with zeroed dirichlet values) to
            //local boundary dofs - This also divides by the mulplicity
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(),
                         tmp.get(), m_map.get(), pLocalIn.get());

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
			    0.0, pLocal.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }
            
            //Assemble local boundary to global non-dirichlet boundary
            m_locToGloMap->AssembleBnd(F_LocBnd,F_HomBnd,nDirBndDofs);
        }

        /**
         * \brief transform the solution vector from low energy back to the
         * original basis.
         *
         * After the conjugate gradient routine the output vector is in the low
         * energy basis and must be trasnformed back to the original basis in
         * order to get the correct solution out. the solution vector
         * i.e. \f$\mathbf{x}=\mathbf{R^{T}}\mathbf{\overline{x}}\f$.
         */
        void PreconditionerLowEnergy::v_DoTransformFromLowEnergy(
            Array<OneD, NekDouble>& pInOut)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            ASSERTL1(pInOut.num_elements() >= nGlobBndDofs,
                     "Output array is greater than the nGlobBndDofs");

            //Block  transformation matrix
            DNekBlkMat &R = *m_RBlk;

            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,
                                              pInOut+nDirBndDofs,
                                              eWrapper);

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, 0.0);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,pLocalIn,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);

            //Global boundary (less dirichlet) to local boundary
            m_locToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd, nDirBndDofs);

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
			    0.0, pLocal.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }

            //Assemble local boundary to global boundary
            Vmath::Assmb(nLocBndDofs, m_locToGloSignMult.get(),pLocal.get(), m_map.get(), tmp.get());

            //Universal assemble across processors
            m_locToGloMap->UniversalAssembleBnd(tmp);

            //copy non-dirichlet boundary values
            Vmath::Vcopy(nGlobBndDofs-nDirBndDofs, tmp.get() + nDirBndDofs, 1, pInOut.get() + nDirBndDofs, 1);
        }

        /**
         * \brief Multiply by the block inverse transformation matrix
         */ 
        void PreconditionerLowEnergy::v_DoMultiplybyInverseTransformationMatrix(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            ASSERTL1(pInput.num_elements() >= nGlobHomBndDofs,
                     "Input array is greater than the nGlobHomBndDofs");
            ASSERTL1(pOutput.num_elements() >= nGlobHomBndDofs,
                     "Output array is greater than the nGlobHomBndDofs");

            //vectors of length number of non-dirichlet boundary dofs
            NekVector<NekDouble> F_GlobBnd(nGlobHomBndDofs,pInput,eWrapper);
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,pOutput,
                                          eWrapper);
            //Block inverse transformation matrix
            DNekBlkMat &invR = *m_InvRBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, 0.0);

            // Allocated array of size number of global boundary dofs and copy
            // the input array to the tmp array offset by Dirichlet boundary
            // conditions.
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobHomBndDofs, pInput.get(), 1, tmp.get() + nDirBndDofs, 1);

            //Global boundary dofs (with zeroed dirichlet values) to
            //local boundary dofs
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(),
                         tmp.get(), m_map.get(), pLocalIn.get());

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
			    0.0, pLocal.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }
            

            //Assemble local boundary to global non-dirichlet boundary
            m_locToGloMap->AssembleBnd(F_LocBnd,F_HomBnd,nDirBndDofs);

	}

        /**
         * \brief Multiply by the block tranposed inverse transformation matrix
         */ 
        void PreconditionerLowEnergy::v_DoMultiplybyInverseTransposedTransformationMatrix(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            ASSERTL1(pInput.num_elements() >= nGlobHomBndDofs,
                     "Input array is greater than the nGlobHomBndDofs");
            ASSERTL1(pOutput.num_elements() >= nGlobHomBndDofs,
                     "Output array is greater than the nGlobHomBndDofs");

            //vectors of length number of non-dirichlet boundary dofs
            NekVector<NekDouble> F_GlobBnd(nGlobHomBndDofs,pInput,eWrapper);
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,pOutput,
                                          eWrapper);
            //Block inverse transformation matrix
            DNekBlkMat &invR = *m_InvRBlk;

            Array<OneD, NekDouble> pLocalIn(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocalIn,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);

            m_locToGloMap->GlobalToLocalBnd(pInput,pLocal, nDirBndDofs);


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
			    0.0, pLocal.get() + cnt, nbndcoeffs);
		cnt  += nbndcoeffs*nexp;
		cnt1 += nexp;
	    }

            
            m_locToGloMap->AssembleBnd(pLocal,pOutput, nDirBndDofs);

            Vmath::Vmul(nGlobHomBndDofs,pOutput,1,m_multiplicity,1,pOutput,1);
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
            const boost::shared_ptr<DNekScalMat > &loc_mat)
	{
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
         
            StdRegions::StdExpansionSharedPtr locExpansion;                
            int eid = expList->GetOffset_Elmt_Id(n);
            locExpansion = expList->GetExp(eid);

            int nbnd=locExpansion->NumBndryCoeffs();
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr pS2 = MemoryManager<DNekMat>::AllocateSharedPtr(nbnd,nbnd,0.0,storage);
            
            //transformation matrices
            DNekMat &R = (*m_RBlk->GetBlock(n,n));

            // Original Schur Complement
            DNekScalMat &S1 = (*loc_mat);

            //create low energy matrix
            DNekMat &S2 = (*pS2);

            S2= R*S1*Transpose(R);

            DNekScalMatSharedPtr return_val;
            return_val = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, pS2);

	    return return_val;
	}

        /**
         * Create the inverse multiplicity map.
         */
        void PreconditionerLowEnergy::CreateMultiplicityMap(void)
        {
            unsigned int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            unsigned int nEntries   = m_locToGloMap->GetNumLocalBndCoeffs();
            unsigned int i;
            
            const Array<OneD, const int> &vMap
                = m_locToGloMap->GetLocalToGlobalBndMap();

            const Array< OneD, const NekDouble > &sign 
                = m_locToGloMap->GetLocalToGlobalBndSign();

            bool m_signChange=m_locToGloMap->GetSignChange();

            // Count the multiplicity of each global DOF on this process
            Array<OneD, NekDouble> vCounts(nGlobalBnd, 0.0);
            for (i = 0; i < nEntries; ++i)
            {
                vCounts[vMap[i]] += 1.0;
            }

            // Get universal multiplicity by globally assembling counts
            m_locToGloMap->UniversalAssembleBnd(vCounts);

            // Construct a map of 1/multiplicity
            m_locToGloSignMult = Array<OneD, NekDouble>(nEntries);
            for (i = 0; i < nEntries; ++i)
            {
                if(m_signChange)
                {
                    m_locToGloSignMult[i] = sign[i]*1.0/vCounts[vMap[i]];
                }
                else
                {
                    m_locToGloSignMult[i] = 1.0/vCounts[vMap[i]];
                }
            }

            int nDirBnd        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBnd    = nGlobalBnd - nDirBnd;
            int nLocBnd        = m_locToGloMap->GetNumLocalBndCoeffs();

            //Set up multiplicity array for inverse transposed transformation matrix
            Array<OneD,NekDouble> tmp(nGlobHomBnd,1.0);
            m_multiplicity = Array<OneD,NekDouble>(nGlobHomBnd,1.0);
            Array<OneD,NekDouble> loc(nLocBnd,1.0);

            m_locToGloMap->GlobalToLocalBnd(tmp,loc, nDirBnd);
            m_locToGloMap->AssembleBnd(loc,m_multiplicity, nDirBnd);
            Vmath::Sdiv(nGlobHomBnd,1.0,m_multiplicity,1,m_multiplicity,1);

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
            
            //boost::shared_ptr<SpatialDomains::PointGeom> verts[6];
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
            const bool   isQuadEdgeFlipped[][4] = { {0,0,1,1}, {0,0,1,1}, {0,0,1,1} };
            // QuadId ordered as 0, 1, 2, otherwise return false
            const int                  quadId[] = { 0,-1,1,-1,2 }; 
            
            //triangle-edge connectivity side-triface-1, side triface-3 
            const int  triEdgeConnectivity[][3] = { {0,5,4}, {2,6,7} };
            const bool    isTriEdgeFlipped[][3] = { {0,0,1}, {0,0,1} };
            // TriId ordered as 0, 1, otherwise return false
            const int                   triId[] = { -1,0,-1,1,-1 }; 
            
            // Populate the list of faces  
            SpatialDomains::Geometry2DSharedPtr faces[nFaces]; 
            for(int f = 0; f < nFaces; ++f){
                if(f == 1 || f == 3) {
                    int i = triId[f];
                    SpatialDomains::SegGeomSharedPtr edgeArray[3];
		    StdRegions::Orientation eorientArray[3];
                    for(int j = 0; j < 3; ++j){
                        edgeArray[j] = edges[triEdgeConnectivity[i][j]];
                        eorientArray[j] = isTriEdgeFlipped[i][j] ? StdRegions::eBackwards : StdRegions::eForwards;
                    }
                    faces[f] = MemoryManager<SpatialDomains::TriGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
                }            
                else {
                    int i = quadId[f];
                    SpatialDomains::SegGeomSharedPtr edgeArray[4];
		    StdRegions::Orientation eorientArray[4]; 
                    for(int j=0; j < 4; ++j){
                        edgeArray[j] = edges[quadEdgeConnectivity[i][j]];
                        eorientArray[j] = isQuadEdgeFlipped[i][j] ? StdRegions::eBackwards : StdRegions::eForwards;
                    }
                    faces[f] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(f, edgeArray, eorientArray);
                }
            } 
            
            SpatialDomains::PrismGeomSharedPtr geom = MemoryManager<SpatialDomains::PrismGeom>::AllocateSharedPtr(faces);

            geom->SetOwnData();

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
            
            boost::shared_ptr<SpatialDomains::PointGeom> verts[4];
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
                boost::shared_ptr<SpatialDomains::PointGeom>
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
            const bool isEdgeFlipped[][3] = {
                {0,0,1}, {0,0,1}, {0,0,1}, {0,0,1}
            };
            
            // Populate the list of faces
            SpatialDomains::TriGeomSharedPtr faces[nFaces];
            for(i=0; i < nFaces; ++i)
            {
                SpatialDomains::SegGeomSharedPtr edgeArray[3];
                StdRegions::Orientation eorientArray[3];
                for(j=0; j < 3; ++j)
                {
                    edgeArray[j] = edges[edgeConnectivity[i][j]];
                    eorientArray[j] = isEdgeFlipped[i][j] ? 
                        StdRegions::eBackwards : StdRegions::eForwards;
                }
                
                
                faces[i] = MemoryManager<SpatialDomains::TriGeom>
                    ::AllocateSharedPtr(i, edgeArray, eorientArray);
            }
            
            SpatialDomains::TetGeomSharedPtr geom =
                MemoryManager<SpatialDomains::TetGeom>::AllocateSharedPtr
                (faces);
            
            geom->SetOwnData();

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
            const bool isEdgeFlipped[][4] = {
                {0,0,0,1}, {0,0,1,1}, {0,0,1,1},
                {0,0,1,1}, {0,0,1,1}, {0,0,0,1}
            };

            // Populate the list of faces
            SpatialDomains::QuadGeomSharedPtr faces[nFaces];
            for( int i = 0; i < nFaces; ++i ) {
                SpatialDomains::SegGeomSharedPtr edgeArray[4];
                StdRegions::Orientation eorientArray[4];
                for( int j = 0; j < 4; ++j ) {
                    edgeArray[j]    = edges[edgeConnectivity[i][j]];
                    eorientArray[j] = isEdgeFlipped[i][j] ? 
                        StdRegions::eBackwards : StdRegions::eForwards;
                }
                faces[i] = MemoryManager<SpatialDomains::QuadGeom>::AllocateSharedPtr(i, edgeArray,
                                                                      eorientArray);
            }

            SpatialDomains::HexGeomSharedPtr geom =
                MemoryManager<SpatialDomains::HexGeom>::AllocateSharedPtr
                (faces);
            
            geom->SetOwnData();

            return geom;
        }


        /**
	 * \brief Loop expansion and determine different variants of the
	 * transformation matrix
	 *
         * Sets up multiple reference elements based on the element expansion. 
	 */
        void PreconditionerLowEnergy::SetUpReferenceElements()
        {
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::StdExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(0);

            DNekScalBlkMatSharedPtr RtetBlk, RprismBlk;
            DNekScalBlkMatSharedPtr RTtetBlk, RTprismBlk;

            DNekScalMatSharedPtr Rprismoriginal;
            DNekScalMatSharedPtr RTprismoriginal;
            DNekMatSharedPtr Rtettmp, RTtettmp, Rhextmp, RThextmp, Rprismtmp, RTprismtmp ;

            /*
             * The building block for all of hte new basis are the Tet
             * and Hex howver if Hex elements are not available then
             * we can use the Tet and Prism expansiosn so set up tet,
             * Prism and Hex elements
             */
            SpatialDomains::TetGeomSharedPtr   tetgeom  = CreateRefTetGeom();
            SpatialDomains::PrismGeomSharedPtr prismgeom= CreateRefPrismGeom();
            SpatialDomains::HexGeomSharedPtr   hexgeom  = CreateRefHexGeom();

            /* Determine the maximum expansion order for all elements */
            m_nummodesmax = 0; 
            for(int n = 0; n < expList->GetNumElmts(); ++n)
            {
                m_nummodesmax = max(m_nummodesmax,expList->GetExp(n)->GetBasisNumModes(0));
                m_nummodesmax = max(m_nummodesmax,expList->GetExp(n)->GetBasisNumModes(1));
                m_nummodesmax = max(m_nummodesmax,expList->GetExp(n)->GetBasisNumModes(2));
            }
            m_comm->AllReduce(m_nummodesmax, LibUtilities::ReduceMax);
            
            /*
             * Set up a transformation matrices for equal max order
             * polynomial meshes
             */

            //Bases for Tetrahedral element
            const LibUtilities::BasisKey TetBa(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey TetBb(
                LibUtilities::eModified_B, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax,
                                        LibUtilities::eGaussRadauMAlpha1Beta0));
            const LibUtilities::BasisKey TetBc(
                LibUtilities::eModified_C, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax,
                                        LibUtilities::eGaussRadauMAlpha2Beta0));

            //Create reference tetrahedral expansion
            LocalRegions::TetExpSharedPtr TetExp;

            TetExp = MemoryManager<LocalRegions::TetExp>
                ::AllocateSharedPtr(TetBa,TetBb,TetBc,
                                    tetgeom);
            
            //Bases for Prism element
            const LibUtilities::BasisKey PrismBa(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));

            const LibUtilities::BasisKey PrismBb(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey PrismBc(
                LibUtilities::eModified_B, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax,
                                        LibUtilities::eGaussRadauMAlpha1Beta0));

            //Create reference prismatic expansion
            LocalRegions::PrismExpSharedPtr PrismExp;

            PrismExp = MemoryManager<LocalRegions::PrismExp>
                ::AllocateSharedPtr(PrismBa,PrismBb,PrismBc,
                                    prismgeom);

            //Bases for Hex element
            const LibUtilities::BasisKey HexBa(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey HexBb(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey HexBc(
                LibUtilities::eModified_A, m_nummodesmax,
                LibUtilities::PointsKey(m_nummodesmax+1,
                                        LibUtilities::eGaussLobattoLegendre));
            
            //Create reference prismatic expansion
            LocalRegions::HexExpSharedPtr HexExp;
            
            HexExp = MemoryManager<LocalRegions::HexExp>
                ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                    hexgeom);
            
            StdRegions::MatrixType PreconR,PreconRT;

            if(m_linSysKey.GetMatrixType() == StdRegions::eMass)
            {
                PreconR  = StdRegions::ePreconRMass;
                PreconRT = StdRegions::ePreconRTMass;
            }
            else
            {
                PreconR  = StdRegions::ePreconR;
                PreconRT = StdRegions::ePreconRT;
            }

            /*
             * Matrix keys - for each element type there are two matrix keys
             * corresponding to the transformation matrix R and its transpose
             */

            //Matrix keys for tetrahedral element transformation matrix
            LocalRegions::MatrixKey TetR
                (PreconR, LibUtilities::eTetrahedron,
                 *TetExp, m_linSysKey.GetConstFactors());

            //Matrix keys for tetrahedral transposed transformation matrix
            LocalRegions::MatrixKey TetRT
                (PreconRT, LibUtilities::eTetrahedron,
                 *TetExp,  m_linSysKey.GetConstFactors());

            //Matrix keys for prismaticl element transformation matrix
            LocalRegions::MatrixKey PrismR
                (PreconR, LibUtilities::ePrism,
                 *PrismExp, m_linSysKey.GetConstFactors());
            
            //Matrix keys for prismatic element transposed transformation
            //matrix
            LocalRegions::MatrixKey PrismRT
                (PreconRT, LibUtilities::ePrism,
                 *PrismExp,  m_linSysKey.GetConstFactors());

            //Matrix keys for hexahedral element transformation matrix
            LocalRegions::MatrixKey HexR
                (PreconR, LibUtilities::eHexahedron,
                 *HexExp, m_linSysKey.GetConstFactors());
            
            //Matrix keys for hexahedral element transposed transformation
            //matrix
            LocalRegions::MatrixKey HexRT
                (PreconRT, LibUtilities::eHexahedron,
                 *HexExp,  m_linSysKey.GetConstFactors());
            
            /*
             * Create transformation matrices for the tetrahedral element
             */
            
            //Get tetrahedral transformation matrix
            m_maxRmat[LibUtilities::eTetrahedron] =
                TetExp->GetLocMatrix(TetR);
            m_edgeMapMaxR[LibUtilities::eTetrahedron] =
                TetExp->GetEdgeInverseBoundaryMap();
            m_faceMapMaxR[LibUtilities::eTetrahedron] =
                TetExp->GetFaceInverseBoundaryMap();
            
            //Get prismatic transformation matrix
            m_maxRmat[LibUtilities::ePrism] =
                PrismExp->GetLocMatrix(PrismR);
            m_edgeMapMaxR[LibUtilities::ePrism] =
                PrismExp->GetEdgeInverseBoundaryMap();
            m_faceMapMaxR[LibUtilities::ePrism] =
                PrismExp->GetFaceInverseBoundaryMap();

            // Note we do not have pyramid here since this has to be
            // constructed from other shapes

            //Get hexahedral transformation matrix
            m_maxRmat[LibUtilities::eHexahedron] =
                HexExp->GetLocMatrix(HexR);
            m_edgeMapMaxR[LibUtilities::eHexahedron] =
                HexExp->GetEdgeInverseBoundaryMap();
            m_faceMapMaxR[LibUtilities::eHexahedron] =
                HexExp->GetFaceInverseBoundaryMap();

#if 0
            -> Need to now have section which takes the maximum matrices on mixed meshes and replaced for example the tet expansion with elements from the hex or if just a tet prism mesh just take elements from the Tet to put in the prism? 
            
            /***********************************/
            /* Setup Tet transformation matrix */
            /***********************************/

            /*
             * Here we take the original R matrix for the lowest
             * polynomial order and replace the values with those from
             * the matrix for the largest polynomail order
             */

            //Get the original number of boundary, edge and face coefficeints
            int nBndCoeffsHet = TetExp->NumBndryCoeffs();
            int nFaceCoeffs = TetExp->GetFaceIntNcoeffs(0);
            int nEdgeCoeffs = TetExp->GetEdgeNcoeffs(0)-2;

            //Allocation the matrix
            DNekMatSharedPtr RtetHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);
            DNekMatSharedPtr RTtetHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);
            DNekMatSharedPtr invRTtetHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);

            //These are the vertex mode locations of R which need to be replaced
            //in the prism element
            int TetVertex0=maxTetExp->GetVertexMap(0);
            int TetVertex1=maxTetExp->GetVertexMap(1);
            int TetVertex2=maxTetExp->GetVertexMap(2);
            int TetVertex3=maxTetExp->GetVertexMap(3);

            Array<OneD, unsigned int> TetEdge0=maxTetExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> TetEdge1=maxTetExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> TetEdge2=maxTetExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> TetEdge3=maxTetExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> TetEdge4=maxTetExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> TetEdge5=maxTetExp->GetEdgeInverseBoundaryMap(5);

            Array<OneD, unsigned int> TetFace0=maxTetExp->GetFaceInverseBoundaryMap(0);
            Array<OneD, unsigned int> TetFace1=maxTetExp->GetFaceInverseBoundaryMap(1);
            Array<OneD, unsigned int> TetFace2=maxTetExp->GetFaceInverseBoundaryMap(2);
            Array<OneD, unsigned int> TetFace3=maxTetExp->GetFaceInverseBoundaryMap(3);

            //Setup R_ve and R_vf

            int nedgemodesconnected = 3*nEdgeCoeffs;
            int nfacemodesconnected = 3*nFaceCoeffs;

            int edgeid, maxvertlocation, vertlocation;

            NekDouble maxVertVertValue;

            //set vertex-vertex values
            for(int vid=0; vid<4; ++vid)
            {
                // Matrix value for each coefficient location
                maxvertlocation=maxTetExp->GetVertexMap(vid);

                maxVertVertValue = (*m_maxRtet)(maxvertlocation,
                                               maxvertlocation);
                vertlocation=TetExp->GetVertexMap(vid);

                // Set the value in the vertex edge/face matrix
                RtetHet->SetValue(vertlocation, vertlocation, maxVertVertValue);
                RTtetHet->SetValue(vertlocation, vertlocation, maxVertVertValue);

                maxVertVertValue = (*m_maxRTinvtet)(maxvertlocation,
                                               maxvertlocation);
                invRTtetHet->SetValue(vertlocation, vertlocation, maxVertVertValue);
            }

            NekDouble maxVertEdgeValue;

            //set vertex edge values
            for(int vid=0; vid<4; ++vid)
            {
                maxvertlocation=maxTetExp->GetVertexMap(vid);
                vertlocation=TetExp->GetVertexMap(vid);

                // Three attached edges
                int edgeid1=tetgeom->GetVertexEdgeMap(vid, 0);
                int edgeid2=tetgeom->GetVertexEdgeMap(vid, 1);
                int edgeid3=tetgeom->GetVertexEdgeMap(vid, 2);

                Array<OneD, unsigned int> maxTetEdge0=maxTetExp->GetEdgeInverseBoundaryMap(edgeid1);
                Array<OneD, unsigned int> maxTetEdge1=maxTetExp->GetEdgeInverseBoundaryMap(edgeid2);
                Array<OneD, unsigned int> maxTetEdge2=maxTetExp->GetEdgeInverseBoundaryMap(edgeid3);

                Array<OneD, unsigned int> TetEdge0=TetExp->GetEdgeInverseBoundaryMap(edgeid1);
                Array<OneD, unsigned int> TetEdge1=TetExp->GetEdgeInverseBoundaryMap(edgeid2);
                Array<OneD, unsigned int> TetEdge2=TetExp->GetEdgeInverseBoundaryMap(edgeid3);

                for(int i=0; i<TetEdge0.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRtet)(maxvertlocation,
                                                    maxTetEdge0[i]);

                    RtetHet->SetValue(vertlocation, TetEdge0[i], maxVertEdgeValue);
                    RTtetHet->SetValue(TetEdge0[i], vertlocation,maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvtet)(maxTetEdge0[i],maxvertlocation);
                    invRTtetHet->SetValue(TetEdge0[i],vertlocation, maxVertEdgeValue);
                }

                for(int i=0; i<TetEdge1.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRtet)(maxvertlocation,
                                                    maxTetEdge1[i]);

                    RtetHet->SetValue(vertlocation, TetEdge1[i], maxVertEdgeValue);
                    RTtetHet->SetValue(TetEdge1[i],vertlocation, maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvtet)(maxTetEdge1[i],maxvertlocation);
                    invRTtetHet->SetValue(TetEdge1[i],vertlocation, maxVertEdgeValue);
                }

                for(int i=0; i<TetEdge2.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRtet)(maxvertlocation,
                                                    maxTetEdge2[i]);

                    RtetHet->SetValue(vertlocation, TetEdge2[i], maxVertEdgeValue);
                    RTtetHet->SetValue(TetEdge2[i],vertlocation, maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvtet)(maxTetEdge2[i],maxvertlocation);
                    invRTtetHet->SetValue(TetEdge2[i],vertlocation, maxVertEdgeValue);
                }


            }

            NekDouble maxVertFaceValue;

            //set vertex face values
            for(int vid=0; vid<4; ++vid)
            {
                maxvertlocation=maxTetExp->GetVertexMap(vid);
                vertlocation=TetExp->GetVertexMap(vid);

                // Three attached edges
                int faceid1=tetgeom->GetVertexFaceMap(vid, 0);
                int faceid2=tetgeom->GetVertexFaceMap(vid, 1);
                int faceid3=tetgeom->GetVertexFaceMap(vid, 2);

                Array<OneD, unsigned int> maxTetFace0=maxTetExp->GetFaceInverseBoundaryMap(faceid1);
                Array<OneD, unsigned int> maxTetFace1=maxTetExp->GetFaceInverseBoundaryMap(faceid2);
                Array<OneD, unsigned int> maxTetFace2=maxTetExp->GetFaceInverseBoundaryMap(faceid3);

                Array<OneD, unsigned int> TetFace0=TetExp->GetFaceInverseBoundaryMap(faceid1);
                Array<OneD, unsigned int> TetFace1=TetExp->GetFaceInverseBoundaryMap(faceid2);
                Array<OneD, unsigned int> TetFace2=TetExp->GetFaceInverseBoundaryMap(faceid3);


                int nfacemodes=TetFace0.num_elements();

                Array<OneD, unsigned int> ExtractTetFace0(nfacemodes);
                Array<OneD, unsigned int> ExtractTetFace1(nfacemodes);
                Array<OneD, unsigned int> ExtractTetFace2(nfacemodes);

                int offset=0;
                int cnt=0;
                for(int k=0; k<nummodes0-3; ++k)
                {
                    for(int j=0; j<nummodes0-3-k; ++j)
                    {
                        ExtractTetFace0[cnt]=maxTetFace0[offset+j];
                        ExtractTetFace1[cnt]=maxTetFace1[offset+j];
                        ExtractTetFace2[cnt]=maxTetFace2[offset+j];
                        cnt++;
                    }
                    offset+=maxnummodes-3-k;
                }

                for(int i=0; i<TetFace0.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRtet)(maxvertlocation,
                                                    ExtractTetFace0[i]);

                    RtetHet->SetValue(vertlocation, TetFace0[i], maxVertFaceValue);
                    RTtetHet->SetValue(TetFace0[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvtet)(ExtractTetFace0[i],maxvertlocation);
                    invRTtetHet->SetValue(TetFace0[i], vertlocation,  maxVertFaceValue);
                }

                for(int i=0; i<TetFace1.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRtet)(maxvertlocation,
                                                    ExtractTetFace1[i]);

                    RtetHet->SetValue(vertlocation, TetFace1[i], maxVertFaceValue);
                    RTtetHet->SetValue(TetFace1[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvtet)(ExtractTetFace1[i],maxvertlocation);
                    invRTtetHet->SetValue(TetFace1[i], vertlocation,  maxVertFaceValue);



                }

                for(int i=0; i<TetFace2.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRtet)(maxvertlocation,
                                                    ExtractTetFace2[i]);

                    RtetHet->SetValue(vertlocation, TetFace2[i], maxVertFaceValue);
                    RTtetHet->SetValue(TetFace2[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvtet)(ExtractTetFace2[i],maxvertlocation);
                    invRTtetHet->SetValue(TetFace2[i], vertlocation,  maxVertFaceValue);

                }
            }

            NekDouble maxEdgeFaceValue;
            //set egde-edge values
            for(int eid=0; eid<6; ++eid)
            {
                int faceid0=tetgeom->GetEdgeFaceMap(eid,0);	
                int faceid1=tetgeom->GetEdgeFaceMap(eid,1);

                Array<OneD, unsigned int> maxTetEdge0=maxTetExp->GetEdgeInverseBoundaryMap(eid);
                Array<OneD, unsigned int> TetEdge0=TetExp->GetEdgeInverseBoundaryMap(eid);

                Array<OneD, unsigned int> maxTetFace0=maxTetExp->GetFaceInverseBoundaryMap(faceid0);
                Array<OneD, unsigned int> maxTetFace1=maxTetExp->GetFaceInverseBoundaryMap(faceid1);

                Array<OneD, unsigned int> TetFace0=TetExp->GetFaceInverseBoundaryMap(faceid0);
                Array<OneD, unsigned int> TetFace1=TetExp->GetFaceInverseBoundaryMap(faceid1);

                int nfacemodes=TetFace0.num_elements();

                Array<OneD, unsigned int> ExtractTetFace0(nfacemodes);
                Array<OneD, unsigned int> ExtractTetFace1(nfacemodes);

                int offset=0;
                int cnt=0;
                for(int k=0; k<nummodes0-3; ++k)
                {
                    for(int j=0; j<nummodes0-3-k; ++j)
                    {
                        ExtractTetFace0[cnt]=maxTetFace0[offset+j];
                        ExtractTetFace1[cnt]=maxTetFace1[offset+j];
                        cnt++;
                    }
                    offset+=maxnummodes-3-k;
                }

                for(int i=0; i<TetEdge0.num_elements(); ++i)
                {
                    for(int j=0; j<TetFace0.num_elements(); ++j)
                    {
                        maxEdgeFaceValue = (*m_maxRtet)(maxTetEdge0[i],
                                                        ExtractTetFace0[j]);
                        RtetHet->SetValue(TetEdge0[i], TetFace0[j], maxEdgeFaceValue);
                        RTtetHet->SetValue(TetFace0[j],TetEdge0[i], maxEdgeFaceValue);

                        maxEdgeFaceValue = (*m_maxRTinvtet)(ExtractTetFace0[j],maxTetEdge0[i]);
                        invRTtetHet->SetValue(TetFace0[j],TetEdge0[i], maxEdgeFaceValue);
                    }
                }

                for(int i=0; i<TetEdge0.num_elements(); ++i)
                {
                    for(int j=0; j<TetFace1.num_elements(); ++j)
                    {
                        maxEdgeFaceValue = (*m_maxRtet)(maxTetEdge0[i],
                                                        ExtractTetFace1[j]);
                        RtetHet->SetValue(TetEdge0[i], TetFace1[j], maxEdgeFaceValue);
                        RTtetHet->SetValue(TetFace1[j],TetEdge0[i], maxEdgeFaceValue);

                        maxEdgeFaceValue = (*m_maxRTinvtet)(ExtractTetFace1[j],maxTetEdge0[i]);
                        invRTtetHet->SetValue(TetFace1[j],TetEdge0[i], maxEdgeFaceValue);

                    }
                }
            }

            for (i = 0; i < RtetHet->GetRows(); ++i)
            {
                RtetHet->SetValue(i, i, 1.0);
                RTtetHet->SetValue(i, i, 1.0);
                invRTtetHet->SetValue(i, i, 1.0);
            }


            m_Rtet = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RtetHet);
            m_RTtet = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RTtetHet);
            m_RTinvtet = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,invRTtetHet);

            /***********************************/
            /* Setup Hex transformation matrix */
            /***********************************/

            //Get the original number of boundary, edge and face coefficeints
            nBndCoeffsHet = HexExp->NumBndryCoeffs();  //change
            nFaceCoeffs = HexExp->GetFaceIntNcoeffs(0); //change
            nEdgeCoeffs = HexExp->GetEdgeNcoeffs(0)-2; //change

            //Allocation the matrix
            DNekMatSharedPtr RhexHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);
            DNekMatSharedPtr RThexHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);
            DNekMatSharedPtr invRThexHet = MemoryManager<DNekMat>::
                AllocateSharedPtr(nBndCoeffsHet,nBndCoeffsHet,zero,eFULL);

            //These are the vertex mode locations of R which need to be replaced
            //in the hex element

            //NekDouble maxVertVertValue;

            //set vertex-vertex values
            for(int vid=0; vid<8; ++vid)
            {
                // Matrix value for each coefficient location
                maxvertlocation=maxHexExp->GetVertexMap(vid);

                maxVertVertValue = (*m_maxRhex)(maxvertlocation,
                                               maxvertlocation);
                vertlocation=HexExp->GetVertexMap(vid);

                // Set the value in the vertex edge/face matrix
                RhexHet->SetValue(vertlocation, vertlocation, maxVertVertValue);
                RThexHet->SetValue(vertlocation, vertlocation, maxVertVertValue);

                maxVertVertValue = (*m_maxRTinvhex)(maxvertlocation,
                                               maxvertlocation);
                invRThexHet->SetValue(vertlocation, vertlocation, maxVertVertValue);
            }

            //NekDouble maxVertEdgeValue;

            //set vertex edge values
            for(int vid=0; vid<8; ++vid)
            {
                maxvertlocation=maxHexExp->GetVertexMap(vid);
                vertlocation=HexExp->GetVertexMap(vid);

                // Three attached edges
                int edgeid1=hexgeom->GetVertexEdgeMap(vid, 0);
                int edgeid2=hexgeom->GetVertexEdgeMap(vid, 1);
                int edgeid3=hexgeom->GetVertexEdgeMap(vid, 2);

                Array<OneD, unsigned int> maxHexEdge0=maxHexExp->GetEdgeInverseBoundaryMap(edgeid1);
                Array<OneD, unsigned int> maxHexEdge1=maxHexExp->GetEdgeInverseBoundaryMap(edgeid2);
                Array<OneD, unsigned int> maxHexEdge2=maxHexExp->GetEdgeInverseBoundaryMap(edgeid3);

                Array<OneD, unsigned int> HexEdge0=HexExp->GetEdgeInverseBoundaryMap(edgeid1);
                Array<OneD, unsigned int> HexEdge1=HexExp->GetEdgeInverseBoundaryMap(edgeid2);
                Array<OneD, unsigned int> HexEdge2=HexExp->GetEdgeInverseBoundaryMap(edgeid3);

                for(int i=0; i<HexEdge0.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRhex)(maxvertlocation,
                                                    maxHexEdge0[i]);

                    RhexHet->SetValue(vertlocation, HexEdge0[i], maxVertEdgeValue);
                    RThexHet->SetValue(HexEdge0[i], vertlocation,maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvhex)(maxHexEdge0[i],maxvertlocation);
                    invRThexHet->SetValue(HexEdge0[i],vertlocation, maxVertEdgeValue);
                }

                for(int i=0; i<HexEdge1.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRhex)(maxvertlocation,
                                                    maxHexEdge1[i]);

                    RhexHet->SetValue(vertlocation, HexEdge1[i], maxVertEdgeValue);
                    RThexHet->SetValue(HexEdge1[i],vertlocation, maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvhex)(maxHexEdge1[i],maxvertlocation);
                    invRThexHet->SetValue(HexEdge1[i],vertlocation, maxVertEdgeValue);
                }

                for(int i=0; i<HexEdge2.num_elements(); ++i)
                {
                    maxVertEdgeValue = (*m_maxRhex)(maxvertlocation,
                                                    maxHexEdge2[i]);

                    RhexHet->SetValue(vertlocation, HexEdge2[i], maxVertEdgeValue);
                    RThexHet->SetValue(HexEdge2[i],vertlocation, maxVertEdgeValue);

                    maxVertEdgeValue = (*m_maxRTinvhex)(maxHexEdge2[i],maxvertlocation);
                    invRThexHet->SetValue(HexEdge2[i],vertlocation, maxVertEdgeValue);
                }


            }

            //set vertex face values
            for(int vid=0; vid<8; ++vid)
            {
                maxvertlocation=maxHexExp->GetVertexMap(vid);
                vertlocation=HexExp->GetVertexMap(vid);

                // Three attached edges
                int faceid1=hexgeom->GetVertexFaceMap(vid, 0);
                int faceid2=hexgeom->GetVertexFaceMap(vid, 1);
                int faceid3=hexgeom->GetVertexFaceMap(vid, 2);

                Array<OneD, unsigned int> maxHexFace0=maxHexExp->GetFaceInverseBoundaryMap(faceid1);
                Array<OneD, unsigned int> maxHexFace1=maxHexExp->GetFaceInverseBoundaryMap(faceid2);
                Array<OneD, unsigned int> maxHexFace2=maxHexExp->GetFaceInverseBoundaryMap(faceid3);

                Array<OneD, unsigned int> HexFace0=HexExp->GetFaceInverseBoundaryMap(faceid1);
                Array<OneD, unsigned int> HexFace1=HexExp->GetFaceInverseBoundaryMap(faceid2);
                Array<OneD, unsigned int> HexFace2=HexExp->GetFaceInverseBoundaryMap(faceid3);
                

                int nfacemodes=HexFace0.num_elements();

                Array<OneD, unsigned int> ExtractHexFace0(nfacemodes);
                Array<OneD, unsigned int> ExtractHexFace1(nfacemodes);
                Array<OneD, unsigned int> ExtractHexFace2(nfacemodes);

                int offset=0;
                int cnt=0;
                for(int k=0; k<nummodes0-2; ++k)
                {
                    for(int j=0; j<nummodes0-2; ++j)
                    {
                        ExtractHexFace0[cnt]=maxHexFace0[offset+j];
                        ExtractHexFace1[cnt]=maxHexFace1[offset+j];
                        ExtractHexFace2[cnt]=maxHexFace2[offset+j];
                        cnt++;
                    }
                    offset+=maxnummodes-2;
                }

                for(int i=0; i<HexFace0.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRhex)(maxvertlocation,
                                                    ExtractHexFace0[i]);

                    RhexHet->SetValue(vertlocation, HexFace0[i], maxVertFaceValue);
                    RThexHet->SetValue(HexFace0[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvhex)(ExtractHexFace0[i],maxvertlocation);
                    invRThexHet->SetValue(HexFace0[i], vertlocation,  maxVertFaceValue);
                }

                for(int i=0; i<HexFace1.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRhex)(maxvertlocation,
                                                    ExtractHexFace1[i]);

                    RhexHet->SetValue(vertlocation, HexFace1[i], maxVertFaceValue);
                    RThexHet->SetValue(HexFace1[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvhex)(ExtractHexFace1[i],maxvertlocation);
                    invRThexHet->SetValue(HexFace1[i], vertlocation,  maxVertFaceValue);
                }

                for(int i=0; i<HexFace2.num_elements(); ++i)
                {
                    maxVertFaceValue = (*m_maxRhex)(maxvertlocation,
                                                    ExtractHexFace2[i]);

                    RhexHet->SetValue(vertlocation, HexFace2[i], maxVertFaceValue);
                    RThexHet->SetValue(HexFace2[i], vertlocation,  maxVertFaceValue);

                    maxVertFaceValue = (*m_maxRTinvhex)(ExtractHexFace2[i],maxvertlocation);
                    invRThexHet->SetValue(HexFace2[i], vertlocation,  maxVertFaceValue);
                }
            }

            //set egde-edge values
            for(int eid=0; eid<12; ++eid)
            {
                int faceid0=hexgeom->GetEdgeFaceMap(eid,0);	
                int faceid1=hexgeom->GetEdgeFaceMap(eid,1);

                Array<OneD, unsigned int> maxHexEdge0=maxHexExp->GetEdgeInverseBoundaryMap(eid);
                Array<OneD, unsigned int> HexEdge0=HexExp->GetEdgeInverseBoundaryMap(eid);

                Array<OneD, unsigned int> maxHexFace0=maxHexExp->GetFaceInverseBoundaryMap(faceid0);
                Array<OneD, unsigned int> maxHexFace1=maxHexExp->GetFaceInverseBoundaryMap(faceid1);

                Array<OneD, unsigned int> HexFace0=HexExp->GetFaceInverseBoundaryMap(faceid0);
                Array<OneD, unsigned int> HexFace1=HexExp->GetFaceInverseBoundaryMap(faceid1);



                int nfacemodes=HexFace0.num_elements();
                Array<OneD, unsigned int> ExtractHexFace0(nfacemodes);
                Array<OneD, unsigned int> ExtractHexFace1(nfacemodes);

                int offset=0;
                int cnt=0;
                for(int k=0; k<nummodes0-2; ++k)
                {
                    for(int j=0; j<nummodes0-2; ++j)
                    {
                        ExtractHexFace0[cnt]=maxHexFace0[offset+j];
                        ExtractHexFace1[cnt]=maxHexFace1[offset+j];
                        cnt++;
                    }
                    offset+=maxnummodes-2;
                }

                for(int i=0; i<HexEdge0.num_elements(); ++i)
                {
                    for(int j=0; j<HexFace0.num_elements(); ++j)
                    {
                        maxEdgeFaceValue = (*m_maxRhex)(maxHexEdge0[i],
                                                        ExtractHexFace0[j]);
                        RhexHet->SetValue(HexEdge0[i], HexFace0[j], maxEdgeFaceValue);
                        RThexHet->SetValue(HexFace0[j],HexEdge0[i], maxEdgeFaceValue);

                        maxEdgeFaceValue = (*m_maxRTinvhex)(ExtractHexFace0[j],maxHexEdge0[i]);
                        invRThexHet->SetValue(HexFace0[j],HexEdge0[i], maxEdgeFaceValue);
                    }
                }

                for(int i=0; i<HexEdge0.num_elements(); ++i)
                {
                    for(int j=0; j<HexFace1.num_elements(); ++j)
                    {
                        maxEdgeFaceValue = (*m_maxRhex)(maxHexEdge0[i],
                                                        ExtractHexFace1[j]);
                        RhexHet->SetValue(HexEdge0[i], HexFace1[j], maxEdgeFaceValue);
                        RThexHet->SetValue(HexFace1[j],HexEdge0[i], maxEdgeFaceValue);

                        maxEdgeFaceValue = (*m_maxRTinvhex)(ExtractHexFace1[j],maxHexEdge0[i]);
                        invRThexHet->SetValue(HexFace1[j],HexEdge0[i], maxEdgeFaceValue);

                    }
                }
            }

            for (i = 0; i < RhexHet->GetRows(); ++i)
            {
                RhexHet->SetValue(i, i, 1.0);
                RThexHet->SetValue(i, i, 1.0);
                invRThexHet->SetValue(i, i, 1.0);
            }

            //LocalTransformToLowEnergy(m_RThex, HexExp);

            m_Rhex = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RhexHet);
            m_RThex = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RThexHet);
            m_RTinvhex = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,invRThexHet);
#endif
        }

        /**
         * \brief transform the solution vector from low energy back to the
         * original basis.
         *
         * After the conjugate gradient routine the output vector is in the low
         * energy basis and must be trasnformed back to the original basis in
         * order to get the correct solution out. the solution vector
         * i.e. \f$\mathbf{x}=\mathbf{R^{T}}\mathbf{\overline{x}}\f$.
         */
        void PreconditionerLowEnergy::LocalTransformToLowEnergy(
            DNekScalMatSharedPtr RTmat,
            LocalRegions::HexExpSharedPtr maxTetExp)
        {
            int mode;

            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();

            LibUtilities::SessionReaderSharedPtr vSession
                                            = expList->GetSession();

            vSession->LoadParameter("mode", mode);
            cout<<"mode: "<<mode<<endl;

            int nRows=RTmat->GetRows();

            /*for(int i=0; i<RTmat->GetRows(); ++i)
            {
                for(int j=0; j<RTmat->GetColumns(); ++j)
                {
                    cout<<(*RTmat)(i,j)<<" ";
                }
                cout<<endl;
            }
            cout<<endl;*/ 


            Array<OneD, NekDouble> u1(nRows,0.0);
            Array<OneD, NekDouble> u2(nRows,0.0);
            Array<OneD, NekDouble> u2phys;

            NekVector<NekDouble> u1vec(nRows,u1,eWrapper);
            NekVector<NekDouble> u2vec(nRows,u2,eWrapper);

            u1[mode]=1.0;

            DNekScalMat RT=(*RTmat);
            //Multiply by the transposed transformation matrix
            u2vec=RT*u1vec;
            //zero the modes that are not in the lower order solution
            u2phys=expList->UpdatePhys();
            expList->BwdTrans(u2,u2phys);

            int npoints = expList->GetNpoints();

            //------------------------------------------------
            //Coordinate arrays
            Array<OneD, NekDouble> x1(npoints,0.0);
            Array<OneD, NekDouble> y1(npoints,0.0);
            Array<OneD, NekDouble> z1(npoints,0.0);
            expList->GetCoords(x1,y1,z1);

            /*cout<<"x,y,z,u"<<endl;
            for(int i=0; i<u2phys.num_elements(); ++i)
            {
                cout<<x1[i]<<","<<y1[i]<<","<<z1[i]<<","<<u2phys[i]<<endl;
                }*/

            int nCoeffs=expList->GetNcoeffs();

            cout<<"ncoeffs: "<<nCoeffs<<" nRows: "<<nRows<<endl;
            //expList->FwdTrans(u2phys, expList->UpdateCoeffs());

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = expList->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

            for(int i = 0; i < FieldDef.size(); ++i)
            {
                FieldDef[i]->m_fields.push_back("u");
                expList->AppendFieldData(FieldDef[i], FieldData[i], u2);
            }
            LibUtilities::Write("out.fld", FieldDef, FieldData);

        }
#if 0 
        /**
         * \brief Modify the prism transformation matrix to align with the
         * tetrahedral modes.
         *
         * This routine replaces the edge and triangular face components of the
         * prismatic vertex transformation matrices \f$\mathbf{R}_{ve}\f$ and
         * \f$\mathbf{R}_{vf}\f$ with the corresponding components from the
         * tetrahedral transformation matrices. Additionally, triangular face
         * components in the prismatic edge transformation matrix
         * \f$\mathbf{R}_{ef}\f$ with the corresponding component from the
         * tetrahedral transformation matrix.
         */
        void PreconditionerLowEnergy::ModifyPrismTransformationMatrix(
            LocalRegions::TetExpSharedPtr TetExp,
            LocalRegions::PrismExpSharedPtr PrismExp,
            DNekMatSharedPtr Rmodprism,
            DNekMatSharedPtr RTmodprism)
        {
            NekDouble Rvalue, RTvalue;
            int i, j;

            //For a tet element the bottom face is made up of the following:
            //vertices: 0, 1 and 2 edges: 0, 1 and 2 face: 0. We first need to
            //determine the mode locations of these vertices, edges and face so
            //we can extract the correct values from the tetrahedral R matrix.

            //These are the vertex mode locations of R which need to be replaced
            //in the prism element
            int TetVertex0=TetExp->GetVertexMap(0);
            int TetVertex1=TetExp->GetVertexMap(1);
            int TetVertex2=TetExp->GetVertexMap(2);
            int TetVertex3=TetExp->GetVertexMap(3);

            cout<<"Tet vertex 0: "<<TetVertex0<<endl;
            cout<<"Tet vertex 1: "<<TetVertex1<<endl;
            cout<<"Tet vertex 2: "<<TetVertex2<<endl;
            cout<<"Tet vertex 3: "<<TetVertex3<<endl;
            cout<<endl;

            //These are the edge mode locations of R which need to be replaced
            //in the prism element
            Array<OneD, unsigned int> TetEdge0=TetExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> TetEdge1=TetExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> TetEdge2=TetExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> TetEdge3=TetExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> TetEdge4=TetExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> TetEdge5=TetExp->GetEdgeInverseBoundaryMap(5);

            cout<<"Tet edge 0: ";
            for(int i=0; i< TetEdge0.num_elements(); ++i)
            {
                cout<<TetEdge0[i]<<" ";
            }
            cout<<endl;

            cout<<"Tet edge 1: ";
            for(int i=0; i< TetEdge1.num_elements(); ++i)
            {
                cout<<TetEdge1[i]<<" ";
            }
            cout<<endl;

            cout<<"Tet edge 2: ";
            for(int i=0; i< TetEdge2.num_elements(); ++i)
            {
                cout<<TetEdge2[i]<<" ";
            }
            cout<<endl;

            cout<<"Tet edge 3: ";
            for(int i=0; i< TetEdge3.num_elements(); ++i)
            {
                cout<<TetEdge3[i]<<" ";
            }
            cout<<endl;

            cout<<"Tet edge 4: ";
            for(int i=0; i< TetEdge4.num_elements(); ++i)
            {
                cout<<TetEdge4[i]<<" ";
            }
            cout<<endl;

            cout<<"Tet edge 5: ";
            for(int i=0; i< TetEdge5.num_elements(); ++i)
            {
                cout<<TetEdge5[i]<<" ";
            }
            cout<<endl;
            cout<<endl;


            //These are the face mode locations of R which need to be replaced
            //in the prism element
            Array<OneD, unsigned int> TetFace=TetExp->GetFaceInverseBoundaryMap(1);

            cout<<"Tet face 1: ";
            for(int i=0; i< TetFace.num_elements(); ++i)
            {
                cout<<TetFace[i]<<" ";
            }
            cout<<endl;
            cout<<endl;



            //Prism vertex modes
            int PrismVertex0=PrismExp->GetVertexMap(0);
            int PrismVertex1=PrismExp->GetVertexMap(1);
            int PrismVertex2=PrismExp->GetVertexMap(2);
            int PrismVertex3=PrismExp->GetVertexMap(3);
            int PrismVertex4=PrismExp->GetVertexMap(4);
            int PrismVertex5=PrismExp->GetVertexMap(5);

            cout<<"Prism vertex 0: "<<PrismVertex0<<endl;
            cout<<"Prism vertex 1: "<<PrismVertex1<<endl;
            cout<<"Prism vertex 2: "<<PrismVertex2<<endl;
            cout<<"Prism vertex 3: "<<PrismVertex3<<endl;
            cout<<"Prism vertex 4: "<<PrismVertex4<<endl;
            cout<<"Prism vertex 5: "<<PrismVertex5<<endl;
            cout<<endl;

            //Prism edge modes
            Array<OneD, unsigned int> PrismEdge0=
                PrismExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> PrismEdge1=
                PrismExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> PrismEdge2=
                PrismExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> PrismEdge3=
                PrismExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> PrismEdge4=
                PrismExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> PrismEdge5=
                PrismExp->GetEdgeInverseBoundaryMap(5);
            Array<OneD, unsigned int> PrismEdge6=
                PrismExp->GetEdgeInverseBoundaryMap(6);
            Array<OneD, unsigned int> PrismEdge7=
                PrismExp->GetEdgeInverseBoundaryMap(7);
            Array<OneD, unsigned int> PrismEdge8=
                PrismExp->GetEdgeInverseBoundaryMap(8);

            cout<<"Prism edge 0: ";
            for(int i=0; i< PrismEdge0.num_elements(); ++i)
            {
                cout<<PrismEdge0[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 1: ";
            for(int i=0; i< PrismEdge1.num_elements(); ++i)
            {
                cout<<PrismEdge1[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 2: ";
            for(int i=0; i< PrismEdge2.num_elements(); ++i)
            {
                cout<<PrismEdge2[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 3: ";
            for(int i=0; i< PrismEdge3.num_elements(); ++i)
            {
                cout<<PrismEdge3[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 4: ";
            for(int i=0; i< PrismEdge4.num_elements(); ++i)
            {
                cout<<PrismEdge4[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 5: ";
            for(int i=0; i< PrismEdge5.num_elements(); ++i)
            {
                cout<<PrismEdge5[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 6: ";
            for(int i=0; i< PrismEdge6.num_elements(); ++i)
            {
                cout<<PrismEdge6[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 7: ";
            for(int i=0; i< PrismEdge7.num_elements(); ++i)
            {
                cout<<PrismEdge7[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism edge 8: ";
            for(int i=0; i< PrismEdge8.num_elements(); ++i)
            {
                cout<<PrismEdge8[i]<<" ";
            }
            cout<<endl;
            cout<<endl;


            //Prism face 1 & 3 face modes
            Array<OneD, unsigned int> PrismFace1=
                PrismExp->GetFaceInverseBoundaryMap(1);
            Array<OneD, unsigned int> PrismFace3=
                PrismExp->GetFaceInverseBoundaryMap(3);
            Array<OneD, unsigned int> PrismFace0=
                PrismExp->GetFaceInverseBoundaryMap(0);
            Array<OneD, unsigned int> PrismFace2=
                PrismExp->GetFaceInverseBoundaryMap(2);
            Array<OneD, unsigned int> PrismFace4=
                PrismExp->GetFaceInverseBoundaryMap(4);

            cout<<"Prism face 0: ";
            for(int i=0; i< PrismFace0.num_elements(); ++i)
            {
                cout<<PrismFace0[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism face 1: ";
            for(int i=0; i< PrismFace1.num_elements(); ++i)
            {
                cout<<PrismFace1[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism face 2: ";
            for(int i=0; i< PrismFace2.num_elements(); ++i)
            {
                cout<<PrismFace2[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism face 3: ";
            for(int i=0; i< PrismFace3.num_elements(); ++i)
            {
                cout<<PrismFace3[i]<<" ";
            }
            cout<<endl;

            cout<<"Prism face 4: ";
            for(int i=0; i< PrismFace4.num_elements(); ++i)
            {
                cout<<PrismFace4[i]<<" ";
            }
            cout<<endl;
            cout<<endl;

            /*DNekScalMat &RTET=(*m_Rtet);
            cout<<"Tet R matrix"<<" rows: "<<RTET.GetRows()<<endl;
            for(i=0; i<RTET.GetRows(); ++i)
            {
                for(j=0; j<RTET.GetColumns(); ++j)
                {
                    cout<<RTET(i,j)<<" ";
                }
                cout<<endl;
            }
            cout<<endl; 


            cout<<"Prism R matrix"<<" rows: "<<Rmodprism->GetRows()<<endl;
            for(i=0; i<Rmodprism->GetRows(); ++i)
            {
                for(j=0; j<Rmodprism->GetColumns(); ++j)
                {
                    cout<<(*Rmodprism)(i,j)<<" ";
                }
                cout<<endl;
            }
            cout<<endl;*/ 


            //vertex 0 edge 0 3 & 4
            for(i=0; i< PrismEdge0.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex0,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge0[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex0,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge3[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex0,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge4[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge0[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge0[i],PrismVertex0,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge2[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge3[i],PrismVertex0,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge3[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge4[i],PrismVertex0,RTvalue);
            }

            //vertex 1 edge 0 1 & 5
            for(i=0; i< PrismEdge1.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex1,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge0[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex1,TetEdge1[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge1[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex1,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge5[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge0[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge0[i],PrismVertex1,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge1[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge1[i],PrismVertex1,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge4[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge5[i],PrismVertex1,RTvalue);
            }

            //vertex 2 edge 1 2 & 6
            for(i=0; i< PrismEdge2.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex2,TetEdge1[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge1[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex1,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge2[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex2,TetEdge5[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge6[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge1[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge1[i],PrismVertex2,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge0[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge2[i],PrismVertex2,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge5[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge6[i],PrismVertex2,RTvalue);
            }

            //vertex 3 edge 3 2 & 7
            for(i=0; i< PrismEdge3.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex2,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge3[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex0,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge2[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex2,TetEdge5[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge7[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge2[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge3[i],PrismVertex3,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge0[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge2[i],PrismVertex3,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge5[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge7[i],PrismVertex3,RTvalue);
            }

            //vertex 4 edge 4 5 & 8
            for(i=0; i< PrismEdge4.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex3,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge4[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex3,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge5[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex0,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge8[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge3[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge4[i],PrismVertex4,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge4[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge5[i],PrismVertex4,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge2[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge8[i],PrismVertex4,RTvalue);
            }

            //vertex 5 edge 6 7 & 8
            for(i=0; i< PrismEdge5.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex3,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge6[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex3,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge7[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex2,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge8[i],Rvalue);

                //transposed values
                RTvalue=(*m_RTtet)(TetEdge3[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge6[i],PrismVertex5,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge4[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge7[i],PrismVertex5,RTvalue);
                RTvalue=(*m_RTtet)(TetEdge2[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge8[i],PrismVertex5,RTvalue);
            }

            // face 1 vertices 0 1 4
            for(i=0; i< PrismFace1.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex0,TetFace[i]);
                Rmodprism->SetValue(PrismVertex0,PrismFace1[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex1,TetFace[i]);
                Rmodprism->SetValue(PrismVertex1,PrismFace1[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex3,TetFace[i]);
                Rmodprism->SetValue(PrismVertex4,PrismFace1[i],Rvalue);
                
                //transposed values
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex0);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex0,RTvalue);
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex1);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex1,RTvalue);
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex3);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex4,RTvalue);
            }

            // face 3 vertices 2, 3 & 5
            for(i=0; i< PrismFace3.num_elements(); ++i)
            {
                Rvalue=(*m_Rtet)(TetVertex1,TetFace[i]);
                Rmodprism->SetValue(PrismVertex2,PrismFace3[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex0,TetFace[i]);
                Rmodprism->SetValue(PrismVertex3,PrismFace3[i],Rvalue);
                Rvalue=(*m_Rtet)(TetVertex3,TetFace[i]);
                Rmodprism->SetValue(PrismVertex5,PrismFace3[i],Rvalue);
                
                //transposed values
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex1);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex2,RTvalue);
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex0);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex3,RTvalue);
                RTvalue=(*m_RTtet)(TetFace[i],TetVertex3);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex5,RTvalue);
            }

            // Face 1 edge 0 4 5
            for(i=0; i< PrismFace1.num_elements(); ++i)
            {
                for(j=0; j<PrismEdge0.num_elements(); ++j)
                {
                    Rvalue=(*m_Rtet)(TetEdge0[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge0[j],PrismFace1[i],Rvalue);
                    Rvalue=(*m_Rtet)(TetEdge3[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge4[j],PrismFace1[i],Rvalue);
                    Rvalue=(*m_Rtet)(TetEdge4[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge5[j],PrismFace1[i],Rvalue);

                    //transposed values
                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge0[j]);
                    RTmodprism->SetValue(PrismFace1[i],PrismEdge0[j],RTvalue);
                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge3[j]);
                    RTmodprism->SetValue(PrismFace1[i],PrismEdge4[j],RTvalue);
                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge4[j]);
                    RTmodprism->SetValue(PrismFace1[i],PrismEdge5[j],RTvalue);
                }
            }
                
            // Face 3 edge 2 6 7
            for(i=0; i< PrismFace3.num_elements(); ++i)
            {
                for(j=0; j<PrismEdge2.num_elements(); ++j)
                {
                    Rvalue=(*m_Rtet)(TetEdge0[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge2[j],PrismFace3[i],Rvalue);
                    Rvalue=(*m_Rtet)(TetEdge4[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge6[j],PrismFace3[i],Rvalue);
                    Rvalue=(*m_Rtet)(TetEdge3[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge7[j],PrismFace3[i],Rvalue);

                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge0[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge2[j],RTvalue);
                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge4[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge6[j],RTvalue);
                    RTvalue=(*m_RTtet)(TetFace[i],TetEdge3[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge7[j],RTvalue);
                }
            }
        }
    
#endif

        DNekMatSharedPtr PreconditionerLowEnergy::
        ExtractLocMat(StdRegions::StdExpansionSharedPtr
                      &locExp)
        {
            LibUtilities::ShapeType eType=locExp->DetShapeType();
            int cnt, cnt1;
            int nverts = locExp->GetNverts();
            int nedges = locExp->GetNedges();
            int nfaces = locExp->GetNfaces();
            NekDouble val;
            NekDouble zero = 0.0;
            
            Array<OneD, Array<OneD, unsigned int> > emap = m_edgeMapMaxR[eType];
            Array<OneD, Array<OneD, unsigned int> > fmap = m_faceMapMaxR[eType];
            
            int nRows = locExp->NumBndryCoeffs();
            DNekMatSharedPtr newmat = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);
            
            
            // fill diagonal
            for(int i = 0; i < nRows; ++i)
            {
                val = 1.0;
                newmat->SetValue(i,i,val);
            }
            
            // fill vertex off diagonal
            for(int v = 0; v < nverts; ++v)
            {
                cnt = 0;
                
                for(int e = 0; e < nedges; ++e)
                {
                    int nEdgeInteriorCoeffs = locExp->GetEdgeNcoeffs(e) -2;
                    for(int i = 0; i < nEdgeInteriorCoeffs; ++i,++cnt)
                    {
                        val = (*m_maxRmat[eType])(v,emap[e][i]);
                        newmat->SetValue(v,nverts+cnt,val);
                    }
                }
                

                for(int f = 0; f < nfaces; ++f)
                {
                    int nFaceInteriorCoeffs = locExp->GetFaceIntNcoeffs(f);
                    for(int i = 0; i < nFaceInteriorCoeffs; ++i,++cnt)
                    {
                        val = (*m_maxRmat[eType])(v,fmap[f][i]);
                        newmat->SetValue(v,nverts+cnt,val);
                    }
                }
            }            
            
            
            // fill edges off diagonal
            int offset = nverts; 
            for(int e = 0; e < nedges; ++e)
            {
                int nEdgeInteriorCoeffs = locExp->GetEdgeNcoeffs(e) - 2;
                offset += nEdgeInteriorCoeffs; 
            }
            
            cnt1 = 0; 
            for(int e = 0; e < nedges; ++e)
            {
                cnt = 0; 
                int nEdgeInteriorCoeffs = locExp->GetEdgeNcoeffs(e) -2;
                
                for(int j = 0; j < nEdgeInteriorCoeffs; ++j, ++cnt1)
                {
                    for(int f = 0; f < nfaces; ++f)
                    {
                        int nFaceInteriorCoeffs =
                            locExp->GetFaceIntNcoeffs(f);
                        
                        for(int i = 0; i < nFaceInteriorCoeffs; ++i,++cnt)
                        {
                            val = (*m_maxRmat[eType])(emap[e][j],fmap[f][i]);
                            newmat->SetValue(nverts+cnt1,offset+cnt,val);
                        }
                    }
                }
            }            
            
            return newmat;
        }
    }
}




