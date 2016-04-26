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
            
            //Sets up reference element and builds transformation matrix
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

            DNekMatSharedPtr    pRS;
            DNekMatSharedPtr    pRSRT;

            //Transformation matrices
            DNekMat R;
            DNekMat RT;
            DNekMat RS;
            DNekMat RSRT;
            
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);

            Array<OneD, NekDouble> vertArray(nNonDirVerts,0.0);
            Array<OneD, long> VertBlockToUniversalMap(nNonDirVerts,-1);

            //maps for different element types
            map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transmatrixmap;
            map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transposedtransmatrixmap;

            //Transformation matrix
            transmatrixmap[LibUtilities::eTetrahedron]= m_Rtet;
            transmatrixmap[LibUtilities::ePrism]      = m_Rprism;
            transmatrixmap[LibUtilities::eHexahedron] = m_Rhex;

            //Transposed transformation matrix
            transposedtransmatrixmap[LibUtilities::eTetrahedron]= m_RTtet;
            transposedtransmatrixmap[LibUtilities::ePrism]      = m_RTprism;
            transposedtransmatrixmap[LibUtilities::eHexahedron] = m_RThex;

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
            
            const Array<OneD,const unsigned int>& nbdry_size
                    = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

            //Variants of R matrices required for low energy preconditioning
            m_RBlk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
            m_RTBlk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
            
            //Here we loop over the expansion and build the block low energy
            //preconditioner as well as the block versions of the transformation
            //matrices.
            for(cnt=n=0; n < n_exp; ++n)
            {
                eid = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(eid);
                nCoeffs=locExpansion->NumBndryCoeffs();
                LibUtilities::ShapeType eType=locExpansion->DetShapeType();

                //Get correct transformation matrix for element type
                R=(*(transmatrixmap[eType]));
                RT=(*(transposedtransmatrixmap[eType]));
                
                pRS = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nCoeffs, nCoeffs, zero, storage);
                RS = (*pRS);
                
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

                //Calculate S*trans(R)  (note R is already transposed)
                RS=R*S;

                //Calculate R*S*trans(R)
                RSRT=RS*RT;

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

                                // Get the face-face value from the low energy matrix (S2)
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
                
                //Here we build the block matrices for R and RT
                m_RBlk->SetBlock(n,n, transmatrixmap[eType]);
                m_RTBlk->SetBlock(n,n, transposedtransmatrixmap[eType]);
            }
            
            if(nNonDirVerts!=0)
            {
                //Exchange vertex data over different processes
                Gs::gs_data *tmp = Gs::Init(VertBlockToUniversalMap, m_comm);
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
            Gs::gs_data *tmp1 = Gs::Init(EdgeBlockToUniversalMap, m_comm);
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
            Gs::gs_data *tmp2 = Gs::Init(FaceBlockToUniversalMap, m_comm);
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

           int n, nel;
 
           const Array<OneD,const unsigned int>& nbdry_size
               = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

           int n_exp=expList->GetNumElmts();

           //maps for different element types
           map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transmatrixmap;
           map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transposedtransmatrixmap;
           map<LibUtilities::ShapeType,DNekScalMatSharedPtr> invtransmatrixmap;
           map<LibUtilities::ShapeType,DNekScalMatSharedPtr> invtransposedtransmatrixmap;

           //Transformation matrix map
           transmatrixmap[LibUtilities::eTetrahedron]=m_Rtet;
           transmatrixmap[LibUtilities::ePrism]=m_Rprism;
           transmatrixmap[LibUtilities::eHexahedron]=m_Rhex;

           //Transposed transformation matrix map
           transposedtransmatrixmap[LibUtilities::eTetrahedron]=m_RTtet;
           transposedtransmatrixmap[LibUtilities::ePrism]=m_RTprism;
           transposedtransmatrixmap[LibUtilities::eHexahedron]=m_RThex;

           //Inverse transfomation map
           invtransmatrixmap[LibUtilities::eTetrahedron]=m_Rinvtet;
           invtransmatrixmap[LibUtilities::ePrism]=m_Rinvprism;
           invtransmatrixmap[LibUtilities::eHexahedron]=m_Rinvhex;

           //Inverse transposed transformation map
           invtransposedtransmatrixmap[LibUtilities::eTetrahedron]=m_RTinvtet;
           invtransposedtransmatrixmap[LibUtilities::ePrism]=m_RTinvprism;
           invtransposedtransmatrixmap[LibUtilities::eHexahedron]=m_RTinvhex;

           MatrixStorage blkmatStorage = eDIAGONAL;
           
           //Variants of R matrices required for low energy preconditioning
           m_RBlk      = MemoryManager<DNekScalBlkMat>
               ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
           m_RTBlk      = MemoryManager<DNekScalBlkMat>
               ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
           m_InvRBlk      = MemoryManager<DNekScalBlkMat>
               ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
           m_InvRTBlk      = MemoryManager<DNekScalBlkMat>
               ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);

           for(n=0; n < n_exp; ++n)
           {
               nel = expList->GetOffset_Elmt_Id(n);
               
               locExpansion = expList->GetExp(nel);
               LibUtilities::ShapeType eType=locExpansion->DetShapeType();

               //Block R matrix
               m_RBlk->SetBlock(n,n, transmatrixmap[eType]);

               //Block RT matrix
               m_RTBlk->SetBlock(n,n, transposedtransmatrixmap[eType]);

               //Block inverse R matrix
               m_InvRBlk->SetBlock(n,n, invtransmatrixmap[eType]);

               //Block inverse RT matrix
               m_InvRTBlk->SetBlock(n,n, invtransposedtransmatrixmap[eType]);
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
            DNekScalBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();

            //Not actually needed but we should only work with the Global boundary dofs
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobBndDofs, pInOut.get(), 1, tmp.get(), 1);

            //Global boundary (with dirichlet values) to local boundary with multiplicity
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(), tmp.get(), m_map.get(), pLocal.get());

            //Multiply by the block transformation matrix
            F_LocBnd=R*F_LocBnd;

            //Assemble local boundary to global non-dirichlet Dofs
            m_locToGloMap->AssembleBnd(F_LocBnd,F_HomBnd, nDirBndDofs);
        }

        /**
         * \brief Transform the solution vector vector to low energy.
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
            DNekScalBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();

            // Allocated array of size number of global boundary dofs and copy
            // the input array to the tmp array offset by Dirichlet boundary
            // conditions.
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobHomBndDofs, pInput.get(), 1, tmp.get() + nDirBndDofs, 1);
            
            //Global boundary dofs (with zeroed dirichlet values) to local boundary dofs
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(), tmp.get(), m_map.get(), pLocal.get());

            //Multiply by the block transformation matrix
            F_LocBnd=R*F_LocBnd;

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

            //Block transposed transformation matrix
            DNekScalBlkMat &RT = *m_RTBlk;

            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,pInOut+nDirBndDofs,
                                              eWrapper);

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);

            //Global boundary (less dirichlet) to local boundary
            m_locToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd, nDirBndDofs);

            //Multiply by the block transposed transformation matrix
            V_LocBnd=RT*V_LocBnd;


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
            DNekScalBlkMat &invR = *m_InvRBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();

            // Allocated array of size number of global boundary dofs and copy
            // the input array to the tmp array offset by Dirichlet boundary
            // conditions.
            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);
            Vmath::Vcopy(nGlobHomBndDofs, pInput.get(), 1, tmp.get() + nDirBndDofs, 1);

            //Global boundary dofs (with zeroed dirichlet values) to local boundary dofs
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(), tmp.get(), m_map.get(), pLocal.get());

            //Multiply by block inverse transformation matrix
            F_LocBnd=invR*F_LocBnd;

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
            DNekScalBlkMat &invRT = *m_InvRTBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            m_map = m_locToGloMap->GetLocalToGlobalBndMap();

            m_locToGloMap->GlobalToLocalBnd(pInput,pLocal, nDirBndDofs);

            //Multiply by the block transposed transformation matrix
            F_LocBnd=invRT*F_LocBnd;

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
        v_TransformedSchurCompl(
            int offset, 
            const boost::shared_ptr<DNekScalMat > &loc_mat)
	{
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
         
            StdRegions::StdExpansionSharedPtr locExpansion;                
            locExpansion = expList->GetExp(offset);
            unsigned int nbnd=locExpansion->NumBndryCoeffs();

            //This is the SC elemental matrix in the orginal basis (S1)
            DNekScalMatSharedPtr pS1=loc_mat;

            //Transformation matrices 
            map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transmatrixmap;
            map<LibUtilities::ShapeType,DNekScalMatSharedPtr> transposedtransmatrixmap;
            transmatrixmap[LibUtilities::eTetrahedron]=m_Rtet;
            transmatrixmap[LibUtilities::ePrism]=m_Rprism;
            transmatrixmap[LibUtilities::eHexahedron]=m_Rhex;
            transposedtransmatrixmap[LibUtilities::eTetrahedron]=m_RTtet;
            transposedtransmatrixmap[LibUtilities::ePrism]=m_RTprism;
            transposedtransmatrixmap[LibUtilities::eHexahedron]=m_RThex;

            DNekScalMat &S1 = (*pS1);
            
            MatrixStorage storage = eFULL;
            
            DNekMatSharedPtr pS2 = MemoryManager<DNekMat>::AllocateSharedPtr(nbnd,nbnd,0.0,storage);
            DNekMatSharedPtr pRS1 = MemoryManager<DNekMat>::AllocateSharedPtr(nbnd,nbnd,0.0,storage);
            
            LibUtilities::ShapeType eType=
                (expList->GetExp(offset))->DetShapeType();
            
            //transformation matrices
            DNekScalMat &R = (*(transmatrixmap[eType]));
            DNekScalMat &RT = (*(transposedtransmatrixmap[eType]));
            
            //create low energy matrix
            DNekMat &RS1 = (*pRS1);
            DNekMat &S2 = (*pS2);
                
            //setup S2
            RS1=R*S1;
            S2=RS1*RT;

            DNekScalMatSharedPtr tmp_mat;
            tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(1.0, pS2);

	    return tmp_mat;
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
	 * \brief Sets up the reference elements needed by the preconditioner
	 *
         * Sets up reference elements which are used to preconditioning the
         * corresponding matrices. Currently we support tetrahedral, prismatic
         * and hexahedral elements
	 */
        void PreconditionerLowEnergy::SetUpReferenceElements()
        {
            int cnt,i,j;
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            StdRegions::StdExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(0);

            DNekScalBlkMatSharedPtr RtetBlk, RprismBlk;
            DNekScalBlkMatSharedPtr RTtetBlk, RTprismBlk;

            DNekScalMatSharedPtr Rprismoriginal;
            DNekScalMatSharedPtr RTprismoriginal;
            DNekMatSharedPtr Rtettmp, RTtettmp, Rhextmp, RThextmp, Rprismtmp, RTprismtmp ;

            /*
             * Set up a Tetrahral & prismatic element which comprises
             * equilateral triangles as all faces for the tet and the end faces
             * for the prism. Using these elements a new expansion is created
             * (which is the same as the expansion specified in the input
             * file).
             */
            SpatialDomains::TetGeomSharedPtr tetgeom=CreateRefTetGeom();
            SpatialDomains::PrismGeomSharedPtr prismgeom=CreateRefPrismGeom();
            SpatialDomains::HexGeomSharedPtr hexgeom=CreateRefHexGeom();

            //Expansion as specified in the input file - here we need to alter
            //this so we can read in different exapansions for different element
            //types
            int nummodes=locExpansion->GetBasisNumModes(0);

            //Bases for Tetrahedral element
            const LibUtilities::BasisKey TetBa(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey TetBb(
                LibUtilities::eModified_B, nummodes,
                LibUtilities::PointsKey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0));
            const LibUtilities::BasisKey TetBc(
                LibUtilities::eModified_C, nummodes,
                LibUtilities::PointsKey(nummodes,LibUtilities::eGaussRadauMAlpha2Beta0));

            //Create reference tetrahedral expansion
            LocalRegions::TetExpSharedPtr TetExp;

            TetExp = MemoryManager<LocalRegions::TetExp>
                ::AllocateSharedPtr(TetBa,TetBb,TetBc,
                                    tetgeom);

            //Bases for prismatic element
            const LibUtilities::BasisKey PrismBa(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey PrismBb(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey PrismBc(
                LibUtilities::eModified_B, nummodes,
                LibUtilities::PointsKey(nummodes,LibUtilities::eGaussRadauMAlpha1Beta0));

            //Create reference prismatic expansion
            LocalRegions::PrismExpSharedPtr PrismExp;

            PrismExp = MemoryManager<LocalRegions::PrismExp>
                ::AllocateSharedPtr(PrismBa,PrismBb,PrismBc,
                                    prismgeom);

            //Bases for prismatic element
            const LibUtilities::BasisKey HexBa(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey HexBb(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey HexBc(
                LibUtilities::eModified_A, nummodes,
                LibUtilities::PointsKey(nummodes+1,LibUtilities::eGaussLobattoLegendre));
            
            //Create reference prismatic expansion
            LocalRegions::HexExpSharedPtr HexExp;
            
            HexExp = MemoryManager<LocalRegions::HexExp>
                ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                    hexgeom);
            

            // retrieve variable coefficient
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                StdRegions::VarCoeffMap::const_iterator x;
                cnt = expList->GetPhys_Offset(0);
                for (x = m_linSysKey.GetVarCoeffs().begin(); 
                     x != m_linSysKey.GetVarCoeffs().end(); ++x)
                {
                    vVarCoeffMap[x->first] = x->second + cnt;
                }
            }

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
                 *TetExp, m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Matrix keys for tetrahedral transposed transformation matrix
            LocalRegions::MatrixKey TetRT
                (PreconRT, LibUtilities::eTetrahedron,
                 *TetExp,  m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Matrix keys for prismatic element transformation matrix
            LocalRegions::MatrixKey PrismR
                (PreconR,   LibUtilities::ePrism,
                 *PrismExp, m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Matrix keys for prismatic element transposed transformation matrix
            LocalRegions::MatrixKey PrismRT
                (PreconRT,  LibUtilities::ePrism,
                 *PrismExp, m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Matrix keys for hexahedral element transformation matrix
            LocalRegions::MatrixKey HexR
                (PreconR, LibUtilities::eHexahedron,
                 *HexExp, m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Matrix keys for hexahedral element transposed transformation
            //matrix
            LocalRegions::MatrixKey HexRT
                (PreconRT, LibUtilities::eHexahedron,
                 *HexExp,  m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            /*
             * Create transformation matrices for the tetrahedral element
             */

            //Get tetrahedral transformation matrix
            m_Rtet = TetExp->GetLocMatrix(TetR);

            //Get tetrahedral transposed transformation matrix
            m_RTtet = TetExp->GetLocMatrix(TetRT);

            // Using the transformation matrix and the inverse transformation
            // matrix create the inverse matrices
            Rtettmp=TetExp->BuildInverseTransformationMatrix(m_Rtet);

            //Inverse transposed transformation matrix
            RTtettmp=TetExp->BuildInverseTransformationMatrix(m_Rtet);
            RTtettmp->Transpose();

            m_Rinvtet = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,Rtettmp);
            m_RTinvtet = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RTtettmp);

            /*
             * Create transformation matrices for the hexahedral element
             */

            //Get hexahedral transformation matrix
            m_Rhex = HexExp->GetLocMatrix(HexR);
            //Get hexahedral transposed transformation matrix
            m_RThex = HexExp->GetLocMatrix(HexRT);

            // Using the transformation matrix and the inverse transformation
            // matrix create the inverse matrices
            Rhextmp=HexExp->BuildInverseTransformationMatrix(m_Rhex);
            //Inverse transposed transformation matrix
            RThextmp=HexExp->BuildInverseTransformationMatrix(m_Rhex);
            RThextmp->Transpose();

            m_Rinvhex = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,Rhextmp);
            m_RTinvhex = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RThextmp);

            /*
             * Create transformation matrices for the prismatic element
             */

            //Get prism transformation matrix
            Rprismoriginal = PrismExp->GetLocMatrix(PrismR);
            //Get prism transposed transformation matrix
            RTprismoriginal = PrismExp->GetLocMatrix(PrismRT);

            unsigned int  nRows=Rprismoriginal->GetRows();
            NekDouble zero=0.0;
            DNekMatSharedPtr Rtmpprism = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);
            DNekMatSharedPtr RTtmpprism = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);
            NekDouble Rvalue, RTvalue;

            //Copy values from the prism transformation matrix
            for(i=0; i<nRows; ++i)
            {
                for(j=0; j<nRows; ++j)
                {
                    Rvalue=(*Rprismoriginal)(i,j);
                    RTvalue=(*RTprismoriginal)(i,j);
                    Rtmpprism->SetValue(i,j,Rvalue);
                    RTtmpprism->SetValue(i,j,RTvalue);
                }
            }

            //Replace triangular faces and edges of the prims transformation
            //matrix with the corresponding values of the tetrahedral
            //transformation matrix.
            ModifyPrismTransformationMatrix(TetExp,PrismExp,Rtmpprism,RTtmpprism);

            m_Rprism = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,Rtmpprism);
            
            m_RTprism = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RTtmpprism);

            //Inverse transformation matrix
            Rprismtmp=PrismExp->BuildInverseTransformationMatrix(m_Rprism);

            //Inverse transposed transformation matrix
            RTprismtmp=PrismExp->BuildInverseTransformationMatrix(m_Rprism);
            RTprismtmp->Transpose();

            m_Rinvprism = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,Rprismtmp);

            m_RTinvprism = MemoryManager<DNekScalMat>
                ::AllocateSharedPtr(1.0,RTprismtmp);
        }

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


            //These are the edge mode locations of R which need to be replaced
            //in the prism element
            Array<OneD, unsigned int> TetEdge0=TetExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> TetEdge1=TetExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> TetEdge2=TetExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> TetEdge3=TetExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> TetEdge4=TetExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> TetEdge5=TetExp->GetEdgeInverseBoundaryMap(5);

            //These are the face mode locations of R which need to be replaced
            //in the prism element
            Array<OneD, unsigned int> TetFace=TetExp->GetFaceInverseBoundaryMap(1);

            //Prism vertex modes
            int PrismVertex0=PrismExp->GetVertexMap(0);
            int PrismVertex1=PrismExp->GetVertexMap(1);
            int PrismVertex2=PrismExp->GetVertexMap(2);
            int PrismVertex3=PrismExp->GetVertexMap(3);
            int PrismVertex4=PrismExp->GetVertexMap(4);
            int PrismVertex5=PrismExp->GetVertexMap(5);

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
        
    }
}






