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
#include <LocalRegions/MatrixKey.h>
#include <math.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */
        string PreconditionerLowEnergy::className1
                = GetPreconFactory().RegisterCreatorFunction(
                    "LowEnergy",
                    PreconditionerLowEnergy::create,
                    "LowEnergy Preconditioning");

        string PreconditionerLowEnergy::className2
                = GetPreconFactory().RegisterCreatorFunction(
                    "Block",
                    PreconditionerLowEnergy::create,
                    "Block Preconditioning");

        string PreconditionerLowEnergy::className3
                = GetPreconFactory().RegisterCreatorFunction(
                    "InverseLinear",
                    PreconditionerLowEnergy::create,
                    "Linear space inverse Preconditioning");

 
       /**
         * @class Preconditioner
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */

         PreconditionerLowEnergy::PreconditionerLowEnergy(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap)
           : Preconditioner(plinsys, pLocToGloMap),
	   m_linsys(plinsys),
           m_locToGloMap(pLocToGloMap),
           m_preconType(pLocToGloMap->GetPreconType())
         {
	 }

        void PreconditionerLowEnergy::v_InitObject()
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(m_preconType)
            {
	    case MultiRegions::eLowEnergy:
	        {
                    if(solvertype != eIterativeStaticCond)
		    {
                        ASSERTL0(0,"Solver type not valid");
		    }
                    SetUpReferenceElements();
                    LowEnergyPreconditioner();
		}
		break;
            case MultiRegions::eBlock:
                {
                    BlockPreconditioner();
		}
		break;
            case MultiRegions::eInverseLinear:
                {
                    if (solvertype == eIterativeFull)
                    {
                        InverseLinearSpacePreconditioner();
                    }
                    else if(solvertype == eIterativeStaticCond)
                    {
                        StaticCondInverseLinearSpacePreconditioner();
                    }
                    else
                    {
                        ASSERTL0(0,"Unsupported solver type");
                    }
                }
                break;
            default:
                ASSERTL0(0,"Unknown preconditioner");
                break;
            }
	}

        /**
         * \brief Inverse of the linear space
	 *
	 * Extracts the linear space and inverts it.
	 *
         */         
        void PreconditionerLowEnergy::InverseLinearSpacePreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector 
                &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;
            
            int localVertId, vMap1, vMap2, nVerts, localCoeff, n, v, m;
            int sign1, sign2, gid1, gid2, i, j;
            int loc_rows, globalrow, globalcol, cnt;
            int  globalLocation, nIntEdgeFace;

            NekDouble globalMatrixValue, MatrixValue, value;
            NekDouble TempValue;
            NekDouble zero=0.0;

            int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
            int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();

            Array<OneD, NekDouble> vOutput(nGlobal,0.0);           
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_S;
            m_S = MemoryManager<DNekMat>::AllocateSharedPtr
                (nNonDirVerts, nNonDirVerts, zero,  storage);
            DNekMat &S = (*m_S);
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr
                (nInt, nInt, zero, storage);
            DNekMat &M = (*m_preconditioner);
 
            DNekScalMatSharedPtr loc_mat;
            int n_exp=expList->GetNumElmts();

            for(cnt=n=0; n < n_exp; ++n)
            {
                //element matrix
                loc_mat = 
                    (m_linsys.lock())->GetBlock(expList->GetOffset_Elmt_Id(n));
                loc_rows = loc_mat->GetRows();
                
                //element expansion
                locExpansion = 
                    boost::dynamic_pointer_cast<StdRegions::StdExpansion>(
                        locExpVector[expList->GetOffset_Elmt_Id(n)]);

                //Get number of vertices
                nVerts=locExpansion->GetGeom()->GetNumVerts();

                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    //Get vertex map
                    vMap1 = locExpansion->GetVertexMap(v);

                    globalrow = m_locToGloMap->
                        GetLocalToGlobalMap(cnt+vMap1)-nDir;
                    
                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2 = locExpansion->GetVertexMap(m);

                            //global matrix location (with offset due to
                            //dirichlet values)
                            globalcol = m_locToGloMap->
                                GetLocalToGlobalMap(cnt+vMap2)-nDir;

                            if(globalcol>=0)
                            {
                                
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalSign(cnt + vMap2);

                                //Global matrix value
                                globalMatrixValue = 
                                    S.GetValue(globalrow,globalcol)
                                    + sign1*sign2*(*loc_mat)(vMap1,vMap2);
                        
                                //build matrix containing the linear finite
                                //element space
                                S.SetValue
                                    (globalrow,globalcol,globalMatrixValue);
                            }
                        }
                    }
                }
                   //move counter down length of loc_rows
                cnt   += loc_rows;
            }
            
            for(n = cnt = 0; n < n_exp; ++n)
            {
                loc_mat = (m_linsys.lock())->
                    GetBlock(expList->GetOffset_Elmt_Id(n));
                loc_rows = loc_mat->GetRows();

                for(i = 0; i < loc_rows; ++i)
                {
                    gid1 = m_locToGloMap->
                        GetLocalToGlobalMap(cnt + i) - nDir-nNonDirVerts;
                    sign1 =  m_locToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_rows; ++j)
                        {
                            gid2 = m_locToGloMap->GetLocalToGlobalMap(cnt + j)
                                 - nDir-nNonDirVerts;
                            sign2 = m_locToGloMap->
                                GetLocalToGlobalSign(cnt + j);
                            if(gid2 == gid1)
                            {
                                value = vOutput[gid1 + nDir + nNonDirVerts]
                                      + sign1*sign2*(*loc_mat)(i,j);
                                vOutput[gid1 + nDir + nNonDirVerts] = value;
                            }
                        }
                    }
                }
                cnt   += loc_rows;
            }

            // Assemble diagonal contributions across processes
            m_locToGloMap->UniversalAssemble(vOutput);

            //Invert vertex space
            if(nNonDirVerts != 0)
            {
                S.Invert();
            }

            //Extract values
            for(int i = 0; i < S.GetRows(); ++i)
            {
                for(int j = 0; j < S.GetColumns(); ++j)
                {
                    MatrixValue=S.GetValue(i,j);
                    M.SetValue(i,j,MatrixValue);
                }
            }

            // Populate preconditioner matrix
            for (unsigned int i = nNonDirVerts; i < M.GetRows(); ++i)
            {
                M.SetValue(i,i,1.0/vOutput[nDir + i]);
            }

        }


        /**
	 * \brief Extracts the entries from a statically condensed matrix
	 * corresponding to the vertex modes.
	 *
	 * This function extracts a static condensed matrix from the nth 
	 * expansion and returns a matrix of the same dimensions containing
	 * only the vertex mode contributions i.e the linear finite element
	 * space.
	 */         
        void PreconditionerLowEnergy::
        StaticCondInverseLinearSpacePreconditioner()
	{
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector 
                &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;
            
            int localVertId, vMap1, vMap2, nVerts, nEdges, nFaces;
            int localCoeff, v, m, n, rows;
            int sign1, sign2;
            int offset, globalrow, globalcol, cnt, cnt1;
            int globalLocation, nDirVertOffset;

            NekDouble globalMatrixValue;
            NekDouble MatrixValue;
            NekDouble localMatrixValue;
            NekDouble zero=0.0;
            DNekMatSharedPtr m_invS;

            int nGlobalBnd    = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd       = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nDirBndLocal  = m_locToGloMap->GetNumLocalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();
            
            //Get Rows of statically condensed matrix
            int gRow = nGlobalBnd - nDirBnd;
            
            //Allocate preconditioner matrix
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_S;
            m_S = 
                MemoryManager<DNekMat>::AllocateSharedPtr
                (nNonDirVerts, nNonDirVerts, zero,  storage);
            DNekMat &S = (*m_S);

            //element expansion
            locExpansion = 
                boost::dynamic_pointer_cast<StdRegions::StdExpansion>(
                    locExpVector[expList->GetOffset_Elmt_Id(0)]);

            //Get total number of vertices
            nVerts=locExpansion->GetGeom()->GetNumVerts();

            m_preconditioner = 
                MemoryManager<DNekMat>::AllocateSharedPtr(
                    nGlobalBnd-nDirBnd, nGlobalBnd-nDirBnd, zero,  storage);
            DNekMat &M = (*m_preconditioner);

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

                bnd_mat=loc_mat->GetBlock(0,0);
		
                //offset by number of rows
                offset = bnd_mat->GetRows();
		
                //loop over vertices of the element and return the vertex map
                //for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    //Get vertex map
                    vMap1 = locExpansion->GetVertexMap(v);
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2 = locExpansion->GetVertexMap(m);

                            //global matrix location (without offset due to
                            //dirichlet values)
                            globalcol = m_locToGloMap->
                                GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol >= 0)
                            {
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->
                                    GetLocalToGlobalBndSign(cnt + vMap2);

                                //Global matrix value
                                globalMatrixValue = 
                                    S.GetValue(globalrow,globalcol)
                                    + sign1*sign2*(*bnd_mat)(vMap1,vMap2);

                                //build matrix containing the linear finite
                                //element space
                                S.SetValue
                                    (globalrow,globalcol,globalMatrixValue);
                            }
                        }
                    }
                }
                cnt   += offset;
            }

            //Invert linear finite element space
            if(nNonDirVerts != 0)
            {
                S.Invert();
            }

            //Extract values
            for(int i = 0; i < S.GetRows(); ++i)
            {
                for(int j = 0; j < S.GetColumns(); ++j)
                {
                    MatrixValue=S.GetValue(i,j);
                    M.SetValue(i,j,MatrixValue);
                }
            }

            Array<OneD, NekDouble> diagonals = 
                AssembleStaticCondGlobalDiagonals();

            // Populate preconditioner matrix
            for (unsigned int i = nNonDirVerts; i < M.GetRows(); ++i)
            {
                  M.SetValue(i,i,1.0/diagonals[i]);
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
            
            //boost::shared_ptr<SpatialDomains::VertexComponent> verts[6];
            SpatialDomains::VertexComponentSharedPtr verts[6];
            for(int i=0; i < nVerts; ++i)
            {
                verts[i] =  MemoryManager<SpatialDomains::VertexComponent>::AllocateSharedPtr
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
                SpatialDomains::VertexComponentSharedPtr vertsArray[2];
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
            /////////////////////////
            // Set up Tetrahedron  //
            /////////////////////////

	    int i,j;
	    const int three=3;
            const int nVerts = 4;
            const double point[][3] = {
                {-1,-1/sqrt(double(3)),-1/sqrt(double(6))},
                {1,-1/sqrt(double(3)),-1/sqrt(double(6))},
                {0,2/sqrt(double(3)),-1/sqrt(double(6))},
                {0,0,3/sqrt(double(6))}};
            
            boost::shared_ptr<SpatialDomains::VertexComponent> verts[4];
	    for(i=0; i < nVerts; ++i)
	    {
	        verts[i] =  
                    MemoryManager<SpatialDomains::VertexComponent>::
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
                boost::shared_ptr<SpatialDomains::VertexComponent> 
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
	 * \brief Sets up the reference elements needed by the preconditioner
	 *
         * Sets up reference elements which are used to preconditioning the
         * corresponding matrices. Currently we support tetrahedral, prismatic
         * and hexahedral elements
	 */         

        void PreconditionerLowEnergy::SetUpReferenceElements()
        {
            int cnt;
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = 
                *(expList->GetExp());
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            StdRegions::StdExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(0);

            DNekScalBlkMatSharedPtr r_mat;
            DNekScalBlkMatSharedPtr rt_mat;

            /*
             * Set up a Tetrahral & prismatic element which comprises
             * equilateral triangles as all faces for the tet and the end faces
             * for the prism. Using these elements a new expansion is created
             * (which is the same as the expansion specified in the input
             * file).*/

            SpatialDomains::TetGeomSharedPtr tetgeom=CreateRefTetGeom();
            SpatialDomains::PrismGeomSharedPtr prismgeom=CreateRefPrismGeom();

            //Expansion as specified in the input file
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

            //Matrix keys for tetrahedral element transformation matrices
            LocalRegions::MatrixKey TetR
                (StdRegions::ePreconR,
                 StdRegions::eTetrahedron,
                 *TetExp,
                 m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);
                
            LocalRegions::MatrixKey TetRT
                (StdRegions::ePreconRT,
                 StdRegions::eTetrahedron,
                 *TetExp,
                 m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);

            //Arrays to store the transformation matrices for each element type;
            m_transformationMatrix = Array<OneD,DNekScalMatSharedPtr>(2);
            m_transposedTransformationMatrix = Array<OneD,DNekScalMatSharedPtr>(2);

            int elmtType=0;
            //Get a LocalRegions static condensed matrix
            r_mat = TetExp->GetLocStaticCondMatrix(TetR);
            rt_mat = TetExp->GetLocStaticCondMatrix(TetRT);
            m_transformationMatrix[elmtType]=r_mat->GetBlock(0,0);
            m_transposedTransformationMatrix[elmtType]=rt_mat->GetBlock(0,0);

            //Matrix keys for Prism element transformation matrices
            LocalRegions::MatrixKey PrismR
                (StdRegions::ePreconR,
                 StdRegions::ePrism,
                 *PrismExp,
                 m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);
                
            LocalRegions::MatrixKey PrismRT
                (StdRegions::ePreconRT,
                 StdRegions::ePrism,
                 *PrismExp,
                 m_linSysKey.GetConstFactors(),
                 vVarCoeffMap);
            elmtType++;
            //Get a LocalRegions static condensed matrix
            r_mat = PrismExp->GetLocStaticCondMatrix(PrismR);
            rt_mat = PrismExp->GetLocStaticCondMatrix(PrismRT);
            m_transformationMatrix[elmtType]=r_mat->GetBlock(0,0);
            m_transposedTransformationMatrix[elmtType]=rt_mat->GetBlock(0,0);
        }


       /**
	 * \brief Construct the low energy preconditioner from
	 * \f$\mathbf{S}_{2}\f$
	 *
	 *\f[\mathbf{M}^{-1}=\left[\begin{array}{ccc}
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
        void PreconditionerLowEnergy::LowEnergyPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = 
                *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;

            int nRow, i, j, nmodes, ntotaledgemodes, ntotalfacemodes, nel;
            int nVerts, nEdges,nFaces; 
            int eid, fid, eid2, fid2, n, cnt, nedgemodes, nfacemodes;
            int nEdgeCoeffs, nFaceCoeffs;
            NekDouble zero = 0.0;
            NekDouble MatrixValue;

            int vMap1, vMap2, sign1, sign2, gid1, gid2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol, bnd_rows, nCoeffs;
            NekDouble globalMatrixValue, globalRValue;

            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;


            DNekMatSharedPtr    m_RS;
            DNekMatSharedPtr    m_RSRT;

            DNekMat R;
            DNekMat RT;
            DNekMat RS;
            DNekMat RSRT;
            
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

            // set up block matrix system
            int nblks=3;
            Array<OneD,unsigned int> exp_size(nblks);

            exp_size[0]=nNonDirVerts;
            exp_size[1]=nNonDirEdges;
            exp_size[2]=nNonDirFaces;

            MatrixStorage blkmatStorage = eDIAGONAL;
            GloBlkMat = MemoryManager<DNekScalBlkMat>::
                AllocateSharedPtr(exp_size,exp_size,blkmatStorage);

	    //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);
            DNekMatSharedPtr EdgeBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirEdges,nNonDirEdges,zero,storage);
            DNekMatSharedPtr FaceBlk = MemoryManager<DNekMat>::
                AllocateSharedPtr(nNonDirFaces,nNonDirFaces,zero,storage);

            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transmatrixmap;
            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transposedtransmatrixmap;
            transmatrixmap[StdRegions::eTetrahedron]=m_transformationMatrix[0];
            transmatrixmap[StdRegions::ePrism]=m_transformationMatrix[1];
            transposedtransmatrixmap[StdRegions::eTetrahedron]=m_transposedTransformationMatrix[0];
            transposedtransmatrixmap[StdRegions::ePrism]=m_transposedTransformationMatrix[1];
            
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {

                nel = expList->GetOffset_Elmt_Id(n);
                
                locExpansion = expList->GetExp(nel);
                nCoeffs=locExpansion->NumBndryCoeffs();
                StdRegions::ExpansionType eType=
                    locExpansion->DetExpansionType();
             
                R=(*(transmatrixmap[eType]));
                RT=(*(transposedtransmatrixmap[eType]));
 
                m_RS = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nCoeffs, nCoeffs, zero, storage);
                RS = (*m_RS);
                m_RSRT = MemoryManager<DNekMat>::AllocateSharedPtr
                    (nCoeffs, nCoeffs, zero, storage);
                RSRT = (*m_RSRT);
                
                nVerts=locExpansion->GetGeom()->GetNumVerts();
                nEdges=locExpansion->GetGeom()->GetNumEdges();
                nFaces=locExpansion->GetGeom()->GetNumFaces();
                Array<OneD, int> 
                    vertModeLocation(nVerts);
                Array<OneD, Array<OneD, unsigned int> > 
                    edgeModeLocation(nEdges);
                Array<OneD, Array<OneD, unsigned int> > 
                    faceModeLocation(nFaces);

                //Change this to array of arrays
                locExpansion->GetModeMappings
                    (vertModeLocation,edgeModeLocation,faceModeLocation);

                
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
                    vMap1=vertModeLocation[v];
                    
                    //Get vertex map
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;
                    
                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=vertModeLocation[m];
                            
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
                                    + sign1*sign2*RSRT(vMap1,vMap2);

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
                    nedgemodes=edgeModeLocation[eid].num_elements();
            
                    for (v=0; v<nedgemodes; ++v)
                    {
                        
                        eMap1=edgeModeLocation[eid][v];

                        globalrow = m_locToGloMap->
                            GetLocalToGlobalBndMap(cnt+eMap1)
                            -nDirBnd-nNonDirVerts;
                        
                        if(globalrow >= 0)
                        {
                            for (m=0; m<nedgemodes; ++m)
                            {
                                eMap2=edgeModeLocation[eid][m];
                                
                                //global matrix location (without offset due to
                                //dirichlet values)
                                globalcol = m_locToGloMap->
                                    GetLocalToGlobalBndMap(cnt+eMap2)
                                    -nDirBnd-nNonDirVerts;
                                
                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + eMap1);
                                    sign2 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + eMap2);

                                    globalMatrixValue = EdgeBlk->
                                        GetValue(globalrow,globalcol)
                                        + sign1*sign2*RSRT(eMap1,eMap2);

                                    //build matrix containing the linear finite
                                    //element space
                                    EdgeBlk->SetValue
                                        (globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
                    }
                }
                
                 //loop over faces of the element and return the face map
                for (fid=0; fid<nFaces; ++fid)
                {
                    nfacemodes=faceModeLocation[fid].num_elements();
            
                    for (v=0; v<nfacemodes; ++v)
                    {
                        fMap1=faceModeLocation[fid][v];
                        
                        globalrow = m_locToGloMap->
                            GetLocalToGlobalBndMap(cnt+fMap1)
                            -nDirBnd-nNonDirVerts-nNonDirEdges;
                        
                        if(globalrow >= 0)
                        {
                            for (m=0; m<nfacemodes; ++m)
                            {
                                fMap2=faceModeLocation[fid][m];

                                //global matrix location (without offset due to
                                //dirichlet values)
                                globalcol = m_locToGloMap->
                                    GetLocalToGlobalBndMap(cnt+fMap2)
                                    -nDirBnd-nNonDirVerts-nNonDirEdges;

                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + fMap1);
                                    sign2 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + fMap2);

                                    globalMatrixValue = FaceBlk->
                                        GetValue(globalrow,globalcol)
                                        + sign1*sign2*RSRT(fMap1,fMap2);
                                    
                                    //build matrix containing the linear finite
                                    //element space
                                    FaceBlk->SetValue
                                        (globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
                    }
                }
                cnt+=offset;
            }

            if (nNonDirVerts != 0)
            {
                VertBlk->Invert();
            }

            if (nNonDirEdges != 0)
            {
                EdgeBlk->Invert();
            }

            if (nNonDirFaces != 0)
            {
                FaceBlk->Invert();
            }

            DNekScalMatSharedPtr     Blktmp;
            NekDouble                one = 1.0;

            GloBlkMat->SetBlock(0,0,Blktmp = 
                                MemoryManager<DNekScalMat>::AllocateSharedPtr
                                (one,VertBlk));
            GloBlkMat->SetBlock(1,1,Blktmp = 
                                MemoryManager<DNekScalMat>::AllocateSharedPtr
                                (one,EdgeBlk));
            GloBlkMat->SetBlock(2,2,Blktmp = 
                                MemoryManager<DNekScalMat>::AllocateSharedPtr
                                (one,FaceBlk));
        }

        /**
	 * \brief Construct the Block preconditioner from \f$\mathbf{S}_{1}\f$
	 *
	 * \f[\mathbf{M}^{-1}=\left[\begin{array}{ccc}
	 *  Diag[(\mathbf{S_{1}})_{vv}] & & \ \ & (\mathbf{S}_{1})_{eb} & \\ & &
	 *  (\mathbf{S}_{1})_{fb} \end{array}\right]\f]
	 *
	 * where \f$\mathbf{R}\f$ is the transformation matrix and
	 * \f$\mathbf{S}_{1}\f$ is the Schur complement of each element.
	 *
	 */
        void PreconditionerLowEnergy::BlockPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector 
                &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;
            
            int nRow, i, j, nmodes, ntotaledgemodes, ntotalfacemodes;
            int nVerts, nEdges,nFaces, eid, fid, eid2, fid2, n, cnt;
            int  nedgemodes, nfacemodes;
            int nEdgeCoeffs, nFaceCoeffs;
            int vMap1, vMap2, sign1, sign2, gid1, gid2;
            int m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;

            NekDouble zero = 0.0;
            NekDouble MatrixValue;
            NekDouble globalMatrixValue, globalRValue;

            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            nRow=vExp->NumBndryCoeffs();

            nVerts=vExp->GetGeom()->GetNumVerts();
            nEdges=vExp->GetGeom()->GetNumEdges();
            nFaces=vExp->GetGeom()->GetNumFaces();

            int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nLocal = m_locToGloMap->GetNumLocalBndCoeffs();
            int nDirBnd    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDir = nGlobal-nDirBnd;
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();
            int nNonDirEdges  = m_locToGloMap->GetNumNonDirEdgeModes();
            int nNonDirFaces  = m_locToGloMap->GetNumNonDirFaceModes();

            // set up block matrix system
            int nblks=3;
            Array<OneD,unsigned int> exp_size(nblks);

            exp_size[0]=nNonDirVerts;
            exp_size[1]=nNonDirEdges;
            exp_size[2]=nNonDirFaces;

            MatrixStorage blkmatStorage = eDIAGONAL;
            GloBlkMat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr
                (exp_size,exp_size,blkmatStorage);

            //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr 
                VertBlk = MemoryManager<DNekMat>::AllocateSharedPtr
                (nNonDirVerts,nNonDirVerts,zero,vertstorage);
            DNekMatSharedPtr 
                EdgeBlk = MemoryManager<DNekMat>::AllocateSharedPtr
                (nNonDirEdges,nNonDirEdges,zero,storage);
            DNekMatSharedPtr 
                FaceBlk = MemoryManager<DNekMat>::AllocateSharedPtr
                (nNonDirFaces,nNonDirFaces,zero,storage);

            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
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
                    vMap1=vertModeLocation[v];
                    
                    //Get vertex map
                    globalrow = m_locToGloMap->
                        GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=vertModeLocation[m];

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
                    nedgemodes=edgeModeLocation[eid].num_elements();
            
                    for (v=0; v<nedgemodes; ++v)
                    {
                        
                        eMap1=edgeModeLocation[eid][v];
                        
                        globalrow = m_locToGloMap->
                            GetLocalToGlobalBndMap(cnt+eMap1)
                            -nDirBnd-nNonDirVerts;

                        if(globalrow >= 0)
                        {
                            for (m=0; m<nedgemodes; ++m)
                            {
                                eMap2=edgeModeLocation[eid][m];

                                //global matrix location (without offset due to
                                //dirichlet values)
                                globalcol = m_locToGloMap->
                                    GetLocalToGlobalBndMap(cnt+eMap2)
                                    -nDirBnd-nNonDirVerts;

                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + eMap1);
                                    sign2 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + eMap2);

                                    globalMatrixValue = 
                                        EdgeBlk->GetValue(globalrow,globalcol)
                                        + sign1*sign2*S(eMap1,eMap2);
		
                                    //build matrix containing the linear finite
                                    //element space
                                    EdgeBlk->SetValue
                                        (globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
                    }
                }
                
                 //loop over faces of the element and return the face map
                for (fid=0; fid<nFaces; ++fid)
                {
                    nfacemodes=faceModeLocation[fid].num_elements();
            
                    for (v=0; v<nfacemodes; ++v)
                    {
                        fMap1=faceModeLocation[fid][v];

                        globalrow = m_locToGloMap->
                            GetLocalToGlobalBndMap(cnt+fMap1)
                            -nDirBnd-nNonDirVerts-nNonDirEdges;
                        
                        if(globalrow >= 0)
                        {
                            for (m=0; m<nfacemodes; ++m)
                            {
                                fMap2=faceModeLocation[fid][m];
                                
                                //global matrix location (without offset due to
                                //dirichlet values)
                                globalcol = m_locToGloMap->
                                    GetLocalToGlobalBndMap(cnt+fMap2)
                                    -nDirBnd-nNonDirVerts-nNonDirEdges;

                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + fMap1);
                                    sign2 = m_locToGloMap->
                                        GetLocalToGlobalBndSign(cnt + fMap2);

                                    globalMatrixValue = 
                                        FaceBlk->GetValue(globalrow,globalcol)
                                        + sign1*sign2*S(fMap1,fMap2);

                                    //build matrix containing the linear finite
                                    //element space
                                    FaceBlk->SetValue
                                        (globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
                    }
                }
                cnt+=offset;
	    }

            if (nNonDirVerts != 0)
            {
                VertBlk->Invert();
            }

            if (nNonDirEdges != 0)
            {
                EdgeBlk->Invert();
            }

            if (nNonDirFaces != 0)
            {
                FaceBlk->Invert();
            }

            DNekScalMatSharedPtr     Blktmp;
            NekDouble                one = 1.0;
            
            GloBlkMat->SetBlock
                (0,0,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr
                 (one,VertBlk));
            GloBlkMat->SetBlock
                (1,1,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr
                 (one,EdgeBlk));
            GloBlkMat->SetBlock
                (2,2,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr
                 (one,FaceBlk));
        }
        
        
        /**
         *
         */
        void PreconditionerLowEnergy::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(m_preconType)
            {
            case MultiRegions::eInverseLinear:
                 {
                    if (solvertype == eIterativeFull)
                    {
                        int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                        int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                        int nNonDir = nGlobal-nDir;
                        DNekMat &M = (*m_preconditioner);

                        NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                        NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                        z = M * r;
                    }
                    else if(solvertype == eIterativeStaticCond)
                    {
                        int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                        int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                        int nNonDir = nGlobal-nDir;
                        DNekMat &M = (*m_preconditioner);

                        NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                        NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                        z = M * r;
                    }
                    else
                    {
                        ASSERTL0(0,"Unsupported solver type");
                    }
                }
                break;
            case MultiRegions::eLowEnergy:
            case MultiRegions::eBlock:
                 {
                    if (solvertype == eIterativeFull)
                    {
                        int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
                        int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                        int nNonDir = nGlobal-nDir;
                        DNekScalBlkMat &M = (*GloBlkMat);

                        NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                        NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                        z = M * r;
                    }
                    else if(solvertype == eIterativeStaticCond)
                    {
                        int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                        int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                        int nNonDir = nGlobal-nDir;
                        DNekScalBlkMat &M = (*GloBlkMat);

                        NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                        NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                        z = M * r;
                    }
                    else
                    {
                        ASSERTL0(0,"Unsupported solver type");
                    }
                }
                break;
            default:
            ASSERTL0(0,"Unknown preconditioner");
            break;
	    }
	}

        /**
         * \brief Get the tetrahedral transformation matrix \f$\mathbf{R}\f$
         */
        const Array<OneD,const DNekScalMatSharedPtr>& PreconditionerLowEnergy::
        v_GetTransformationMatrix() const
	{
	    return m_transformationMatrix;
	}

        /**
         * \brief Get the transposed transformation matrix \f$\mathbf{R}^{T}\f$
         */
        const Array<OneD,const DNekScalMatSharedPtr>& PreconditionerLowEnergy::
        v_GetTransposedTransformationMatrix() const
	{
	    return m_transposedTransformationMatrix;
	}
    }
}






