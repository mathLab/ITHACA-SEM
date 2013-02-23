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
                    "InverseLinear",
                    PreconditionerLowEnergy::create,
                    "Linear space inverse Preconditioning");

 
       /**
         * @class PreconditionerLowEnergy
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */
        
        PreconditionerLowEnergy::PreconditionerLowEnergy(
            const boost::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap),
              m_linsys(plinsys),
              m_preconType(pLocToGloMap->GetPreconType()),
              m_locToGloMap(pLocToGloMap)
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

            CreateMultiplicityMap();

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

            int vMap1, vMap2, nVerts, n, v, m;
            int sign1, sign2, gid1, gid2, i, j;
            int loc_rows, globalrow, globalcol, cnt;

            NekDouble globalMatrixValue, MatrixValue, value;
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
            
            int vMap1, vMap2, nVerts;
            int v, m, n;
            int sign1, sign2;
            int offset, globalrow, globalcol, cnt;

            NekDouble globalMatrixValue;
            NekDouble MatrixValue;
            NekDouble zero=0.0;
            DNekMatSharedPtr m_invS;

            int nGlobalBnd    = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd       = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();
            
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
            int cnt,i,j;
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            StdRegions::StdExpansionSharedPtr locExpansion;
            locExpansion = expList->GetExp(0);

            DNekScalBlkMatSharedPtr RtetBlk, RprismBlk;
            DNekScalBlkMatSharedPtr RTtetBlk, RTprismBlk;

            DNekScalMatSharedPtr Rtet, Rprism, Rmodified;
            DNekScalMatSharedPtr RTtet, RTprism, RTmodified;

            /*
             * Set up a Tetrahral & prismatic element which comprises
             * equilateral triangles as all faces for the tet and the end faces
             * for the prism. Using these elements a new expansion is created
             * (which is the same as the expansion specified in the input
             * file).*/

            SpatialDomains::TetGeomSharedPtr tetgeom=CreateRefTetGeom();
            SpatialDomains::PrismGeomSharedPtr prismgeom=CreateRefPrismGeom();

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
            RtetBlk = TetExp->GetLocStaticCondMatrix(TetR);
            RTtetBlk = TetExp->GetLocStaticCondMatrix(TetRT);
            Rtet=RtetBlk->GetBlock(0,0);
            RTtet=RTtetBlk->GetBlock(0,0);
            m_transformationMatrix[elmtType]=Rtet;
            m_transposedTransformationMatrix[elmtType]=RTtet;

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


            //For a tet element the bottom face is made up of the following:
            //vertices: 0, 1 and 2 edges: 0, 1 and 2 face: 0. We first need to
            //determine the mode locations of these vertices, edges and face so
            //we can extract the correct values from the tetrahedral R matrix.

            //These are the vertex mode locations of R which need to be replaced
            //in the prism element
            int TetVertex0;
            TetVertex0=TetExp->GetVertexMap(0);
            int TetVertex1;
            TetVertex1=TetExp->GetVertexMap(1);
            int TetVertex2;
            TetVertex2=TetExp->GetVertexMap(2);
            int TetVertex3;
            TetVertex3=TetExp->GetVertexMap(3);


            //These are the edge mode locations of R which need to be replaced
            //in the prism element - THESE ARE WRONG
            Array<OneD, unsigned int> TetEdge0;
            TetEdge0=TetExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> TetEdge1;
            TetEdge1=TetExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> TetEdge2;
            TetEdge2=TetExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> TetEdge3;
            TetEdge3=TetExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> TetEdge4;
            TetEdge4=TetExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> TetEdge5;
            TetEdge5=TetExp->GetEdgeInverseBoundaryMap(5);

            //These are the face mode locations of R which need to be replaced
            //in the prism element
            Array<OneD, unsigned int> TetFace;
            TetFace=TetExp->GetFaceInverseBoundaryMap(1);

            elmtType++;
            //Get a LocalRegions static condensed matrix
            RprismBlk = PrismExp->GetLocStaticCondMatrix(PrismR);
            RTprismBlk = PrismExp->GetLocStaticCondMatrix(PrismRT);
            Rprism=RprismBlk->GetBlock(0,0);
            RTprism=RTprismBlk->GetBlock(0,0);

            unsigned int  nRows=Rprism->GetRows();
            NekDouble zero=0.0;
            DNekMatSharedPtr Rmodprism = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);
            DNekMatSharedPtr RTmodprism = MemoryManager<DNekMat>::
                AllocateSharedPtr(nRows,nRows,zero,eFULL);
            NekDouble Rvalue, RTvalue;

            //Modified Prism - copy values from R
            for(i=0; i<nRows; ++i)
            {
                for(j=0; j<nRows; ++j)
                {
                    Rvalue=(*Rprism)(i,j);
                    RTvalue=(*RTprism)(i,j);
                    Rmodprism->SetValue(i,j,Rvalue);
                    RTmodprism->SetValue(i,j,RTvalue);
                }
            }

            //Prism vertex modes
            int PrismVertex0;
            PrismVertex0=PrismExp->GetVertexMap(0);
            int PrismVertex1;
            PrismVertex1=PrismExp->GetVertexMap(1);
            int PrismVertex2;
            PrismVertex2=PrismExp->GetVertexMap(2);
            int PrismVertex3;
            PrismVertex3=PrismExp->GetVertexMap(3);
            int PrismVertex4;
            PrismVertex4=PrismExp->GetVertexMap(4);
            int PrismVertex5;
            PrismVertex5=PrismExp->GetVertexMap(5);

            //Prism edge modes
            Array<OneD, unsigned int> PrismEdge0;
            PrismEdge0=PrismExp->GetEdgeInverseBoundaryMap(0);
            Array<OneD, unsigned int> PrismEdge1;
            PrismEdge1=PrismExp->GetEdgeInverseBoundaryMap(1);
            Array<OneD, unsigned int> PrismEdge2;
            PrismEdge2=PrismExp->GetEdgeInverseBoundaryMap(2);
            Array<OneD, unsigned int> PrismEdge3;
            PrismEdge3=PrismExp->GetEdgeInverseBoundaryMap(3);
            Array<OneD, unsigned int> PrismEdge4;
            PrismEdge4=PrismExp->GetEdgeInverseBoundaryMap(4);
            Array<OneD, unsigned int> PrismEdge5;
            PrismEdge5=PrismExp->GetEdgeInverseBoundaryMap(5);
            Array<OneD, unsigned int> PrismEdge6;
            PrismEdge6=PrismExp->GetEdgeInverseBoundaryMap(6);
            Array<OneD, unsigned int> PrismEdge7;
            PrismEdge7=PrismExp->GetEdgeInverseBoundaryMap(7);
            Array<OneD, unsigned int> PrismEdge8;
            PrismEdge8=PrismExp->GetEdgeInverseBoundaryMap(8);

            //Prism face 1 & 3 face modes
            Array<OneD, unsigned int> PrismFace1;
            PrismFace1=PrismExp->GetFaceInverseBoundaryMap(1);
            Array<OneD, unsigned int> PrismFace3;
            PrismFace3=PrismExp->GetFaceInverseBoundaryMap(3);

            //vertex 0 edge 0 3 & 4
            for(i=0; i< PrismEdge0.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex0,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge0[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex0,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge3[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex0,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex0,PrismEdge4[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge0[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge0[i],PrismVertex0,RTvalue);
                RTvalue=(*RTtet)(TetEdge2[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge3[i],PrismVertex0,RTvalue);
                RTvalue=(*RTtet)(TetEdge3[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge4[i],PrismVertex0,RTvalue);
            }

            //vertex 1 edge 0 1 & 5
            for(i=0; i< PrismEdge1.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex1,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge0[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex1,TetEdge1[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge1[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex1,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex1,PrismEdge5[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge0[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge0[i],PrismVertex1,RTvalue);
                RTvalue=(*RTtet)(TetEdge1[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge1[i],PrismVertex1,RTvalue);
                RTvalue=(*RTtet)(TetEdge4[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge5[i],PrismVertex1,RTvalue);
            }

            //vertex 2 edge 1 2 & 6
            for(i=0; i< PrismEdge2.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex2,TetEdge1[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge1[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex1,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge2[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex2,TetEdge5[i]);
                Rmodprism->SetValue(PrismVertex2,PrismEdge6[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge1[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge1[i],PrismVertex2,RTvalue);
                RTvalue=(*RTtet)(TetEdge0[i],TetVertex1);
                RTmodprism->SetValue(PrismEdge2[i],PrismVertex2,RTvalue);
                RTvalue=(*RTtet)(TetEdge5[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge6[i],PrismVertex2,RTvalue);
            }

            //vertex 3 edge 3 2 & 7
            for(i=0; i< PrismEdge3.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex2,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge3[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex0,TetEdge0[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge2[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex2,TetEdge5[i]);
                Rmodprism->SetValue(PrismVertex3,PrismEdge7[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge2[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge3[i],PrismVertex3,RTvalue);
                RTvalue=(*RTtet)(TetEdge0[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge2[i],PrismVertex3,RTvalue);
                RTvalue=(*RTtet)(TetEdge5[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge7[i],PrismVertex3,RTvalue);
            }

            //vertex 4 edge 4 5 & 8
            for(i=0; i< PrismEdge4.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex3,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge4[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex3,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge5[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex0,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex4,PrismEdge8[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge3[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge4[i],PrismVertex4,RTvalue);
                RTvalue=(*RTtet)(TetEdge4[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge5[i],PrismVertex4,RTvalue);
                RTvalue=(*RTtet)(TetEdge2[i],TetVertex0);
                RTmodprism->SetValue(PrismEdge8[i],PrismVertex4,RTvalue);
            }

            //vertex 5 edge 6 7 & 8
            for(i=0; i< PrismEdge5.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex3,TetEdge3[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge6[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex3,TetEdge4[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge7[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex2,TetEdge2[i]);
                Rmodprism->SetValue(PrismVertex5,PrismEdge8[i],Rvalue);

                //transposed values
                RTvalue=(*RTtet)(TetEdge3[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge6[i],PrismVertex5,RTvalue);
                RTvalue=(*RTtet)(TetEdge4[i],TetVertex3);
                RTmodprism->SetValue(PrismEdge7[i],PrismVertex5,RTvalue);
                RTvalue=(*RTtet)(TetEdge2[i],TetVertex2);
                RTmodprism->SetValue(PrismEdge8[i],PrismVertex5,RTvalue);
            }

            // face 1 vertices 0 1 4
            for(i=0; i< PrismFace1.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex0,TetFace[i]);
                Rmodprism->SetValue(PrismVertex0,PrismFace1[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex1,TetFace[i]);
                Rmodprism->SetValue(PrismVertex1,PrismFace1[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex3,TetFace[i]);
                Rmodprism->SetValue(PrismVertex4,PrismFace1[i],Rvalue);
                
                //transposed values
                RTvalue=(*RTtet)(TetFace[i],TetVertex0);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex0,RTvalue);
                RTvalue=(*RTtet)(TetFace[i],TetVertex1);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex1,RTvalue);
                RTvalue=(*RTtet)(TetFace[i],TetVertex3);
                RTmodprism->SetValue(PrismFace1[i],PrismVertex4,RTvalue);
            }

            // face 3 vertices 2, 3 & 5
            for(i=0; i< PrismFace3.num_elements(); ++i)
            {
                Rvalue=(*Rtet)(TetVertex1,TetFace[i]);
                Rmodprism->SetValue(PrismVertex2,PrismFace3[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex0,TetFace[i]);
                Rmodprism->SetValue(PrismVertex3,PrismFace3[i],Rvalue);
                Rvalue=(*Rtet)(TetVertex3,TetFace[i]);
                Rmodprism->SetValue(PrismVertex5,PrismFace3[i],Rvalue);
                
                //transposed values
                RTvalue=(*RTtet)(TetFace[i],TetVertex1);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex2,RTvalue);
                RTvalue=(*RTtet)(TetFace[i],TetVertex0);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex3,RTvalue);
                RTvalue=(*RTtet)(TetFace[i],TetVertex3);
                RTmodprism->SetValue(PrismFace3[i],PrismVertex5,RTvalue);
            }

            // Face 1 edge 0 4 5
            for(i=0; i< PrismFace1.num_elements(); ++i)
            {
                for(j=0; j<PrismEdge0.num_elements(); ++j)
                {
                    Rvalue=(*Rtet)(TetEdge0[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge0[j],PrismFace1[i],Rvalue);
                    Rvalue=(*Rtet)(TetEdge3[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge4[j],PrismFace1[i],Rvalue);
                    Rvalue=(*Rtet)(TetEdge4[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge5[j],PrismFace1[i],Rvalue);

                    //transposed values
                    RTvalue=(*RTtet)(TetEdge0[j],TetFace[i]);
                    RTmodprism->SetValue(PrismEdge0[j],PrismFace1[i],RTvalue);
                    RTvalue=(*RTtet)(TetEdge3[j],TetFace[i]);
                    RTmodprism->SetValue(PrismEdge4[j],PrismFace1[i],RTvalue);
                    RTvalue=(*RTtet)(TetEdge4[j],TetFace[i]);
                    RTmodprism->SetValue(PrismEdge5[j],PrismFace1[i],RTvalue);
                }
            }
                
            // Face 3 edge 2 6 7
            for(i=0; i< PrismFace3.num_elements(); ++i)
            {
                for(j=0; j<PrismEdge2.num_elements(); ++j)
                {
                    Rvalue=(*Rtet)(TetEdge0[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge2[j],PrismFace3[i],Rvalue);
                    Rvalue=(*Rtet)(TetEdge4[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge6[j],PrismFace3[i],Rvalue);
                    Rvalue=(*Rtet)(TetEdge3[j],TetFace[i]);
                    Rmodprism->SetValue(PrismEdge7[j],PrismFace3[i],Rvalue);

                    RTvalue=(*RTtet)(TetFace[i],TetEdge0[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge2[j],RTvalue);
                    RTvalue=(*RTtet)(TetFace[i],TetEdge4[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge6[j],RTvalue);
                    RTvalue=(*RTtet)(TetFace[i],TetEdge3[j]);
                    RTmodprism->SetValue(PrismFace3[i],PrismEdge7[j],RTvalue);
                }
            }

            DNekScalMatSharedPtr returnvalR = MemoryManager<DNekScalMat>
                                            ::AllocateSharedPtr(1.0,Rmodprism);

            DNekScalMatSharedPtr returnvalRT = MemoryManager<DNekScalMat>
              ::AllocateSharedPtr(1.0,RTmodprism);

            m_transformationMatrix[elmtType]=Rprism;
            m_transposedTransformationMatrix[elmtType]=RTprism;
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

            DNekMatSharedPtr    m_RS;
            DNekMatSharedPtr    m_RSRT;

            //Transformation matrices
            DNekMat R;
            DNekMat RT;
            DNekMat RS;
            DNekMat RSRT;
            
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

            //maps for different element types
            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transmatrixmap;
            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transposedtransmatrixmap;
            transmatrixmap[StdRegions::eTetrahedron]=m_transformationMatrix[0];
            transmatrixmap[StdRegions::ePrism]=m_transformationMatrix[1];
            transposedtransmatrixmap[StdRegions::eTetrahedron]=m_transposedTransformationMatrix[0];
            transposedtransmatrixmap[StdRegions::ePrism]=m_transposedTransformationMatrix[1];

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
            Array<OneD, int> m_edgemodeoffset(nNonDirEdgeIDs,0);
            Array<OneD, int> m_facemodeoffset(nNonDirFaceIDs,0);

            Array<OneD, int> m_edgeglobaloffset(nNonDirEdgeIDs,0);
            Array<OneD, int> m_faceglobaloffset(nNonDirFaceIDs,0);

            const Array<OneD, const ExpListSharedPtr>& bndCondExp = expList->GetBndCondExpansions();
            StdRegions::StdExpansion2DSharedPtr bndCondFaceExp;
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>& bndConditions = expList->GetBndConditions();

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

            const Array<OneD,const unsigned int>& nbdry_size
                    = m_locToGloMap->GetNumLocalBndCoeffsPerPatch();

            //Variants of R matrices required for low energy preconditioning
            m_RBlk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
            m_RTBlk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);
            m_S1Blk      = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nbdry_size, nbdry_size , blkmatStorage);

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

                //Get correct transformation matrix for element type
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
                                      + sign1*sign2*RSRT(vMap1,vMap2);

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

                                m_EdgeBlockArray[edgematrixoffset+v*nedgemodes+m]=globalEdgeValue;
                            }
                        }
                        edgematrixoffset+=maxEdgeDof*maxEdgeDof;
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
                                NekDouble globalFaceValue = sign1*sign2*RSRT(fMap1,fMap2);

                                //test here with local to global map
                                m_FaceBlockArray[facematrixoffset+v*nfacemodes+m]=globalFaceValue;
                            }
                        }
                        facematrixoffset+=nfacemodes*nfacemodes;
                        //facematrixoffset+=maxFaceDof*maxFaceDof;
                    }
                }

                //offset for the expansion
                cnt+=offset;

                //Here we build the block matrices for R and RT
                m_RBlk->SetBlock(n,n, transmatrixmap[eType]);
                m_RTBlk->SetBlock(n,n, transposedtransmatrixmap[eType]);
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
            if(nNonDirVerts!=0)
            {
                Gs::gs_data *tmp = Gs::Init(m_VertBlockToUniversalMap, m_comm);
                Gs::Gather(vertArray, Gs::gs_add, tmp);
            }

            //Exchange edge data over different processes
            Gs::gs_data *tmp1 = Gs::Init(m_EdgeBlockToUniversalMap, m_comm);
            Gs::Gather(m_GlobalEdgeBlock, Gs::gs_add, tmp1);

            //Exchange face data over different processes
            Gs::gs_data *tmp2 = Gs::Init(m_FaceBlockToUniversalMap, m_comm);
            Gs::Gather(m_GlobalFaceBlock, Gs::gs_add, tmp2);

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

                tmp_mat->Invert();
                BlkMat->SetBlock(i,i,tmp_mat);
            }
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
                {
                    if(solvertype == eIterativeStaticCond)
                    {
                        int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                        int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                        int nNonDir = nGlobal-nDir;
                        DNekBlkMat &M = (*BlkMat);
                         
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

        /**
         * \brief transform the solution vector vector to low energy
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

            NekVector<NekDouble> F_GlobBnd(nGlobBndDofs,pInput,eWrapper);
            NekVector<NekDouble> F_HomBnd(nGlobHomBndDofs,pOutput,
                                          eWrapper);
            
            DNekScalBlkMat &R = *m_RBlk;

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> F_LocBnd(nLocBndDofs,pLocal,eWrapper);
            Array<OneD, int> m_map = m_locToGloMap->GetLocalToGlobalBndMap();
            
            Vmath::Gathr(m_map.num_elements(), m_locToGloSignMult.get(), pInput.get(), m_map.get(), pLocal.get());
                        
            F_LocBnd=R*F_LocBnd;
            
            m_locToGloMap->AssembleBnd(F_LocBnd,F_HomBnd, nDirBndDofs);            
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
            Array<OneD, NekDouble>& pInput)
        {
            int nGlobBndDofs       = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBndDofs        = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs    = nGlobBndDofs - nDirBndDofs;
            int nLocBndDofs        = m_locToGloMap->GetNumLocalBndCoeffs();

            DNekScalBlkMat &RT = *m_RTBlk;
            NekVector<NekDouble> V_GlobHomBnd(nGlobHomBndDofs,pInput+nDirBndDofs,
              eWrapper);

            Array<OneD, NekDouble> pLocal(nLocBndDofs, 0.0);
            NekVector<NekDouble> V_LocBnd(nLocBndDofs,pLocal,eWrapper);

            Array<OneD, int> m_map = m_locToGloMap->GetLocalToGlobalBndMap();

            Array<OneD,NekDouble> tmp(nGlobBndDofs,0.0);

            //only want to map non-dirichlet dofs
            m_locToGloMap->GlobalToLocalBnd(V_GlobHomBnd,V_LocBnd, nDirBndDofs);

            V_LocBnd=RT*V_LocBnd;
            
            Vmath::Assmb(nLocBndDofs, m_locToGloSignMult.get(),pLocal.get(), m_map.get(), tmp.get());

            m_locToGloMap->UniversalAssembleBnd(tmp);

            Vmath::Vcopy(nGlobBndDofs-nDirBndDofs, tmp.get() + nDirBndDofs, 1, pInput.get() + nDirBndDofs, 1);
        }


        /**
         * \brief Set up the transformed block  matrix system
         *
         * Sets up a block elemental matrix in which each of the block matrix is
         * the low energy equivalent
         * i.e. \f$\mathbf{S}_{2}=\mathbf{R}\mathbf{S}_{1}\mathbf{R}^{T}\f$
         */     
        DNekScalBlkMatSharedPtr PreconditionerLowEnergy::
        v_TransformedSchurCompl(int offset)
	{
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
         
            StdRegions::StdExpansionSharedPtr locExpansion;                
            locExpansion = expList->GetExp(offset);
            int nbnd=locExpansion->NumBndryCoeffs();   
            int ncoeffs=locExpansion->GetNcoeffs();
            int nint=ncoeffs-nbnd;
            
            //This is the SC elemental matrix in the orginal basis (S1)
            DNekScalBlkMatSharedPtr loc_mat = (m_linsys.lock())->GetStaticCondBlock(expList->GetOffset_Elmt_Id(offset));
            DNekScalMatSharedPtr m_S1=loc_mat->GetBlock(0,0);

            //Transformation matrices 
            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transmatrixmap;
            map<StdRegions::ExpansionType,DNekScalMatSharedPtr> transposedtransmatrixmap;
            transmatrixmap[StdRegions::eTetrahedron]=m_transformationMatrix[0];
            transmatrixmap[StdRegions::ePrism]=m_transformationMatrix[1];
            transposedtransmatrixmap[StdRegions::eTetrahedron]=m_transposedTransformationMatrix[0];
            transposedtransmatrixmap[StdRegions::ePrism]=m_transposedTransformationMatrix[1];

            DNekScalMat &S1 = (*m_S1);
            
            NekDouble zero = 0.0;
            NekDouble one  = 1.0;
            MatrixStorage storage = eFULL;
            
            DNekMatSharedPtr m_S2 = MemoryManager<DNekMat>::AllocateSharedPtr(nbnd,nbnd,zero,storage);
            DNekMatSharedPtr m_RS1 = MemoryManager<DNekMat>::AllocateSharedPtr(nbnd,nbnd,zero,storage);
            
            StdRegions::ExpansionType eType=
                (expList->GetExp(offset))->DetExpansionType();
            
            //transformation matrices
            DNekScalMat &R = (*(transmatrixmap[eType]));
            DNekScalMat &RT = (*(transposedtransmatrixmap[eType]));
            
            //create low energy matrix
            DNekMat &RS1 = (*m_RS1);
            DNekMat &S2 = (*m_S2);
                
            //setup S2
            RS1=R*S1;
            S2=RS1*RT;

            DNekScalBlkMatSharedPtr returnval;
            DNekScalMatSharedPtr tmp_mat;
            unsigned int exp_size[] = {nbnd, nint};
            int nblks = 1;
            returnval = MemoryManager<DNekScalBlkMat>::
                AllocateSharedPtr(nblks, nblks, exp_size, exp_size);

            returnval->SetBlock(0,0,tmp_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,m_S2));

	    return returnval;
	}

        /**
         * Create the inverse multiplicity map.
         */
        void PreconditionerLowEnergy::CreateMultiplicityMap(void)
        {
            const Array<OneD, const int> &vMap
                                    = m_locToGloMap->GetLocalToGlobalBndMap();

            const Array< OneD, const NekDouble > &sign = m_locToGloMap->GetLocalToGlobalBndSign();

            unsigned int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            unsigned int nEntries   = m_locToGloMap->GetNumLocalBndCoeffs();
            unsigned int i;

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
                m_locToGloSignMult[i] = sign[i]*1.0/vCounts[vMap[i]];
            }

        }

    }
}






