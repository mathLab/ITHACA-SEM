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

                    //CreateReferenceGeometryAndMatrix(StdRegions::eTetrahedron);

                    //SetUpLowEnergyBasis();
                    LowEnergyPreconditioner();
		}
		break;
            case MultiRegions::eBlock:
                {
                    CreateReferenceGeometryAndMatrix(StdRegions::eTetrahedron);
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
	 *
	 *
         */         
        void PreconditionerLowEnergy::InverseLinearSpacePreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;

            int localVertId, vMap1, vMap2, nVerts, localCoeff, n, v, m;
            int sign1, sign2, gid1, gid2, i, j;
            int loc_rows, globalrow, globalcol, cnt, globalLocation, nIntEdgeFace;
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
            m_S = MemoryManager<DNekMat>::AllocateSharedPtr(nNonDirVerts, nNonDirVerts, zero,  storage);
            DNekMat &S = (*m_S);
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, zero, storage);
            DNekMat &M = (*m_preconditioner);
 
            DNekScalMatSharedPtr loc_mat;

            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
                //element matrix
                loc_mat = (m_linsys.lock())->GetBlock(expList->GetOffset_Elmt_Id(n));
                loc_rows = loc_mat->GetRows();

                //element expansion
                locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(
                                      locExpVector[expList->GetOffset_Elmt_Id(n)]);

                //Get number of vertices
                nVerts=locExpansion->GetGeom()->GetNumVerts();

                //loop over vertices of the element and return the vertex map for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    //Get vertex map
                    vMap1 = locExpansion->GetVertexMap(v);

                    globalrow = m_locToGloMap->GetLocalToGlobalMap(cnt+vMap1)-nDir;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2 = locExpansion->GetVertexMap(m);

                            //global matrix location (with offset due to dirichlet values)
                            globalcol = m_locToGloMap->GetLocalToGlobalMap(cnt+vMap2)-nDir;

                            if(globalcol>=0)
                            {

                                //modal connectivity between elements
                                sign1 = m_locToGloMap->GetLocalToGlobalSign(cnt + vMap1);
                                sign2 = m_locToGloMap->GetLocalToGlobalSign(cnt + vMap2);

                                //Global matrix value
                                globalMatrixValue = S.GetValue(globalrow,globalcol)
                                                  + sign1*sign2*(*loc_mat)(vMap1,vMap2);
                        
                                //build matrix containing the linear finite element space
                                S.SetValue(globalrow,globalcol,globalMatrixValue);
                            }
                        }
                    }
                }

                //move counter down length of loc_rows
                cnt   += loc_rows;
            }

            int loc_lda;
            for(n = cnt = 0; n < expList->GetNumElmts(); ++n)
            {
                loc_mat = (m_linsys.lock())->GetBlock(expList->GetOffset_Elmt_Id(n));
                loc_lda = loc_mat->GetRows();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = m_locToGloMap->GetLocalToGlobalMap(cnt + i) - nDir-nNonDirVerts;
                    sign1 =  m_locToGloMap->GetLocalToGlobalSign(cnt + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = m_locToGloMap->GetLocalToGlobalMap(cnt + j)
                                 - nDir-nNonDirVerts;
                            sign2 = m_locToGloMap->GetLocalToGlobalSign(cnt + j);
                            if(gid2 == gid1)
                            {
                                // When global matrix is symmetric,
                                // only add the value for the upper
                                // triangular part in order to avoid
                                // entries to be entered twice
                                value = vOutput[gid1 + nDir + nNonDirVerts]
                                      + sign1*sign2*(*loc_mat)(i,j);
                                vOutput[gid1 + nDir + nNonDirVerts] = value;
                            }
                        }
                    }
                }
                cnt   += loc_lda;
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
        void PreconditionerLowEnergy::StaticCondInverseLinearSpacePreconditioner()
	{
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;

            int localVertId, vMap1, vMap2, nVerts, nEdges, nFaces, localCoeff, v, m, n, rows;
            int sign1, sign2;
            int offset, globalrow, globalcol, cnt, cnt1, globalLocation, nDirVertOffset;
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
            m_S = MemoryManager<DNekMat>::AllocateSharedPtr(nNonDirVerts, nNonDirVerts, zero,  storage);
            DNekMat &S = (*m_S);

            //element expansion
            locExpansion = boost::dynamic_pointer_cast<StdRegions::StdExpansion>(
                                  locExpVector[expList->GetOffset_Elmt_Id(0)]);

            //Get total number of vertices
            nVerts=locExpansion->GetGeom()->GetNumVerts();

            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(
                               nGlobalBnd-nDirBnd, nGlobalBnd-nDirBnd, zero,  storage);
            DNekMat &M = (*m_preconditioner);

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);
		
                //Extract boundary block
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                offset = bnd_mat->GetRows();
		
                //loop over vertices of the element and return the vertex map for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    //Get vertex map
                    vMap1 = locExpansion->GetVertexMap(v);
                    globalrow = m_locToGloMap->GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2 = locExpansion->GetVertexMap(m);

                            //global matrix location (without offset due to dirichlet values)
                            globalcol = m_locToGloMap->GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol >= 0)
                            {
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + vMap2);

                                //Global matrix value
                                globalMatrixValue = S.GetValue(globalrow,globalcol)
                                                  + sign1*sign2*(*bnd_mat)(vMap1,vMap2);

                                //build matrix containing the linear finite element space
                                S.SetValue(globalrow,globalcol,globalMatrixValue);
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

            Array<OneD, NekDouble> diagonals = AssembleStaticCondGlobalDiagonals();

            // Populate preconditioner matrix
            for (unsigned int i = nNonDirVerts; i < M.GetRows(); ++i)
            {
                  M.SetValue(i,i,1.0/diagonals[i]);
            }
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
            int offset, globalrow, globalcol, bnd_rows;
            NekDouble globalMatrixValue, globalRValue;

            vExp = expList->GetExp(0);
            int nCoeffs=vExp->NumBndryCoeffs();

            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekScalBlkMatSharedPtr r_mat;
            DNekScalBlkMatSharedPtr rt_mat;

            DNekMatSharedPtr    m_RS;
            DNekMatSharedPtr    m_RSRT;

            DNekMat R;
            DNekMat RT;
            
            DNekScalMatSharedPtr m_transformationmatrix;
            DNekScalMatSharedPtr m_transposedtransformationmatrix;


            m_RS = MemoryManager<DNekMat>::AllocateSharedPtr
                (nCoeffs, nCoeffs, zero, storage);
            DNekMat &RS = (*m_RS);
            m_RSRT = MemoryManager<DNekMat>::AllocateSharedPtr
                (nCoeffs, nCoeffs, zero, storage);
            DNekMat &RSRT = (*m_RSRT);
            
            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            nVerts=vExp->GetGeom()->GetNumVerts();
            nEdges=vExp->GetGeom()->GetNumEdges();
            nFaces=vExp->GetGeom()->GetNumFaces();

            Array<OneD, int>                           vertModeLocation(nVerts);
            Array<OneD, Array<OneD, unsigned int> >    edgeModeLocation(nEdges);
            Array<OneD, Array<OneD, unsigned int> >    faceModeLocation(nFaces);

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

            
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {

                nel = expList->GetOffset_Elmt_Id(n);

                locExpansion = expList->GetExp(nel);

                StdRegions::ExpansionType eType=
                    locExpansion->DetExpansionType();



                if(eType==StdRegions::eTetrahedron && !m_transformationmatrix)
                {
                    // retrieve variable coefficient
                    if(m_linSysKey.GetNVarCoeffs() > 0)
                    {
                        StdRegions::VarCoeffMap::const_iterator x;
                        cnt = expList->GetPhys_Offset(nel);
                        for (x = m_linSysKey.GetVarCoeffs().begin(); 
                             x != m_linSysKey.GetVarCoeffs().end(); ++x)
                        {
                            vVarCoeffMap[x->first] = x->second + cnt;
                        }
                    }

                    LocalRegions::MatrixKey r_matkey
                        (StdRegions::ePreconR,
                         eType,
                         *locExpansion,
                         m_linSysKey.GetConstFactors(),
                         vVarCoeffMap);
                    
                    LocalRegions::MatrixKey rt_matkey
                        (StdRegions::ePreconRT,
                         eType,
                         *locExpansion,
                         m_linSysKey.GetConstFactors(),
                         vVarCoeffMap);
                    
                    
                    //Get a LocalRegions static condensed matrix
                    r_mat = locExpansion->GetLocStaticCondMatrix(r_matkey);
                    rt_mat = locExpansion->GetLocStaticCondMatrix(rt_matkey);
                    
                    m_transformationmatrix=r_mat->GetBlock(0,0);
                    m_transposedtransformationmatrix=rt_mat->GetBlock(0,0);
                    
                    R=(*m_transformationmatrix);
                    RT=(*m_transposedtransformationmatrix);
                    
                    locExpansion->GetModeMappings
                        (vertModeLocation,edgeModeLocation,faceModeLocation);

                    //number of rows=columns of the schur complement
                    bnd_rows=m_transformationmatrix->GetRows();
                }
                
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
v                        
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
                for(i=0; i<EdgeBlk->GetRows(); ++i)
                {
                    for (j=0; j<EdgeBlk->GetRows(); ++j)
                    {
                        cout<<(*EdgeBlk)(i,j)<<" ";
                    }
                    cout<<endl;
                }
                EdgeBlk->Invert();
            }

            if (nNonDirFaces != 0)
            {
                for(i=0; i<FaceBlk->GetRows(); ++i)
                {
                    for (j=0; j<FaceBlk->GetRows(); ++j)
                    {
                        cout<<(*FaceBlk)(i,j)<<" ";
                    }
                    cout<<endl;
                }

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
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;

            int nRow, i, j, nmodes, ntotaledgemodes, ntotalfacemodes;
            int nVerts, nEdges,nFaces, eid, fid, eid2, fid2, n, cnt, nedgemodes, nfacemodes;
            int nEdgeCoeffs, nFaceCoeffs;
            NekDouble zero = 0.0;
            NekDouble MatrixValue;

            int vMap1, vMap2, sign1, sign2, gid1, gid2, m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;
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
            GloBlkMat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size,blkmatStorage);

            //Vertex, edge and face preconditioner matrices
            DNekMatSharedPtr VertBlk = MemoryManager<DNekMat>::AllocateSharedPtr(nNonDirVerts,nNonDirVerts,zero,vertstorage);
            DNekMatSharedPtr EdgeBlk = MemoryManager<DNekMat>::AllocateSharedPtr(nNonDirEdges,nNonDirEdges,zero,storage);
            DNekMatSharedPtr FaceBlk = MemoryManager<DNekMat>::AllocateSharedPtr(nNonDirFaces,nNonDirFaces,zero,storage);

            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
                //Get statically condensed matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);
		
                //Extract boundary block (elemental S1)
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                offset = bnd_mat->GetRows();

                DNekScalMat &S=(*bnd_mat);

                //loop over vertices of the element and return the vertex map for each vertex
                for (v=0; v<nVerts; ++v)
                {
                    vMap1=vertModeLocation[v];
                    
                    //Get vertex map
                    globalrow = m_locToGloMap->GetLocalToGlobalBndMap(cnt+vMap1)-nDirBnd;

                    if(globalrow >= 0)
                    {
                        for (m=0; m<nVerts; ++m)
                        {
                            vMap2=vertModeLocation[m];

                            //global matrix location (without offset due to dirichlet values)
                            globalcol = m_locToGloMap->GetLocalToGlobalBndMap(cnt+vMap2)-nDirBnd;

                            //offset for dirichlet conditions
                            if (globalcol == globalrow)
                            {
                                //modal connectivity between elements
                                sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + vMap1);
                                sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + vMap2);

                                //Global matrix value
                                globalMatrixValue = VertBlk->GetValue(globalrow,globalcol)
                                                  + sign1*sign2*S(vMap1,vMap2);

                                //build matrix containing the linear finite element space
                                VertBlk->SetValue(globalrow,globalcol,globalMatrixValue);
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

                        globalrow = m_locToGloMap->GetLocalToGlobalBndMap(cnt+eMap1)-nDirBnd-nNonDirVerts;

                        if(globalrow >= 0)
                        {
                            for (m=0; m<nedgemodes; ++m)
                            {
                                eMap2=edgeModeLocation[eid][m];

                                //global matrix location (without offset due to dirichlet values)
                                globalcol = m_locToGloMap->GetLocalToGlobalBndMap(cnt+eMap2)-nDirBnd-nNonDirVerts;

                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + eMap1);
                                    sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + eMap2);

                                    globalMatrixValue = EdgeBlk->GetValue(globalrow,globalcol)
                                                      + sign1*sign2*S(eMap1,eMap2);
		
                                    //build matrix containing the linear finite element space
                                    EdgeBlk->SetValue(globalrow,globalcol,globalMatrixValue);
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

                        globalrow = m_locToGloMap->GetLocalToGlobalBndMap(cnt+fMap1)-nDirBnd-nNonDirVerts-nNonDirEdges;

                        if(globalrow >= 0)
                        {
                            for (m=0; m<nfacemodes; ++m)
                            {
                                fMap2=faceModeLocation[fid][m];

                                //global matrix location (without offset due to dirichlet values)
                                globalcol = m_locToGloMap->GetLocalToGlobalBndMap(cnt+fMap2)-nDirBnd-nNonDirVerts-nNonDirEdges;

                                //offset for dirichlet conditions
                                if (globalcol >= 0)
                                {
                                    //modal connectivity between elements
                                    sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + fMap1);
                                    sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + fMap2);

                                    globalMatrixValue = FaceBlk->GetValue(globalrow,globalcol)
                                                      + sign1*sign2*S(fMap1,fMap2);

                                    //build matrix containing the linear finite element space
                                    FaceBlk->SetValue(globalrow,globalcol,globalMatrixValue);
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

            GloBlkMat->SetBlock(0,0,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,VertBlk));
            GloBlkMat->SetBlock(1,1,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,EdgeBlk));
            GloBlkMat->SetBlock(2,2,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,FaceBlk));
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


    }
}






