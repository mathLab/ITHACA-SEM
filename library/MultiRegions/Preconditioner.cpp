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
#include <MultiRegions/Preconditioner.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <math.h>


namespace Nektar
{
    namespace MultiRegions
    {
        std::string Preconditioner::lookupIds[3] = {
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","Diagonal",MultiRegions::eDiagonal),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","InverseLinear",MultiRegions::eInverseLinear),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","LowEnergy",MultiRegions::eLowEnergy),
        };
        std::string Preconditioner::def = LibUtilities::SessionReader::RegisterDefaultSolverInfo("Preconditioner","Diagonal");

        /**
         * @class Preconditioner
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */

         Preconditioner::Preconditioner(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
                         const AssemblyMapSharedPtr &pLocToGloMap):
	   m_linsys(plinsys),
           m_locToGloMap(pLocToGloMap),
           m_preconType(pLocToGloMap->GetPreconType())
         {
             GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
             switch(m_preconType)
             {
             case MultiRegions::eDiagonal:
                 {
                     if (solvertype == eIterativeFull)
                     {
                         DiagonalPreconditionerSum();
                     }
                     else if(solvertype == eIterativeStaticCond)
                     {
                         StaticCondDiagonalPreconditionerSum();
                     }
                     else
                     {
                         ASSERTL0(0,"Unsupported solver type");
                     }
                 }
	         break;
             case MultiRegions::eInverseLinear:
                 {
                     if (solvertype == eIterativeFull)
                     {
                         LinearInversePreconditioner();
                     }
                     else if(solvertype == eIterativeStaticCond)
                     {
                         StaticCondLinearInversePreconditioner();
                     }
                     else
                     {
                         ASSERTL0(0,"Unsupported solver type");
                     }
                 }
                 break;
             case MultiRegions::eLowEnergy:
                 {
		   //ASSERTL0(0,"Not yet implemented");
                     StaticCondLinearInversePreconditioner();
                     LowEnergyPreconditioner();
                 }
                 break;
             default:
                 ASSERTL0(0,"Unknown preconditioner");
             break;
             }
        }

        /**
         * Populates preconditioner matrix with the identity i.e no
         * preconditioning.
         * @param   pLocToGloMap    Local to Global mapping.
         */
         void Preconditioner::NullPreconditioner(
                        const boost::weak_ptr<GlobalLinSys> &plinsys,
                        const AssemblyMapSharedPtr &pLocToGloMap)
        {
            int nGlobal = pLocToGloMap->GetNumGlobalCoeffs();
            int nDir    = pLocToGloMap->GetNumGlobalDirBndCoeffs();
            int nInt = nGlobal - nDir;
            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
            DNekMat &M = (*m_preconditioner);

            for (unsigned int i = 0; i < nInt; ++i)
            {
                M.SetValue(i,i,1.0);
            }
        }

        /**
         * Diagonal preconditioner computed by summing the relevant elements of
         * the local matrix system.
         * @param   pLocToGloMap    Local to global mapping.
         */
         void Preconditioner::DiagonalPreconditionerSum()
         {
             boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();

             const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
             StdRegions::StdExpansionSharedPtr locExpansion;

             int i,j,n,cnt,gid1,gid2;
             NekDouble sign1,sign2,value;
             int nGlobal = m_locToGloMap->GetNumGlobalCoeffs();
             int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
             int nInt    = nGlobal - nDir;

             NekDouble zero = 0.0;

             // fill global matrix
             DNekScalMatSharedPtr loc_mat;
             Array<OneD, NekDouble> vOutput(nGlobal,0.0);
             MatrixStorage storage = eDIAGONAL;
             m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(nInt, nInt, storage);
             DNekMat &M = (*m_preconditioner);

             int loc_lda;
             for(n = cnt = 0; n < expList->GetNumElmts(); ++n)
             {
                 loc_mat = (m_linsys.lock())->GetBlock(expList->GetOffset_Elmt_Id(n));
                 loc_lda = loc_mat->GetRows();

                 for(i = 0; i < loc_lda; ++i)
                 {
                     gid1 = m_locToGloMap->GetLocalToGlobalMap(cnt + i) - nDir;
                     sign1 =  m_locToGloMap->GetLocalToGlobalSign(cnt + i);
                     if(gid1 >= 0)
                     {
                         for(j = 0; j < loc_lda; ++j)
                         {
                             gid2 = m_locToGloMap->GetLocalToGlobalMap(cnt + j)
                                                                    - nDir;
                             sign2 = m_locToGloMap->GetLocalToGlobalSign(cnt + j);
                             if(gid2 == gid1)
                             {
                                 // When global matrix is symmetric,
                                 // only add the value for the upper
                                 // triangular part in order to avoid
                                 // entries to be entered twice
                                 value = vOutput[gid1 + nDir]
                                            + sign1*sign2*(*loc_mat)(i,j);
                                 vOutput[gid1 + nDir] = value;
                             }
                         }
                     }
                 }
                 cnt   += loc_lda;
             }

             // Assemble diagonal contributions across processes
             m_locToGloMap->UniversalAssemble(vOutput);

             // Populate preconditioner with reciprocal of diagonal elements
             for (unsigned int i = 0; i < nInt; ++i)
             {
                 M.SetValue(i,i,1.0/vOutput[i + nDir]);
             }
         }


        /**
         *  Performs global assembly of diagonal entries
         *  to global schur complement matrix.
         */
        Array<OneD, NekDouble> Preconditioner::AssembleStaticCondGlobalDiagonals()
        {
            int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int rows = nGlobalBnd - nDirBnd;

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;
            int sign1, sign2, gid1, gid2, i, j, n, cnt;
            Array<OneD, NekDouble> diagonals(rows,0.0);

            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();

            // Extract diagonal contributions of globally assembled
            // schur complement matrix
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
                //Get statically condensed local matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

                //Extract boundary block
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                int bnd_lda = bnd_mat->GetRows();

                for(i = 0; i < bnd_lda; ++i)
                {
                    gid1  = m_locToGloMap->GetLocalToGlobalBndMap (cnt + i) - nDirBnd;
                    sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < bnd_lda; ++j)
                        {
                            gid2  = m_locToGloMap->GetLocalToGlobalBndMap(cnt+j) - nDirBnd;
                            sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if(gid2 == gid1)
                            {
                                //vOutput[gid1 + nDirBnd] += sign1*sign2*(*bnd_mat)(i,j);
                                diagonals[gid1] += sign1*sign2*(*bnd_mat)(i,j);
                            }
                        }
                   }
                }
                cnt += bnd_lda;
            }
            return diagonals;
        }



        /**
         * Diagonal preconditioner defined as the inverse of the main
	 * diagonal of the Schur complement
	 *
         * @param   pLocToGloMap    Local to global mapping.
         */
        void Preconditioner::StaticCondDiagonalPreconditionerSum()
        {
            int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int rows = nGlobalBnd - nDirBnd;

            MatrixStorage storage = eDIAGONAL;
            m_preconditioner = MemoryManager<DNekMat>::AllocateSharedPtr(rows, rows, storage);
            DNekMat &M = (*m_preconditioner);

            Array<OneD, NekDouble> vOutput(nGlobalBnd,0.0);

            // Extract diagonal contributions
            Array<OneD, NekDouble> diagonals = AssembleStaticCondGlobalDiagonals();
            for (unsigned int i = 0; i < rows; ++i)
            {
                vOutput[nDirBnd + i] = diagonals[i];
            }

            // Assemble diagonal contributions across processes
            m_locToGloMap->UniversalAssembleBnd(vOutput);

            // Populate preconditioner matrix
            for (unsigned int i = 0; i < rows; ++i)
            {
                M.SetValue(i,i,1.0/vOutput[nDirBnd + i]);
            }
        }



        /**
         *
         */         
        void Preconditioner::LinearInversePreconditioner()
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
            S.Invert();

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
        void Preconditioner::StaticCondLinearInversePreconditioner()
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
            
            //Get Global statically condensed matrix
            //DNekMatSharedPtr m_gmat=(m_linsys.lock())->GetGmat();
            int gRow = nGlobalBnd - nDirBnd; //m_gmat->GetRows();
            
            
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

            for(int i = 0; i < S.GetRows(); ++i)
            {
                for(int j = 0; j < S.GetColumns(); ++j)
                {
                    MatrixValue=S.GetValue(i,j);
		    cout<<MatrixValue<<" ";
                }
		cout<<endl;
            }

            //Invert linear finite element space
            S.Invert();

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
//                M.SetValue(i,i,1.0/(*m_gmat)(i,i));
                  M.SetValue(i,i,1.0/diagonals[i]);
            }
        }

        void Preconditioner::LowEnergyPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
            StdRegions::StdExpansionSharedPtr locExpansion;

            //Local regions matrix and geometrical info
            CreateReferenceGeometryAndMatrix();

            //Determine the low energy modes
            SetLowEnergyModes_Rv();

            SetLowEnergyModes_Ref();

            //Construct the low energy basis and preconditioner
            BuildPreconditioner();
        }


        /**
	 * \brief Create reference element and statically condensed matrix
	 *
	 **/
        void Preconditioner::CreateReferenceGeometryAndMatrix()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();

            int nVerts, nEdges, nFaces, bndry_rows;
            int vMap, eid, fid, vid, cnt, n, i, j;
            int nEdgeCoeffs, nFaceCoeffs;

            DNekScalBlkMatSharedPtr loc_mat;

            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();

            //offset of preconditioning element
            int nel = expList->GetOffset_Elmt_Id(0);

            //only need a single local matrix for this method; this will
            //be changed in future versions.
            vExp = expList->GetExp(nel);

            const Array<OneD, const LibUtilities::BasisSharedPtr>& m_base=vExp->GetBase();

            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            // retrieve variable coefficients
            if(m_linSysKey.GetNVarCoeffs() > 0)
            {
                StdRegions::VarCoeffMap::const_iterator x;
                cnt = expList->GetPhys_Offset(n);
                for (x = m_linSysKey.GetVarCoeffs().begin(); x != m_linSysKey.GetVarCoeffs().end(); ++x)
                {
                    vVarCoeffMap[x->first] = x->second + cnt;
                }
            }

            //Local matrix key - the matrix key "ePreconditioner"
            //is a helmholz matrix constructed from a equalateral
            //Tetrahedron or Hexahedron
            LocalRegions::MatrixKey matkey(StdRegions::ePreconditioner,
                                           vExp->DetExpansionType(),
                                           *vExp,
                                           m_linSysKey.GetConstFactors(),
                                           vVarCoeffMap);

            //Get a LocalRegions static condensed matrix
            loc_mat = vExp->GetLocStaticCondMatrix(matkey);

            //local schur complement (boundary-boundary block)
            bndry_mat = loc_mat->GetBlock(0,0);
	    
            //number of rows=columns of the schur complement
            bndry_rows=bndry_mat->GetRows();

            int nCoeffs=vExp->GetNcoeffs();
            int nint=nCoeffs-bndry_rows;

            Array<OneD,unsigned int> bmap(bndry_rows);
            vExp->GetBoundaryMap(bmap);
	    
            Array<OneD,unsigned int> imap(nint);
            vExp->GetInteriorMap(imap);

            //map from full system to statically condensed system
            //i.e reverse GetBoundaryMap

	    map<int,int> invmap;
            for(j = 0; j < bmap.num_elements(); ++j)
            {
                invmap[bmap[j]] = j;
            }

            //Get geometric information about this element
            nVerts=vExp->GetGeom()->GetNumVerts();
            nEdges=vExp->GetGeom()->GetNumEdges();
            nFaces=vExp->GetGeom()->GetNumFaces();

            //Set up map between element vertex, edge or face on the reference
            //element and modes in the matrix
            vertModeLocation = Array<OneD, int > (nVerts);
            edgeModeLocation = Array<OneD, Array<OneD, unsigned int> > (nEdges);
            faceModeLocation = Array<OneD, Array<OneD, unsigned int> > (nFaces);
	    
            //Get expansion type
            StdRegions::ExpansionType eType=vExp->DetExpansionType();
	    
            //loop over vertices and determine the location of vertex coefficients in the storage array
            for (vid=0; vid<nVerts; ++vid)
            {
                //location in matrix
                vMap = vExp->GetVertexMap(vid);
                vertModeLocation[vid]=vMap;
            }

            //loop over edges and determine location of edge coefficients in the storage array
            for (eid=0; eid<nEdges; ++eid)
            {
                //Number of interior edge coefficients
                nEdgeCoeffs=vExp->GetEdgeNcoeffs(eid)-2;

                StdRegions::Orientation eOrient=vExp->GetGeom()->GetEorient(eid);
                Array< OneD, unsigned int > maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
                Array< OneD, int > signarray = Array<OneD, int>(nEdgeCoeffs,1);

                //maparray is the location of the edge within the matrix
                vExp->GetEdgeInteriorMap(eid,eOrient,maparray,signarray);

                for (n=0; n<maparray.num_elements(); ++n)
                {
                    maparray[n]=invmap[maparray[n]];
                }
                edgeModeLocation[eid]=maparray;
            }

            int nTotFaceCoeffs=vExp->GetTotalFaceIntNcoeffs();
            Array< OneD, unsigned int > teststorage = Array<OneD, unsigned int>(nTotFaceCoeffs);
            
            //loop over edges and determine location of face coefficients in the storage array
            for (cnt=fid=0; fid<nFaces; ++fid)
            {
                //Number of interior edge coefficients
                nFaceCoeffs=vExp->GetFaceIntNcoeffs(fid);
 
                StdRegions::Orientation fOrient=vExp->GetFaceOrient(fid);
                Array< OneD, unsigned int > maparray = Array<OneD, unsigned int>(nFaceCoeffs);
                Array< OneD, int > signarray = Array<OneD, int>(nFaceCoeffs,1);

                //maparray is the location of the edge within the matrix
                vExp->GetFaceInteriorMap(fid,fOrient,maparray,signarray);
		
                for (n=0; n<maparray.num_elements(); ++n)
                {
                    maparray[n]=invmap[maparray[n]];
                }
                faceModeLocation[fid]=maparray;
            }
        }


        /**
	 * \brief Build vertex transformation matrix \f$\mathbf{R_{v}}\f$
	 *
	 * The matrix component of \f$\mathbf{R}\f$ is given by
	 *\f[
	 *  \mathbf{R^{T}_{v}}=-\mathbf{S}^{-1}_{ef,ef}\mathbf{S}^{T}_{v,ef}\f]
	 *
	 * For every vertex mode we extract the submatrices from statically condensed 
	 * matrix \f$\mathbf{S}\f$ corresponding to the coupling between the attached 
	 * edges and faces of a vertex (\f$\mathbf{S_{ef,ef}}\f$). This matrix is then
	 * inverted and multiplied by the submatrix representing the coupling between
	 * a vertex and the attached edges and faces (\f$\mathbf{S_{v,ef}}\f$). 
	 */
        void Preconditioner::SetLowEnergyModes_Rv()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            int nVerts, nEdges, nFaces, bndry_rows, nmodes;
            int vMap, eid, eid2, fid, fid2, vid, cnt, cnt2, n, m;
            int FaceTotNCoeffs, EdgeTotNCoeffs;
            NekDouble MatrixValue, VertexEdgeFaceValue;
            NekDouble zero = 0.0;

            //The number of connected edges/faces is 3 (for all elements)
            int nConnectedEdges=3;
            int nConnectedFaces=3;

            //location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedEdges);
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedFaces);

            int nCoeffs=vExp->NumBndryCoeffs();

            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            m_transformationmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(
                                     nCoeffs, nCoeffs, zero, storage);
            DNekMat &R = (*m_transformationmatrix);

            //Build the vertex-edge/face transform matrix: This matrix is constructed
            //from the submatrices corresponding to the couping between each vertex
            //and the attached edges/faces
            for(vid=0; vid<vExp->GetGeom()->GetNumVerts(); ++vid)
            {
                //row and column size of the vertex-edge/face matrix
                int efRow = vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,0)) +
                            vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,1)) +
                            vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,2)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,0)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,1)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,2)) - 6;

                //Allocation of vector to store vertex-edge/face coupling
                Array<OneD, NekDouble> vef_array(efRow, 0.0);

                int nedgemodesconnected=nConnectedEdges * (vExp->GetEdgeNcoeffs(
                                                           vExp->GetGeom()->GetVertexEdgeMap(vid,0))-2);
                Array<OneD, unsigned int> edgemodearray(nedgemodesconnected);

                int nfacemodesconnected=nConnectedFaces * (vExp->GetFaceIntNcoeffs(
                                                           vExp->GetGeom()->GetVertexFaceMap(vid,0)));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);

                //Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> Vvef(efRow,vef_array,eWrapper);

                //create array of edge modes
                for(eid=0; eid < nConnectedEdges; ++eid)
                {
                    MatEdgeLocation[eid]=edgeModeLocation[vExp->GetGeom()->GetVertexEdgeMap(vid,eid)];
                    nmodes=MatEdgeLocation[eid].num_elements();
                    Vmath::Vcopy(nmodes, &MatEdgeLocation[eid][0], 1, &edgemodearray[eid*nmodes], 1);
                }

                //create array of face modes
                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    MatFaceLocation[fid]=faceModeLocation[vExp->GetGeom()->GetVertexFaceMap(vid,fid)];
                    nmodes=MatFaceLocation[fid].num_elements();
                    Vmath::Vcopy(nmodes, &MatFaceLocation[fid][0], 1, &facemodearray[fid*nmodes], 1);
                }

                //vertex-edge coupling
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=(*bndry_mat)(edgemodearray[n],vertModeLocation[vid]);

                    //Set the value in the vertex edge/face matrix
                    Vvef[n]=VertexEdgeFaceValue;
                }

                //vertex-face coupling
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=(*bndry_mat)(facemodearray[n],vertModeLocation[vid]);

                    //Set the value in the vertex edge/face matrix
                    Vvef[n+nedgemodesconnected]=VertexEdgeFaceValue;
                }


                /*Build the edge-face transform matrix: This matrix is constructed
                  from the submatrices corresponding to the couping between the edges
                  and faces on the attached faces/edges of a vertex*/

                //Allocation of matrix to store edge/face-edge/face coupling
                m_vertexedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                           efRow, efRow,zero, storage);
                DNekMat &Mvef = (*m_vertexedgefacecoupling);
                NekDouble EdgeEdgeValue, FaceFaceValue;

                //edge-edge coupling (S_{ee})
                for (m=0; m<nedgemodesconnected; ++m)
                {
                    for (n=0; n<nedgemodesconnected; ++n)
                    {
                        //Matrix value for each coefficient location
                        EdgeEdgeValue=(*bndry_mat)(edgemodearray[m],edgemodearray[n]);

                        //Set the value in the vertex edge/face matrix
                        Mvef.SetValue(m,n,EdgeEdgeValue);
                    }
                }

                //face-face coupling (S_{ff})
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bndry_mat)(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Mvef.SetValue(nedgemodesconnected+n,nedgemodesconnected+m,FaceFaceValue);
                    }
                }

                //edge-face coupling (S_{ef} and trans(S_{ef}))
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bndry_mat)(edgemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix (and transpose)
                        Mvef.SetValue(n,nedgemodesconnected+m,FaceFaceValue);

                        //and transpose
                        Mvef.SetValue(nedgemodesconnected+m,n,FaceFaceValue);
                    }
                }                

                // Invert edge-face coupling matrix
                Mvef.Invert();

                // Allocate array storage
                Array<OneD, NekDouble> rv_array(efRow, 0.0);

                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> Rv_component(efRow,rv_array,eWrapper);

                //trans(R_{v})=-inv(S_{ef,ef})*trans(S_{v,ef})
                Rv_component=-(Mvef*Vvef);

                // Populate R with R_{ve} components
                for(n=0; n<edgemodearray.num_elements(); ++n)
                {
                    R.SetValue(edgemodearray[n],vertModeLocation[vid],Rv_component[n]);
                }

                // Populate R with R_{vf} components
                for(n=0; n<facemodearray.num_elements(); ++n)
                {
                    R.SetValue(facemodearray[n],vertModeLocation[vid],Rv_component[n+nedgemodesconnected]);
                }
            }
        }


        /**
	 * \brief Build edge-face transformation matrix (\f$\mathbf{R_{ef}}\f$)
	 *
	 * The matrix component of \f$\mathbf{R}\f$ is given by
	 *\f[
	 *  \mathbf{R^{T}_{ef}}=-\mathbf{S}^{-1}_{ff}\mathbf{S}^{T}_{ef}\f]
	 *
	 * For each edge extract the submatrices from statically condensed 
	 * matrix \f$\mathbf{S}\f$ corresponding to inner products of modes on the two
	 * attached faces within themselves as well as the coupling matrix 
	 * between the two faces (\f$\mathbf{S}_{ff}\f$). This matrix of face coupling is then inverted 
	 * and multiplied by the submatrices of corresponding to the 
	 * coupling between the edge and attached faces (\f$\mathbf{S}_{ef}\f$).
	 *
	 */
        void Preconditioner::SetLowEnergyModes_Ref()
        {
            int eid, fid, fid2, cnt, i, j, n, m, nmodes, nedgemodes;
            int efRow, efCol, FaceTotNCoeffs, EdgeTotNCoeffs;
            NekDouble zero = 0.0;
	    
            NekDouble EdgeFaceValue, FaceFaceValue, Rvalue;
	    
            //number of attached faces is always 2
            int nConnectedFaces=2;

            //location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> > (vExp->GetGeom()->GetNumEdges());
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedFaces);

            FaceTotNCoeffs=vExp->GetTotalFaceIntNcoeffs();
            EdgeTotNCoeffs=vExp->GetTotalEdgeIntNcoeffs();

            // Define storage for vertex transpose matrix
            MatrixStorage storage = eFULL;
            DNekMat &R = (*m_transformationmatrix);

            //Build the edge/face transform matrix: This matrix is constructed
            //from the submatrices corresponding to the couping between a specific
            //edge and the two attached faces.
            for (cnt=eid=0; eid<vExp->GetGeom()->GetNumEdges(); ++eid)
            {
                //row and column size of the vertex-edge/face matrix
                efRow=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,0))+
                      vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,1));
                efCol=vExp->GetEdgeNcoeffs(eid)-2;

                // Edge-face coupling matrix
                m_efedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                       efRow, efCol, zero, storage);
                DNekMat &Mef = (*m_efedgefacecoupling);

                // Face-face coupling matrix
                m_effacefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                       efRow, efRow, zero, storage);
                DNekMat &Meff = (*m_effacefacecoupling);

                // Edge-face transformation matrix
                m_edgefacetransformmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(
                                            efRow, efCol, zero, storage);
                DNekMat &Meft = (*m_edgefacetransformmatrix);

                int nfacemodesconnected=nConnectedFaces * (vExp->GetFaceIntNcoeffs(
                                                           vExp->GetGeom()->GetEdgeFaceMap(eid,0)));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);

                //create array of edge modes
                nedgemodes=edgeModeLocation[eid].num_elements();
                Array<OneD, unsigned int> edgemodearray(nedgemodes);
                Vmath::Vcopy(nedgemodes, &edgeModeLocation[eid][0], 1, &edgemodearray[0], 1);

                //create array of face modes
                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    MatFaceLocation[fid]=faceModeLocation[vExp->GetGeom()->GetEdgeFaceMap(eid,fid)];
                    nmodes=MatFaceLocation[fid].num_elements();
                    Vmath::Vcopy(nmodes, &MatFaceLocation[fid][0], 1, &facemodearray[fid*nmodes], 1);
                }

                //edge-face coupling
                for (n=0; n<nedgemodes; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        EdgeFaceValue=(*bndry_mat)(facemodearray[m],edgemodearray[n]);

                        //Set the value in the edge/face matrix
                        Mef.SetValue(m,n,EdgeFaceValue);
                    }
                }

                //face-face coupling
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bndry_mat)(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Meff.SetValue(n,m,FaceFaceValue);
                    }
                }

                // Invert edge-face coupling matrix
                Meff.Invert();

                // trans(R_{ef})=-inv(S_{ef})*trans(S_{ef})
                Meft=-Meff*Mef;

                //Populate transformation matrix with Meft
                for(n=0; n<Meft.GetRows(); ++n)
                {
                    for(m=0; m<Meft.GetColumns(); ++m)
                    {
                        Rvalue=Meft(n,m);
                        R.SetValue(facemodearray[n],edgemodearray[m],Rvalue);
                    }
                }
            }

            for (i = 0; i < R.GetRows(); ++i)
            {
                R.SetValue(i,i,1.0);
            }
        }

        /**
	 * \brief Construct the low energy preconditioner
	 *
	 *\f[\mathbf{R}^{T}\left[\begin{array}{ccc} Diag[(\mathbf{S_{2}})_{vv}] & & \\
	 *  & (\mathbf{S}_{2})_{eb} & \\
	 *  &  & (\mathbf{S}_{2})_{fb} \end{array}\right] \mathbf{R} \f]
	 *
	 * where \f$\mathbf{R}\f$ is the transformation matrix and \f$\mathbf{S}_{2}\f$
	 * the Schur complement of the modified basis, given by
	 *
	 * \f[\mathbf{S}_{2}=\mathbf{R}\mathbf{S}_{1}\mathbf{R}^{T}\f]
	 *
	 */
        void Preconditioner::BuildPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
             const StdRegions::StdExpansionVector &locExpVector = *(expList->GetExp());
             StdRegions::StdExpansionSharedPtr locExpansion;

            int nRow, i, j, nmodes, ntotaledgemodes, ntotalfacemodes;
            int nVerts, nEdges,nFaces, eid, fid, n, cnt, nedgemodes, nfacemodes;
	    int nEdgeCoeffs, nFaceCoeffs;
            NekDouble zero = 0.0;
            NekDouble MatrixValue;

            int vMap1, vMap2, sign1, sign2, gid1, gid2, m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;
            NekDouble globalMatrixValue, globalRValue;

            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            DNekMat &R = (*m_transformationmatrix);
            DNekScalMat &S=(*bndry_mat);

            DNekMatSharedPtr    m_SR;
            DNekMatSharedPtr    m_RSRT;

            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            nRow=R.GetRows();

            m_SR = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &SR = (*m_SR);
            m_RSRT = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &RSRT = (*m_RSRT);
            m_SP = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &SP = (*m_SP);

            //Calculate S*trans(R)  (note R is already transposed)
            SR=S*R;

            //Transpose R i.e R->trans(R)
            R.Transpose();
	    
            //Calculate R*S*trans(R)
            RSRT=R*SR;

            /*for (i=0; i<R.GetRows(); ++i)
            {
                for (int j=0; j<R.GetColumns(); ++j)
                {
                    cout<<RSRT(i,j)<<" ";
                }
                cout<<endl;
            }

	    cout<<endl;
	    cout<<endl;

            for (i=0; i<S.GetRows(); ++i)
            {
                for (int j=0; j<S.GetColumns(); ++j)
                {
                    cout<<S(i,j)<<" ";
                }
                cout<<endl;
		}*/

            nVerts=vExp->GetGeom()->GetNumVerts();
            nEdges=vExp->GetGeom()->GetNumEdges();
            nFaces=vExp->GetGeom()->GetNumFaces();

            int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int nNonDir = nGlobal-nDirBnd;
            int nNonDirVerts  = m_locToGloMap->GetNumNonDirVertexModes();
            int nNonDirEdges  = m_locToGloMap->GetNumNonDirEdgeModes();
            int nNonDirFaces  = m_locToGloMap->GetNumNonDirFaceModes();

            // set up block matrix system
            int nblks = 3;
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


            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

	    // Set vertex block
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
		offset=vExp->NumBndryCoeffs();

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
				  + sign1*sign2*RSRT(vMap1,vMap2);

                                //build matrix containing the linear finite element space
                                VertBlk->SetValue(globalrow,globalcol,globalMatrixValue);
                            }
                        }
                    }
                }
                cnt   += offset;
            }

 	    //set edge block
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
		offset=vExp->NumBndryCoeffs();	

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

                                    //Global matrix value
                                    globalMatrixValue = EdgeBlk->GetValue(globalrow,globalcol)
				      + sign1*sign2*RSRT(eMap1,eMap2);

                                    //build matrix containing the linear finite element space
                                    EdgeBlk->SetValue(globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
		    }
                }
                cnt   += offset;
            }

	    //set face block
            for(cnt=n=0; n < expList->GetNumElmts(); ++n)
            {
		offset=vExp->NumBndryCoeffs();	

                //loop over edges of the element and return the edge map
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

                                    //Global matrix value
                                    globalMatrixValue = FaceBlk->GetValue(globalrow,globalcol)
				      + sign1*sign2*RSRT(fMap1,fMap2);

                                    //build matrix containing the linear finite element space
                                    FaceBlk->SetValue(globalrow,globalcol,globalMatrixValue);
                                }
                            }
                        }
		    }
                }
                cnt   += offset;
            }

	    VertBlk->Invert();
	    EdgeBlk->Invert();
	    FaceBlk->Invert();

            DNekScalMatSharedPtr     Blktmp;
            NekDouble                one = 1.0;

            GloBlkMat->SetBlock(0,0,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,VertBlk));
            GloBlkMat->SetBlock(1,1,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,EdgeBlk));
            GloBlkMat->SetBlock(2,2,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,FaceBlk));
        }


        void Preconditioner::BuildPreconditioner_Reordered()
        {
            int nRow, i, j, nmodes, ntotaledgemodes, ntotalfacemodes;
            int nVerts, nEdges,nFaces, eid, fid;
            NekDouble zero = 0.0;
            NekDouble MatrixValue;

            MatrixStorage storage = eFULL;
            DNekMat &R = (*m_transformationmatrix);
            DNekMat &S=(*m_om);

            DNekMatSharedPtr    m_SR;
            DNekMatSharedPtr    m_RSRT;

            nRow=R.GetRows();

            m_SR = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &SR = (*m_SR);
            m_RSRT = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &RSRT = (*m_RSRT);
            m_SP = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &SP = (*m_SP);

            DNekMatSharedPtr  m_SPV;
            DNekMatSharedPtr  m_SPE;
            DNekMatSharedPtr  m_SPF;

            //Calculate S*trans(R)  (note R is already transposed)
            SR=S*R;

            //Transpose R i.e R->trans(R)
            R.Transpose();
	    
            //Calculate R*S*trans(R)
            RSRT=R*SR;

        }

        /**
	 * \brief Reorder the schur complement matrix into vertex, edge and face storage
	 *
	 * \f[\mathbf{S} = \left[ \begin{array}{ccc}
              \mathbf{S_{vv}} & \mathbf{S_{ve}} & \mathbf{S_{vf}} \\
              \mathbf{S^{T}_{ve}} & \mathbf{S^{T}_{ee}} & \mathbf{S_{ef}} \\
              \mathbf{S^{T}_{vf}} & \mathbf{S^{T}_{ef}} & \mathbf{S_{ff}} \end{array} \right]\f]
	 */
        void Preconditioner::VertexEdgeFaceMatrix()
	{
            int i,j, eid, fid;
            int nCoeffs=vExp->NumBndryCoeffs();
            NekDouble MatrixValue;
 
            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            m_om = MemoryManager<DNekMat>::AllocateSharedPtr(nCoeffs, nCoeffs, 0, storage);
            DNekMat &OM = (*m_om);

            int nVerts=vExp->GetGeom()->GetNumVerts();
            int nEdges=vExp->GetGeom()->GetNumEdges();
            int nFaces=vExp->GetGeom()->GetNumFaces();
            int nedgemodes, nfacemodes;

            Array<OneD, unsigned int> edgemodearray(nEdges*edgeModeLocation[0].num_elements());
            Array<OneD, unsigned int> facemodearray(nFaces*faceModeLocation[0].num_elements());

            //create array of edge modes
            for(eid=0; eid < nEdges; ++eid)
            {
                nedgemodes=edgeModeLocation[eid].num_elements();
                Vmath::Vcopy(nedgemodes, &edgeModeLocation[eid][0], 1, &edgemodearray[eid*nedgemodes], 1);
            }

            //create array of face modes
            for(fid=0; fid < nFaces; ++fid)
            {
                nfacemodes=faceModeLocation[fid].num_elements();
                Vmath::Vcopy(nfacemodes, &faceModeLocation[fid][0], 1, &facemodearray[fid*nfacemodes], 1);
            }
 
            int nedgemodestotal=nedgemodes*nEdges;
            int nfacemodestotal=nfacemodes*nFaces;

            //vertex-vertex/edge/face
            for (i=0; i<nVerts; ++i)
            {
                for(j=0; j<nVerts; ++j)
                {
                    MatrixValue=(*bndry_mat)(vertModeLocation[i],vertModeLocation[j]);

                    OM.SetValue(i,j,MatrixValue);
                }
 
                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(vertModeLocation[i],edgemodearray[j]);
                    OM.SetValue(i,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(vertModeLocation[i],facemodearray[j]);
                    OM.SetValue(i,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }

            //edge-vertex/edge/face
            for (i=0; i<nedgemodestotal; ++i)
            {
                for(j=0; j<nVerts; ++j)
                {
                    MatrixValue=(*bndry_mat)(edgemodearray[i],vertModeLocation[j]);
                    OM.SetValue(i+nVerts,j,MatrixValue);
                }

                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(edgemodearray[i],edgemodearray[j]);
                    OM.SetValue(i+nVerts,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(edgemodearray[i],facemodearray[j]);
                    OM.SetValue(i+nVerts,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }

            //face-vertex/edge/face
            for (i=0; i<nfacemodestotal; ++i)
            {
                for(j=0; j<nVerts; ++j)
                {
                    MatrixValue=(*bndry_mat)(facemodearray[i],vertModeLocation[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j,MatrixValue);
                }
 
                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(facemodearray[i],edgemodearray[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=(*bndry_mat)(facemodearray[i],facemodearray[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }
        }


        /**
	 *
	 */
        void Preconditioner::SetLowEnergyModes_Rordered()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            int nVerts, nEdges, nFaces, bndry_rows, nmodes;
            int vMap, eid, eid2, fid, fid2, vid, cnt, cnt2, n, m;
            int FaceTotNCoeffs, EdgeTotNCoeffs;
            int edge,face, nedgemodes, nfacemodes;
            NekDouble MatrixValue, VertexEdgeFaceValue;
            NekDouble zero = 0.0;

            //The number of connected edges/faces is 3 (for all elements)
            int nConnectedEdges=3;
            int nConnectedFaces=3;

            //location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedEdges);
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedFaces);

            nVerts=vExp->GetGeom()->GetNumVerts();
            int nCoeffs=vExp->NumBndryCoeffs();

            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            DNekMat &OM = (*m_om);

            m_transformationmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(
                                     nCoeffs, nCoeffs, zero, storage);
            DNekMat &R = (*m_transformationmatrix);

            //Build the vertex-edge/face transform matrix: This matrix is constructed
            //from the submatrices corresponding to the couping between each vertex
            //and the attached edges/faces
            for(vid=0; vid<nVerts; ++vid)
            {
                //row and column size of the vertex-edge/face matrix
                int efRow = vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,0)) +
                            vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,1)) +
                            vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,2)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,0)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,1)) +
                            vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,2)) - 6;

                //Allocation of vector to store vertex-edge/face coupling
                Array<OneD, NekDouble> vef_array(efRow, 0.0);

                int nedgemodesconnected=nConnectedEdges*(vExp->GetEdgeNcoeffs(
                                                         vExp->GetGeom()->GetVertexEdgeMap(vid,0))-2);
                Array<OneD, unsigned int> edgemodearray(nedgemodesconnected);
 
                int nfacemodesconnected=nConnectedFaces*(vExp->GetFaceIntNcoeffs(
                                                         vExp->GetGeom()->GetVertexFaceMap(vid,0)));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);

                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> Vvef(efRow,vef_array,eWrapper);

                //create array of edge modes
                for(eid=0; eid < nConnectedEdges; ++eid)
                {
                    nedgemodes=vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,eid))-2;
                    edge=vExp->GetGeom()->GetVertexEdgeMap(vid,eid);
                    for(n=0; n<nedgemodes; ++n)
                    {
                        edgemodearray[nedgemodes*eid+n]=nVerts+edge*nedgemodes+n;
                    }
                }

                nEdges=vExp->GetGeom()->GetNumEdges();
                nFaces=vExp->GetGeom()->GetNumFaces();

                //create array of face modes
                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    nfacemodes=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,fid));
                    face=vExp->GetGeom()->GetVertexFaceMap(vid,fid);
                    for(n=0; n<nfacemodes; ++n)
                    {
                        facemodearray[nfacemodes*fid+n]=nVerts+nEdges*nedgemodes+face*nfacemodes+n;
                    }
                }


                //Get rows where edge modes, of a specific vertex-edge coupling, are located
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=OM(edgemodearray[n],vid);
 
                    //Set the value in the vertex edge/face matrix
                    Vvef[n]=VertexEdgeFaceValue;
                }

                //Get rows where face modes, of a specific vertex vertex-face coupling, are located
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=OM(facemodearray[n],vid);
 
                    //Set the value in the vertex edge/face matrix
                    Vvef[n+nedgemodesconnected]=VertexEdgeFaceValue;
                }


                /*Build the edge-face transform matrix: This matrix is constructed
                  from the submatrices corresponding to the couping between the edges
                  and faces on the attached faces/edges of a vertex*/

                //Allocation of matrix to store edge/face-edge/face coupling
                m_vertexedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(efRow, efRow, storage);
                DNekMat &Mvef = (*m_vertexedgefacecoupling);
                NekDouble EdgeEdgeValue, FaceFaceValue;

                //Get rows where edge modes, of a specific edge-edge coupling, are located
                for (m=0; m<nedgemodesconnected; ++m)
                {
                    for (n=0; n<nedgemodesconnected; ++n)
                    {
                        //Matrix value for each coefficient location
                        EdgeEdgeValue=OM(edgemodearray[m],edgemodearray[n]);
 
                        //Set the value in the vertex edge/face matrix
                        Mvef.SetValue(m,n,EdgeEdgeValue);
                    }
                }


                //face-face coupling (S_{ff})
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=OM(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Mvef.SetValue(nedgemodesconnected+n,nedgemodesconnected+m,FaceFaceValue);
                    }
                }



                //edge-face coupling (S_{ef} and trans(S_{ef}))
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=OM(edgemodearray[n],facemodearray[m]);
 
                        //Set the value in the vertex edge/face matrix
                        Mvef.SetValue(n,nedgemodesconnected+m,FaceFaceValue);
 
                        //Set transposed values
                        Mvef.SetValue(nedgemodesconnected+m,n,FaceFaceValue);
                    }
                }

                // Invert edge-face coupling matrix
                Mvef.Invert();

                // Allocate array storage
                Array<OneD, NekDouble> rv_array(efRow, 0.0);

                // Create NekVector wrappers for linear algebra operations
                NekVector<NekDouble> Rv_component(efRow,rv_array,eWrapper);

                // Multiply the inverse edge-face (Mef) coupling matrix
                // with the vertex-face/edge vector (Vvef)
                Rv_component=-(Mvef*Vvef);

                for(eid=0; eid < nConnectedEdges; ++eid)
                {
                    nedgemodes=vExp->GetEdgeNcoeffs(vExp->GetGeom()->GetVertexEdgeMap(vid,eid))-2;
                    edge=vExp->GetGeom()->GetVertexEdgeMap(vid,eid);
                    for(n=0; n<nedgemodes; ++n)
                    {
                        R.SetValue(nVerts+nedgemodes*edge+n,vid,Rv_component[nedgemodes*eid+n]);
                    }
                }

                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    nfacemodes=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetVertexFaceMap(vid,fid));
                    face=vExp->GetGeom()->GetVertexFaceMap(vid,fid);
                    for(n=0; n<nfacemodes; ++n)
                    {
                        R.SetValue(nVerts+nedgemodes*nEdges+face*nfacemodes+n,
                                   vid,Rv_component[nedgemodes*nConnectedEdges+nfacemodes*fid+n]);
                    }
                }
            }
        }

        /**
	 *
	 */
        void Preconditioner::SetLowEnergyModes_Ref_Reordered()
        {
            int eid, fid, fid2, cnt, i, n, m, nmodes, nedgemodes;
            int efRow, efCol, FaceTotNCoeffs, EdgeTotNCoeffs;
	    
            NekDouble EdgeFaceValue, FaceFaceValue, Rvalue;
	    
            //number of attached faces is always 2
            int nConnectedFaces=2;
                
            int nEdges=vExp->GetGeom()->GetNumEdges();
            int nVerts=vExp->GetGeom()->GetNumVerts();

            //location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> > (vExp->GetGeom()->GetNumEdges());
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedFaces);

            FaceTotNCoeffs=vExp->GetTotalFaceIntNcoeffs();
            EdgeTotNCoeffs=vExp->GetTotalEdgeIntNcoeffs();

            // Define storage for vertex transpose matrix
            MatrixStorage storage = eFULL;
            DNekMat &R = (*m_transformationmatrix);
            // Define storage for vertex transpose matrix and zero all entries

            DNekMat &OM = (*m_om);

            //Get rows where face modes, of a specific edge-face coupling, are located
            for (cnt=eid=0; eid<vExp->GetGeom()->GetNumEdges(); ++eid)
            {
                //row and column size of the vertex-edge/face matrix
                efRow=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,0))+
                      vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,1));
                efCol=vExp->GetEdgeNcoeffs(eid)-2;

                //Edge-face coupling matrix
                m_efedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(efRow, efCol, storage);
                DNekMat &Mef = (*m_efedgefacecoupling);

                //Face-face coupling matrix
                m_effacefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(efRow, efRow, storage);
                DNekMat &Meff = (*m_effacefacecoupling);

                //Edge-face transformation matrix
                m_edgefacetransformmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(efRow, efCol, storage);
                DNekMat &Meft = (*m_edgefacetransformmatrix);

                int nfacemodesconnected=nConnectedFaces * (vExp->GetFaceIntNcoeffs(
                                                           vExp->GetGeom()->GetEdgeFaceMap(eid,0)));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);

                //create array of edge modes
                nedgemodes=edgeModeLocation[eid].num_elements();
                Array<OneD, unsigned int> edgemodearray(nedgemodes);

                //create array of edge modes
                for(n=0; n<nedgemodes; ++n)
                {
                    edgemodearray[n]=nVerts+eid*nedgemodes+n;
                }

                int face;
                int nfacemodes;

                //create array of edge modes
                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    nfacemodes=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,fid));
                    face=vExp->GetGeom()->GetEdgeFaceMap(eid,fid);
                    for(n=0; n<nfacemodes; ++n)
                    { 
                        facemodearray[nfacemodes*fid+n]=nVerts+nEdges*nedgemodes+face*nfacemodes+n;
                    }
                }


                //Get rows where face modes, of a specific edge-face coupling, are located
                for (n=0; n<nedgemodes; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        EdgeFaceValue=OM(facemodearray[m],edgemodearray[n]);
 
                        //Set the value in the edge/face matrix
                        Mef.SetValue(m,n,EdgeFaceValue);
                    }
                }

                //Get rows where face modes, of a specific face-face coupling, are located
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=OM(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Meff.SetValue(n,m,FaceFaceValue);
                    }
                }

                //Invert edge-face coupling matrix
                Meff.Invert();

                //
                Meft=Meff*Mef;


                for(fid=0; fid < nConnectedFaces; ++fid)
                {
                    nfacemodes=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,fid));
                    face=vExp->GetGeom()->GetEdgeFaceMap(eid,fid);
                    for(n=0; n<nfacemodes; ++n)
                    {
                        for(m=0; m<nedgemodes; ++m)
                        {
                            Rvalue=Meft(fid*nfacemodes+n,m);
                            R.SetValue(nVerts+nedgemodes*nEdges+face*nfacemodes+n,
                                       nVerts+nedgemodes*eid+m,-Rvalue);
                        }
                    }
                }
            }

            for (i = 0; i < R.GetRows(); ++i)
            {
                R.SetValue(i,i,1.0);
            }

        }


        /**
         *
         */
        void Preconditioner::DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(m_preconType)
            {
            case MultiRegions::eDiagonal:
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
		    int i;
                    DNekScalMatSharedPtr     Blktmp;
                    NekDouble                one = 1.0;
		    boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
                    int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                    int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
		    int nLocal  = m_locToGloMap->GetNumLocalBndCoeffs();
                    int nNonDir = nGlobal-nDir;
                    DNekScalBlkMat &M = (*GloBlkMat);
		    DNekMat &R = (*m_transformationmatrix);

		    int nblks = expList->GetExpSize();
		    Array<OneD,unsigned int> exp_size(nblks);
		    
		    for(i=0; i<expList->GetExpSize(); ++i)
		    {
		        exp_size[i]=R.GetRows();
		    }

                    MatrixStorage blkmatStorage = eDIAGONAL;
                    DNekScalBlkMatSharedPtr RBlkMat = MemoryManager<DNekScalBlkMat>::AllocateSharedPtr(exp_size,exp_size,blkmatStorage);

		    for(i=0; i<expList->GetExpSize(); ++i)
		    {
		      RBlkMat->SetBlock(i,i,Blktmp = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,m_transformationmatrix));
		    }

                    DNekScalBlkMat &RL = (*RBlkMat);

		    NekVector<NekDouble> r1(nLocal,0.0);
		    NekVector<NekDouble> z1(nLocal,0.0);
                    NekVector<NekDouble> Globalin(nNonDir,pInput,eWrapper);

		    m_locToGloMap->GlobalToLocalBnd(Globalin,r1,nDir);

		    r1 = RL * r1;

		    NekVector<NekDouble> z1glo(nNonDir,0.0);
		    m_locToGloMap->AssembleBnd(r1,z1glo,nDir);

                    NekVector<NekDouble> z2glo(nNonDir,0.0);

		    z2glo = M * z1glo;

		    NekVector<NekDouble> z2(nLocal,0.0);
		    m_locToGloMap->GlobalToLocalBnd(z2glo,z2,nDir);
		    RL.Transpose();

		    z2 = RL * z2;

                    NekVector<NekDouble> Globalout(nNonDir,pOutput,eWrapper);

		    m_locToGloMap->AssembleBnd(z2,Globalout,nDir);

		}
		break;
            default:
            ASSERTL0(0,"Unknown preconditioner");
            break;
	    }
	}
    }
}








