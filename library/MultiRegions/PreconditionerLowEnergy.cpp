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

                    SetUpLowEnergyBasis();
                    LowEnergyPreconditioner();
		}
		break;
            case MultiRegions::eBlock:
                {
                    CreateReferenceGeometryAndMatrix();
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

            int vMap1, vMap2, nVerts, v, m, n;
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

        void PreconditionerLowEnergy::SetUpLowEnergyBasis()
        {
            //Local regions matrix and geometrical info
            CreateReferenceGeometryAndMatrix();

            //Determine the low energy modes
            SetLowEnergyModes_Rv();

            SetLowEnergyModes_Ref();

            SetUpInverseTransformationMatrix();

        }

        /**
	 * \brief Create reference element and statically condensed matrix
	 *
	 **/
        void PreconditionerLowEnergy::CreateReferenceGeometryAndMatrix()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();

            int nVerts, nEdges, nFaces, bndry_rows;
            int vMap, eid, fid, vid, cnt, n,  j;
            int nEdgeCoeffs, nFaceCoeffs;

            DNekScalBlkMatSharedPtr loc_mat;

            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();

            //offset of preconditioning element
             int nel = expList->GetOffset_Elmt_Id(0);

            //only need a single local matrix for this method.
            vExp = expList->GetExp(nel);

            // need to be initialised with zero size for non variable
            // coefficient case
            StdRegions::VarCoeffMap vVarCoeffMap;

            // retrieve variable coefficient1
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

            if(m_preconType == MultiRegions::eLowEnergy || m_preconType == MultiRegions::eLocalLowEnergy)
            {
                LocalRegions::MatrixKey matkey(StdRegions::ePreconditioner,
                                               vExp->DetShapeType(),
                                               *vExp,
                                               m_linSysKey.GetConstFactors(),
                                               vVarCoeffMap);

                //Get a LocalRegions static condensed matrix
                loc_mat = vExp->GetLocStaticCondMatrix(matkey);
            }
            else if(m_preconType == MultiRegions::eBlock)
            {
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(nel);
            }
            else
            {
                ASSERTL0(0,"Unknown preconditiner for this method");
            }

            //local schur complement (boundary-boundary block)
            bnd_mat = loc_mat->GetBlock(0,0);

            //number of rows=columns of the schur complement
            bndry_rows=bnd_mat->GetRows();

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

            //loop over faces and determine location of face coefficients in the storage array
            for (cnt=fid=0; fid<nFaces; ++fid)
            {
                //Number of interior edge coefficients
                nFaceCoeffs=vExp->GetFaceIntNcoeffs(fid);
 
                StdRegions::Orientation fOrient=vExp->GetFaceOrient(fid);
                Array< OneD, unsigned int > maparray = Array<OneD, unsigned int>(nFaceCoeffs);
                Array< OneD, int > signarray = Array<OneD, int>(nFaceCoeffs,1);

                //maparray is the location of the face within the matrix
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
        void PreconditionerLowEnergy::SetLowEnergyModes_Rv()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            int nmodes;
            int eid, fid, vid, n, m;
            NekDouble VertexEdgeFaceValue;
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
            m_transposedtransformationmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(
                                     nCoeffs, nCoeffs, zero, storage);
            DNekMat &R = (*m_transformationmatrix);
            DNekMat &RT = (*m_transposedtransformationmatrix);

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

                int nedgemodesconnected=nConnectedEdges * (vExp->GetEdgeNcoeffs(
                                                           vExp->GetGeom()->GetVertexEdgeMap(vid,0))-2);
                Array<OneD, unsigned int> edgemodearray(nedgemodesconnected);

                int nfacemodesconnected=nConnectedFaces * (vExp->GetFaceIntNcoeffs(
                                                           vExp->GetGeom()->GetVertexFaceMap(vid,0)));
                Array<OneD, unsigned int> facemodearray(nfacemodesconnected);


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
                
                m_vertexedgefacetransformmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(
											    1, efRow, zero, storage);
                DNekMat &Sveft = (*m_vertexedgefacetransformmatrix);

                m_vertexedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
										     1, efRow, zero, storage);
                DNekMat &Svef = (*m_vertexedgefacecoupling);

                //vertex-edge coupling
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=(*bnd_mat)(vertModeLocation[vid], edgemodearray[n]);

                    //Set the value in the vertex edge/face matrix
                    Svef.SetValue(0,n,VertexEdgeFaceValue);
                }

                //vertex-face coupling
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    //Matrix value for each coefficient location
                    VertexEdgeFaceValue=(*bnd_mat)(vertModeLocation[vid],facemodearray[n]);

                    //Set the value in the vertex edge/face matrix
                    //Svef.SetValue(vid,n+nedgemodesconnected,VertexEdgeFaceValue);
                    Svef.SetValue(0,n+nedgemodesconnected,VertexEdgeFaceValue);
                }


                /*Build the edge-face transform matrix: This matrix is constructed
                  from the submatrices corresponding to the couping between the edges
                  and faces on the attached faces/edges of a vertex*/

                //Allocation of matrix to store edge/face-edge/face coupling
                m_edgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                           efRow, efRow,zero, storage);
                DNekMat &Sefef = (*m_edgefacecoupling);


                NekDouble EdgeEdgeValue, FaceFaceValue;

                //edge-edge coupling (S_{ee})
                for (m=0; m<nedgemodesconnected; ++m)
                {
                    for (n=0; n<nedgemodesconnected; ++n)
                    {
                        //Matrix value for each coefficient location
                        EdgeEdgeValue=(*bnd_mat)(edgemodearray[n],edgemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Sefef.SetValue(n,m,EdgeEdgeValue);
                    }
                }

                //face-face coupling (S_{ff})
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bnd_mat)(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Sefef.SetValue(nedgemodesconnected+n,nedgemodesconnected+m,FaceFaceValue);
                    }
                }

                //edge-face coupling (S_{ef} and trans(S_{ef}))
                for (n=0; n<nedgemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bnd_mat)(edgemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix (and transpose)
                        Sefef.SetValue(n,nedgemodesconnected+m,FaceFaceValue);

                        //and transpose
                        Sefef.SetValue(nedgemodesconnected+m,n,FaceFaceValue);
                    }
                }                

                // Invert edge-face coupling matrix
                Sefef.Invert();

                //R_{v}=-S_{v,ef}inv(S_{ef,ef})
                Sveft=-Svef*Sefef;

                // Populate R with R_{ve} components
                for(n=0; n<edgemodearray.num_elements(); ++n)
                {
                    RT.SetValue(edgemodearray[n], vertModeLocation[vid], Sveft(0,n));
                    R.SetValue(vertModeLocation[vid], edgemodearray[n], Sveft(0,n));
                }

                // Populate R with R_{vf} components
                for(n=0; n<facemodearray.num_elements(); ++n)
                {
                    RT.SetValue(facemodearray[n], vertModeLocation[vid], Sveft(0,n+nedgemodesconnected));
                    R.SetValue(vertModeLocation[vid], facemodearray[n], Sveft(0,n+nedgemodesconnected));
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
        void PreconditionerLowEnergy::SetLowEnergyModes_Ref()
        {
            int eid, fid,  cnt, i, n, m, nmodes, nedgemodes;
            int efRow, efCol, FaceTotNCoeffs, EdgeTotNCoeffs;
            NekDouble zero = 0.0;
	    
            NekDouble EdgeFaceValue, FaceFaceValue;
	    
            //number of attached faces is always 2
            int nConnectedFaces=2;
            int nEdges=vExp->GetGeom()->GetNumEdges();

            //location in the matrix
            MatEdgeLocation = Array<OneD, Array<OneD, unsigned int> > (vExp->GetGeom()->GetNumEdges());
            MatFaceLocation = Array<OneD, Array<OneD, unsigned int> > (nConnectedFaces);

            FaceTotNCoeffs=vExp->GetTotalFaceIntNcoeffs();
            EdgeTotNCoeffs=vExp->GetTotalEdgeIntNcoeffs();

            // Define storage for vertex transpose matrix
            MatrixStorage storage = eFULL;
            DNekMat &R = (*m_transformationmatrix);
            DNekMat &RT = (*m_transposedtransformationmatrix);

            //Build the edge/face transform matrix: This matrix is constructed
            //from the submatrices corresponding to the couping between a specific
            //edge and the two attached faces.
            for (cnt=eid=0; eid<nEdges; ++eid)
            {
                //row and column size of the vertex-edge/face matrix
                efCol=vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,0))+
                      vExp->GetFaceIntNcoeffs(vExp->GetGeom()->GetEdgeFaceMap(eid,1));
                efRow=vExp->GetEdgeNcoeffs(eid)-2;

                // Edge-face coupling matrix
                m_efedgefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                       efRow, efCol, zero, storage);
                DNekMat &Mef = (*m_efedgefacecoupling);

                // Face-face coupling matrix
                m_effacefacecoupling = MemoryManager<DNekMat>::AllocateSharedPtr(
                                       efCol, efCol, zero, storage);
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
                        EdgeFaceValue=(*bnd_mat)(edgemodearray[n],facemodearray[m]);

                        //Set the value in the edge/face matrix
                        Mef.SetValue(n,m,EdgeFaceValue);
                    }
                }

                //face-face coupling
                for (n=0; n<nfacemodesconnected; ++n)
                {
                    for (m=0; m<nfacemodesconnected; ++m)
                    {
                        //Matrix value for each coefficient location
                        FaceFaceValue=(*bnd_mat)(facemodearray[n],facemodearray[m]);

                        //Set the value in the vertex edge/face matrix
                        Meff.SetValue(n,m,FaceFaceValue);
                    }
                }

                // Invert edge-face coupling matrix
                Meff.Invert();

                // trans(R_{ef})=-S_{ef}*(inv(S_{ff})
                Meft=-Mef*Meff;

                //Populate transformation matrix with Meft
                for(n=0; n<Meft.GetRows(); ++n)
                {
                    for(m=0; m<Meft.GetColumns(); ++m)
                    {
                        R.SetValue(edgemodearray[n], facemodearray[m], Meft(n,m));
                        RT.SetValue(facemodearray[m], edgemodearray[n], Meft(n,m));
                    }
                }
            }

            for (i = 0; i < R.GetRows(); ++i)
            {
                R.SetValue(i,i,1.0);
                RT.SetValue(i,i,1.0);
            }
        }

        /**
	 * \brief Build inverse and inverse transposed transformation matrix: \f$\mathbf{R^{-1}}\f$ and \f$\mathbf{R^{-T}}\f$
	 *
	 * \f\mathbf{R^{-T}}=[\left[\begin{array}{ccc} \mathbf{I} & -\mathbf{R}_{ef} & -\mathbf{R}_{ve}+\mathbf{R}_{ve}\mathbf{R}_{vf} \\
	 *  0 & \mathbf{I} & \mathbf{R}_{ef} \\
	 *  0 & 0 & \mathbf{I}} \end{array}\right]\f]
	 *
	 */
        void PreconditionerLowEnergy::SetUpInverseTransformationMatrix()
	{
            int i,j,n, eid, fid;
            int nCoeffs=vExp->NumBndryCoeffs();
            NekDouble MatrixValue;
            NekDouble zero=0.0;
            DNekMat &R = (*m_transformationmatrix);
            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            m_inversetransformationmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(nCoeffs, nCoeffs, 
                                                                                        zero, storage);
            DNekMat &InvR = (*m_inversetransformationmatrix);
            //transposed inverse transformation matrix
            m_inversetransposedtransformationmatrix = MemoryManager<DNekMat>::AllocateSharedPtr(nCoeffs, 
                                                                                 nCoeffs,zero, storage);
            DNekMat &InvRT = (*m_inversetransposedtransformationmatrix);

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

            //vertex-edge/face
            for (i=0; i<nVerts; ++i)
            {
                for(j=0; j<nedgemodestotal; ++j)
                {
                    InvR.SetValue(vertModeLocation[i],edgemodearray[j],
                             -R(vertModeLocation[i],edgemodearray[j]));
                    InvRT.SetValue(edgemodearray[j],vertModeLocation[i],
                              -R(vertModeLocation[i],edgemodearray[j]));
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    InvR.SetValue(vertModeLocation[i],facemodearray[j],
                             -R(vertModeLocation[i],facemodearray[j]));
                    InvRT.SetValue(facemodearray[j],vertModeLocation[i],
                              -R(vertModeLocation[i],facemodearray[j]));
                    for(n=0; n<nedgemodestotal; ++n)
                    {
                        MatrixValue=InvR.GetValue(vertModeLocation[i],facemodearray[j])
                                               +R(vertModeLocation[i],edgemodearray[n])
                                                 *R(edgemodearray[n],facemodearray[j]);
                        InvR.SetValue(vertModeLocation[i],facemodearray[j],MatrixValue);
                        InvRT.SetValue(facemodearray[j],vertModeLocation[i],MatrixValue);
                    }
                }
            }

            //edge-face contributions
            for (i=0; i<nedgemodestotal; ++i)
            {
                for(j=0; j<nfacemodestotal; ++j)
                {
                    InvR.SetValue(edgemodearray[i],facemodearray[j],-R(edgemodearray[i],facemodearray[j]));
                    InvRT.SetValue(facemodearray[j],edgemodearray[i],-R(edgemodearray[i],facemodearray[j]));
                }
            }

            for (i = 0; i < nCoeffs; ++i)
            {
                InvR.SetValue(i,i,1.0);
                InvRT.SetValue(i,i,1.0);
            }
        }


        /**
         *
         */
        void PreconditionerLowEnergy::VertexEdgeFaceMatrix()
	{
            int i,j, eid, fid;
            int nCoeffs=vExp->NumBndryCoeffs();
            NekDouble MatrixValue;
            DNekMat &R = (*m_transformationmatrix);
            // Define storage for vertex transpose matrix and zero all entries
            MatrixStorage storage = eFULL;
            DNekMatSharedPtr m_om = MemoryManager<DNekMat>::AllocateSharedPtr(nCoeffs, nCoeffs, 0, storage);
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
                    MatrixValue=R(vertModeLocation[i],vertModeLocation[j]);

                    OM.SetValue(i,j,MatrixValue);
                }
 
                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=R(vertModeLocation[i],edgemodearray[j]);
                    OM.SetValue(i,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=R(vertModeLocation[i],facemodearray[j]);
                    OM.SetValue(i,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }

            //edge-vertex/edge/face
            for (i=0; i<nedgemodestotal; ++i)
            {
                for(j=0; j<nVerts; ++j)
                {
                    MatrixValue=R(edgemodearray[i],vertModeLocation[j]);
                    OM.SetValue(i+nVerts,j,MatrixValue);
                }

                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=R(edgemodearray[i],edgemodearray[j]);
                    OM.SetValue(i+nVerts,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=R(edgemodearray[i],facemodearray[j]);
                    OM.SetValue(i+nVerts,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }

            //face-vertex/edge/face
            for (i=0; i<nfacemodestotal; ++i)
            {
                for(j=0; j<nVerts; ++j)
                {
                    MatrixValue=R(facemodearray[i],vertModeLocation[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j,MatrixValue);
                }
 
                for(j=0; j<nedgemodestotal; ++j)
                {
                    MatrixValue=R(facemodearray[i],edgemodearray[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j+nVerts,MatrixValue);
                }

                for(j=0; j<nfacemodestotal; ++j)
                {
                    MatrixValue=R(facemodearray[i],facemodearray[j]);
                    OM.SetValue(i+nVerts+nedgemodestotal,j+nVerts+nedgemodestotal,MatrixValue);
                }
            }
        }

       /**
	 * \brief Construct the low energy preconditioner from \f$\mathbf{S}_{2}\f$
	 *
	 *\f[\mathbf{M}^{-1}=\left[\begin{array}{ccc} Diag[(\mathbf{S_{2}})_{vv}] & & \\
	 *  & (\mathbf{S}_{2})_{eb} & \\
	 *  &  & (\mathbf{S}_{2})_{fb} \end{array}\right] \f]
	 *
	 * where \f$\mathbf{R}\f$ is the transformation matrix and \f$\mathbf{S}_{2}\f$
	 * the Schur complement of the modified basis, given by
	 *
	 * \f[\mathbf{S}_{2}=\mathbf{R}\mathbf{S}_{1}\mathbf{R}^{T}\f]
	 *
	 * where \f$\mathbf{S}_{1}\f$ is the local schur complement matrix for each
	 * element.
	 */
        void PreconditionerLowEnergy::LowEnergyPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            StdRegions::StdExpansionSharedPtr locExpansion;

            int nRow;
            int nVerts, nEdges,nFaces, eid, fid, n, cnt, nedgemodes, nfacemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2, m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;
            NekDouble globalMatrixValue;

            MatrixStorage storage = eFULL;
            MatrixStorage vertstorage = eDIAGONAL;
            DNekMat &R = (*m_transformationmatrix);
            DNekMat &RT = (*m_transposedtransformationmatrix);

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;

            DNekMatSharedPtr    m_RS;
            DNekMatSharedPtr    m_RSRT;

            DNekMatSharedPtr m_VertBlk;
            DNekMatSharedPtr m_EdgeBlk;
            DNekMatSharedPtr m_FaceBlk;

            nRow=vExp->NumBndryCoeffs();

            m_RS = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &RS = (*m_RS);
            m_RSRT = MemoryManager<DNekMat>::AllocateSharedPtr(nRow, nRow, zero, storage);
            DNekMat &RSRT = (*m_RSRT);

            nVerts=vExp->GetGeom()->GetNumVerts();
            nEdges=vExp->GetGeom()->GetNumEdges();
            nFaces=vExp->GetGeom()->GetNumFaces();

            int nDirBnd    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
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

                //Calculate S*trans(R)  (note R is already transposed)
                RS=R*S;

                //Calculate R*S*trans(R)
                RSRT=RS*RT;

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
                                                      + sign1*sign2*RSRT(eMap1,eMap2);

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
                                                      + sign1*sign2*RSRT(fMap1,fMap2);

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
	 * \brief Construct the Block preconditioner from \f$\mathbf{S}_{1}\f$
	 *
	 * \f[\mathbf{M}^{-1}=\left[\begin{array}{ccc} Diag[(\mathbf{S_{1}})_{vv}] & & \\
	 *  & (\mathbf{S}_{1})_{eb} & \\
	 *  &  & (\mathbf{S}_{1})_{fb} \end{array}\right]\f]
	 *
	 * where \f$\mathbf{R}\f$ is the transformation matrix and \f$\mathbf{S}_{1}\f$
	 * is the Schur complement of each element.
	 *
	 */
        void PreconditionerLowEnergy::BlockPreconditioner()
        {
            boost::shared_ptr<MultiRegions::ExpList> expList=((m_linsys.lock())->GetLocMat()).lock();
            StdRegions::StdExpansionSharedPtr locExpansion;

            int nRow;
            int nVerts, nEdges,nFaces, eid, fid, n, cnt, nedgemodes, nfacemodes;
            NekDouble zero = 0.0;

            int vMap1, vMap2, sign1, sign2, m, v, eMap1, eMap2, fMap1, fMap2;
            int offset, globalrow, globalcol;
            NekDouble globalMatrixValue;

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

            int nDirBnd    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
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

        /**
         * \brief Get the transformation matrix \f$\mathbf{R}\f$
         */
        const DNekMatSharedPtr& PreconditionerLowEnergy::v_GetTransformationMatrix() const
	{
	    return m_transformationmatrix;
	}

        /**
         * \brief Get the transposed transformation matrix \f$\mathbf{R}^{T}\f$
         */
        const DNekMatSharedPtr& PreconditionerLowEnergy::v_GetTransposedTransformationMatrix() const
	{
	    return m_transposedtransformationmatrix;
	}

        /**
         * \brief Get the inverse transformation matrix \f$\mathbf{R}^{-1}\f$
         */
        const DNekMatSharedPtr& PreconditionerLowEnergy::v_GetInverseTransformationMatrix() const
	{
	    return m_inversetransformationmatrix;
	}

        /**
         * \brief Get the inverse of the transposed transformation matrix \f$\mathbf{R}^{-T}\f$
         */
        const DNekMatSharedPtr& PreconditionerLowEnergy::v_GetInverseTransposedTransformationMatrix() const
	{
	    return m_inversetransposedtransformationmatrix;
	}

    }
}






