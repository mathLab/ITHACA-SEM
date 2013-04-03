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
#include <MultiRegions/PreconditionerLinear.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysIterativeStaticCond.h>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/GlobalLinSysXxtFull.h>
#include <LocalRegions/MatrixKey.h>
#include <math.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */

        string PreconditionerLinear::className1
                = GetPreconFactory().RegisterCreatorFunction(
                    "Linear",
                    PreconditionerLinear::create,
                    "Linear space inverse Preconditioning");
 
       /**
         * @class PreconditionerLinear
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */
        
        PreconditionerLinear::PreconditionerLinear(
            const boost::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }
        
        void PreconditionerLinear::v_InitObject()
        {
            //CreateMultiplicityMap();
        }

        void PreconditionerLinear::v_BuildPreconditioner()
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(solvertype)
            {
                case MultiRegions::eIterativeFull:
                {
                    InverseLinearSpacePreconditioner();
                }
                break;
                case MultiRegions::eIterativeStaticCond:
                {
#if 1
                    boost::shared_ptr<MultiRegions::ExpList> 
                        expList=((m_linsys.lock())->GetLocMat()).lock();
                    m_vertLocToGloMap = m_locToGloMap->XxtLinearSpaceMap(*expList);

                    // Generate XXT system. 
                    GlobalLinSysKey preconKey(StdRegions::ePreconLinearSpace,
                                              m_vertLocToGloMap,
                                              (m_linsys.lock())->GetKey().GetConstFactors());


                    m_vertLinsys = MemoryManager<GlobalLinSysXxtFull>::
                        AllocateSharedPtr(preconKey,expList,m_vertLocToGloMap);
#else
                    StaticCondInverseLinearSpacePreconditioner();
#endif
                }
                break;
                default:
                    ASSERTL0(0,"This type of preconditioning is not implemented for this solver");
                    break;
            }
	}

        /**
	 * \brief Extracts the entries from a statically condensed matrix
	 * corresponding to the vertex modes.
	 *
	 * This function extracts a static condensed matrix from the nth 
	 * expansion and returns a matrix containing
	 * only the vertex mode contributions i.e the linear finite element
	 * space.
	 */         
        void PreconditionerLinear::
        StaticCondInverseLinearSpacePreconditioner()
	{
            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            GlobalLinSysKey m_linSysKey=(m_linsys.lock())->GetKey();
            StdRegions::VarCoeffMap vVarCoeffMap;
            StdRegions::StdExpansionSharedPtr locExpansion;

            int cnt, n, nel, nVerts;

            int n_exp=expList->GetNumElmts();

            //Allocate Linear space matrix
            DNekScalBlkMatSharedPtr vertMatBlk;
            DNekScalMatSharedPtr vertMat;


           for(cnt=n=0; n < n_exp; ++n)
           {
               nel = expList->GetOffset_Elmt_Id(n);
               
               locExpansion = expList->GetExp(nel);

               //Get total of vertices
               nVerts=locExpansion->GetNverts();

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
               LocalRegions::MatrixKey vertexspace(
                   StdRegions::ePreconLinearSpace,
                   locExpansion->DetShapeType(),
                   *locExpansion, m_linSysKey.GetConstFactors(),
                   vVarCoeffMap);

               vertMatBlk = locExpansion->GetLocStaticCondMatrix(vertexspace);
               vertMat=vertMatBlk->GetBlock(0,0);
           }

        }


        /**
         * \brief Inverse of the linear space
	 *
	 * Extracts the linear space and inverts it.
	 *
         */         
        void PreconditionerLinear::InverseLinearSpacePreconditioner()
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
         *
         */
        void PreconditionerLinear::v_DoPreconditioner(
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
                    DNekMat &M = (*m_preconditioner);
                    
                    NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                    NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                    z = M * r;
                }
                break;
            case MultiRegions::eXxtFullMatrix:
                {
                }
                break;
            case MultiRegions::eIterativeStaticCond:
                {
#if 1
                    int i,val;
                    int nloc = m_vertLocToGloMap->GetNumLocalCoeffs();
                    int nglo = m_vertLocToGloMap->GetNumGlobalCoeffs();
                    // mapping from full space to vertices
                    Array<OneD, int> LocToGloBnd = m_vertLocToGloMap->GetLocalToGlobalBndMap();

                    // Global to local for linear solver (different from above)
                    Array<OneD, int> LocToGlo = m_vertLocToGloMap->GetLocalToGlobalMap();

                    // number of Dir coeffs in full system. 
                    int nDirFull = m_locToGloMap->GetNumGlobalDirBndCoeffs();

                    Array<OneD,NekDouble> In(nglo,0.0);
                    Array<OneD,NekDouble> Out(nglo,0.0);
                    
                    // Gather rhs
                    for(i = 0; i < nloc; ++i)
                    {
                        val = LocToGloBnd[i];
                        if(val >= nDirFull)
                        {
                            In[LocToGlo[i]] = pInput[val-nDirFull];
                        }
                    }

                    // Do solve without enforcing any boundary conditions. 
                    m_vertLinsys->SolveLinearSystem(m_vertLocToGloMap->GetNumLocalCoeffs(),
                                                In,Out,m_vertLocToGloMap);
                    

                    //Copy input to output as a unit preconditioner on any other value
                    Vmath::Vcopy(pInput.num_elements(),pInput,1,pOutput,1);


                    // Scatter back soln from linear solve
                    for(i = 0; i < nloc; ++i)
                    {
                        val = LocToGloBnd[i];
                        if(val >= nDirFull)
                        {
                            pOutput[val-nDirFull] = Out[LocToGlo[i]];
                        }
                    }
#else
                    int nDir    = m_locToGloMap->GetNumGlobalDirBndCoeffs();
                    int nGlobal = m_locToGloMap->GetNumGlobalBndCoeffs();
                    int nNonDir = nGlobal-nDir;
                    DNekMat &M = (*m_preconditioner);
                    
                    NekVector<NekDouble> r(nNonDir,pInput,eWrapper);
                    NekVector<NekDouble> z(nNonDir,pOutput,eWrapper);
                    z = M * r;
#endif
                }
                break;
                default:
                    ASSERTL0(0,"Unsupported solver type");
                    break;
	    }
        }


        /**
         * Create the inverse multiplicity map.
         */
        void PreconditionerLinear::CreateMultiplicityMap(void)
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






