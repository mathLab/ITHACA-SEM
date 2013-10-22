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
                    "FullLinearSpace",
                    PreconditionerLinear::create,
                    "Full Linear space inverse Preconditioning");
 
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
        }

        void PreconditionerLinear::v_BuildPreconditioner()
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            if(solvertype != eIterativeStaticCond)
            {
                ASSERTL0(0,"This type of preconditioning is not implemented for this solver");
            }

            boost::shared_ptr<MultiRegions::ExpList> 
                expList=((m_linsys.lock())->GetLocMat()).lock();
            m_vertLocToGloMap = m_locToGloMap->XxtLinearSpaceMap(*expList);

            // Generate XXT system. 
            GlobalLinSysKey preconKey(StdRegions::ePreconLinearSpace,
                                      m_vertLocToGloMap,
                                      (m_linsys.lock())->GetKey().GetConstFactors());

            m_vertLinsys = MemoryManager<GlobalLinSysXxtFull>::
                AllocateSharedPtr(preconKey,expList,m_vertLocToGloMap);

	}

        /**
         *
         */
        void PreconditionerLinear::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput)
        {
            v_DoPreconditionerWithNonVertOutput(pInput,pOutput,NullNekDouble1DArray);
        }

        /**
         *
         */
        void PreconditionerLinear::v_DoPreconditionerWithNonVertOutput(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput,
            const Array<OneD, NekDouble>& pNonVertOutput)
        {
            GlobalSysSolnType solvertype=m_locToGloMap->GetGlobalSysSolnType();
            switch(solvertype)
            {
                case MultiRegions::eIterativeStaticCond:
                {
                    int i,val;
                    int nloc = m_vertLocToGloMap->GetNumLocalCoeffs();
                    int nglo = m_vertLocToGloMap->GetNumGlobalCoeffs();
                    // mapping from full space to vertices
                    Array<OneD, int> LocToGloBnd = m_vertLocToGloMap->GetLocalToGlobalBndMap();

                    // Global to local for linear solver (different from above)
                    Array<OneD, int> LocToGlo = m_vertLocToGloMap->GetLocalToGlobalMap();
                    
                    // number of Dir coeffs in from full problem
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
                    
                    
                    if(pNonVertOutput != NullNekDouble1DArray)
                    {
                        ASSERTL1(pNonVertOutput.num_elements() >= pOutput.num_elements(),"Non Vert output is not of sufficient length");
                        Vmath::Vcopy(pOutput.num_elements(),pNonVertOutput,1,pOutput,1);
                    }
                    else
                    {
                        //Copy input to output as a unit preconditioner on
                        //any other value
                        Vmath::Vcopy(pInput.num_elements(),pInput,1,pOutput,1);
                    }

                    // Scatter back soln from linear solve
                    for(i = 0; i < nloc; ++i)
                    {
                        val = LocToGloBnd[i];
                        if(val >= nDirFull)
                        {
                            pOutput[val-nDirFull] = Out[LocToGlo[i]];
                        }
                    }

                }
                break;
            default:
                ASSERTL0(0,"Unsupported solver type");
                break;
	    }
        }
    }
}






