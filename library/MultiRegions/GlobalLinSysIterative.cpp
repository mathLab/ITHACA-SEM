///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysIterative.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: GlobalLinSysIterative definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysIterative.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        std::string GlobalLinSysIterative::IteratSolverlookupIds[2] =
        {
            LibUtilities::SessionReader::RegisterEnumValue(
                "LinSysIterSovler", "ConjugateGradient", MultiRegions::
                eConjugateGradient),
            LibUtilities::SessionReader::RegisterEnumValue(
                "LinSysIterSovler", "GMRES", MultiRegions::eGMRES),
        };

        std::string GlobalLinSysIterative::IteratSolverdef =
            LibUtilities::SessionReader::RegisterDefaultSolverInfo(
                "LinSysIterSovler", "ConjugateGradient");

        /**
         * @class GlobalLinSysIterative
         *
         * Solves a linear system using iterative methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysIterative::GlobalLinSysIterative(
                const GlobalLinSysKey &pKey,
                const std::weak_ptr<ExpList> &pExpList,
                const std::shared_ptr<AssemblyMap>
                &pLocToGloMap)
                : GlobalLinSys(pKey, pExpList, pLocToGloMap),
                  m_rhs_magnitude(NekConstants::kNekUnsetDouble),
                  m_rhs_mag_sm(0.9),
                  m_precon(NullPreconditionerSharedPtr),
                  m_totalIterations(0),
                  m_useProjection(false),
                  m_numPrevSols(0)
        {
            m_tolerance = pLocToGloMap->GetIterativeTolerance();
            m_maxiter   = pLocToGloMap->GetMaxIterations();

            LibUtilities::CommSharedPtr vComm = m_expList.lock()->GetComm()->GetRowComm();
            m_root    = (vComm->GetRank())? false : true;

            int successiveRHS;

            if((successiveRHS = pLocToGloMap->GetSuccessiveRHS()))
            {
                m_prevLinSol.set_capacity(successiveRHS);
                m_useProjection = true;
            }
            else
            {
                m_useProjection = false;
            }
        }

        GlobalLinSysIterative::~GlobalLinSysIterative()
        {
        }

        /**
         *
         */
        void GlobalLinSysIterative::v_SolveLinearSystem(
                    const int nGlobal,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const AssemblyMapSharedPtr &plocToGloMap,
                    const int nDir)
        {
            if (!m_linsol)
            {
                LibUtilities::CommSharedPtr v_Comm = 
                m_expList.lock()->GetComm()->GetRowComm();
                LibUtilities::SessionReaderSharedPtr pSession =
                m_expList.lock()->GetSession();

                std::vector<std::string>  variables(1);
                variables[0] =  pSession->GetVariable(0);
                string variable = variables[0];

                std::string LinSysIterSovlerType = "ConjugateGradient";
                if (pSession->DefinesGlobalSysSolnInfo(variable,
                                                    "LinSysIterSovler"))
                {
                    LinSysIterSovlerType = pSession->GetGlobalSysSolnInfo(
                                          variable,"LinSysIterSovler");
                }
                else
                {
                    if (pSession->DefinesSolverInfo("LinSysIterSovler"))
                    {
                        LinSysIterSovlerType = pSession->GetSolverInfo(
                        "LinSysIterSovler");
                    }
                }
                
                // Check such a module exists for this equation.
                ASSERTL0(LibUtilities::GetNekLinSysIterFactory().
                ModuleExists(LinSysIterSovlerType),
                    "NekLinSysIter '" + LinSysIterSovlerType + 
                    "' is not defined.\n");
                m_linsol = LibUtilities::GetNekLinSysIterFactory().
                           CreateInstance(LinSysIterSovlerType, pSession,
                                          v_Comm, nGlobal - nDir);

                m_NekSysOp.DefineNekSysLhsEval(
                    &GlobalLinSysIterative::DoMatrixMultiplyFlag, this);
                m_NekSysOp.DefineNekSysPrecond(
                    &GlobalLinSysIterative::DoPreconditionerFlag, this);
                m_linsol->SetSysOperators(m_NekSysOp);
                    v_UniqueMap();
                m_linsol->setUniversalUniqueMap(m_map);
            }
            
            if (!m_precon)
            {
                m_precon = CreatePrecon(plocToGloMap);
                m_precon->BuildPreconditioner();
            }

            m_linsol->setRhsMagnitude(m_rhs_magnitude);
            int ntmpGMRESIts = m_linsol->SolveSystem(
            nGlobal, pInput, pOutput, nDir, m_tolerance);
            boost::ignore_unused(ntmpGMRESIts);
        }
        void GlobalLinSysIterative::Set_Rhs_Magnitude(
            const NekVector<NekDouble> &pIn)
        {
            Array<OneD, NekDouble> vExchange(1, 0.0);
            if (m_map.size() > 0)
            {
                vExchange[0] = Vmath::Dot2(pIn.GetDimension(),
                                        &pIn[0],&pIn[0],&m_map[0]);
            }

            m_expList.lock()->GetComm()->GetRowComm()->AllReduce(
                vExchange, Nektar::LibUtilities::ReduceSum);

            // To ensure that very different rhs values are not being
            // used in subsequent solvers such as the velocit solve in
            // INC NS. If this works we then need to work out a better
            // way to control this.
            NekDouble new_rhs_mag = (vExchange[0] > 1e-6)? vExchange[0] : 1.0;

            if(m_rhs_magnitude == NekConstants::kNekUnsetDouble)
            {
                m_rhs_magnitude = new_rhs_mag;
            }
            else
            {
                m_rhs_magnitude = (m_rhs_mag_sm*(m_rhs_magnitude) +
                                   (1.0-m_rhs_mag_sm)*new_rhs_mag);
            }
        }

    }
}
