///////////////////////////////////////////////////////////////////////////////
//
// File:  NekNonlinSysNewton.cpp
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
// Description:  NekNonlinSysNewton definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysNewton.h>
#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{  
    namespace LibUtilities
    {     
        /**
         * @class  NekNonlinSysNewton
         *
         * Solves a nonlinear system using iterative methods.
         */
		 string NekNonlinSysNewton::className =
        LibUtilities::GetNekNonlinSysFactory().RegisterCreatorFunction(
            "Newton", NekNonlinSysNewton::create,
            "NekNonlinSysNewton solver.");

        NekNonlinSysNewton::NekNonlinSysNewton(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nscale)
            : NekNonlinSys(pSession, vComm, nscale)
        {
            std::vector<std::string>  variables(1);
            variables[0] =  pSession->GetVariable(0);
            string variable = variables[0];

            if (pSession->DefinesGlobalSysSolnInfo(variable,
                                                "NonlinIteTolRelativeL2"))
            {
                m_NonlinIteTolRelativeL2 = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "NonlinIteTolRelativeL2").c_str());
            }
            else
            {
                pSession->LoadParameter("NonlinIteTolRelativeL2",
                                        m_NonlinIteTolRelativeL2,
                                        1.0E-10);
            }

            if (pSession->DefinesGlobalSysSolnInfo(variable,
                                                   "NonlinIteTolLinRelatTol"))
            {
                m_NonlinIteTolLinRelatTol = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "NonlinIteTolLinRelatTol").c_str());
            }
            else
            {
                pSession->LoadParameter("NonlinIteTolLinRelatTol",
                                        m_NonlinIteTolLinRelatTol,
                                        m_NonlinIteTolRelativeL2);
            }

            m_LinSysIterSovlerType = "GMRES";
            if (pSession->DefinesGlobalSysSolnInfo(variable,
                                                "LinSysIterSovler"))
            {
                m_LinSysIterSovlerType = pSession->GetGlobalSysSolnInfo(
                                        variable,
                                        "LinSysIterSovler");
            }
            else
            {
                if (pSession->DefinesSolverInfo("LinSysIterSovler"))
                {
                    m_LinSysIterSovlerType = pSession->GetSolverInfo(
                    "LinSysIterSovler");
                }
            }
        }

        void NekNonlinSysNewton::v_InitObject()
        {
            NekSys::v_InitObject();
            m_Residual = Array<OneD, NekDouble> (m_SysDimen, 0.0);
            m_DeltSltn = Array<OneD, NekDouble> (m_SysDimen, 0.0);

            ASSERTL0(LibUtilities::GetNekLinSysIterFactory().
                    ModuleExists(m_LinSysIterSovlerType),
                    "NekLinSysIter '" + m_LinSysIterSovlerType + 
                    "' is not defined.\n");

            m_linsol = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
                                m_LinSysIterSovlerType, 
                                m_session.lock(), m_Comm, m_SysDimen);

        }

        NekNonlinSysNewton::~NekNonlinSysNewton()
        {
        }

        /**
         *
         */
        int NekNonlinSysNewton::v_SolveSystem(
            const int                           nGlobal,
            const Array<OneD, const NekDouble>  &pInput,
            Array<OneD,      NekDouble>         &pOutput,
            const int                           nDir,
            const NekDouble                     tol,
            const NekDouble                     factor)
        {
            int nwidthcolm = 12;
            m_linsol->setSysOperators(m_operator);

            ASSERTL0(0 == nDir,"0 != nDir not tested");
            ASSERTL0(m_SysDimen == nGlobal, "m_SysDimen!=nGlobal");

            boost::ignore_unused(factor);

            int ntotal = nGlobal - nDir;
            int NtotLinSysIts = 0;
            
            int NttlNonlinIte    = 0;
            m_converged = false;

            m_Solution = pOutput;
            Vmath::Vcopy(ntotal, pInput, 1, m_Solution, 1);
            for (int k = 0; k < m_maxiter; ++k)
            {
                m_operator.DoNekSysRhsEval(m_Solution, m_Residual);
                
                m_converged = v_ConvergenceCheck(k,m_Residual, tol);
                if (m_converged) break;

                NekDouble   LinSysTol = m_NonlinIteTolLinRelatTol * 
                                        sqrt(m_SysResNorm);
                int ntmpGMRESIts =  m_linsol->SolveSystem(ntotal, 
                                                          m_Residual,
                                                          m_DeltSltn,
                                                          0, 
                                                          LinSysTol);
                NtotLinSysIts   +=  ntmpGMRESIts;
                Vmath::Vsub(ntotal, m_Solution, 1, 
                            m_DeltSltn, 1, m_Solution, 1);
                NttlNonlinIte++;
            }

            if ((m_root || (!m_converged)) && m_verbose)
            {
                WARNINGL0(m_converged,
                    "     # Nonlinear solver not converge in DoImplicitSolve");
               cout << right << scientific << setw(nwidthcolm) 
                    << setprecision(nwidthcolm-6)
                    <<"     * Newton-Its converged (RES="
                    << sqrt(m_SysResNorm) << " Res/(DtRHS): "
                    << sqrt(m_SysResNorm / m_SysResNorm0)
                    << " with " << setw(3) << NttlNonlinIte
                    << " Non-Its)" <<endl;
            }
            return NttlNonlinIte;
        }

        bool NekNonlinSysNewton::v_ConvergenceCheck(
                const int                           nIteration,
                const Array<OneD, const NekDouble>  &Residual,
                const NekDouble                     tol         )
        {
            bool converged = false;
            // int nwidthcolm = 12;
            // NekDouble resmaxm  = 0.0;
            NekDouble resratio = 1.0;
            int ntotal = Residual.size();

            m_SysResNorm = Vmath::Dot(ntotal,Residual,Residual);
            m_Comm->AllReduce(m_SysResNorm, Nektar::LibUtilities::ReduceSum);

            if (0 == nIteration)
            {
                m_SysResNorm0 = m_SysResNorm;
                resratio = 1.0;
            }
            else
            {
                resratio = m_SysResNorm / m_SysResNorm0;
            }

            if (resratio < (m_NonlinIteTolRelativeL2 * m_NonlinIteTolRelativeL2) 
                || m_SysResNorm < tol)
            {
                converged = true;
            }
            
            return converged;
        }
    }
}

