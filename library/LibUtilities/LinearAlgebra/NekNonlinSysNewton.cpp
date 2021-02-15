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

#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/LinearAlgebra/NekNonlinSysNewton.h>

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
        "Newton", NekNonlinSysNewton::create, "NekNonlinSysNewton solver.");

NekNonlinSysNewton::NekNonlinSysNewton(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vComm, const int nscale,
    const NekSysKey &pKey)
    : NekNonlinSys(pSession, vComm, nscale, pKey)
{
    
}

void NekNonlinSysNewton::v_InitObject()
{
    NekSys::v_InitObject();
    m_Residual = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_DeltSltn = Array<OneD, NekDouble>(m_SysDimen, 0.0);
    m_SourceVec = Array<OneD, NekDouble>(m_SysDimen, 0.0);
}

NekNonlinSysNewton::~NekNonlinSysNewton()
{
}

/**
 *
 **/
int NekNonlinSysNewton::v_SolveSystem(
    const int nGlobal, const TensorOfArray1D<NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, const int nDir, const NekDouble tol,
    const NekDouble factor)
{
    boost::ignore_unused(factor);

    int nwidthcolm = 12;

    v_SetupNekNonlinSystem(nGlobal, pInput, pInput, nDir);

    m_Solution = pOutput;
    Vmath::Vcopy(nGlobal-nDir, pInput, 1, m_Solution, 1);

    int ntotal        = nGlobal - nDir;
    m_NtotLinSysIts = 0;

    int NttlNonlinIte = 0;
    m_converged       = false;

    NekDouble resnormOld = 0.0;
    for (int k = 0; k < m_maxiter; ++k)
    {
        m_converged = v_ConvergenceCheck(k, m_Residual, tol);
        if (m_converged)
            break;

        NekDouble LinSysRelativeIteTol;
        CalcInexactNewtonForcing(k, resnormOld, m_SysResNorm, 
            LinSysRelativeIteTol);
        
        resnormOld = m_SysResNorm;

        NekDouble LinSysTol = LinSysRelativeIteTol * sqrt(m_SysResNorm);
        int ntmpLinSysIts =
            m_linsol->SolveSystem(ntotal, m_Residual, m_DeltSltn, 0, LinSysTol);
        m_NtotLinSysIts += ntmpLinSysIts;
        Vmath::Vsub(ntotal, m_Solution, 1, m_DeltSltn, 1, m_Solution, 1);
        NttlNonlinIte++;
        m_operator.DoNekSysResEval(m_Solution, m_Residual);
    }

    if ( ((!m_converged) || m_verbose) && m_root)
    {
        WARNINGL0(m_converged,
                  "     # Nonlinear solver not converge in DoImplicitSolve");
        cout << right << scientific << setw(nwidthcolm)
             << setprecision(nwidthcolm - 6)
             << "     * Newton-Its converged (RES=" << sqrt(m_SysResNorm)
             << " Res/(DtRHS): " << sqrt(m_SysResNorm / m_SysResNorm0)
             << " with " << setw(3) << NttlNonlinIte << " Non-Its)" << endl;
    }

    m_ResidualUpdated = false;
    return NttlNonlinIte;
}

bool NekNonlinSysNewton::v_ConvergenceCheck(
    const int nIteration, const TensorOfArray1D<NekDouble> &Residual,
    const NekDouble tol)
{
    bool converged     = false;
    NekDouble resratio = 1.0;
    int ntotal         = Residual.size();

    m_SysResNorm = Vmath::Dot(ntotal, Residual, Residual);
    m_Comm->AllReduce(m_SysResNorm, Nektar::LibUtilities::ReduceSum);

    if (0 == nIteration)
    {
        m_SysResNorm0 = m_SysResNorm;
        resratio      = 1.0;
    }
    else
    {
        resratio = m_SysResNorm / m_SysResNorm0;
    }

    if (resratio < (m_NonlinIterTolRelativeL2 * m_NonlinIterTolRelativeL2) ||
        m_SysResNorm < tol)
    {
        converged = true;
    }

    return converged;
}

void NekNonlinSysNewton::CalcInexactNewtonForcing(
    const int       &k,
    NekDouble       &resnormOld,
    const NekDouble &resnorm,
    NekDouble       &forcing)
{
    if (0 == k)
    {
        forcing = m_LinSysRelativeTolInNonlin;
        resnormOld = resnorm;
    }
    else
    {
        switch(m_InexactNewtonForcing)
        {
        case 0:
            {
                forcing = m_LinSysRelativeTolInNonlin;
                break;
            }
        case 1: 
            {
                NekDouble tmpForc = m_ForcingGama * 
                    pow((resnorm / resnormOld), m_ForcingAlpha); 
                NekDouble tmp = m_ForcingGama * 
                    pow(forcing, m_ForcingAlpha);
                if (tmp > 0.1)
                {
                    forcing = min(m_LinSysRelativeTolInNonlin, 
                        max(tmp, tmpForc));
                }
                else
                {
                    forcing = min(m_LinSysRelativeTolInNonlin, tmpForc);
                }

                forcing = max(forcing,  1.0E-6);
                break;
            }
        }
    }
}

void NekNonlinSysNewton::v_SetupNekNonlinSystem(
    const int nGlobal, const TensorOfArray1D<NekDouble> &pInput,
    const TensorOfArray1D<NekDouble> &pSource,
    const int nDir)
{
    boost::ignore_unused(nGlobal, nDir);

    ASSERTL0(0 == nDir, "0 != nDir not tested");
    ASSERTL0(m_SysDimen == nGlobal, "m_SysDimen!=nGlobal");

    m_SourceVec = pSource;
    Vmath::Vcopy(nGlobal-nDir, pSource, 1, m_SourceVec, 1);

    if (!m_ResidualUpdated)
    {
        m_operator.DoNekSysResEval(pInput, m_Residual);
        m_ResidualUpdated = true;
    }
    m_linsol->SetSysOperators(m_operator);
}


} // namespace LibUtilities
} // namespace Nektar
