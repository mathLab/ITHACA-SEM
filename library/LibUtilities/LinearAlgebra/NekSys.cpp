///////////////////////////////////////////////////////////////////////////////
//
// File:  NekSys.cpp
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
// Description:  NekSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekSys.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{
/**
 * @class  NekSys
 *
 * Solves a nonlinear or linear system.
 */

NekSys::NekSys(const LibUtilities::SessionReaderSharedPtr &pSession,
               const LibUtilities::CommSharedPtr &vComm, const int nDimen,
               const NekSysKey &pKey)
{
    m_tolerance = pKey.m_Tolerance;
    m_verbose   = false;
    m_root      = false;
    m_Comm      = vComm;

    m_FlagWarnings = true;

    if (0 == m_Comm->GetRank())
    {
        m_root = true;
    }
    m_verbose = pSession->DefinesCmdLineArgument("verbose");

    m_converged = false;

    m_SysDimen = nDimen;
}

NekSys::~NekSys()
{
}

bool NekSys::v_ConvergenceCheck(const int nIteration,
                                const Array<OneD, const NekDouble> &Residual,
                                const NekDouble tol)
{
    bool converged = false;
    int ntotal     = Residual.size();
    boost::ignore_unused(nIteration);

    NekDouble SysResNorm = Vmath::Dot(ntotal, Residual, Residual);
    m_Comm->AllReduce(SysResNorm, Nektar::LibUtilities::ReduceSum);

    if (SysResNorm < tol)
    {
        converged = true;
    }
    return converged;
}

/**
 * Natural guess
**/
void NekSys::v_NekSysInitialGuess(
        const Array<OneD, const NekDouble> &pInput,
        Array<OneD, NekDouble> &pguess)
{
    size_t ndim = pInput.size();
    if(pguess.size() != ndim)
    {
        pguess = Array<OneD, NekDouble> {ndim};
    }

    Vmath::Vcopy(ndim, pInput, 1, pguess, 1);
}

} // namespace LibUtilities
} // namespace Nektar
