///////////////////////////////////////////////////////////////////////////////
//
// File  NekNonlinSysNewton.h
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
// Description: NekNonlinSysNewton header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_NEWTON_H

#include <LibUtilities/LinearAlgebra/NekNonlinSys.h>

namespace Nektar
{
namespace LibUtilities
{

class NekNonlinSysNewton;

class NekNonlinSysNewton : public NekNonlinSys
{
public:
    /// Constructor for full direct matrix solve.
    friend class MemoryManager<NekNonlinSysNewton>;

    LIB_UTILITIES_EXPORT static NekNonlinSysSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekNonlinSysSharedPtr p =
            MemoryManager<NekNonlinSysNewton>::AllocateSharedPtr(pSession,
                vComm, nDimen, pKey);
        p->InitObject();
        return p;
    }

    static std::string className;
    LIB_UTILITIES_EXPORT NekNonlinSysNewton(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nscale,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekNonlinSysNewton();

protected:

    NekDouble m_SysResNorm0;
    NekDouble m_SysResNorm;

    int m_InexactNewtonForcing = 0;
    NekDouble   m_ForcingGama  = 1.0;
    NekDouble   m_ForcingAlpha = 0.5 * (1.0 + sqrt(5));

    virtual void v_InitObject();

    virtual int v_SolveSystem(
        const int nGlobal,
        const Array<OneD, const NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput, const int nDir,
        const NekDouble tol, const NekDouble factor);

    virtual bool v_ConvergenceCheck(
        const int nIteration, const Array<OneD, const NekDouble> &Residual,
        const NekDouble tol);
    void CalcInexactNewtonForcing(
        const int       &k,
        NekDouble       &resnormOld,
        const NekDouble &resnorm,
        NekDouble       &forcing);

    virtual void v_SetupNekNonlinSystem(
        const int nGlobal, const Array<OneD, const NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput, const int nDir);


private:
};
} // namespace LibUtilities
} // namespace Nektar

#endif
