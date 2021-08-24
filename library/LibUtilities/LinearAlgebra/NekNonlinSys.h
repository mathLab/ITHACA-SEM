///////////////////////////////////////////////////////////////////////////////
//
// File  NekNonlinSys.h
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
// Description: NekNonlinSys header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_NONLINSYS_H

#include <LibUtilities/LinearAlgebra/NekSys.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>

namespace Nektar
{
namespace LibUtilities
{
class NekNonlinSys;

typedef std::shared_ptr<NekNonlinSys> NekNonlinSysSharedPtr;

typedef LibUtilities::NekFactory<std::string, NekNonlinSys,
                                 const LibUtilities::SessionReaderSharedPtr &,
                                 const LibUtilities::CommSharedPtr &, const int,
                                 const NekSysKey &>
    NekNonlinSysFactory;
LIB_UTILITIES_EXPORT NekNonlinSysFactory &GetNekNonlinSysFactory();

class NekNonlinSys : public NekSys
{
public:
    friend class MemoryManager<NekNonlinSys>;
    LIB_UTILITIES_EXPORT static NekNonlinSysSharedPtr CreateInstance(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey)
    {
        NekNonlinSysSharedPtr p =
            MemoryManager<NekNonlinSys>::AllocateSharedPtr(pSession, vComm,
                                                           nDimen, pKey);
        return p;
    }
    LIB_UTILITIES_EXPORT NekNonlinSys(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm, const int nDimen,
        const NekSysKey &pKey);
    LIB_UTILITIES_EXPORT ~NekNonlinSys();

    LIB_UTILITIES_EXPORT virtual void v_SetupNekNonlinSystem(
        const int nGlobal, const Array<OneD, const NekDouble> &pInput,
        const Array<OneD, const NekDouble> &pSource,
        const int nDir);

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefSolution()
        const
    {
        return m_Solution;
    }

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefResidual()
        const
    {
        return m_Residual;
    }

    LIB_UTILITIES_EXPORT const Array<OneD, const NekDouble> &GetRefSourceVec()
        const
    {
        return m_SourceVec;
    }

    LIB_UTILITIES_EXPORT void SetRefResidual(
        const Array<OneD, const NekDouble> &in)
    {
        ASSERTL0(in.size() == m_SysDimen, 
            "SetRefResidual dimension not correct");
        Vmath::Vcopy(m_SysDimen, in, 1, m_Residual, 1);
        
        m_ResidualUpdated = true;
    }

    LIB_UTILITIES_EXPORT void SetNekNonlinSysTolerance(const NekDouble in)
    {
        m_tolerance = in;
    }

    LIB_UTILITIES_EXPORT void SetNekNonlinSysMaxIterations(
        const unsigned int in)
    {
        m_maxiter = in;
    }

    LIB_UTILITIES_EXPORT const NekLinSysIterSharedPtr &
        GetLinSys()
    {
        return m_linsol;
    }

    LIB_UTILITIES_EXPORT void SetNonlinIterTolRelativeL2(const NekDouble in)
    {
        m_NonlinIterTolRelativeL2 = in;
    }

    LIB_UTILITIES_EXPORT void SetLinSysRelativeTolInNonlin(const NekDouble in)
    {
        m_LinSysRelativeTolInNonlin = in;
    }
    
    LIB_UTILITIES_EXPORT int GetNtotLinSysIts()
    {
        return m_NtotLinSysIts;
    }

protected:
    NekLinSysIterSharedPtr m_linsol;

    NekDouble m_NonlinIterTolRelativeL2;
    NekDouble m_LinSysRelativeTolInNonlin;

    std::string m_LinSysIterSolverType;

    int m_totalIterations = 0;
    int m_NtotLinSysIts = 0;


    Array<OneD, NekDouble> m_Solution;
    Array<OneD, NekDouble> m_Residual;
    Array<OneD, NekDouble> m_DeltSltn;
    Array<OneD, NekDouble> m_SourceVec;

    bool m_ResidualUpdated = false;

    virtual void v_InitObject();

private:
};
} // namespace LibUtilities
} // namespace Nektar
#endif
