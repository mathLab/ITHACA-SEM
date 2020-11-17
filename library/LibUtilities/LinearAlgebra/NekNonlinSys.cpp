///////////////////////////////////////////////////////////////////////////////
//
// File:  NekNonlinSys.cpp
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
// Description:  NekNonlinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSys.h>

using namespace std;

namespace Nektar
{
namespace LibUtilities
{
/**
 * @class  NekNonlinSys
 *
 * Solves a nonlinear system using iterative methods.
 */
NekNonlinSysFactory &GetNekNonlinSysFactory()
{
    static NekNonlinSysFactory instance;
    return instance;
}

NekNonlinSys::NekNonlinSys(const LibUtilities::SessionReaderSharedPtr &pSession,
                           const LibUtilities::CommSharedPtr &vComm,
                           const int nDimen,
                           const NekSysKey &pKey)
    : NekSys(pSession, vComm, nDimen, pKey)
{
    std::vector<std::string> variables(1);
    variables[0]    = pSession->GetVariable(0);
    string variable = variables[0];

    if (pSession->DefinesGlobalSysSolnInfo(variable, "NekNonlinSysTolerance"))
    {
        m_tolerance = boost::lexical_cast<NekDouble>(
            pSession->GetGlobalSysSolnInfo(variable, "NekNonlinSysTolerance")
                .c_str());
    }
    else
    {
        pSession->LoadParameter("NekNonlinSysTolerance", m_tolerance,
            pKey.m_DefaultNekNonlinSysTolerance);
    }

    if (pSession->DefinesGlobalSysSolnInfo(variable,
                                           "NekNonlinSysMaxIterations"))
    {
        m_maxiter = boost::lexical_cast<int>(
            pSession
                ->GetGlobalSysSolnInfo(variable, "NekNonlinSysMaxIterations")
                .c_str());
    }
    else
    {
        pSession->LoadParameter("NekNonlinSysMaxIterations", m_maxiter, 
            pKey.m_DefaultNekNonlinSysMaxIterations);
    }

    if (pSession->DefinesGlobalSysSolnInfo(variable, "NonlinIterTolRelativeL2"))
    {
        m_NonlinIterTolRelativeL2 = boost::lexical_cast<int>(
            pSession->GetGlobalSysSolnInfo(variable, "NonlinIterTolRelativeL2")
                .c_str());
    }
    else
    {
        pSession->LoadParameter("NonlinIterTolRelativeL2",
            m_NonlinIterTolRelativeL2, pKey.m_DefaultNonlinIterTolRelativeL2);
    }

    if (pSession->DefinesGlobalSysSolnInfo(variable,
                                           "LinSysRelativeTolInNonlin"))
    {
        m_LinSysRelativeTolInNonlin = boost::lexical_cast<int>(
            pSession
                ->GetGlobalSysSolnInfo(variable, "LinSysRelativeTolInNonlin")
                .c_str());
    }
    else
    {
        pSession->LoadParameter("LinSysRelativeTolInNonlin",
            m_LinSysRelativeTolInNonlin,
            pKey.m_DefaultLinSysRelativeTolInNonlin);
    }

    // cout << " m_LinSysRelativeTolInNonlin = " << m_LinSysRelativeTolInNonlin << endl;

    m_LinSysIterSovlerType = pKey.m_DefaultLinSysIterSovlerTypeInNonlin;
    if (pSession->DefinesGlobalSysSolnInfo(variable, 
            "LinSysIterSovlerTypeInNonlin"))
    {
        m_LinSysIterSovlerType =
            pSession->GetGlobalSysSolnInfo(variable, 
            "LinSysIterSovlerTypeInNonlin");
    }
    else
    {
        if (pSession->DefinesSolverInfo("LinSysIterSovlerTypeInNonlin"))
        {
            m_LinSysIterSovlerType =
                pSession->GetSolverInfo("LinSysIterSovlerTypeInNonlin");
        }
    }

    ASSERTL0(LibUtilities::GetNekLinSysIterFactory().ModuleExists(
                 m_LinSysIterSovlerType),
             "NekLinSysIter '" + m_LinSysIterSovlerType +
                 "' is not defined.\n");
                 
    m_linsol = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
        m_LinSysIterSovlerType, pSession, m_Comm, m_SysDimen, pKey);
}

void NekNonlinSys::v_InitObject()
{
    NekSys::v_InitObject();
}

NekNonlinSys::~NekNonlinSys()
{
}

void NekNonlinSys::v_SetupNekNonlinSystem(
    const int nGlobal, const TensorOfArray1D<NekDouble> &pInput,
    const TensorOfArray1D<NekDouble> &pSource,
    const int nDir)
{
    boost::ignore_unused(nGlobal, pInput, pSource, nDir);
    ASSERTL0(false, "v_SetupNekNonlinSystem not defined");
}

} // namespace LibUtilities
} // namespace Nektar
