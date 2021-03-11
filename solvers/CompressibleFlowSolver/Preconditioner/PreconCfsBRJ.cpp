///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfsBRJ.cpp
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
// Description:  PreconCfsBRJ definition
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Preconditioner/PreconCfsBRJ.h>

using namespace std;

namespace Nektar
{
    /**
     * @class  PreconCfsBRJ
     *
     * Solves a linear system using iterative methods.
     */
    std::string PreconCfsBRJ::className = GetPreconCfsOpFactory().
        RegisterCreatorFunction("PreconCfsBRJ",
            PreconCfsBRJ::create,
            "Block Relaxed Jacobi Preconditioner for CFS.");

    PreconCfsBRJ::PreconCfsBRJ(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm)
        : PreconCfsOp(pSession, vComm)
    {
        pSession->LoadParameter("JFNKPreconStep",      
            m_PreconItsStep, 7);
        pSession->LoadParameter("BRJRelaxParam",        
            m_BRJRelaxParam, 1.0);
    }

    void PreconCfsBRJ::v_InitObject()
    {
        PreconCfsOp::v_InitObject();
    }

    void PreconCfsBRJ::v_DoPreconCfs(
        const Array<OneD, NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput,
        const bool &flag)
    {
        NEKERROR(ErrorUtil::efatal, "v_DoPreconCfs not defined");
    }

    void PreconCfsBRJ::v_BuildPreconCfs()
    {
        NEKERROR(ErrorUtil::efatal, "v_BuildPreconCfs not defined");
    }

} // namespace Nektar
