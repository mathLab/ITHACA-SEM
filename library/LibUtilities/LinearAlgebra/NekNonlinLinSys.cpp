///////////////////////////////////////////////////////////////////////////////
//
// File:  NekNonlinLinSys.cpp
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
// Description:  NekNonlinLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinLinSys.h>

using namespace std;

namespace Nektar
{       
    namespace LibUtilities
    {
    /**
     * @class  NekNonlinLinSys
     *
     * Solves a linear system using iterative methods.
     */

    /// Constructor for full direct matrix solve.
    NekNonlinLinSys::NekNonlinLinSys(
        const LibUtilities::SessionReaderSharedPtr  &pSession,
        const LibUtilities::CommSharedPtr           &vComm,
        const int                                   nDimen)
    {
        m_session = pSession;
        m_tolerance         =   1.0E-15;
        m_verbose           =   false;
        m_root              =   false;
        m_Comm              =   vComm;

        if (0==m_Comm->GetRank())
        {
            m_root              =   true;
        }
        m_verbose   =   pSession->DefinesCmdLineArgument("verbose");

        m_converged = false;

        m_SysDimen = nDimen;
    }

    NekNonlinLinSys::~NekNonlinLinSys()
    {
    }

    }

}

