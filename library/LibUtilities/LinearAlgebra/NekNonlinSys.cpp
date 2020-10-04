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
         * Solves a linear system using iterative methods.
         */
		NekNonlinSysFactory & GetNekNonlinSysFactory()
        {
            static NekNonlinSysFactory instance;
            return instance;
        }

        /// Constructor for full direct matrix solve.
        NekNonlinSys::NekNonlinSys(
            const LibUtilities::SessionReaderSharedPtr  &pSession,
            const LibUtilities::CommSharedPtr           &vComm,
            const int                                   nDimen)
            : NekSys(pSession, vComm, nDimen)
        {
            std::vector<std::string>  variables(1);
            variables[0] =  pSession->GetVariable(0);
            string variable = variables[0];

            if (pSession->DefinesGlobalSysSolnInfo(variable,
                                              "NonlinIteratSolverTolerance"))
            {
                m_tolerance = boost::lexical_cast<NekDouble>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "NonlinIteratSolverTolerance").c_str());
            }
            else
            {
                pSession->LoadParameter("NonlinIteratSolverTolerance",
                                        m_tolerance,
                                        NekConstants::kNekIterativeTol);
            }

            if (pSession->DefinesGlobalSysSolnInfo(variable,
                                                  "NonlinIteratMaxIterations"))
            {
                m_maxiter = boost::lexical_cast<int>(
                        pSession->GetGlobalSysSolnInfo(variable,
                                "NonlinIteratMaxIterations").c_str());
            }
            else
            {
                pSession->LoadParameter("NonlinIteratMaxIterations",
                                        m_maxiter,
                                        5000);
            }
        }

        void NekNonlinSys::v_InitObject()
        {
            NekSys::v_InitObject();
        }

        NekNonlinSys::~NekNonlinSys()
        {
        }
    }
}

