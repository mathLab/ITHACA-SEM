///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionSystem.h
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
// Description: Advection system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTIONSYSTEM_H
#define NEKTAR_SOLVERUTILS_ADVECTIONSYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/Advection/Advection.h>

namespace Nektar {
namespace SolverUtils {

/// A base class for PDEs which include an advection component
class AdvectionSystem: virtual public UnsteadySystem
{
public:
    SOLVER_UTILS_EXPORT AdvectionSystem(
            const LibUtilities::SessionReaderSharedPtr &pSession);

    SOLVER_UTILS_EXPORT virtual ~AdvectionSystem();

    SOLVER_UTILS_EXPORT virtual void v_InitObject();

    /// Returns the advection object held by this instance.
    AdvectionSharedPtr GetAdvObject()
    {
        return m_advObject;
    }

protected:
    /// Advection term
    SolverUtils::AdvectionSharedPtr m_advObject;
};

/// Shared pointer to an AdvectionSystem class
typedef boost::shared_ptr<AdvectionSystem> AdvectionSystemSharedPtr;

}
}

#endif
