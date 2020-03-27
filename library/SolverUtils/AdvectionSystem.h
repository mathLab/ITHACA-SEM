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
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph);

    SOLVER_UTILS_EXPORT virtual ~AdvectionSystem();

    SOLVER_UTILS_EXPORT virtual void v_InitObject();

    /// Returns the advection object held by this instance.
    SOLVER_UTILS_EXPORT AdvectionSharedPtr GetAdvObject()
    {
        return m_advObject;
    }

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble>  GetElmtCFLVals(void);
    SOLVER_UTILS_EXPORT NekDouble               GetCFLEstimate(int &elmtid);

protected:
    /// Advection term
    SolverUtils::AdvectionSharedPtr m_advObject;

    SOLVER_UTILS_EXPORT virtual bool v_PostIntegrate(int step);

    SOLVER_UTILS_EXPORT virtual Array<OneD, NekDouble> v_GetMaxStdVelocity()
    {
        ASSERTL0(false,
            "v_GetMaxStdVelocity is not implemented by the base class.");
        Array<OneD, NekDouble> dummy(1);
        return dummy;
    }

private:
    /// dump cfl estimate
    int m_cflsteps;
    /// Write field if cfl is higher than IO_CFLWriteFld treshold
    NekDouble m_cflWriteFld;
    /// Number of timesteps after which IO_CFLWriteFld is activated
    int m_cflWriteFldWaitSteps;    

};

/// Shared pointer to an AdvectionSystem class
typedef std::shared_ptr<AdvectionSystem> AdvectionSystemSharedPtr;

}
}

#endif
