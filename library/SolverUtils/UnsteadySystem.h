///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySystem.h
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
// Description: Generic timestepping for Unsteady solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H
#define NEKTAR_SOLVERUTILS_UNSTEADYSYSTEM_H

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSchemeOperators.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Filters/Filter.h>

namespace Nektar
{
namespace SolverUtils
{
/// Base class for unsteady solvers.
class UnsteadySystem : public EquationSystem
{
public:
    /// Destructor
    SOLVER_UTILS_EXPORT virtual ~UnsteadySystem();

    /// Calculate the larger time-step mantaining the problem stable.
    SOLVER_UTILS_EXPORT NekDouble
    GetTimeStep(const Array<OneD, const Array<OneD, NekDouble>> &inarray);

    /// CFL safety factor (comprise between 0 to 1).
    NekDouble m_cflSafetyFactor;
    /// CFL growth rate
    NekDouble m_CFLGrowth;
    /// maximun cfl in cfl growth
    NekDouble m_CFLEnd;

protected:
    /// Number of time steps between outputting status information.
    int m_infosteps;
    /// Number of steps between checks for abort conditions.
    int m_abortSteps;
    /// Number of time steps between outputting filters information.
    int m_filtersInfosteps;
    int m_nanSteps;
    /// Wrapper to the time integration scheme
    LibUtilities::TimeIntegrationSchemeSharedPtr    m_intScheme;
    /// The time integration scheme operators to use.
    LibUtilities::TimeIntegrationSchemeOperators    m_ode;
    ///
    LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr  m_intSoln;
    ///
    NekDouble m_epsilon;
    /// Indicates if explicit or implicit treatment of diffusion is used.
    bool m_explicitDiffusion;
    /// Indicates if explicit or implicit treatment of advection is used.
    bool m_explicitAdvection;
    /// Indicates if explicit or implicit treatment of reaction is used.
    bool m_explicitReaction;
    /// Flag to determine if simulation should start in homogeneous
    /// forward transformed state.
    bool m_homoInitialFwd;

    /// Tolerance to which steady state should be evaluated at
    NekDouble m_steadyStateTol;
    /// Check for steady state at step interval
    int m_steadyStateSteps;
    /// Storage for previous solution for steady-state check
    Array<OneD, Array<OneD, NekDouble>> m_previousSolution;
    // Steady-state residual file
    std::ofstream m_errFile;

    std::vector<int> m_intVariables;

    std::vector<std::pair<std::string, FilterSharedPtr>> m_filters;

    /// Number of time steps between outputting status information.
    NekDouble m_filterTimeWarning;

    /// Initialises UnsteadySystem class members.
    SOLVER_UTILS_EXPORT UnsteadySystem(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    /// Init object for UnsteadySystem class.
    SOLVER_UTILS_EXPORT virtual void v_InitObject();

    /// Get the maximum timestep estimator for cfl control.
    SOLVER_UTILS_EXPORT NekDouble MaxTimeStepEstimator();

    /// Solves an unsteady problem.
    SOLVER_UTILS_EXPORT virtual void v_DoSolve();

    /// Sets up initial conditions.
    SOLVER_UTILS_EXPORT virtual void v_DoInitialise();

    /// Print a summary of time stepping parameters.
    SOLVER_UTILS_EXPORT virtual void v_GenerateSummary(SummaryList &s);

    /// Print the solution at each solution point in a txt file
    SOLVER_UTILS_EXPORT virtual void v_AppendOutput1D(
        Array<OneD, Array<OneD, NekDouble>> &solution1D);

    SOLVER_UTILS_EXPORT virtual NekDouble v_GetTimeStep(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray);

    SOLVER_UTILS_EXPORT virtual bool v_PreIntegrate(int step);
    SOLVER_UTILS_EXPORT virtual bool v_PostIntegrate(int step);

    SOLVER_UTILS_EXPORT virtual bool v_RequireFwdTrans()
    {
        return true;
    }

    SOLVER_UTILS_EXPORT void CheckForRestartTime(NekDouble &time, int &nchk);

    /// \brief Evaluate the SVV diffusion coefficient
    /// according to Moura's paper where it should
    /// proportional to h time velocity
    SOLVER_UTILS_EXPORT void SVVVarDiffCoeff(
        const Array<OneD, Array<OneD, NekDouble>> vel,
        StdRegions::VarCoeffMap &varCoeffMap);

private:
    void InitializeSteadyState();

    bool CheckSteadyState(int step);
};

} // namespace SolverUtils
} // namespace Nektar

#endif
