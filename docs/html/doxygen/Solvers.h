namespace Nektar {
/**
 * \page pageSolvers Solvers
 *
 * A number of solvers are provided with the Nektar++ to demonstrate how to
 * use the toolkit to solve a number of common (as well as less-common)
 * problems. These include:
 * - @subpage pageADRSolver ADRSolver - a general solver for solving a number of linear PDE systems
 *   which are in this form using either a Continuous Galerkin or Discontinuous
 *   Galerkin approach.
 * - ADR2DManifold - solves mono-domain models such as the
 *   FitzHugh-Nagumo equations or the Aliev-Panfilov model on a 2D manifold
 *   embedded in 3-space.
 * - IncNavierStokes - An incompressible Navier-Stokes solver.
 *
 * All these solvers utilise a general timestepping scheme:
 * - @subpage pageGeneralLinearMethods.
 */

/**
 * \page pageADRSolver ADRSolver
 *
 * This solves a range of problems using the spectral/hp discretisation and is
 * a good starting point for using the library. A number of example session
 * files are included in the \c solvers/ADRSolver/Examples directory.
 *
 * To run the solver provide an XML session file defining the problem to solve:
 * @verbatim ./ADRSolver <session>.xml @endverbatim
 *
 * A range of solvers are already implemented in the framework:
 * <table border>
 * <tr><td>Helmholtz</td><td></td></tr>
 * <tr><td>Laplace</td><td></td></tr>
 * <tr>
 * <td>MonodomainAlievPanfilov</td>
 * <td>Phenomological model of canine cardiac electrophysiology.</td>
 * </tr>
 * <tr>
 * <td>MonodomainFitzHughNagumo</td>
 * <td>Phenomological model of nerve cell electrophysiology.</td>
 * </tr>
 * <tr><td>Poisson</td><td></td></tr>
 * <tr><td>SteadyAdvectionDiffusion</td><td></td></tr>
 * <tr><td>SteadyAdvectionDiffusionReaction</td><td></td></tr>
 * <tr><td>SteadyDiffusion</td>
 *     <td>Alias for Poisson</td></tr>
 * <tr><td>SteadyDiffusionReaction</td>
 *     <td>Alias for Helmholtz</td></tr>
 * <tr><td>UnsteadyAdvection</td>
 *     <td>Unsteady Advection equation.</td></tr>
 * <tr><td>UnsteadyAdvectionDiffusion</td><td></td></tr>
 * <tr><td>UnsteadyDiffusion</td><td></td></tr>
 * <tr><td>UnsteadyInviscidBurger</td><td></td></tr>
 * </table>
 *
 * \section sectionADRSolverAdvection Example: Advection Equation
 *
 * As an example, consider the advection equation
 * \f[ \frac{\partial u}{\partial t} = \frac{\partial u}{\partial x} \f]
 * on a 2D domain with mixed time-dependent boundary conditions.
 * A corresponding session file is given in
 * \c solvers/ADRSolver/Examples/Advection.xml.
 * The file consists of three sections:
 * - GEOMETRY: defining the mesh
 * - EXPANSIONS: defining the polynomial basis to use on each mesh element
 * - CONDITIONS: defining the parameters of the problem to solve
 * General information on the XML file format used by Nektar++ can be found
 * @ref pageXML here.
 *
 * A number of components of the \<CONDITIONS> section are of importance here.
 * First, we set the \<PARAMETERS> for our time integration by setting the
 * FinTime, NumSteps and TimeStep values.
 * @code
 *    <P> FinTime        = 0.75             </P>
 *    <P> NumSteps       = 500              </P>
 *    <P> TimeStep       = FinTime/NumSteps </P>
 * @endcode
 * The frequency of checkpoint files and status updates can be controlled
 * through a further two parameters:
 * @code
 *    <P> IO_CheckSteps  = NumSteps         </P>
 *    <P> IO_InfoSteps   = 1000             </P>
 * @endcode
 * The final two parameters indicate the advection velocity. Although these are
 * not directly used by the solver, they are embedded in the boundary
 * conditions defined later.
 * @code
 *    <P> advx           = 1                </P>
 *    <P> advy           = 1                </P>
 * @endcode
 *
 * The \<SOLVERINFO> section allows us to specify the solver, the type
 * of projection (continuous or discontinuous) and the time integration scheme
 * to use.
 * @code
 * <SOLVERINFO>
 *    <I PROPERTY="EQTYPE" VALUE="UnsteadyAdvection"/>
 *    <I PROPERTY="Projection" VALUE="Continuous"/>
 *    <I PROPERTY="TimeIntegrationMethod"  VALUE="ClassicalRungeKutta4"/>
 * </SOLVERINFO>
 * @endcode
 *
 * For the advection equation we must specify the advection velocity for the
 * solver. This is done in the \<USERDEFINEDEQNS> subsection:
 * @code
 * <USERDEFINEDEQNS>
 *    <F LHS="Vx" VALUE="1"  />
 *    <F LHS="Vy" VALUE="1" />
 * </USERDEFINEDEQNS>
 * @endcode
 *
 * \section sectionADRSolverNewModule Creating new modules
 * Details on how to implement a new module can be found
 * \subpage pageADRSolverImplement here.
 */

/**
 * \page pageADRSolverImplement ADRSolver Framework
 *
 * The ADRSolver provides a general-purpose framework for implementing solvers
 * for both \ref subsectionADRSolverModuleImplementationSteady steady-state
 * and \ref subsectionADRSolverModuleImplementationUnsteady time-dependent
 * problems. It is designed in a modular fashion to allow the implementation
 * of additional solvers with maximum ease and minimal complexity. Furthermore,
 * each module automatically registers itself with the LibUtilities::NekFactory
 * class, allowing the ADRSolver to automatically load the correct module based
 * on the information in the supplied session file without any need to change
 * the underlying framework code. A complete list of available modules can be
 * printed through a call to the static
 * LibUtilities::NekFactory::PrintAvailableClasses:
 * \code
 * EquationSystemFactory::PrintAvailableClasses();
 * \endcode
 *
 * Some of the solvers currently included in the framework are:
 * - Laplace, Poisson, Helmholtz - solvers for these steady-state equations.
 * - UnsteadyAdvection - a solver for time-dependent advection problems.
 * - UnsteadyDiffusion - a solver for time-dependent diffusion problems.
 * - AlievPanfilov - a solver for cardiac electrophysiology problems.
 *
 * \section sectionADRSolverModules Modules
 * Modules are added to an inheritence tree rooted at the EquationSystem class.
 * This class provides the underlying functionality and means by which to
 * generalise the solver to accommodate the broad range of problems.
 * Each distinct problem is implemented in a unique module, derived in one way
 * or another from the EquationSystem class. Consequently, only the aspects
 * of the problem which makes it unique need to be implemented in the module.
 *
 *
 * \section sectionADRSolverModuleImplementation Implementing a new Module
 * To implement a new module requires implementation of a number of functions
 * which will be called by the solver framework or time integration scheme.
 * As a first step in implementing a new solver, it is necessary to add the
 * common elements needed to interface with the framework.
 * -# Create a new class header in the ADRSolver/EquationSystems/MySolver.h:
 *    \code
 *    // MySolver.h
 *    #ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MYSOLVER_H
 *    #define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MYSOLVER_H
 *
 *    #include <ADRSolver/EquationSystem.h>
 *
 *    namespace Nektar
 *    {
 *        class MySolver : public EquationSystem
 *        {
 *        public:
 *            MySolver(SessionReaderSharedPtr& pSession);
 *            virtual ~MySolver;
 *        }
 *    }
 *
 *    #endif
 *    \endcode
 *    The base class depends on the type of solver:
 *    - For steady-state solvers, use EquationSystem (as above).
 *    - For time-dependent problems, use UnsteadySystem.
 *      \code
 *      #include <ADRSolver/EquationSystems/UnsteadySystem.h>
 *      \endcode
 * -# Add the ability for the class to register with LibUtilities::NekFactory.
 *    This is essential for the module to fit into the ADRSolver framework.
 *    \code
 *    // MySolver.h
 *    #ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MYSOLVER_H
 *    #define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MYSOLVER_H
 *
 *    #include <ADRSolver/EquationSystem.h>
 *
 *    namespace Nektar
 *    {
 *        class MySolver : public EquationSystem
 *        {
 *        public:
 *            static EquationSystemSharedPtr create(SessionReaderSharedPtr& pSession)
 *            {
 *                return MemoryManager<MySolver>::AllocateSharedPtr(pSession);
 *            }
 *
 *            static std::string className;
 *
 *            MySolver(SessionReaderSharedPtr& pSession);
 *            virtual ~MySolver;
 *        }
 *    }
 *
 *    #endif
 *    \endcode
 * -# Create an implementation file ADRSolver/EquationSystems/MySolver.cpp and
 *    implement the registration of the new class with the NekFactory, along
 *    with implementation of the constructor and destructor.
 *    \code
 *    // MySolver.cpp
 *    #include <ADRSolver/EquationSystems/MySolver.h>
 *
 *    namespace Nektar
 *    {
 *        string MySolver::className = EquationSystemFactory::RegisterCreatorFunction("MySolver", MySolver::create, "Description of solver");
 *
 *        MySolver::MySolver(SessionReaderSharedPtr& pSession)
 *                : EquationSystem(pSession)
 *        {
 *
 *        }
 *
 *        MySolver::~MySolver()
 *        {
 *
 *        }
 *    }
 *    \endcode
 *
 * \subsection subsectionADRSolverModuleImplementationSteady Steady-state Solvers
 * To create a steady-state solver, we provide the problem-specific
 * implementation. This is achieved by adding the routines which will actually
 * perform the solve,
 * \code
 * protected:
 *     virtual void v_DoSolve();
 *     virtual void v_PrintSummary(std::ostream &out);
 * \endcode
 *
 * The first of these two routines is where the code to solve the problem
 * should be placed. For example, in the Laplace problem, we simply call the
 * HelmSolve routine on each of the solution variables:
 * \code
 * void Laplace::v_DoSolve()
 * {
 *     for(int i = 0; i < m_fields.num_elements(); ++i)
 *     {
 *         m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
 *                                m_fields[i]->UpdateCoeffs(),
 *                                m_lambda);
 *         m_fields[i]->SetPhysState(false);
 *     }
 * }
 * \endcode
 *
 *
 * \subsection subsectionADRSolverModuleImplementationUnsteady Time-dependent Solvers
 * In the case of a time-dependent problem, the new module should inherit the
 * class UnsteadySystem (which itself inherits EquationSystem). This provides
 * the time-stepping algorithms for unsteady problems.
 *
 * Depending on the nature of the problem, and the time integration scheme
 * chosen, one or more of the following routines will need to be implemented.
 * \code
 * protected:
 *     void DoOdeRhs(
 *              const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
 *                    Array<OneD,        Array<OneD, NekDouble> >&outarray,
 *              const NekDouble time);
 *
 *     void DoImplicitSolve(
 *              const Array<OneD, const Array<OneD,      NekDouble> >&inarray,
 *                    Array<OneD, Array<OneD, NekDouble> >&outarray,
 *                    NekDouble time,
 *                    NekDouble lambda);
 *
 *     void DoOdeProjection(
 *              const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
 *                    Array<OneD,  Array<OneD, NekDouble> > &outarray,
 *              const NekDouble time);
 * \endcode
 *
 * Finally, to configure the time integration scheme to use these routines,
 * we must identify them with the TimeIntegrationSchemeOperators object defined
 * in UnsteadySystem. In the constructor for our new module, we add suitable
 * statements, for example
 * \code
 * m_ode.DefineOdeRhs        (&MySolver::DoOdeRhs,        this);
 * m_ode.DefineImplicitSolve (&MySolver::DoImplicitSolve, this);
 * m_ode.DefineProjection    (&MySolver::DoOdeProjection, this);
 * \endcode
 */
}
