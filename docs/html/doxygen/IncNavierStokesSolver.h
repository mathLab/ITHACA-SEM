namespace Nektar {
/**
* \page pageIncNavierStokesSolver Incompressible Navier-Stokes Equations Solver
 
 Two very popular approaches to numerically solving the discrete incompressible Navier-Stokes equations
 are to either directly solve the matrix problem arising from the discretisation of the Stokes problem or 
 to use a splitting/projection method where the velocity matrix system and pressure systems are typically decoupled.
 
 \f[
 \frac{\partial V}{\partial t}+(V\cdot \nabla)V=-\nabla p + \nu \Delta V + f
 \f]
 \f[
 \nabla \cdot V=0
 \f]
 
 The direct solution of the Stokes system introduces the problem of appropriate space for the velocity and the pressure systems
 to satisfy the inf-sup condition [2] and requires the solution of the full velocity-pressure system.
 However, if a discontinuous pressure space is used then all but the constant mode of the pressure system can be decoupled
 from the velocity. Futhermore, when implementing this approach with a spectral/hp element dicretisation, the remaining velocity
 system may then be statically condensed to decouple the so called interior elemental degrees of freedom, reducing the Stokes problem
 to a smaller system expressed on the elemental boundaries [3,4]. The direct solution of the Stokes problem provides a very natural
 setting for the solution of the pressure system which is not easily dealt with in a splitting scheme. Futher, the solution
 of the full coupled velocity systems means that the introduction of a spatially varying viscosity, which arise for non-Newtonian flows,
 is only a minor modification.
 
 Splitting schemes are typically favoured for their numerical efficiency since for a Newtonian fluid the velocity and pressure
 are handled indipendently requiring the solution of three (in two dimensions) elliptic systems of rank N as opposed to a single system
 of rank 3N solved in the Stokes problem. However, a drawback of this approach is the splitting scheme error which is introduced
 when decoupling the pressure and the velocity system, although this can be made consistent wuth the overall temporal accuracy of
 the scheme by appropriate discretisation of the pressure boundary conditions [5].
 
 (The previous overview has been directly taken from [1]).
 
 Both of these approaches have been implemented in Nektar++, here futher details about the formulations:
 * - @subpage pageVelocityCorrectionScheme - a splitting scheme for the incompressible Navier-Stokes equations.
 * - @subpage pageCoupledSolver - a direct solution of the Navier-Stokes equations (coupled solver).
 
 The Incommpressible Navier-Stokes Equations solver contains other features:
 * - @subpage pageStability - stability analysis of the Navier-Stokes equations.

 [1] S.J. Sherwin and M. Ainsworth: <i>Unsteady Navier-Stokes Solvers Using Hybrid Spectral/hp Element Methods</i>, Conference Paper, 2000.<br />
 [2] R. Stenberg and M. Suri: <i>Mixed hp finite element methods for problems in elasticity and Stokes flows</i>, Numer. Math, 72, 367-389, 1996.<br />
 [3] M. Ainsworth and S.J. Sherwin: <i>Domain decomposition preconditioners for p and hp finite element approximation of the Stokes equations</i>, Comput. Methods Appl. Mech. Engrg., 175, 243-266, 1999. <br />
 [4] P. La Tallec and A. Patra: <i>Non-overlapping domain decomposition methods for adaptive hp approximations of the Stokes problem with discontiunous pressure field</i>, Comput. Methods Appl. Mech. Engrg., 145, 361-379, 1997.<br />
 [5] G. E. Karniadakis, M. Israeli, and S. A. Orszag: <i>High-order splitting methodsfor the incompressible Navier-Stokes equations</i>, J. Comput. Phys., 97, 414-443, 1991.<br />
 */	
}
