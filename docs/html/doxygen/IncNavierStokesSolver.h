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
 to satisfy the inf-sup condition [] and requires the solution of the full velocity-pressure system.
 However, if a discontinuous pressure space is used then all but the constant mode of the pressure system can be decoupled
 from the velocity. Futhermore, when implementing this approach with a spectral/hp element dicretisation, the remaining velocity
 system may then be statically condensed to decouple the so called interior elemental degrees of freedom, reducing the Stokes problem
 to a smaller system expressed on the elemental boundaries []. The direct solution of the Stokes problem provides a very natural
 setting for the solution of the pressure system which is not easily dealt with in a splitting scheme. Futher, the solution
 of the full coupled velocity systems means that the introduction of a spatially varying viscosity, which arise for non-Newtonian flows,
 is only a minor modification.
 
 Splitting schemes are typically favoured for their numerical efficiency since for a Newtonian fluid the velocity and pressure
 are handled indipendently requiring the solution of three (in two dimensions) elliptic systems of rank N as opposed to a single system
 of rank 3N solved in the Stokes problem. However, a drawback of this approach is the splitting scheme error which is introduced
 when decoupling the pressure and the velocity system, although this can be made consistent wuth the overall temporal accuracy of
 the scheme by appropriate discretisation of the pressure boundary conditions [].
 
 (The previous overview has been directly taken from []).
 
 Both of these approaches have been implemented in Nektar++, here futher detail about the formulations:
 * - @subpage pageVelocityCorrectionScheme - a splitting scheme for the incompressible Navier-Stokes equations.
 * - @subpage pageCoupledSolver - a direct solution of the Navier-Stokes equations (coupled solver).
 
 The Incommpressible Navier-Stokes Equations solver contain other features:
 * - @subpage pageAdvection - a collection of methods to deal with the advection term calculation inside the equations.
 * - @subpage pageStability - stability analysis of the Navier-Stokes equations.
*/	
	
}
