namespace Nektar {
/**
 * \page pageShallowWaterSolver Shallow Water Equation Solver
 *
 * The <code>ShallowWaterSolver</code> is a solver for depth-integrated wave 
 * equations of shallow water type. Presently the following equations 
 * are supported: 
 * <table border>
 * <tr><td>LinearSWE</td><td>Linearized SWE solver in primitive variables
 * (constant still water depth)</td></tr>
 * <tr><td>NonlinearSWE</td><td>Nonlinear SWE solver in conservative variables
 * (constant still water depth)</td></tr>
 * </table>
 *
 * \section The Shallow Water Equations
 * The shallow water equations (SWE) is a two-dimensional system of nonlinear partial
 * differential equations of hyperbolic type that are fundamental in hydraulic, coastal
 * and environmental engineering. In deriving the SWE the vertical velocity is 
 * considered negligible and the horizontal velocities are assumed uniform with depth. 
 * The SWE are hence valid when the water depth can be considered small compared to the 
 * characteristic length scale of the problem, as typical for flows in rivers and shallow 
 * coastal areas. Despite the limiting restrictions the SWE can be used to describe many 
 * important phenomena, for example storm surges, tsunamis and river flooding. 
 *
 * The two-dimensional SWE is stated in conservation form as
 * 
 * \f[
 * \frac{\partial {\bf U}}{\partial t} + \nabla \cdot {\bf F(U)} = {\bf
 * S(U)}\,, \label{eq:SWE}
 * \f]
 *
 * where \f${\bf F(U)} = \left[{\bf E(U)}\,, {\bf G(U)} \right]\f$ is the
 * flux vector and the vector of conserved variables read
 * \f${\bf U}=\left[H\,,Hu\,,Hv \right]^\mathrm{T}\,.\f$ Here 
 * \f$H({\bf x},t)=\zeta({\bf x},t) + d({\bf x})\f$ is the total
 * water depth, \f$\zeta({\bf x},t)\f$ is the free surface elevation and
 * \f$d({\bf x})\f$ is the still water depth. The depth-averaged velocity
 * is denoted by \f${\bf u} = \left[u({\bf x},t)\,, v({\bf x},t)\right]^\mathrm{T}\f$,
 * where \f$u\f$ and \f$v\f$ are the velocities in the \f$x\f$- and \f$y\f$-directions,
 * respectively. The content of the flux vector is
 * \f[
 * {\bf E(U)} = \left[ \begin{array}{c} Hu\\Hu^2 +
 * gH^2/2\\Huv\end{array}\right]\,, \qquad {\bf G(U)} = \left[
 * \begin{array}{c} Hv\\Hvu\\Hv^2 + gH^2/2\end{array}\right]\,,
 * \f]
 *
 * in which \f$g\f$ is the acceleration due to gravity. The source term
 * \f${\bf S(U)}\f$ accounts for, e.g., forcing due to bed
 * friction, bed slope, Coriolis force and higher-order dispersive effects 
 * (Boussinesq terms). In the distributed version of the <code>ShallowWaterSolver</code>
 * only the Coriolis force is included. 
 *
 * \section sectionSWEExample An example: the Rossby modon case
 * In the following we will look closer on the file
 * @code
 * Nektar++/solvers/ShallowWaterSolver/Examples/Rossby_Nonlinear_DG.xml
 * @endcode
 * which gives the input data for performing a discontinuous Galerkin 
 * simulation of the westward propagation of an equatorial Rossby modon. 
 *
 * \subsection subsectionSWEInput Input Options
 * A detailed description of the input file for Nektar++ can be found in 
 * @ref pageXML. For what concern the <code>ShallowWaterSolver</code> the
 * <SOLVERINFO> section allows us to specify the solver, the type of projection 
 * (continuous or discontinuous), the explicit time integration scheme to
 * use and (in the case the discontinuous Galerkin method is used) 
 * the choice of numerical flux. A typical example would be:
 * @code
 * <SOLVERINFO>      
 *   <I PROPERTY="EqType" VALUE="NonlinearSWE">
 *   <I PROPERTY="Projection" VALUE="DisContinuous">
 *   <I PROPERTY="TimeIntegrationMethod" VALUE="ClassicalRungeKutta4">
 *   <I PROPERTY="UpwindType" VALUE="HLLC">
 * </SOLVERINFO>
 * @endcode
 * 
 * In the <PARAMETERS> section we, in addition to the normal setting
 * of time step etc., also define the acceleration of gravity by 
 * setting the parameter "Gravity": 
 * @code
 * <PARAMETERS>
 *    <P> TimeStep       = 0.04             </P>
 *    <P> NumSteps       = 1000             </P>
 *    <P> IO_CheckSteps  = 100              </P>
 *    <P> IO_InfoSteps   = 100              </P>
 *    <P> Gravity        = 1.0              </P>
 * </PARAMETERS>
 * @endcode
 *  
 * In the <USERDEFINEDEQNS> we specify "f" which is the Coriolis
 * parameter and "d" denoting the still water depth:
 * @code
 * <USERDEFINEDEQNS>
 *    <F LHS="f" VALUE="0+1*y" />
 *    <F LHS="d" VALUE="1"     />
 * </USERDEFINEDEQNS>   
 * @endcode 
 *
 * Initial values and boundary conditions are given in terms of primitive variables
 * (please note that also the output files are given in terms of primitive
 * variables). For the discontinuous Galerkin we typically enforce any slip wall boundaries
 * weakly using symmetry technique. This is given by the USERDEFINEDTYPE="Wall" choice in the
 * <BOUNDARYCONDITIONS> section:
 * @code 
 * <BOUNDARYCONDITIONS>
 *   <REGION REF="0">
 *    <D VAR="eta"  USERDEFINEDTYPE="Wall"  VALUE="0" />
 *	<D VAR="u"    USERDEFINEDTYPE="Wall"  VALUE="0" />
 *	<D VAR="v"    USERDEFINEDTYPE="Wall"  VALUE="0" />
 *   </REGION>
 * </BOUNDARYCONDITIONS>
 * @endcode
 *
 * \subsection subsectionExecutiongSWE Running the code
 * After the input file has been copied to the build directory
 * of the <code>ShallowWaterSolver</code> the code can be executed by:
 * @code
 * ./ShallowWaterSolver Rossby_Nonlinear_DG.xml
 * @endcode
 *
 * \subsection subsectionSWEPostProcess Post-proceesing
 * After the final time step the solver will write an output file 
 * <code>RossbyModon_Nonlinear_DG.fld</code>. We can convert it to tecplot
 * format by using the <code>FldToTecplot</code> utility. Thus we execute the following command:
 * @code
 * ../../../utilities/builds/PostProcessing/FldToTecplot RossbyModon_Nonlinear_DG.xml RossbyModon_Nonlinear_DG.fld 
 * @endcode
 * This will generate a file called <code>RossbyModon_Nonlinear_DG.dat</code> that
 * can be loaded directly into tecplot:
 *
 *  \image html RossbyModon.png "Rossby modon after 40 time units." 
 */
}
