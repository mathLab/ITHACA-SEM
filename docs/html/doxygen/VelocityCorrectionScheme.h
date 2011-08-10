namespace Nektar {
/**
* \page pageVelocityCorrectionScheme Velocity Correction Scheme
 
 Here, we will explain how to use the incompressible Navier-Stokes solver of Nektar++ using a splitting method (Velocity Correction Scheme). 
 Its source code is located at <code>Nektar++/solvers/IncNavierStokesSolver/EquationSystems</code>.
 We will start by showing the formulation and demonstrating how to run a first example with the incompressible Navier-Stokes 
 solver and briefly explain the options of the solver specified in the input file. 
 These explanations will be based on the example file <code>QuadMesh.xml</code> 
 that can be found in <code>Nektar++/solvers/IncNavierStokesSolver/Examples</code>. 
 The  <code>QuadMesh.xml</code> example file describes a channel flow in x-direction.
 After these introductory explanations, we will present some numerical results 
 in order to demonstrate the capabilities and the accuracy of the solver.
 - @ref sectionFormulation
 - @ref sectionRunning1example
 - @ref sectionInputOptions
 - @ref sectionResult
        - @ref subsectionKovasznay
        - @ref subsectionTaylor
        - @ref subsectionShedding
 - @ref sectionRef
 
 
 
 \section sectionFormulation Formulation
 
 Briefly, the Navier-Stokes solver computes the solution of the Navier-Stokes equation 
 in two space dimensions in terms of the primitive variables velocity
 \f$V=(u,v)\f$ and pressure \f$p\f$ for a fluid with kinematic viscosity \f$\nu\f$
 subject to appropriate initial and boundary conditions:
 \f[
 \frac{\partial V}{\partial t}+(V\cdot \nabla)V=-\nabla p + \nu \Delta V + f
 \f]
 \f[
 \nabla \cdot V=0
 \f]
 
 For example, a typical set of boundary conditions would be given by:
 \f[
 V=V_{in} \qquad at\ the\ inflow\ boundary,
 \f]
 \f[
 V=0 \qquad at\ walls,
 \f]
 \f[
 \frac{\partial V}{\partial n}=0 \qquad at\ the\ outflow\ boundary,
 \f]
 
 where \f$n\f$ is the unit outward normal vector on the boundary. 
 For simplicity, we assume that the density of the fluid is \f$\rho=1\f$. 
 The Navier-Stokes equations describe an incompressible, isothermal flow with constant properties.
 
 The incompressible Navier-Stokes equations are solved in time using a <strong>velocity correction splitting scheme</strong> combined 
 with <strong>stiffly-stable time integration</strong> which is derived from the work of Karniadakis, Israeli and Orszag [3]. 
 Briefly, the time integration scheme consists of the following steps: 
 
 - Calculate a first intermediate velocity field by calculating the advection term explicitly and 
 combining it with the solution at previous time-steps
 \f[
 \frac{\tilde{V}-\sum_{q=0}^{J-1}\frac{\alpha_q}{\gamma_0}V^{n-q}}{\Delta t}=\sum_{q=0}^{J-1}\frac{\beta_q}{\gamma_0}[-(V\cdot\nabla)V]^{n-q} + f^{n+1}
 \f]

 - Solve a Poisson equation to obtain the pressure solution at the new time level
 \f[
 \Delta p^{n+1}=\Bigl(\frac{\gamma_0}{\Delta t}\Bigr)\nabla\cdot\tilde{V}
 \f] 
 with consistent boundary conditions
 \f[
 \frac{\partial p}{\partial n}^{n+1}= - \Bigl[ \frac{\partial V}{\partial t}^{n+1} + \nu \sum_{q=0}^{J-1}\beta_q(\nabla \times \nabla \times V)^{n-q} + \sum_{q=0}^{J-1}\beta_q [(V \cdot \nabla)V]^{n-q}\Bigr]\cdot n
 \f]

 - Calculate a second intermediate velocity field 
 \f[
 \tilde{\tilde{V}}=\tilde{V}-\Bigl(\frac{\Delta t}{\gamma_0}\Bigr)\nabla p^{n+1}
 \f]

 - Use the second intermediate velocity field as a forcing term in a Helmholtz problem to obtain the velocity field at the new time level
 \f[
 \Bigl(\Delta-\frac{\gamma_0}{\Delta t \nu}\Bigr)V^{n+1}=-\Bigl(\frac{\gamma_0}{\Delta t \nu}\Bigr)\tilde{\tilde{V}}
 \f]
 
 Here, \f$J\f$ is the integration order and \f$\gamma_0\f$, \f$\alpha_q\f$, \f$\beta_q\f$ 
 are the stiffly stable time integration coefficients given in the table below.
 
\f{displaymath}
 \begin{array}{cccc}
          & 1st\ order & 2nd\ order & 3rd\ order \\
 \gamma_0 & 1 & 3/2 & 11/6\\
 \alpha_0 & 1 & 2 & 3\\
 \alpha_1 & 0 & -1/2 & -3/2\\
 \alpha_2 & 0 & 0 & 1/3\\
 \beta_0 & 1 & 2 & 3\\
 \beta_1 & 0 & -1 & -3\\
 \beta_2 & 0 & 0 & 1\\
 \end{array}
\f}
 
This multi-step implicit-explicit splitting scheme decouples the velocity field \f$V\f$ from the pressure \f$p\f$, leading to an
explicit treatment of the advection term and an implicit treatment of the pressure and the diffusion term.
For details on the velocity correction scheme and the stiffly stable time integrators see Chapter 8.3.2.5 and 8.3.2.6 in the book of Karniadakis and Sherwin [1]. 
We can choose between the first-order, the second-order or the third-order scheme by setting <code class="source">TimeIntegrationMethod</code> to <code class="source">
IMEXOrder1,IMEXOrder2</code> or <code class="source">IMEXOrder3</code>. For more details about the input file see @ref pageXML.
 
 \section sectionRunning1example Running a first example
 
 The folder <code>Nektar++/solvers/IncNavierStokesSolver/Examples</code> contains several <code>*.xml</code> files.
 (Further examples can be found in <code>Nektar++/regressionTests/Solvers/IncNavierStokesSolver/InputFiles</code>)
 These are input files for the Navier-Stokes solver specifying the geometry (i.e. the mesh and 
 the spectal/hp expansion), the parameters and boundary conditions. Further details on the structure 
 of the input file can be found in @ref pageXML. 
  
 Now, lets try to run the <code>QuadMesh.xml</code> example that can be found in 
 <code>Nektar++/solvers/IncNavierStokesSolver</code>:

 - Copy the input file, <code>QuadMesh.xml</code> to the directory where the solver is compiled, e.g. 
   <code>Nektar++/solvers/IncNavierStokesSolver/build/IncNavierStokesSolver</code>.
 - Run the code by typing
 @code
 ./IncNavierStokesSolver QuadMesh.xml
 @endcode

 The solution should now have been written to the file <code>QuadMesh.fld</code>. 
 This file is formatted in the Nektar++ output format.
 To visualise the solution, we can convert the fld-file into Vtk, Gmsh or Tecplot file formats using 
 the Post-processing tools in <code>Nektar++/utilities/builds/PostProcessing/</code>. 
 Here, we will demonstrate how to convert the <code>fld</code>-file into Vtk-file format. 
 The <a href="http://www.vtk.org/">Visualization Toolkit (VTK)</a> is an open-source visualisation library for 3D computer graphics. 
 There is a wide range of open source GUIs available for VTK, e.g. <a href="http://www.paraview.org/">paraview</a>.
 We convert the <code>fld</code>-file into Vtk-file format by typing
 @code
 ../../../utilities/builds/PostProcessing/FldToVtk QuadMesh.xml QuadMesh.fld
 @endcode
 
 (Make sure you are still at the directory <code>Nektar++/solvers/IncNavierStokesSolver/build/IncNavierStokesSolver</code>.)
 This should have been created the file <code>QuadMesh.vtu</code> in the current directory. 
 The <code>*.vtu</code> files can now be opened with paraview, for example. 
 
 \image html QuadMeshu.png "Channel Flow (u-velocity component)" 
 
\section sectionInputOptions Input Options
 
 A detailed descriptions of the input file for Nektar++ can be found in @ref pageXML. 
 For what concern the incompressible Navier-Stokes solver a typical set is:
 @code
 <SOLVERINFO>      
   <I PROPERTY="EQTYPE" VALUE="UnsteadyNavierStokes"/>
   <I PROPERTY="SolverType"  VALUE="VelocityCorrectionScheme"/>
   <I PROPERTY="AdvectionForm" VALUE="Convective"/>
   <I PROPERTY="Projection" VALUE="Galerkin"/>
   <I PROPERTY="TimeIntegrationMethod" VALUE="IMEXOrder2"/>
 </SOLVERINFO>
 @endcode
 
 @code
 <PARAMETERS>
   <P> TimeStep      = 0.001 <P>
   <P> NumSteps      = 100   <P>
   <P> IO_CheckSteps = 10    <P>
   <P> IO_InfoSteps  = 10    <P>
   <P> Kinvis        = 1.0   <P>
 </PARAMETERS>
 @endcode
 
 The only material property we have to specify for the incompressible Navier-Stokes equation is the kinematic viscosity, 
 which is set by assigning a floating point number to <code class="source">Kinvis</code>. 
 The parameter <code class="source">TimeStep</code> defines \f$\Delta t\f$ in the time integration scheme and <code class="source">NumSteps</code> 
 is the total number of steps the time integration scheme will perform. 
 That means the physical endtime will be <code>TimeStep*NumSteps</code>. <code class="source">IO_CheckSteps</code> specifies the number of time steps 
 that will be performed between writing the solution to output files. 
 These output files are named <code>&lt;INPUTFILENAME&gt;_&lt;NUMBER&gt;.chk</code>. 
 They can be converted into Vtk, Tecplot or Gmsh file formats using the same Post-Processing tools as for the <code>fld</code>-files. 
 The conversion tools can be found in <code>Nektar++/utilities/builds/PostProcessing/</code>.
 <code class="source">IO_InfoSteps</code> specifies the number of time steps between displaying information on the console.
 
 @code
 <VARIABLES>
   <V ID="0"> u </V> 
   <V ID="1"> v </V> 
   <V ID="2"> p </V> 
 </VARIABLES>
 @endcode
 
 @code
 <BOUNDARYREGIONS>
   <B ID="0"> C[1] </B>
   <B ID="1"> C[2] </B>
   <B ID="2"> C[3] </B>
 </BOUNDARYREGIONS>
 @endcode
 
 @code
 <BOUNDARYCONDITIONS>
  <REGION REF="0">
    <D VAR="u" VALUE="0" />
    <D VAR="v" VALUE="0" />
	<N VAR="p" USERDEFINEDTYPE="H" VALUE="0"/>  
  </REGION>
  <REGION REF="1">
    <D VAR="u" VALUE="y*(1-y)" />
    <D VAR="v" VALUE="0" />
    <N VAR="p" USERDEFINEDTYPE="H" VALUE="0"/>  
  </REGION>
  <REGION REF="2">
    <N VAR="u" VALUE="0" />
    <N VAR="v" VALUE="0" />
    <D VAR="p" VALUE="0" />
  </REGION>
 </BOUNDARYCONDITIONS>
 @endcode
 
 In this example, Region 0 shows typical boundary conditions for walls, Region 1 for inlet boundaries and Region 2 for outlet boundaries.
 The options <code class="source">N VAR="p" USERDEFINEDTYPE="H" VALUE="0"</code> specifies high-order boundary conditions for the pressure which
 are calculated as shown in @ref sectionFormulation.
 Additionally, <strong>time-dependent</strong> boundary conditions can be specified in the input file by setting 
 <code class="source">USERDEFINEDTYPE="TimeDependent"</code> in the boundary region (see @ref pageXML). 
 For example, see the  <code>Test_TaylorVor_dt1.xml</code> file in the folder  
 <code>Nektar++/regressionTests/Solvers/IncNavierStokesSolver/InputFiles/</code> describing a Taylor Decaying Vortex.
 
 In case of a forcing term, the term can be specified as all the others Analitic functions
 @code
 <FUNCTION NAME="BodyForce">
   <E VAR="u" VALUE="0" />
   <E VAR="v" VALUE="0" />
   <E VAR="p" VALUE="0" />
 </FUNCTION>
 @endcode
 
 For the set of the initial conditions, exact solutions, special bounary conditons and the Quasi-3D approach see @ref pageXML.
 
 \section sectionResult Numerical Results
 
 \subsection subsectionKovasznay Kovasznay Flow 2D
 
 To investigate the convergence rate of the Navier-Stokes solver, 
 we compare the computed solution with the analytical solution for a Kovasznay flow. 
 The input file <code>Test_KovFlow.xml</code> for this flow configuration can be found at 
 <code>Nektar++/solvers/IncNavierStokesSolver/Examples</code>.
 In 1948, Kovasznay solved the problem of steady, laminar flow behind a two-dimensional grid. 
 The exact solution to the Navier-Stokes equations is given by:
 \f[
 u=1-e^{\lambda x}\cos{2 \pi y}
 \f]
 
 \f[
 v=\frac{\lambda}{2\pi}e^{\lambda x}\sin{2\pi y}
 \f]
 
 \f[
 p=\frac{1}{2}(1-e^{2 \lambda x})
 \f]
 
 \f[
 \lambda=\frac{1}{2\nu}-\Bigl [ \frac{1}{4\nu^2}+4\pi^2\Bigr]^{\frac{1}{2}}
 \f]
 
 The figure below shows the streamlines of the 2D solution obtained with the Navier-Stokes solver.
 A 12 quadrilaterals mesh has been used in a rectangular domain defined as \f$x\in[-0.5,1]\times[-0.5,1.5]\f$.
 The solution looks similar to the low-speed flow of a viscous fluid past an array of cylinders. 
 The Reynolds number has been set to \f$Re=1/\nu=40\f$. The analytical solution for the velocity and the pressure have 
 been used to set Dirichlet boundary conditions. 
 Knowing the exact solution, we can calculate the convergence rate with respect to the expansion order. 
 We observe exponential convergence of the error up to error values of \f$10^{-12}\f$.
 
 \image html kov.png "Streamlines of the 2D solution for a Kovasznay flow computed with the Navier-Stokes solver"
 
 \image html kov_dt0001.png "Convergence rate in the L2-norm"
 
 
 \subsection subsectionTaylor Taylor Decaying Vortex 2D
 
 To evaluate the time accuracy of the Navier-Stokes solver, we compute a Taylor decaying vortex described by the equations
 
 \f[
 u=-\cos(x)\sin(y)e^{-2t/Re}
 \f]
 
 \f[
 v=\sin(x)\cos(y)e^{-2t/Re}
 \f]
 
 \f[
 p=-\frac{1}{4}\Bigl(\cos(2x)+\cos(2y)\Bigr)e^{-4t/Re}
 \f]
 
 in a square domain \f$x\in[-\pi/2,\pi/2]\times[-\pi/2,\pi/2]\f$. 
 The input file <code>Test_TaylorVor_dt1.xml</code> for this flow configuration can be found at 
 <code>Nektar++/regressionTests/Solvers/IncNavierStokesSolver/InputFiles/</code>.
 The figures below show the numerical results for the  velocity components and the vorticity using a 16 quadrilaterals mesh, 
 an expansion order of P=11 and the third order time-integration scheme. 
 The exact solution has been used to set Dirichlet boundary conditions on the edges. Using a time-step \f$\Delta t=0.001\f$
 the L2-error obtained is \f$10^{-13}\f$ for the velocity field and \f$10^{-12}\f$ for the pressure.
 
 \image html DecVort_u.png 
 \image html DecVort_v.png 
 \image html DecVort_w.png 
 
 \subsection subsectionShedding Flow Past a Cylinder 2D
 
 As an illustrative example to demonstrate the capability of the solver to solve more complex fluid dynamics problems, 
 we compute the two-dimensional flow past a circular cylinder in a free stream. 
 The figure below shows the vorticity field of a flow past a cylinder with a Reynolds number 
 \f$Re=100\f$ using the third order time integration scheme. 
 
 \image html vxshed.png
 
 Now, consider the following configuration:
 
 \image html sphereschematics.png
 
 Consider an infinitely long cylinder of diameter \f$D\f$ placed midway between two parallel plates 
 which are a distance \f$H\f$ apart, as shown in the figure above.
 A no-slip condition is imposed on the plates and the cylinder boundary, 
 a parabolic inflow profile is specified as the inflow condition with a maximum inlet speed of \f$U_{max}=1\f$ and 
 at the outflow boundary homogeneous Neumann boundary conditions are imposed.  
 The Reynolds number for the flow around a cylinder is defined by
 \f$Re=U_{max}D/\nu\f$,  where \f$\nu\f$ is the kinematic viscosity. The figures below show the mesh configuration and 
 contour plots of the computed velocity components for a Reynolds number of 100 for a  
 mesh consisting of 856 quadrilateral elements and a polynomial order P=4 after a physical time of 150s. 
 
 \image html vortexsheddingnewmesh.png
 
 \image html vortexsheddingnewvelu.png "Contour plot of velocity component u"
 
 \image html vortexsheddingnewvelv.png "Contour plot of velocity component v"
 
 We chose a kinematic viscosity of 0.01 and the first order time integration scheme with a timestep of 0.001. 
 We can observe vortices that are shedding periodically from the cylinder. 
 The vortex shedding behaviour can be observed nicely in the following two videos.
 The first video shows the evolution of the velocity component \f$u\f$ and \f$v\f$ and the velocity vector length.
 The second video shows a comparison of the vortex shedding behaviour for different Reynolds number. 
 For a Reynolds number of 50 no vortex shedding is observed. The vortex shedding frequency increases with increasing Reynolds number.
 
 \image html vortexstreetvelocitycomponentsoneperiodsmall.gif "Video of contour plots for the velocity components and the velocity vector lenght"
 
 \image html vortexstreetreoneperiodsmall.gif "Video of contour plots for the velocity componen u for a range of Reynolds numbers"
 
\section sectionRef References
 
 [1] G.E. Karniadakis and S.J. Sherwin: <i>Spectral/hp Element Methods for Computational Fluid Dynamics</i>, Oxford Science Publications, 2005<br />
 [2] M. Dubiner: <i>Spectral methods on triangles and other domains</i>, J. Sci. Comp., 6, 345-390, 1991 <br />
 [3] G. E. Karniadakis, M. Israeli, and S. A. Orszag: <i>High-order splitting methodsfor the incompressible Navier-Stokes equations</i>, J. Comput. Phys., 97, 414-443, 1991.
 
 

*/	
	
}
