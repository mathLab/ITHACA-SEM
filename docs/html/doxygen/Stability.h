namespace Nektar {
/**
* \page pageStability Stability Analysis
 
 The starting point for the stability analysis of the incompressible Navier-Stokes equations is 
 a proper base flow \f$\boldsymbol{U}\f$ (steady or unsteady) which a solution of the Navier-Stokes equation.
 The base flow is commonly obtained as a DNS solution of the incompressible Navier-Stokes equations (associate with pressure field \f$P\f$).
 
 To analyse the system stability, we study the evolution of an infinitesimal perturbation of the base flow (\f$\boldsymbol{u'}\f$).
 The equations which govern the perturbation evolution are simply the Navier-Stokes equations where \f$\boldsymbol{u} = \boldsymbol{U} +\epsilon \boldsymbol{u'}\f$
 and \f$p = P +\epsilon p'\f$ (\f$p'\f$ is the pressure perturbation).
 If we keep the lower order terms respect to \f$\epsilon\f$, i.e. the linear terms, we obtain the linearized Navier-Stokes equations, which can be written
 as:
 
 \f[
 \frac{\partial \boldsymbol{u'}}{\partial t} = -(\boldsymbol{U}\cdot\nabla)\boldsymbol{u'} - (\boldsymbol{u'}\cdot\nabla)\boldsymbol{U}-\nabla p' + Re^{-1}\Delta \boldsymbol{u'}
 \f]
 \f[
 \nabla \cdot \boldsymbol{u'}=0
 \f]
 
 The linear stability is studied starting from a linear evolution operator \f$\mathcal{A}\f$  which evolves the perturbation forward in time:
 
 \f[
 \boldsymbol{u'}(t)=\mathcal{A}\boldsymbol{u'}(0)
 \f]
 
 The eigenvalue problem is
 
 \f[
\mathcal{A}(T)\boldsymbol{\tilde{u}_j}=\mu_j\boldsymbol{\tilde{u}_j}
 \f]
 
 with
 
 \f[
 \mu_j=exp(\lambda_j T)
 \f]
 
 Inside Nektar++ the stability analysis can be perfomed using both the direct solver and the splitting scheme.
 We have two numerical methods to work out the eigenvalues/eigenvectors of the problem, which are:
 - The Implicitly Restarted Arnoldi Method (see <a href="http://www.caam.rice.edu/software/ARPACK/UG/ug.html">ARPACK Users' Guide</a>).
 - A Modified Arnoldi Method (see [1]).
 
 The switch between those two approaches can be done in Nektar++ from the session file and it is managed by the Driver class in 
 the folder <code>Nektar++/solvers/Auxiliary</code>.
 

 A typical set of parameters for using for exaple ARPACK and the Velocity Correction Scheme is:
 @code
 <SOLVERINFO>
   <I PROPERTY="EQTYPE" VALUE="UnsteadyNavierStokes"/>
   <I PROPERTY="SolverType"  VALUE="VelocityCorrectionScheme"/>
   <I PROPERTY="AdvectionForm" VALUE="Linearised"/>
   <I PROPERTY="Projection" VALUE="Galerkin"/>
   <I PROPERTY="TimeIntegrationMethod" VALUE="IMEXOrder3"/>
   <I PROPERTY="InitialVector" VALUE="Random"/>
   <I PROPERTY="Driver" VALUE="Arpack"/>
 </SOLVERINFO>
 @endcode
 if we want to use for example the Coupled Solver and the Modified Arnoldi Method:
 @code
 <SOLVERINFO>
 <I PROPERTY="EQTYPE" VALUE="UnsteadyNavierStokes"/>
 <I PROPERTY="SolverType"  VALUE="CoupledLinearisedNS"/>
 <I PROPERTY="AdvectionForm" VALUE="Linearised"/>
 <I PROPERTY="Projection" VALUE="Galerkin"/>
 <I PROPERTY="TimeIntegrationMethod" VALUE="IMEXOrder3"/>
 <I PROPERTY="InitialVector" VALUE="Random"/>
 <I PROPERTY="Driver" VALUE="ModifiedArnoldi"/>
 </SOLVERINFO>
 @endcode
 
 Here some input parameters:
 @code
 <PARAMETERS>
   <P> TimeStep      = 0.002  </P>
   <P> NumSteps      = 70     </P>
   <P> IO_CheckSteps = 1000   </P>
   <P> IO_InfoSteps  = 20     </P>
   <P> Re            = 7500   </P>
   <P> Kinvis        = 1.0/Re </P>
   <P> kdim          = 16     </P> 
   <P> nvec          = 2      </P> 
   <P> evtol         = 1e-6   </P> 
 </PARAMETERS>
 @endcode
 
 As mentioned before, the Driver class manage the eigenvales calculations:
 @code
 <I PROPERTY="Driver" VALUE="ModifiedArnoldi"/>
 @endcode
 or
 @code
 <I PROPERTY="Driver" VALUE="Arpack"/>
 @endcode
 
 The initial vector, to start the Arnoldi Method, in both case can be random
 @code
 <I PROPERTY="InitialVector" VALUE="Random"/>
 @endcode
 or can be specified by the user removing the previous tag and specifying in the initail condition
 the values of the vector.
 
 The base flow does not need to be specified in the session file but it need to be located 
 in the same folder where you are running the simulation.
 It does have to have the same name of the .xml file, but the extension must be .bse (it is as usually an fld-type file).
 
 In the parameter list, further option can be specified like:
 
 The Krylov space dimension:
 @code
  <P> kdim = 16     </P>
 @endcode
 
 The number of eigenvalues/eigenvectors you want:
 @code
  <P> nvec = 2      </P>
 @endcode
 
 The eigenvalues/eigenvector tollerance:
 @code
 <P> evtol = 1e-6   </P>
 @endcode
 
 The Arnoldi Method maximum iterations number (defualt is 500):
 @code
 <P> nit   = 1000   </P> 
 @endcode
 
 Here some results for the session file <code>ChanStability.xml</code> which can be found
 in <code>Nektar++/regressionTests/Solvers/IncNavierStokesSolver/InputFiles</code>.
 
 The file can be run as usually
 @code
 ./IncNavierStokesSolver ChanStability.xml
 @endcode
 taking care of copying the .xml, .rst and .bse file in the folder from where you are running the executable.
 
  \image html eigen_u.png "Eigen values (u-velocity component)" 
 
  \image html eigen_v.png "Eigen values (v-velocity component)" 
 
 [1] D. Barkley, M.H.Blackburn and S.J. Sherwin <i>Direct optimal growth analysis for timesteppers</i>, Int. J. Numer. Meth Fluids, 57, 1435-1458, 2008<br />
*/	
	
}
