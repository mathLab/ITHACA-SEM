/**
 * \page pageGeneralLinearMethods General Linear Methods
 *
 *  The implementation of time integration schemes in Nektar++
 *
 * \section sectionGeneralLinearMethodsIntroduction Introduction
 *
 * General linear methods can be considered as the
 * generalisation of a broad range of different numerical
 * methods for ordinary differential equations.  They were
 * introduced by Butcher and they provide a unified
 * formulation for traditional methods such as the Runge-Kutta
 * methods and the linear multi-step methods.  From an
 * implementational point of view, this means that all these
 * numerical methods can be abstracted in a similar way. As
 * this allows a high level of generality, it is chosen in
 * Nektar++ to cast all time integration schemes in the
 * framework of general linear methods.
 *
 *
 * For background information about general linear methods, please consult the
 * following references:<BR>
 * [1] Butcher, J.C. (2006) <em>General linear methods</em>, Acta Numerica 15, 157-256 <BR>
 * [2] http://www.math.auckland.ac.nz/~butcher/conferences.html
 *
 *
 * \section  sectionGeneralLinearMethods General linear methods
 *
 * The standard initial value problem can written in the form
 * \f[
 * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y}),\quad \boldsymbol{y}(t_0)=\boldsymbol{y}_0
 * \f]
 * where \f$\boldsymbol{y}\f$ is a vector containing the variable (or an array of array contianing the variables). <BR>
 * In the formulation of general linear methods, it is more convenient to consider the ODE in autonomous form, i.e.
 * \f[
 * \frac{d\boldsymbol{\hat{y}}}{dt}=\boldsymbol{\hat{f}}(\boldsymbol{\hat{y}}),\quad \boldsymbol{\hat{y}}(t_0)=
 * \boldsymbol{\hat{y}}_0.
 * \f]
 *
 * \subsection subsectionGeneralLinearMethodsFomulation Formulation
 * Suppose the governing differential equation is given in autonomous form,
 * the \f$n^{th}\f$ step of the general linear method comprising
 *  -  \f$r\f$ steps (as in a multi-step method)
 *  -  \f$s\f$ stages (as in a Runge-Kutta method)
 *
 * is formulated as:
 * \f{eqnarray*}
 * \boldsymbol{Y}_i = \Delta t\sum_{j=0}^{s-1}a_{ij}\boldsymbol{F}_j+\sum_{j=0}^{r-1}u_{ij}
 * \boldsymbol{\hat{y}}_{j}^{[n-1]}, \qquad i=0,1,\ldots,s-1 \\
 * \boldsymbol{\hat{y}}_{i}^{[n]}=\Delta t\sum_{j=0}^{s-1}b_{ij}\boldsymbol{F}_j+\sum_{j=0}^{r-1}v_{ij}
 * \boldsymbol{\hat{y}}_{j}^{[n-1]}, \qquad i=0,1,\ldots,r-1
 * \f}
 * where \f$\boldsymbol{Y}_i\f$ are referred to as the stage values and \f$\boldsymbol{F}_i\f$ as
 * the stage derivatives. Both quantities are related by the differential equation:
 * \f{displaymath}
 * \boldsymbol{F}_i = \boldsymbol{\hat{f}}(\boldsymbol{Y}_i).
 * \f}
 * The matrices \f$A=[a_{ij}]\f$, \f$U=[u_{ij}]\f$, \f$B=[b_{ij}]\f$, \f$V=[v_{ij}]\f$ are
 * characteristic of a specific method. Each scheme can then be uniquely defined by a partioned 
 * \f$(s+r)\times(s+r)\f$ matrix
 * \f{displaymath}
 * \left[ \begin{array}{cc}
 * A & U\\
 * B & V
 * \end{array}\right]
 * \f}
 *
 *
 * \subsection subsectionGeneralLinearMethodsMatrixNotation Matrix notation
 *
 * Adopting the notation:
 * \f{displaymath}
 * \boldsymbol{\hat{y}}^{[n-1]}=
 * \left[\begin{array}{c}
 * \boldsymbol{\hat{y}}^{[n-1]}_0\\
 * \boldsymbol{\hat{y}}^{[n-1]}_1\\
 * \vdots\\
 * \boldsymbol{\hat{y}}^{[n-1]}_{r-1}
 * \end{array}\right],\quad
 * \boldsymbol{\hat{y}}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{\hat{y}}^{[n]}_0\\
 * \boldsymbol{\hat{y}}^{[n]}_1\\
 * \vdots\\
 * \boldsymbol{\hat{y}}^{[n]}_{r-1}
 * \end{array}\right],\quad
 * \boldsymbol{Y}=
 * \left[\begin{array}{c}
 * \boldsymbol{Y}_0\\
 * \boldsymbol{Y}_1\\
 * \vdots\\
 * \boldsymbol{Y}_{s-1}
 * \end{array}\right],\quad
 * \boldsymbol{F}=
 * \left[\begin{array}{c}
 * \boldsymbol{F}_0\\
 * \boldsymbol{F}_1\\
 * \vdots\\
 * \boldsymbol{F}_{s-1}
 * \end{array}\right]\quad
 * \f}
 * the general linear method can be written more compactly in the following form:
 * \f{displaymath}
 * \left[ \begin{array}{c}
 * \boldsymbol{Y}\\
 * \boldsymbol{\hat{y}}^{[n]}
 * \end{array}\right] =
 * \left[ \begin{array}{cc}
 * A\otimes I_N & U\otimes I_N \\
 * B\otimes I_N & V\otimes I_N
 * \end{array}\right]
 * \left[ \begin{array}{c}
 * \Delta t\boldsymbol{F}\\
 * \boldsymbol{\hat{y}}^{[n-1]}
 * \end{array}\right]
 * \f}
 * where \f$I_N\f$ is the identity matrix of dimension \f$N\times N\f$.
 *
 *
 * \subsection subsectionGeneralLinearMethodsEvaluation Evaluation of an General Linear Method
 *
 * Although the General linear method is essentially presented for ODE's in its autonomous form, in
 * Nektar++ it will be used to solve ODE's formulated in non-autonomous form.
 * Given the ODE,
 * \f[
 * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y}),\quad \boldsymbol{y}(t_0)=\boldsymbol{y}_0
 * \f]
 * a single step of General Linear Method can then be evaluated in the following way:<BR><BR>
 * - <b>input</b>
 *   - \f$\boldsymbol{y}^{[n-1]}\f$, i.e. the \f$r\f$ subvectors comprising \f$\boldsymbol{y}^{[n-1]}_i\f$
 *   - \f$t^{[n-1]}\f$, i.e. the equivalent of \f$\boldsymbol{y}^{[n-1]}\f$ for the time variable \f$t\f$
 * - <b>step 1:</b> The stage values \f$\boldsymbol{Y}_i\f$, \f$\boldsymbol{T}_i\f$
 *   and the stage derivatives \f$\boldsymbol{F}_i\f$ are calculated through the relations:
 * \f{eqnarray*}
 * \boldsymbol{Y}_i &=& \Delta t\sum_{j=0}^{s-1}a_{ij}\boldsymbol{F}_j+\sum_{j=0}^{r-1}u_{ij}\boldsymbol{y}_{j}^{[n-1]}, \qquad i=0,1,\ldots,s-1\\
 * T_i &=& \Delta t\sum_{j=0}^{s-1}a_{ij}+\sum_{j=0}^{r-1}u_{ij}
 * t_{j}^{[n-1]}, \qquad i=0,1,\ldots,s-1\\
 * \boldsymbol{F}_i &=& f(T_i,\boldsymbol{Y}_i), \qquad i=0,1,\ldots,s-1
 * \f}
 * - <b>step 2:</b> The approximation at the new time level \f$\boldsymbol{y}^{[n]}\f$ is calculated
 * as a linear combination of the stage derivatives \f$\boldsymbol{F}_i\f$ and the input vector \f$\boldsymbol{y}^{[n-1]}\f$.
 * In addition, the time vector \f$t^{[n]}\f$ is also updated
 * \f{eqnarray*}
 * \boldsymbol{y}_{i}^{[n]}&=&\Delta t\sum_{j=0}^{s-1}b_{ij}\boldsymbol{F}_j+
 * \sum_{j=0}^{r-1}v_{ij}\boldsymbol{y}_{j}^{[n-1]}, \qquad i=0,1,\ldots,r-1\\
 * t_{i}^{[n]}&=&\Delta t\sum_{j=0}^{s-1}b_{ij}+
 * \sum_{j=0}^{r-1}v_{ij}t_{j}^{[n-1]}, \qquad i=0,1,\ldots,r-1\\
 * \f}
 * - <b>output</b>
 *   - \f$\boldsymbol{y}^{[n]}\f$, i.e. the \f$r\f$ subvectors comprising \f$\boldsymbol{y}^{[n]}_i\f$.
 *     \f$\boldsymbol{y}^{[n-1]}_0\f$ corresponds to the actual approximation at the new time level.
 *   - \f$t^{[n]}\f$ where \f$t^{[n]}_0\f$ is equal to the new time level \f$t+\Delta t\f$
 *
 *
 * \section sectionGeneralLinearMethodsNektar General Linear Methods in Nektar++
 *
 * In Nektar++, we do not use the standard General Linear Methods formulation. As mentioned above,
 * we will use the non-autonomous form as the explicit treatment of the time variable \f$t\f$
 * allows for more flexibility. 
 *
 * For a detailed describtion of the formulation and a deeper insight of the numerical method see:
 * <BR><BR>
 * Peter E.J. Vos, Claes Eskilsson, Alessandro Bolis, Sehun Chun,Robert M. Kirby & Spencer J. Sherwin, (2011),<BR> 
 * <em>A generic framework for time-stepping partial differential equations (PDEs): <BR>
 * general linear methods, object-oriented implementation and application to fluid problems</em>, <BR>
 * International Journal of Computational Fluid Dynamics, 25:3, 107-125
 * <BR>
 * 
 * \subsection sectionGeneralLinearMethodsNektarSchemes Type of Time Integration Schemes
 *
 * Nektar++ contains various classes and methods which implement the concept of the General
 * Linear Methods. This toolbox is capable of numerically solving the generalised ODE using a broad range
 * of different time-stepping methods. We distinguish the following types of General Linear Methods:
 * - <b>Formally Explicit Methods</b><BR>
 *   These types of methods are considered explicit from an ODE point of view. They are characterised by a lower
 *   triangular coefficient matrix \f$A\f$, i.e. \f$a_{ij} = 0\f$ for \f$j\geq i\f$. To avoid confusion, we
 *   make a further distinction:
 *   - <i>direct explicit method</i>: Only forward operators are required.
 *   - <i>indirect explicit method</i>: The inverse operator is required.
 * - <b>Diagonally Implicit Methods</b><BR>
 *   Compared to the explicit methods, the coefficient matrix \f$A\f$ has now non-zero entries on the diagonal.
 *   This means that each stage value depend on the stage derivative at the same stage, requiring an implicit
 *   step. However, the calculation of the different stage values is still uncoupled. Best known are the
 *   DIRK schemes.
 * - <b>IMEX schemes</b><BR>
 * These schemes support the concept of being able to split
 * right hand forcing term into an explicit and
 * implicit component. This is useful in advection diffusion
 * type problems where the advection is handled explicity and
 * the diffusion is handled implicit.
 * - <b>Fully Implicit Methods Methods</b><BR>
 *   The coefficient matrix has a non-zero upper triangular part. The calculation of all stages values is fully coupled.
 *
 * The aim in Nektar++ is to fully support the first three
 * types of General Linear Methods.  Fully implicit methods
 * are currently not implemented.
 *
 * \subsection sectionGeneralLinearMethodsNektarSchemesHowTo How to use
 * The goal of abstracting the concept of General Linear Methods is to provide users with a single interface
 * for time-stepping, independent of the chosen method. The TimeIntegrationScheme class allow the user to numerically
 * integrate ODE's using high-order complex schemes, as if it were done using the forward euler method.
 * Switching between time-stepping schemes should be as easy as changing a parameter in an input file.
 * The only thing the user should provide, is an implementation of the left and right hand side operation
 * of the generalised ODE to be solved.
 *
 * To introduce the implementation of time stepping schemes in Nektar++, consider the following example:
 * \code
 * NekDouble timestepsize = 0.1;
 * NekDouble time         = 0.0;
 * 
 * Array<OneD, Array<OneD, NekDouble> >  y;
 * 
 * LibUtilities::TimeIntegrationSchemeOperators ode;
 *
 * ode.DefineImplicitSolve (&function1, &object1);
 * ode.DefineOdeRhs        (&function2, &object2);
 *
 *
 * LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eClassicalRungeKutta4);
 * LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
 *
 * LibUtilities::TimeIntegrationSolutionSharedPtr y_0;
 *
 * y_0 = IntScheme->InitializeScheme(timestepsize,y,time,ode);
 *
 * for(int n = 0; n < nsteps; ++n)
 * {
 *  y = IntScheme->TimeIntegrate(timestep,y_0,ode);
 *  time += timestepsize;
 * }
 * \endcode
 * <BR>
 * We can distinguish four different sections in the code above:
 * - <b>Initialisation</b>
 * \code
 * NekDouble timestepsize = 0.1;
 * NekDouble time         = 0.0;
 * Array<OneD, Array<OneD, NekDouble> >  y;
 * \endcode
 * The vector <code>y</code> will be used to store the solution at the end of each time-step.
 * It corresponds to the vector \f$\boldsymbol{y}^{[n]}_0\f$ defined above.
 * It is an Array of Arrays where the first dimension
 * corresponds to the number of variables (eg. u,v,w) and
 * the second dimension corresponds to the number length of the variables
 * (e.g. the number of modes or the number of physical points).<BR>
 * \code
 * LibUtilities::TimeIntegrationSchemeOperators ode;
 *
 * ode.DefineImplicitSolve (&function1, &object1);
 * ode.DefineOdeRhs        (&function2, &object2);
 * \endcode
 * The variable <code>ode</code> is an object containig the methods. A class representing a PDE eqaution (or system of equations)
 * should have some a series of functions representing the implicit/explicit part of the method, which represents the reduction of the PDE/s to a system of ODEs. 
 * The user should make sure this class contains a proper implementation of these necessary methods. 
 * <code>&function1</code> is a functor, i.e. a pointer
 * to a function where the method is implemented. <code>&object1</code> is a pointer to the object, i.e. the class, where the function/method (<code>&function1</code>) is implemented.
 * For more information about these methods, click \ref sectionGeneralLinearMethodsNektarSchemesRequiredMethods "here", were the hypotical class where the methods/functions
 * are implemented is called <code>Foo</code>.
 
 * - <b>Loading the time integration scheme</b>
 * \code
 * LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eClassicalRungeKutta4);
 * LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
 * \endcode
 * First, a key which uniquely defines the desired time integration scheme is created. Then, Using this key,
 * the TimeIntegrationSchemeManager is called which will return you the requested time integration scheme.
 * The concept of a time integration scheme is abstracted as the
 * class TimeIntegrationScheme. Objects of this class will hold the proper
 * integration coefficients and know how to perform integration.
 * - <b>Initialising the time integration scheme</b>
 * \code
 * LibUtilities::TimeIntegrationSolutionSharedPtr y_0;
 *
 * y_0 = IntScheme->InitializeScheme(timestepsize,y,time,ode);
 * \endcode
 * Given the initial solution (stored in <code>y</code>) and some additional parameters, this method constructs and returns an object <code>y_0</code> which is
 * the abstraction of the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$. This abstraction is done through the class
 * <code>TimeIntegrationSolution</code>. This initialisation is essential in case of multi-step schemes, where the vector
 * \f$\boldsymbol{y}^{[n]}\f$ consist of more than one entry. The object <code>y_0</code> can than later be passed to the actual
 * integration method, as it contains all necessary input.
 * - <b>Perform the time integration</b>
 * \code
 * for(int n = 0; n < nsteps; ++n)
 * {
 *  y = IntScheme->TimeIntegrate(timestep,y_0,ode);
 *  time += timestepsize;
 * } 
 * \endcode
 * The actual time integration can now be done by looping over the functions call above.
 * This function updates the object <code>y_0</code> every time step to hold the vectors
 * \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$ at every time level <code>n</code>. In addition, it also
 * returns the actual solution \f$\boldsymbol{y}^{[n]}_0\f$ (which in fact is also embedded in the
 * object <code>y_0</code>)
 *
 *
 * \subsection sectionGeneralLinearMethodsNektarSchemesRequiredMethods Required Methods
 * An object <code>bar</code> of the class <code>Foo</code> should be passed to the TimeIntegrationScheme
 * class routines. This class <code>Foo</code> should have the following public member functions:
 * - <code>Foo::DoOdeRhs</code>
 *   \code
 *   void Foo::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
 *                          Array<OneD,       Array<OneD, NekDouble> >& outarray,
 *                    const NekDouble time);
 *   \endcode
 *   This function should which evaluate the explict part of the right hand side operator
 *   of the generalised ODE.
 *   - Input Parameters
 *     - <code>inarray</code>: the vector \f$\boldsymbol{y}\f$
 *     - <code>time</code>: the time \f$t\f$
 *   - Output Parameters
 *     - <code>outarray</code>: the result of the right hand side operator
 *     .
 *   .
 *   Both <code>inarray</code> and <code>outarray</code> are Array of Arrays where the first dimension
 *   corresponds to the number of variables (eg. u,v,w) and
 *   the second dimension corresponds to the number length of the variables
 *   (e.g. the number of modes).
 *
 * - <code>Foo::DoImplicitSolve</code>
 *   \code
 *   void Foo::DoImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
 *                                Array<OneD,       Array<OneD, NekDouble> >& outarray,
 *                          const NekDouble time,
 *                          const NekDouble lambda);
 *   \endcode
 * This function should which evaluate the implict part of the right hand side operator
 * of the generalised ODE.
 *   - Input Parameters
 *     - <code>inarray</code>: the vector \f$\boldsymbol{b}\f$
 *     - <code>time</code>: the time \f$t\f$
 *     - <code>lambda</code>: the coefficient \f$\lambda\f$
 *   - Output Parameters
 *     - <code>outarray</code>: the result \f$\boldsymbol{y}\f$
 *     .
 *   .
 *   Both <code>inarray</code> and <code>outarray</code> are Array of Arrays where the first dimension
 *   corresponds to the number of variables (eg. u,v,w) and
 *   the second dimension corresponds to the number length of the variables
 *   (e.g. the number of modes).
 *
 * \subsection sectionGeneralLinearMethodsImplementationEssentialBC Strongly imposed essential boundary conditions.
 
 * Dirichlet boundary conditions can be strongly imposed by lifting the known Dirichlet solution.
 * This is equivalent to decompose the approximate solution \f$y\f$ into an known lifted function,
 * \f$y^{\mathcal{D}}\f$, which satisfies the Dirichlet boundary conditions, and an unknown
 * homogeneous function, \f$y^{\mathcal{D}}\f$, which
 * is zero on the Dirichlet boundaries, i.e.
 * \f[
 * y = y^{\mathcal{D}} + y^{\mathcal{H}}
 * \f]
 * In a Finite Element discretisation, this corresponds to splitting the solution vector of coefficients \f$\boldsymbol{y}\f$
 * into the known Dirichlet degrees of freedom \f$\boldsymbol{y}^{\mathcal{D}}\f$ and the unknown homogeneous
 * degrees of freedom \f$\boldsymbol{y}^{\mathcal{H}}\f$. If ordering the known coefficients first, this corresponds to:
 * \f[
 * \boldsymbol{y} = \left[ \begin{array}{c}
 * \boldsymbol{y}^{\mathcal{D}} \\
 * \boldsymbol{y}^{\mathcal{H}} \end{array} \right]
 * \f]
 * The generalised formulation of the General Linear Method (i.e. the introduction of a left hand side operator)
 * allows for an easier treatment of these types of boundary conditions. To better appreciate this, consider the
 * equation for the stage values for an explicit general linear method where both the left and right hand side operator
 * are linear operators,
 * i.e. they can be represented by a matrix.
 * \f[
 * \boldsymbol{M}\boldsymbol{Y}_i = \Delta t\sum_{j=0}^{i-1}a_{ij}\boldsymbol{L}\boldsymbol{Y}_j+\sum_{j=0}^{r-1}u_{ij}\boldsymbol{y}_{j}^{[n-1]}, \qquad i=0,1,\ldots,s-1
 * \f]
 * In case of a lifted known solution, this can be written as:
 * \f[
 * \left[ \begin{array}{cc}
 * \boldsymbol{M}^{\mathcal{DD}} & \boldsymbol{M}^{\mathcal{DH}} \\
 * \boldsymbol{M}^{\mathcal{HD}} & \boldsymbol{M}^{\mathcal{HH}} \end{array} \right]
 * \left[ \begin{array}{c}
 * \boldsymbol{Y}^{\mathcal{D}}_i \\
 * \boldsymbol{Y}^{\mathcal{H}}_i \end{array} \right]
 * = \Delta t\sum_{j=0}^{i-1}a_{ij}
 * \left[ \begin{array}{cc}
 * \boldsymbol{L}^{\mathcal{DD}} & \boldsymbol{L}^{\mathcal{DH}} \\
 * \boldsymbol{L}^{\mathcal{HD}} & \boldsymbol{L}^{\mathcal{HH}} \end{array} \right]
 * \left[ \begin{array}{c}
 * \boldsymbol{Y}^{\mathcal{D}}_j \\
 * \boldsymbol{Y}^{\mathcal{H}}_j \end{array} \right]
 * +\sum_{j=0}^{r-1}u_{ij}
 * \left[ \begin{array}{c}
 * \boldsymbol{y}^{\mathcal{D}[n-1]}_j \\
 * \boldsymbol{y}^{\mathcal{H}[n-1]}_j \end{array} \right], \qquad i=0,1,\ldots,s-1
 * \f]
 *
 * In order to calculate the stage values correctly, the \ref sectionGeneralLinearMethodsNektarSchemesRequiredMethods
 *  should now be implemented to do the following:
 *
 * - <code>Foo::DoOdeRhs</code>
 *   \code
 *   void Foo::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
 *                          Array<OneD,       Array<OneD, NekDouble> >& outarray,
 *                    const NekDouble time);
 *   \endcode
 *   This function should now evaluate the right hand side operator in order to calculate:
 *   \f[
 *   \left[ \begin{array}{c}
 *   \boldsymbol{b}^{\mathcal{D}} \\
 *   \boldsymbol{b}^{\mathcal{H}} \end{array} \right]
 *    =
 *   \left[ \begin{array}{cc}
 *   \boldsymbol{L}^{\mathcal{DD}} & \boldsymbol{L}^{\mathcal{DH}} \\
 *   \boldsymbol{L}^{\mathcal{HD}} & \boldsymbol{L}^{\mathcal{HH}} \end{array} \right]
 *   \left[ \begin{array}{c}
 *   \boldsymbol{y}^{\mathcal{D}} \\
 *   \boldsymbol{y}^{\mathcal{H}} \end{array} \right]
 *   \f]
 *   Note that only the homogeneous part \f$\boldsymbol{b}^{\mathcal{H}}\f$ will be used to
 *   calculate the stage values, as seen in the <code>Foo::ODElhsSolve</code> method. This means
 *   essentially, only the bottom part of the operation above, i.e.
 *   \f[ \boldsymbol{L}^{\mathcal{HD}}\boldsymbol{y}^{\mathcal{D}} + \boldsymbol{L}^{\mathcal{HH}}\boldsymbol{y}^{\mathcal{H}} \f]
 *   is required. However, sometimes it might be more convenient to use/implement routines for the
 *   <code>Foo::DoOdeRhs</code> method that also calculate \f$\boldsymbol{b}^{\mathcal{D}}\f$.
 *
 * - <code>Foo::DoImplicitSolve</code>
 *   \code
 *   void Foo::DoImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
 *                                Array<OneD,       Array<OneD, NekDouble> >& outarray,
 *                          const NekDouble time,
 *                          const NekDouble lambda);
 *   \endcode
 *   This method, needed for diagonally implicit methods, should now solve the system:
 *   \f[
 *   \left(\left[ \begin{array}{cc}
 *   \boldsymbol{M}^{\mathcal{DD}} & \boldsymbol{M}^{\mathcal{DH}} \\
 *   \boldsymbol{M}^{\mathcal{HD}} & \boldsymbol{M}^{\mathcal{HH}} \end{array} \right]
 *   - \lambda \left[ \begin{array}{cc}
 *   \boldsymbol{L}^{\mathcal{DD}} & \boldsymbol{L}^{\mathcal{DH}} \\
 *   \boldsymbol{L}^{\mathcal{HD}} & \boldsymbol{L}^{\mathcal{HH}} \end{array} \right]\right)
 *   \left[ \begin{array}{c}
 *   \boldsymbol{y}^{\mathcal{D}} \\
 *   \boldsymbol{y}^{\mathcal{H}} \end{array} \right]
 *   =
 *   \left[ \begin{array}{cc}
 *   \boldsymbol{H}^{\mathcal{DD}} & \boldsymbol{H}^{\mathcal{DH}} \\
 *   \boldsymbol{H}^{\mathcal{HD}} & \boldsymbol{H}^{\mathcal{HH}} \end{array} \right]
 *   \left[ \begin{array}{c}
 *   \boldsymbol{y}^{\mathcal{D}} \\
 *   \boldsymbol{y}^{\mathcal{H}} \end{array} \right]
 *   =
 *   \left[ \begin{array}{c}
 *   \boldsymbol{b}^{\mathcal{D}} \\
 *   \boldsymbol{b}^{\mathcal{H}} \end{array} \right]
 *   \f]
 *   for the unknown vector \f$\boldsymbol{y}\f$. This can be done in three steps:
 *   -# Set the known solution \f$\boldsymbol{y}^{\mathcal{D}}\f$
 *   -# Calculate the modified right hand side term:
 *      \f[  \boldsymbol{b}^{\mathcal{H}} - \boldsymbol{H}^{\mathcal{HD}}\boldsymbol{y}^{\mathcal{D}} \f]
 *   -# Solve the system below for the unknown \f$\boldsymbol{y}^{\mathcal{H}}\f$,
 *      \f[ \boldsymbol{H}^{\mathcal{HH}}\boldsymbol{y}^{\mathcal{H}} = \boldsymbol{b}^{\mathcal{H}} - \boldsymbol{H}^{\mathcal{HD}}\boldsymbol{y}^{\mathcal{D}} \f]
 *
 *
 * \subsection sectionGeneralLinearMethodsImplementationSchemes Implemented integration schemes
 * Currently the time integration schemes below are implemented in Nektar++. We will list their coefficients here
 * and we will also indicate what the auxiliary parameters \f$\boldsymbol{\hat{y}}_{j}^{[n]}\f$ (\f$j=1,\ldots,r-1\f$)
 * of the multistep methods represent:
 * - Forward Euler (or 1st order Adams Bashfort)
 *   - enum value: <code>eForwardEuler</code>, <code>eAdamsBashforthOrder1</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{c|c}
 * 0 & 1 \\
 * \hline
 * 1 & 1
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - Backward Euler (or 1st order Adams Moulton)
 *   - enum value: <code>eBackwardEuler</code>, <code>eAdamsMoultonOrder1</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{c|c}
 * 1 & 1 \\
 * \hline
 * 1 & 1
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - 2nd order Adams Bashforth
 *   - enum value: <code>eAdamsBashforthOrder2</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{c|cc}
 * 0 & 1 & 0 \\
 * \hline
 * \frac{3}{2} & 1 & \frac{-1}{2}  \\
 * 1 & 0 & 0
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0\\
 * \boldsymbol{y}^{[n]}_1\\
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)\\
 * \Delta t \boldsymbol{l}(t^{n-1},\boldsymbol{y}^{n-1})
 * \end{array}\right]
 * \f]
 * - 1st order IMEX Euler backwards/ Euler Forwards
 *   - enum value: <code>eIMEXOrder1</code>
 *   - coefficients and parameters:
 * \f[
 * \left[ \begin{array}{cc|c}
 * A^{\mathrm{IM}} & A^{\mathrm{EM}} & U \\
 * \hline
 * B^{\mathrm{IM}} & B^{\mathrm{EM}} & V
 * \end{array} \right ] =
 * \left[\begin{array}{cc|c}
 *  \left [ \begin{array}{c} 1 \end{array} \right ] & \left [ \begin{array}{c} 0 \end{array} \right ] &
 *  \left [ \begin{array}{cc} 1 & 1 \end{array} \right ] \\
 * \hline
 *  \left [ \begin{array}{c} 1 \\ 0 \end{array} \right ] & \left [ \begin{array}{c} 0 \\ 1 \end{array} \right ]&
 *  \left [ \begin{array}{cc} 1 & 1 \\ 0 & 0 \end{array} \right ]
 * \end{array}\right] \quad \mathrm{with}\quad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0\\
 * \boldsymbol{y}^{[n]}_1
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^n, \boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - 2nd order IMEX Backward Different Formula & Extrapolation
 *   - enum value: <code>eIMEXOrder2</code>
 *   - coefficients and parameters:
 * \f[
 * \left[ \begin{array}{cc|c}
 * A^{\mathrm{IM}} & A^{\mathrm{EM}} & U \\
 * \hline
 * B^{\mathrm{IM}} & B^{\mathrm{EM}} & V
 * \end{array} \right ] =
 * \left[\begin{array}{cc|c}
 *  \left [ \begin{array}{c} \frac{2}{3} \end{array} \right ] & \left [ \begin{array}{c} 0 \end{array} \right ] &
 *  \left [ \begin{array}{cccc} \frac{4}{3} & -\frac{1}{3} & \frac{4}{3} &  -\frac{2}{3} \end{array} \right ] \\
 * \hline
 *  \left [ \begin{array}{c} \frac{2}{3} \\ 0 \\ 0 \\ 0 \end{array} \right ] & \left [ \begin{array}{c} 0 \\0 \\ 1 \\ 0 \end{array} \right ]&
 *  \left [ \begin{array}{cccc} \frac{4}{3} & -\frac{1}{3} & \frac{4}{3} & -\frac{2}{3} \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0\\ 0 & 0 & 1 & 0 \end{array} \right ]
 * \end{array}\right] \quad \mathrm{with}\quad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0\\
 * \boldsymbol{y}^{[n]}_1 \\
 * \boldsymbol{y}^{[n]}_2 \\
 * \boldsymbol{y}^{[n]}_3
 * \end{array}\right]=
 * \left[\begin{array}{l}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)\\
 * \boldsymbol{m}\left(t^{n-1},\boldsymbol{y}^{n-1}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^n, \boldsymbol{y}^{n}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^{n-1}, \boldsymbol{y}^{n-1}\right)
 * \end{array}\right]
 * \f]
 * (Note: the first two rows are normalised so the coefficient on \f$ \boldsymbol{y}_n^{[n+1]} \f$ is one. In the standard formulation it is 3/2)
 * - 3rdorder IMEX Backward Different Formula & Extrapolation
 *   - enum value: <code>eIMEXOrder3</code>
 *   - coefficients and parameters:
 * \f[
 * \left[ \begin{array}{cc|c}
 * A^{\mathrm{IM}} & A^{\mathrm{EM}} & U \\
 * \hline
 * B^{\mathrm{IM}} & B^{\mathrm{EM}} & V
 * \end{array} \right ] =
 * \left[\begin{array}{cc|c}
 *  \left [ \begin{array}{c} \frac{6}{11} \end{array} \right ] & \left [ \begin{array}{c} 0 \end{array} \right ] &
 *  \left [ \begin{array}{cccccc} \frac{18}{11} & -\frac{9}{11} & \frac{2}{11} & \frac{18}{11} &  -\frac{18}{11} & \frac{6}{11} \end{array} \right ] \\
 * \hline
 *  \left [ \begin{array}{c} \frac{6}{11}\\ 0 \\ 0 \\ 0 \\0 \\0 \end{array} \right ] & \left [ \begin{array}{c} 0 \\0 \\ 0 \\ 1 \\ 0 \\ 0 \end{array} \right ]&
 *  \left [ \begin{array}{cccccc} \frac{18}{11} & -\frac{9}{11} & \frac{2}{11} & \frac{18}{11} & -\frac{18}{11} & \frac{6}{11} \\ 1 & 0 & 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0\\ 0 & 0 & 0 & 1 & 0 & 0 \\ 0 & 0 &  0 & 0 & 1 & 0 \end{array} \right ]
 * \end{array}\right] \quad \mathrm{with}\quad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{l}
 * \boldsymbol{y}^{[n]}_0 \\
 * \boldsymbol{y}^{[n]}_1 \\
 * \boldsymbol{y}^{[n]}_2 \\
 * \boldsymbol{y}^{[n]}_3 \\
 * \boldsymbol{y}^{[n]}_4 \\
 * \boldsymbol{y}^{[n]}_5
 * \end{array}\right]=
 * \left[\begin{array}{l}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)\\
 * \boldsymbol{m}\left(t^{n-1},\boldsymbol{y}^{n-1}\right)\\
 * \boldsymbol{m}\left(t^{n-2},\boldsymbol{y}^{n-2}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^n, \boldsymbol{y}^{n}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^{n-1}, \boldsymbol{y}^{n-1}\right)\\
 * \Delta t \boldsymbol{l}\left ( t^{n-2}, \boldsymbol{y}^{n-2}\right)
 * \end{array}\right]
 * \f]
 * (Note: the first two rows are normalised so the coefficient on \f$ \boldsymbol{y}_n^{[n+1]} \f$ is one. In the standard formulation it is 11/6)
 * - 2nd order Adams Moulton
 *   - enum value: <code>eAdamsMoultonOrder2</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{c|cc}
 * \frac{1}{2} & 1 & \frac{1}{2} \\
 * \hline
 * \frac{1}{2} & 1 & \frac{1}{2} \\
 * 1 & 0 & 0
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0\\
 * \boldsymbol{y}^{[n]}_1\\
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)\\
 * \Delta t \boldsymbol{l}(t^{n},\boldsymbol{y}^{n})
 * \end{array}\right]
 * \f]
 * - the midpoint method
 *   - enum value: <code>eMidpoint</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{cc|c}
 * 0 & 0 & 1 \\
 * \frac{1}{2} & 0 & 1 \\
 * \hline
 * 0 & 1 & 1
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - RK4: the standard fourth order Runge-Kutta scheme
 *   - enum value: <code>eClassicalRungeKutta4</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{cccc|c}
 * 0 & 0 & 0 & 0 & 1 \\
 * \frac{1}{2} & 0 & 0 & 0 & 1 \\
 * 0 & \frac{1}{2} & 0 & 0 & 1 \\
 * 0 & 0 & 1 & 0 & 1 \\
 * \hline
 * \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6} & 1
 * \end{array}\right],\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - 2nd order Diagonally Implicit Runge Kutta (DIRK)
 *   - enum value: <code>eDIRKOrder2</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{cc|c}
 * \lambda & 0 & 1 \\
 * \left(1-\lambda\right) & \lambda & 1 \\
 * \hline
 * \left(1-\lambda\right) & \lambda & 1
 * \end{array}\right]\quad \mathrm{with}\quad \lambda=\frac{2-\sqrt{2}}{2},\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - 3rd order Diagonally Implicit Runge Kutta (DIRK)
 *   - enum value: <code>eDIRKOrder3</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{c|c}
 * A & U \\
 * \hline
 * B & V
 * \end{array}\right] =
 * \left[\begin{array}{ccc|c}
 * \lambda & 0 & 0 & 1 \\
 * \frac{1}{2}\left(1-\lambda\right) & \lambda & 0 & 1 \\
 * \frac{1}{4}\left(-6\lambda^2+16\lambda-1\right) & \frac{1}{4}\left(6\lambda^2-20\lambda+5\right) & \lambda & 1 \\
 * \hline
 * \frac{1}{4}\left(-6\lambda^2+16\lambda-1\right) & \frac{1}{4}\left(6\lambda^2-20\lambda+5\right) & \lambda & 1
 * \end{array}\right]\quad \mathrm{with}\quad \lambda=0.4358665215,\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * - 3rd order L-stable, three stage IMEX DIRK(3,4,3)
 *   - enum value: <code>eIMEXdirk_3_4_3</code>
 *   - coefficients and parameters:
 * \f[
 * \left[\begin{array}{cc|c}
 * A^{\mathrm{IM}} & A^{\mathrm{EM}} & U \\
 * \hline
 * B^{\mathrm{IM}} & B^{\mathrm{EM}} & V
 * \end{array}\right] =
 * \left[\begin{array}{cc|c}
 * \left[\begin{array}{cccc}
 * 0 & 0 & 0 & 0 \\
 * 0 & \lambda & 0 & 0 \\
 * 0 & \frac{1}{2}\left(1-\lambda\right) & \lambda & 0 \\
 * 0 & \frac{1}{4}\left(-6\lambda^2+16\lambda-1\right) & \frac{1}{4}\left(6\lambda^2-20\lambda+5\right) & \lambda
 * \end{array}\right] &
 * \left[\begin{array}{cccc}
 * 0 & 0 & 0 & 0 \\
 * \lambda & 0 & 0 & 0 \\
 * 0.3212788860 & 0.3966543747 & 0 & 0 \\
 * -0.105858296 & 0.5529291479 & 0.5529291479 & 0
 * \end{array}\right] &
 * \left[\begin{array}{c}
 * 1\\
 * 1\\
 * 1\\
 * 1
 * \end{array}\right] \\
 * \hline
 * \left[\begin{array}{cccc}
 * 0 & \frac{1}{4}\left(-6\lambda^2+16\lambda-1\right) & \frac{1}{4}\left(6\lambda^2-20\lambda+5\right) & \lambda
 * \end{array}\right] &
 * \left[\begin{array}{cccc}
 * 0 & \frac{1}{4}\left(-6\lambda^2+16\lambda-1\right) & \frac{1}{4}\left(6\lambda^2-20\lambda+5\right) & \lambda
 * \end{array}\right] &
 * \left[\begin{array}{c}
 * 1
 * \end{array}\right]
 * \end{array}\right]\quad \mathrm{with}\quad \lambda=0.4358665215,\qquad
 * \boldsymbol{y}^{[n]}=
 * \left[\begin{array}{c}
 * \boldsymbol{y}^{[n]}_0
 * \end{array}\right]=
 * \left[\begin{array}{c}
 * \boldsymbol{m}\left(t^n,\boldsymbol{y}^{n}\right)
 * \end{array}\right]
 * \f]
 * \subsection sectionGeneralLinearMethodsImplementationAddMethods How to add a method
 * To add a new time integration scheme, follow the steps below:
 * - Choose a name for the method and add it to the ::TimeIntegrationMethod enum list.
 * - Populate the switch statement in the TimeIntegrationScheme constructor
 *   with the coefficients of the new method.
 * - Use ( or modify) the function TimeIntegrationScheme::InitializeScheme
 *   to select (or implement) a proper initialisation strategy for the method.
 * - Add documentation for the method (especially indicating what the auxiliary parameters
 *   of the input and output vectors of the multi-step method represent)
 */
