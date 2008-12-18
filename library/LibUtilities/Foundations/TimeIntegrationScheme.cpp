///////////////////////////////////////////////////////////////////////////////
//
// File TimeIntegrationScheme.cpp
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
// Description: implementation of time integration key class 
//
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/LibUtilities.h>
#include <iostream>
#include <LibUtilities/Foundations/TimeIntegrationScheme.h>
#include <LibUtilities/Foundations/Foundations.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {        
        /**
         * \page pageGeneralLinearMethods General Linear Methods 
         *  the implementation of time integration schemes in Nektar++
         *
         * \section sectionGeneralLinearMethodsIntroduction Introduction
         * General linear methods can be considered as the generalisation of a broad 
         * range of different numerical methods for ordinary differential equations. 
         * They were introduced by Butcher and they provide a unified formulation for 
         * traditional methods such as the Runge-Kutta methods and the linear multi-step methods. 
         * From an implementational point of view, this means that all these numerical methods 
         * can be abstracted in a similar way. As this allows a high level of generality, 
         * it is chosen in Nektar++ to cast all time integration schemes in the framework of 
         * general linear methods.
         *
         *
         * For background information about general linear methods, please consult the 
         * following references:<BR>
         * [1] Butcher, J.C. (2006) <em>General linear methods</em>, Acta Numerica 15, 157-256 <BR>
         * [2] http://www.math.auckland.ac.nz/~butcher/conferences.html
         *
         * \section  sectionGeneralLinearMethods General linear methods
         * The standard initial value problem can written in the form
         * \f[
         * \frac{d\boldsymbol{y}}{dt}=\boldsymbol{f}(t,\boldsymbol{y}),\quad \boldsymbol{y}(t_0)=\boldsymbol{y}_0
         * \f]
         * where \f$\boldsymbol{y}=[y^0\quad y^1\quad  \cdots \quad y^{N-1}]^T\f$. <BR>
         * In the formulation of general linear methods, it is more convenient to consider the ODE in autonomous form, i.e.
         * \f[
         * \frac{d\boldsymbol{\hat{y}}}{dt}=\boldsymbol{\hat{f}}(\boldsymbol{\hat{y}}),\quad \boldsymbol{\hat{y}}(t_0)=
         * \boldsymbol{\hat{y}}_0.
         * \f]
         *
         * When the original problem is presented in non-autonomous form, it can be written as an 
         * autonomous system in the following way:
         * \f[
         * \boldsymbol{\hat{y}}=
         * \left [ \begin{array}{c}
         * t \\
         * \boldsymbol{y}
         * \end{array}\right],\quad
         * \boldsymbol{\hat{f}}(\boldsymbol{\hat{y}})=
         * \left [ \begin{array}{c}
         * 1 \\
         * \boldsymbol{f}(\boldsymbol{\hat{y}})
         * \end{array}\right].
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
         * characteristic of a specific method.
         * Furthermore, the following two quantities serve as input and output of step \f$n\f$:
         * - <b>output</b>: \f$\boldsymbol{\hat{y}}_{j}^{[n]}\f$ the new set of approximations.
         *    - \f$\boldsymbol{\hat{y}}_{0}^{[n]}\f$ refers to the approximation \f$\boldsymbol{\hat{y}}(t+n\Delta t)\f$. 
         *      This is the solution were are looking for.
         *    - \f$\boldsymbol{\hat{y}}_{j}^{[n]}\f$ (\f$j=1,\ldots,r-1\f$) refers to the approximation at
         *       the current time-level of an auxiliary set of parameters inherent to the method.
         * - <b>input</b>: \f$\boldsymbol{\hat{y}}_{j}^{[n-1]}\f$ the set of approximations calculated at step \f$n-1\f$.
         *    - \f$\boldsymbol{\hat{y}}_{0}^{[n-1]}\f$ refers to the approximation of the solution at time 
         *      level \f$n-1\f$ calculated during the step \f$n-1\f$.
         *    - \f$\boldsymbol{\hat{y}}_{j}^{[n-1]}\f$ (\f$j=1,\ldots,r-1\f$) refers to the approximation 
         *      of the auxiliary parameters calculated during the step \f$n-1\f$.
         *
         * \subsection subsectionGeneralLinearMethodsMatrixNotation Matrix notation
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
         * A\otimes I & U\otimes I \\
         * B\otimes I & V\otimes I
         * \end{array}\right]
         * \left[ \begin{array}{c}
         * \Delta t\boldsymbol{F}\\
         * \boldsymbol{\hat{y}}^{[n-1]}
         * \end{array}\right]
         * \f}
         * where \f$I\f$ is the identity matrix of dimension \f$N\times N\f$.
         *
         * \section sectionGeneralLinearMethodsEvaluation Evaluation of an General Linear Method
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
         * \section sectionGeneralLinearMethodsImplementation Implementation in Nektar++
         * \subsection sectionGeneralLinearMethodsImplementationHowTo How to use
         * Nektar++ contains various classes and methods which implement the concept of the General
         * Linear Method. As an introduction of how to use this functionality, consider the following example which implements
         * the explicit integration of a particular ODE.       
         * \code
         * NekDouble timestepsize = 0.1;
         * NekDouble time         = 0.0;
         * Array<OneD, Array<OneD, NekDouble> >  y_t0;
         * Array<OneD, Array<OneD, NekDouble> >  y_0;
         * SomeClass object;
         *        
         * LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eClassicalRungeKutta4);
         * LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
         *  
         * int nInitSteps;
         * LibUtilities::TimeIntegrationSolutionSharedPtr y = IntScheme->InitializeScheme(timestepsize,time,nInitSteps,object,y_t0);
         *
         * for(n = nInitSteps; n < nsteps; ++n)
         * {
         *  y0 = IntScheme->ExplicitIntegration(m_timestep,object,y);
         *  time += timestepsize;
         * }
         * \endcode
         * <BR>
         * We can distinguish four different sections in the code above:
         * - <b>Initialisation</b>
         * \code
         * NekDouble timestepsize = 0.1;
         * NekDouble time         = 0.0;
         * Array<OneD, Array<OneD, NekDouble> >  y_t0;
         * Array<OneD, Array<OneD, NekDouble> >  y_0;
         * Foo bar;
         * \endcode
         * The object <code>y_t0</code> should contain the solution at the initial time t0.
         * The object <code>y_0</code> will later be used to store the solution at every time level.
         * It corresponds to the vector \f$\boldsymbol{y}^{[n]}_0\f$ defined above.
         * Both <code>y_t0</code> and <code>y_0</code> are Array of Arrays where the first dimension
         * corresponds to the number of variables (eg. u,v,w) and
         * the second dimension corresponds to the number of points
         * every variable is defined (eg. the quadrature points).<BR>
         * The variable <code>bar</code> is the object of some class <code>Foo</code>. This class
         * Foo should have some method called <code>ODEforcing</code>:
         * \code    
         * void Foo::ODEforcing(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time)
         * \endcode
         * which calculates the right hand side \f$\boldsymbol{f}(t,\boldsymbol{y})\f$ of the ODE.
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
         * int nInitSteps;
         * LibUtilities::TimeIntegrationSolutionSharedPtr y = IntScheme->InitializeScheme(timestepsize,time,nInitSteps,bar,y_t0);
         * \endcode
         * Given the initial solution and some additional parameters, this method constructs and returns an object <code>y</code> which is
         * the abstraction of the vectors \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$. This abstraction is done through the class
         * <code>TimeIntegrationSolution</code>. This initialisation is essential in case of multi-step schemes, where the vector
         * \f$\boldsymbol{y}^{[n]}\f$ consist of more than one entry. The object <code>y</code> can than later be passed to the actual 
         * integration method, as it contains all necessary input.
         * - <b>Perform the time integration</b>         
         * \code
         * for(n = nInitSteps; n < nsteps; ++n)
         * {
         *  y0 = IntScheme->ExplicitIntegration(m_timestep,object,y);
         *  time += timestepsize;
         * }
         * \endcode
         * The actual explicit time integration can now be done by looping over the function call above.
         * This function updates the object <code>y</code> every time step to hold the vectors
         * \f$\boldsymbol{y}^{[n]}\f$ and \f$t^{[n]}\f$ at every time level <code>n</code>. In addition, it also 
         * returns the actual solution \f$\boldsymbol{y}^{[n]}_0\f$ (which in fact is also embedded in the 
         * object <code>y</code>)
         * \subsection sectionGeneralLinearMethodsImplementationMethods Implemented methods
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
         * \boldsymbol{y}^{n}
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
         * \boldsymbol{y}^{n}
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
         * \boldsymbol{y}^{n}\\
         * \Delta t f(t^{n-1},\boldsymbol{y}^{n-1})
         * \end{array}\right]
         * \f]
         * - 2nd order Adams Moulton
         *   - enum value: <code>eAdamsMoultonOrder2</code>
         *   - coefficients and parameters:
         * \f[
         * \left[\begin{array}{c|c}
         * A & U \\
         * \hline 
         * B & V
         * \end{array}\right] = 
         * \left[\begin{array}{c|c}
         * \frac{1}{2} & 1 \\
         * \hline
         * 1 & 1 
         * \end{array}\right],\qquad
         * \boldsymbol{y}^{[n]}=
         * \left[\begin{array}{c}
         * \boldsymbol{y}^{[n]}_0
         * \end{array}\right]=
         * \left[\begin{array}{c}
         * \boldsymbol{y}^{n}
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
         * \boldsymbol{y}^{n}
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
         * \boldsymbol{y}^{n}
         * \end{array}\right]
         * \f]
         * \subsection sectionGeneralLinearMethodsImplementationAddMethods How to add a method
         * To add a new time integration scheme, follow the steps below:
         * - Choose a name for the method and add it to the ::TimeIntegrationType enum list.
         * - Add a RegisterCreator for the ::TimeIntegrationSchemeManager
         * - Populate the switch statement in the TimeIntegrationScheme constructor 
         *   with the coefficients of the new method.
         * - Populate the switch statement in the function TimeIntegrationScheme::InitializeScheme  
         *   with a proper initialisation strategy for the method.
         * - Add documentation for the method (especially indicating what the auxiliary parameters
         *   of the input and output vectors of the multi-step method represent)
         */


        namespace
        {            
            const bool AdamsBashforthOrder1_Inited = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eAdamsBashforthOrder1), TimeIntegrationScheme::Create);
            const bool AdamsBashforthOrder2_Inited   = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eAdamsBashforthOrder2), TimeIntegrationScheme::Create);
            const bool AdamsMoultonOrder1_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eAdamsMoultonOrder1), TimeIntegrationScheme::Create);
            const bool AdamsMoultonOrder2_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eAdamsMoultonOrder2), TimeIntegrationScheme::Create);
            const bool ClassicalRungeKutta4_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eClassicalRungeKutta4), TimeIntegrationScheme::Create);
            const bool ForwardEuler_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eForwardEuler), TimeIntegrationScheme::Create);
            const bool BackwardEuler_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eBackwardEuler), TimeIntegrationScheme::Create);
            const bool Midpoint_Inited     = TimeIntegrationSchemeManager().RegisterCreator(TimeIntegrationSchemeKey(eMidpoint), TimeIntegrationScheme::Create);
        };

        TimeIntegrationSchemeManagerT &TimeIntegrationSchemeManager(void)
        {
            return Loki::SingletonHolder<TimeIntegrationSchemeManagerT>::Instance();
        }
        

        bool operator==(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {
            return (lhs.m_integrationtype == rhs.m_integrationtype);
        }

        bool operator<(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs)
        {            
            return (lhs.m_integrationtype < rhs.m_integrationtype);
        }
        
        bool TimeIntegrationSchemeKey::opLess::operator()(const TimeIntegrationSchemeKey &lhs, const TimeIntegrationSchemeKey &rhs) const
        {
            return (lhs.m_integrationtype < rhs.m_integrationtype);
        }

        std::ostream& operator<<(std::ostream& os, const TimeIntegrationSchemeKey& rhs)
        {
            os << "Time Integration Scheme: " << TimeIntegrationTypeMap[rhs.GetIntegrationSchemeType()] << endl;

            return os;
        }

        TimeIntegrationSolution::TimeIntegrationSolution(TimeIntegrationType schemeType, const DoubleArray& y, NekDouble t):
            m_schemeType(schemeType),
            m_y(1),
            m_t(1)
        {        
            m_y[0] = y;
            m_t[0] = t;    
        }

        TimeIntegrationSolution::TimeIntegrationSolution(TimeIntegrationType schemeType, const TripleArray& y, const Array<OneD, NekDouble>& t):
            m_schemeType(schemeType),
            m_y(y),
            m_t(t)
        {
        }

        TimeIntegrationSolution::TimeIntegrationSolution(TimeIntegrationType schemeType, 
                                                         unsigned int nsteps,
                                                         unsigned int nvar,
                                                         unsigned int npoints):
            m_schemeType(schemeType),
            m_y(nsteps),
            m_t(nsteps)
        {
            for(int i = 0; i < nsteps; i++)
            {
                m_y[i] = Array<OneD, Array<OneD, NekDouble> >(nvar);
                for(int j = 0; j < nvar; j++)
                {
                    m_y[i][j] = Array<OneD,NekDouble>(npoints);
                }
            }            
        }
        
        TimeIntegrationSchemeSharedPtr TimeIntegrationScheme::Create(const TimeIntegrationSchemeKey &key)
        {
            TimeIntegrationSchemeSharedPtr returnval(new TimeIntegrationScheme(key));
            return returnval;
        }

        TimeIntegrationScheme::TimeIntegrationScheme(const TimeIntegrationSchemeKey &key):
            m_schemeKey(key)
        {
            switch(key.GetIntegrationSchemeType())
            {
            case eForwardEuler:
            case eAdamsBashforthOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eAdamsBashforthOrder2:
                {
                    m_numsteps = 2;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,0.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,0.0);

                    m_B[0][0] = 3.0/2.0;
                    m_B[1][0] = 1.0;

                    m_U[0][0] = 1.0;

                    m_V[0][0] = 1.0;
                    m_V[0][1] = -0.5;
                }
                break;
            case eBackwardEuler:
            case eAdamsMoultonOrder1:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,1.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eAdamsMoultonOrder2:
                {
                    m_numsteps = 1;
                    m_numstages = 1;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.5);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,1.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);
                }
                break;
            case eMidpoint:
                {
                    m_numsteps = 1;
                    m_numstages = 2;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,0.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);

                    m_A[1][0] = 0.5;
                    m_B[0][1] = 1.0;
                }
                break;
            case eClassicalRungeKutta4:
                {
                    m_numsteps = 1;
                    m_numstages = 4;
                    m_A = Array<TwoD,NekDouble>(m_numstages,m_numstages,0.0);
                    m_B = Array<TwoD,NekDouble>(m_numsteps,m_numstages,0.0);
                    m_U = Array<TwoD,NekDouble>(m_numstages,m_numsteps,1.0);
                    m_V = Array<TwoD,NekDouble>(m_numsteps,m_numsteps,1.0);

                    m_A[1][0] = 0.5;
                    m_A[2][1] = 0.5;
                    m_A[3][2] = 1.0;

                    m_B[0][0] = 1.0/6.0;
                    m_B[0][1] = 1.0/3.0;
                    m_B[0][2] = 1.0/3.0;
                    m_B[0][3] = 1.0/6.0;
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,"Invalid Time Integration Scheme Type");
                }
            }
        }

        bool TimeIntegrationScheme::IsExplicitMethod() const
        {
            int i,j;
            
            for(i = 0; i < m_numstages; i++)
            {
                for(j = i; j < m_numstages; j++)
                {
                    if( fabs(m_A[i][j]) > NekConstants::kNekZeroTol )
                    {
                        return false;
                    }
                }
            }
            
            return true;
        }
        
        std::ostream& operator<<(std::ostream& os, const TimeIntegrationScheme& rhs)
        {
            int i,j;
            int r = rhs.GetNsteps();
            int s = rhs.GetNstages();

            int oswidth = 8;

            os << "Time Integration Scheme: " << TimeIntegrationTypeMap[rhs.GetIntegrationSchemeType()] << endl;
            os << "- number of steps:  " << r << endl;
            os << "- number of stages: " << s << endl;
            os << "General linear method tableau: " << endl;

            for(i = 0; i < s; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetA())[i][j] << " ";
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetU())[i][j] << " ";
                }
                os << endl;
            }
            for(int i = 0; i < (r+s)*oswidth+2; i++)
            {
                os << "-";
            }
            os << endl;
            for(i = 0; i < r; i++)
            {
                for(j = 0; j < s; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetB())[i][j] << " ";
                }
                os << " |"; 

                for(j = 0; j < r; j++)
                {
                    os.width(oswidth);
                    os << right << (rhs.GetV())[i][j] << " ";
                }
                os << endl;
            }
            return os;
        }


    }
}

