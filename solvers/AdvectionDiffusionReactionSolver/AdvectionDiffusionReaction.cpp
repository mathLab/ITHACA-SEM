///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.cpp
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
// Description: Advection Diffusion Reaction class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <AdvectionDiffusionReactionSolver/AdvectionDiffusionReaction.h>
#include <cstdio>
#include <cstdlib>
namespace Nektar
{
    /**
     * @class AdvectionDiffusionReaction
     * This is a solver class for solving the ADR class of problems.
     * - Laplace equation: \f$ \nabla^2 u=0 \f$
     * - Poisson equation: \f$ \nabla^2 u=f \f$
     * - Helmholtz equation: \f$ \nabla^2 u + k^2u = 0 \f$
     * - Steady Advection equation: 
     *   \f$ c \cdot \nabla u = 0\f$
     * - Steady Diffusion equation:
     *   \f$ -\nabla \cdot (k\nabla u) = f(x)\f$
     * - Steady Diffusion-Reaction equation:
     *   \f$ -\nabla \cdot (k\nabla u) = f(u,x)\f$
     * - Unsteady Advection equation:
     *   \f$ d_t u + c \cdot \nabla u = 0\f$
     * - Unsteady Inviscid Burger equation:
     *   \f$ d_t u + u \nabla u = 0\f$
     * - Unsteady Diffusion equation:
     *   \f$ d_t u -\nabla \cdot (k\nabla u) = f(x)\f$
     * - Unsteady Diffusion-Reaction equation:
     *   \f$ d_t u -\nabla \cdot (k\nabla u) = f(u,x)\f$
     */

    /**
     * Performs default initialisation of members.
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction():
            ADRBase(),
            m_infosteps(100),
            m_explicitAdvection(true),
            m_explicitDiffusion(true),
            m_explicitReaction(true)
    {
    }

    /**
     * Initialises members and loads a given session. Extracts the equation
     * type and then performs equation-specific initialisations.
     * @todo Missing equation types!
     * @param   fileNameString      Filename of session file.
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(
                                                string &fileNameString):
            ADRBase(fileNameString,true),
            m_infosteps(10),
            m_explicitDiffusion(true),
            m_explicitReaction(true)
    {
        
        int i;

        if(m_boundaryConditions->CheckForParameter("epsilon") == true)
        {
            m_epsilon = m_boundaryConditions->GetParameter("epsilon");
        }
        else
        {
            m_epsilon  = 0;
        }

        if(m_boundaryConditions->CheckForParameter("wavefreq") == true)
        {
            m_wavefreq = m_boundaryConditions->GetParameter("wavefreq");
        }
        else
        {
            m_wavefreq  = 0.0;
        }

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(NoCaseStringCompare(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i;
                break;
            }
        }

        ASSERTL0(i != (int) eEquationTypeSize, "Invalid expansion type.");

        // Velocity vectors
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
        for(int i = 0; i < m_spacedim; ++i)
        {
            m_velocity[i] = Array<OneD, NekDouble> (GetNpoints(),0.0);
        }

        // Equation specific Setups
        switch(m_equationType)
        {
        case eHelmholtz:
        case eLaplace:
        case eSteadyDiffusion:
        case ePoisson:
        case eSteadyDiffusionReaction:
            break;
        case eSteadyAdvection:

            EvaluateAdvectionVelocity();
            break;

        case eUnsteadyAdvection:

            EvaluateAdvectionVelocity();

            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;

            goto UnsteadySetup;
            break;
        case eUnsteadyInviscidBurger:
            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;
            goto UnsteadySetup;
            break;
        case eUnsteadyDiffusion:
            m_timeIntMethod = LibUtilities::eClassicalRungeKutta4;
            if(fabs(m_epsilon)<0.000000001)
            {
                m_epsilon = 1.0;
            }
            goto UnsteadySetup;
        case eUnsteadyAdvectionDiffusion:

            EvaluateAdvectionVelocity();

            m_timeIntMethod = LibUtilities::eIMEXdirk_3_4_3;

            goto UnsteadySetup;
            break;

        UnsteadySetup:
        {
            std::string Implicit = "Implicit"; 
            if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
            {
                m_infosteps =  m_boundaryConditions->GetParameter(
                                                            "IO_InfoSteps");
            }

            // check that any user defined boundary condition is indeed 
            // implemented
            for(int n = 0;
                    n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {
                // Time Dependent Boundary Condition (if no use defined then 
                // this is empty)
                if (m_fields[0]->GetBndConditions()[n]
                            ->GetUserDefined().GetEquation() != "")
                {
                    if (m_fields[0]->GetBndConditions()[n]
                            ->GetUserDefined().GetEquation() != "TimeDependent")
                    {
                        ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary "
                                       "condition");
                    }
                }
            }

            // Check for definition of Implicit/Explicit terms in solverinfo
            if(m_boundaryConditions->SolverInfoExists("ADVECTIONADVANCEMENT"))
            {
                std::string AdvStr = m_boundaryConditions
                                        ->GetSolverInfo("ADVECTIONADVANCEMENT");

                if(NoCaseStringCompare(AdvStr,Implicit) == 0)
                {
                    m_explicitAdvection = false;
                }
                else
                {
                    m_explicitAdvection = true;
                }
            }
            else
            {
                m_explicitAdvection = true;
            }

            if(m_boundaryConditions->SolverInfoExists("DIFFUSIONADVANCEMENT"))
            {
                std::string AdvStr = m_boundaryConditions
                                        ->GetSolverInfo("DIFFUSIONADVANCEMENT");

                if(NoCaseStringCompare(AdvStr,Implicit) == 0 )
                {
                    m_explicitDiffusion = false;
                    // Reset default for implicit diffusion
                    if(m_equationType == eUnsteadyDiffusion)
                    {
                        if(m_wavefreq>1000)
                        {
                            m_timeIntMethod = LibUtilities::eIMEXdirk_3_4_3;
                        }

                        else
                        {
                            m_timeIntMethod = LibUtilities::eDIRKOrder3;
                        }
                    }
                }
                else
                {
                    m_explicitDiffusion = true;
                }
            }
            else
            {
                m_explicitDiffusion = true;
            }

            if(m_boundaryConditions->SolverInfoExists("REACTIONADVANCEMENT"))
            {
                std::string AdvStr = m_boundaryConditions
                                        ->GetSolverInfo("REACTIONADVANCEMENT");

                if(NoCaseStringCompare(AdvStr,Implicit) == 0)
                {
                    m_explicitReaction = false;
                }
                else
                {
                    m_explicitReaction = true;
                }
            }
            else
            {
                m_explicitReaction = true;
            }

            // check to see if time stepping has been reset
            if(m_boundaryConditions->SolverInfoExists("TIMEINTEGRATIONMETHOD"))
            {
                std::string TimeIntStr = m_boundaryConditions
                                    ->GetSolverInfo("TIMEINTEGRATIONMETHOD");
                int i;
                for(i = 0;
                    i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
                {
                    if(NoCaseStringCompare(LibUtilities
                            ::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
                    {
                        m_timeIntMethod
                            = (LibUtilities::TimeIntegrationMethod)i;
                        break;
                    }
                }

                ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
                                            "Invalid time integration type.");
            }

            break;
        }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
    }


    /**
     * Populates the m_velocity fields using the user-defined equation in the
     * session file. This defines the velocity field for the SteadyAdvection
     * equation \f$ c\nabla u = 0\f$.
     */
    void AdvectionDiffusionReaction::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetNpoints();

        std::string velStr[3] = {"Vx","Vy","Vz"};

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_velocity.num_elements(); i++)
        {
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc
                        = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);

            for(int j = 0; j < nq; j++)
            {
                m_velocity[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
            }
        }
    }


    /**
     *
     */
    void AdvectionDiffusionReaction::ODErhs(const Array<OneD,
                        const Array<OneD, NekDouble> >&inarray,
                        Array<OneD,       Array<OneD, NekDouble> >&outarray,
                        const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

        Array<OneD, Array<OneD, NekDouble> > Forcing(1);

        switch(m_projectionType)
        {
            case eDiscontinuousGalerkin:
            {
                switch(m_equationType)
                {
                    case eUnsteadyAdvection:
		    case eUnsteadyAdvectionDiffusion:
                    case eUnsteadyInviscidBurger:
                    {

                        // cout << "time = " << time << endl;

                       SetBoundaryConditions(time);
                        WeakDGAdvection(inarray, outarray);
                        for(i = 0; i < nvariables; ++i)
                        {
                            m_fields[i]->MultiplyByElmtInvMass(outarray[i],
                                                                outarray[i]);
                            Vmath::Neg(ncoeffs,outarray[i],1);
                        }

                        if(m_wavefreq>0)
                        {
                            Forcing[0] = Array<OneD, NekDouble> (ncoeffs);
                            ODEeReaction(outarray,Forcing,time);
                            Vmath::Vadd(ncoeffs, Forcing[0], 1, outarray[0], 1, outarray[0], 1);
                        }

                        break;
                    }
                    case eUnsteadyDiffusion:
                    {
                        // BoundaryConditions are imposed weakly at the
                        // Diffusion operator

                        WeakDGDiffusion(inarray,outarray);

                        for(i = 0; i < nvariables; ++i)
                        {
                            m_fields[i]->MultiplyByElmtInvMass(outarray[i],
                                                                outarray[i]);
                        }

                        if(m_wavefreq>0)
                        {
                            Forcing[0] = Array<OneD, NekDouble> (ncoeffs);
                            ODEeReaction(outarray,Forcing,time);
                            Vmath::Vadd(ncoeffs, Forcing[0], 1, outarray[0], 1, outarray[0], 1);
                        }

                        break;
                    }
                }
                break;
            }
            case eGalerkin:
            {
                switch(m_equationType)
                {
                    case eUnsteadyAdvection:
		    case eUnsteadyAdvectionDiffusion:
                    {
                        SetBoundaryConditions(time);
                        Array<OneD, NekDouble> physfield(GetNpoints());

                        for(i = 0; i < nvariables; ++i)
                        {
                            m_fields[i]->MultiplyByInvMassMatrix(inarray[i],
                                                            outarray[i], false);
                            // Calculate -(\phi, V\cdot Grad(u))
                            m_fields[i]->BwdTrans_IterPerExp(outarray[i],
                                                                physfield);

                            WeakAdvectionNonConservativeForm(m_velocity,
                                                        physfield, outarray[i]);

                            Vmath::Neg(ncoeffs,outarray[i],1);
                        }
                        break;
                    }

                }
                break;
            }
            default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }


    /**
     *
     */
    void AdvectionDiffusionReaction::ODEeReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                  Array<OneD, Array<OneD, NekDouble> >&outarray,
                                                  const NekDouble time)
    {
        int i,k;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

	// PI*PI*exp(-1.0*PI*PI*FinTime)*sin(PI*x)*sin(PI*y)
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        Array<OneD, NekDouble> physfield(nq);

        NekDouble kt, kx, ky;
	for (i=0; i<nq; ++i)
	  {
              kt = m_wavefreq*time;
              kx = m_wavefreq*x0[i];
              ky = m_wavefreq*x1[i];

              // F(x,y,t) = du/dt + V \cdot \nabla u - \varepsilon \nabla^2 u 
              physfield[i] = (2.0*m_epsilon*m_wavefreq*m_wavefreq + m_wavefreq*cos(kt))*exp(sin(kt))*sin(kx)*sin(ky) 
                  + m_wavefreq*exp(sin(kt))*( m_velocity[0][i]*cos(kx)*sin(ky) + m_velocity[1][i]*sin(kx)*cos(ky) );
	  }
          m_fields[0]->FwdTrans(physfield, outarray[0]);
    }


    /**
     *
     */
    void AdvectionDiffusionReaction::ODElhs(const Array<OneD,
                        const Array<OneD, NekDouble> >&inarray,
                        Array<OneD,       Array<OneD, NekDouble> >&outarray,
                        const NekDouble time)
    {
        int nvariables = inarray.num_elements();
        MultiRegions::GlobalMatrixKey key(StdRegions::eMass);

        for(int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,inarray[i],
                                                                   outarray[i]);
        }
    }


    /**
     *
     */
    void AdvectionDiffusionReaction::ODElhsSolve(const Array<OneD,
                        const Array<OneD, NekDouble> >&inarray,
                        Array<OneD,       Array<OneD, NekDouble> >&outarray,
                        const NekDouble time)
    {
        SetBoundaryConditions(time);
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
            case eDiscontinuousGalerkin:

                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByElmtInvMass(inarray[i],outarray[i]);
                }
                break;
            case eGalerkin:
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(inarray[i],outarray[i],
                                                                    false);
                }
                break;
            }
            default:
                ASSERTL0(false,"Unknown projection scheme");
                break;
        }
    }


    /**
     *
     */
    void AdvectionDiffusionReaction::ODEhelmSolve(const Array<OneD,
                        const Array<OneD, NekDouble> >&inarray,
                        Array<OneD, Array<OneD, NekDouble> >&outarray,
                        NekDouble time,
                        NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        int nq = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > Forcing(1);
        Forcing[0] = Array<OneD, NekDouble>(ncoeffs,0.0);

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        MultiRegions::GlobalMatrixKey key(StdRegions::eMass);

        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply by inverse of mass matrix
            if(m_projectionType==eGalerkin)
            {
                //   m_fields[i]->MultiplyByInvMassMatrix(inarray[i],outarray[i],false);
            }

            // Multiply 1.0/timestep/lambda
            Vmath::Smul(ncoeffs, -1.0/lambda, inarray[i], 1, outarray[i], 1);

            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];

            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
            // m_fields[i]->SetPhysState(true);

            NekDouble kappa = 1.0/lambda/m_epsilon;

            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                            m_fields[i]->UpdateCoeffs(),kappa);
            m_fields[i]->SetPhysState(false);

            // The solution is Y[i]
             outarray[i] = m_fields[i]->GetCoeffs();

            // Multiply back by mass matrix
            if(m_projectionType==eGalerkin)
            {
                //   m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,
                //                                                  outarray[i],outarray[i]);
            }
        }
    }


    /**
     * Solve the Helmholtz problem
     * \f[ \nabla^2 \boldsymbol{u} + \lambda \boldsymbol{u} = \boldsymbol{f} \f]
     * for constant \f$\lambda\f$ and forcing term \f$\boldsymbol{f}\f$. This
     * is achieved by solving the separate Helmholtz problems for each
     * dependent variable. After solving, the result is in transformed space.
     * @param   lambda          Parameter.
     */
    void AdvectionDiffusionReaction::SolveHelmholtz(NekDouble lambda)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(),
                                   lambda);
            m_fields[i]->SetPhysState(false);
        }
    }

    // For Continuous Galerkin projections with time-dependent dirichlet boundary conditions,
    // the time integration can be done as follows:
    // The ODE resulting from the PDE can be formulated as:
    // 
    // M du/dt = F(u)  or du/dt = M^(-1) F(u)
    //
    // Now suppose that M does not depend of time, the ODE can than be written as:
    //
    // d(Mu)/dt = F(u)
    //
    // Introducing the variable u* = Mu, this yields
    //
    // du*/dt = F( M^(-1) u*  ) = F*(u*)
    //
    // So rather than solving the initial ODE, it is advised to solve this new ODE for u*
    // as this allows for an easier treatment of the dirichlet boundary conditions.
    // However, note that at the end of every time step, the actual solution u can
    // be calculated as:
    // 
    // u = M^(-1) u*;
    //
    // This can be viewed as projecting the solution u* onto the known boundary conditions.
    // Note that this step is also done inside the ODE rhs function F*.
    //
    // In order for all of this to work appropriately, make sure that the operator M^(-1)
    // does include the enforcment of the dirichlet boundary conditionst

  void AdvectionDiffusionReaction::GeneralTimeIntegration(int nsteps, 
                              LibUtilities::TimeIntegrationMethod IntMethod,
                              LibUtilities::TimeIntegrationSchemeOperators ode)
  {
    int i,n,nchk = 0;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nvariables = m_fields.num_elements();
 
    // Set up wrapper to fields data storage. 
    Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
    Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
    
    for(i = 0; i < nvariables; ++i)
      {
    m_fields[i]->SetPhysState(false);
    fields[i]  = m_fields[i]->UpdateCoeffs();
      }
  
    if(m_projectionType==eGalerkin)
      {
    // calculate the variable u* = Mu
    // we are going to TimeIntegrate this new variable u*
    MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
    for(int i = 0; i < nvariables; ++i)
      {
        tmp[i] = Array<OneD, NekDouble>(ncoeffs);
        //    m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
      }
      }
  
    // Declare an array of TimeIntegrationSchemes
    // For multi-stage methods, this array will have just one entry containing
    // the actual multi-stage method...
    // For multi-steps method, this can have multiple entries
    //  - the first scheme will used for the first timestep (this is an initialization scheme)
    //  - the second scheme will used for the first timestep (this is an initialization scheme)
    //  - ...
    //  - the last scheme will be used for all other time-steps (this will be the actual scheme)
    Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
    LibUtilities::TimeIntegrationSolutionSharedPtr u;
    int numMultiSteps;
    
        switch(IntMethod)
        {
        case LibUtilities::eIMEXdirk_2_3_2:
        case LibUtilities::eIMEXdirk_3_4_3:
        case LibUtilities::eDIRKOrder2:
        case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eForwardEuler:      
        case LibUtilities::eClassicalRungeKutta4:
            {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
            }
            break;
        case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);
                
                // Used for all other time steps 
                LibUtilities::TimeIntegrationSchemeKey IntKey1(IntMethod); 
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,ode);
            }
            break;
        default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }
                   
        for(n = 0; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
            }

            m_time += m_timestep;

            if(m_projectionType==eGalerkin)
            {
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
          m_fields[i]->SetPhysState(false);

          //      m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i],false);
          // fields[i] = tmp[i];                
                }
            }

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%m_infosteps))
            {
          cout << "Steps: " << n+1 << "\t Time: " << m_time << endl;
            }
            
            if(n&&(!((n+1)%m_checksteps)))
            {
          for(i = 0; i < nvariables; ++i)
        {
          (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
          Checkpoint_Output(nchk++);
            }
        }
        

        for(i = 0; i < nvariables; ++i)
        {
            (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }
    
    
  //----------------------------------------------------
  void AdvectionDiffusionReaction::SetBoundaryConditions(NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->EvaluateBoundaryConditions(time);
    }
  }
  
  // Evaulate flux = m_fields*ivel for i th component of Vu 
  // alt
  //          flux = 0.5(m_fields*m_fields) for the 
  void AdvectionDiffusionReaction::GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
                         Array<OneD, Array<OneD, NekDouble> > &flux)
  {

      switch(m_equationType)
      {
      case eUnsteadyAdvection:
      case eUnsteadyDiffusion:
      case eUnsteadyAdvectionDiffusion:
    {
      ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");
      
      for(int j = 0; j < flux.num_elements(); ++j)
        {
          Vmath::Vmul(GetNpoints(),physfield[i],1,
              m_velocity[j],1,flux[j],1);
        }
    }
    break;
      case eUnsteadyInviscidBurger:
    {
      for(int j = 0; j < flux.num_elements(); ++j)
        {
          Vmath::Vmul(GetNpoints(),physfield[i],1,
              physfield[i],1,flux[j],1);
          Vmath::Smul(GetNpoints(),0.5,flux[j],1,flux[j],1);
        }
    }
    break;
      default:
    ASSERTL0(false,"unknown equationType");
      }
  }
    
 // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
  void AdvectionDiffusionReaction::GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield, 
                         Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    switch(m_equationType)
      {
      case eUnsteadyAdvection:
      case eUnsteadyDiffusion:
      case eUnsteadyAdvectionDiffusion:
    {
      ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");
      
      for(int k = 0; k < flux.num_elements(); ++k)
        {
          Vmath::Zero(GetNpoints(),flux[k],1);
        }
      Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }
    break;
      case eUnsteadyInviscidBurger:
    {
      ASSERTL0(false,"should never arrive here ...");
    }
    break;
      default:
    ASSERTL0(false,"unknown equationType");
      }
  }
  
  
    void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                           Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

    switch(m_equationType)
      {
      case eUnsteadyAdvection:
      case eUnsteadyAdvectionDiffusion:
        {

          // Get Edge Velocity - Could be stored if time independent
          for(i = 0; i < nvel; ++i)
        {
          m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
          Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }
          
          for(i = 0; i < numflux.num_elements(); ++i)
        {
          m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
          //evaulate upwinded m_fields[i]
          m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
          // calculate m_fields[i]*Vn
          Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
        }
        break;
      case eUnsteadyInviscidBurger:
        {
          // Get Edge Velocity - Could be stored if time independent
          m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
          for(i = 0; i < nvel; ++i)
        {
          // m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
          Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }
          
          for(i = 0; i < numflux.num_elements(); ++i)
        {
          m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
          //evaulate upwinded m_fields[i]
          m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
          // calculate m_fields[i]*Vn
          Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
          Vmath::Smul(nTraceNumPoints,0.5,numflux[i],1,numflux[i],1);
        }
        }
        break;
      default:
        ASSERTL0(false,"unknown equationType");
      }
    }

  void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                         Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
                         Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
    Array<OneD, NekDouble > tmp(nTraceNumPoints,0.0);

    Array<OneD, Array<OneD, NekDouble > > traceVelocity(2);

    traceVelocity[0] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);
    traceVelocity[1] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);

    switch(m_equationType)
      {
      case eUnsteadyAdvection:
      case eUnsteadyAdvectionDiffusion:
        {

          // Get Edge Velocity - Could be stored if time independent
          m_fields[0]->ExtractTracePhys(m_velocity[0], traceVelocity[0]);
          m_fields[0]->ExtractTracePhys(m_velocity[1], traceVelocity[1]);
          
          m_fields[0]->GetFwdBwdTracePhys(physfield[0],Fwd,Bwd);
          
          m_fields[0]->GetTrace()->Upwind(traceVelocity,Fwd,Bwd,tmp);
          
          Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[0],1,numfluxX[0],1);
          Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[1],1,numfluxY[0],1);
        }
        break;
      default:
        ASSERTL0(false,"unknown equationType");
      }
    }


  // Compute the fluxes of q from the scalar functin u.
  // Input:   ufield (1 by nTraceNumPoints) - Should be in physical field
  // Output:  ufluxFwd  (2 by nTraceNumPoints) - Flux values for forward edges
  //          ufluxBwd  (2 by nTraceNumPoints) - Flux values for backward edges

    void AdvectionDiffusionReaction::NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield, 
                              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
    {
        int i,j;
        int nTraceNumPoints = GetTraceNpoints();
    int nvariables = m_fields.num_elements();
        int nqvar = uflux.num_elements();
    
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
    Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
    Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);
                      
    // Get the sign of (v \cdot n), v = an arbitrary vector

    // Vn = V \cdot n, where n is tracenormal for eForward edges. Set V = (1,0)
    // Vmath::Vcopy(nTraceNumPoints,m_traceNormals_tbasis[0],1,Vn,1);

        //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
        for (j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables ; ++i)
            {
                //  Compute Forward and Backward value of ufield of i direction

                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd
                
                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],Fwd,Bwd,fluxtemp);  
    
                // Imposing weak boundary condition with flux
                // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd
                
                // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

	        if(m_fields[0]->GetBndCondExpansions().num_elements())
		  {
                    WeakPenaltyforScalar(i,ufield[i],fluxtemp);
		  }

                // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )
                
                // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

        Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[j][i],1);

            }
    }
    }

    // Compute the fluxes of q and u vector fields for discontinuous diffusion term
    // Input:   qfield : 2 by # of total trace points
    // Output:  qflux  : 2 by # of total trace points
    void AdvectionDiffusionReaction::NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                         Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                                         Array<OneD, Array<OneD, NekDouble> >  &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
    int nvariables = m_fields.num_elements();
        int nqvar = qfield.num_elements();

    NekDouble C11 = 1.0;            
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
    Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
            
    Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
    Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

    Array<OneD, NekDouble > uterm(nTraceNumPoints);
              
        // Get the sign of (v \cdot n), v = an arbitrary vector
    // Vn = V \cdot n, where n is tracenormal for eForward edges
    // Vmath::Vcopy(nTraceNumPoints,m_traceNormals[0],1,Vn,1);
        
        // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)                 
        for(int i = 0; i < nvariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
            for(int j = 0; j < nqvar; ++j)
            {
        //  Compute Forward and Backward value of ufield of jth direction
        m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);
                
                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-
                
                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd =q+

        m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],qBwd,qFwd,qfluxtemp);
        Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);

        // Generate Stability term = - C11 ( u- - u+ )
        m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
        Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);                     
        Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

        //  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
        Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);
        
        // Imposing weak boundary condition with flux
        if(m_fields[0]->GetBndCondExpansions().num_elements())
          {
            WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
          }
        
        // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
        // n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
        // n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
        }
        }
    }
    

    
  // Diffusion: Imposing weak boundary condition for u with flux 
  //  uflux = g_D  on Dirichlet boundary condition
  //  uflux = u_Fwd  on Neumann boundary condition
  void AdvectionDiffusionReaction::WeakPenaltyforScalar(const int var,
							const Array<OneD, const NekDouble> &physfield, 
							Array<OneD, NekDouble> &penaltyflux,
							NekDouble time)
  {
    int i, j, e, npoints, id1, id2;
    // Number of boundary regions
    int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
    int Nfps, numBDEdge;
    int nTraceNumPoints = GetTraceNpoints();
    int cnt = 0;

    Array<OneD, NekDouble > uplus(nTraceNumPoints);
    
    m_fields[var]->ExtractTracePhys(physfield,uplus);   
    for(i = 0; i < nbnd; ++i)
      {     
	// Number of boundary expansion related to that region
	numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
	// Evaluate boundary values g_D or g_N from input files
	SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(0);
	npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
	Array<OneD,NekDouble> BDphysics(npoints);
	Array<OneD,NekDouble> x0(npoints,0.0);
	Array<OneD,NekDouble> x1(npoints,0.0);
	Array<OneD,NekDouble> x2(npoints,0.0);  
        
	m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
	for(j = 0; j < npoints; j++)
	  {
	    BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
	  }
	
	// Weakly impose boundary conditions by modifying flux values
	for (e = 0; e < numBDEdge ; ++e)
	  {
	    // Number of points on the expansion
	    Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0) ;
	    id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
	    id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

	    // For Dirichlet boundary condition: uflux = g_D
	    if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
	      {
		Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
	      }
	    
	    // For Neumann boundary condition: uflux = u+
	    else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
	      {
		Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
	      }
	  }
      }
  }
  
  // Diffusion: Imposing weak boundary condition for q with flux 
  //  uflux = g_D  on Dirichlet boundary condition
  //  uflux = u_Fwd  on Neumann boundary condition
  void AdvectionDiffusionReaction::WeakPenaltyforVector(const int var,
							const int dir,
							const Array<OneD, const NekDouble> &physfield,
							Array<OneD, NekDouble> &penaltyflux,
							NekDouble C11,
							NekDouble time)
  {
    int i, j, e, npoints, id1, id2;
    int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
    int numBDEdge, Nfps;
    int nTraceNumPoints = GetTraceNpoints();
    Array<OneD, NekDouble > uterm(nTraceNumPoints);
    Array<OneD, NekDouble > qtemp(nTraceNumPoints);
    int cnt = 0;

    m_fields[var]->ExtractTracePhys(physfield,qtemp);            
    
    for(i = 0; i < nbnd; ++i)
      {    
        numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();     
        // Evaluate boundary values g_D or g_N from input files
	SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(0);
	npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
	
	Array<OneD,NekDouble> BDphysics(npoints);
	Array<OneD,NekDouble> x0(npoints,0.0);
	Array<OneD,NekDouble> x1(npoints,0.0);
	Array<OneD,NekDouble> x2(npoints,0.0);  
        
	m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
	for(j = 0; j < npoints; j++)
	  {
	    BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
	  }
	
	// Weakly impose boundary conditions by modifying flux values
	for (e = 0; e < numBDEdge ; ++e)
	  {
	    Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);
     
	    id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
	    id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

	    // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ - g_D) (nx, ny)
	    if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
	      {
		Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&qtemp[id2],1,&penaltyflux[id2],1);
		
		// Vmath::Vsub(Nfps,&Fwd[id2],1,&BDphysics[id1],1,&uterm[id2],1);
		//    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&uterm[id2],1,&uterm[id2],1);
		// Vmath::Svtvp(Nfps,-1.0*C11,&uterm[id2],1,&qFwd[id2],1,&penaltyflux[id2],1);
	      }
	    
	    // For Neumann boundary condition: qflux = g_N
	    else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
	      {
		Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
	      }
	  }
      }       
  }
  
  
  void AdvectionDiffusionReaction::Summary(std::ostream &out)
  {   
    cout << "=======================================================================" << endl;
    cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
    ADRBase::SessionSummary(out);
    switch(m_equationType)
      {
      case eSteadyDiffusion: case eSteadyDiffusionReaction:
      case eHelmholtz: case eLaplace: case ePoisson:
    out << "\tLambda          : " << m_boundaryConditions->GetParameter("Lambda") << endl;
    for(int i = 0; i < m_fields.num_elements(); ++i)
          {
        out << "\tForcing (field " << i << ") : " << m_boundaryConditions->GetForcingFunction(i)->GetEquation() << endl;
          }
    
    break;
      case eUnsteadyAdvection: 
    if(m_explicitAdvection)
          {
        out << "\tAdvection Advancement   : Explicit" <<endl;
          }
    else
          {
        out << "\tAdvection Advancement   : Implicit" <<endl;
          }
    out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
    ADRBase::TimeParamSummary(out);
    break;
      case eUnsteadyDiffusion:
          out << "\tepsilon          : " << m_epsilon << endl;
    if(m_explicitDiffusion)
          {
        out << "\tDiffusion Advancement   : Explicit" <<endl;
          }
    else
          {
              out << "\tDiffusion Advancement   : Implicit" <<endl;
          }
          out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
          ADRBase::TimeParamSummary(out);
          break;
      case eUnsteadyAdvectionDiffusion:
          out << "\tepsilon          : " << m_epsilon << endl;
          if(m_explicitDiffusion)
          {
            out << "\tDiffusion Advancement   : Explicit" <<endl;
          }
          else
          {
              out << "\tDiffusion Advancement   : Implicit" <<endl;
          }
          if(m_explicitReaction)
          {
              out << "\tReaction Advancement    : Explicit" <<endl;
          }
          else
          {
              out << "\tReaction Advancement    : Implicit" <<endl;
          }
          out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
          ADRBase::TimeParamSummary(out);
          break;
      }
      cout << "=======================================================================" << endl;

    }
    
} //end of namespace

/**
* $Log: AdvectionDiffusionReaction.cpp,v $
* Revision 1.25  2010/03/11 15:31:29  sehunchun
* Changes for Neumann boundary test problems
*
* Revision 1.24  2010/03/10 23:40:25  sehunchun
* Add UnsteadyAdvectionDiffusion and new variables, m_epsilon and m_wavefreq
*
* Revision 1.23  2009/11/20 11:11:58  cbiotto
* Fixing the Neumann boundary condition.
*
* Revision 1.22  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.21  2009/07/23 05:32:28  sehunchun
* Implicit and Explicit diffusion debugging
*
* Revision 1.20  2009/07/09 21:24:57  sehunchun
* Upwind function is update
*
* Revision 1.19  2009/07/01 21:51:58  sehunchun
* Modification of ExDiffusion Solver and tabbing
*
* Revision 1.18  2009/06/11 02:20:43  claes
* *** empty log message ***
*
* Revision 1.17  2009/06/11 01:54:08  claes
* Added Inviscid Burger
*
* Revision 1.16  2009/04/29 20:45:09  sherwin
* Update for new eNum definition of EQTYPE
*
* Revision 1.15  2009/04/27 21:37:14  sherwin
* Updated to dump .fld and .chk file in compressed coefficient format
*
* Revision 1.14  2009/03/06 12:00:10  sehunchun
* Some minor changes on nomenclatures and tabbing errors
*
* Revision 1.13  2009/03/05 14:02:38  pvos
* Fixed bug
*
* Revision 1.12  2009/03/05 11:50:32  sehunchun
* Implicit scheme and IMEX scheme are now implemented
*
* Revision 1.11  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.10  2009/03/03 16:11:26  pvos
* New version of TimeIntegrator classes
*
* Revision 1.9  2009/02/28 21:59:09  sehunchun
* Explicit Diffusion solver is added
*
* Revision 1.8  2009/02/16 16:07:03  pvos
* Update of TimeIntegration classes
*
* Revision 1.7  2009/02/10 16:39:35  sherwin
* Added new format of SolverInfo reader to identify EQTYPE
*
* Revision 1.6  2009/02/08 09:13:08  sherwin
* Updates to go with Multiple matrix/variable solve
*
* Revision 1.6  2009/01/06 21:10:34  sherwin
* Updates for virtual calls to IProductWRTBase and introduced reader to handle SOLVERINFO section to specify different solvers
*
* Revision 1.5  2008/11/19 10:53:51  pvos
* Made 2D CG version working
*
* Revision 1.4  2008/11/17 08:20:14  claes
* Temporary fix for CG schemes. 1D CG working (but not for userdefined BC). 1D DG not working
*
* Revision 1.3  2008/11/12 12:12:26  pvos
* Time Integration update
*
* Revision 1.2  2008/11/02 22:38:51  sherwin
* Updated parameter naming convention
*
* Revision 1.1  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
