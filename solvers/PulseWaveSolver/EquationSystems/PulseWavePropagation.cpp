///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePropagation.cpp
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
// Description: Pulse Wave Propagation solve routines based on the weak formulation (1): 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <PulseWaveSolver/EquationSystems/PulseWavePropagation.h>

namespace Nektar
{
    string PulseWavePropagation::className = GetEquationSystemFactory().RegisterCreatorFunction("PulseWavePropagation", PulseWavePropagation
																								::create, "Pulse Wave Propagation equation.");
	/**
     *  @class PulseWavePropagation 
     *
	 *  Set up the routines based on the weak formulation from "Computational Modelling of 1D blood
	 *  flow with variable mechanical properties" by S. J. Sherwin et al. The weak formulation (1) 
	 *  reads:
	 *  \f$ \sum_{e=1}^{N_{el}} \left[ \left( \frac{\partial \mathbf{U}^{\delta} }{\partial t} , 
	 *    \mathbf{\psi}^{\delta} \right)_{\Omega_e} - \left( \frac{\partial \mathbf{F(\mathbf{U})}^{\delta} }
	 *    {\partial x}, \mathbf{\psi}^{\delta}  \right)_{\Omega_e} + \left[ \mathbf{\psi}^{\delta} 
	 *    \cdot \{ \mathbf{F}^u - \mathbf{F}(\mathbf{U}^{\delta}) \} \right]_{x_e^l}^{x_eû} \right] = 0 \f$
     */ 
    PulseWavePropagation::PulseWavePropagation(const LibUtilities::SessionReaderSharedPtr& pSession)
	: PulseWaveSystem(pSession)
    {
    }

    void PulseWavePropagation::v_InitObject()
    {
        PulseWaveSystem::v_InitObject();
		
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs        (&PulseWavePropagation::DoOdeRhs, this);
            m_ode.DefineProjection	  (&PulseWavePropagation::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit Pulse Wave Propagation not set up.");
        }
    }

    PulseWavePropagation::~PulseWavePropagation()
    {
    }
	
	
	/**
	 *  Computes the right hand side of (1). The RHS is everything except the term that contains
	 *  the time derivative \f$\frac{\partial \mathbf{U}}{\partial t}\f$. In case of a Discontinuous
	 *  Galerkin projection, the routine WeakDGAdvection will be called which then calls 
	 *  v_GetFluxVector and v_NumericalFlux implemented in the PulseWavePropagation class. The 
	 *  continuous Galerkin projection calls this two routines directly.
	 */
    void PulseWavePropagation::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
										Array<OneD,        Array<OneD, NekDouble> >&outarray,
										const NekDouble time)
    {
        int i;
		int ndim = m_spacedim;
        int nvariables = inarray.num_elements();
        int nq = GetNpoints();
		int ncoeffs = GetNcoeffs();

        switch (m_projectionType)
        {
			case MultiRegions::eDiscontinuous:
            {
				Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
				Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);
				
				for (i = 0; i < nvariables; ++i)
				{
					physarray[i] = Array<OneD, NekDouble>(nq);
					Vmath::Vcopy(nq,inarray[i],1,physarray[i],1);
					
					modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
				}
				
				WeakDGAdvection(physarray, modarray, true, true);
				
				for(i = 0; i < nvariables; ++i)
				{
					Vmath::Neg(ncoeffs,modarray[i],1);
				}	  
				
				for(i = 0; i < nvariables; ++i)
				{
					m_fields[i]->MultiplyByElmtInvMass(modarray[i],modarray[i]);
					m_fields[i]->BwdTrans(modarray[i],outarray[i]);
				}
			}
			break;
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
				Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
                for(i = 0; i < nvariables; ++i)
                {
                    physarray[i] = Array<OneD, NekDouble> (nq);
					Vmath::Vcopy(nq,inarray[i],1,physarray[i],1);
                }
				
				Array<OneD, Array<OneD, NekDouble> > fluxvector(ndim);
				for(i = 0; i < ndim; ++i)
                {
                    fluxvector[i] = Array<OneD, NekDouble> (nq);
                }
				
				// Calculate second term in weak formulation (1) 
                for(i = 0; i < nvariables; ++i)
                {
					// Get the ith component of the flux vector \f$\maththbf{F} = [Au,p_t]^T\f$
					PulseWavePropagation::GetFluxVector(i,physarray,fluxvector);
					
					// Calculate Derivative \f$ \partial \mathbf{F} / \partial x \f$
					m_fields[0]->PhysDeriv(0,fluxvector[0],outarray[i]);
                    
					// Negate as shifted on right hand side
					Vmath::Neg(nq,outarray[i],1);
                }
			}	
			break;
        }
    }


    /**
     *	Does the projection between ... space and the ... space. Also checks for Q-inflow boundary 
	 *  conditions at the inflow of the current arterial segment and applies the Q-inflow if specified
     */
    void PulseWavePropagation::DoOdeProjection(const Array<OneD,const Array<OneD, NekDouble> >&inarray,
											   Array<OneD, Array<OneD, NekDouble> >&outarray,
											   const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
	NekDouble Q, A_r, u_r;
	NekDouble A_u, u_u;
        NekDouble R_t, A_l, u_l, u_0, c_0, c_l;
	if (time == 0)
	  {pc = 0.0;
	  }
		
		SetBoundaryConditions(time);
				
		// Loop over Boundary Regions to find the Q-inflow type and the terminal resistance type
		for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{					
			if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eQinflow && m_omega == 0)
			{
				// Note: The Q value is contained in A in the inputfile, the value in u has to be 1.0 
				ASSERTL0((m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] == 1.0,
						 "For the Q-inflow BC the value in u must be 1.0");
				
				// Get the values of all variables needed for the Riemann problem
				Q = (m_fields[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];
				A_r = m_fields[0]->GetCoeffs()[0];
				u_r = m_fields[1]->GetCoeffs()[0];
				
				// Call the Q-inflow Riemann solver
				Q_inflowRiemannSolver(Q,A_r,u_r,m_A_0[0],m_beta[0],A_u,u_u);

				// Set the boundary conditions to  prescribe
				A_l=A_r;
				u_l=2*u_u-u_r;

				// Store the updated values in the boundary condition
				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = A_l;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = u_l;
			}
			/* Find the terminal resistance boundary condition and calculate the reflection. We assume 
			 * A_r = A_l and apply the reflection in u_r after paper "Computational Modelling of 1D 
			 * blood flow"*/
			else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTerminal)
			{
			  // Note: The R_t value is contained in A in the inputfile
				R_t = (m_fields[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];
				ASSERTL0((-1<=R_t && R_t<=1),
						 "R_t must be comprised between -1 and 1");
								
				// Get the left values A_l and u_l needed for Eq. 37
				A_l = m_fields[0]->GetCoeffs()[1];
				u_l = m_fields[1]->GetCoeffs()[1];
				
				// Get the values at initial state u_0, c_0
				u_0 = 0.0; //for all vessels start from initial condition 0
				c_0 = sqrt(m_betaglobal[m_omega][0]/(2*m_rho))*sqrt(sqrt(m_A_0global[m_omega][0]));				
				
				// Calculate the boundary values
				A_r = A_l;
				c_l = sqrt(m_beta[GetTotPoints()-1]/(2*m_rho))*sqrt(sqrt(A_l));
				u_r = (1-R_t)*((u_l-u_0) + 4*(c_l-c_0)) - u_l;
				
				// Store the new values in the boundary condition
				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = A_r;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = u_r;
			}

			/* Find the terminal R boundary condition and calculates the updated velocity and area as well as the updated boundary conditions */
			else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eRterminal)
			{
			  // cout << "R boundary condition found" << endl;
			        NekDouble RT = m_RT;
				NekDouble pout = m_pout;

				// Get the values of all variables needed for the Riemann problem
				A_l = m_fields[0]->GetCoeffs()[1];
				u_l = m_fields[1]->GetCoeffs()[1];

				// Goes through the resistance			

				// Call the R RiemannSolver
				R_RiemannSolver(RT,A_l,u_l,m_A_0[0],m_beta[GetTotPoints()-1],pout,A_u,u_u);

				// Calculates the new boundary conditions
				A_r=A_l;
				u_r=2*u_u-u_l;

				// Store the updated values in the boundary condition

				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = A_r;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = u_r;
			}

			/* Find the terminal CR boundary condition and calculates the updated velocity and area as well as the updated boundary conditions */
			else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eCRterminal)
			{
			        // cout << "CR boundary condition found" << endl;
			        NekDouble C = m_C;
			        NekDouble R = m_RT;
				NekDouble pout = m_pout;
				
				// Get the values of all variables needed for the Riemann problem
				A_l = m_fields[0]->GetCoeffs()[1];
				u_l = m_fields[1]->GetCoeffs()[1];
				
				// Call the CR Riemann solver
				CR_RiemannSolver(C,R,A_l,u_l,m_A_0[0],m_beta[GetTotPoints()-1],pout,A_u,u_u);

				// Calculates the boundary conditions
				u_r=u_l;
				A_r=4*sqrt(m_beta[GetTotPoints()-1]/(2*m_rho))*(2*sqrt(sqrt(A_u))-sqrt(sqrt(A_l)));
				A_r=A_r*A_r;
				A_r=A_r*A_r;
				
				// Store the updated values in the boundary condition

				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = A_r;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = u_r;
			}

			/* Find the terminal RCR boundary condition and calculates the updated velocity and area as well as the updated boundary conditions */
			else if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eRCRterminal)
			{

			        NekDouble C = m_C;
			        NekDouble RT = m_RT;
				NekDouble R1;
				NekDouble R2;
				NekDouble pout = m_pout;
				NekDouble rho = m_rho;

				// Get the values of all variables needed for the Riemann problem
				A_l = m_fields[0]->GetCoeffs()[1];
				u_l = m_fields[1]->GetCoeffs()[1];

				// Goes through the first resistance

				// Calculate c_0

				c_0 = sqrt(m_betaglobal[m_omega][0]/(2*m_rho))*sqrt(sqrt(m_A_0global[m_omega][0]));			

				// Calculate R1 and R2, R1 being calculated so as to eliminate reflections in the vessel
				R1 = rho*c_0/m_A_0[0];
				R2 = RT-R1;

				// Call the R RiemannSolver
				R_RiemannSolver(R1,A_l,u_l,m_A_0[0],m_beta[GetTotPoints()-1],pc,A_u,u_u);
				A_r = A_l;
				u_r = 2*u_u-u_l;

				// Goes through the CR system, it consists in updating the pressure pc

				pc = pc + m_timestep/C*(A_u*u_u-(pc-pout)/R2);

				// Store the updated values in the boundary condition

				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = A_r;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = u_r;
			}

			
		}
			
		// Do actual projection
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();
                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
            }
            break;
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
				
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i],coeffs);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
                }
            }
			break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

	
	/**
	 *  Calculates the second term of the weak form (1):
	 *  \f$ \left( \frac{\partial \mathbf{F(\mathbf{U})}^{\delta} }{\partial x}, \mathbf{\psi}^{\delta} \right)_{\Omega_e} \f$
	 *  The variables of the system are $\mathbf{U} = [A,u]^T$
	 *  physfield[0] = A        physfield[1] = u
	 *  flux[0] = F[0] = A*u    flux[1] = F[1] = u^2/2 + p/rho
	 *  p-A-relationship: p = p_ext + beta*(sqrt(A)-sqrt(A_0))
	 */
    void PulseWavePropagation::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
											   Array<OneD, Array<OneD, NekDouble> > &flux)
    {
		int nq = m_fields[0]->GetTotPoints();
		NekDouble rho = m_rho; 
		NekDouble pext = m_pext; 
		NekDouble p = 0.0;
		NekDouble p_t = 0.0;
		
        switch (i)
		{
			case 0:   // Flux for A equation
			{
				for (int j = 0; j < nq; j++)
				{
					flux[0][j] = physfield[0][j]*physfield[1][j]; 
				}
			}
			break;
			case 1:  // Flux for u equation
 			{
				for (int j = 0; j < nq; j++)
				{
					ASSERTL0(physfield[0][j]>=0,"Negative A not allowed.");
					p = pext + m_beta[j]*(sqrt(physfield[0][j]) - sqrt(m_A_0[j]));
					p_t = (physfield[1][j]*physfield[1][j])/2 + p/rho;
					flux[0][j] =  p_t;
				}
			}
			break;
			default:
				ASSERTL0(false,"GetFluxVector: illegal vector index");
			break;
		}
    }

	
	/**
	 *  Calculates the third term of the weak form (1): numerical flux at boundary
	 *  \f$ \left[ \mathbf{\psi}^{\delta} \cdot \{ \mathbf{F}^u - \mathbf{F}(\mathbf{U}^{\delta})
	 *  \} \right]_{x_e^l}^{x_eû} \f$
	 */
    void PulseWavePropagation::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
											   Array<OneD, Array<OneD, NekDouble> > &numflux)
    {		
        int i;
		int nTraceNumPoints = GetTraceTotPoints();
		int nvariables      = 2; //(A,u)
		
		Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
		Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
		for (i = 0; i < nvariables; ++i)
		{
			Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
			Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		}
		
		// Get the physical values at the trace
		for (i = 0; i < nvariables; ++i)
		{
			m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
		}
				
		// Get A_0 at the trace
		Array<OneD, NekDouble> A_trace(GetTraceTotPoints());
		m_fields[0]->ExtractTracePhys(m_A_0,A_trace);
		A_trace[GetTraceTotPoints()-1] = m_A_0[GetTotPoints()-1];
		
		// Get the material properties at the trace
		Array<OneD, NekDouble> beta_trace(GetTraceTotPoints());
		m_fields[0]->ExtractTracePhys(m_beta,beta_trace);
		beta_trace[GetTraceTotPoints()-1] = m_beta[GetTotPoints()-1];

        // Solve the upwinding Riemann problem within one arterial segment by calling the upwinding
		// Riemann solver implemented in this file
        NekDouble Aflux, uflux;
        for (i = 0; i < nTraceNumPoints; ++i)
		{
			switch(m_upwindTypePulse)
			{
				case eUpwindPulse:
				{
					RiemannSolverUpwind(Fwd[0][i],Fwd[1][i],Bwd[0][i],Bwd[1][i],
										Aflux, uflux,i, A_trace[i], beta_trace[i]);
				}
					break;
				default:
				{
					ASSERTL0(false,"populate switch statement for upwind flux");
				}
					break;
			}
			numflux[0][i]  = Aflux;
			numflux[1][i] = uflux;
		}
    }
	
	
	/**
	 *  Riemann solver for upwinding at an interface between two elements. Uses the characteristic variables 
	 *  for calculating the upwinded state \f$(A_u,u_u)\f$ from the left \f$(A_L,u_L)\f$ and right state \f$(A_R,u_R)\f$. 
	 *  Returns the upwinded flux $\mathbf{F}^u$ needed for the weak formulation (1). Details can be found in "Pulse 
	 *  wave propagation in the human vascular system", section 3.3 
	 */
	void PulseWavePropagation::RiemannSolverUpwind(NekDouble AL,NekDouble uL,NekDouble AR,NekDouble uR, 
												   NekDouble &Aflux, NekDouble &uflux, int i, NekDouble A_0, NekDouble beta)
	{
		int nvariables      = 2;
		Array<OneD, NekDouble> characteristic(4);
		Array<OneD, NekDouble> W(2);
		Array<OneD, NekDouble> lambda(nvariables);
		Array<OneD, NekDouble> upwindedphysfield(2);
		NekDouble cL = 0.0;
		NekDouble cR = 0.0;
		NekDouble c_Roe = 0.0;
		NekDouble u_Roe =0.0;
		NekDouble rho = m_rho; 
		NekDouble pext = m_pext; 
		NekDouble p = 0.0;
		NekDouble p_t = 0.0;
		
		// Compute the wave speeds
		cL = sqrt(beta*sqrt(AL)/(2*rho));
		cR = sqrt(beta*sqrt(AR)/(2*rho));
		c_Roe =(cL+cR)/2;		
		u_Roe = (uL+uR)/2;
		lambda[0]= u_Roe + c_Roe;
		lambda[1]= u_Roe - c_Roe;
		
		// Calculate the caracteristic variables 
		// Left characteristics \f$W_1^l, W_2^l\f$
		characteristic[0] = uL + 4*cL; //sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		characteristic[1] = uL - 4*cL; //sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		// Right characteristics \f$W_1^r, W_2^r\f$
		characteristic[2] = uR + 4*cR; //sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		characteristic[3] = uR - 4*cR; //sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		
		// Take left or right value of characteristic variable
		for (int j=0; j<2; j++)
		{
			if (lambda[j]>=0.0)
			{	 
				W[j]=characteristic[j];
			}
			if(lambda[j]<0.0)
			{
				W[j]=characteristic[j+2];
			}
		}

		// Calculate conservative variables from characteristics
		upwindedphysfield[0]= ((W[0]-W[1])/4)*((W[0]-W[1])/4)*((W[0]-W[1])/4)*((W[0]-W[1])/4)*(rho/(2*beta))*(rho/(2*beta));
		upwindedphysfield[1]= (W[0] + W[1])/2;
		
		// Compute the fluxes
		Aflux = upwindedphysfield[0] * upwindedphysfield[1];
		p = pext + beta*(sqrt(upwindedphysfield[0]) - sqrt(A_0));
		p_t = (upwindedphysfield[1]*upwindedphysfield[1])/2 + p/rho;				
		uflux =  p_t;
		
	}

	
	/**
	 *  Q-inflow Riemann solver for pulse wave propagation. This Riemann solver is called
	 *  by DoOdeProjection in case of a Q-inflow boundary condition. It is based on the 
	 *  conservation of mass and total pressure and on the characteristic information. For
	 *  further details see "Pulse wave propagation in the human vascular system", section 3.4.1
	 *  Returns the upwinded quantities \f$(A_u,u_u)\f$ and stores them into the boundary values
	 */
	void PulseWavePropagation::Q_inflowRiemannSolver(NekDouble Q,NekDouble A_r,NekDouble u_r,NekDouble A_0, NekDouble beta,
													 NekDouble &Au,NekDouble &uu)
	{		
		NekDouble W2 = 0.0;
		NekDouble A_calc = 0.0;
		NekDouble fa = 0.0;
		NekDouble dfa = 0.0;
		NekDouble delta_A_calc = 0.0;
		NekDouble rho = m_rho; 
	 
		int proceed = 1;
		int iter = 0;
		int MAX_ITER = 200;
	 
		// Tolerances for the algorithm
		NekDouble Tol = 1.0e-10;
	 
		// Riemann invariant \f$W_2(Ar,ur)\f$
		W2 = u_r - 4*sqrt(beta/(2*rho))*sqrt(sqrt(A_r));
	 
		// Newton Iteration (Area only)
		A_calc = A_r;
		while ((proceed) && (iter < MAX_ITER))
		{	
			iter =iter+1;
	 
			fa = Q - W2*A_calc - A_calc*4*sqrt(beta/(2*rho))*sqrt(sqrt(A_calc));
			dfa = -W2 - 5*sqrt(beta/(2*rho))*sqrt(sqrt(A_calc));
			delta_A_calc = fa/dfa;
			A_calc = A_calc - delta_A_calc;
	 
			if (sqrt(delta_A_calc*delta_A_calc) < Tol)
				proceed = 0;
		}
		
		// Obtain u_u and A_u
		uu = W2+4*sqrt(beta/(2*rho))*sqrt(sqrt(A_calc)); 
		Au = A_calc;
		
		//cout << "-----------------------------------------------------"<<endl;
		//cout << "| Q_inflow Riemann solver; number of iterations: "<<iter<<"  |"<<endl;
		//cout << "| A_u = "<<Au<<"\tu_u = "<<uu<<"\tQ = "<<Au*uu<<"\t\t    |"<<endl;
		//cout << "----------------------------------------------------"<< endl;
		
	 }

  void PulseWavePropagation::R_RiemannSolver(NekDouble R,NekDouble A_l,NekDouble u_l,NekDouble A_0, NekDouble beta, NekDouble pout,
													 NekDouble &A_u,NekDouble &u_u)
	{		
		NekDouble W1 = 0.0;
		NekDouble c_l = 0.0;
		NekDouble pext = m_pext;
		NekDouble A_calc = 0.0;
		NekDouble fa = 0.0;
		NekDouble dfa = 0.0;
		NekDouble delta_A_calc = 0.0;
		NekDouble rho = m_rho;

		int proceed = 1;
		int iter = 0;
		int MAX_ITER = 200;

		// Tolerances for the algorithm
		NekDouble Tol = 1.0e-10;

		// Calculate the wave speed
		c_l = sqrt(beta/(2*rho))*sqrt(sqrt(A_l));
		// cout << "c_l=" << c_l << endl;
	 
		// Riemann invariant \f$W_1(Al,ul)\f$
		W1 = u_l + 4*c_l;	 

		// Newton Iteration (Area only)
		A_calc = A_l;
		while ((proceed) && (iter < MAX_ITER))
		{	
			iter =iter+1;
	 
			fa = R*W1*A_calc-4*R*sqrt(beta/(2*rho))*A_calc*sqrt(sqrt(A_calc))-pext-beta*(sqrt(A_calc)-sqrt(A_0))+pout;
			dfa = R*W1-5*R*sqrt(beta/(2*rho))*sqrt(sqrt(A_calc))-beta/(2*sqrt(A_calc));
			delta_A_calc = fa/dfa;
			A_calc = A_calc - delta_A_calc;
	 
			if (sqrt(delta_A_calc*delta_A_calc) < Tol)
				proceed = 0;
		}

		// Obtain u_u and A_u
		//u_u = W1 - 4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_calc))); 
		u_u=(pext+beta*(sqrt(A_calc)-sqrt(A_0))-pout)/(R*A_calc);
		A_u = A_calc;

		//cout << "----------------------------------------------------"<< endl;
		//cout << "| A_u = "<<Au<<"\tu_u = "<<uu<<"\tQ = "<<Au*uu<<"\t\t    |"<<endl;
		//cout << "----------------------------------------------------"<< endl;
		
	 }

	/**
	 *  CR Riemann solver for pulse wave propagation. 
	 */
  void PulseWavePropagation::CR_RiemannSolver(NekDouble C,NekDouble R,NekDouble A_l,NekDouble u_l,NekDouble A_0, NekDouble beta, NekDouble pout,
													 NekDouble &A_u,NekDouble &u_u)
	{		
 	  //cout << "Entering CR_RiemannSolver" << endl;
		NekDouble pext = m_pext;
		NekDouble A_calc = 0.0;
		// to modify
		NekDouble delta_t = m_timestep;

		// First order finite difference scheme
		
		A_calc = sqrt(A_l)+delta_t/(C*beta)*(A_l*u_l+1/R*(pout-pext-beta*(sqrt(A_l)-sqrt(A_0))));
		A_u=A_calc*A_calc;

		// u_u is assumed to be equal to u_l
		u_u = u_l; 
		
		//cout << "----------------------------------------------------"<< endl;
		//cout << "| A_u = "<<Au<<"\tu_u = "<<uu<<"\tQ = "<<Au*uu<<"\t\t    |"<<endl;
		//cout << "----------------------------------------------------"<< endl;
		
	 }
	
	
	/**
	 *  Print summary routine, calls virtual routine reimplemented in UnsteadySystem
	 */
    void PulseWavePropagation::v_PrintSummary(std::ostream &out)
    {
        PulseWaveSystem::v_PrintSummary(out);
    }

}
