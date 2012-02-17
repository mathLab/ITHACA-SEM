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
// Description: Pulse Wave Propagation solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <PulseWaveSolver/EquationSystems/PulseWavePropagation.h>

namespace Nektar
{
    string PulseWavePropagation::className = GetEquationSystemFactory().RegisterCreatorFunction("PulseWavePropagation", PulseWavePropagation
																								::create, "Pulse Wave Propagation equation.");

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
	 * DoOdeRhs
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
			case MultiRegions::eDiscontinuousGalerkin:
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
				
				// Calculate second term in weak formulation (17) 
                for(i = 0; i < nvariables; ++i)
                {
					// Get the ith component of the flux vector F=(Q;p_t)
					PulseWavePropagation::GetFluxVector(i,physarray,fluxvector);
					
					// Calculate Derivative dF/dx
					m_fields[0]->PhysDeriv(0,fluxvector[0],outarray[i]);
                    
					// Negate as shifted on right hand side
					Vmath::Neg(nq,outarray[i],1);
                }
			}	
			break;
        }
    }


    /**
     *	DoOdeProjection
     */
    void PulseWavePropagation::DoOdeProjection(const Array<OneD,const Array<OneD, NekDouble> >&inarray,
											   Array<OneD, Array<OneD, NekDouble> >&outarray,
											   const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
		NekDouble Q, A_r, u_r, Au, uu;
        
		SetBoundaryConditions(time);
				
		// Loop over Boundary Regions to find the Q-inflow type
		for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
		{					
			if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eQinflow)
			{
				// Note: The Q value is contained in A in the inputfile, the value in u has to be 1.0 
				ASSERTL0((m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] == 1.0,
						 "For the Q-inflow BC the value in u must be 1.0");
				
				// Get the values of all variables needed for the Riemann problem
				Q = (m_fields[0]->UpdateBndCondExpansion(n))->GetCoeffs()[0];
				A_r = m_fields[0]->GetCoeffs()[0];
				u_r = m_fields[1]->GetCoeffs()[0];
				
				// Call the Q-inflow Riemann solver
				Q_inflowRiemannSolver(Q,A_r,u_r,m_A_0[0],m_beta[0],Au,uu);
				
				// Store the upwinded values in the boundary condition
				(m_fields[0]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = Au;
				(m_fields[1]->UpdateBndCondExpansion(n))->UpdateCoeffs()[0] = uu;
			}
		}
			
	
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
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
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
				
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i],coeffs,false);
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
	 * Calculates the second term of the weak form: dF/dx 
	 * The variables ot the system are (A;u)
	 * physfield[0] = A
	 * physfield[1] = u
	 * flux[0] = F[0] = A*u
	 * flux[1] = F[1] = u^2/2 + p/rho
	 * p-A-relationship: p = p_ext + beta*(sqrt(A)-sqrt(A_0))
	 */
    void PulseWavePropagation::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
											   Array<OneD, Array<OneD, NekDouble> > &flux)
    {
		int nq = m_fields[0]->GetTotPoints();
		NekDouble rho = m_rho; 
		NekDouble pext = m_pext; 
		NekDouble p = 0.0;
		NekDouble p_t = 0.0;
		NekDouble h0 = m_h0; 
		NekDouble nue = m_nue; 
		
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
	 * Calculates the third term of the weak form: numerical flux at boundary
	 *
	 */
    void PulseWavePropagation::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
											   Array<OneD, Array<OneD, NekDouble> > &numflux)
    {		
        int i;
		int nTraceNumPoints = GetTraceTotPoints();
		int nvariables      = 2; //(A,u)
		int nq = m_fields[0]->GetNpoints();		
		NekDouble rho = m_rho; 
		NekDouble pext = m_pext; 
		NekDouble h0 = m_h0; 
		NekDouble nue = m_nue; 
		
		Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
		Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
		for (i = 0; i < nvariables; ++i)
		{
			Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
			Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		}
		
		/*/ Print the values in physfield
		for (int i=0; i<physfield[0].num_elements(); i++)
		{
			cout << "physfield[A]["<<i<<"] = "<<physfield[0][i]<<"\t";
			cout << "physfield[u]["<<i<<"] = "<<physfield[1][i]<<endl;
		}*/
		
		// Get the physical values at the trace
		for (i = 0; i < nvariables; ++i)
		{
			m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
		}
				
		/*/ Print the values in Fwd and Bwd
		 for (int i=0; i<Fwd[0].num_elements(); i++)
		 {
		 cout << "Fwd[A]["<<i<<"] = "<<Fwd[0][i]<<"\t";
		 cout << "Bwd[A]["<<i<<"] = "<<Bwd[0][i]<<"\t\t";
		 
		 cout << "Fwd[u]["<<i<<"] = "<<Fwd[1][i]<<"\t";
		 cout << "Bwd[u]["<<i<<"] = "<<Bwd[1][i]<<endl;
		 }*/
		
		// Get A_0 at the trace
		Array<OneD, NekDouble> A_trace(GetTraceTotPoints());
		m_fields[0]->ExtractTracePhys(m_A_0,A_trace);
		A_trace[GetTraceTotPoints()-1] = m_A_0[GetTotPoints()-1];
		
		// Get the material properties at the trace
		Array<OneD, NekDouble> beta_trace(GetTraceTotPoints());
		m_fields[0]->ExtractTracePhys(m_beta,beta_trace);
		beta_trace[GetTraceTotPoints()-1] = m_beta[GetTotPoints()-1];

        // Solve the upwinding Riemann problem within one arterial segment
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
	 * Upwinding Riemann solver for pulse wave propagaiton
	 * 
	 */
	void PulseWavePropagation::RiemannSolverUpwind(NekDouble AL,NekDouble uL,NekDouble AR,NekDouble uR, 
												   NekDouble &Aflux, NekDouble &uflux, int i, NekDouble A_0, NekDouble beta)
	{
		int nvariables      = 2;
		int nq = m_fields[0]->GetNpoints();
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
		NekDouble h0 = m_h0; 
		NekDouble nue = m_nue; 
		
		// Compute the wave speeds
		cL = sqrt(beta*sqrt(AL)/(2*rho));
		cR = sqrt(beta*sqrt(AR)/(2*rho));
		c_Roe =(cL+cR)/2;
		//cout << "cL = "<<cL<<"\tcR = "<<cR<<"\tc_Roe = "<<c_Roe<<endl;
		
		u_Roe = (uL+uR)/2;
		//cout << "uL = "<<uL<<"\tuR = "<<uR<<"\tu_Roe = "<<u_Roe<<endl;
		
		lambda[0]= u_Roe + c_Roe;
		lambda[1]= u_Roe - c_Roe;
		//cout << "lambda[0] = "<<lambda[0]<<"\tlambda[1] = "<<lambda[1]<<endl;
		
		// Calculate the caracteristic variables 
		// Left characteristics
		characteristic[0] = uL + 4*sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		characteristic[1] = uL - 4*sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		// Right characteristics
		characteristic[2] = uR + 4*sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		characteristic[3] = uR - 4*sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		
		for (int k=0; k<4; k++)
		{
			//cout << "characteristic["<<k<<"] = "<<characteristic[k]<<endl;
		}
		
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
		
		for (int i=0; i<2; i++)
		{
			//cout << "upwinded W["<<i<<"] = "<<W[i]<<endl;
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
	 * Q-inflow Riemann solver for pulse wave propagation.
	 * This Riemann solver is called by SetBoundaryCondition_Pulse()
	 * in case of the inflow boundary condition is "Q_INFLOW" type.
	 * Returns the upwinded quantities and stores them in the bc's
	 */
	void PulseWavePropagation::Q_inflowRiemannSolver(NekDouble Q,NekDouble A_r,NekDouble u_r,NekDouble A_0, NekDouble beta,
													 NekDouble &Au,NekDouble &uu)
	{		
		NekDouble W2 = 0.0;
		NekDouble c = 0.0;
		NekDouble A_calc = 0.0;
		NekDouble fa = 0.0;
		NekDouble dfa = 0.0;
		NekDouble delta_A_calc = 0.0;
		NekDouble p = 0.0;
		NekDouble pext = 0.0;
		NekDouble p_t = 0.0;
		NekDouble rho = m_rho; 
	 
		int proceed = 1;
		int iter = 0;
		int MAX_ITER = 200;
	 
		// Tolerances for the algorithm
		NekDouble Tol = 1.0e-10;
	 
		// Riemann invariant W2(Ar,ur)
		W2 = u_r - 4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_r)) - sqrt(sqrt(A_0)));
	 
		// Calculate the wave speed
		c = sqrt(beta/(2*rho))*sqrt(sqrt(A_r));
		
		// Newton Iteration (Area only)
		A_calc = A_r;
		while ((proceed) && (iter < MAX_ITER))
		{	
			iter =iter+1;
	 
			fa = Q - W2*A_calc - A_calc*4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_calc)) - sqrt(sqrt(A_0)));
			dfa = -W2 - A_calc*4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_calc)) - sqrt(sqrt(A_0))) - sqrt(beta/(2*rho))*sqrt(sqrt(A_calc));
			delta_A_calc = fa/dfa;
			A_calc = A_calc - delta_A_calc;
	 
			if (sqrt(delta_A_calc*delta_A_calc) < Tol)
				proceed = 0;
		}
		
		// Obtain u from W2 and A_calc
		uu = W2 + 4*sqrt(beta/(2*rho))*(sqrt(sqrt(A_calc)) - sqrt(sqrt(A_0))); 
		Au = A_calc;
		
		cout << "-----------------------------------------------------"<<endl;
		cout << "| Q_inflow Riemann solver; number of iterations: "<<iter<<"  |"<<endl;
		cout << "| A_u = "<<Au<<"\tu_u = "<<uu<<"\tQ = "<<Au*uu<<"\t\t    |"<<endl;
		cout << "----------------------------------------------------"<< endl;
		
	 }
	
	
	/**
	 * Print summary routine
	 */
    void PulseWavePropagation::v_PrintSummary(std::ostream &out)
    {
        PulseWaveSystem::v_PrintSummary(out);
    }

}
