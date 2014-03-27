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
    string PulseWavePropagation::className = GetEquationSystemFactory().RegisterCreatorFunction("PulseWavePropagation", PulseWavePropagation::create, "Pulse Wave Propagation equation.");
    /**
     *  @class PulseWavePropagation 
     *
     *  Set up the routines based on the weak formulation from
     *  "Computational Modelling of 1D blood flow with variable
     *  mechanical properties" by S. J. Sherwin et al. The weak
     *  formulation (1) reads:
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
	
        m_pressureArea=GetPressureAreaFactory().CreateInstance("Lymphatic",m_vessels,m_session);
        m_pressureArea->DoPressure();
	
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs       (&PulseWavePropagation::DoOdeRhs, this);
            m_ode.DefineProjection   (&PulseWavePropagation::DoOdeProjection, this);
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
     *  v_GetFluxVector and v_NumericalFlux implemented in the PulseWavePropagation class. 
     *
     */
    void PulseWavePropagation::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                        Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                       const NekDouble time)
    {
        int i;
            
        Array<OneD, Array<OneD, NekDouble> > physarray(m_nVariables);
        Array<OneD, Array<OneD, NekDouble> > modarray (m_nVariables);
	
	Array<OneD, NekDouble> tmpArray;

        int cnt = 0;


        // Set up Inflow and Outflow boundary conditions. 
        SetPulseWaveBoundaryConditions(inarray, outarray, time);
    
        // Set up any interface conditions and write into boundary condition
        EnforceInterfaceConditions(inarray);
        
        // do advection evauation in all domains
        for(int omega=0; omega < m_nDomains; ++omega)
        {
            m_currentDomain = omega;
            int nq = m_vessels[omega*m_nVariables]->GetTotPoints();
            int ncoeffs = m_vessels[omega*m_nVariables]->GetNcoeffs();
            
            for (i = 0; i < m_nVariables; ++i)
            {
                physarray[i] = inarray[i]+cnt;
                modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
            }

            for(i = 0; i < m_nVariables; ++i)
            {
                m_fields[i] = m_vessels[omega*m_nVariables+ i];
            }
            
            WeakDGAdvection(physarray, modarray, true, true);
            
            for(i = 0; i < m_nVariables; ++i)
            {
                Vmath::Neg(ncoeffs,modarray[i],1);
            }	  
            
            for(i = 0; i < m_nVariables; ++i)
            {
                m_vessels[omega*m_nVariables+i]->MultiplyByElmtInvMass(modarray[i],modarray[i]);
                m_vessels[omega*m_nVariables+i]->BwdTrans(modarray[i],tmpArray = outarray[i]+cnt);
            }
            cnt += nq;
        }
    }

    void PulseWavePropagation::DoOdeProjection(const Array<OneD,const Array<OneD, NekDouble> >&inarray,
                                               Array<OneD, Array<OneD, NekDouble> >&outarray,
                                               const NekDouble time)
    {
        // Just copy over array
        for(int i = 0; i < m_nVariables; ++i)
        {
            Vmath::Vcopy(inarray[i].num_elements(),inarray[i],1,outarray[i],1);
        }
    }


	
    /**
     *	Does the projection between ... space and the ... space. Also checks for Q-inflow boundary 
     *  conditions at the inflow of the current arterial segment and applies the Q-inflow if specified
     */
    void PulseWavePropagation::SetPulseWaveBoundaryConditions(
        const Array<OneD,const Array<OneD, NekDouble> >&inarray, 
        Array<OneD, Array<OneD, NekDouble> >&outarray, 
        const NekDouble time)
        
    {
        int omega;
	NekDouble Q, A_r, u_r;
	NekDouble A_u, u_u;
        NekDouble R_t, A_l, u_l, u_0, c_0, c_l;
        
        Array<OneD, MultiRegions::ExpListSharedPtr>     vessel(2);

        int offset=0; 

        //This will be moved to the RCR boundary condition once factory is setup
	if (time == 0)
        {
            m_Boundary = Array<OneD,PulseWaveBoundarySharedPtr>(2*m_nDomains);

            for(omega = 0; omega < m_nDomains; ++omega)
            {
                vessel[0] = m_vessels[2*omega];

                for(int j = 0; j < 2; ++j)
                {	
                    std::string BCType = vessel[0]->GetBndConditions()[j]->
                        GetBndTypeAsString(vessel[0]->GetBndConditions()[j]->GetUserDefined());
                    m_Boundary[2*omega+j]=GetBoundaryFactory().CreateInstance(BCType,m_vessels,m_session,m_pressureArea);
                }
            }

        }

        SetBoundaryConditions(time);

        // Loop over all vessesls and set boundary conditions
        for(omega = 0; omega < m_nDomains; ++omega)
        {
            for(int n = 0; n < 2; ++n)
            {	
                m_Boundary[2*omega+n]->DoBoundary(inarray,m_A_0,m_beta,time,omega,offset,n);
            }
            offset += m_vessels[2*omega]->GetTotPoints();
        }

    }

	
	/**
	 *  Calculates the second term of the weak form (1): \f$
	 *  \left( \frac{\partial \mathbf{F(\mathbf{U})}^{\delta}
	 *  }{\partial x}, \mathbf{\psi}^{\delta} \right)_{\Omega_e}
	 *  \f$
	 *  The variables of the system are $\mathbf{U} = [A,u]^T$
	 *  physfield[0] = A        physfield[1] = u
	 *  flux[0] = F[0] = A*u    flux[1] = F[1] = u^2/2 + p/rho
	 *  p-A-relationship: p = p_ext + beta*(sqrt(A)-sqrt(A_0))
	 */
    void PulseWavePropagation::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                                               Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        int nq = m_vessels[m_currentDomain*m_nVariables]->GetTotPoints();
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

                    p = m_pext + m_beta[m_currentDomain][j]*
                        (sqrt(physfield[0][j]) - sqrt(m_A_0[m_currentDomain][j]));

                    p_t = (physfield[1][j]*physfield[1][j])/2 + p/m_rho;
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
        int nTracePts = GetTraceTotPoints();
        
        Array<OneD, Array<OneD, NekDouble> > Fwd(m_nVariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd(m_nVariables);
        
        for (i = 0; i < m_nVariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            Bwd[i] = Array<OneD, NekDouble>(nTracePts);
        }
	
        // Get the physical values at the trace
        for (i = 0; i < m_nVariables; ++i)
        {
            m_vessels[m_currentDomain*m_nVariables+ i]->
                GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
        }
        
        // Solve the upwinding Riemann problem within one arterial
        // segment by calling the upwinding Riemann solver implemented
        // in this file
        NekDouble Aflux, uflux;
        for (i = 0; i < nTracePts; ++i)
        {
            switch(m_upwindTypePulse)
            {
            case eUpwindPulse:
                {
                    RiemannSolverUpwind(Fwd[0][i],Fwd[1][i],Bwd[0][i],Bwd[1][i],
                                        Aflux, uflux, m_A_0_trace[m_currentDomain][i],
                                        m_beta_trace[m_currentDomain][i],
                                        m_trace_fwd_normal[m_currentDomain][i]);
                }
                break;
            default:
                {
                    ASSERTL0(false,"populate switch statement for upwind flux");
                }
                break;
            }
            numflux[0][i] = Aflux;
            numflux[1][i] = uflux;
        }
    }
    
    
    /**
     *  Riemann solver for upwinding at an interface between two
     *  elements. Uses the characteristic variables for calculating
     *  the upwinded state \f$(A_u,u_u)\f$ from the left
     *  \f$(A_L,u_L)\f$ and right state \f$(A_R,u_R)\f$.  Returns the
     *  upwinded flux $\mathbf{F}^u$ needed for the weak formulation
     *  (1). Details can be found in "Pulse wave propagation in the
     *  human vascular system", section 3.3
     *
     */
    void PulseWavePropagation::RiemannSolverUpwind(NekDouble AL,NekDouble uL,
                                                   NekDouble AR,NekDouble uR, 
                                                   NekDouble &Aflux, 
                                                   NekDouble &uflux, 
                                                   NekDouble A_0, 
                                                   NekDouble beta,
                                                   NekDouble n)
    {
        Array<OneD, NekDouble> W(2);
        Array<OneD, NekDouble> upwindedphysfield(2);
        NekDouble cL = 0.0;
        NekDouble cR = 0.0;
        NekDouble rho = m_rho; 
        NekDouble pext = m_pext; 
        NekDouble p = 0.0;
        NekDouble p_t = 0.0;
        
#if 0
        // Compute the wave speeds in the normal direction according
        // to the definition of Fwd and Bwd and indicated by n 
        cL = sqrt(beta*sqrt(AL)/(2*rho));
        cR = sqrt(beta*sqrt(AR)/(2*rho));

        NekDouble c_Roe = 0.0;
        NekDouble u_Roe =0.0;
        Array<OneD, NekDouble> lambda(2);
        Array<OneD, NekDouble> characteristic(4);
        
        c_Roe = (cL+cR)/2;		
        u_Roe = (uL+uR)/2;
        lambda[0]= u_Roe + c_Roe*n;
        lambda[1]= u_Roe - c_Roe*n;
        
        // Calculate the caracteristic variables 
        // Left characteristics \f$W_1^l, W_2^l\f$
        characteristic[0] = uL + 4*cL;
        characteristic[1] = uL - 4*cL;
        // Right characteristics \f$W_1^r, W_2^r\f$
        characteristic[2] = uR + 4*cR; 
        characteristic[3] = uR - 4*cR; 
        
        // Take left or right value of characteristic variable
        for (int j=0; j<2; j++)
        {
            if (lambda[j]>=0.0)
            {	 
                W[j]=characteristic[j];
            }
            else 
            {
                W[j]=characteristic[j+2];
            }
        }
#else
        // Compute the wave speeds in the normal direction according
        // to the definition of Fwd and Bwd and indicated by n 
        cL = sqrt(beta*sqrt(AL)/(2*rho))*n;
        cR = sqrt(beta*sqrt(AR)/(2*rho))*n;

        ASSERTL1(fabs(cL+cR) > fabs(uL+uR),"Conditions are not sub-sonic");

        // If upwinding from left and right for subsonic domain
        // then know characteristics immediately
        W[0] = uL + 4*cL;
        W[1] = uR - 4*cR;
#endif

        // Calculate conservative variables from characteristics
        NekDouble w0mw1 = 0.25*(W[0]-W[1]);
        NekDouble fac = rho/(2*beta);
        w0mw1 *= w0mw1; // squared
        w0mw1 *= w0mw1; // fourth power
        fac *= fac;     // squared
        upwindedphysfield[0]= w0mw1*fac;
        upwindedphysfield[1]= 0.5*(W[0] + W[1]);

        //fprintf(stdout," Upwind A: %16.14lf\n",upwindedphysfield[0]);
        //fprintf(stdout," Upwind u: %16.14lf\n",upwindedphysfield[1]);
        
        // Compute the fluxes
        Aflux = upwindedphysfield[0] * upwindedphysfield[1];
        p = pext + beta*(sqrt(upwindedphysfield[0]) - sqrt(A_0));
        p_t = 0.5*(upwindedphysfield[1]*upwindedphysfield[1]) + p/rho;				
        uflux =  p_t;
    }    
    
    /**
     *  Print summary routine, calls virtual routine reimplemented in UnsteadySystem
     */
    void PulseWavePropagation::v_PrintSummary(std::ostream &out)
    {
        PulseWaveSystem::v_PrintSummary(out);
    }

}
