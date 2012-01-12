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
    string PulseWavePropagation::className = GetEquationSystemFactory().RegisterCreatorFunction("PulseWavePropagation", PulseWavePropagation::create, "Pulse Wave Propagation equation.");

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
				//inarray in physical space
				Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
				Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);
				
				for (i = 0; i < nvariables; ++i)
				{
					physarray[i] = Array<OneD, NekDouble>(nq);
					Vmath::Vcopy(nq,inarray[i],1,physarray[i],1);
					
					modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
				}
				
				WeakDGAdvection(physarray, modarray, true, true);
				
				// negate the outarray since moving terms to the rhs
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
				// Initialise variables
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
     *
     */
    void PulseWavePropagation::DoOdeProjection(const Array<OneD,const Array<OneD, NekDouble> >&inarray,
											   Array<OneD, Array<OneD, NekDouble> >&outarray,
											   const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();
        SetBoundaryConditions(time);
		
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

	
	/*Calculates the second term of the weak form: dF/dx 
	 *The variables ot the system are (A;u)
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

		//cout << "pext = "<<pext<<"\th0 = "<<h0<<"\tnue = "<<nue<<"\trho = "<<rho<<endl;
		
		//Get A_0 at equilibrium state
		Array<OneD, NekDouble> A_0(nq);
		StaticArea(A_0,0.0);
		for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
		{
			//cout << "A_0 = "<<A_0[j]<<endl;
		}
		

		//Get the material properties of the artery from the inputfile and calculate the beta
		Array<OneD, NekDouble> YoungsModulus(nq);
		Array<OneD, NekDouble> beta(nq);
		MaterialProperties(YoungsModulus,0.0);
		for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
		{
			beta[j] = sqrt(3.1415)*h0*YoungsModulus[j]/((1-nue*nue)*A_0[j]);
			//cout << "beta = "<<beta[j]<<endl;
		}
		
		
        switch (i)
		{
			//flux for A equation	
			case 0:
			{
				for (int j = 0; j < nq; j++)
				{
					flux[0][j] = physfield[0][j]*physfield[1][j]; 
				}
			}
			break;
				
			//flux for u equation	
			case 1:
			{
				for (int j = 0; j < nq; j++)
				{
					ASSERTL0(physfield[0][j]>=0,"Negative A not allowed.");
					p = pext + beta[j]*(sqrt(physfield[0][j]) - sqrt(A_0[j]));
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

	
	/*Calculates the third term of the weak form: numerical flux at boundary
	 *
	 */
    void PulseWavePropagation::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {		
        int i;
		int nTraceNumPoints = GetTraceTotPoints();
		int nvariables      = 2; //(A,u)
		int nq = m_fields[0]->GetNpoints();
		
		NekDouble rho = m_rho; 
		NekDouble pext = m_pext; 
		NekDouble p = 0.0;
		NekDouble p_t = 0.0;
		NekDouble h0 = m_h0; 
		NekDouble nue = m_nue; 
		
		// get temporary arrays
		Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
		Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
		for (i = 0; i < nvariables; ++i)
		{
			Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
			Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		}
		
		
		// get the physical values at the trace
		for (i = 0; i < nvariables; ++i)
		{
			m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
		}
				
		
		/*/print the values in Fwd and Bwd
		 for (int i=0; i<Fwd[0].num_elements(); i++)
		 {
		 cout << "Fwd[A]["<<i<<"] = "<<Fwd[0][i]<<"\t";
		 cout << "Bwd[A]["<<i<<"] = "<<Bwd[0][i]<<"\t\t";
		 
		 cout << "Fwd[u]["<<i<<"] = "<<Fwd[1][i]<<"\t";
		 cout << "Bwd[u]["<<i<<"] = "<<Bwd[1][i]<<endl;
		 }*/
		
		
		//Get A_0 at equilibrium state, !Is hard coded by A_0[0]
		Array<OneD, NekDouble> A_0(nq);
		Array<OneD, NekDouble> A_trace(GetTraceTotPoints());
		StaticArea(A_0,0.0);
		for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
		{
			//cout << "A_0 = "<<A_0[j]<<endl;
		}
		//m_fields[0]->GetFwdBwdTracePhys(A_0,A_trace,A_trace2);
		m_fields[0]->ExtractTracePhys(A_0,A_trace);
		A_trace[GetTraceTotPoints()-1]=A_0[GetTotPoints()-1];
		
		for (int j = 0; j < GetTraceTotPoints(); ++j)
		{
			//cout << "A_trace["<<j<<"] = "<<A_trace[j]<<endl;
		}
		
		
		
		//Get the material properties of the artery from the inputfile and calculate the beta
		Array<OneD, NekDouble> YoungsModulus(nq);
		Array<OneD, NekDouble> beta(nq);
		Array<OneD, NekDouble> beta_trace(GetTraceTotPoints());
		MaterialProperties(YoungsModulus,0.0);
		for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
		{
			beta[j] = sqrt(3.1415)*h0*YoungsModulus[j]/((1-nue*nue)*A_0[j]);
			//cout << j<<"\tbeta = "<<beta[j]<<endl;
		}
		
		m_fields[0]->ExtractTracePhys(beta,beta_trace);
		beta_trace[GetTraceTotPoints()-1]=beta[GetTotPoints()-1];

		for (int j = 0; j < GetTraceTotPoints(); ++j)
		{
			//cout << "beta_trace["<<j<<"] = "<<beta_trace[j]<<endl;
		}
		
		
        // Solve the Riemann problem
        NekDouble Aflux, uflux;
        for (i = 0; i < nTraceNumPoints; ++i)
		{
			switch(m_upwindTypePulse)
			{
				case eUpwindPulse:
				{
					RiemannSolverUpwind(Fwd[0][i],Fwd[1][i],
										Bwd[0][i],Bwd[1][i],
										Aflux, uflux,i, A_trace[i], beta_trace[i]);
				}
					break;
				default:
				{
					ASSERTL0(false,"populate switch statement for upwind flux");
				}
					break;
			}
			
			// rotate back to Cartesian
			numflux[0][i]  = Aflux;
			numflux[1][i] = uflux;
		}
    }
	
	
	/*Upwinding Riemann solver for pulse wave propagaiton
	 * 
	 */
	void PulseWavePropagation::RiemannSolverUpwind(NekDouble AL,NekDouble uL,
												   NekDouble AR,NekDouble uR, 
												   NekDouble &Aflux, NekDouble &uflux, int i, NekDouble A_0, NekDouble beta) //(croth)
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
		
		
		// compute the wave speeds
		cL = sqrt(beta*sqrt(AL)/(2*rho));
		cR = sqrt(beta*sqrt(AR)/(2*rho));
		c_Roe =(cL+cR)/2;
		//cout << "cL = "<<cL<<"\tcR = "<<cR<<"\tc_Roe = "<<c_Roe<<endl;
		
		u_Roe = (uL+uR)/2;
		//cout << "uL = "<<uL<<"\tuR = "<<uR<<"\tu_Roe = "<<u_Roe<<endl;
		
		lambda[0]= u_Roe + c_Roe;
		lambda[1]= u_Roe - c_Roe;
		//cout << "lambda[0] = "<<lambda[0]<<"\tlambda[1] = "<<lambda[1]<<endl;
		
		// calculate the caracteristic variables 
		//left characteristics
		characteristic[0] = uL + 4*sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		characteristic[1] = uL - 4*sqrt(sqrt(AL))*sqrt(beta/(2*rho));
		//right characteristics
		characteristic[2] = uR + 4*sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		characteristic[3] = uR - 4*sqrt(sqrt(AR))*sqrt(beta/(2*rho));
		
		for (int k=0; k<4; k++)
		{
			//cout << "characteristic["<<k<<"] = "<<characteristic[k]<<endl;
		}
		
		//take left or right value of characteristic variable
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
		
		//calculate conservative variables from characteristics
		upwindedphysfield[0]= ((W[0]-W[1])/4)*((W[0]-W[1])/4)*((W[0]-W[1])/4)*((W[0]-W[1])/4)*(rho/(2*beta))*(rho/(2*beta));
		upwindedphysfield[1]= (W[0] + W[1])/2;
		
		
		// compute the fluxes
		Aflux = upwindedphysfield[0] * upwindedphysfield[1];
		p = pext + beta*(sqrt(upwindedphysfield[0]) - sqrt(A_0));
		p_t = (upwindedphysfield[1]*upwindedphysfield[1])/2 + p/rho;				
		uflux =  p_t;
		
	}
	
	
	
	
	/*Gets the Material Properties of the artery
	 * specified in the inputfile
	 */
	void PulseWavePropagation::MaterialProperties(Array<OneD, NekDouble> &YoungsModulus, const NekDouble time)
    {
		int nq = m_fields[0]->GetNpoints();
		std::string velStr[1] = {"E0"};        
		
		LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("MaterialProperties",velStr[0]);
		
		EvaluateFunction(YoungsModulus,ifunc,time);
    }
	
	
	/*Gets the Area at static equilibrium
	 * specified in the inputfile
	 */
	void PulseWavePropagation::StaticArea(Array<OneD, NekDouble> &A_0, const NekDouble time)
    {
		int nq = m_fields[0]->GetNpoints();
		std::string velStr[1] = {"A"};        
		
		LibUtilities::EquationSharedPtr ifunc = m_session->GetFunction("A_0",velStr[0]);
		
		EvaluateFunction(A_0,ifunc,time);
    }
	
	
	/*Handle the pressure boundary condition
	 * if a boundary condition is set for
	 * the inflow pressure
	 */
	void PulseWavePropagation::SetBoundaryConditions_new(NekDouble time)
	{
		
	}
	
	
	/*Print summary routine
	 *
	 */
    void PulseWavePropagation::v_PrintSummary(std::ostream &out)
    {
        PulseWaveSystem::v_PrintSummary(out);
    }


}
