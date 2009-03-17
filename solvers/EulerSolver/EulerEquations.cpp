///////////////////////////////////////////////////////////////////////////////
//
// File EulerEquations.cpp
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
// Description: Euler Equations class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <EulerSolver/EulerEquations.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
  /**
   * Basic construnctor
   */
  EulerEquations::EulerEquations(void):
    ADRBase(),
    m_infosteps(100)
  {     
  }
  
  /**
   * Constructor. Creates ... of #DisContField2D fields
   *
   * \param 
   * \param
   */
  EulerEquations::EulerEquations(string &fileNameString):
    ADRBase(fileNameString,true),
    m_infosteps(10)
  {
    
    if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
      {
	m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
      }

    if(m_boundaryConditions->CheckForParameter("Gamma") == true)
      {
	m_gamma =  m_boundaryConditions->GetParameter("Gamma");
      }
    else
      {
	ASSERTL0(false,"Gamma not specified");
      }

    // check that any user defined boundary condition is implemented
    
    
  }


  void EulerEquations::ODElhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
				     Array<OneD,       Array<OneD, NekDouble> >&outarray, 
				     const NekDouble time) 
  {
    int nvariables = inarray.num_elements();
    MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
    for(int i = 0; i < nvariables; ++i)
      {
	m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,inarray[i],outarray[i]);
      }
  }
  
  
  void EulerEquations::ODElhsSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
					  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                                 const NekDouble time)   
    {
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
                  m_fields[i]->MultiplyByInvMassMatrix(inarray[i],  
                                                       outarray[i],
                                                       false);
              }
          }
          break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

  void EulerEquations::ODEdirkSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
						  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                                  const NekDouble lambda,
                                                  const NekDouble time) 
  {
    ASSERTL0(false, "this routine needs implementation");
  }
  
  
  void EulerEquations::GeneralTimeIntegration(int nsteps, 
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
	fields[i]  = m_fields[i]->UpdateCoeffs();
      }
    
    if(m_projectionType==eGalerkin)
      {
	// calculate the variable u* = Mu
	// we are going to TimeIntegrate this new variable u*
	MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
	for(int i = 0; i < nvariables; ++i)
	  {
	    tmp[i] = Array<OneD, NekDouble>(ncoeffs);
	    m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
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
      case LibUtilities::eIMEXdirk_3_4_3:
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
	    ASSERTL0(false,"CG not implemented for Euler");
              //   // Project the solution u* onto the boundary conditions to
//                 // obtain the actual solution
//                 SetBoundaryConditions(m_time);
//                 for(i = 0; i < nvariables; ++i)
//                 {
//                     m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i],false);
//                     fields[i] = tmp[i];	   		    
//                 }
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
  
  
  
  
  void EulerEquations::ODErhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			            Array<OneD,       Array<OneD, NekDouble> >&outarray, 
			      const NekDouble time) 
  {
    int i;
    int nvariables = inarray.num_elements();
    int ncoeffs    = inarray[0].num_elements();
    int nq         = GetTotPoints();
    
    //-------------------------------------------------------
    // go to physical space
    
    Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	physarray[i] = Array<OneD, NekDouble>(nq);
	m_fields[i]->BwdTrans(inarray[i],physarray[i]);
      }
    //-------------------------------------------------------
    
    SetBoundaryConditions(physarray, time);

    switch(m_projectionType)
      {
      case eDiscontinuousGalerkin:
	{
	  
	  WeakDGAdvection(physarray, outarray, false, true);

	  // any source terms should be added here
	  
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
	      Vmath::Neg(ncoeffs,outarray[i],1);
	    }
	}
	break;
      case eGalerkin:
	ASSERTL0(false,"Continouos scheme not implemented for Euler");
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme");
	break;
      }
  }


    
  //   void EulerEquations::ExplicitlyIntegrateAdvection(int nsteps)
//     {
//         int i,n,nchk = 0;
//         int ncoeffs = m_fields[0]->GetNcoeffs();
//         int nvariables = m_fields.num_elements();

//         // Get Integration scheme details
//         LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eForwardEuler);
//         LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];


// 	// HACK!!! 
// 	//hardcoded initial conditions for Isentropic vortex
// 	 SetIsenTropicVortex();

//         // Set up wrapper to fields data storage. 
//         Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        
//         for(i = 0; i < nvariables; ++i)
//         {
//             fields[i] = m_fields[i]->UpdateCoeffs();
//         }
                
//          int nInitSteps;
// 	 LibUtilities::TimeIntegrationSolutionSharedPtr u = IntScheme->InitializeScheme(m_timestep,m_time,nInitSteps,*this,fields);


// 	 for(n = nInitSteps; n < nsteps; ++n)
// 	  {
	    
// 	    //----------------------------------------------
//             // Perform time step integration
//             //----------------------------------------------

//             fields = IntScheme->ExplicitIntegration(m_timestep,*this,u);
	   
// 	    m_time += m_timestep;
//             //----------------------------------------------

//             //----------------------------------------------
//             // Dump analyser information
//             //----------------------------------------------
//             if(!((n+1)%m_infosteps))
//             {
// 	      cout << "Steps: " << n+1 << "\t Time: " <<m_time<< endl;
//             }
            
//             if(n&&(!((n+1)%m_checksteps)))
//             {
// 	      for(i = 0; i < nvariables; ++i)
// 		{
// 		  (m_fields[i]->UpdateCoeffs()) = fields[i];
// 		}
// 	      Checkpoint_Output(nchk++);
//             }
// 	  }
        
// 	for(i = 0; i < nvariables; ++i)
// 	  {
// 	    (m_fields[i]->UpdateCoeffs()) = fields[i];
// 	  }
//     }
    
  
  //----------------------------------------------------
  void EulerEquations::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time)
  {
    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Wall Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
	  {
	    WallBoundary(n,physarray);
	  }
	
	// Time dependent Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
	  {
	    for (int i = 0; m_fields.num_elements(); ++i)
	      {
		m_fields[i]->EvaluateBoundaryConditions(time);
	      }
	  }
      }
  }
  
  //----------------------------------------------------
 
  void EulerEquations::WallBoundary(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  { 

     // get physical values of rho, rho u, rho v and E for the forward trace
    Array<OneD, NekDouble> rho(GetTraceTotPoints());
    Array<OneD, NekDouble> rhou(GetTraceTotPoints());
    Array<OneD, NekDouble> rhov(GetTraceTotPoints());
    Array<OneD, NekDouble> E(GetTraceTotPoints());
    m_fields[0]->ExtractTracePhys(physarray[0],rho);
    m_fields[1]->ExtractTracePhys(physarray[1],rhou);
    m_fields[2]->ExtractTracePhys(physarray[2],rhov);
    m_fields[3]->ExtractTracePhys(physarray[3],E);

    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts, cnt = 0; 
    
    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
	
	Array<OneD, NekDouble> tmp_n(npts);
	Array<OneD, NekDouble> tmp_t(npts);
	  
	// rotate to compute the normal and tangential flux components
	Vmath::Vmul(npts,&rhou[id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
	Vmath::Vvtvp(npts,&rhov[id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	
	Vmath::Vmul(npts,&rhou[id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
	Vmath::Vvtvm(npts,&rhov[id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	  
	// negate the normal flux
	Vmath::Neg(npts,tmp_n,1);		      
	
	// rotate back to Cartesian
	Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&rhou[id2],1);
	Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&rhou[id2],1,&rhou[id2],1);
	
	Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&rhov[id2],1);
	Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&rhov[id2],1,&rhov[id2],1);
	
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&rho[id2], 1,&(m_fields[0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	Vmath::Vcopy(npts,&rhou[id2],1,&(m_fields[1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	Vmath::Vcopy(npts,&rhov[id2],1,&(m_fields[2]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	Vmath::Vcopy(npts,&E[id2],1,&(m_fields[3]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
      }
    cnt +=e;
  }
  

  void EulerEquations::GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux)
    {
      NekDouble gamma = m_gamma;
        
      Array<OneD, NekDouble>pressure(m_fields[0]->GetTotPoints());

      switch(i){
	
	// flux function for the \rho equation
      case 0:
	
	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
	    flux[0][j]  =  physfield[1][j];
	    flux[1][j]  =  physfield[2][j];
	  }
	break;
	
	// flux function for the \rho u equation
      case 1:
	
	GetPressure(physfield,pressure);


	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
	    flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] + pressure[j];
	    flux[1][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	  }
	break;
	
	// flux function for the \rho v equation
      case 2:

	GetPressure(physfield,pressure);

	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
	    flux[0][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	    flux[1][j] = physfield[2][j]*physfield[2][j]/physfield[0][j] + pressure[j];
	  }
	break;

	// flux function for the E equation
      case 3:
	
	GetPressure(physfield,pressure);
	
	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
	    flux[0][j] = (physfield[1][j]/physfield[0][j])*(physfield[3][j] + pressure[j]);
	    flux[1][j] = (physfield[2][j]/physfield[0][j])*(physfield[3][j] + pressure[j]);
	  }
	break;
	
      default:
	ASSERTL0(false,"GetFluxVector: illegal vector index");
      }
  }
  	
  void EulerEquations::GetPressure(Array<OneD, Array<OneD, NekDouble> > &physfield,
				   Array<OneD, NekDouble> &pressure)
  {
    NekDouble gamma = m_gamma;

    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	pressure[i] = (gamma - 1.0)*(physfield[3][i] - 0.5*(physfield[1][i]*physfield[1][i]/physfield[0][i] +
							    physfield[2][i]*physfield[2][i]/physfield[0][i]));
      }
  }
  
  void EulerEquations::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numfluxX, Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int i;

        int nTraceNumPoints = GetTraceTotPoints();
	int nvariables      = m_fields.num_elements();

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
	
	// rotate the values to the normal direction
        NekDouble tmpX, tmpY;
        for (i = 0; i < nTraceNumPoints; ++i)
	  {
            tmpX =  Fwd[1][i]*m_traceNormals[0][i]+Fwd[2][i]*m_traceNormals[1][i];
            tmpY = -Fwd[1][i]*m_traceNormals[1][i]+Fwd[2][i]*m_traceNormals[0][i];
            Fwd[1][i] = tmpX;
            Fwd[2][i] = tmpY;
	    
            tmpX =  Bwd[1][i]*m_traceNormals[0][i]+Bwd[2][i]*m_traceNormals[1][i];
            tmpY = -Bwd[1][i]*m_traceNormals[1][i]+Bwd[2][i]*m_traceNormals[0][i];
            Bwd[1][i] = tmpX;
            Bwd[2][i] = tmpY;
	  }
	
      // Solve the Riemann problem
        NekDouble rhoflux, rhouflux, rhovflux, Eflux;
      
        for (int i = 0; i < nTraceNumPoints; ++i)
	{
            RiemannSolver(Fwd[0][i],Fwd[1][i],Fwd[2][i],Fwd[3][i],
                          Bwd[0][i],Bwd[1][i],Bwd[2][i],Bwd[3][i],
                          rhoflux, rhouflux, rhovflux, Eflux );
	  
            // rotate back to Cartesian
            numfluxX[0][i] =  rhoflux*m_traceNormals[0][i];
            numfluxY[0][i] =  rhoflux*m_traceNormals[1][i];
            numfluxX[1][i] = (rhouflux*m_traceNormals[0][i] - rhovflux*m_traceNormals[1][i]) * m_traceNormals[0][i];
            numfluxY[1][i] = (rhouflux*m_traceNormals[0][i] - rhovflux*m_traceNormals[1][i]) * m_traceNormals[1][i];
            numfluxX[2][i] = (rhouflux*m_traceNormals[1][i] + rhovflux*m_traceNormals[0][i]) * m_traceNormals[0][i];
            numfluxY[2][i] = (rhouflux*m_traceNormals[1][i] + rhovflux*m_traceNormals[0][i]) * m_traceNormals[1][i];
            numfluxX[3][i] =  Eflux*m_traceNormals[0][i];
            numfluxY[3][i] =  Eflux*m_traceNormals[1][i];
	}
    }
  
   void EulerEquations::RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
				       NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
				       NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux)
    {
  
        NekDouble gamma = m_gamma;
      
        NekDouble uL = rhouL/rhoL;
        NekDouble vL = rhovL/rhoL;
        NekDouble uR = rhouR/rhoR;
        NekDouble vR = rhovR/rhoR;
        NekDouble pL = (gamma - 1.0) * (EL - 0.5 * (rhouL*uL + rhovL*vL));
        NekDouble pR = (gamma - 1.0) * (ER - 0.5 * (rhouR*uR + rhovR*vR));
        NekDouble cL = sqrt(gamma * pL / rhoL);
        NekDouble cR = sqrt(gamma * pR / rhoR);
        NekDouble hL = (EL+pL)/rhoL;
        NekDouble hR = (ER+pR)/rhoR;
 
        // compute Roe averages
        NekDouble rhoRoe = sqrt(rhoL) * sqrt(rhoR);
        NekDouble uRoe   = (sqrt(rhoL)*uL + sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble vRoe   = (sqrt(rhoL)*vL + sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble hRoe   = (sqrt(rhoL)*hL + sqrt(rhoR)*hR)/(sqrt(rhoL)+sqrt(rhoR));
        NekDouble cRoe   = sqrt( (gamma - 1.0)*(hRoe - 0.5*(uRoe*uRoe + vRoe*vRoe)) );

        // compute the wave speeds
        NekDouble SL = min(uL-cL, uRoe-cRoe);
        NekDouble SR = max(uR+uR, uRoe+cRoe);

        // compute the HLL flux
        if (SL >= 0)
	{
            rhoflux  = rhoL*uL;
            rhouflux = rhoL*uL*uL + pL;
            rhovflux = rhoL*uL*vL;
            Eflux    = uL*(EL + pL);
	}
        else if (SR <= 0)
	{
            rhoflux  = rhoR*uR;
            rhouflux = rhoR*uR*uR + pR;
            rhovflux = rhoR*uR*vR;
            Eflux    = uR*(ER + pR);
	}
        else
	{
            rhoflux  = (( SR*(rhoL*uL) - SL*(rhoR*uR)+SR*SL*(rhoR-rhoL))/(SR-SL) );
            rhouflux = (( SR*(rhoL*uL*uL + pL) - SL*(rhoR*uR*uR + pR) +
                          SR*SL*(rhouR-rhouL)) / (SR-SL) );
            rhovflux = (( SR*(rhoL*uL*vL) - SL*(rhoR*uR*vR) +
                          SR*SL*(rhovR-rhovL)) / (SR-SL) );
            Eflux    = (( SR*(uL*EL+uL*pL) - SL*(uR*ER+uR*pR) +
                          SR*SL*(ER-EL)) / (SR-SL) );
	}

    }

    void EulerEquations::Summary(std::ostream &out)
    {
      cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : Compressible Euler Equations" << endl;
      ADRBase::Summary(out);
      cout << "=======================================================================" << endl;
      cout << endl;
    
    }


  void EulerEquations::SetIsenTropicVortex(void)
  {
    int nTotQuadPoints  = GetTotPoints();
    
    Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
    Array<OneD, NekDouble> rhou(nTotQuadPoints);
    Array<OneD, NekDouble> rhov(nTotQuadPoints);
    Array<OneD, NekDouble> E(nTotQuadPoints);
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);
  
    m_fields[0]->GetCoords(x,y,z);

    //---------------------------------
    // flow parameters
    NekDouble x0   = 5.0;
    NekDouble y0   = 0.0;
    NekDouble beta  = 5.0;
    NekDouble u0    = 1.0;
    NekDouble v0    = 0.0;
    NekDouble gamma = m_gamma;
    NekDouble time  = m_time;
    NekDouble r;

    for (int i = 0; i < nTotQuadPoints; ++i)
      {
        r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
        rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
        rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
        rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
        E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
    }

    m_fields[0]->SetPhys(rho);
    m_fields[1]->SetPhys(rhou);
    m_fields[2]->SetPhys(rhov);
    m_fields[3]->SetPhys(E);

    // forward transform to fill the modal coeffs
    for(int i = 0; i < m_fields.num_elements(); ++i)
      {
	m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
      }
  }
  
  void EulerEquations::GetExactIsenTropicVortex(Array<OneD, NekDouble> &outarray, int field)
  {
    int nTotQuadPoints  = GetTotPoints();
  
    Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
    Array<OneD, NekDouble> rhou(nTotQuadPoints);
    Array<OneD, NekDouble> rhov(nTotQuadPoints);
    Array<OneD, NekDouble> E(nTotQuadPoints);
    Array<OneD, NekDouble> x(nTotQuadPoints);
    Array<OneD, NekDouble> y(nTotQuadPoints);
    Array<OneD, NekDouble> z(nTotQuadPoints);
  
    m_fields[0]->GetCoords(x,y,z);
  
    //---------------------------------
    // flow parameters
    NekDouble x0   = 5.0;
    NekDouble y0   = 0.0;
    NekDouble beta  = 5.0;
    NekDouble u0    = 1.0;
    NekDouble v0    = 0.0;
    NekDouble gamma = m_gamma;
    NekDouble time  = m_time;
    NekDouble r;

    for (int i = 0; i < nTotQuadPoints; ++i)
    {
        r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
        rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
        rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
        rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
        E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
    }

    switch (field){
    case 0:
        outarray = rho;
        break;
    case 1:
        outarray = rhou;
        break;
    case 2:
        outarray = rhov;
        break;
    case 3:
        outarray = E;
        break;
    }
    
  }
  
} //end of namespace

/**
* $Log: EulerEquations.cpp,v $
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.1  2008/11/17 08:42:06  claes
* Initial commit of restructured Euler Solver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
