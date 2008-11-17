///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterEquations.cpp
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
// Description: Shallow Water Equations class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <ShallowWaterSolver/ShallowWaterEquations.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
  /**
   * Basic construnctor
   */
  ShallowWaterEquations::ShallowWaterEquations(void):
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
  ShallowWaterEquations::ShallowWaterEquations(string &fileNameString):
    ADRBase(fileNameString,true),
    m_infosteps(10)
  {
    
    // Set Coriolis forcing if specified
    if(m_boundaryConditions->CheckForParameter("Coriolis") == true)
      {
	if (m_expdim == 2)
	  {
	    m_coriolis = Array<OneD, NekDouble> (GetPointsTot());
	    
	    EvaluateCoriolis();
	  }
	else
	  {
	    ASSERTL0(false,"Coriolis defined for 1D run");
	  }
      }
    
    
    if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
      {
	m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
      }
    
    
    if(m_boundaryConditions->CheckForParameter("Gravity") == true)
      {
	m_g =  m_boundaryConditions->GetParameter("Gravity");
      }
    else
      {
	ASSERTL0(false,"Gravity not specified");
      }

    // check that any user defined boundary condition is indeed implemented
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	// If no User Defined then this entry is empty
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "")
	  {
	    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "TimeDependent" &&
		m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "Wall"	)
	      {
		ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition.Implemented are:\n TimeDependent \n Wall \n Please check the spelling in the input file");
	      }
	  }
      }
  }
	  
  void ShallowWaterEquations::EvaluateCoriolis(void)
  {
    int nq = m_fields[0]->GetPointsTot();
    
    std::string coriolisStr[1] = {"f"};
    
    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);
    
    // get the coordinates (assuming all fields have the same
    // discretisation)
    m_fields[0]->GetCoords(x0,x1,x2);
    
    SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(coriolisStr[0]);
    
    for(int j = 0; j < nq; j++)
      {
	m_coriolis[j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
      }
  }
  
  
  void ShallowWaterEquations::ODEforcing(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
					 Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time) 
  {
    int i;
    int nVelDim    = m_spacedim;
    int nvariables = inarray.num_elements();
    int ncoeffs    = inarray[0].num_elements();
    int nq         = GetPointsTot();
    
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
	  
	  //--------------------------------------------
	  // Add source terms here

	  //coriolis forcing
	  if (m_coriolis[0])
	    AddCoriolis(physarray,outarray);
	  //--------------------------------------------
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
	      
	      Vmath::Neg(ncoeffs,outarray[i],1);
	    }
	}
	break;
      case eGalerkin:
	{
	
	  Array<OneD, Array<OneD, NekDouble> > flux(nVelDim);
	  
	  for(i = 0; i < nVelDim; ++i)
	    {
	      flux[i] = Array<OneD, NekDouble>(nq);
	    }
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      GetFluxVector(i, physarray, flux);
	      WeakAdvectionDivergenceForm(flux,outarray[i]);
	      
	    }
	  // --------------------------------------------
	  // Add source terms here
	  
	  //coriolis forcing
	  if (m_coriolis[0])
	    AddCoriolis(physarray,outarray);
	  //--------------------------------------------
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      // Multiply by inverse of mass matrix to get forcing term
	    //   m_fields[i]->MultiplyByInvMassMatrix(outarray[i],  
// 						   outarray[i],
// 						   false, true);
	      Vmath::Neg(ncoeffs,outarray[i],1);
	      
	    }
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme");
	break;
      }
  }
  
    void ShallowWaterEquations::ExplicitlyIntegrateAdvection(int nsteps)
    {
      int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Get Integration scheme details
        LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eForwardEuler);
        LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
	Array<OneD, Array<OneD, NekDouble> >   in(nvariables);
	Array<OneD, Array<OneD, NekDouble> >   out(nvariables);
	Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
	Array<OneD, Array<OneD, NekDouble> >   phys(nvariables);
        
        for(i = 0; i < nvariables; ++i)
        {
            fields[i]  = m_fields[i]->UpdateCoeffs();
	    in[i] = Array<OneD, NekDouble >(ncoeffs);
	    out[i] = Array<OneD, NekDouble >(ncoeffs);
	    tmp[i] = Array<OneD, NekDouble >(ncoeffs);
	    phys[i] = Array<OneD, NekDouble>(m_fields[0]->GetPointsTot());
	    Vmath::Vcopy(ncoeffs,m_fields[i]->GetCoeffs(),1,in[i],1);
        }
                
        int nInitSteps;
        LibUtilities::TimeIntegrationSolutionSharedPtr u = IntScheme->InitializeScheme(m_timestep,m_time,nInitSteps,*this,fields);

        for(n = nInitSteps; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
 
	  switch(m_projectionType)
	    {
	    case eDiscontinuousGalerkin:
	      fields = IntScheme->ExplicitIntegration(m_timestep,*this,u);
	      break;
	    case eGalerkin:
	      {
		//---------------------------------------------------------
		// this is just a forward Euler to illustate that CG works
		 
		// get -D u^n
		ODEforcing(in,out,m_time); // note that MultiplyByInvMassMatrix is not performed inside ODEforcing
	  
		// compute M u^n
		for (i = 0; i < nvariables; ++i)
		  {
		    m_fields[0]->BwdTrans(in[i],phys[i]);
		    m_fields[0]->IProductWRTBase(phys[i],tmp[i]);
		    
		    // f = M u^n - Dt D u^n
		    Vmath::Svtvp(ncoeffs,m_timestep,out[i],1,tmp[i],1,tmp[i],1);
		    
		    // u^{n+1} = M^{-1} f
		    m_fields[i]->MultiplyByInvMassMatrix(tmp[i],out[i],false,false);
		    
		    // fill results
		    Vmath::Vcopy(ncoeffs,out[i],1,in[i],1);
		    Vmath::Vcopy(ncoeffs,out[i],1,fields[i],1);
		  }
		//---------------------------------------------------------
	      }
	      break;
	    }
	  m_time += m_timestep;
            //----------------------------------------------

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
  void ShallowWaterEquations::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time)
  {
    
    int nvariables = m_fields.num_elements();
    
    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Wall Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
	  {
	    if (m_expdim == 1)
	      {
		WallBoundary1D(n,inarray);
	      }
	    else if (m_expdim == 2)
	      {
		WallBoundary2D(n,inarray);
	      }
	  }
	
	// Time Dependent Boundary Condition (specified in meshfile)
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
	  {
	    for (int i = 0; i < nvariables; ++i)
	      {
		m_fields[i]->EvaluateBoundaryConditions(time);
	      }
	  }
      }
  }
  
  //----------------------------------------------------
 
  void ShallowWaterEquations::WallBoundary2D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  { 

    int i;
    int nTraceNumPoints = GetTracePointsTot();
    int nvariables      = m_fields.num_elements();
    
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
      }
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts, cnt = 0; 
    
    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
	
	switch(m_expdim)
	  {
	  case 1:
	    {
	      // negate the forward flux
	      Vmath::Neg(npts,&Fwd[1][id2],1);	
	      
	      //ASSERTL0(false,"1D not yet implemented for SWE");
	    }
	    break;
	  case 2:
	    {
	      Array<OneD, NekDouble> tmp_n(npts);
	      Array<OneD, NekDouble> tmp_t(npts);
	      
	      Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
	      Vmath::Vvtvp(npts,&Fwd[2][id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	      
	      Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
	      Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	      
	      // negate the normal flux
	      Vmath::Neg(npts,tmp_n,1);		      
	      
	      // rotate back to Cartesian
	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
	      Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);
	      
	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
	      Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
	    }
	    break;
	  case 3:
	    ASSERTL0(false,"3D not implemented for Shallow Water Equations");
	    break;
	  default:
	    ASSERTL0(false,"Illegal expansion dimension");
	  }

	// copy boundary adjusted values into the boundary expansion
	for (i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(npts,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	  }
      }
    cnt +=e;
  }
  
    void ShallowWaterEquations::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  { 
    ASSERTL0(false,"1D not working with user defined BC");
  }


  void ShallowWaterEquations::AddCoriolis(Array<OneD, Array<OneD, NekDouble> > &physarray,
					  Array<OneD, Array<OneD, NekDouble> > &outarray)
  {
    
    int ncoeffs = outarray[0].num_elements();
    int nq      = physarray[0].num_elements();
    
    Array<OneD, NekDouble> tmp(nq);
    
    // add to hu equation
    Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
    m_fields[0]->IProductWRTBase(tmp,tmp);
    Vmath::Neg(nq,tmp,1); 
    Vmath::Vadd(ncoeffs,tmp,1,outarray[1],1,outarray[1],1);
    
    // add to hv equation
    Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
    //Vmath::Neg(nq,tmp1,1);
    m_fields[0]->IProductWRTBase(tmp,tmp);
    Vmath::Vadd(ncoeffs,tmp,1,outarray[2],1,outarray[2],1);
    
  }


  void ShallowWaterEquations::GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    
    NekDouble g = m_g;
    
    switch(i){
      
      // flux function for the h equation
    case 0:
      for (int j = 0; j < m_fields[0]->GetPointsTot(); ++j)
	{
	  flux[0][j]  =  physfield[1][j];
	}
      break;
      
      // flux function for the hu equation
    case 1:
      for (int j = 0; j < m_fields[0]->GetPointsTot(); ++j)
	{
	  flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] +
	    0.5*g*physfield[0][j]*physfield[0][j];
	}
      break;
    default:
      ASSERTL0(false,"GetFluxVector1D: illegal vector index");
    }
  }
  
  
  void ShallowWaterEquations::GetFluxVector2D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    
    NekDouble g = m_g;
    
    switch(i){
      
      // flux function for the h equation
    case 0:
      for (int j = 0; j < m_fields[0]->GetPointsTot(); ++j)
	{
	  flux[0][j]  =  physfield[1][j];
	  flux[1][j]  =  physfield[2][j];
	}
      break;
      
      // flux function for the hu equation
      case 1:
	for (int j = 0; j < m_fields[0]->GetPointsTot(); ++j)
	  {
	    flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] +
	      0.5*g*physfield[0][j]*physfield[0][j];
	    flux[1][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	  }
	break;
	
	// flux function for the hv equation
    case 2:
      for (int j = 0; j < m_fields[0]->GetPointsTot(); ++j)
	{
	  flux[0][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	  flux[1][j] = physfield[2][j]*physfield[2][j]/physfield[0][j] +
	    0.5*g*physfield[0][j]*physfield[0][j];
	}
      break;
      
    default:
      ASSERTL0(false,"GetFluxVector2D: illegal vector index");
    }
  }
  
  void ShallowWaterEquations::NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxX)
    {
      ASSERTL0(false,"1D DG not working");
    }
  
  void ShallowWaterEquations::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int i;

        int nTraceNumPoints = GetTracePointsTot();
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
        NekDouble hflux, huflux, hvflux;
	
        for (i = 0; i < nTraceNumPoints; ++i)
	  {
            RiemannSolver(Fwd[0][i],Fwd[1][i],Fwd[2][i],
                          Bwd[0][i],Bwd[1][i],Bwd[2][i],
                          hflux, huflux, hvflux );
	    
            // rotate back to Cartesian
            numfluxX[0][i]  = hflux*m_traceNormals[0][i];
            numfluxY[0][i]  = hflux*m_traceNormals[1][i];
            numfluxX[1][i] = (huflux*m_traceNormals[0][i] - hvflux*m_traceNormals[1][i]) * m_traceNormals[0][i];
            numfluxY[1][i] = (huflux*m_traceNormals[0][i] - hvflux*m_traceNormals[1][i]) * m_traceNormals[1][i];
            numfluxX[2][i] = (huflux*m_traceNormals[1][i] + hvflux*m_traceNormals[0][i]) * m_traceNormals[0][i];
            numfluxY[2][i] = (huflux*m_traceNormals[1][i] + hvflux*m_traceNormals[0][i]) * m_traceNormals[1][i];
	  }
    }
  
  void ShallowWaterEquations::RiemannSolver(NekDouble hL,NekDouble huL,NekDouble hvL,NekDouble hR,NekDouble huR, 
					    NekDouble hvR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux )
  {
    
        NekDouble hC,huC,hvC,SL,SR,hstar,Sstar;
      
        NekDouble g = m_g;
      
        NekDouble uL = huL/hL;
        NekDouble vL = hvL/hL;
        NekDouble uR = huR/hR;
        NekDouble vR = hvR/hR;
        NekDouble cL = sqrt(g * hL);
        NekDouble cR = sqrt(g * hR);
      
        // the two-rarefaction wave assumption
        hstar = 0.5*(cL + cR) + 0.25*(uL - uR);
        hstar *= hstar;
        hstar *= (1.0/g);
      
        // Compute SL
        if (hstar > hL)
            SL = uL - cL * sqrt(0.5*((hstar*hstar + hstar*hL)/(hL*hL)));
        else
            SL = uL - cL;
      
        // Compute SR
        if (hstar > hR)
            SR = uR + cR * sqrt(0.5*((hstar*hstar + hstar*hR)/(hR*hR)));
        else
            SR = uR + cR;
      
      
        if (fabs(hR*(uR-SR)-hL*(uL-SL)) <= 1.0e-15)
            Sstar = 0.0;
        else
            Sstar = (SL*hR*(uR-SR)-SR*hL*(uL-SL))/(hR*(uR-SR)-hL*(uL-SL));
      
        if (SL >= 0)
	{
            hflux  = hL * uL;
            huflux  = uL * uL * hL + 0.5 * g * hL * hL;
            hvflux  = hL * uL * vL;
	}
        else if (SR <= 0)
	{
            hflux  = hR * uR;
            huflux  = uR * uR * hR + 0.5 * g * hR * hR;
            hvflux  = hR * uR *vR;
	}
        else
	{
            if ((SL < 0) && (Sstar >= 0))
	    {
                hC  = hL * ((SL - uL) / (SL - Sstar));
                huC = hC * Sstar;
                hvC = hC * vL;
	      
                hflux = hL*uL + SL * (hC - hL);
                huflux = (uL*uL*hL+0.5*g*hL*hL) + SL * (huC - hL*uL);
                hvflux = (uL*vL*hL) + SL * (hvC - hL*vL);
	    }
            else
	    {
                hC  = hR * ((SR - uR) / (SR - Sstar));
                huC = hC * Sstar;
                hvC = hC * vR;
	      
                hflux = hR*uR + SR * (hC - hR);
                huflux = (uR*uR*hR+0.5*g*hR*hR) + SR * (huC - hR*uR);
                hvflux = (uR*vR*hR) + SR * (hvC - hR*vR);
	    }
	}
      
    }

    void ShallowWaterEquations::Summary(std::ostream &out)
    {
      cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : Shallow Water Equations" << endl;
      ADRBase::Summary(out);
      cout << "=======================================================================" << endl;
      cout << endl;
    
    }


} //end of namespace

/**
* $Log: ShallowWaterEquations.cpp,v $
* Revision 1.2  2008/09/15 14:54:15  claes
* Small changes associated with the BoussinesqSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
