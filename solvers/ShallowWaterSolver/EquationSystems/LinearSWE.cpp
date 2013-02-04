///////////////////////////////////////////////////////////////////////////////
//
// File LinearSWE.cpp
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
// Description: Linearized Shallow water equations in primitive variables
//              valid for a constant water depth
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/LinearSWE.h>

namespace Nektar
{
  string LinearSWE::className = GetEquationSystemFactory().RegisterCreatorFunction("LinearSWE", LinearSWE::create, "Linear shallow water equation in primitive variables.");
  
  LinearSWE::LinearSWE(const LibUtilities::SessionReaderSharedPtr& pSession)
    : ShallowWaterSystem(pSession)
  {
  }

  void LinearSWE::v_InitObject()
  {
      ShallowWaterSystem::v_InitObject();
      
    if (m_explicitAdvection)
      {
	m_ode.DefineOdeRhs     (&LinearSWE::DoOdeRhs,        this);
	m_ode.DefineProjection (&LinearSWE::DoOdeProjection, this);
      }
    else
      {
	ASSERTL0(false, "Implicit SWE not set up.");
      }

  }
  
  LinearSWE::~LinearSWE()
  {
    
  }
 
  void LinearSWE::AddCoriolis(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
			            Array<OneD,       Array<OneD, NekDouble> > &outarray)
  {
   
    int ncoeffs = outarray[0].num_elements();
    int nq      = physarray[0].num_elements();
    
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> h(nq);
    
    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  // add to u equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,tmp);
	  Vmath::Vadd(ncoeffs,tmp,1,outarray[1],1,outarray[1],1);
	  
	  // add to v equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,tmp);
	  Vmath::Vadd(ncoeffs,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	  // add to u equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[1],1,outarray[1],1);
	  
	  // add to v equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the LinearSWE");
	break;
      }
  }
  void LinearSWE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			        Array<OneD,       Array<OneD, NekDouble> >&outarray, 
			  const NekDouble time) 
  {
    int i;
    int ndim    = m_spacedim;
    int nvariables = inarray.num_elements();
    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();
    
    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  //-------------------------------------------------------
	  //inarray in physical space
	  
	  Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);
	  for (i = 0; i < nvariables; ++i)
	    {
	      modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
	    }
	  //-------------------------------------------------------


	  //-------------------------------------------------
	  // get the advection part
	  // input: physical space
	  // output: modal space 
	  
	  // straighforward DG
	  WeakDGAdvection(inarray, modarray, false, true);
	  //-------------------------------------------------

	  
	  //-------------------------------------------------------
	  // negate the outarray since moving terms to the rhs
	  for(i = 0; i < nvariables; ++i)
	    {
	      Vmath::Neg(ncoeffs,modarray[i],1);
	    }
	  //-------------------------------------------------------
	  
	  
	  //-------------------------------------------------
 	  // Add "source terms"
 	  // input: physical space
 	  // output: modal space
	  
 	  // coriolis forcing
 	  if (m_coriolis.num_elements() != 0)
 	    {
 	      AddCoriolis(inarray,modarray);
 	    }
 	  //------------------------------------------------- 
	   
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
	  
	  Array<OneD,NekDouble> tmp(nq);
	  Array<OneD, NekDouble>tmp1(nq);   
	  Array<OneD, Array<OneD, NekDouble> > fluxvector(ndim);
	  for(i = 0; i < ndim; ++i)
	    {
	      fluxvector[i]    = Array<OneD, NekDouble>(nq);
	    }
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      // Get the ith component of the  flux vector in (physical space)
	      LinearSWE::GetFluxVector2D(i, inarray, fluxvector);
	      
	      //m_fields[0]->PhysDeriv(0,fluxvector[0],tmp);
	      //m_fields[0]->PhysDeriv(1,fluxvector[1],tmp1);
	      m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],fluxvector[0],tmp);
	      m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],fluxvector[1],tmp1);
	      Vmath::Vadd(nq,tmp,1,tmp1,1,outarray[i],1);
	      Vmath::Neg(nq,outarray[i],1);
	    }
	  

	  //-------------------------------------------------
 	  // Add "source terms"
	 	  
 	  // coriolis forcing
 	  if (m_coriolis.num_elements() != 0)
 	    {
 	      AddCoriolis(inarray,outarray);
 	    }
 	  //------------------------------------------------- 

	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the LinearSWE");
	break;
      }
  }
 

  void LinearSWE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
				       Array<OneD,       Array<OneD, NekDouble> >&outarray,
				 const NekDouble time)
  {
    int i;
    int nvariables = inarray.num_elements();
   
    
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
	   SetBoundaryConditions(outarray,time);
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	  EquationSystem::SetBoundaryConditions(time);
	  Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
	  
	  for(i = 0; i < nvariables; ++i)
          {
              m_fields[i]->FwdTrans(inarray[i],coeffs);
              m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
          }
	  break;
	}
      default:
	ASSERTL0(false,"Unknown projection scheme");
	break;
      }
  }
  

   //----------------------------------------------------
  void LinearSWE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time)
  {
    
    int nvariables = m_fields.num_elements();
    int cnt = 0;

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Wall Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWall)
	  {
	    if (m_expdim == 1)
	      {
		WallBoundary1D(n,inarray);
	      }
	    else if (m_expdim == 2)
	      {
		WallBoundary2D(n,cnt,inarray);
	      }
	  }
	
	// Time Dependent Boundary Condition (specified in meshfile)
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
	  {
	    for (int i = 0; i < nvariables; ++i)
	      {
		m_fields[i]->EvaluateBoundaryConditions(time);
	      }
	  }
	cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }
  
  //----------------------------------------------------
 
  void LinearSWE::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
  { 

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = physarray.num_elements();
    
    // get physical values of the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
      }
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts;
    
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
  }
  
    void LinearSWE::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {     
    ASSERTL0(false,"1D not working with user defined BC");
  }


    void LinearSWE::GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    
    NekDouble g = m_g;
    
    switch(i){
      
      // flux function for the h equation
    case 0:
      for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	{
	  flux[0][j]  =  physfield[1][j];
	}
      break;
      
      // flux function for the hu equation
    case 1:
      for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	{
	  flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] +
	    0.5*g*physfield[0][j]*physfield[0][j];
	}
      break;
    default:
      ASSERTL0(false,"GetFluxVector1D: illegal vector index");
    }
  }
  
  
  void LinearSWE::GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
				  Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    int nq = m_fields[0]->GetTotPoints();
    
    NekDouble g = m_g;
    
    switch(i){
      
      // flux function for the eta equation
    case 0:
      {
	Vmath::Vmul(nq, m_depth, 1, physfield[1], 1, flux[0], 1);
	Vmath::Vmul(nq, m_depth, 1, physfield[2], 1, flux[1], 1);
      }
      break;
      
      // flux function for the u equation
    case 1:
      {
	Vmath::Smul(nq, g, physfield[0], 1, flux[0], 1);
	Vmath::Zero(nq, flux[1], 1);
      }
      break;
      
      // flux function for the v equation
    case 2:
      {
	Vmath::Zero(nq, flux[0], 1);
	Vmath::Smul(nq, g, physfield[0], 1, flux[1], 1);
      }
      break;
      
    default:
      ASSERTL0(false,"GetFluxVector2D: illegal vector index");
    }
  }
  
  
  void LinearSWE::NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxX)
    {
      ASSERTL0(false,"1D DG not working");
    }
  
  void LinearSWE::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {

    switch(m_upwindType)
      {
      case eNotSet:
	{
	  ASSERTL0(false,"No upwind flux set in the input file");
	  break;
	}
      case eAverage:
	{
	  AverageLinearNumericalFlux2D(physfield, numfluxX, numfluxY);
	  break;
	}
      case eHLL:
	{
	RiemannLinearNumericalFlux2D(physfield, numfluxX, numfluxY);
	break;
	}
      default:
	{
	  ASSERTL0(false,"populate switch statement for upwind flux");
	}
	break;
      }
  }

  void LinearSWE::RiemannLinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					       Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
					       Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {

    int i;
    
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; // only the dependent variables 
    
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
  
  void LinearSWE::AverageLinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						     Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
						     Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {

    int i;
    
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3;
    
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
    
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
      }
    
    // get the physical values at the trace from the dependent variables
    for (i = 0; i < nvariables; ++i)
      {
	m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
      }
    
    NekDouble eta, u, v, d;
    NekDouble g = m_g;
	
	
	// averaging
	for (i = 0; i < nTraceNumPoints; ++i)
	  {
	    eta = 0.5*(Fwd[0][i] + Bwd[0][i]);
	    u   = 0.5*(Fwd[1][i] + Bwd[1][i]);
	    v   = 0.5*(Fwd[2][i] + Bwd[2][i]);
	    d   = m_depth[0]; // constant depth only!
	    
	    numfluxX[0][i]  = d * u;
	    numfluxY[0][i]  = d * v;
	    numfluxX[1][i]  = g * eta;
	    numfluxY[1][i]  = 0.0;
	    numfluxX[2][i]  = 0.0;
	    numfluxY[2][i]  = g * eta;
	  }
		
  }
  

  /***
      
   */
  void LinearSWE::RiemannSolver(NekDouble etaL,NekDouble uL,NekDouble vL,NekDouble etaR,NekDouble uR, 
					    NekDouble vR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux )
  {
    NekDouble g = m_g;
    NekDouble d = m_depth[0]; // only valid for constant depth!

    NekDouble SL,SR;
    NekDouble cL = sqrt(g * d);
    NekDouble cR = sqrt(g * d);
    
    // Compute SL
    SL = uL - cL;
    
    // Compute SR
    SR = uR + cR;
        
    if (SL >= 0)
      {
	hflux  = d * uL;
	huflux  = g * etaL;
	hvflux  = 0.0;
      }
    else if (SR <= 0)
      {
	hflux  = d * uR;
	huflux  = g * etaR;
	hvflux  = 0.0;
      }
    else
      {
	hflux  = (SR*d*uL-SL*d*uR+SL*SR*(etaR-etaL))/(SR-SL);
	huflux = (SR*g*etaL-SL*g*etaR+SL*SR*(uR-uL))/(SR-SL);
	hvflux = 0.0;
      }
  }
    
} //end of namespace

