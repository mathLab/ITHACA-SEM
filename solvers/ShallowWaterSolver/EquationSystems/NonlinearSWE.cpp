///////////////////////////////////////////////////////////////////////////////
//
// File NonlinearSWE.cpp
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
// Description: Nonlinear Shallow water equations in conservative variables
//              valid for a constant water depth
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/NonlinearSWE.h>

namespace Nektar
{
  string NonlinearSWE::className = GetEquationSystemFactory().RegisterCreatorFunction("NonlinearSWE", NonlinearSWE::create, "Nonlinear shallow water equation in conservative variables.");
  
  NonlinearSWE::NonlinearSWE(
          const LibUtilities::SessionReaderSharedPtr& pSession)
    : ShallowWaterSystem(pSession)
  {
  }

  void NonlinearSWE::v_InitObject()
  {
      ShallowWaterSystem::v_InitObject();
      
    if (m_explicitAdvection)
      {
	m_ode.DefineOdeRhs     (&NonlinearSWE::DoOdeRhs,        this);
	m_ode.DefineProjection (&NonlinearSWE::DoOdeProjection, this);
      }
    else
      {
	ASSERTL0(false, "Implicit SWE not set up.");
      }
  }
  
  NonlinearSWE::~NonlinearSWE()
  {
    
  }
 

  void NonlinearSWE::AddCoriolis(Array<OneD, Array<OneD, NekDouble> > &physarray,
			      Array<OneD, Array<OneD, NekDouble> > &outarray)
  {
     // routine works for both primitive and conservative formulations
    int ncoeffs = outarray[0].num_elements();
    int nq      = physarray[0].num_elements();
    
    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> h(nq);
	
    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  
	  // conservative formulation compute h
	  // h = \eta + d
	  Vmath::Vadd(nq,physarray[0],1,m_depth,1,h,1);
	  
	  // add to hu equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  Vmath::Vmul(nq,h,1,tmp,1,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,tmp);
	  Vmath::Vadd(ncoeffs,tmp,1,outarray[1],1,outarray[1],1);
	  
	  // add to hv equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Vmul(nq,h,1,tmp,1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,tmp);
	  Vmath::Vadd(ncoeffs,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	    // conservative formulation compute h
	  // h = \eta + d
	  Vmath::Vadd(nq,physarray[0],1,m_depth,1,h,1);
	  
	  // add to hu equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  Vmath::Vmul(nq,h,1,tmp,1,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[1],1,outarray[1],1);
	  
	  // add to hv equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Vmul(nq,h,1,tmp,1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the NonlinearSWE");
	break;
      }
  }
  
  void NonlinearSWE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
	  
	  Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
	  Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);

	  for (i = 0; i < nvariables; ++i)
	    {
	      physarray[i] = Array<OneD, NekDouble>(nq);
	      modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
	    }

	  ConservativeToPrimitive(inarray,physarray);
	  //-------------------------------------------------------
	  	  
	  
	  //-------------------------------------------------
	  // get the advection part
	  // input: physical space
	  // output: modal space 
	  
	  // straighforward DG
	  WeakDGAdvection(physarray, modarray, false, true);
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
	      AddCoriolis(physarray,modarray);
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
	  Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
	  Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);

	  for (i = 0; i < nvariables; ++i)
	    {
	      physarray[i] = Array<OneD, NekDouble>(nq);
	      modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
	    }

	  ConservativeToPrimitive(inarray,physarray);
	  
	  Array<OneD, Array<OneD, NekDouble> > fluxvector(ndim);
	  for(i = 0; i < ndim; ++i)
	    {
	      fluxvector[i]    = Array<OneD, NekDouble>(nq);
	    }
	  
	  Array<OneD,NekDouble> tmp(nq);
	  Array<OneD, NekDouble>tmp1(nq);           
	  
	  for(i = 0; i < nvariables; ++i)
	    {
	      // Get the ith component of the  flux vector in (physical space)
	      NonlinearSWE::GetFluxVector2D(i, physarray, fluxvector);
	      
	      //m_fields[0]->PhysDeriv(0,fluxvector[0],tmp);
	      //m_fields[0]->PhysDeriv(1,fluxvector[1],tmp1);
	      m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],fluxvector[0],tmp);
	      m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],fluxvector[1],tmp1);
	      Vmath::Vadd(nq,tmp,1,tmp1,1,outarray[i],1);
	      Vmath::Neg(nq,outarray[i],1);
	    }
	  
	  //-------------------------------------------------
	  // Add "source terms"
	  // input: physical space
	  // output: modal space
	  
	  // coriolis forcing
	  if (m_coriolis.num_elements() != 0)
	    {
	      AddCoriolis(physarray,outarray);
	    }
	  //------------------------------------------------- 
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the NonlinearSWE");
	break;
      }
  }
 

  void NonlinearSWE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
				           Array<OneD,       Array<OneD, NekDouble> >&outarray,
				     const NekDouble time)
  {
    int i;
    int nvariables = inarray.num_elements();
   
    
    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  ConservativeToPrimitive(inarray,outarray);
	  SetBoundaryConditions(outarray,time);
	  PrimitiveToConservative(outarray,outarray);
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	  ConservativeToPrimitive(inarray,outarray);
	  EquationSystem::SetBoundaryConditions(time);
	  Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
	  
	  for(i = 0; i < nvariables; ++i)
          {
              m_fields[i]->FwdTrans(outarray[i],coeffs);
	      m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
          }
	  PrimitiveToConservative(outarray,outarray);
	  break;
	}
      default:
	ASSERTL0(false,"Unknown projection scheme");
	break;
      }
  }
  

   //----------------------------------------------------
  void NonlinearSWE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time)
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
 
  void NonlinearSWE::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
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
  
    void NonlinearSWE::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {     
    ASSERTL0(false,"1D not working with user defined BC");
  }


    void NonlinearSWE::GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    
    NekDouble g = m_g;
    int nq = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> tmp(nq);
    
    switch(i){
      
      // flux function for the h equation
    case 0:
     {
       // h in tmp
       Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, tmp, 1);
       
       // hu in flux 0
       Vmath::Vmul(nq, tmp, 1, physfield[1], 1, flux[0], 1);
     }
     break;
     
     // flux function for the hu equation
    case 1:
      {
	// h in tmp
	Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, tmp, 1);
	    
	// hu in flux 0
	Vmath::Vmul(nq, tmp, 1, physfield[1], 1, flux[0], 1);
	
	// huu in flux 0
	Vmath::Vmul(nq, flux[0], 1, physfield[1], 1, flux[0], 1);
	
	//  hh in tmp
	Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
	    
	// huu + 0.5 g hh in flux 0
	Blas::Daxpy(nq, 0.5*g, tmp, 1, flux[0], 1);
      }
      break;
    default:
      ASSERTL0(false,"GetFluxVector1D: illegal vector index");
    }
  }
  
  
  void NonlinearSWE::GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    int nq = m_fields[0]->GetTotPoints();
    
    NekDouble g = m_g;
    
    switch(i){
      
      //----------------------------------------------------
      // flux function for the h equation
      
    case 0:
      {
	// h in flux 1
	Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, flux[1], 1);
	
	// hu in flux 0
	Vmath::Vmul(nq, flux[1], 1, physfield[1], 1, flux[0], 1);
	
	// hv in flux 1
	Vmath::Vmul(nq, flux[1], 1, physfield[2], 1, flux[1], 1);
      }
      break;
      //----------------------------------------------------
      
      
      //----------------------------------------------------
      // flux function for the hu equation
      
    case 1:
      {
	Array<OneD, NekDouble> tmp(nq);
	
	// h in tmp
	Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, tmp, 1);
	
	// hu in flux 1
	Vmath::Vmul(nq, tmp, 1, physfield[1], 1, flux[1], 1);
	
	// huu in flux 0
	Vmath::Vmul(nq, flux[1], 1, physfield[1], 1, flux[0], 1);
	
	//  hh in tmp
	Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
	
	// huu + 0.5 g hh in flux 0
	// Daxpy overwrites flux[0] on exit
	Blas::Daxpy(nq, 0.5*g, tmp, 1, flux[0], 1);
	
	// huv in flux 1
	Vmath::Vmul(nq, flux[1], 1, physfield[2], 1, flux[1], 1);
      }
      break;
      //----------------------------------------------------
      
      
      //----------------------------------------------------
      // flux function for the hv equation
      
    case 2:
      {
	Array<OneD, NekDouble> tmp(nq);
	
	// h in tmp
	Vmath::Vadd(nq, physfield[0], 1, m_depth, 1, tmp, 1);
	
	// hv in flux 0
	Vmath::Vmul(nq, tmp, 1, physfield[2], 1, flux[0], 1);
	
	// hvv in flux 1
	Vmath::Vmul(nq, flux[0], 1, physfield[2], 1, flux[1], 1);
	
	//  hh in tmp
	Vmath::Vmul(nq, tmp, 1, tmp, 1, tmp, 1);
	
	// huv in flux 0
	Vmath::Vmul(nq, flux[0], 1, physfield[1], 1, flux[0], 1);
	
	// hvv + 0.5 g hh in flux 1
	Blas::Daxpy(nq, 0.5*g, tmp, 1, flux[1], 1);
      }
      break;
      //----------------------------------------------------
      
    default:
      ASSERTL0(false,"GetFluxVector2D: illegal vector index");
    }
  }
  
  
  void NonlinearSWE::NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxX)
    {
      ASSERTL0(false,"1D DG not working");
    }
  
  void NonlinearSWE::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; // only the dependent variables 
    
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
    Array<OneD, NekDouble>  DepthFwd(nTraceNumPoints);
    Array<OneD, NekDouble>  DepthBwd(nTraceNumPoints);
    
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
	m_fields[0]->GetFwdBwdTracePhys(m_depth,DepthFwd,DepthBwd);
	
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
	    // note that we are using the same depth - i.e. the depth is assumed continuous...

	    switch(m_upwindType)
	      {
	      case eHLLC:
		{
		  RiemannSolverHLLC(Fwd[0][i]+DepthFwd[i],Fwd[1][i],Fwd[2][i],
				    Bwd[0][i]+DepthFwd[i],Bwd[1][i],Bwd[2][i],
				    hflux, huflux, hvflux );
		}
		break;
		
	      default:
		{
		  ASSERTL0(false,"populate switch statement for upwind flux");
		}
		break;
	      }
	    
            // rotate back to Cartesian
            numfluxX[0][i]  = hflux*m_traceNormals[0][i];
            numfluxY[0][i]  = hflux*m_traceNormals[1][i];
            numfluxX[1][i] = (huflux*m_traceNormals[0][i] - hvflux*m_traceNormals[1][i]) * m_traceNormals[0][i];
            numfluxY[1][i] = (huflux*m_traceNormals[0][i] - hvflux*m_traceNormals[1][i]) * m_traceNormals[1][i];
            numfluxX[2][i] = (huflux*m_traceNormals[1][i] + hvflux*m_traceNormals[0][i]) * m_traceNormals[0][i];
            numfluxY[2][i] = (huflux*m_traceNormals[1][i] + hvflux*m_traceNormals[0][i]) * m_traceNormals[1][i];
	  }
  }
  
  void NonlinearSWE::RiemannSolverHLLC(NekDouble hL,NekDouble uL,NekDouble vL,NekDouble hR,NekDouble uR, 
					    NekDouble vR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux )
  {
    NekDouble g = m_g;
    NekDouble hC,huC,hvC,SL,SR,hstar,Sstar;
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
  
    void NonlinearSWE::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					      Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    int nq = GetTotPoints();
      
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(3);
	for (int i = 0; i < 3; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	
	// \eta = h - d
	Vmath::Vsub(nq,tmp[0],1,m_depth,1,physout[0],1);
	
	// u = hu/h
	Vmath::Vdiv(nq,tmp[1],1,tmp[0],1,physout[1],1);
	
	// v = hv/ v
	Vmath::Vdiv(nq,tmp[2],1,tmp[0],1,physout[2],1);
      }
    else
      {
	// \eta = h - d
	Vmath::Vsub(nq,physin[0],1,m_depth,1,physout[0],1);
	
	// u = hu/h
	Vmath::Vdiv(nq,physin[1],1,physin[0],1,physout[1],1);
	
	// v = hv/ v
	Vmath::Vdiv(nq,physin[2],1,physin[0],1,physout[2],1);
      }
  }


   void NonlinearSWE::v_ConservativeToPrimitive( )
  {
    int nq = GetTotPoints();
    
    // u = hu/h
    Vmath::Vdiv(nq,m_fields[1]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);
	
    // v = hv/ v
    Vmath::Vdiv(nq,m_fields[2]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);

    // \eta = h - d
    Vmath::Vsub(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);
  }

  void NonlinearSWE::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					     Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    
    int nq = GetTotPoints();
    
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(3);
	for (int i = 0; i < 3; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	
	// h = \eta + d
	Vmath::Vadd(nq,tmp[0],1,m_depth,1,physout[0],1);
	
	// hu = h * u
	Vmath::Vmul(nq,physout[0],1,tmp[1],1,physout[1],1);
	
	// hv = h * v
	Vmath::Vmul(nq,physout[0],1,tmp[2],1,physout[2],1);
      
      }
    else
      {
	// h = \eta + d
	Vmath::Vadd(nq,physin[0],1,m_depth,1,physout[0],1);
	
	// hu = h * u
	Vmath::Vmul(nq,physout[0],1,physin[1],1,physout[1],1);
	
	// hv = h * v
	Vmath::Vmul(nq,physout[0],1,physin[2],1,physout[2],1);
	
      }
     
  }

  void NonlinearSWE::v_PrimitiveToConservative( )
  {
    int nq = GetTotPoints();
    
    // h = \eta + d
    Vmath::Vadd(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);
    
    // hu = h * u
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[1]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);
    
    // hv = h * v
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[2]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);
  }

} //end of namespace

