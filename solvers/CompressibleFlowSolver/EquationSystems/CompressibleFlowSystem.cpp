///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <MultiRegions/LocalToGlobalDGMap.h>
#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

namespace Nektar
{

  string CompressibleFlowSystem::className = GetEquationSystemFactory().RegisterCreatorFunction("CompressibleFlowSystem", CompressibleFlowSystem::create, "Auxiliary functions for the compressible flow system.");
  
  CompressibleFlowSystem::CompressibleFlowSystem(
          LibUtilities::CommSharedPtr& pComm,
          LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pComm, pSession)
  {

  }

  void CompressibleFlowSystem::v_InitObject()
  {
      UnsteadySystem::v_InitObject();
  }

  CompressibleFlowSystem::~CompressibleFlowSystem()
  {
    
  }

  void CompressibleFlowSystem::v_PrintSummary(std::ostream &out)
  {
    UnsteadySystem::v_PrintSummary(out);
  }

  //----------------------------------------------------
 
  void CompressibleFlowSystem::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
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
	      ASSERTL0(false,"1D not yet implemented for the Compressible Flow Equations");	
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
	    ASSERTL0(false,"3D not implemented for Compressible Flow Equations");
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
  
  void CompressibleFlowSystem::SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {  
    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = physarray.num_elements();
    
    // get physical values of the forward trace (from exp to phys)
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
      }
    
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
	      ASSERTL0(false,"1D not yet implemented for the Compressible Flow Equations");
	    }
	    break;
	  case 2:
	    {
	      Array<OneD, NekDouble> tmp_t(npts);
	      
	      Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
	      Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	      
	      Array<OneD, NekDouble> tmp_n(npts,0.0);

	      // rotate back to Cartesian
	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
	      Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);
	      
	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
	      Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
	    }
	    break;
	  case 3:
	    ASSERTL0(false,"3D not implemented for the Compressible Flow Equations");
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

  void CompressibleFlowSystem::GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
				  Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    int nq = m_fields[0]->GetTotPoints();

    NekDouble gamma = m_gamma;
    Array<OneD, NekDouble> pressure(nq);
    
    switch(i){
      
      // flux function for the rho equation
    case 0:
      {
	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
 	    flux[0][j]  =  physfield[1][j];
 	    flux[1][j]  =  physfield[2][j];
 	  }
      }
      break;
      
      // flux function for the rhou equation
    case 1:
      {
	GetPressure(physfield,pressure);


 	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
 	  {
 	    flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] + pressure[j];
 	    flux[1][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
 	  }
      }
      break;
      
      // flux function for the rhov equation
    case 2:
      {
	GetPressure(physfield,pressure);

 	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
 	  {
 	    flux[0][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
 	    flux[1][j] = physfield[2][j]*physfield[2][j]/physfield[0][j] + pressure[j];
 	  }
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
      ASSERTL0(false,"GetFluxVector2D: illegal vector index");
    }
  }

  void CompressibleFlowSystem::GetPressure(const Array<OneD, const Array<OneD, NekDouble> > &physfield,
			     Array<OneD, NekDouble> &pressure)
  {
    NekDouble gamma = m_gamma;
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
 	pressure[i] = (gamma - 1.0)*(physfield[3][i] - 0.5*(physfield[1][i]*physfield[1][i]/physfield[0][i] +
 							    physfield[2][i]*physfield[2][i]/physfield[0][i]));
      }
  }
  
  void CompressibleFlowSystem::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
				 Array<OneD, Array<OneD, NekDouble> > &numfluxY)
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
  
  void CompressibleFlowSystem::RiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
			       NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
			       NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux)
  {
    switch(m_upwindType)
      {
      case eNotSet:
	{
	  ASSERTL0(false,"No upwind flux set in the input file");
	}
	break;
      case eAverage:
	{
	  AverageRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
      case eLLF:
	{
	  LFRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      case eHLL:
	{
	  HLLRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      case eHLLC:
	{
	  HLLCRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      case eExact:
	{
	  ExactRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      case eAUSM:
      case eAUSMPlus:
      case eAUSMPlusUp:
      case eAUSMPlusUpAllSpeed:
	{
	  AUSMRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      default:
	{
	  ASSERTL0(false,"populate switch statement for upwind flux");
	}
	break;
      }
  }
  
  
  void CompressibleFlowSystem::HLLRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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
    NekDouble SR = max(uR+cR, uRoe+cRoe);
    
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

  void CompressibleFlowSystem::HLLCRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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
    NekDouble SR = max(uR+cR, uRoe+cRoe);
    NekDouble SM = (pR-pL+rhouL*(SL-uL)-rhouR*(SR-uR))/(rhoL*(SL-uL)-rhoR*(SR-uR));

    // compute the HLLC flux
    if (SL >= 0.0)
      {
	rhoflux  = rhoL*uL;
	rhouflux = rhoL*uL*uL + pL;
	rhovflux = rhoL*uL*vL;
	Eflux    = uL*(EL + pL);
      }
    else if (SR <= 0.0)
      {
	rhoflux  = rhoR*uR;
	rhouflux = rhoR*uR*uR + pR;
	rhovflux = rhoR*uR*vR;
	Eflux    = uR*(ER + pR);
      }
    else
      {
	
	NekDouble rhoML = rhoL*(SL-uL)/(SL-SM);
	NekDouble rhouML = rhoML*SM;
	NekDouble rhovML = rhoML*vL;
	NekDouble EML = rhoML*(EL/rhoL+(SM-uL)*(SM+pL/(rhoL*(SL-uL))));
	
	NekDouble rhoMR = rhoR*(SR-uR)/(SR-SM);
	NekDouble rhouMR = rhoMR*SM;
	NekDouble rhovMR = rhoMR*vR;
	NekDouble EMR = rhoMR*(ER/rhoR+(SM-uR)*(SM+pR/(rhoR*(SR-uR))));

	if (SL < 0.0 && SM >= 0.0)
	  {
	    rhoflux  = rhoL*uL         + SL*(rhoML - rhoL);
	    rhouflux = rhoL*uL*uL + pL + SL*(rhouML - rhouL);
	    rhovflux = rhoL*uL*vL      + SL*(rhovML - rhovL);
	    Eflux    = uL*(EL + pL)    + SL*(EML - EL);
	  }
	else if(SM < 0.0 && SR > 0.0)
	  {
	    rhoflux  = rhoR*uR         + SR*(rhoMR - rhoR);
	    rhouflux = rhoR*uR*uR + pR + SR*(rhouMR - rhouR);
	    rhovflux = rhoR*uR*vR      + SR*(rhovMR - rhovR);
	    Eflux    = uR*(ER + pR)    + SR*(EMR - ER);
	  }

      }
    
  }


  void CompressibleFlowSystem::LFRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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
    
    // compute the wave speeds
    NekDouble S = max(uR+cR, -uL+cL);
    NekDouble sign = 1;
    if(S == -uL+cL)
      sign = -1;

    // compute the Lax-Friedrichs flux
    rhoflux  = 0.5*((rhouL+rhouR)                 - sign*S*(rhoR -rhoL));
    rhouflux = 0.5*((rhoL*uL*uL+pL+rhoR*uR*uR+pR) - sign*S*(rhouR-rhouL));
    rhovflux = 0.5*((rhoL*uL*vL+rhoR*uR*vR)       - sign*S*(rhovR-rhovL));
    Eflux    = 0.5*((uL*(EL+pL)+uR*(ER+pR))       - sign*S*(ER   -EL))   ;
    
  }

  void CompressibleFlowSystem::AverageRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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

    // compute the Average flux
    rhoflux  = 0.5*(rhouL+rhouR);
    rhouflux = 0.5*(rhoL*uL*uL+pL+rhoR*uR*uR+pR);
    rhovflux = 0.5*(rhoL*uL*vL+rhoR*uR*vR);
    Eflux    = 0.5*(uL*(EL+pL)+uR*(ER+pR));
    
  }

  void CompressibleFlowSystem::AUSMRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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
    
    // average sound speed
    NekDouble cA = 0.5 * (cL + cR);
    
    // local Mach numbers
    NekDouble ML = uL/cA;
    NekDouble MR = uR/cA;

    // parameters that specifies the upwinding...
    NekDouble Mbar, pbar;

    switch(m_upwindType)
      {
      case eAUSM:
	{
	  NekDouble beta  = 0.0;
	  NekDouble alpha = 0.0;
	  Mbar = M4Function(0, beta, ML) + M4Function(1, beta, MR);
	  pbar = pL*P5Function(0, alpha, ML) + pR*P5Function(1, alpha, MR);
	}
	break;
      case eAUSMPlus:
	{
	  NekDouble beta  = 1.0/8.0;
	  NekDouble alpha = 3.0/16.0;
	  Mbar = M4Function(0, beta, ML) + M4Function(1, beta, MR);
	  pbar = pL*P5Function(0, alpha, ML) + pR*P5Function(1, alpha, MR);
	}
	break;
      case eAUSMPlusUp:
	{
	  NekDouble beta  = 1.0/8.0;
	  NekDouble alpha = 3.0/16.0;
	  NekDouble sigma = 1.0;
	  NekDouble Kp    = 0.25;
	  NekDouble Ku    = 0.75;
	  NekDouble Mtilde = 0.5 * (ML*ML + MR*MR);
	  NekDouble rhoA   = 0.5 * (rhoL + rhoR);
	  NekDouble Mp     = -Kp * ((pR-pL)/(rhoA*cA*cA)) * max(1.0 - sigma * Mtilde, 0.0);
	  Mbar =  M4Function(0, beta, ML) + M4Function(1, beta, MR) + Mp;
	  NekDouble pu     = -2.0*Ku*rhoA*cA*cA*(MR-ML)*P5Function(0, alpha, ML)*P5Function(1, alpha, MR);
	  pbar = pL*P5Function(0, alpha, ML) + pR*P5Function(1, alpha, MR) + pu;
	}
	break;
      case eAUSMPlusUpAllSpeed:
	{
	  // if fa = 1 then AUSMPlusUpAllSpeed equal to AUSMPlusUp...
	  NekDouble Mco    = 0.01;
	  NekDouble Mtilde = 0.5 * (ML*ML + MR*MR);
	  NekDouble Mo     = min(1.0, max(Mtilde,Mco*Mco));
	  NekDouble fa    = Mo*(2.0-Mo);
	  NekDouble beta  = 1.0/8.0;
	  NekDouble alpha = 3.0/16.0;
	  NekDouble sigma = 1.0;
	  NekDouble Kp    = 0.25;
	  NekDouble Ku    = 0.75;
	  NekDouble rhoA   = 0.5 * (rhoL + rhoR);
	  NekDouble Mp     = -(Kp/fa) * ((pR-pL)/(rhoA*cA*cA)) * max(1.0 - sigma * Mtilde, 0.0);
	  Mbar =  M4Function(0, beta, ML) + M4Function(1, beta, MR) + Mp;
	  NekDouble pu     = -2.0*Ku*rhoA*cA*cA*(MR-ML)*P5Function(0, alpha, ML)*P5Function(1, alpha, MR);
	  pbar = pL*P5Function(0, alpha, ML) + pR*P5Function(1, alpha, MR) + pu;
	}
	break;
      default:
	{
	  ASSERTL0(false,"AUSM family chosen but specific scheme not implemented");
	}
      }
    
    if (Mbar >= 0.0)
      {
	rhoflux  = cA * Mbar * rhoL;
	rhouflux = cA * Mbar * rhoL * uL + pbar;
	rhovflux = cA * Mbar * rhoL * vL;
	Eflux    = cA * Mbar * (EL + pL);
      }
    else
      {
	rhoflux  = cA * Mbar * rhoR;
	rhouflux = cA * Mbar * rhoR * uR + pbar;
	rhovflux = cA * Mbar * rhoR * vR;
	Eflux    = cA * Mbar * (ER + pR);
      }
  }

  NekDouble CompressibleFlowSystem::M1Function(int A, NekDouble M)
  {
    NekDouble out;
    
    if (A == 0) // i.e. "plus"
      {
	out = 0.5 * (M + fabs(M));
      }
    else 
      {
	out = 0.5 * (M - fabs(M));
      }
    
    return out; 
  }

  NekDouble CompressibleFlowSystem::M2Function(int A, NekDouble M)
  {
    NekDouble out;
    
    if (A == 0) // i.e. "plus"
      {
	out = 0.25 * (M + 1.0) * (M + 1.0);
      }
    else
      {
	out = -0.25 * (M - 1.0) * (M - 1.0);
      }
    
    return out;
  }

  NekDouble CompressibleFlowSystem::M4Function(int A, NekDouble beta, NekDouble M)
  {
    NekDouble out;
    
    if (fabs(M) >= 1.0)
      {
	out = M1Function(A,M);
      }
    else
      {
	out = M2Function(A,M);
	
	if (A == 0)
	  {
	    out *= 1.0 - 16.0*beta*M2Function(1,M);
	  }
	else
	  {
	    out *= 1.0 + 16.0*beta*M2Function(0,M);
	  }
      }

    return out;
  }

  NekDouble CompressibleFlowSystem::P5Function(int A, NekDouble alpha, NekDouble M)
  {
    NekDouble out;

    if (fabs(M) >= 1.0)
      {
	out = (1.0/M) * M1Function(A,M);
      }
    else
      {
	out = M2Function(A,M);
	
	if (A == 0)
	  {
	    out *= (2.0-M) - 16.0*alpha*M*M2Function(1,M);
	  }
	else
	  {
	    out *= (-2.0-M) + 16.0*alpha*M*M2Function(0,M);
	  }
      }
    
    return out;
  }

  void CompressibleFlowSystem::ExactRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
				    NekDouble rhoR, NekDouble rhouR, NekDouble rhovR, NekDouble ER,
				    NekDouble &rhoflux, NekDouble &rhouflux, NekDouble &rhovflux, NekDouble &Eflux)
  {

    // Exact Riemann Solver (GOTTLIEB AND GROTH - 1987; TORO - 1998) 

    NekDouble gamma = m_gamma;
    // gamma operation
    NekDouble f1=(gamma-1.0)/(2.0*gamma);
    NekDouble f2=2.0/(gamma-1.0);
    NekDouble f3=(gamma+1.0)/(2.0*gamma);
    NekDouble f4=(gamma+1.0)/4.0;

    // Right and left variables
    NekDouble uL = rhouL/rhoL;
    NekDouble vL = rhovL/rhoL;
    NekDouble uR = rhouR/rhoR;
    NekDouble vR = rhovR/rhoR;
    NekDouble pL = (gamma - 1.0) * (EL - 0.5 * (rhouL*uL + rhovL*vL));
    NekDouble pR = (gamma - 1.0) * (ER - 0.5 * (rhouR*uR + rhovR*vR));
    NekDouble cL = sqrt(gamma * pL / rhoL);
    NekDouble cR = sqrt(gamma * pR / rhoR);

    // Five possible configuration: rcs=1,scr=2,scs=3,rcr=4,rcvcr=5
    // r = rarefaction
    // c = contact discontinuity
    // s = shock
    // v = vacuum

    NekDouble pratio,z;
    NekDouble u_ncr,u_rcn,u_scn,u_ncs;
    NekDouble u_rcvr = uL +(f2*cL)+(f2*cR);
    int pattern = -1;
    
    if(uR<u_rcvr)
      {
	if(uR>=uL)
	  {
	    if(pR>=pL)
	      { 
		pratio=pL/pR;
		u_ncr=uL+(f2*cR)*(1.0-pow(pratio,f1));
		if(uR>=u_ncr)
		  pattern=4;
		else
		  pattern=2;
	      }
	    else if(pR<pL)
	      {
		pratio=pR/pL;
		u_rcn=uL+(f2*cL)*(1.0-pow(pratio,f1));
		if(uR>=u_rcn)
		  pattern=4;
		else
		  pattern=1;
	      }
	  }
	else if(uR<uL)
	  {
	    if(pR>=pL)
	      {
		pratio=pR/pL;
		u_scn=uL-((cL/gamma)*(pratio-1.0)/sqrt((f3*pratio)-f1));
		if(uR>=u_scn)
		  pattern=2;
		else
		  pattern=3;
	      }
	    else if(pR<pL)
	      {
		pratio=pL/pR;
		u_ncs=uL-((cR/gamma)*(pratio-1.0)/sqrt((f3*pratio)-f1));
		if(uR>=u_ncs)
		  pattern=1;
		else
		  pattern=3;
	      }
	  }
      }
    else if(uR>=u_rcvr)
      pattern=5;
    
    // Initial Guess for u_int
    NekDouble aL = (gamma*pL)/cL;
    NekDouble aR = (gamma*pR)/cR;
    if(pL>=pR)
      z = (f2/f2)*(cR/cL)*pow((pL/pR),f1);
    else if(pL<pR)
      z = (f2/f2)*(cR/cL)*pow((pL/pR),f1);
    NekDouble u_int = ((z*(uL+(f2*cL)))+(uR-f2*cR))/(1.0+z); 

    /*   PATTERN RAREFACTION WAVE - CONTACT SURFACE - SHOCK  */
    
    NekDouble c_intR,c_intL,derp_intR,derp_intL,wR,wL,p_int,p_intL,p_intR,u_intL,u_intR;
    NekDouble swave1,swave2,swave3,swave4,swave5; 
    NekDouble u0,u1,u2,u3;
    NekDouble uaux,aaux,paux;
    NekDouble chi = 0.0;
    NekDouble EPSILON = 1.0e-6;
    
    if(pattern==1)
      {
	p_intR = 1.0;
	p_intL = p_intR*(1.0+10.0*EPSILON);
	while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
	  { 
	    c_intL    = cL-((u_int-uL)/f2);
	    p_intL    = pL*pow((c_intL/cL),(1.0/f1));
	    derp_intL = (-gamma*p_intL)/c_intL;
	    wR        = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
	    p_intR    = pR+(aR*wR*(u_int-uR));
	    derp_intR = (2.0*aR*pow(wR,3.0))/(1.0+(wR*wR));
	    u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
	  }
	p_int  = (p_intL+p_intR)/2.0;
	c_intR = cR*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pR))/(gamma+1.0+(gamma-1.0)*(pR/p_int)));
	c_intL = cL-((u_int-uL)/f2);
	wR     = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
	swave1 = (wR*cR)+uR;
	swave3 = u_int;
	swave4 = u_int-c_intL;
	swave5 = uL-cL;
	if(chi>=swave1)
	  {
	    u0 = gamma*(pR/(cR*cR));
	    u1 = u0*uR;
	    u2 = u0*vR;
	    u3 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR + vR*vR)/(cR*cR)));
	  }
	else if((chi<swave1)&&(chi>=swave3))
	  {
	    u0 = gamma*(p_int/(c_intR*c_intR));
	    u1 = u0*u_int;
	    u2 = u0*vR;
	    u3 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int + vR*vR)/(c_intR*c_intR)));
	  }
	else if((chi<swave3)&&(chi>=swave4))
	  {
	    u0 = gamma*(p_int/(c_intL*c_intL));
	    u1 = u0*u_int;
	    u2 = u0*vL;
	    u3 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int + vL*vL)/(c_intL*c_intL)));
	  }
	else if((chi<swave4)&&(chi>=swave5))
	  {
	    uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
	    aaux = cL+((1.0/f2)*(uL-uaux));
	    paux = pL*pow((aaux/cL),(1.0/f1));
	    u0   = (gamma*paux)/(aaux*aaux);
	    u1   = u0*uaux;
	    u2   = u0*vL;
	    u3   = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL)*0.5);        
	  }
	else if(chi<swave5)
	  {
	    u0 = gamma*(pL/(cL*cL));
	    u1 = u0*uL;
	    u2 = u0*vL;
	    u3 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL)/(cL*cL)));  
	  }    
      }
    
    /*        PATTERN:   SHOCK - CONTACT SURFACE - RAREFACTION WAVE    */
    
    else if(pattern==2)
      {
	p_intR = 1.0;
	p_intL = p_intR*(1.0+10.0*EPSILON);
	while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
	  {    
	    wL        = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
	    p_intL    = pL+(aL*wL*(u_int-uL));
	    derp_intL = (2.0*aL*pow(wL,3.0))/(1.0+(wL*wL));
	    c_intR    = cR+((u_int-uR)/f2);
	    p_intR    = pR*pow((c_intR/cR),(1.0/f1));
	    derp_intR = (gamma*p_intR)/c_intR;
	    u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
	  }
	p_int  = (p_intL+p_intR)/2.0;
	c_intL = cL*sqrt((gamma+1.0+((gamma-1.0)*(p_int/pL)))/(gamma+1.0+((gamma-1.0)*(pL/p_int))));
	wL     = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
	swave1 = uR+cR;
	swave2 = u_int+c_intR;
	swave3 = u_int;
	swave5 = (wL*cL)+uL;
	if(chi>=swave1)
	  {
	    u0 = gamma*(pR/(cR*cR));
	    u1 = u0*uR;
	    u2 = u0*vR;
	    u3 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR)/(cR*cR)));
	  }
	else if((chi<swave1)&&(chi>=swave2))
	  { 
	    uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
	    aaux = cR+((1.0/f2)*(uaux-uR));
	    paux = pR*pow((aaux/cR),(1.0/f1));
	    u0 = (gamma*paux)/(aaux*aaux);
	    u1 = u0*uaux;
	    u2 = u0*vR;
	    u3 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR)*0.5);  
	  }
	else if((chi<swave2)&&(chi>=swave3))
	  { 
	    u0 = gamma*(p_int/(c_intR*c_intR));
	    u1 = u0*u_int;
	    u2 = u0*vR;
	    u3 = (p_int/(gamma-1.0))+(u0*(u_int*u_int+vR*vR)*0.5);
	  }
	else if((chi<swave3)&&(chi>=swave5))
	  { 
	    u0 = gamma*(p_int/(c_intL*c_intL));
	    u1 = u0*u_int;
	    u2 = u0*vL;
	    u3 = (p_int/(gamma-1.0))+(u0*(u_int*u_int+vL*vL)*0.5);
	  }
	else if(chi<swave5)
	  { 
	    u0 = gamma*(pL/(cL*cL));
	    u1 = u0*uL;
	    u2 = u0*vL;
	    u3 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL)/(cL*cL)));
	  }
      }

    /*       PATTERN: SHOCK - CONTACT SURFACE -SHOCK       */
    
    else if(pattern==3)
      {
	p_intR = 1.0;
	p_intL = p_intR*(1.0+10.0*EPSILON);
	while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
	  {      
	    wL        = f4*((u_int-uL)/cL)-sqrt(1.0+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
	    p_intL    = pL+(aL*wL*(u_int-uL));
	    derp_intL = (2.0*aL*pow(wL,3.0))/(1.0+(wL*wL));
	    wR        = f4*((u_int-uR)/cR)+sqrt(1.0+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
	    p_intR    = pR+(aR*wR*(u_int-uR));
	    derp_intR = (2.0*aR*pow(wR,3.0))/(1.0+(wR*wR));
	    u_int     = u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
	  }
	p_int  = (p_intL+p_intR)/2.0;
	c_intL = cL*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pL))/(gamma+1.0+(gamma-1.0)*(pL/p_int)));
	c_intR = cR*sqrt((gamma+1.0+(gamma-1.0)*(p_int/pR))/(gamma+1.0+(gamma-1.0)*(pR/p_int)));
	wR     = f4*((u_int-uR)/cR)+sqrt(1+(f4*((u_int-uR)/cR))*(f4*((u_int-uR)/cR)));
	wL     = f4*((u_int-uL)/cL)-sqrt(1+(f4*((u_int-uL)/cL))*(f4*((u_int-uL)/cL)));
	swave1 = (wR*cR)+uR;
	swave3 = u_int;
	swave5 = (wL*cL)+uL;
	if(chi>=swave1)
	  {
	    u0 = gamma*(pR/(cR*cR));
	    u1 = u0*uR;
	    u2 = u0*vR;
	    u3 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR)/(cR*cR)));
	  }
	else if((chi<swave1)&&(chi>=swave3))
	  {
	    u0 = gamma*(p_int/(c_intR*c_intR));
	    u1 = u0*u_int;
	    u2 = u0*vR;
	    u3 = (p_int/(gamma-1.0))+(0.5*u0*(u_int*u_int+vR*vR));
	  }
	else if((chi<swave3)&&(chi>=swave5))
	  {
	    u0 = gamma*(p_int/(c_intL*c_intL));
	    u1 = u0*u_int;
	    u2 = u0*vL;
	    u3 = (p_int/(gamma-1.0))+(0.5*u0*(u_int*u_int+vL*vL));
	  }
	else if(chi<swave5)
	  {
	    u0 = gamma*(pL/(cL*cL));
	    u1 = u0*uL;
	    u2 = u0*vL;
	    u3 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL)/(cL*cL+vL*vL)));
	  }
      }

    /*    PATTERN:  RAREFACTION WAVE - CONTAC SURFACE - RAREFACTION WAVE   */
    
    else if(pattern==4)
      { 
	p_intR=1.0;
	p_intL=p_intR*(1.0+10.0*EPSILON);
	while(fabs(1.0-(p_intL/p_intR))>=EPSILON)
	  {
	    c_intL=cL-((u_int-uL)/f2);
	    p_intL=pL*pow((c_intL/cL),(1.0/f1));
	    derp_intL=(-gamma*p_intL)/c_intL;
	    c_intR=cR+((u_int-uR)/f2);
	    p_intR=pR*pow((c_intR/cR),(1.0/f1));
	    derp_intR=(gamma*p_intR)/c_intR;
	    u_int=u_int-((p_intL-p_intR)/(derp_intL-derp_intR));
	  }
	p_int=(p_intL+p_intR)/2.0;
	c_intL=cL-((u_int-uL)/f2);
	c_intR=cR+((u_int-uR)/f2);
	swave1=uR+cR;
	swave2=u_int+c_intR;
	swave3=u_int;
	swave4=u_int-c_intL;
	swave5=uL-cL;
	  
	  if(chi>=swave1)
	    {
	      u0 = gamma*(pR/(cR*cR));
	      u1 = u0*uR;
	      u2 = u0*vR;
	      u3 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR)/(cR*cR)));
	    }
	  else if((chi<swave1)&&(chi>=swave2))
	    {
	      uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
	      aaux = cR+((1.0/f2)*(uaux-uR));
	      paux = pR*pow((aaux/cR),(1.0/f1));
	      u0 = (gamma*paux)/(aaux*aaux);
	      u1 = u0*uaux;
	      u2 = u0*vR;
	      u3 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR)*.5);  
	    }
	  else if((chi<swave2)&&(chi>=swave3))
	    {
	      u0 = gamma*(p_int/(c_intR*c_intR));
	      u1 = u0*u_int;
	      u2 = u0*vR;
	      u3 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int+vR*vR)/(c_intR*c_intR)));
	    }
	  else if((chi<swave3)&&(chi>=swave4))
	    {
	      u0 = gamma*(p_int/(c_intL*c_intL));
	      u1 = u0*u_int;
	      u2 = u0*vL;
	      u3 = (p_int/(gamma-1.0))+(gamma*p_int*0.5*((u_int*u_int+vL*vL)/(c_intL*c_intL)));       
	    }
	  else if((chi<swave4)&&(chi>=swave5))
	    {
	      uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
	      aaux = cL+((1.0/f2)*(uL-uaux));
	      paux = pL*pow((aaux/cL),(1.0/f1));
	      u0 = (gamma*paux)/(aaux*aaux);
	      u1 = u0*uaux;
	      u2 = u0*vL;
	      u3 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL)*.5);  
	    }
	  else if(chi<swave5)
	    {
	      u0 = gamma*(pL/(cL*cL));
	      u1 = u0*uL;
	      u2 = u0*vL;
	      u3 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL)/(cL*cL)));
	    }
      }
    
    /*PATTERN:   RAREFACTION WAVE - CONTACT SURFACE - VACUUM -
      CONTACT SURFACE - RAREFACTION WAVE    */
    else if(pattern==5)
      {
	p_int  = 0.0;
	u_intR = uR-(f2*cR);
	u_intL = uL+(f2*cL);
	c_intR = 0.0;
	c_intL = 0.0;
	swave1 = uR+cR;
	swave2 = u_intR+c_intR;
	swave4 = u_intL-c_intL;
	swave5 = uL-cL;
	if(chi>=swave1)
	  {
	    u0 = gamma*(pR/(cR*cR));
	    u1 = u0*uR;
	    u2 = u0*vR;
	    u3 = (pR/(gamma-1.0))+(gamma*pR*0.5*((uR*uR+vR*vR)/(cR*cR)));
	  }
	else if((chi<swave1)&&(chi>=swave2))
	  {
	    uaux = (2.0/(gamma+1.0))*(chi-cR+(uR/f2));
	    aaux = cR+((1.0/f2)*(uaux-uR));
	    paux = pR*pow((aaux/cR),(1.0/f1));
	    u0 = (gamma*paux)/(aaux*aaux);
	    u1 = u0*uaux;
	    u2 = u0*vR;
	    u3 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vR*vR)*.5);  
	  }
	else if((chi<swave2)&&(chi>=swave4))
	  {
	    u0 = 0.0;
	    u1 = 0.0;
	    u2 = 0.0;
	    u3 = 0.0;
	  }
	else if((chi<swave4)&&(chi>=swave5))
	  {
	    uaux = (2.0/(gamma+1.0))*(chi+cL+(uL/f2));
	    aaux = cL+((1.0/f2)*(uL-uaux));
	    paux = pL*pow((aaux/cL),(1.0/f1));
	    u0 = (gamma*paux)/(aaux*aaux);
	    u1 = u0*uaux;
	    u2 = u0*vL;
	    u3 = (paux/(gamma-1.0))+(u0*(uaux*uaux+vL*vL)*.5);  
	  }
	else if(chi<swave5)
	  {
	    u0 = gamma*(pL/(cL*cL));
	    u1 = u0*uL;
	    u2 = u0*vL;
	    u3 = (pL/(gamma-1.0))+(gamma*pL*0.5*((uL*uL+vL*vL)/(cL*cL)));
	  }
      }
    
    rhoflux  = u1;
    rhouflux = u1*u1/u0 + (gamma-1.0)*(u3 - 0.5*(u1*u1/u0+u2*u2/u0));
    rhovflux = u1*u2/u0;
    Eflux    = u1/u0 * (gamma*u3 - (gamma-1.0)*0.5*(u1*u1/u0+u2*u2/u0));
   
  }

  void CompressibleFlowSystem::GetVelocityVector(const Array<OneD, Array<OneD, NekDouble> > &physfield,
					 Array<OneD, Array<OneD, NekDouble> > &velocity)
  {
    NekDouble invRho;
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	invRho = 1/physfield[0][i];
	velocity[0][i] = invRho*physfield[1][i];
	velocity[1][i] = invRho*physfield[2][i];
      }
  }
  
  void CompressibleFlowSystem::GetTemperature(Array<OneD, Array<OneD, NekDouble> > &physfield,
				      Array<OneD, NekDouble> &pressure,
				      Array<OneD, NekDouble> &temperature)
  {
    NekDouble gamma = m_gamma;
    
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	temperature[i] = pressure[i]/(physfield[0][i]*m_GasConstant);
      }
  }
  
  void CompressibleFlowSystem::GetSoundSpeed(const Array<OneD, Array<OneD, NekDouble> > &physfield,
				     Array<OneD, NekDouble> &pressure,
				     Array<OneD, NekDouble> &soundspeed)
  {
    NekDouble gamma = m_gamma;
    
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	soundspeed[i] = sqrt(gamma*pressure[i]/physfield[0][i]);
      }
  }
  
  void CompressibleFlowSystem::GetMach(Array<OneD, Array<OneD, NekDouble> > &physfield,
			       Array<OneD, NekDouble> &soundspeed,
			       Array<OneD, NekDouble> &mach)
  {
    NekDouble velocity;
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	velocity = sqrt(physfield[1][i]/physfield[0][i]*physfield[1][i]/physfield[0][i]+physfield[2][i]/physfield[0][i]*physfield[2][i]/physfield[0][i]);
	mach[i] = velocity/soundspeed[i];
      }
  }
  
  NekDouble CompressibleFlowSystem::v_GetTimeStep(const Array<OneD, Array<OneD,NekDouble> > physarray, const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFL, NekDouble timeCFL)
  { 

    int nvariables = m_fields.num_elements();               // Number of variables in the mesh
    int nTotQuadPoints  = GetTotPoints();
    int n_element  = m_fields[0]->GetExpSize(); 
    Array<OneD, NekDouble> tstep(n_element,0.0);
    const NekDouble minLengthStdTri  = 0.7072*0.5;
    const NekDouble minLengthStdQuad = 0.5;
    const NekDouble cLambda = 0.2; // Spencer book pag. 317
    Array<OneD, NekDouble> stdVelocity(n_element,0.0);
    stdVelocity = GetStdVelocity(physarray);

    for(int el = 0; el < n_element; ++el)
      {
	int npoints = m_fields[0]->GetExp(el)->GetTotPoints();
	Array<OneD, NekDouble> one2D(npoints, 1.0);
	NekDouble Area = m_fields[0]->GetExp(el)->Integral(one2D);
	if(boost::dynamic_pointer_cast<LocalRegions::TriExp>(m_fields[0]->GetExp(el)))
	  {
	    //tstep[el] =  timeCFL*minLengthStdTri/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
	    tstep[el] = CFL[el]*minLengthStdTri/(stdVelocity[el]);
	  }
	else if(boost::dynamic_pointer_cast<LocalRegions::QuadExp>(m_fields[0]->GetExp(el)))
	  { 
	    //tstep[el] =  timeCFL*minLengthStdQuad/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
	    tstep[el] = CFL[el]*minLengthStdQuad/(stdVelocity[el]);
	  }
      }
    
    NekDouble TimeStep = Vmath::Vmin(n_element,tstep,1);

    return TimeStep;
  }

  
  Array<OneD,NekDouble> CompressibleFlowSystem::GetStdVelocity(const Array<OneD, Array<OneD,NekDouble> > inarray)
  {
    // Checking if the problem is 2D
    ASSERTL0(m_expdim==2,"Method not implemented for 1D and 3D");
    
    int nTotQuadPoints  = GetTotPoints();
    int n_element  = m_fields[0]->GetExpSize();// number of element in the mesh
    int npts = 0;

    // Getting the velocity vector on the 2D normal space
    Array<OneD, Array<OneD, NekDouble> > velocity(2);
    Array<OneD, Array<OneD, NekDouble> > stdVelocity(2);
    Array<OneD, NekDouble> pressure(nTotQuadPoints);
    Array<OneD, NekDouble> soundspeed(nTotQuadPoints);

    Array<OneD, NekDouble> stdV(n_element,0.0);
    for (int i = 0; i < 2; ++i)
      {
	velocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
	stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
      }
    GetVelocityVector(inarray,velocity);
    GetPressure(inarray,pressure);
    GetSoundSpeed(inarray,pressure,soundspeed);

    for(int el = 0; el < n_element; ++el)
      { 
	Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom2D()->GetJac();
	Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom2D()->GetGmat();

	int n_points = m_fields[0]->GetExp(el)->GetTotPoints();

	if(m_fields[0]->GetExp(el)->GetGeom2D()->GetGtype() == SpatialDomains::eDeformed)
	  {
	    // d xi/ dx = gmat = 1/J * d x/d xi
	    for(int i=0; i<n_points; i++)
	      {
		stdVelocity[0][i] = gmat[0][i]*velocity[0][i] + gmat[2][i]*velocity[1][i];
		stdVelocity[1][i] = gmat[1][i]*velocity[0][i] + gmat[3][i]*velocity[1][i];
	      }
	  }
	else
	  {
	    for(int i=0; i<n_points; i++)
	      {
		stdVelocity[0][i] = gmat[0][0]*velocity[0][i] + gmat[2][0]*velocity[1][i];
		stdVelocity[1][i] = gmat[1][0]*velocity[0][i] + gmat[3][0]*velocity[1][i];
	      }
	  }
	NekDouble pntVelocity;
	for(int i=0; i<n_points; i++)
	  {
	    pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] + stdVelocity[1][i]*stdVelocity[1][i]) + soundspeed[npts];
	    if(pntVelocity>stdV[el])
		stdV[el] = pntVelocity;
	    npts++;
	  }
      }

    return stdV;
  }
}
