	///////////////////////////////////////////////////////////////////////////////
//
// File APE.cpp
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
// Description: Acoustic perturbation equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <APESolver/EquationSystems/APE.h>

namespace Nektar
{
  string APE::className = GetEquationSystemFactory().RegisterCreatorFunction("APE", APE::create, "Acoustic perturbation equations in conservative variables.");
  
  APE::APE(
           const LibUtilities::SessionReaderSharedPtr& pSession)
    : APESystem(pSession)
  {
  }

  void APE::v_InitObject()
  {
	  APESystem::v_InitObject();
	  
    if (m_explicitAdvection)
      {

	m_ode.DefineOdeRhs     (&APE::DoOdeRhs,        this);
	m_ode.DefineProjection (&APE::DoOdeProjection, this);
      }
    else
      {

	ASSERTL0(false, "Implicit APE not set up.");
      }
  }
  
  APE::~APE()
  {
    
  }
 

  
  void APE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
	  // Add "source term"
	  // input: physical space
	  // output: modal space (JOSEF)
	  AddSource(physarray,modarray);
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
	      APE::GetFluxVector2D(i, physarray, fluxvector);
	      
	      m_fields[0]->PhysDeriv(0,fluxvector[0],tmp);
	      m_fields[0]->PhysDeriv(1,fluxvector[1],tmp1);
	      Vmath::Vadd(nq,tmp,1,tmp1,1,outarray[i],1);
	      Vmath::Neg(nq,outarray[i],1);
	    }
	  
	  //-------------------------------------------------
	  // Add "source term"
	  // input: physical space
	  // output: modal space
	  AddSource(physarray,outarray);
	  //------------------------------------------------- 
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the APE");
	break;
      }
  }
 

  void APE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
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
  void APE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time)
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
  void APE::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
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
	    ASSERTL0(false,"3D not implemented for Acoustic perturbation equations");
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
  
    void APE::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {     
    ASSERTL0(false,"1D not yet working for APE");
  }


    void APE::GetFluxVector1D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    ASSERTL0(false,"1D not yet working for APE");
  }
  
  
  void APE::GetFluxVector2D(const int i, const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    NekDouble Rho0 = m_Rho0;
    NekDouble gamma = m_gamma;

    
	basefield =Array<OneD, Array<OneD, NekDouble> >(m_spacedim+1);
	/*basefield[0]=U0
    basefield[1]=V0
    basefield[2]=P0

    physfield[0]=p'
    physfield[1]=u'
    physfield[2]=v'*/

    InitialiseBaseFlowAnalytical(basefield, m_time);
	

    switch(i){
      
      // flux function for the p' equation
    case 0:
      {
	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	  {
	  flux[0][j] = basefield[0][j]*physfield[0][j] + gamma*basefield[2][j]*physfield[1][j];
	  flux[1][j] = basefield[1][j]*physfield[0][j] + gamma*basefield[2][j]*physfield[2][j];
 	  }
      }
      break;
      
      // flux function for the u' equation
    case 1:
      {

 	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
 	  {
 	    flux[0][j] = basefield[0][j]*physfield[1][j] + basefield[1][j]*physfield[2][j] + physfield[0][j]/Rho0;
 	    flux[1][j] = 0;
 	  }
      }
      break;
      
      // flux function for the v' equation
    case 2:
      {

 	for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
 	  {
 	    flux[0][j] = 0;
 	    flux[1][j] = basefield[0][j]*physfield[1][j] + basefield[1][j]*physfield[2][j] + physfield[0][j]/Rho0;
 	  }
      }
      break;
    default:
      ASSERTL0(false,"GetFluxVector2D: illegal vector index");
    }
  }
  
  
  void APE::NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &numfluxX)
    {
      ASSERTL0(false,"1D DG not yet working for APE");
    }
  
   //Evaluation of the upwinded DG fluxes	
   void APE::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
				  Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {

    int i;
	//Number the points of the "shared" edges of the elements
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; //p', u', v' 
	  
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
	Array<OneD, Array<OneD, NekDouble> > rotbasefield(nvariables);
	Array<OneD, Array<OneD, NekDouble> > rotbasefieldBwd(nvariables);
	
	int nq = m_fields[0]->GetNpoints();
	Array<OneD,NekDouble> x0(nq);
	Array<OneD,NekDouble> x1(nq);
	Array<OneD,NekDouble> x2(nq);
	  
	m_fields[0]->GetCoords(x0,x1,x2);
	  
	  
	for (i = 0; i < nvariables; ++i)
    {
		Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		rotbasefield[i] = Array<OneD, NekDouble>(nTraceNumPoints);
		rotbasefieldBwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	}
    
	// get the physical values at the trace
	for (i = 0; i < nvariables; ++i)
	{
		m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
		m_fields[i]->GetFwdBwdTracePhys(basefield[i],rotbasefield[i],rotbasefieldBwd[i]);  
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
		  
		//rotates the baseflow
		tmpX =  rotbasefield[0][i]*m_traceNormals[0][i]+rotbasefield[1][i]*m_traceNormals[1][i];
		tmpY = -rotbasefield[0][i]*m_traceNormals[1][i]+rotbasefield[1][i]*m_traceNormals[0][i];
		rotbasefield[0][i] = rotbasefield[2][i];
		rotbasefield[1][i] = tmpX;
		rotbasefield[2][i] = tmpY;
	}
	  
	
        // Solve the Riemann problem
        NekDouble pflux, uflux, vflux;
	
        for (i = 0; i < nTraceNumPoints; ++i)
	  {
		  //cout << "Tracepoint: "<<i<<" of "<<nTraceNumPoints<<endl;

	    switch(m_upwindType)
	      {
	      case eUpwind:
		{
		  RiemannSolverUpwind(Fwd[0][i],Fwd[1][i],Fwd[2][i],
				    Bwd[0][i],Bwd[1][i],Bwd[2][i],rotbasefield[0][i],rotbasefield[1][i],rotbasefield[2][i],
				    pflux, uflux, vflux );
		}
		break;
		
	      default:
		{
		  ASSERTL0(false,"populate switch statement for upwind flux");
		}
		break;
	      }
	    
            // rotate back to Cartesian
			numfluxX[0][i]  = pflux*m_traceNormals[0][i];
			numfluxY[0][i]  = pflux*m_traceNormals[1][i];
            numfluxX[1][i] = (uflux*m_traceNormals[0][i] - vflux*m_traceNormals[1][i]) * m_traceNormals[0][i];
            numfluxY[1][i] = (uflux*m_traceNormals[0][i] - vflux*m_traceNormals[1][i]) * m_traceNormals[1][i];
            numfluxX[2][i] = (uflux*m_traceNormals[1][i] + vflux*m_traceNormals[0][i]) * m_traceNormals[0][i];
            numfluxY[2][i] = (uflux*m_traceNormals[1][i] + vflux*m_traceNormals[0][i]) * m_traceNormals[1][i];
	  }
	  
  }
  
  void APE::RiemannSolverUpwind(NekDouble pL,NekDouble uL,NekDouble vL,NekDouble pR,NekDouble uR, NekDouble vR, NekDouble P0, NekDouble U0, NekDouble V0,
					     NekDouble &pflux, NekDouble &uflux, NekDouble &vflux ) //(croth)
  {
	  int nvariables      = 2;
	  Array<OneD, NekDouble> characteristic(4);
	  Array<OneD, NekDouble> W(2);
	  Array<OneD, NekDouble> lambda(nvariables);
	  Array<OneD, NekDouble> upphysfield(3);
	  
	  // compute the wave speeds
	  lambda[0]=U0 + sqrt(P0*m_gamma*m_Rho0)/m_Rho0;
	  lambda[1]=U0 - sqrt(P0*m_gamma*m_Rho0)/m_Rho0;
	  
	  
	  // calculate the caracteristic variables 
	  //left characteristics
	  characteristic[0] = pL/2 + uL*sqrt(P0*m_gamma*m_Rho0)/2;
	  characteristic[1] = pL/2 - uL*sqrt(P0*m_gamma*m_Rho0)/2;
	  //right characteristics
	  characteristic[2] = pR/2 + uR*sqrt(P0*m_gamma*m_Rho0)/2;
	  characteristic[3] = pR/2 - uR*sqrt(P0*m_gamma*m_Rho0)/2;

	  //cout << endl;
	  for (int i=0; i<4; i++)
	  {
		  //cout << "characteristic["<<i<<"] = "<<characteristic[i]<<endl;
	  }
	  
	  
	  //take left or right value of characteristic variable
	  for (int j=0; j<nvariables; j++)
	  {
		  if (lambda[j]>=0)
		  {	  
			  W[j]=characteristic[j];
		  }
		  if(lambda[j]<0)
		  {
			  W[j]=characteristic[j+2];
		  }
	  }
	  
	  //cout << "lambda[0] = "<<lambda[0]<<"\tlambda[1] = "<<lambda[1]<<endl;

	  for (int i=0; i<2; i++)
	  {
		 // cout << "upwinded W["<<i<<"] = "<<W[i]<<endl;
	  }
	  
	  
	  //calculate conservative variables from characteristics
	  upphysfield[0]= W[0]+W[1];
	  upphysfield[1]= (W[0]-W[1])/sqrt(P0*m_gamma*m_Rho0);
	  upphysfield[2]= vL; 
		  
	  // compute the fluxes
	  pflux = U0*upphysfield[0] + m_gamma*P0*upphysfield[1];   	  
	  uflux = U0*upphysfield[1]+V0*upphysfield[2] + upphysfield[0]/m_Rho0;
	  vflux = 0.0;
	  
	  //cout << "pflux = "<<pflux<< "\tuflux = "<<uflux<< "\tvflux = "<<vflux<<endl;
	  
  }
  
    void APE::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					      Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    int i;
    int nq = GetTotPoints();
    int nvariables = physin.num_elements();
    
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(nvariables);
	for(i = 0; i < nvariables; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	for(i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(nq,tmp[i],1,physout[i],1);
	  }

      }
    else
      {
	for(i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(nq,physin[i],1,physout[i],1);
	  }
     
      }

  }


   void APE::v_ConservativeToPrimitive( )
  {
    
  }

  void APE::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					     Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    
    int i;
    int nq = GetTotPoints();
    int nvariables = physin.num_elements();
    
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(nvariables);
	for(i = 0; i < nvariables; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	for(i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(nq,tmp[i],1,physout[i],1);
	  }

      }
    else
      {
	for(i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(nq,physin[i],1,physout[i],1);
	  }
     
      }
     
  }

  void APE::v_PrimitiveToConservative( )
  {

  }

  // Initialise baseflow from the inputfile
  void APE::InitialiseBaseFlowAnalytical(Array<OneD, Array<OneD, NekDouble> > &base, const NekDouble time)
    {
        base = Array<OneD, Array<OneD, NekDouble> >(m_spacedim+1);
    	int nq = m_fields[0]->GetNpoints();
        std::string velStr[3] = {"U0","V0","P0"};        

  	        for(int i = 0; i <= m_spacedim; ++i)
        	{	
        		base[i] = Array<OneD, NekDouble> (nq,0.0);
        		EvaluateFunction(velStr[i], base[i], "Baseflow", time);
           	}
           
    }
  

	// Get sourceterm for p' equation from the inputfile
    void APE::GetSource(Array<OneD, NekDouble> &source, const NekDouble time)
    {
    	EvaluateFunction("S", source, "Source", time);
    }

	// Add sourceterm for p' equation obtained from GetSource
    void APE::AddSource(const Array< OneD, Array< OneD, NekDouble > > &inarray,
  		Array< OneD, Array< OneD, NekDouble > > &outarray)
    {
       
        int ncoeffs = outarray[0].num_elements();
        int nq      = inarray[0].num_elements();
    
        Array<OneD, NekDouble> source(nq);
    
        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {                   
                GetSource(source,m_time);
                
                //Source term solely for the p' equation (outarray[0])
                m_fields[0]->IProductWRTBase(source,source); 
                Vmath::Vadd(ncoeffs,source,1,outarray[0],1,outarray[0],1);
                
            }
            break;
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            {
                GetSource(source,m_time);
                
                //Source term solely for the p' equation (outarray[0])
                Vmath::Vadd(ncoeffs,source,1,outarray[0],1,outarray[0],1);	               }
            break;
            
        default:
            ASSERTL0(false,"Unknown projection scheme for the APE");
            break;
        }
    }

    void APE::v_PrintSummary(std::ostream &out)
    {
        EquationSystem::v_PrintSummary(out);
	out << "\tUpwind Type     : " << UpwindTypeMap[m_upwindType] << endl;
        out << "\tAdvection       : " << (m_explicitAdvection ? "explicit" : "implicit") << endl;
	out << "\tIntegration Type: " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
        out << "\tTime Step       : " << m_timestep << endl;
        out << "\tNo. of Steps    : " << m_steps << endl;
        out << "\tCheckpoints     : " << m_checksteps << " steps" << endl;

    }


} //end of namespace

