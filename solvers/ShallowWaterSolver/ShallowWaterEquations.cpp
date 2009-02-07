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
     
    //--------------------------------------------
    // Set conservative or primitive varibles
    if(m_boundaryConditions->CheckForParameter("Conservative") == true)
      {
	if((int) m_boundaryConditions->GetParameter("Conservative") == 0)
	  {
	    m_variableType = ePrimitive;
	  }
	else if ((int) m_boundaryConditions->GetParameter("Conservative") == 1)
	  {
	    m_variableType = eConservative;
	  }
	else
	  {
	    ASSERTL0(false,"Illegal value for variableType");
	  }
      }
    else
      {
	ASSERTL0(false,"variableType undefined");
      }
    //--------------------------------------------

    // HACK!!! the m_depth should be removed and replaced by d field[3]
    
    //--------------------------------------------
    // Set linear or nonlinear scheme
    if(m_boundaryConditions->CheckForParameter("NonLinear") == true)
      {
	if((int) m_boundaryConditions->GetParameter("NonLinear") == 0)
	  {
	    m_linearType = eLinear;

	    if(m_boundaryConditions->CheckForParameter("Depth") == true)
	      {
		m_depth = m_boundaryConditions->GetParameter("Depth");
	      }
	    else
	      {
		ASSERTL0(false,"Depth not specified");
	      }
	  }
	else if ((int) m_boundaryConditions->GetParameter("NonLinear") == 1)
	  {
	    m_linearType = eNonLinear;

	    if(m_boundaryConditions->CheckForParameter("Depth") == true)
	      {
		m_depth = m_boundaryConditions->GetParameter("Depth");
	      }
	  }
	else
	  {
	    ASSERTL0(false,"Illegal value for linearType");
	  }
      }
    else
      {
	ASSERTL0(false,"linearType undefined");
      }
    //--------------------------------------------

    
    //--------------------------------------------
    // Set Coriolis forcing if specified
    if(m_boundaryConditions->CheckForParameter("Coriolis") == true)
      {
	if (m_expdim == 2)
	  {
	    m_coriolis = Array<OneD, NekDouble> (GetTotPoints());
	    
	    EvaluateCoriolis();
	  }
	else
	  {
	    ASSERTL0(false,"Coriolis defined for 1D run");
	  }
      }
    //--------------------------------------------


    
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

    // HACK!!! Depth should be removed asap...
    // check if depth 
    if(m_boundaryConditions->CheckForParameter("Depth") == true)
      {
	m_depth = m_boundaryConditions->GetParameter("Depth");
      }
    else
      {
	ASSERTL0(false,"Depth not specified");
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
    int nq = m_fields[0]->GetTotPoints();
    
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
  
  
  //  void ShallowWaterEquations::ODEforcing(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
  // 					 Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time) 
//   {
//     int i;
//     int nVelDim    = m_spacedim;
//     int nvariables = inarray.num_elements();
//     int ncoeffs    = inarray[0].num_elements();
//     int nq         = GetTotPoints();
    
//     //-------------------------------------------------------
//     // go to physical space
    
//     Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
//     for (i = 0; i < nvariables; ++i)
//       {
// 	physarray[i] = Array<OneD, NekDouble>(nq);
// 	m_fields[i]->BwdTrans(inarray[i],physarray[i]);
//       }
//     //-------------------------------------------------------
    
//     SetBoundaryConditions(physarray, time);
    
//     switch(m_projectionType)
//       {
//       case eDiscontinuousGalerkin:
// 	{
	 
// 	  //-------------------------------------------------
// 	  // get the advection part
// 	  // input: physical space
// 	  // output: modal space 
// 	  WeakDGAdvection(physarray, outarray, false, true);

// 	  // negate the outarray since moving to the rhs
// 	  for(i = 0; i < nvariables; ++i)
// 	    {
// 	      Vmath::Neg(ncoeffs,outarray[i],1);
// 	    }
// 	  //-------------------------------------------------------
	  
	  
// 	  //-------------------------------------------------
// 	  // Add "source terms"
// 	  // input: physical space
// 	  // output: modal space
	  
// 	  //coriolis forcing
// 	  if (m_coriolis[0])
// 	    AddCoriolis(physarray,outarray);
// 	  //------------------------------------------------- 
	  
	  
// 	  //--------------------------------------------
// 	  // solve the block-diagonal system
	  
// 	  for(i = 0; i < nvariables; ++i)
// 	    {
// 	      m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
// 	    }
// 	  //--------------------------------------------

// 	}
// 	break;
//       case eGalerkin:
// 	{
// 	  ASSERTL0(false,"Continouos scheme not implemented for SWE");
// 	}
// 	break;
//       default:
// 	ASSERTL0(false,"Unknown projection scheme for the SWE");
// 	break;
//       }
//   }
  
  void ShallowWaterEquations::ODErhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
				           Array<OneD,       Array<OneD, NekDouble> >&outarray, 
				     const NekDouble time) 
  {
 int i;
    int nVelDim    = m_spacedim;
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
    
    
    //-------------------------------------------------------
    // set time dependent boundary conditions
    
    SetBoundaryConditions(physarray, time);
    //------------------------------------------------------

    switch(m_projectionType)
      {
      case eDiscontinuousGalerkin:
	{
	  
	  //-------------------------------------------------
	  // get the advection part
	  // input: physical space
	  // output: modal space 
	  
	  if (m_variableType == eConservative)
	    {
	      
	      // straighforward DG
	      // (note - that we only use the 3 dependent variables)
	      WeakDGAdvection(physarray, outarray, false, true,3);
	    }
	  else if (m_variableType == ePrimitive)
	    {
	      if (m_linearType == eLinear)
		{
		  
		  // for linear scheme we we see it as a conservation law
		  // providing the depth is constant
		  // (note - that we only use the 3 dependent variables)
		  WeakDGAdvection(physarray, outarray, false, true,3);
		}
	      else
		{
		  
		  // general nonconservative case
		  // can also be used for the linearized form (should be standard choice ??)
		  SWEAdvectionNonConservativeForm(physarray,outarray); 
		}
	    }
	  else
	    {
	      ASSERTL0(false,"Illegal variableType");
	    }
	  //-------------------------------------------------

	  
	  //-------------------------------------------------------
	  // negate the outarray since moving terms to the rhs
	  for(i = 0; i < nvariables; ++i)
	    {
	      Vmath::Neg(ncoeffs,outarray[i],1);
	    }
	  //-------------------------------------------------------
	  
	  
	  //-------------------------------------------------
	  // Add "source terms"
	  // input: physical space
	  // output: modal space
	  
	  // coriolis forcing
	  if (m_coriolis[0])
	    AddCoriolis(physarray,outarray);
	  //------------------------------------------------- 
	  
	}
	break;
      case eGalerkin:
	{
	  ASSERTL0(false,"Continouos scheme not implemented for SWE");
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the SWE");
	break;
      }
  }
 
void ShallowWaterEquations::ODElhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
    
  void ShallowWaterEquations::ODElhsSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
                                                       false,true);
              }
          }
          break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

    void ShallowWaterEquations::ODEdirkSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
						  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                                  const NekDouble lambda,
                                                  const NekDouble time) 
    {
        ASSERTL0(false, "this routine needs implementation");
    }
  
  
    void ShallowWaterEquations::ExplicitlyIntegrateAdvection(int nsteps)
    {
      int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Get Integration scheme details
        LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eClassicalRungeKutta4);
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
	    phys[i] = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints());
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
	// 	// this is just a forward Euler to illustate that CG works
		 
// 		// get -D u^n
// 		ODEforcing(in,out,m_time); // note that MultiplyByInvMassMatrix is not performed inside ODEforcing
	  
// 		// compute M u^n
// 		for (i = 0; i < nvariables; ++i)
// 		  {
// 		    m_fields[0]->BwdTrans(in[i],phys[i]);
// 		    m_fields[0]->IProductWRTBase(phys[i],tmp[i]);
		    
// 		    // f = M u^n - Dt D u^n
// 		    Vmath::Svtvp(ncoeffs,m_timestep,out[i],1,tmp[i],1,tmp[i],1);
		    
// 		    // u^{n+1} = M^{-1} f
// 		    m_fields[i]->MultiplyByInvMassMatrix(tmp[i],out[i],false,false);
		    
// 		    // fill results
// 		    Vmath::Vcopy(ncoeffs,out[i],1,in[i],1);
// 		    Vmath::Vcopy(ncoeffs,out[i],1,fields[i],1);
// 		  }
// 		//--------------------------------------------------------
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
	      
	      if (m_variableType == eConservative)
		{
		  ConservativeToPrimitive();
		  Checkpoint_Output(nchk++);
		  PrimitiveToConservative();
		}
	      else
		{
		  Checkpoint_Output(nchk++);
		}
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
    int cnt = 0;

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
		WallBoundary2D(n,cnt,inarray);
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
	cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }
  
  //----------------------------------------------------
 
  void ShallowWaterEquations::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
  { 

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; // only the dependent variables 
    
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
      }
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts;// cnt = 0; 
    
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
	// note that h(eta) in [0], hu(u) in [1] and hv(v) in [2]
	for (i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(npts,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	  }
      }
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
    // Vmath::Neg(nq,tmp,1); 
    Vmath::Vadd(ncoeffs,tmp,1,outarray[1],1,outarray[1],1);
    
    // add to hv equation
    Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
    Vmath::Neg(nq,tmp,1);
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
  
  
  void ShallowWaterEquations::GetFluxVector2D(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
					      Array<OneD, Array<OneD, NekDouble> > &flux)
  {
    
    // since this function can be called
    // from Boussinesq with i > 2
    if (i > 2)
      {
	return;
      }
  
    NekDouble g = m_g;

    switch(m_linearType)
      {
      case eLinear:
	
	switch(i){
	  
	  // flux function for the eta equation
	case 0:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j]  =  m_depth * physfield[1][j];
	      flux[1][j]  =  m_depth * physfield[2][j];
	    }
	  break;
	  
	  // flux function for the u equation
	case 1:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j] = g*physfield[0][j];
	      flux[1][j] = 0.0;
	    }
	  break;
	  
	  // flux function for the v equation
	case 2:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j] = 0.0;
	      flux[1][j] = g*physfield[0][j];
	    }
	  break;
	  
	default:
	  ASSERTL0(false,"GetFluxVector2D: illegal vector index");
	}
	break;
	
      case eNonLinear:
     	switch(i){
	  
	  // flux function for the eta equation
	case 0:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j]  =  physfield[1][j];
	      flux[1][j]  =  physfield[2][j];
	    }
	  break;
	  
	  // flux function for the hu equation
	case 1:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] +
		0.5*g*physfield[0][j]*physfield[0][j];
	      flux[1][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	    }
	  break;
	  
	  // flux function for the hv equation
	case 2:
	  for (int j = 0; j < m_fields[0]->GetTotPoints(); ++j)
	    {
	      flux[0][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	      flux[1][j] = physfield[2][j]*physfield[2][j]/physfield[0][j] +
		0.5*g*physfield[0][j]*physfield[0][j];
	    }
	  break;
	  
	default:
	  ASSERTL0(false,"GetFluxVector2D: illegal vector index");
	}
	break;
	
      default:
	ASSERTL0(false,"GetFluxVector2D: illegal lineartype");
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
     
      switch(m_linearType)
	{
	case eLinear:
	  LinearNumericalFlux2D(physfield, numfluxX, numfluxY);
	  break;
	case eNonLinear:
	  NonlinearNumericalFlux2D(physfield, numfluxX, numfluxY);
	  break;
	}
    }

  void ShallowWaterEquations::NonlinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
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
  
  void ShallowWaterEquations::LinearNumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						    Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
						    Array<OneD, Array<OneD, NekDouble> > &numfluxY)
  {
    int i;
    
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3;//m_fields.num_elements();
    
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
	
	
	NekDouble eta, u, v;
	NekDouble g = m_g;
	NekDouble d = m_depth;
	
	// averaging
	for (i = 0; i < nTraceNumPoints; ++i)
	  {
	    eta = 0.5*(Fwd[0][i] + Bwd[0][i]);
	    u   = 0.5*(Fwd[1][i] + Bwd[1][i]);
	    v   = 0.5*(Fwd[2][i] + Bwd[2][i]);
	    
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
  void ShallowWaterEquations::RiemannSolver(NekDouble hL,NekDouble huL,NekDouble hvL,NekDouble hR,NekDouble huR, 
					    NekDouble hvR, NekDouble &hflux, NekDouble &huflux,NekDouble &hvflux )
  {
    NekDouble g = m_g;
    
    NekDouble hC,huC,hvC,SL,SR,hstar,Sstar;
    
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
      
  void ShallowWaterEquations::SWEAdvectionNonConservativeForm(const Array<OneD, const Array<OneD, NekDouble> >&physarray,
							            Array<OneD,       Array<OneD, NekDouble> >&outarray)
  {
   
    //   ASSERTL0(false,"not implemented");
    //--------------------------------------
    // local parameters
    
    int ncoeffs         = outarray[0].num_elements();
    int nq              = GetTotPoints();
    int nTraceNumPoints = GetTraceTotPoints();
    
    NekDouble g        = m_g;
    NekDouble d        = m_depth; 
    //------------------------------------
  
    Array<OneD, NekDouble> tmp0(ncoeffs);
    Array<OneD, NekDouble> tmp1(ncoeffs);
    Array<OneD, NekDouble> tmp2(ncoeffs);
	
    Array<OneD, NekDouble> phys0(nq);
    Array<OneD, NekDouble> phys1(nq);
    Array<OneD, NekDouble> phys2(nq);

    //----------------------------------------
    // compute {\bf e}_1 = \nabla \eta
    
    Array<OneD, Array<OneD, NekDouble> > e1(2);
    e1[0] = Array<OneD, NekDouble> (nq); 
    e1[1] = Array<OneD, NekDouble> (nq);
    // boundary conditions already up to date
    //SetGradientBoundary(physarray[0],m_time,0);
    GradientFluxTerms(physarray[0],e1,0);
    //----------------------------------------
  

    //----------------------------------------
    // compute {\bf m}_1 = \nabla u
    
    Array<OneD, Array<OneD, NekDouble> > m1(2);
    m1[0] = Array<OneD, NekDouble> (nq); 
    m1[1] = Array<OneD, NekDouble> (nq);
    // boundary conditions already up to date
    //SetGradientBoundary(physarray[1],m_time,1);
    GradientFluxTerms(physarray[1],m1,1);
    //----------------------------------------
    

    //----------------------------------------
    // compute {\bf m}_2 = \nabla v
    
    Array<OneD, Array<OneD, NekDouble> > m2(2);
    m2[0] = Array<OneD, NekDouble> (nq); 
    m2[1] = Array<OneD, NekDouble> (nq);
    // boundary conditions already up to date
    //SetGradientBoundary(physarray[2],m_time,2);
    GradientFluxTerms(physarray[2],m2,2);
    //----------------------------------------
    
    
    //----------------------------------------
    // compute a_3 = \nabla \cdot (h{\bf u})
   

    // HACK! only correct for constant depth
    // (d*a3 below)
    // boundary conditions already up to date

    Array<OneD, NekDouble > a3(nq);
    {
      Array<OneD, Array<OneD, NekDouble> >in(2);

      // here we make sure we get deep copies
      in[0] = Array<OneD, NekDouble>(nq);
      Vmath::Vcopy(nq,physarray[1],1,in[0],1);
      in[1] = Array<OneD, NekDouble>(nq);
      Vmath::Vcopy(nq,physarray[2],1,in[1],1);
      
      
      //Vmath::Vmul(nq,physarray[3],1,in[0],1,in[0],1);
      //Vmath::Smul(nq,d,in[0],1,in[0],1);
      //Vmath::Vmul(nq,physarray[3],1,in[1],1,in[1],1);
      //Vmath::Smul(nq,d,in[1],1,in[1],1);
      //SetDivergenceBoundary(in,m_time,1,2);
      
      DivergenceFluxTerms(in,a3,1,2);
    } 
    //----------------------------------------
    
    
    //----------------------------------------
    // compute a_7 = \nabla \cdot (\eta {\bf u})
    
    Array<OneD, NekDouble > a7(nq);
    {
      Array<OneD, Array<OneD, NekDouble> >in(2);
      
      // here we make sure we get deep copies
      in[0] = Array<OneD, NekDouble>(nq);
      Vmath::Vcopy(nq,physarray[1],1,in[0],1);
      in[1] = Array<OneD, NekDouble>(nq);
      Vmath::Vcopy(nq,physarray[2],1,in[1],1);
      
      Vmath::Vmul(nq,physarray[0],1,in[0],1,in[0],1);
      Vmath::Vmul(nq,physarray[0],1,in[1],1,in[1],1);
      
      // needs to update BC
      SetDivergenceBoundary(in,m_time,1,2);
      
      DivergenceFluxTerms(in,a7,1,2);
    } 
    //----------------------------------------
    

    for (int i = 0; i < nq; ++i)
      {
	// linear part
	phys0[i] = d*a3[i];
  	phys1[i] = g*e1[0][i];
 	phys2[i] = g*e1[1][i];
	
	if (m_linearType == eNonLinear)
	  {
	    // nonlinear part
	    phys0[i] += a7[i];
	    phys1[i] += physarray[1][i]*m1[0][i] + physarray[2][i]*m1[1][i];
	    phys2[i] += physarray[1][i]*m2[0][i] + physarray[2][i]*m2[1][i];
	  }
      }
    
    // Get modal values
    m_fields[0]->IProductWRTBase(phys0,outarray[0]);
    m_fields[0]->IProductWRTBase(phys1,outarray[1]);
    m_fields[0]->IProductWRTBase(phys2,outarray[2]);

  }
  
  void ShallowWaterEquations::SetGradientBoundary(Array<OneD, NekDouble> &inarray,  
						  NekDouble time, int field_0)
 {
   
   int cnt = 0;
   
   // loop over Boundary Regions
   for (int n = 0; n < m_fields[field_0]->GetBndConditions().num_elements(); ++n)
     {	
       
       // Wall Boundary Condition
       if (m_fields[field_0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
	 {
	   GradientWallBoundary2D(n,cnt,inarray,field_0);
	 }
       
       cnt +=m_fields[field_0]->GetBndCondExpansions()[n]->GetExpSize();
     }
 }
  
  //----------------------------------------------------
  
  void ShallowWaterEquations::GradientWallBoundary2D(int bcRegion, int cnt, Array<OneD, NekDouble> &physarray, int field_0)
  { 

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
        
    // get physical values of h, hu, hv for the forward trace
    Array<OneD, NekDouble> Fwd(nTraceNumPoints);
    m_fields[field_0]->ExtractTracePhys(physarray,Fwd);
    
    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts;
    
    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
	
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&Fwd[id2], 1,&(m_fields[field_0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
      }
  }
  //----------------------------------------------------

  
  void ShallowWaterEquations::SetDivergenceBoundary(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time, int field_0, int field_1)
  {
    
    int cnt = 0;

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[field_0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Wall Boundary Condition
	if (m_fields[field_0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
	  {
	    DivergenceWallBoundary2D(n,cnt,inarray,field_0,field_1);
	  }
	
	cnt +=m_fields[field_0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }
  
  //----------------------------------------------------
 
  void ShallowWaterEquations::DivergenceWallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray, int field_0, int field_1)
  { 

    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 2;
    
    // get physical values for the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
      }
    
    m_fields[field_0]->ExtractTracePhys(physarray[0],Fwd[0]);
    m_fields[field_1]->ExtractTracePhys(physarray[1],Fwd[1]);


    // Adjust the physical values of the trace to take 
    // user defined boundaries into account
    int e, id1, id2, npts;
    
    for(e = 0; e < m_fields[field_0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
	
	
	Array<OneD, NekDouble> tmp_n(npts);
	Array<OneD, NekDouble> tmp_t(npts);
	      
	Vmath::Vmul(npts,&Fwd[0][id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
	Vmath::Vvtvp(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);
	      
	Vmath::Vmul(npts,&Fwd[0][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
	Vmath::Vvtvm(npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
	
	// negate the normal flux
	Vmath::Neg(npts,tmp_n,1);		      
	      
	// rotate back to Cartesian
	Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[0][id2],1);
	Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[0][id2],1,&Fwd[0][id2],1);
	
	Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1);
	Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);
      
	// copy boundary adjusted values into the boundary expansion
	Vmath::Vcopy(npts,&Fwd[0][id2], 1,&(m_fields[field_0]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	Vmath::Vcopy(npts,&Fwd[1][id2], 1,&(m_fields[field_1]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);

	
      }
  }

  // in and out in physical space
  void ShallowWaterEquations::GradientFluxTerms(const Array<OneD, NekDouble> &in, 
						      Array<OneD, Array<OneD, NekDouble> > &out,
						int field_0)
  {
    int ncoeffs         = GetNcoeffs();
    int nTraceNumPoints = GetTraceTotPoints();
    
    Array<OneD, NekDouble> tmpX(ncoeffs);
    Array<OneD, NekDouble> tmpY(ncoeffs);
    
    Array<OneD, NekDouble> upwindX(nTraceNumPoints);
    Array<OneD, NekDouble> upwindY(nTraceNumPoints);
    Array<OneD, NekDouble> upwindZero(nTraceNumPoints,0.0);
    
    m_fields[0]->IProductWRTDerivBase(0,in,tmpX);
    m_fields[0]->IProductWRTDerivBase(1,in,tmpY);
    
    Vmath::Neg(ncoeffs,tmpX,1);
    Vmath::Neg(ncoeffs,tmpY,1);

    NumericalFluxGradient(in,upwindX,upwindY,field_0);
    
    m_fields[0]->AddTraceIntegral(upwindX,upwindZero,tmpX);
    m_fields[0]->AddTraceIntegral(upwindZero,upwindY,tmpY);
    
    m_fields[0]->MultiplyByElmtInvMass(tmpX,tmpX);
    m_fields[0]->MultiplyByElmtInvMass(tmpY,tmpY);

    m_fields[0]->BwdTrans(tmpX,out[0]);
    m_fields[0]->BwdTrans(tmpY,out[1]);
  }
  

  // in and out in physical space
  void ShallowWaterEquations::DivergenceFluxTerms(const Array<OneD, const Array<OneD, NekDouble> >&in, 
						  Array<OneD, NekDouble> &out,
						  int field_0, int field_1)
  {
    int ncoeffs         = GetNcoeffs();
    int nTraceNumPoints = GetTraceTotPoints();
    
    Array<OneD, NekDouble> tmpX(ncoeffs);
    Array<OneD, NekDouble> tmpY(ncoeffs);
    
    Array<OneD, NekDouble> upwindX(nTraceNumPoints);
    Array<OneD, NekDouble> upwindY(nTraceNumPoints);
 
    m_fields[0]->IProductWRTDerivBase(0,in[0],tmpX);
    m_fields[0]->IProductWRTDerivBase(1,in[1],tmpY);
    
    Vmath::Vadd(ncoeffs,tmpX,1,tmpY,1,tmpX,1);
    Vmath::Neg(ncoeffs,tmpX,1);
    
    NumericalFluxDivergence(in,upwindX,upwindY,field_0,field_1);
    
    m_fields[0]->AddTraceIntegral(upwindX,upwindY,tmpX);
    
    m_fields[0]->MultiplyByElmtInvMass(tmpX,tmpX);

    m_fields[0]->BwdTrans(tmpX,out);
    
  }
  
  
  /**
   * Computes the \hat{f} term in the  \int_{\partial \Omega^e} \phi \hat{f} {\bf n} dS 
   * integral. Using averaged fluxes.
   **/
  void ShallowWaterEquations::NumericalFluxGradient(const Array <OneD, NekDouble> &physarray,
						    Array<OneD, NekDouble> &outX,
						    Array<OneD, NekDouble> &outY,
						    int field_0)
  {
    int nTraceNumPoints = GetTraceTotPoints();

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(1);
    Array<OneD, Array<OneD, NekDouble> > Bwd(1);
    
    Fwd[0] = Array<OneD, NekDouble> (nTraceNumPoints,0.0); 
    Bwd[0] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
    
    // get the physical values at the trace
    m_fields[field_0]->GetFwdBwdTracePhys(physarray,Fwd[0],Bwd[0]);
    
    for (int i = 0; i < nTraceNumPoints; ++i)
      {
 	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
      }

  }

  
  /**
   * Computes the \hat{\bf f} term in the  \int_{\partial \Omega^e} \phi \hat{\bf f} \cdot {\bf n} dS 
   * integral. Using averaged fluxes.
   **/
  void ShallowWaterEquations::NumericalFluxDivergence(const Array <OneD, const Array<OneD, NekDouble> > &physarray,
						      Array <OneD, NekDouble> &outX, 
						      Array<OneD, NekDouble> &outY,
						      int field_0, int field_1)
  {
    int nTraceNumPoints = GetTraceTotPoints();
    
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(2);
    Array<OneD, Array<OneD, NekDouble> > Bwd(2);
    
    for (int i = 0; i < 2; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
      }
    
    // get the physical values at the trace
    m_fields[field_0]->GetFwdBwdTracePhys(physarray[0],Fwd[0],Bwd[0]);
    m_fields[field_1]->GetFwdBwdTracePhys(physarray[1],Fwd[1],Bwd[1]);
    
    for (int i = 0; i < nTraceNumPoints; ++i)
      {
	outX[i]  = 0.5*(Fwd[0][i] + Bwd[0][i]);
	outY[i]  = 0.5*(Fwd[1][i] + Bwd[1][i]);
      }
  }
  

  void ShallowWaterEquations::PrimitiveToConservative(void)
  {
    
    int nq = GetTotPoints();
    
    // physical space
    for (int i = 0; i < 4; ++i)
      {
	m_fields[i]->BwdTrans(*m_fields[i]);
      }
    
    // h = eta + d
    Vmath::Vadd(nq,m_fields[0]->GetPhys(),1,m_fields[3]->GetPhys(),1,m_fields[0]->UpdatePhys(),1);

    // hu = h * u
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[1]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // hv = h * v
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[2]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);

    // modal space
    for (int i = 0; i < 4; ++i)
      {
	m_fields[i]->FwdTrans(*m_fields[i]);
      }
  }
  
  void ShallowWaterEquations::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
						            Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    
    int nq = GetTotPoints();
    
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(4);
	for (int i = 0; i < 4; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	
	// h = \eta + d
	Vmath::Vadd(nq,tmp[0],1,tmp[3],1,physout[0],1);
	
	// hu = h * u
	Vmath::Vmul(nq,tmp[0],1,tmp[1],1,physout[1],1);
	
	// hv = h * v
	Vmath::Vmul(nq,tmp[0],1,tmp[2],1,physout[2],1);
      
      }
    else
      {
	// h = \eta + d
	Vmath::Vadd(nq,physin[0],1,physin[3],1,physout[0],1);
	
	// hu = h * u
	Vmath::Vmul(nq,physin[0],1,physin[1],1,physout[1],1);
	
	// hv = h * v
	Vmath::Vmul(nq,physin[0],1,physin[2],1,physout[2],1);
	
      }
  }
  
  void ShallowWaterEquations::ConservativeToPrimitive(void)
  {
    
    int nq = GetTotPoints();
    
    // physical space
    for (int i = 0; i < 4; ++i)
      {
	m_fields[i]->BwdTrans(*m_fields[i]);
      }
    
    // u = hu/h
    Vmath::Vdiv(nq,m_fields[1]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // v = hv/h
    Vmath::Vdiv(nq,m_fields[2]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);
    
    // eta = h - d
    Vmath::Vsub(nq,m_fields[0]->GetPhys(),1,m_fields[3]->GetPhys(),1,m_fields[0]->UpdatePhys(),1);
    
    // modal space
    for (int i = 0; i < 4; ++i)
      {
	m_fields[i]->FwdTrans(*m_fields[i]);
      }
  }

  void ShallowWaterEquations::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
						            Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    int nq = GetTotPoints();
      
    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(4);
	for (int i = 0; i < 4; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }
	
	// \eta = h - d
	Vmath::Vsub(nq,tmp[0],1,tmp[3],1,physout[0],1);
	
	// u = hu/h
	Vmath::Vdiv(nq,tmp[1],1,tmp[0],1,physout[1],1);
	
	// v = hv/ v
	Vmath::Vdiv(nq,tmp[2],1,tmp[0],1,physout[2],1);
      }
    else
      {
	// \eta = h - d
	Vmath::Vsub(nq,physin[0],1,physin[3],1,physout[0],1);
	
	// u = hu/h
	Vmath::Vdiv(nq,physin[1],1,physin[0],1,physout[1],1);
	
	// v = hv/ v
	Vmath::Vdiv(nq,physin[2],1,physin[0],1,physout[2],1);
      }
  }
    
  void ShallowWaterEquations::Summary(std::ostream &out)
  {
    cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : Shallow Water Equations" << endl;
      cout << "\t                  eta should be in field[0]" <<endl;
      cout << "\t                  u   should be in field[1]" <<endl;
      cout << "\t                  v   should be in field[2]" <<endl;
      cout << "\t                  d   should be in field[3]" <<endl;
      ADRBase::Summary(out);
      cout << "=======================================================================" << endl;
      cout << endl;
    
    }


} //end of namespace

/**
* $Log: ShallowWaterEquations.cpp,v $
* Revision 1.3  2009/02/06 16:38:23  claes
* Added primitive formulation
*
* Revision 1.2  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.1  2008/11/17 08:37:00  claes
* Restructured Shallow Water Solver. Only 2D DG is properly working
*
* Revision 1.2  2008/09/15 14:54:15  claes
* Small changes associated with the BoussinesqSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
