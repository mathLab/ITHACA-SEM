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

      
    //--------------------------------------------
    // Set linear or nonlinear scheme
    if(m_boundaryConditions->CheckForParameter("NonLinear") == true)
      {
	if((int) m_boundaryConditions->GetParameter("NonLinear") == 0)
	  {
	    m_linearType = eLinear;
	  }
	else if ((int) m_boundaryConditions->GetParameter("NonLinear") == 1)
	  {
	    m_linearType = eNonLinear;
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
    // Set still water depth
    EvaluateDepth();

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
    else
      {
	m_coriolis = Array<OneD, NekDouble>();
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

  void ShallowWaterEquations::EvaluateDepth(void)
  {
    int nq = m_fields[0]->GetTotPoints();
    
    m_depth      = Array<OneD, NekDouble >(nq);
    m_d_depth    = Array<OneD, Array<OneD, NekDouble> >(2);
    m_d_depth[0] = Array<OneD, NekDouble>(nq);
    m_d_depth[1] = Array<OneD, NekDouble>(nq);

    std::string coriolisStr[1] = {"d"};
    
    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);
    
    // get the coordinates (assuming all fields have the same
    // discretisation)
    m_fields[0]->GetCoords(x0,x1,x2);
    
    SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(coriolisStr[0]);
    
    for(int j = 0; j < nq; j++)
      {
	m_depth[j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
      }

    
    // compute d_x and d_y
    // this is needed for Boussinesq equations to run smoothly...
    // SetGradientBoundary(m_depth,0.0,1);
    // GradientFluxTerms(m_depth,m_d_depth,1);
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
  
  
  void ShallowWaterEquations::ODErhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
				           Array<OneD,       Array<OneD, NekDouble> >&outarray, 
				     const NekDouble time) 
  {
    int i;
    int ndim    = m_spacedim;
    int nvariables = inarray.num_elements();
    int ncoeffs    = inarray[0].num_elements();
    int nq         = GetTotPoints();
    
    switch(m_projectionType)
      {
      case eDiscontinuousGalerkin:
	{
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
	  
	  //-------------------------------------------------
	  // get the advection part
	  // input: physical space
	  // output: modal space 
	  
	  if (m_variableType == eConservative)
	    {
	      
	      // straighforward DG
	      WeakDGAdvection(physarray, outarray, false, true);
	    }
	  else if (m_variableType == ePrimitive)
	    {
	      if (m_linearType == eLinear)
		{
		  
		  // for linear scheme we we see it as a conservation law
		  // providing the depth is constant
		  WeakDGAdvection(physarray, outarray, false, true);
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
	  if (m_coriolis.num_elements() != 0)
	    {
	      AddCoriolis(physarray,outarray);
	    }
	  //------------------------------------------------- 
	   
	  for(i = 0; i < nvariables; ++i)
	    {
	      m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
	    }
	}
	break;
      case eGalerkin:
	{
	  if (m_variableType == ePrimitive)
	    {
	      ASSERTL0(false,"primitive scheme not implemented for CG SWE");
	    }
	  
	  //SetCGBoundaryConditions(time);
	  
	  // Go to physical space 
	  Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
	  Array<OneD, Array<OneD, NekDouble> > fluxvector(ndim);
	  for(i = 0; i < ndim; ++i)
	    {
	      fluxvector[i]    = Array<OneD, NekDouble>(nq);
	    }

	  for (i = 0; i < nvariables; ++i)
	    {
	      physarray[i] = Array<OneD, NekDouble>(nq);
	      m_fields[i]->MultiplyByInvMassMatrix(inarray[i],  
						   outarray[i],
						   false);
	      m_fields[i]->BwdTrans_IterPerExp(outarray[i],physarray[i]);
	    }

	  
	  Array<OneD,NekDouble> tmp(nq);
	  Array<OneD, NekDouble>tmp1(nq);           
	  // tmp = Array<OneD, NekDouble> (ndim*nq);
// 	  tmp1 = tmp + nq;

	  for(i = 0; i < nvariables; ++i)
	    {
	      
	      // Get the ith component of the  flux vector in (physical space)
	      ShallowWaterEquations::GetFluxVector2D(i, physarray, fluxvector);
         
	      m_fields[0]->PhysDeriv(0,fluxvector[0],tmp);
	      m_fields[0]->PhysDeriv(1,fluxvector[1],tmp1);
	      Vmath::Vadd(nq,tmp,1,tmp1,1,tmp,1);
	      m_fields[0]->IProductWRTBase_IterPerExp(tmp,outarray[i]);
	      Vmath::Neg(ncoeffs,outarray[i],1);
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

	  //ASSERTL0(false,"Continouos scheme not implemented for SWE");
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
   MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
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
                                                       false);
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
  
  void ShallowWaterEquations::GeneralTimeIntegration(int nsteps, 
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
            MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
            for(int i = 0; i < nvariables; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(ncoeffs);
                m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
            }

	    for(int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->SetPhysState(false);
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
	      // ASSERTL0(false,"CG not implemented for SWE");
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
	      //SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->MultiplyByInvMassMatrix(fields[i],tmp[i],false);
                    fields[i] = tmp[i];	   		    
                }
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
	      
	      if (m_variableType == eConservative)
		{
		  ConservativeToPrimitive();
		  Checkpoint_Output(nchk++);
		  PrimitiveToConservative();
		  //dump depth
		  //Array_Output(nchk,"depth",m_depth,true);
		}
	      else
		{
		  Checkpoint_Output(nchk++);
		  //dump depth
		  //Array_Output(nchk,"depth",m_depth,true);
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
  }
  
    void ShallowWaterEquations::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {     
    ASSERTL0(false,"1D not working with user defined BC");
  }


  void ShallowWaterEquations::AddCoriolis(Array<OneD, Array<OneD, NekDouble> > &physarray,
					  Array<OneD, Array<OneD, NekDouble> > &outarray)
  {
   
    // routine works for both primitive and conservative formulations
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
    int nq = m_fields[0]->GetTotPoints();
    
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
	  for (int j = 0; j < nq; ++j)
	    {
	      flux[0][j]  =  m_depth[j] * physfield[1][j];
	      flux[1][j]  =  m_depth[j] * physfield[2][j];
	    }
	  break;
	  
	  // flux function for the u equation
	case 1:
	  for (int j = 0; j < nq; ++j)
	    {
	      flux[0][j] = g*physfield[0][j];
	      flux[1][j] = 0.0;
	    }
	  break;
	  
	  // flux function for the v equation
	case 2:
	  for (int j = 0; j < nq; ++j)
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
	  for (int j = 0; j < nq; ++j)
	    {
	      flux[0][j]  =  physfield[1][j];
	      flux[1][j]  =  physfield[2][j];
	    }
	  break;
	  
	  // flux function for the hu equation
	case 1:
	  for (int j = 0; j < nq; ++j)
	    {
	      flux[0][j] = physfield[1][j]*physfield[1][j]/physfield[0][j] +
		0.5*g*physfield[0][j]*physfield[0][j];
	      flux[1][j] = physfield[1][j]*physfield[2][j]/physfield[0][j];
	    }
	  break;
	  
	  // flux function for the hv equation
	case 2:
	  for (int j = 0; j < nq; ++j)
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
    int nvariables      = 4;// we need the boundary values of the depth//m_fields.num_elements();
    
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
    
    for (i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
      }
    
    // get the physical values at the trace from the dependent variables
    for (i = 0; i < nvariables-1; ++i)
      {
	m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
      }
    
    // we also needs the values of d...
    // need to fill the in m_fields[0]
    SetGradientBoundary(m_depth,m_time,0);
    m_fields[0]->GetFwdBwdTracePhys(m_depth,Fwd[3],Bwd[3]);

    

    NekDouble eta, u, v, d;
    NekDouble g = m_g;
	
	
	// averaging
	for (i = 0; i < nTraceNumPoints; ++i)
	  {
	    eta = 0.5*(Fwd[0][i] + Bwd[0][i]);
	    u   = 0.5*(Fwd[1][i] + Bwd[1][i]);
	    v   = 0.5*(Fwd[2][i] + Bwd[2][i]);
	    d   = 0.5*(Fwd[3][i] + Bwd[3][i]);

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
    
    //--------------------------------------
    // local parameters
    
    int ncoeffs         = outarray[0].num_elements();
    int nq              = GetTotPoints();
    int nTraceNumPoints = GetTraceTotPoints();
    
    NekDouble g        = m_g;
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
    // compute a_3 = \nabla \cdot (d{\bf u})
   

    Array<OneD, NekDouble > a3(nq);
    {
      Array<OneD, Array<OneD, NekDouble> >in(2);

      // here we make sure we get deep copies
      in[0] = Array<OneD, NekDouble>(nq);
      Vmath::Vmul(nq,physarray[1],1,m_depth,1,in[0],1);
      in[1] = Array<OneD, NekDouble>(nq);
      Vmath::Vmul(nq,physarray[2],1,m_depth,1,in[1],1);
      
      SetDivergenceBoundary(in,m_time,1,2);
      
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
	phys0[i] = a3[i];
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
 // in and out in physical space
  void ShallowWaterEquations::GradientFluxUpwindTerms(Array<OneD, NekDouble> &in, 
						      Array<OneD, Array<OneD, NekDouble> > &out,
						      int field_0,
						      const Array<OneD, Array<OneD, NekDouble> > &phys)
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

    NumericalFluxUpwindGradient(in,upwindX,upwindY,field_0,phys);
    
    m_fields[0]->AddTraceIntegral(upwindX,upwindZero,tmpX);
    m_fields[0]->AddTraceIntegral(upwindZero,upwindY,tmpY);
    
    m_fields[0]->MultiplyByElmtInvMass(tmpX,tmpX);
    m_fields[0]->MultiplyByElmtInvMass(tmpY,tmpY);

    m_fields[0]->BwdTrans(tmpX,out[0]);
    m_fields[0]->BwdTrans(tmpY,out[1]);
  }
  void ShallowWaterEquations::NumericalFluxUpwindGradient(Array<OneD, NekDouble> &in, 
							  Array<OneD, NekDouble> &outX,
							  Array<OneD, NekDouble> &outY,
							  int field_0,
							  const Array<OneD, Array<OneD, NekDouble> > &phys)
  {
    
    int nTraceNumPoints = GetTraceTotPoints();

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(1);
    Array<OneD, Array<OneD, NekDouble> > Bwd(1);
    
    Fwd[0] = Array<OneD, NekDouble> (nTraceNumPoints,0.0); 
    Bwd[0] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
    Array<OneD, NekDouble>uFwd (nTraceNumPoints,0.0);
    Array<OneD, NekDouble>uBwd (nTraceNumPoints,0.0);
    
    Array<OneD, NekDouble>vFwd (nTraceNumPoints,0.0);
    Array<OneD, NekDouble>vBwd (nTraceNumPoints,0.0);


    // get the physical values at the trace
    m_fields[field_0]->GetFwdBwdTracePhys(in,Fwd[0],Bwd[0]);
    m_fields[1]->GetFwdBwdTracePhys(phys[1],uFwd,uBwd);
    m_fields[2]->GetFwdBwdTracePhys(phys[2],vFwd,vBwd);

    
    
    for (int i = 0; i < nTraceNumPoints; ++i)
      {
	if (uFwd[i]*m_traceNormals[0][i]+
	    vFwd[i]*m_traceNormals[1][i] >= 0.0)
	  {
	    outX[i]  = Fwd[0][i];
	    outY[i]  = Fwd[0][i];
	  }
	else
	  {
	    outX[i]  = Bwd[0][i];
	    outY[i]  = Bwd[0][i];
	  }
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
    for (int i = 0; i < 3; ++i)
      {
	m_fields[i]->BwdTrans((m_fields[i]->GetCoeffs()),(m_fields[i]->UpdatePhys()));
      }
    
    // h = eta + d
    Vmath::Vadd(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);

    // hu = h * u
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[1]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // hv = h * v
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[2]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);

    // modal space
    for (int i = 0; i < 3; ++i)
      {
	m_fields[i]->FwdTrans((m_fields[i]->GetPhys()),(m_fields[i]->UpdateCoeffs()));
      }
  }
  
  void ShallowWaterEquations::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
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
	Vmath::Vmul(nq,tmp[0],1,tmp[1],1,physout[1],1);
	
	// hv = h * v
	Vmath::Vmul(nq,tmp[0],1,tmp[2],1,physout[2],1);
      
      }
    else
      {
	// h = \eta + d
	Vmath::Vadd(nq,physin[0],1,m_depth,1,physout[0],1);
	
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
    for (int i = 0; i < 3; ++i)
      {
	m_fields[i]->BwdTrans((m_fields[i]->GetCoeffs()),(m_fields[i]->UpdatePhys()));
      }
    
    // u = hu/h
    Vmath::Vdiv(nq,m_fields[1]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // v = hv/h
    Vmath::Vdiv(nq,m_fields[2]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);
    
    // eta = h - d
    Vmath::Vsub(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);
    
    // modal space
    for (int i = 0; i < 3; ++i)
      {
	m_fields[i]->FwdTrans((m_fields[i]->GetPhys()),(m_fields[i]->UpdateCoeffs()));
      }
  }

  void ShallowWaterEquations::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
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
    
  void ShallowWaterEquations::Summary(std::ostream &out)
  {
    cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : Shallow Water Equations" << endl;
      if (m_variableType == ePrimitive)
	{
	  cout << "\t                  Primitive varibles" << endl;
	}
      else if (m_variableType == eConservative)
	{
	  cout << "\t                  Conservative variables" << endl;
	}
      if (m_linearType == eNonLinear)
	{
	  cout << "\t                  Nonlinear equations " << endl;
	}
      else if (m_linearType == eLinear)
	{
	  cout << "\t                  Linearized equations" << endl;
	}
      cout << "\t                  eta should be in field[0]" <<endl;
      cout << "\t                  u   should be in field[1]" <<endl;
      cout << "\t                  v   should be in field[2]" <<endl;
      ADRBase::Summary(out);
      cout << "=======================================================================" << endl;
      cout << endl;
    
    }


} //end of namespace

/**
* $Log: ShallowWaterEquations.cpp,v $
* Revision 1.6  2009/04/28 10:17:41  pvos
* Some updates to make the solvers compile properly with the newly added sparse matrix library
*
* Revision 1.5  2009/03/10 23:37:14  claes
* Updated the ShallowWaterSolver to work with the general timestepping scheme
*
* Revision 1.4  2009/02/07 23:58:08  claes
* Changed so I/O always are in terms of primitive variables
*
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
