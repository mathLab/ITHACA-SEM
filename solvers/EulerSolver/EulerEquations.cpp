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

    if(m_boundaryConditions->SolverInfoExists("UPWINDTYPE"))
      {

	std::string UpwindTypeStr = m_boundaryConditions->GetSolverInfo("UPWINDTYPE");
	int i;
	for(i = 0; i < (int) SIZE_UpwindType; ++i)
	  {
	    if(NoCaseStringCompare(UpwindTypeMap[i],UpwindTypeStr) == 0)
	      {
		m_upwindType = (UpwindType)i;
		break;
	      }
	  }
      }
    else
      {
	m_upwindType = (UpwindType)0;  // Upwind flux scheme not set
      }

    if(m_boundaryConditions->SolverInfoExists("PROBLEMTYPE"))
      {

	std::string ProblemTypeStr = m_boundaryConditions->GetSolverInfo("PROBLEMTYPE");
	int i;
	for(i = 0; i < (int) SIZE_ProblemType; ++i)
	  {
	    if(NoCaseStringCompare(ProblemTypeMap[i],ProblemTypeStr) == 0)
	      {
		m_problemType = (ProblemType)i;
		break;
	      }
	  }
      }
    else
      {
	m_problemType = (ProblemType)0;
      }
    
    // check that any user defined boundary condition is implemented
    
  }


  void EulerEquations::ODElhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
	MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
	for(int i = 0; i < nvariables; ++i)
	  {
	    tmp[i] = Array<OneD, NekDouble>(ncoeffs);
	    m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,fields[i],fields[i]);
	  }
      }

    // Derived Variables:
    Derived_field[0] = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*m_fields[0]);

    // Declaring the DG-CFL limit
    int num_mods = GetNumExpModes(); //ADRBase.h
    NekDouble CFLDG = GetCFLNumber(num_mods);
    NekDouble CFL;
    Array<OneD, NekDouble> MinLength(m_fields[0]->GetExpSize());
    GetMinLength(MinLength);

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
	  
	  // Maybe it is better to save it somewhere in the integration class
	  NekDouble CFLRK = 2.784;
	  CFL = CFLRK/CFLDG;
	  
	  IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
	  
	  LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
	  IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
	  
	  u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
	}
	break;
      case LibUtilities::eAdamsBashforthOrder2:
	{

	  // Maybe it is better to save it somewhere in the integration class
	  NekDouble CFLRK = 1.0;
	  CFL = CFLRK/CFLDG;

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
 
	
    n = 0;
    double check = m_timecheck;

    while(m_time<m_fintime)//while(n<nsteps)
      {
	//----------------------------------------------
	// Perform time step integration
	//---------------------------------------------- 
	GetTimeStep(CFL,MinLength,fields,m_timestep);

	if( n < numMultiSteps-1)
	  {
	    // Use initialisation schemes
	    fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
	  }
	else
	  {
	    fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
	  }
	
	if(m_time+m_timestep>m_fintime)
	  m_timestep = m_fintime - m_time;
	
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
	    cout << "Steps: " << n+1 << "\t Time: " << m_time <<  "\t TimeStep: " << m_timestep << endl;
	  }
	
	//if(n&&(!((n+1)%m_checksteps)))
	if(m_time>=m_timecheck || m_time>=m_fintime)
	  { 
	    nchk++;
	    cout << "Printing file: "<< nchk << ".chk" << endl;
	    for(i = 0; i < nvariables; ++i)
	      {
		(m_fields[i]->UpdateCoeffs()) = fields[i];
	      }
	    
	    // Extracting primitive variables on the wall and Mach in the field
	    int nq         = m_fields[0]->GetTotPoints(); 
	    Array<OneD, NekDouble> mach(nq);
	    Array<OneD, NekDouble> pressure(nq);
	    Array<OneD, NekDouble> soundspeed(nq);
	    Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);

	    for (i = 0; i < nvariables; ++i)
	      {
		physarray[i] = Array<OneD, NekDouble>(nq);
		m_fields[i]->BwdTrans(fields[i],physarray[i]);
	      }
	    GetPressure(physarray,pressure);
	    GetSoundSpeed(physarray,pressure,soundspeed);
	    GetMach(physarray,soundspeed,mach);
	    
	    // Writing the binary Mach file 
	    WriteVar(nchk,Derived_field,mach,"Mach");
	    
	    // Writing Wall file in gnuplot
	    ExtractWall(nchk,physarray);
	    
	    // Writing binary file of the solution
	    Checkpoint_Output(nchk);
	    
	    m_timecheck += check;
	  }    
	n++;
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
    
    
    //-------------------------------------------------------
    // set time dependent boundary conditions
    
    SetBoundaryConditions(physarray, time);
    //------------------------------------------------------
    
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


  //----------------------------------------------------
  void EulerEquations::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    int cnt = 0;

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Wall Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Wall")
	  {
	    if (m_expdim == 2)
	      {
		WallBoundary(n,cnt,physarray);
	      }
	    else
	      {
		ASSERTL0(false,"1D, 3D not yet implemented");
	      }
	  }

	// Symmetric Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "Symmetry")
	  {
	    if (m_expdim == 2)
	      {
		SymmetryBoundary(n,cnt,physarray);
	      }
	    else
	      {
		ASSERTL0(false,"1D, 3D not yet implemented");
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

	// Ringleb Flow Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "RinglebFlow")
	  {
	    SetBoundaryRinglebFlow(n,cnt,physarray);
	  }
	
	cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }
 
  
  //----------------------------------------------------
 
  void EulerEquations::WallBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
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
	      // negate the forward flux
	      Vmath::Neg(npts,&Fwd[1][id2],1);	
	      
	      ASSERTL0(false,"1D not yet implemented for the Euler");
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
	    ASSERTL0(false,"3D not implemented for the Euler Equations");
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

  void EulerEquations::ExtractWall(int nchk, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {
    int cnt = 0;
    for(int bcRegion = 0; bcRegion < m_fields[0]->GetBndConditions().num_elements(); ++bcRegion)
      {	
	// Wall Boundary Condition
	if (m_fields[0]->GetBndConditions()[bcRegion]->GetUserDefined().GetEquation() == "Wall")
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
	    
	    int e, id2, npts;
	    
	    NekDouble L = 0.0;
	    char buffer[1024];
	    sprintf(buffer,"WallRegion%d_%d.data",bcRegion,nchk);
	    FILE* out=fopen(buffer,"w");
	    
	    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
	      {
		npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
		id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e)); 
		
		// Get quadrature points
		Array<OneD,NekDouble> x0(npts,0.0);
		Array<OneD,NekDouble> x1(npts,0.0);
		Array<OneD,NekDouble> x2(npts,0.0); 
		m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);
		
		// Tangential velocity
		Array<OneD, NekDouble> tmp_t(npts);	
		Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
		Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
		
		
		fprintf(out, "%f\t%f\t%f\t%f\t%f\t%f\n",L,Fwd[0][id2],Fwd[1][id2],Fwd[2][id2],Fwd[3][id2],fabs(tmp_t[0])/Fwd[0][id2]);
		// write the wall file
		for (i = 1; i < npts; ++i)
		  {
		    NekDouble dx = x0[i]-x0[i-1];
		    NekDouble dy = x1[i]-x1[i-1];
		    NekDouble dl  = sqrt(1.0+dy/dx*dy/dx)*fabs(dx);
		    L = L + dl;
		    fprintf(out, "%f\t%f\t%f\t%f\t%f\t%f\n",L,Fwd[0][id2+i],Fwd[1][id2+i],Fwd[2][id2+i],Fwd[3][id2+i],tmp_t[i]/Fwd[0][id2+1]);
		  }
	      }
	    fclose(out);
	  } 
	cnt +=m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
      } 
  }

  //----------------------------------------------------
  
  void EulerEquations::SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
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
	      ASSERTL0(false,"1D not yet implemented for the Euler");
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
	    ASSERTL0(false,"3D not implemented for the Euler Equations");
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

  void EulerEquations::GetSoundSpeed(Array<OneD, Array<OneD, NekDouble> > &physfield,
				     Array<OneD, NekDouble> &pressure,
				     Array<OneD, NekDouble> &soundspeed)
  {
    NekDouble gamma = m_gamma;
    
    for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
      {
	soundspeed[i] = sqrt(gamma*pressure[i]/physfield[0][i]);
      }
  }

  void EulerEquations::GetMach(Array<OneD, Array<OneD, NekDouble> > &physfield,
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
      case eUpwind:
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
      case eRoe:
      case eExact:
	{
	  ExactRiemannSolver(rhoL,rhouL,rhovL,EL,rhoR,rhouR,rhovR,ER,rhoflux,rhouflux,rhovflux,Eflux);
	}
	break;
      default:
	{
	  ASSERTL0(false,"populate switch statement for upwind flux");
	}
	break;
      }
  }


  void EulerEquations::HLLRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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

  void EulerEquations::HLLCRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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

  void EulerEquations::LFRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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

  void EulerEquations::AverageRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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

  void EulerEquations::ExactRiemannSolver(NekDouble rhoL, NekDouble rhouL, NekDouble rhovL, NekDouble EL,
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
  
  
  void EulerEquations::GetMinLength(Array<OneD, NekDouble> &MinLength)
  {
    int n_element  = m_fields[0]->GetExpSize();                // number of element in the mesh
    for(int el = 0; el < n_element; ++el)
      {
	int n_edges = m_fields[0]->GetExp(el)->GetNedges();     // Nedges for that element
	int n_points = m_fields[0]->GetExp(el)->GetTotPoints(); // Nquadrature points for that element
	Array<OneD, NekDouble> one2D(n_points, 1.0);
	NekDouble Area = m_fields[0]->GetExp(el)->Integral(one2D);
	Array<OneD, NekDouble> Lengths(n_edges);
	for(int edge = 0; edge< n_edges; ++edge)
	  {	
	    int n_pointsEdge = m_fields[0]->GetExp(el)->GetEdgeExp(edge,false)->GetTotPoints(); // Number of Quadrature Points on each edge
	    Array<OneD, NekDouble> one1D(n_pointsEdge, 1.0);
	    NekDouble L1 = m_fields[0]->GetExp(el)->GetEdgeExp(edge,false)->Integral(one1D);
	    NekDouble L2 = 2.0*Area/L1;
	    Lengths[edge] = (L1<L2) ? L1 : L2;
	  }
	MinLength[el] = Vmath::Vmin(n_edges,&Lengths[0],1);
      }
  }

  void EulerEquations::GetTimeStep(const NekDouble CFL,Array<OneD,NekDouble> &MinLength, Array<OneD, Array<OneD,NekDouble> > inarray ,NekDouble &TimeStep)
  {
    int n_elements  = m_fields[0]->GetExpSize();            // number of element in the mesh
    int nvariables = m_fields.num_elements();               // Number of variables in the mesh
    int nTraceNumTotPoints = GetTraceTotPoints();           // Total Number of Quadrature Points on the edges
    NekDouble gamma = m_gamma;
    int nTotQuadPoints  = GetTotPoints();

    // get temporary array of phys variables
    Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
    //Array<OneD, Array<OneD, NekDouble> > inarray(nvariables);
    Array<OneD, Array<OneD, NekDouble> > PhysVar(nvariables);
    Array<OneD, NekDouble> MinParameter(n_elements);
    NekDouble tmpV, tmpP, tmpR;

    // Fill physical arrays on the traces
    for (int i = 0; i < nvariables; ++i)
      {
	PhysVar[i] = Array<OneD, NekDouble>(nTraceNumTotPoints);
	physarray[i] = Array<OneD, NekDouble>(nTotQuadPoints);;
	//inarray[i] = m_fields[i]->UpdateCoeffs();	
	m_fields[i]->BwdTrans(inarray[i],physarray[i]);
	m_fields[i]->ExtractTracePhys(physarray[i],PhysVar[i]); //QuadExp.cpp GetEdgePhysVals
      }
    
    // Loop on the elements
    for(int el = 0; el < n_elements; ++el)
      {
	int n_edges = m_fields[0]->GetExp(el)->GetNedges();     // Nedges for that element
	Array<OneD, NekDouble> MaxEigenvalue(n_edges);
	for(int edge = 0; edge< n_edges; ++edge)                // Loop over the edges of that el
	  {
	    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
	    int n_pointsEdge = m_fields[0]->GetExp(el)->GetEdgeNumPoints(edge);
	    for (int i = 0; i < nvariables; ++i)
	      {
		Fwd[i] = Array<OneD, NekDouble>(n_pointsEdge);
	      }

	    int traceId = m_fields[0]->GetTraceMap()->GetElmtToTrace()[el][edge]->GetElmtId();
	    int offset = m_fields[0]->GetTrace()->GetPhys_Offset(traceId);
	    Array<OneD, NekDouble> pressure(n_pointsEdge);

	    // Fill physical array for the considered edge
	    for (int i = 0; i < nvariables; ++i)
	      {
		for(int j=0; j<n_pointsEdge; j++)
		  {
		    Fwd[i][j] = PhysVar[i][offset+j];
		  }
	      } 
	    for (int i = 0; i< n_pointsEdge; ++i)
	      {
		pressure[i] = (gamma - 1.0)*(Fwd[3][i] - 0.5*(Fwd[1][i]*Fwd[1][i]/Fwd[0][i] + Fwd[2][i]*Fwd[2][i]/Fwd[0][i]));
	      }
	    Array<OneD, NekDouble> normals(2*n_pointsEdge); // Nx(n_pointsEdge) - Ny(n_pointsEdge)
	    normals = m_fields[0]->GetExp(el)->GetEdgeExp(edge,true)->GetPhysNormals();
	    
	    tmpV = 0.0;
	    tmpP = 0.0;
	    tmpR = 0.0;
	    for (int j = 0; j < n_pointsEdge; ++j)
	      {
		tmpV += (-Fwd[1][j]*normals[j+n_pointsEdge]+Fwd[2][j]*normals[j])/Fwd[0][j];
		tmpP += pressure[j];
		tmpR += Fwd[0][j];
	      }
	    double EdgeNormVelocity = ((tmpV > 0)? tmpV: -tmpV)/n_pointsEdge;
	    double tmp = tmpP/tmpR;
	    NekDouble EdgeSpeedSound;
	    if(tmp >= 0.0)
	      EdgeSpeedSound = sqrt(gamma * tmp);
	    else
	      {
		cout << tmpR/n_pointsEdge << " " << tmpV/n_pointsEdge << " " << tmpP/n_pointsEdge << endl;
		ASSERTL0(false,"Negative Pressure in EulerSolver::GetTimeStep");
		break;
	      }
	    MaxEigenvalue[edge] = EdgeNormVelocity+EdgeSpeedSound;
	  }
	
	MinParameter[el] = MinLength[el]/(Vmath::Vmax(n_edges,&MaxEigenvalue[0],1));
      }
    NekDouble CFLParameter = Vmath::Vmin(n_elements,&MinParameter[0],1);
    TimeStep = CFL*CFLParameter;
  }
  

  
  void EulerEquations::Summary(std::ostream &out)
  {
    cout << "=======================================================================" << endl;
    cout << "\tEquation Type   : Compressible Euler Equations" << endl;
    ADRBase::Summary(out);
    cout << "\tTime max        : " << m_fintime   << endl;
    cout << "\tChecktime       : " << m_timecheck << endl;
    cout << "\tUpwind Flux     : " << UpwindTypeMap[m_upwindType] << endl;
    cout << "\tProblem Type    : " << ProblemTypeMap[m_problemType] << endl;
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
  
  void EulerEquations::GetExactRinglebFlow(Array<OneD, NekDouble> &outarray, int field)
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
    NekDouble theta    = M_PI/4.0;
    NekDouble kExt     = 0.7;
    NekDouble V        = kExt*sin(theta);
    NekDouble dV,dtheta;
    NekDouble toll     = 1.0e-8;
    NekDouble errV     = 1.0;
    NekDouble errTheta = 1.0;
    NekDouble gamma = m_gamma;
    NekDouble gamma_1_2 = (gamma-1.0)/2.0;
    NekDouble c,k,phi,r,J,VV,pp,sint,P,ss;
    NekDouble J11,J12,J21,J22,det;
    NekDouble Fx, Fy;
    NekDouble xi,yi;
    NekDouble par1;
    
    for (int i = 0; i < nTotQuadPoints; ++i)
      {
	while((abs(errV) > toll) || (abs(errTheta)>toll))
	  {
	    VV = V*V;
	    sint = sin(theta);
	    c = sqrt(1.0-gamma_1_2*VV);
	    k = V/sint;
	    phi = 1.0/k;
	    pp = phi*phi;
	    J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
	    r = pow(c,1.0/gamma_1_2);
	    xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
	    yi = phi/(r*V)*sqrt(1.0-VV*pp);
	    par1 = 25.0-5.0*VV;
	    ss = sint*sint;

	    Fx = xi - x[i];
	    Fy = yi - y[i];
	    
	    J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
	    J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
	    J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);
	    // the matrix is singular when theta = pi/2
	    if(abs(y[i])<toll && abs(cos(theta))<toll)
	      {
		J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));

		// dV = -dV/dx * Fx
		dV = -1.0/J22*Fx;
		dtheta = 0.0;
		theta = M_PI/2;
	      }
	    else
	      {
		J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
		det = -1.0/(J11*J22-J12*J21);
		
		// [dV dtheta]' = -[invJ]*[Fx Fy]'
		dV     = det*( J22*Fx-J12*Fy);
		dtheta = det*(-J21*Fx+J11*Fy);
	      }
	    
	    V = V + dV;
	    theta = theta + dtheta;
	    
	    errV     = abs(dV);
	    errTheta = abs(dtheta);

	  }

	c = sqrt(1.0-gamma_1_2*VV);
	r = pow(c,1.0/gamma_1_2);
	
	rho[i] = r;
	rhou[i] = rho[i] * V * cos(theta);
	rhov[i] = rho[i] * V * sin(theta);
	P = (c*c)*rho[i]/gamma;
	E[i]    = P/(gamma-1.0) + 0.5*(rhou[i]*rhou[i]/rho[i]+rhov[i]*rhov[i]/rho[i]);

	// Resetting the guess value
	errV = 1.0;
	errTheta = 1.0;
	theta    = M_PI/4.0;
	V        = kExt*sin(theta);
      }

    /*  */

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
    case 4:
      {
	m_fields[0]->SetPhys(rho);
	m_fields[1]->SetPhys(rhou);
	m_fields[2]->SetPhys(rhov);
	m_fields[3]->SetPhys(E);
	
	// forward transform to fill the modal coeffs
	for(int i = 0; i < m_fields.num_elements(); ++i)
	  {
	    m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
	  }
	
	// dump analytical solution conditions to file
	std::string outname = m_sessionName + "_RinglebFlow_analytical.chk";
	ofstream outfile(outname.c_str());
	WriteFld(outfile);
	
	// Extracting primitive variables on the wall and Mach in the field
	int nq         = m_fields[0]->GetTotPoints();
	int nvariables = m_fields.num_elements();
	Array<OneD, NekDouble> mach(nq);
	Array<OneD, NekDouble> pressure(nq);
	Array<OneD, NekDouble> soundspeed(nq);
	Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
	
	for (int i = 0; i < nvariables; ++i)
	  {
	    physarray[i] = Array<OneD, NekDouble>(nq);
	    m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),physarray[i]);
	  }
	GetPressure(physarray,pressure);
	GetSoundSpeed(physarray,pressure,soundspeed);
	GetMach(physarray,soundspeed,mach);
	Derived_field[0] = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*m_fields[0]);
	// Writing the binary Mach file 
	WriteVar(-1,Derived_field,mach,"Mach");
      }
      break;
    }
  }

  void EulerEquations::SetBoundaryRinglebFlow(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {
    int nvariables      = physarray.num_elements();
    int nTraceNumPoints = GetTraceTotPoints();
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);

    // get physical values of the forward trace (from exp to phys)
    for (int i = 0; i < nvariables; ++i)
      {
	Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
	m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
      }

    for(int e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	
	int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	int id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	int id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));


	Array<OneD,NekDouble> x0(npoints,0.0);
	Array<OneD,NekDouble> x1(npoints,0.0);
	Array<OneD,NekDouble> x2(npoints,0.0); 
	
	m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);

	//---------------------------------
	// flow parameters
	NekDouble theta    = M_PI/4.0;
	NekDouble kExt     = 0.7;
	NekDouble V        = kExt*sin(theta);
	NekDouble dV,dtheta;
	NekDouble toll     = 1.0e-8;
	NekDouble errV     = 1.0;
	NekDouble errTheta = 1.0;
	NekDouble gamma = m_gamma;
	NekDouble gamma_1_2 = (gamma-1.0)/2.0;
	NekDouble c,k,phi,J,VV,pp,sint,P,ss;
	NekDouble J11,J12,J21,J22,det;
	NekDouble Fx, Fy;
	NekDouble xi,yi;
	NekDouble par1;
	NekDouble r;
	
	// Loop on all the points of that edge
	for(int j = 0; j < npoints; j++)
	  {
	    
	    while((abs(errV) > toll) || (abs(errTheta)>toll))
	      {
		VV = V*V;
		sint = sin(theta);
		c = sqrt(1.0-gamma_1_2*VV);
		k = V/sint;
		phi = 1.0/k;
		pp = phi*phi;
		J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
		r = pow(c,1.0/gamma_1_2);
		xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
		yi = phi/(r*V)*sqrt(1.0-VV*pp);
		par1 = 25.0-5.0*VV;
		ss = sint*sint;
		
		Fx = xi - x0[j];
		Fy = yi - x1[j];
		
		J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
		J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
		J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);
		// the matrix is singular when theta = pi/2
		if(abs(x1[j])<toll && abs(cos(theta))<toll)
		  {
		    J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
		    
		    // dV = -dV/dx * Fx
		    dV = -1.0/J22*Fx;
		    dtheta = 0.0;
		    theta = M_PI/2;
		  }
		else
		  {
		    J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
		    det = -1.0/(J11*J22-J12*J21);
		    
		    // [dV dtheta]' = -[invJ]*[Fx Fy]'
		    dV     = det*( J22*Fx-J12*Fy);
		    dtheta = det*(-J21*Fx+J11*Fy);
		  }
		
		V = V + dV;
		theta = theta + dtheta;
		
		errV     = abs(dV);
		errTheta = abs(dtheta);
		
	      }
	    
	    c  = sqrt(1.0-gamma_1_2*VV);
	    int kk = id2+j;
	    Fwd[0][kk]  = pow(c,1.0/gamma_1_2);
	    Fwd[1][kk] = Fwd[0][kk] * V * cos(theta);
	    Fwd[2][kk] = Fwd[0][kk] * V * sin(theta);
	    P  = (c*c)*Fwd[0][kk]/gamma;
	    Fwd[3][kk]  = P/(gamma-1.0) + 0.5*(Fwd[1][kk]*Fwd[1][kk]/Fwd[0][kk]+Fwd[2][kk]*Fwd[2][kk]/Fwd[0][kk]);

	    errV = 1.0;
	    errTheta = 1.0;
	    theta    = M_PI/4.0;
	    V        = kExt*sin(theta);
	    
	  }
	for (int i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(npoints,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);	
	  }
      }
  }   

  void EulerEquations::SetInitialRinglebFlow(void)
  {

    int nbnd    = m_fields[0]->GetBndConditions().num_elements();// Get number of different boundaries in the input file

    // Loop on all the edges of the input file
    for(int bcRegion=0; bcRegion < nbnd; ++bcRegion)
      {

	    int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetNpoints();
	    Array<OneD,NekDouble> x0(npoints,0.0);
	    Array<OneD,NekDouble> x1(npoints,0.0);
	    Array<OneD,NekDouble> x2(npoints,0.0); 

	    Array<OneD, NekDouble> rho(npoints,0.0);
	    Array<OneD, NekDouble> rhou(npoints,0.0);
	    Array<OneD, NekDouble> rhov(npoints,0.0);
	    Array<OneD, NekDouble> E(npoints,0.0);
	
	    m_fields[0]->GetBndCondExpansions()[bcRegion]->GetCoords(x0,x1,x2);
	    //---------------------------------
	    // flow parameters
	    NekDouble theta    = M_PI/4.0;
	    NekDouble kExt     = 0.7;
	    NekDouble V        = kExt*sin(theta);
	    NekDouble dV,dtheta;
	    NekDouble toll     = 1.0e-8;
	    NekDouble errV     = 1.0;
	    NekDouble errTheta = 1.0;
	    NekDouble gamma = m_gamma;
	    NekDouble gamma_1_2 = (gamma-1.0)/2.0;
	    NekDouble c,k,phi,J,VV,pp,sint,P,ss;
	    NekDouble J11,J12,J21,J22,det;
	    NekDouble Fx, Fy;
	    NekDouble xi,yi;
	    NekDouble par1;
	    NekDouble r;
	    
	    // Loop on all the points of that edge
	    for(int j = 0; j < npoints; j++)
	      {
		
		while((abs(errV) > toll) || (abs(errTheta)>toll))
		  {
		    
		    VV = V*V;
		    sint = sin(theta);
		    c = sqrt(1.0-gamma_1_2*VV);
		    k = V/sint;
		    phi = 1.0/k;
		    pp = phi*phi;
		    J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
		    r = pow(c,1.0/gamma_1_2);
		    xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
		    yi = phi/(r*V)*sqrt(1.0-VV*pp);
		    par1 = 25.0-5.0*VV;
		    ss = sint*sint;
		    
		    Fx = xi - x0[j];
		    Fy = yi - x1[j];
		    
		    J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
		    J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
		    J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);
		    
		    // the matrix is singular when theta = pi/2
		    if(abs(x1[j])<toll && abs(cos(theta))<toll)
		      {
			
			J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
			
			// dV = -dV/dx * Fx
			dV = -1.0/J22*Fx;
			dtheta = 0.0;
			theta = M_PI/2;
		      }
		    else
		      {
			
			J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
			det = -1.0/(J11*J22-J12*J21);
			
			// [dV dtheta]' = -[invJ]*[Fx Fy]'
			dV     = det*( J22*Fx-J12*Fy);
			dtheta = det*(-J21*Fx+J11*Fy);
		      }
		    
		    V = V + dV;
		    theta = theta + dtheta;
		    
		    errV     = abs(dV);
		    errTheta = abs(dtheta);
		    
		  }
		
		c  = sqrt(1.0-gamma_1_2*VV);
		rho[j]  = pow(c,1.0/gamma_1_2);
		rhou[j] = rho[j] * V * cos(theta);
		rhov[j] = rho[j] * V * sin(theta);
		P  = (c*c)*rho[j]/gamma;
		E[j]  = P/(gamma-1.0) + 0.5*(rhou[j]*rhou[j]/rho[j]+rhov[j]*rhov[j]/rho[j]);
		
		errV = 1.0;
		errTheta = 1.0;
		theta    = M_PI/4.0;
		V        = kExt*sin(theta);
		
	      }

	    m_fields[0]->GetBndCondExpansions()[bcRegion]->SetPhys(rho);
	    m_fields[1]->GetBndCondExpansions()[bcRegion]->SetPhys(rhou);
	    m_fields[2]->GetBndCondExpansions()[bcRegion]->SetPhys(rhov);
	    m_fields[3]->GetBndCondExpansions()[bcRegion]->SetPhys(E);

	// forward transform to fill the modal coeffs
	for(int i = 0; i < m_fields.num_elements(); ++i)
	  {
	    m_fields[i]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[i]->GetBndCondExpansions()[bcRegion]->GetPhys(),m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
	  }

      }

    // Calculation of the initial internal values as a weighted average over the distance from the boundaries
    int nq = m_fields[0]->GetNpoints();
    
    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);
      
    // get the coordinates (assuming all fields have the same discretisation)
    m_fields[0]->GetCoords(x0,x1,x2);
    
    for(int j = 0; j < nq; j++)
      {
	NekDouble Dist = 0.0; 
	NekDouble rho  = 0.0;
	NekDouble rhou = 0.0;
	NekDouble rhov = 0.0;
	NekDouble E    = 0.0;	   
	NekDouble SumDist = 0.0;
	
	// Calculation of all the distances
	// Loop on all the edges of the input file
	for(int bcRegion=0; bcRegion < nbnd; ++bcRegion)
	  {
	    
	    int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetNpoints(); // get quadrature points on the edge
	    
	    Array<OneD,NekDouble> xb0(npoints,0.0);
	    Array<OneD,NekDouble> xb1(npoints,0.0);
	    Array<OneD,NekDouble> xb2(npoints,0.0); 
	    
	    m_fields[0]->GetBndCondExpansions()[bcRegion]->GetCoords(xb0,xb1,xb2);
	    
	    for(int k=0; k<npoints; k++)
	      {
		Dist = sqrt((xb0[k]-x0[j])*(xb0[k]-x0[j])+(xb1[k]-x1[j])*(xb1[k]-x1[j])+(xb2[k]-x2[j])*(xb2[k]-x2[j]));
		SumDist += Dist;
		rho  += Dist*(m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
		rhou += Dist*(m_fields[1]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
		rhov += Dist*(m_fields[2]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
		E    += Dist*(m_fields[3]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
	      }
	  }
	
	rho  = rho/SumDist;
	rhou = rhou/SumDist;
	rhov = rhov/SumDist;
	E    = E/SumDist;

	(m_fields[0]->UpdatePhys())[j] = rho;
	(m_fields[1]->UpdatePhys())[j] = rhou;
	(m_fields[2]->UpdatePhys())[j] = rhov;
	(m_fields[3]->UpdatePhys())[j] = E;

      }

    for(int i = 0 ; i < m_fields.num_elements(); i++)
      {
	m_fields[i]->SetPhysState(true);
	m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
      }
    
    // dump initial conditions to file
    std::string outname = m_sessionName + "_initialRingleb.chk";
    ofstream outfile(outname.c_str());
    WriteFld(outfile);

  }
  
} //end of namespace

/**
* $Log: EulerEquations.cpp,v $
* Revision 1.5  2009/06/29 07:47:33  claes
* bug fix in WallBoundary
*
* Revision 1.4  2009/04/28 10:17:41  pvos
* Some updates to make the solvers compile properly with the newly added sparse matrix library
*
* Revision 1.3  2009/03/17 12:32:01  claes
* Updates to get the EulerSolver to work with the latest version of the TimeIntegrationScheme
*
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
