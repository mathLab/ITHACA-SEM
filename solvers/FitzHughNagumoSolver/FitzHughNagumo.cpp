///////////////////////////////////////////////////////////////////////////////
//
// File FitzHughNagumo.cpp
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
// Description: Advection Diffusion Reaction class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <FitzHughNagumoSolver/FitzHughNagumo.h>
#include <cstdio>
#include <cstdlib>
namespace Nektar
{
    /**
     * Basic construnctor
     */
    FitzHughNagumo::FitzHughNagumo(void):
        ADRBase(),
        m_infosteps(100),
        m_explicitAdvection(true),
        m_explicitDiffusion(true),
        m_explicitReaction(true)
    {     
    }

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    FitzHughNagumo::FitzHughNagumo(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10),
        m_explicitDiffusion(true),
        m_explicitReaction(true)
    {

        int i;

        if(m_boundaryConditions->CheckForParameter("epsilon") == true)
        {
            m_epsilon = m_boundaryConditions->GetParameter("epsilon");
        }
        else
        {
            m_epsilon = 0.20;
        }

        if(m_boundaryConditions->CheckForParameter("beta") == true)
        {
            m_beta = m_boundaryConditions->GetParameter("beta");
        }
        else
        {
            m_beta  = 0.77;
        }

        // initialwave type
        // 0 = plane wave propagation from the left
        // 1 = plane wave propagation from the bottom
        // 2 = A circular wave from the corner (topright)
        // 3 = Point initialization from the center

        if(m_boundaryConditions->CheckForParameter("initialwavetype") == true)
        {
            m_initialwavetype = m_boundaryConditions->GetParameter("initialwavetype");
        }
        else
        {
            m_initialwavetype  = 0;
        }

        if(m_boundaryConditions->CheckForParameter("timedelay") == true)
        {
            m_timedelay = m_boundaryConditions->GetParameter("timedelay");
        }
        else
        {
            m_timedelay  = 0.0;
        }

        if(m_boundaryConditions->CheckForParameter("duration") == true)
        {
            m_duration = m_boundaryConditions->GetParameter("duration");
        }
        else
        {
            m_duration  = 3.5;
        }

        if(m_boundaryConditions->CheckForParameter("frequency1") == true)
        {
            m_frequency1 = m_boundaryConditions->GetParameter("frequency1");
        }
        else
        {
            m_frequency1  = 5.0;
        }

        if(m_boundaryConditions->CheckForParameter("frequency2") == true)
        {
            m_frequency2 = m_boundaryConditions->GetParameter("frequency2");
        }
        else
        {
            m_frequency2  = 5.0;
        }

        if(m_boundaryConditions->CheckForParameter("x1center") == true)
        {
            m_x1center = m_boundaryConditions->GetParameter("x1center");
        }
        else
        {
            m_x1center  = 0.0;
        }

        if(m_boundaryConditions->CheckForParameter("y1center") == true)
        {
            m_y1center = m_boundaryConditions->GetParameter("y1center");
        }
        else
        {
            m_y1center  = 0.0;
        }

        if(m_boundaryConditions->CheckForParameter("x2center") == true)
        {
            m_x2center = m_boundaryConditions->GetParameter("x2center");
        }
        else
        {
            m_x2center  = 0.0;
        }

        if(m_boundaryConditions->CheckForParameter("y2center") == true)
        {
            m_y2center = m_boundaryConditions->GetParameter("y2center");
        }
        else
        {
            m_y2center  = 0.0;
        }

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(NoCaseStringCompare(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }

        ASSERTL0(i != (int) eEquationTypeSize, "Invalid expansion type.");
 
	std::string Implicit = "Implicit"; 
	if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
	  {
	    m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
	  }
	
	// check that any user defined boundary condition is indeed implemented
	for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
	  {	
	    // Time Dependent Boundary Condition (if no use defined then this is empty)
	    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "")
	      {
		if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "TimeDependent")
		  {
		    ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
		  }
	      }
	  }
	
	// Check for definition of Implicit/Explicit terms in solverinfo
	if(m_boundaryConditions->SolverInfoExists("ADVECTIONADVANCEMENT"))
	  {
	    std::string AdvStr = m_boundaryConditions->GetSolverInfo("ADVECTIONADVANCEMENT");
            
	    if(NoCaseStringCompare(AdvStr,Implicit) == 0)
	      {
		m_explicitAdvection = false;
	      }
	    else
	      {
		m_explicitAdvection = true;
	      }
	  }
	else
	  {
	    m_explicitAdvection = true;
	  }
	
        
	if(m_boundaryConditions->SolverInfoExists("DIFFUSIONADVANCEMENT"))
	  {
	    std::string AdvStr = m_boundaryConditions->GetSolverInfo("DIFFUSIONADVANCEMENT");
            
	    if(NoCaseStringCompare(AdvStr,Implicit) == 0 )
	      {
		m_explicitDiffusion = false;
		// Reset default for implicit diffusion
	      }
	    else
	      {
		m_explicitDiffusion = true;
	      }
	  }
	else
	  {
	    m_explicitDiffusion = true;
	  }
	
	if(m_boundaryConditions->SolverInfoExists("REACTIONADVANCEMENT"))
	  {
	    std::string AdvStr = m_boundaryConditions->GetSolverInfo("REACTIONADVANCEMENT");
            
	    if(NoCaseStringCompare(AdvStr,Implicit) == 0)
	      {
		m_explicitReaction = false;
	      }
	    else
	      {
		m_explicitReaction = true;
	      }
	  }
	else
	  {
	    m_explicitReaction = true;
	  }

        m_timeIntMethod = LibUtilities::eIMEXdirk_3_4_3;	


	// check to see if time stepping has been reset
	if(m_boundaryConditions->SolverInfoExists("TIMEINTEGRATIONMETHOD"))
	  {
	    std::string TimeIntStr = m_boundaryConditions->GetSolverInfo("TIMEINTEGRATIONMETHOD");
	    int i;
	    for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
	      {
		if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
		  {
		    m_timeIntMethod = (LibUtilities::TimeIntegrationMethod)i; 
		    break;
		  }
	      }
	    
	    ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
	  }
    }

    void FitzHughNagumo::SetFHNInitialConditions(const int initialwavetype, NekDouble initialtime)
    {
        int nq = m_fields[0]->GetNpoints();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        NekDouble unew=100.0, uinit=0.0, vinit, f, fd, fval, gval,Tol=0.00000001;

        // get the coordinates (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        // Get the initial constant u and v
        for(int j = 0; j<1000; j++)
        {
            // f = fvalue(uinit);
            f = uinit*uinit*uinit + 3.0*uinit + 6.0*m_beta;

            // fd = fderiv(uinit);
            fd = 3.0*uinit*uinit + 3.0;

            unew = uinit - f/fd;

            if(abs(unew-uinit) < Tol)
            {
                break;
            }
            uinit = unew;
        }
        vinit = uinit - (1.0/3.0)*uinit*uinit*uinit;

        cout << "Static u0 = " << uinit << ", Static v0 = " << vinit << endl;

        NekDouble xmin = Vmath::Vmin(nq,x0,1);
        NekDouble xmax = Vmath::Vmax(nq,x0,1);
        NekDouble ymin = Vmath::Vmin(nq,x1,1);
        NekDouble ymax = Vmath::Vmax(nq,x1,1);
        NekDouble rad;
        NekDouble xc = 0.5*(xmin+xmax);
        NekDouble yc = 0.5*(ymin+ymax);

        switch(initialwavetype)
        {
        // plane wave from the bottom
        case(1):
            {
                for(int j = 0; j < nq; j++)
                {
                    if( x1[j] <= (ymin+m_duration) )
                    {
                        (m_fields[0]->UpdatePhys())[j] = 2.0;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                    else
                    {
                        (m_fields[0]->UpdatePhys())[j] = uinit;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                }   
            }            
            
            // Circular wave from the top-right corner
        case(2):
            {
                for(int j = 0; j < nq; j++)
                {
                    rad = sqrt( (x0[j]-xmax)*(x0[j]-xmax) + (x1[j]-ymax)*(x1[j]-ymax) );
                    
                    if( rad <= m_duration )
                    {
                        (m_fields[0]->UpdatePhys())[j] = 2.0;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                    
                    else
                    {
                        (m_fields[0]->UpdatePhys())[j] = uinit;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                }
            }

            // Point initialization at 
        case(3):
            {
                for(int j = 0; j < nq; j++)
                {
                    rad = sqrt( (x0[j]-m_x1center)*(x0[j]-m_x1center) + (x1[j]-m_y1center)*(x1[j]-m_y1center) );
                    
                    if( rad <= m_duration )
                    {
                        (m_fields[0]->UpdatePhys())[j] = 2.0;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                    
                    else
                    {
                        (m_fields[0]->UpdatePhys())[j] = uinit;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                }
            }
                       
            // otherwise planar wave from the left
        default:
            {
                for(int j = 0; j < nq; j++)
                {
                    if( x0[j] <= (xmin+m_duration) )
                    {
                        (m_fields[0]->UpdatePhys())[j] = 2.0;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                    else
                    {
                        (m_fields[0]->UpdatePhys())[j] = uinit;
                        (m_fields[1]->UpdatePhys())[j] = vinit;
                    }
                }   
            }
            
        }
        
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	  {
              m_fields[i]->SetPhysState(true);
              m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
          }

	// dump initial conditions to file
	for(int i = 0; i < m_fields.num_elements(); ++i)
	  {
            // dump initial conditions to file
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
	  }       
    }

    void FitzHughNagumo::Generatesecondstimulus(const int initialwavetype, 
                                                Array<OneD, NekDouble>&outarray)
    {
        int nvar = m_fields.num_elements();
        int nq = m_fields[0]->GetNpoints();
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        NekDouble unew=100.0, uinit=0.0, vinit, f, fd, fval, gval,Tol=0.00000001;

        // get the coordinates (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        Array<OneD, NekDouble> physfield(nq,0.0); 

        NekDouble xmin = Vmath::Vmin(nq,x0,1);
        NekDouble xmax = Vmath::Vmax(nq,x0,1);
        NekDouble ymin = Vmath::Vmin(nq,x1,1);
        NekDouble ymax = Vmath::Vmax(nq,x1,1);
        NekDouble rad;

        switch(initialwavetype)
        {
            // Plane wave from the left
        case(0):
            {
                for(int j = 0; j < nq; j++)
                {
                    if( x0[j] <= (xmin+m_duration) )
                    {
                        physfield[j] = 2.0;
                    }
                }   
            }

            // Plane wave from the bottom
        case(1):
            {
                for(int j = 0; j < nq; j++)
                {
                    if( x1[j] <= (ymin+m_duration) )
                    {
                        physfield[j] = 2.0;
                    }
                }   
            }
            
            // Circular wave from the corner
        case(2):
            {
                for(int j = 0; j < nq; j++)
                {
                    rad = sqrt( (x0[j]-xmax)*(x0[j]-xmax) + (x1[j]-ymax)*(x1[j]-ymax) );
                    
                    if( rad <= m_duration )
                    {
                        physfield[j] = 2.0;
                    }
                }   
            }
            
            // Point initialization at (m_x2center, m_y2center)
        case(4):
            {
                for(int j = 0; j < nq; j++)
                {
                    rad = sqrt( (x0[j]-m_x2center)*(x0[j]-m_x2center) + (x1[j]-m_y2center)*(x1[j]-m_y2center) );
                    
                    if( rad <= m_duration )
                    {
                        physfield[j] = 2.0;
                    }
                }   
            }
        }
        
        m_fields[0]->FwdTrans(physfield,outarray);
    }

    void FitzHughNagumo::ODEeReactionIMEXtest(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
                                              Array<OneD, Array<OneD, NekDouble> >&outarray, 
                                              const NekDouble time)
	
    {
        NekDouble PI = 3.141592653589793238462;
        
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();        
        int npoints = m_fields[0]->GetNpoints();
        
        Array<OneD,NekDouble> x0(npoints,0.0);
        Array<OneD,NekDouble> x1(npoints,0.0);
        Array<OneD,NekDouble> x2(npoints,0.0);  
        
        Array<OneD,NekDouble> physfield(npoints);
        Array<OneD,NekDouble> ft(npoints,0.0);
        
        NekDouble csx,csy;
       
        //Compute f(u) in physical field
        m_fields[0]->BwdTrans(inarray[0],physfield);
	m_fields[0]->SetPhysState(true);        
        
        Vmath::Vmul(npoints,physfield, 1, physfield, 1, physfield, 1);
        
        m_fields[0]->GetCoords(x0,x1,x2);
        
        for (int i =0; i< npoints; ++i)
        {
            csx = cos(PI*x0[i]);
            csy = cos(PI*x1[i]);
            ft[i] = exp(-2.0*time)*csx*csx*csy*csy + (2*PI*PI-1.0)*exp(-1.0*time)*csx*csy; 
        }
        
        Vmath::Vsub(npoints, ft, 1, physfield, 1, physfield, 1);
        
        // Back to modal field
        m_fields[0]->FwdTrans(physfield,outarray[0]);
      	m_fields[0]->SetPhysState(false);        
    }
    
    void FitzHughNagumo::ODEeReactiontest1(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
					   Array<OneD, Array<OneD, NekDouble> >&outarray, 
					   const NekDouble time)
  
    {
      int nvariables = inarray.num_elements();
      int ncoeffs    = inarray[0].num_elements();
      int npoints = m_fields[0]->GetNpoints();
      
      const NekDouble coeff = 2.0/m_epsilon;
      
      Array<OneD, NekDouble> physfield(npoints);
      Array<OneD, NekDouble> utemp2(npoints,0.0);
      Array<OneD, NekDouble> utemp3(npoints,0.0);

      for (int i = 0; i < nvariables; ++i)
	{  
	  m_fields[i]->BwdTrans(inarray[i],physfield);
          m_fields[i]->SetPhysState(true);        
	  
	  // temp3 = (2/epsilon)*(u*u - u*u*u)
	  Vmath::Vmul(npoints, physfield, 1, physfield, 1, utemp2, 1);
	  Vmath::Vmul(npoints, physfield, 1, utemp2, 1, utemp3, 1);
	  Vmath::Vsub(npoints, utemp2, 1, utemp3, 1, physfield, 1);
	  Vmath::Smul(npoints, coeff, physfield, 1, physfield, 1);
	  
	  m_fields[i]->FwdTrans(physfield,outarray[i]);
          m_fields[i]->SetPhysState(false);        
        }
    }

    void FitzHughNagumo::ODEeReactiontest2(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
					   Array<OneD, Array<OneD, NekDouble> >&outarray, 
					   const NekDouble time)
	
    {
      int nvariables = inarray.num_elements();
      int ncoeffs    = inarray[0].num_elements();
      int npoints = m_fields[0]->GetNpoints();
  
      Array<OneD,NekDouble> x0(npoints,0.0);
      Array<OneD,NekDouble> x1(npoints,0.0);
      Array<OneD,NekDouble> x2(npoints,0.0);  

      const NekDouble coeff = 2.0/m_epsilon;
      const NekDouble Tol = 0.0000001;
      
      Array<OneD, NekDouble> physfield(npoints);
      Array<OneD, NekDouble> dist(npoints,0.0);
      Array<OneD, NekDouble> temp(npoints,0.0);
      
      for (int i = 0; i < nvariables; ++i)
	{  
	  m_fields[i]->BwdTrans(inarray[i],physfield);
          m_fields[0]->SetPhysState(true);        
	  
	  // temp = u - u*u
	  Vmath::Vmul(npoints, physfield, 1, physfield, 1, temp, 1);
	  Vmath::Vsub(npoints, physfield, 1, temp, 1, temp, 1);
	  
	  m_fields[0]->GetCoords(x0,x1,x2);
	  
	  for (int j =0; j< npoints; ++j)
	    {  
	      dist[j] = 1.0/sqrt(x0[j]*x0[j]+x1[j]*x1[j]+Tol);
	    }
	  
	  // v = (2/epsilon)*u + 1/|x|
	  Vmath::Svtvp(npoints, coeff, physfield, 1, dist, 1, physfield, 1);
	  
	  // f(u) = temp*v = u(1-u)*( (2/epsilon)*u + 1/|x| )
	  Vmath::Vmul(npoints, temp, 1, physfield, 1, physfield, 1);
	  
	  m_fields[i]->FwdTrans(physfield,outarray[i]);
          m_fields[0]->SetPhysState(false);        
	}
    }

    void FitzHughNagumo::ODEeReactionmono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
					  Array<OneD, Array<OneD, NekDouble> >&outarray, 
					  const NekDouble time)
	
    {

      NekDouble m_gamma = 0.5;

      int nvariables = inarray.num_elements();
      int ncoeffs    = inarray[0].num_elements();
      int npoints    = m_fields[0]->GetNpoints();
      
      Array<OneD, NekDouble> physfield(npoints);
      Array<OneD, NekDouble> temp2(npoints,0.0);
      Array<OneD, NekDouble> temp3(npoints,0.0);
      
      Array<OneD, NekDouble> temp(ncoeffs, 0.0);
      
      m_fields[0]->BwdTrans(inarray[0],physfield);
      m_fields[0]->SetPhysState(true);        

      // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
      // physfield = u - (1.0/3.0)*u*u*u
      Vmath::Vmul(npoints, physfield, 1, physfield, 1, temp2, 1);
      Vmath::Vmul(npoints, physfield, 1, temp2, 1, temp3, 1);
      Vmath::Svtvp(npoints, (-1.0/3.0), temp3, 1, physfield, 1, physfield, 1);
      
      m_fields[0]->FwdTrans(physfield,outarray[0]);
      m_fields[0]->SetPhysState(false);        
     
      Vmath::Vsub(ncoeffs, inarray[1], 1, outarray[0], 1, outarray[0], 1);
      Vmath::Smul(ncoeffs, -1.0/m_epsilon, outarray[0], 1, outarray[0], 1);
      
      // For v: m_epsilon*( u + m_beta - m_gamma*v )
      Vmath::Svtvp(ncoeffs, -1.0*m_gamma, inarray[1], 1, inarray[0], 1, outarray[1], 1);
      Vmath::Sadd(ncoeffs, m_beta, outarray[1], 1, outarray[1], 1);
      Vmath::Smul(ncoeffs, m_epsilon, outarray[1], 1, outarray[1], 1);
    }
  
  void FitzHughNagumo::ODEhelmSolvetest(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, NekDouble> >&outarray,
                                        NekDouble time, 
                                        NekDouble lambda)
  {
    int nvariables = inarray.num_elements();
    int ncoeffs    = inarray[0].num_elements();
    
    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}    
    // outarray = output: nabla^2 \hat{Y}       
    // where \hat = modal coeffs

    NekDouble kappa = 1.0/(lambda*m_epsilon);

    for (int i = 0; i < nvariables; ++i)
      {	
	// Multiply rhs[i] with -1.0/gamma/timestep
	Vmath::Smul(ncoeffs, -1.0*kappa, inarray[i], 1, outarray[i], 1);
	
	// Update coeffs to m_fields
	m_fields[i]->UpdateCoeffs() = outarray[i];
			
	// Backward Transformation to nodal coefficients
	m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());

	// Solve a system of equations with Helmholtz solver
	m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),kappa);
	m_fields[i]->SetPhysState(false);        
	
	// The solution is Y[i]
	outarray[i] = m_fields[i]->GetCoeffs();	  
      }
  }


  void FitzHughNagumo::ODEhelmSolvemono(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, NekDouble> >&outarray,
                                        NekDouble time, 
                                        NekDouble lambda)
  {
    int nvariables = inarray.num_elements();
    int ncoeffs    = inarray[0].num_elements();
    
    NekDouble kappa = 1.0/lambda;

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}    
    // outarray = output: nabla^2 \hat{Y}       
    // where \hat = modal coeffs
    
    // Multiply rhs[i] with -1.0/gamma/timestep
    Vmath::Smul(ncoeffs, -1.0*kappa, inarray[0], 1, outarray[0], 1);
    
    // Update coeffs to m_fields
    m_fields[0]->UpdateCoeffs() = outarray[0];
    
    // Backward Transformation to nodal coefficients
    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(), m_fields[0]->UpdatePhys());
    
    // Solve a system of equations with Helmholtz solver
    m_fields[0]->HelmSolve(m_fields[0]->GetPhys(),m_fields[0]->UpdateCoeffs(),kappa);
    m_fields[0]->SetPhysState(false);        
    
    // The solution is Y[i]
    outarray[0] = m_fields[0]->GetCoeffs();	  
        
    // For q: No helmholtz solver is needed=============================
    Vmath::Vcopy(ncoeffs, inarray[1], 1, outarray[1], 1);
  }
  
  
    void FitzHughNagumo::SolveHelmholtz(NekDouble lambda)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),lambda);
            m_fields[i]->SetPhysState(false);
        }
    }

    // For Continuous Galerkin projections with time-dependent dirichlet boundary conditions,
    // the time integration can be done as follows:
    // The ODE resulting from the PDE can be formulated as:
    // 
    // M du/dt = F(u)  or du/dt = M^(-1) F(u)
    //
    // Now suppose that M does not depend of time, the ODE can than be written as:
    //
    // d(Mu)/dt = F(u)
    //
    // Introducing the variable u* = Mu, this yields
    //
    // du*/dt = F( M^(-1) u*  ) = F*(u*)
    //
    // So rather than solving the initial ODE, it is advised to solve this new ODE for u*
    // as this allows for an easier treatment of the dirichlet boundary conditions.
    // However, note that at the end of every time step, the actual solution u can
    // be calculated as:
    // 
    // u = M^(-1) u*;
    //
    // This can be viewed as projecting the solution u* onto the known boundary conditions.
    // Note that this step is also done inside the ODE rhs function F*.
    //
    // In order for all of this to work appropriately, make sure that the operator M^(-1)
    // does include the enforcment of the dirichlet boundary conditionst

  void FitzHughNagumo::GeneralTimeIntegration(int nsteps, 
					      LibUtilities::TimeIntegrationMethod IntMethod,
					      LibUtilities::TimeIntegrationSchemeOperators ode)
  {
    int i,n,nchk = 0;
    int ncoeffs = m_fields[0]->GetNcoeffs();
    int nvariables = m_fields.num_elements();
    
    // Set up wrapper to fields data storage. 
    Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
    Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);
    Array<OneD, NekDouble>   initialimpulse(ncoeffs);
    Array<OneD, NekDouble>   secondimpulse(ncoeffs);
    
    for(i = 0; i < nvariables; ++i)
      {
	m_fields[i]->SetPhysState(false);
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
					          
        NekDouble timenow1, timenow2,Tol=0.00001;
        int nofreq1=0, nofreq2=0;
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

                timenow1 = abs(1.0*m_timestep*n - ( (nofreq1+1)*m_frequency1) );
                if(timenow1<Tol)
                {
                    Generatesecondstimulus(m_secondwavetype,initialimpulse);
                    Vmath::Vadd(ncoeffs, &initialimpulse[0], 1, &fields[0][0], 1, &fields[0][0], 1);
                    nofreq1++;
                }

                timenow2 = abs(1.0*m_timestep*n - ( (nofreq2*m_frequency2)+ m_timedelay) );
                if(timenow2<Tol)
                {
                    Generatesecondstimulus(m_secondwavetype,secondimpulse);
                    Vmath::Vadd(ncoeffs, &secondimpulse[0], 1, &fields[0][0], 1, &fields[0][0], 1);
                    nofreq2++;
                }
            }

            m_time += m_timestep;

            if(m_projectionType==eGalerkin)
            {
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                SetBoundaryConditions(m_time);
                for(i = 0; i < nvariables; ++i)
                {
		  m_fields[i]->SetPhysState(false);

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
	      Checkpoint_Output(nchk++);
            }
        }
        

        for(i = 0; i < nvariables; ++i)
        {
	  (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }
	
    
  //----------------------------------------------------
  void FitzHughNagumo::SetBoundaryConditions(NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->EvaluateBoundaryConditions(time);
    }
  }

  // Compute the fluxes of q from the scalar functin u.
  // Input:   ufield (1 by nTraceNumPoints) - Should be in physical field
  // Output:  ufluxFwd  (2 by nTraceNumPoints) - Flux values for forward edges
  //          ufluxBwd  (2 by nTraceNumPoints) - Flux values for backward edges

    void FitzHughNagumo::NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield, 
						      Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
    {
        int i,j;
        int nTraceNumPoints = GetTraceNpoints();
	int nvariables = m_fields.num_elements();
        int nqvar = uflux.num_elements();
	
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
	Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);
        	  		  
	// Get the sign of (v \cdot n), v = an arbitrary vector

	// Vn = V \cdot n, where n is tracenormal for eForward edges. Set V = (1,0)
	// Vmath::Vcopy(nTraceNumPoints,m_traceNormals_tbasis[0],1,Vn,1);

        //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
        for (j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables ; ++i)
            {
                //  Compute Forward and Backward value of ufield of i direction
                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
                
                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd
                
                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

		m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],Fwd,Bwd,fluxtemp);  
                
                // Imposing weak boundary condition with flux
                // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd
                
                // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd
                if(m_fields[0]->GetBndCondExpansions().num_elements())
		{
                    WeakPenaltyforScalar(i,ufield[i],fluxtemp);
		}
                
                // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )
                
                // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

		Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[j][i],1);
            }
	}
    }

    // Compute the fluxes of q and u vector fields for discontinuous diffusion term
    // Input:   qfield : 2 by # of total trace points
    // Output:  qflux  : 2 by # of total trace points
    void FitzHughNagumo::NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                         Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                                         Array<OneD, Array<OneD, NekDouble> >  &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
	int nvariables = m_fields.num_elements();
        int nqvar = qfield.num_elements();

	NekDouble C11 = 1.0;			
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
			
	Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
	Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

	Array<OneD, NekDouble > uterm(nTraceNumPoints);
	  		  
        // Get the sign of (v \cdot n), v = an arbitrary vector
	// Vn = V \cdot n, where n is tracenormal for eForward edges
	// Vmath::Vcopy(nTraceNumPoints,m_traceNormals[0],1,Vn,1);
        
        // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)		   			
        for(int i = 0; i < nvariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
            for(int j = 0; j < nqvar; ++j)
            {
		//  Compute Forward and Backward value of ufield of jth direction
		m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);
                
                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-
                
                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd =q+

		m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],qBwd,qFwd,qfluxtemp);
		Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);

		// Generate Stability term = - C11 ( u- - u+ )
		m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
		Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);					  
		Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

		//  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
		Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);
		
		// Imposing weak boundary condition with flux
		if(m_fields[0]->GetBndCondExpansions().num_elements())
		  {
		    WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
		  }
		
		// q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
		// n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
		// n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
	    }
        }
    }
    

	
  // Diffusion: Imposing weak boundary condition for u with flux 
  //  uflux = g_D  on Dirichlet boundary condition
  //  uflux = u_Fwd  on Neumann boundary condition
  void FitzHughNagumo::WeakPenaltyforScalar(const int var,
					    const Array<OneD, const NekDouble> &physfield, 
					    Array<OneD, NekDouble> &penaltyflux,
					    NekDouble time)
  {
    int i, j, e, npoints, id1, id2;
    int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
    int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
    int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
    int nTraceNumPoints = GetTraceNpoints();

    Array<OneD, NekDouble > uplus(nTraceNumPoints);

    m_fields[var]->ExtractTracePhys(physfield,uplus);            
    for(i = 0; i < nbnd; ++i)
      {                 
	// Evaluate boundary values g_D or g_N from input files
	SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(i);
	npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();
	
	Array<OneD,NekDouble> BDphysics(npoints);
	Array<OneD,NekDouble> x0(npoints,0.0);
	Array<OneD,NekDouble> x1(npoints,0.0);
	Array<OneD,NekDouble> x2(npoints,0.0);  
        
	m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
	for(j = 0; j < npoints; j++)
	  {
	    BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
	  }
	
	// Weakly impose boundary conditions by modifying flux values
	for (e = 0; e < numBDEdge ; ++e)
	  {
	    id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
	    id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));
	    
	    // For Dirichlet boundary condition: uflux = g_D
	    if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
	      {
                   Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
	      }
	    
	    // For Neumann boundary condition: uflux = u+
	    else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
	      {
		Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
	      }
	  }
      }           
  }
  
  // Diffusion: Imposing weak boundary condition for q with flux 
  //  uflux = g_D  on Dirichlet boundary condition
  //  uflux = u_Fwd  on Neumann boundary condition
    void FitzHughNagumo::WeakPenaltyforVector(const int var,
                                             const int dir,
                                             const Array<OneD, const NekDouble> &physfield,
                                             Array<OneD, NekDouble> &penaltyflux,
                                             NekDouble C11,
                                             NekDouble time)
  {
    int i, j, e, npoints, id1, id2;
    int nbnd = m_fields[0]->GetBndCondExpansions().num_elements();
    int numBDEdge = m_fields[0]->GetBndCondExpansions()[0]->GetExpSize();
    int Nfps = m_fields[0]->GetBndCondExpansions()[0]->GetExp(0)->GetNumPoints(0) ;
    int nTraceNumPoints = GetTraceNpoints();
    Array<OneD, NekDouble > uterm(nTraceNumPoints);
    Array<OneD, NekDouble > qtemp(nTraceNumPoints);
    
    m_fields[var]->ExtractTracePhys(physfield,qtemp);            

    for(i = 0; i < nbnd; ++i)
      {                 
      	// Evaluate boundary values g_D or g_N from input files
	SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(i);
	npoints = m_fields[0]->GetBndCondExpansions()[i]->GetNpoints();
	
	Array<OneD,NekDouble> BDphysics(npoints);
	Array<OneD,NekDouble> x0(npoints,0.0);
	Array<OneD,NekDouble> x1(npoints,0.0);
	Array<OneD,NekDouble> x2(npoints,0.0);  
        
	m_fields[0]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
	for(j = 0; j < npoints; j++)
	  {
	    BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
	  }
	
	// Weakly impose boundary conditions by modifying flux values
	for (e = 0; e < numBDEdge ; ++e)
	  {
	    id1 = m_fields[i]->GetBndCondExpansions()[0]->GetPhys_Offset(e);
	    id2 = m_fields[i]->GetTrace()->GetPhys_Offset(m_fields[i]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(e));
	    
	    // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ - g_D) (nx, ny)
	    if(m_fields[0]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
	      {
                  Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&qtemp[id2],1,&penaltyflux[id2],1);
                  
                  //  Vmath::Vsub(Nfps,&Fwd[id2],1,&BDphysics[id1],1,&uterm[id2],1);
                  //Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&uterm[id2],1,&uterm[id2],1);
                  //Vmath::Svtvp(Nfps,-1.0*C11,&uterm[id2],1,&qFwd[id2],1,&penaltyflux[id2],1);
	      }
	   
	    // For Neumann boundary condition: qflux = g_N
	    else if((m_fields[0]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
	      {
		Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
	      }
	  }
      }       
  }
  

   void FitzHughNagumo::MassMultiply(const Array<OneD, NekDouble> &inarray, 
				     Array<OneD, NekDouble> &outarray, 
				     const int direction, bool turnon )
   {
     int ncoeffs = inarray.num_elements();
     
     if(turnon)
       {
	 if(direction == -1)
	   {
	     m_fields[0]->MultiplyByInvMassMatrix(inarray,outarray,false);
	   }
	 else if(direction == 1)
	   {
	     
	     MultiRegions::GlobalMatrixKey key(StdRegions::eMass);
	     m_fields[0]->MultiRegions::ExpList::GeneralMatrixOp(key,inarray,outarray);
	   }
       }
     
     else
       {
	 Vmath::Vcopy(ncoeffs, inarray, 1, outarray, 1);
       }
   }  

  void FitzHughNagumo::Summary(std::ostream &out)
  {   
    cout << "=======================================================================" << endl;
    cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
    ADRBase::SessionSummary(out);
    if(m_explicitDiffusion)
      {
	out << "\tDiffusion Advancement   : Explicit" <<endl;
      }
    else
      {
	out << "\tDiffusion Advancement   : Implicit" <<endl;
      }
    if(m_explicitReaction)
      {
	out << "\tReaction Advancement    : Explicit" <<endl;
      }
    else
      {
	out << "\tReaction Advancement    : Implicit" <<endl;
      }

    out << "\tTime Integration Method : " << LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod] << endl;
    out << "\tEpsilon    : " << m_epsilon << endl;
    out << "\tBeta    : " << m_beta << endl;
    out << "\tinitialwavetype : " << m_initialwavetype << endl;
    out << "\tTimedelay : " << m_timedelay << endl;
    out << "\tDuration : " << m_duration << endl;

    ADRBase::TimeParamSummary(out);

    cout << "=======================================================================" << endl;
  }
    

} //end of namespace
