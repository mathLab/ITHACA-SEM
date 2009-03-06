///////////////////////////////////////////////////////////////////////////////
//
// File FHN.cpp
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
// Description: FHN class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <FitzHugh-Nagumo/FHN.h>
#include <cstdio>
#include <cstdlib>
namespace Nektar
{
    /**
     * Basic construnctor
     */
    FHN::FHN(void):
        ADRBase(),
        m_infosteps(100)
    {     
    }
    
    int nocase_cmp(const string & s1, const string& s2);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    FHN::FHN(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10)
    {

        int i;

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        // for(i = 0; i < (int) eEquationTypeSize; ++i)
        for(i = 0; i < 1; ++i)
        {
            if(nocase_cmp(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }        
        
        // Equation Setups 
  
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (GetNpoints());
            }
            
            EvaluateAdvectionVelocity();
            
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
    }

    void FHN::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetNpoints();
        
        std::string velStr[3] = {"Vx","Vy","Vz"};

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_velocity.num_elements(); i++)
	{
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);
            
            for(int j = 0; j < nq; j++)
	    {
                m_velocity[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
	    }
	}
    }
    
    void FHN::ODEeReaction(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
			   Array<OneD, Array<OneD, NekDouble> >&outarray, 
			   const NekDouble time)
	
	{
		int i,k;
                int nvariables = inarray.num_elements();
                int ncoeffs    = inarray[0].num_elements();
		const NekDouble coeff = 3.14159265*3.14159265;
					
		for (i = 0; i < nvariables; ++i)
		{
		     Vmath::Smul(ncoeffs, coeff, inarray[i], 1, outarray[i], 1);
		}
	}
	   

    void FHN::ODEhelmSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			   Array<OneD, Array<OneD, NekDouble> >&outarray,
			   NekDouble time, NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
							
		// We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
		// inarray = input: \hat{rhs} -> output: \hat{Y}    
		// outarray = output: nabla^2 \hat{Y}       
		// where \hat = modal coeffs
		
		MultiRegions::GlobalLinSysKey key(StdRegions::eMass);
		
		for (int i = 0; i < nvariables; ++i)
		{
		        // Multiply by inverse of mass matrix
                         m_fields[i]->MultiplyByInvMassMatrix(inarray[i],outarray[i],false);
				
			// Multiply rhs[i] with -1.0/gamma/timestep
                         Vmath::Smul(ncoeffs, -1.0/lambda, outarray[i], 1, outarray[i], 1);
			
			// Update coeffs to m_fields
			 m_fields[i]->UpdateCoeffs() = outarray[i];
			
			// Backward Transformation to nodal coefficients
			  m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());

                          NekDouble kappa = 1.0/lambda;
			  			
	    	        // Solve a system of equations with Helmholtz solver
                          m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs(),kappa);
			
			// The solution is Y[i]
			  outarray[i] = m_fields[i]->GetCoeffs();	  
						  
			// Multiply back by mass matrix
                           m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,outarray[i],outarray[i]);
		}
	}
	

    void FHN::GeneralTimeIntegration(int nsteps, 
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
                // Project the solution u* onto the boundary conditions to
                // obtain the actual solution
                SetBoundaryConditions(m_time);
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
	      Checkpoint_Output(nchk++);
            }
        }
        
        for(i = 0; i < nvariables; ++i)
        {
	  (m_fields[i]->UpdateCoeffs()) = fields[i];
        }
    }
	
    
  //----------------------------------------------------
  void FHN::SetBoundaryConditions(NekDouble time)
  {
    int nvariables = m_fields.num_elements();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->EvaluateBoundaryConditions(time);
    }
    
//     // loop over Boundary Regions
//     for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
//       {	
	
// 	// Time Dependent Boundary Condition
// 	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
// 	  {
// 	    for (int i = 0; i < nvariables; ++i)
// 	      {
// 		m_fields[i]->EvaluateBoundaryConditions(time);
// 	      }
// 	  }
//       }
  }
  
 

    void FHN::Summary(std::ostream &out)
    {   
      cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
      ADRBase::SessionSummary(out);
      ADRBase::TimeParamSummary(out);
      cout << "=======================================================================" << endl;

    }
    
    /*
    // case insensitive string comparison from web
    int nocase_cmp(const string & s1, const string& s2) 
    {
        string::const_iterator it1=s1.begin();
        string::const_iterator it2=s2.begin();
        
        //stop when either string's end has been reached
        while ( (it1!=s1.end()) && (it2!=s2.end()) ) 
        { 
            if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
            {
                // return -1 to indicate smaller than, 1 otherwise
                return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1; 
            }
            //proceed to the next character in each string
            ++it1;
            ++it2;
        }
        size_t size1=s1.size(), size2=s2.size();// cache lengths
        
        //return -1,0 or 1 according to strings' lengths
        if (size1==size2) 
        {
            return 0;
        }
        return (size1 < size2) ? -1 : 1;
    }
    */
    
} //end of namespace

/**
* $Log: FHN.cpp,v $
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Sehun' FHN solver
*
**/
