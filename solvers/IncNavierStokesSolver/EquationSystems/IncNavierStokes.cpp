///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokesSolver.cpp
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
// Description: Incompressible Navier Stokes class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    
    IncNavierStokes::IncNavierStokes(SessionReaderSharedPtr& pSession, string &globoptfile):
        EquationSystem(pSession),
        m_infosteps(10)
    {
        int i,j;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};
                
        // Set up Velocity field to point to the first m_expdim of m_fields; 
        m_velocity = Array<OneD,int>(m_spacedim);
        
        for(i = 0; i < m_expdim; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(NoCaseStringCompare(velids[i],var) == 0)
                {
                    m_velocity[i] = j;
                    break;
                }
                
                if(j == numfields)
                {
                    std::string error = "Failed to find field: " + var; 
                    ASSERTL0(false,error.c_str());
                }
            }
        }
        
         // Set up equation type enum using kEquationTypeStr
         for(i = 0; i < (int) eEquationTypeSize; ++i)
         {
             bool match;
             pSession->MatchSolverInfo("EQTYPE",kEquationTypeStr[i],match,false);
             if(match)
             {
                 m_equationType = (EquationType)i; 
                 break;
             }
         }
         ASSERTL0(i != eEquationTypeSize,"EQTYPE not found in SOLVERINFO section");

         m_advectionForm = eNoAdvectionForm;
         
         // This probably should to into specific implementations 
         // Equation specific Setups 
         switch(m_equationType)
         {
         case eSteadyStokes: 
         case eSteadyOseen: 
             break;
         case eUnsteadyNavierStokes:
             {
                 // Set up advection form
                 for(i = 0; i < (int) eAdvectionFormSize; ++i)
                 {
                     bool match;
                     pSession->MatchSolverInfo("ADVECTIONFORM",kAdvectionFormStr[i],match,false);
                     
                     if(match)
                     {
                         m_advectionForm = (AdvectionForm)i; 
                         break;
                     }
                 }
                 ASSERTL0(i != eAdvectionFormSize,"ADVECTIONFORM not found in SOLVERINFO section");
             }
         case eUnsteadyStokes:
             
             if(m_boundaryConditions->CheckForParameter("IO_InfoSteps") == true)
             {
                 m_infosteps =  m_boundaryConditions->GetParameter("IO_InfoSteps");
             }

             // check to see if any user defined boundary condition is
             // indeed implemented
             
             for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
             {	
                 // Time Dependent Boundary Condition (if no user
                 // defined then this is empty)
                 if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "")
                 {
                     if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "TimeDependent")
                     {
                         ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                     }
                 }
             }
             break;
         case eNoEquationType:
         default:
             ASSERTL0(false,"Unknown or undefined equation type");
         }
         
        if(m_boundaryConditions->CheckForParameter("Kinvis") == true)
        {
            m_kinvis = m_boundaryConditions->GetParameter("Kinvis");
        }
        else
        {
            ASSERTL0(false,"Kinvis is not specified");
        }
    }

    IncNavierStokes::~IncNavierStokes(void)
    {

    }

    // Evaluation -N(V) for all fields except pressure using m_velocity
    void IncNavierStokes::EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                 Array<OneD, NekDouble> &wk)
    {
        int i,j;
        int nvariables = inarray.num_elements();
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, NekDouble > Deriv;
        
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]]; 
        }

        // Set up Derivative work space; 
        if(wk.num_elements())
        {
            ASSERTL0(wk.num_elements() > nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }


        switch(m_advectionForm)
        {
        case eConvective: case eNonConservative:
            {
                int i;
                for(i = 0; i < m_nConvectiveFields; ++i)
                {
                    AdvectionNonConservativeForm(velocity,inarray[i],outarray[i],Deriv);
                    Vmath::Neg(nqtot,outarray[i],1);
                }
            }
            break; 
        default:
            ASSERTL0(false,"Advection form not known");
            break;
        }        
    }
	
	//time dependent boundary conditions updating
	
	void IncNavierStokes::SetBoundaryConditions(NekDouble time)
	{
		int  nvariables = m_fields.num_elements();
		
		for (int i = 0; i < nvariables; ++i)
		{
			for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
			{	
			   if(m_fields[i]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
			   {
				  m_fields[i]->EvaluateBoundaryConditions(time);
			   }
			}
		}
	}
    

    // case insensitive string comparison from web
} //end of namespace

/**
* $Log: IncNavierStokes.cpp,v $
* Revision 1.3  2010/01/28 15:16:03  abolis
* Time-Dependent boundary conditions
*
* Revision 1.2  2009/09/06 22:31:15  sherwin
* First working version of Navier-Stokes solver and input files
*
*
**/
