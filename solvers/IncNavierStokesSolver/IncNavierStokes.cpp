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

#include <IncNavierStokesSolver/IncNavierStokes.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    /**
     * Basic construnctor
     */
    IncNavierStokes::IncNavierStokes(void):
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

    IncNavierStokes::IncNavierStokes(string &fileNameString):
        ADRBase(fileNameString,false,true),
        m_infosteps(10)
    {
	switch(m_expdim)
	  {
	  case 2:
	    {
	      SpatialDomains::MeshGraph2DSharedPtr mesh2D;
              
	      if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)))
              {
		  ASSERTL0(false,"Dynamics cast failed");
              }
	      
              m_pressure =  MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*mesh2D,*m_boundaryConditions/??);
              
	      break;
	    }
	  case 3:
	    ASSERTL0(false,"3 D not set up");
	  default:
	    ASSERTL0(false,"Expansion dimension not recognised");
	    break;
	  }

        // Set up equation type enum using kEquationTypeStr
        
        const std::string typeStr = m_boundaryConditions->GetEquationTypeStr();

        for(int i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(nocase_cmp(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }

        
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eUnsteadyStokes: case eUnsteadyNavierStokes:
            
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
    }

    void IncNavierStokes::EvaluateAdvectionTerms(Array<OneD, Array<OneD, const NekDouble> > &Velocity,
                                                 Array<OneD, Array<OneD, const NekDouble> > &inarray, 
                                                 Array<OneD, Array< OneD, NekDouble> > &outarray, Array<OneD, NekDouble> &wk)
    {
        int j;
        int nvariables = inarray.num_elements();
        int nqtot      = m_fields[i]->GetTotPoints();
        int VelDim     = Velocity.num_elements();
        Array<OneD, OneD, NekDouble > Deriv;
        
        // Set up Derivative work space; 
        if(wk.num_elements())
        {
            ASSERTL0(wk.num_elements() > nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDoube> (nqtot*VelDim);;p
        }


        switch(m_advectionForm)
        {
        case eConvective: case eNonConservative:
            {
                AdvectionnNonConservativeForm(Velocity,inarray,outarray,Deriv);
            }
            break;
        default:
            ASSERTL0(false,"Advection form not known");
            break;
        }
        
    }
    

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
    
} //end of namespace

/**
* $Log: AdvectionDiffusionReaction.cpp,v $
* Revision 1.6  2009/01/06 21:10:34  sherwin
* Updates for virtual calls to IProductWRTBase and introduced reader to handle SOLVERINFO section to specify different solvers
*
* Revision 1.5  2008/11/19 10:53:51  pvos
* Made 2D CG version working
*
* Revision 1.4  2008/11/17 08:20:14  claes
* Temporary fix for CG schemes. 1D CG working (but not for userdefined BC). 1D DG not working
*
* Revision 1.3  2008/11/12 12:12:26  pvos
* Time Integration update
*
* Revision 1.2  2008/11/02 22:38:51  sherwin
* Updated parameter naming convention
*
* Revision 1.1  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
