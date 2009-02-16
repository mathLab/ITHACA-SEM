///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.cpp
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

#include <AdvectionDiffusionReactionSolver/AdvectionDiffusionReaction.h>
#include <cstdio>
#include <cstdlib>
namespace Nektar
{
    /**
     * Basic construnctor
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(void):
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
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10)
    {

        int i;

        // Set up equation type enum using kEquationTypeStr
        std::string typeStr = m_boundaryConditions->GetSolverInfo("EQTYPE");

        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            if(nocase_cmp(kEquationTypeStr[i],typeStr) == 0 )
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }

        ASSERTL0(i != (int) eEquationTypeSize, "Invalid expansion type.");
        
        
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eHelmholtz:
            break;
        case eAdvection:
            m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
            for(int i = 0; i < m_spacedim; ++i)
            {
                m_velocity[i] = Array<OneD, NekDouble> (GetNpoints());
            }
            
            EvaluateAdvectionVelocity();
            goto UnsteadySetup;
            break;

        UnsteadySetup:
            
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
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
    }

    void AdvectionDiffusionReaction::EvaluateAdvectionVelocity()
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
    
    void AdvectionDiffusionReaction::ODErhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
                                                  Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                            const NekDouble time) 
    {
        int i;
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();

	SetBoundaryConditions(time);
	
        switch(m_projectionType)
        {
        case eDiscontinuousGalerkin:
	  
	  WeakDGAdvection(inarray, outarray);
	  for(i = 0; i < nvariables; ++i)
            {
		Vmath::Neg(ncoeffs,outarray[i],1);
            }
	  break;
        case eGalerkin:
	  {
                Array<OneD, NekDouble> physfield(GetNpoints());
		
                for(i = 0; i < nvariables; ++i)
		  {
                    // Calculate -(\phi, V\cdot Grad(u))
                    m_fields[i]->BwdTrans_IterPerExp(inarray[i],physfield);
		    
		    WeakAdvectionNonConservativeForm(m_velocity,
						     physfield, outarray[i]);
		    
                    Vmath::Neg(ncoeffs,outarray[i],1);		   		    
                }
            }
            break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }
    
    void AdvectionDiffusionReaction::ODElhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
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
    
    void AdvectionDiffusionReaction::ODElhsSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
                                                       Array<OneD,       Array<OneD, NekDouble> >&outarray, 
                                                 const NekDouble time)   
    {
	SetBoundaryConditions(time);
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

    void AdvectionDiffusionReaction::SolveHelmholtz(NekDouble lambda)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->HelmSolve(*(m_fields[i]),lambda);
        }
    }


    void AdvectionDiffusionReaction::ExplicitlyIntegrateAdvection(int nsteps)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        LibUtilities::TimeIntegrationSchemeOperators ode;
        ode.DefineOdeLhs       (&AdvectionDiffusionReaction::ODElhs,      this);
        ode.DefineOdeLhsSolve  (&AdvectionDiffusionReaction::ODElhsSolve, this);
        ode.DefineOdeRhs       (&AdvectionDiffusionReaction::ODErhs,      this);

        // Get Integration scheme details
        LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eAdamsBashforthOrder2);
        LibUtilities::TimeIntegrationSchemeSharedPtr IntScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        
        for(i = 0; i < nvariables; ++i)
        {
            fields[i]  = m_fields[i]->UpdateCoeffs();
        }
                
        int nInitSteps;
        LibUtilities::TimeIntegrationSolutionSharedPtr u = IntScheme->InitializeScheme(m_timestep,m_time,nInitSteps,fields,ode);

        for(n = nInitSteps; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------

            fields = IntScheme->TimeIntegrate(m_timestep,u,ode);
            m_time += m_timestep;

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
  void AdvectionDiffusionReaction::SetBoundaryConditions(NekDouble time)
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
  
  // Evaulate flux = m_fields*ivel for i th component of Vu 
  void AdvectionDiffusionReaction::GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &flux)
  {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                        m_velocity[j],1,flux[j],1);
        }
    }

    void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						   Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

	// Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
	  {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
	  }

        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            //evaulate upwinded m_fields[i]
	    m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
	    // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }

  void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &numfluxX, 
						 Array<OneD, Array<OneD, NekDouble> > &numfluxY)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
	Array<OneD, NekDouble > tmp(nTraceNumPoints,0.0);

	Array<OneD, Array<OneD, NekDouble > > traceVelocity(2);

	traceVelocity[0] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);
	traceVelocity[1] = Array<OneD,NekDouble>(nTraceNumPoints,0.0);

	// Get Edge Velocity - Could be stored if time independent
	m_fields[0]->ExtractTracePhys(m_velocity[0], traceVelocity[0]);
	m_fields[0]->ExtractTracePhys(m_velocity[1], traceVelocity[1]);

	m_fields[0]->GetFwdBwdTracePhys(physfield[0],Fwd,Bwd);
	
	m_fields[0]->GetTrace()->Upwind(traceVelocity,Fwd,Bwd,tmp);
	
	Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[0],1,numfluxX[0],1);
	Vmath::Vmul(nTraceNumPoints,tmp,1,traceVelocity[1],1,numfluxY[0],1);
  
    }

    void AdvectionDiffusionReaction::Summary(std::ostream &out)
    {   
      cout << "=======================================================================" << endl;
      cout << "\tEquation Type   : "<< kEquationTypeStr[m_equationType] << endl;
      ADRBase::SessionSummary(out);
      switch(m_equationType)
      {
      case eSteadyDiffusion: case eSteadyDiffusionReaction:
      case eHelmholtz: case eLaplace: case ePoisson:
          out << "\tLambda          : " << m_boundaryConditions->GetParameter("Lambda") << endl;
          
          break;
      case eAdvection:
          ADRBase::TimeParamSummary(out);
          break;
      }
      cout << "=======================================================================" << endl;

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
* Revision 1.7  2009/02/10 16:39:35  sherwin
* Added new format of SolverInfo reader to identify EQTYPE
*
* Revision 1.6  2009/02/08 09:13:08  sherwin
* Updates to go with Multiple matrix/variable solve
*
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
