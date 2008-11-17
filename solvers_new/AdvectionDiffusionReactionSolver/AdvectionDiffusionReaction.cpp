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
    
    /**
     * Constructor. Creates ... of #DisContField2D fields
     *
     * \param 
     * \param
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(string &fileNameString):
        ADRBase(fileNameString,true),
        m_infosteps(10)
    {
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
        for(int i = 0; i < m_spacedim; ++i)
        {
            m_velocity[i] = Array<OneD, NekDouble> (GetPointsTot());
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

    void AdvectionDiffusionReaction::EvaluateAdvectionVelocity()
    {
        int nq = m_fields[0]->GetPointsTot();
        
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
    

    void AdvectionDiffusionReaction::ODEforcing(const Array<OneD, const Array<OneD, NekDouble> >&inarray,  
						Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time) 
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
                m_fields[i]->MultiplyByElmtInvMass(outarray[i],outarray[i]);
		Vmath::Neg(ncoeffs,outarray[i],1);
            }
	  break;
        case eGalerkin:
	  {
                Array<OneD, NekDouble> physfield(GetPointsTot());
		
                for(i = 0; i < nvariables; ++i)
		  {
                    // Calculate -(\phi, V\cdot Grad(u))
                    m_fields[i]->BwdTrans(inarray[i],physfield);
		    
		    WeakAdvectionNonConservativeForm(m_velocity,
						     physfield, outarray[i]);
		    
                    Vmath::Neg(ncoeffs,outarray[i],1);
                    
                    // Multiply by inverse of mass matrix to get forcing term
		    // m_fields[i]->MultiplyByInvMassMatrix(outarray[i],  
                    //                                     outarray[i],
                    //                                     false, true);
		   		    
                }
            }
            break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

    void AdvectionDiffusionReaction::ExplicitlyIntegrateAdvection(int nsteps)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Get Integration scheme details
        LibUtilities::TimeIntegrationSchemeKey       IntKey(LibUtilities::eForwardEuler);
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
	    phys[i] = Array<OneD, NekDouble>(m_fields[0]->GetPointsTot());
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
		// this is just a forward Euler to illustate that CG works
		 
		// get -D u^n
		ODEforcing(in,out,m_time); // note that MultiplyByInvMassMatrix is not performed inside ODEforcing
	  
		// compute M u^n
		for (i = 0; i < nvariables; ++i)
		  {
		    m_fields[0]->BwdTrans(in[i],phys[i]);
		    m_fields[0]->IProductWRTBase(phys[i],tmp[i]);
		    
		    // f = M u^n - Dt D u^n
		    Vmath::Svtvp(ncoeffs,m_timestep,out[i],1,tmp[i],1,tmp[i],1);
		    
		    // u^{n+1} = M^{-1} f
		    m_fields[i]->MultiplyByInvMassMatrix(tmp[i],out[i],false,false);
		    
		    // fill results
		    Vmath::Vcopy(ncoeffs,out[i],1,in[i],1);
		    Vmath::Vcopy(ncoeffs,out[i],1,fields[i],1);
		  }
		//---------------------------------------------------------
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
    
    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
      {	
	
	// Time Dependent Boundary Condition
	if (m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() == "TimeDependent")
	  {
	    for (int i = 0; i < nvariables; ++i)
	      {
		m_fields[i]->EvaluateBoundaryConditions(time);
	      }
	  }
      }
  }
  
  // Evaulate flux = m_fields*ivel for i th component of Vu 
  void AdvectionDiffusionReaction::GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, 
						 Array<OneD, Array<OneD, NekDouble> > &flux)
  {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetPointsTot(),physfield[i],1,
                        m_velocity[j],1,flux[j],1);
        }
    }

    void AdvectionDiffusionReaction::NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
						   Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTracePointsTot();
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

        int nTraceNumPoints = GetTracePointsTot();
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
      cout << "\tEquation Type   : Advection Equation" << endl;
      ADRBase::Summary(out);
      cout << "=======================================================================" << endl;

    }
} //end of namespace

/**
* $Log: AdvectionDiffusionReaction.cpp,v $
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
