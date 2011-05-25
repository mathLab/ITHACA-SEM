///////////////////////////////////////////////////////////////////////////////
//
// File EquationSystem.cpp
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
// Description: Main wrapper class for Advection Diffusion Reaction Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <string>
using std::string;

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    /**
     * @class EquationSystem
     *
     * This class is a base class for all solver implementations. It provides
     * the underlying generic functionality and interface for solving equations.
     *
     * To solve a steady-state equation, create a derived class from this class
     * and reimplement the virtual functions to provide custom implementation
     * for the problem.
     *
     * To solve unsteady problems, derive from the UnsteadySystem class instead
     * which provides general time integration.
     */

    /**
     * This constructor is protected as the objects of this class are never
     * instantiated directly.
     * @param   pSession        The session reader holding problem parameters.
     */
    EquationSystem::EquationSystem(SessionReaderSharedPtr& pSession)
        : ADRBase(pSession->GetFilename(), true),
          m_session(pSession)
    {
        // If a tangent vector policy is defined then the local tangent vectors
        // on each element need to be generated.
        if (pSession->DefinesGeometricInfo("TANGENTDIR"))
        {
            m_fields[0]->SetUpTangents();
        }

        // Zero all physical fields initially.
        ZeroPhysFields();
        filename = pSession->GetFilename();
    }

    /**
     *
     */
    EquationSystem::~EquationSystem()
    {

    }

    /**
     * This allows initialisation of the solver which cannot be completed
     * during object construction (such as setting of initial conditions).
     *
     * Public interface routine to virtual function implementation.
     */
    void EquationSystem::DoInitialise(void)
    {
        v_DoInitialise();
    }


    /**
     * Performs the actual solve.
     *
     * Public interface routine to virtual function implementation.
     */
    void EquationSystem::DoSolve(void)
    {
        v_DoSolve();
    }


    void EquationSystem::Output(void)
    {
        v_Output();
    }


    /**
     * Prints a summary of variables and problem parameters.
     *
     * Public interface routine to virtual function implementation.
     *
     * @param   out             The ostream object to write to.
     */
    void EquationSystem::PrintSummary(std::ostream &out)
    {
        out << "=======================================================================" << endl;
        out << "\tEquation Type   : " << m_session->GetSolverInfo("EQTYPE") << endl;
        ADRBase::SessionSummary(out);

        v_PrintSummary(out);

        out << "=======================================================================" << endl;
    }


    /**
     * Evaluates a physical function at each quadrature point in the domain.
     *
     * @param   pArray          The array into which to write the values.
     * @param   pEqn            The equation to evaluate.
     */
    void EquationSystem::EvaluateFunction(Array<OneD, NekDouble>& pArray,
                              SpatialDomains::ConstUserDefinedEqnShPtr pEqn)
    {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        if (pArray.num_elements() != nq)
        {
            pArray = Array<OneD, NekDouble>(nq);
        }
        for(int i = 0; i < nq; i++)
        {
            pArray[i] = pEqn->Evaluate(x0[i],x1[i],x2[i]);
        }

    }


    /**
     * If boundary conditions are time-dependent, they will be evaluated at the
     * time specified.
     *
     * @param   time            The time at which to evaluate the BCs
     */
    void EquationSystem::SetBoundaryConditions(NekDouble time)
    {
      int nvariables = m_fields.num_elements();
      for (int i = 0; i < nvariables; ++i)
      {
          m_fields[i]->EvaluateBoundaryConditions(time);
      }
    }


    /**
     * By default, nothing needs initialising at the EquationSystem level.
     */
    void EquationSystem::v_DoInitialise()
    {

    }
    
    void EquationSystem::InitialiseForce()
    {
       int nq = m_fields[0]->GetNpoints();
       ADRBase::SetInitialForce(0.0);
       //copy force in m_fields if it exists otherwise set m_fields to zero 
       if(m_bforce)
       {	     
    	  for(int i=0; i<m_FDim; ++i)
    	  {
    	       Vmath::Vcopy(nq,(m_forces[i]->GetPhys()),1,(m_fields[i]
        	->UpdatePhys()),1);
          }
       }
       else
       {     
          for(int i=0; i<m_FDim; ++i)
    	  {
    	       Vmath::Zero(nq,(m_fields[i]
        	->UpdatePhys()),1);
          }    
       }		
    }
    
    void EquationSystem::InitialiseBaseFlow(Array<OneD, Array<OneD, NekDouble> > &base)
    {
    	base = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
    	int nq = m_fields[0]->GetNpoints();
        std::string velStr[3] = {"Vx","Vy","Vz"};        
        if(m_boundaryConditions->SolverInfoExists("BASEFLOWFILE"))
        {
        	std::string baseyn= m_boundaryConditions->GetSolverInfo("BASEFLOWFILE");
        	//capitalise the string baseyn
        	std::transform(baseyn.begin(), baseyn.end(), baseyn.begin()
        		, ::toupper);
        	if(baseyn=="YES")
        	{
        		SetUpBaseFields(m_graph);
        		string basename = filename;
        		basename= (filename).substr(0
        		, basename.find_last_of("."));
        		ImportFldBase(basename + "-Base.fld",m_graph);
        		cout<<"Base flow from file:  "<<basename<<"-Base.fld"<<endl;
        		for(int i=0; i<m_spacedim; ++i)
        		{
        			base[i] = Array<OneD, NekDouble> (nq,0.0);
        			Vmath::Vcopy(nq,m_base[i]->GetPhys(),1,base[i],1);	
        		}
        	}
        	else
                {
                	for(int i = 0; i < m_spacedim; ++i)
                	{	
        		   base[i] = Array<OneD, NekDouble> (nq,0.0);
        		   SpatialDomains::ConstUserDefinedEqnShPtr ifunc
        		   = m_boundaryConditions->GetUserDefinedEqn(velStr[i]);
        		   EvaluateFunction(base[i],ifunc);
        		}
        	}
	  }        
	  else
	  {
        	for(int i = 0; i < m_spacedim; ++i)
        	{	
        		base[i] = Array<OneD, NekDouble> (nq,0.0);
        		SpatialDomains::ConstUserDefinedEqnShPtr ifunc
        		= m_boundaryConditions->GetUserDefinedEqn(velStr[i]);
               		EvaluateFunction(base[i],ifunc);
              	}
           }
    }    	    
    
    void
    EquationSystem::SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr
    	    &mesh)
    {
    	   int i;
       	   //NUM VARIABLES CAN BE DIFFERENT FROM THE DIMENSION OF THE BASE FLOW
    	   m_base =Array<OneD, MultiRegions::ExpListSharedPtr> (m_spacedim);
    	   if (m_projectionType = ADRBase::eGalerkin)
	   {
	        switch (m_expdim)
	        {
	           case 1:
	           {
                       SpatialDomains::MeshGraph1DSharedPtr mesh1D;
                       if( !(mesh1D = boost::dynamic_pointer_cast<
                        SpatialDomains::MeshGraph1D>(mesh)) )
                       {
                           ASSERTL0(false,"Dynamics cast failed");
                       }
                       for(i = 0 ; i < m_base.num_elements(); i++)
                       {
                           m_base[i] = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(*mesh1D,
                                                    *m_boundaryConditions,i);
                       }
	           }
                   break;
                   case 2:
	           {
                       SpatialDomains::MeshGraph2DSharedPtr mesh2D;
                       if(!(mesh2D = boost::dynamic_pointer_cast<
                         SpatialDomains::MeshGraph2D>(mesh)))
                       {
                           ASSERTL0(false,"Dynamics cast failed");
                       }
                       i = 0;
                       MultiRegions::ContField2DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*mesh2D,
                                                    *m_boundaryConditions,i);
                       m_base[0]=firstbase;
                       for(i = 1 ; i < m_base.num_elements(); i++)
                       {
                            m_base[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*firstbase,*mesh2D,
                                                    *m_boundaryConditions,i);
                       }
	           }
                   break;
	           case 3:
	           {
	               SpatialDomains::MeshGraph3DSharedPtr mesh3D;
                       if(!(mesh3D = boost::dynamic_pointer_cast<
                          SpatialDomains::MeshGraph3D>(mesh)))
                       {
                            ASSERTL0(false,"Dynamics cast failed");
                       }
                       MultiRegions::ContField3DSharedPtr firstbase =
                            MemoryManager<MultiRegions::ContField3D>
                                ::AllocateSharedPtr(*mesh3D,
                                                    *m_boundaryConditions,i);
                       m_base[0] = firstbase;
                       for(i = 1 ; i < m_base.num_elements(); i++)
                       {
                           m_base[i] = MemoryManager<MultiRegions::ContField3D>
                           ::AllocateSharedPtr(*firstbase,
                                        *mesh3D,*m_boundaryConditions,i);
                       }
	           }
                   break;
                   default:
                     ASSERTL0(false,"Expansion dimension not recognised");
                   break;
	       }
        }
        else
	{
               switch(m_expdim)
               {
                   case 1:
                   {
                       SpatialDomains::MeshGraph1DSharedPtr mesh1D;
                       if(!(mesh1D = boost::dynamic_pointer_cast<SpatialDomains
                         ::MeshGraph1D>(mesh)))
                       {
                           ASSERTL0(false,"Dynamics cast failed");
                       }
                       for(i = 0 ; i < m_base.num_elements(); i++)
                       {
                           m_base[i] = MemoryManager<MultiRegions
                            ::DisContField1D>::AllocateSharedPtr(*mesh1D,
                                                             *m_boundaryConditions,i);
                       }
                       break;
                   }
                   case 2:
                   {
                         SpatialDomains::MeshGraph2DSharedPtr mesh2D;
                         if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains
                         ::MeshGraph2D>(mesh)))
                         {
                               ASSERTL0(false,"Dynamics cast failed");
                         }
                         for(i = 0 ; i < m_base.num_elements(); i++)
                         {
                                 m_base[i] = MemoryManager<MultiRegions
                                 ::DisContField2D>::AllocateSharedPtr(*mesh2D,
                                                             *m_boundaryConditions,i);
                         }
                         break;
                    }
                    case 3:
                         ASSERTL0(false,"3 D not set up");
                    default:
                          ASSERTL0(false,"Expansion dimension not recognised");
                    break;
                    }
	  }
       }    	    
 	
    //Import base flow from file and load in m_base    	
    void EquationSystem::ImportFldBase(std::string pInfile, 
    	    SpatialDomains::MeshGraphSharedPtr pGraph)
    {
    	    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
    	    std::vector<std::vector<NekDouble>   > FieldData;
    	    pGraph->Import(pInfile, FieldDef,FieldData);
       	    int nvar= m_spacedim;
      	    //copy data to m_velocity
    	    for(int j=0; j <nvar; ++j)
    	    {
    	       for(int i=0; i<FieldDef.size(); ++i)
    	       {
    	    	  bool flag = FieldDef[i]->m_fields[j]
    	             ==m_boundaryConditions->GetVariable(j);
    	             ASSERTL1(flag, (std::string("Order of ") +pInfile
    	             	     + std::string("  data and that defined in "
    	             	     	     "m_boundaryconditions differs")).c_str());   
    	          m_base[j]->ExtractDataToCoeffs(FieldDef[i]
    	            	    , FieldData[i], FieldDef[i]->m_fields[j]);
    	       }
    	    	m_base[j]->BwdTrans(m_base[j]->GetCoeffs(), m_base[j]
    	    		->UpdatePhys());    
   	    }	    
     }	      	    

    /**
     *
     */
    void EquationSystem::v_DoSolve()
    {

    }


    /**
     * By default, there are no further parameters to display.
     */
    void EquationSystem::v_PrintSummary(std::ostream &out)
    {

    }

    /** 
     *
     */
    void EquationSystem::v_Output(void)
    {
        ADRBase::Output();
    }

}
