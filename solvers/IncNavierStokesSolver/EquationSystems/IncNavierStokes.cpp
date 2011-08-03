///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.cpp
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
#include <Auxiliary/EquationSystem.h>
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
    IncNavierStokes::IncNavierStokes(LibUtilities::CommSharedPtr& pComm,
            LibUtilities::SessionReaderSharedPtr& pSession):
        EquationSystem(pComm, pSession),
        m_infosteps(10)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        EquationSystem::v_InitObject();

        int i,j,expdim;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};
        
        // Set up Velocity field to point to the first m_expdim of m_fields; 
        m_velocity = Array<OneD,int>(m_spacedim);
        
        for(i = 0; i < m_spacedim; ++i)
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
            m_session->MatchSolverInfo("EQTYPE",kEquationTypeStr[i],match,false);
            if(match)
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }
        ASSERTL0(i != eEquationTypeSize,"EQTYPE not found in SOLVERINFO section");
        
        // This probably should to into specific implementations 
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eSteadyStokes: 
        case eSteadyOseen: 
        case eSteadyLinearisedNS: 
            break;
        case eUnsteadyNavierStokes:
        case eUnsteadyStokes:
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
            
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
                        if(m_fields[0]->GetBndConditions()[n]->GetUserDefined().GetEquation() != "I")
                        {  	 	 
                            ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                        }
                    }
                }
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
        
        m_session->LoadParameter("Kinvis", m_kinvis);
        
        if (m_equationType == eUnsteadyNavierStokes)
        {
            std::string vConvectiveType = "Convective";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                vConvectiveType = m_session->GetTag("AdvectiveType");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_comm, m_session, m_graph, m_boundaryConditions);
        }
    }

    IncNavierStokes::~IncNavierStokes(void)
    {

    }
    
    void IncNavierStokes::AdvanceInTime(int nsteps)
    {
        int i,n;
        int phystot = m_fields[0]->GetTotPoints();
        static int nchk = 0;
        
        int n_fields = m_fields.num_elements();
	
        //if(m_HomogeneousType != eNotHomogeneous) //Homogeneous case semi-phys integration
        //{
        //    for(i = 0; i < n_fields; ++i)
        //    {
        //        m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
        //        m_fields[i]->SetFourierSpace(MultiRegions::eCoef);
        //        m_fields[i]->SetPhysState(false);
        //    }
        //}
	
        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >   fields(m_nConvectiveFields);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields[i]  = m_fields[i]->UpdatePhys();
        }
        
        // Initialise NS solver which is set up to use a GLM method
        // with calls to EvaluateAdvection_SetPressureBCs and
        // SolveUnsteadyStokesSystem
        LibUtilities::TimeIntegrationSolutionSharedPtr 
            IntegrationSoln = m_integrationScheme[m_intSteps-1]->InitializeScheme(m_timestep, fields, m_time, m_integrationOps);
        
        std::string   outname = m_session->GetFilename();
        outname = (m_session->GetFilename()).substr(0, outname.find_last_of(".")) + ".his";
        std::ofstream hisFile (outname.c_str());

        //Time advance
        for(n = 0; n < nsteps; ++n)
        {
            // Advance velocity fields
            fields = m_integrationScheme[min(n,m_intSteps-1)]->TimeIntegrate(m_timestep, IntegrationSoln, m_integrationOps);
            
            m_time += m_timestep;
       		
            if(!((n+1)%m_infosteps))
            {
                cout << "Step: " << n+1 << "  Time: " << m_time << endl;
                WriteHistoryData(hisFile);
            }
            
            // dump data in m_fields->m_coeffs to file. 
            if(n&&(!((n+1)%m_checksteps)))
            {
                //if(m_HomogeneousType != eNotHomogeneous)
                //{
                //    for(i = 0; i < n_fields; ++i)
                //    {
                //        m_fields[i]->SetFourierSpace(MultiRegions::ePhys);
                //        m_fields[i]->SetPhysState(false);
                //    }
                //    
                //    Checkpoint_Output(nchk++);
                //    
                //    for(i = 0; i < n_fields; ++i)
                //    {
                //        m_fields[i]->SetFourierSpace(MultiRegions::eCoef);
                //        m_fields[i]->SetPhysState(false);
                //    }
                //}
                //else 
                //{
                    for(i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_fields[i]->SetPhys(fields[i]);
                        m_fields[i]->SetPhysState(true);
                    }
                    Checkpoint_Output(nchk++);
                    
                //}
            }
        }
        
        //updating physical space
        //if(m_HomogeneousType != eNotHomogeneous)
        //{
        //    for(i = 0; i < n_fields; ++i)
        //    {
        //        m_fields[i]->SetFourierSpace(MultiRegions::ePhys);
        //        m_fields[i]->SetPhysState(false);
        //    }			
        //}
        //else 
        //{
            for(i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[i]->SetPhys(fields[i]);
                m_fields[i]->SetPhysState(true);
            }
        //}
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
            ASSERTL0(wk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }
        
        m_advObject->DoAdvection(m_fields,m_nConvectiveFields, 
                                 m_velocity,inarray,outarray,Deriv);
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
