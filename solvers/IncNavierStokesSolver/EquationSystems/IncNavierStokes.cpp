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
    IncNavierStokes::IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession):
        EquationSystem(pSession),
        m_infosteps(10),
        m_steadyStateSteps(0)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        EquationSystem::v_InitObject();
        
        int i,j;
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
            m_session->LoadParameter("IO_EnergySteps", m_energysteps, 0);
            m_session->LoadParameter("SteadyStateSteps", m_steadyStateSteps, 0);
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 1e-6);
            
            // check to see if any user defined boundary condition is
            // indeed implemented
            
            for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {	
                // Time Dependent Boundary Condition (if no user
                // defined then this is empty)
                if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eNoUserDefined)
                {
                    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eTimeDependent)
                    {                     	     
                        if(m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eI)
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
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
		
		if (m_equationType == eUnsteadyLinearisedNS)
        {
            std::string vConvectiveType = "Linearised";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                vConvectiveType = m_session->GetTag("Linearised");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
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
		
		bool integrate_in_wave_space = false;
		
		if(m_session->DefinesSolverInfo("INTEGRATIONSPACE"))
		{
			if(m_HomogeneousType == eNotHomogeneous)
			{
				ASSERTL0(false,"INTEGRATIONSPACE type is meant to be for homogeneous cases");
			}
			
			std::string IntegrationSpaceStr = m_session->GetSolverInfo("INTEGRATIONSPACE");
			
			if((IntegrationSpaceStr == "WaveSpace") || (IntegrationSpaceStr == "WAVESPACE"))
			{
				integrate_in_wave_space = true;
			}
			else if((IntegrationSpaceStr == "RealSpace") || (IntegrationSpaceStr == "REALSPACE"))
			{
				integrate_in_wave_space = false;
			}
			else 
			{
				ASSERTL0(false,"INTEGRATIONSPACE type not allowed, try WaveSpace or RealSpace");
			}
		}

		//SingleMode integration in wave space
		if(m_session->DefinesSolverInfo("SingleMode")==true && m_session->GetSolverInfo("SingleMode")=="ModifiedBasis")
		{
			integrate_in_wave_space = true;
		}
		
        int n_fields = m_fields.num_elements();
		

		if(integrate_in_wave_space)
		{
			for(i = 0; i < n_fields; ++i)
			{
				m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
				m_fields[i]->SetWaveSpace(true);
				m_fields[i]->SetPhysState(false);
			}
		}
	
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
		
        std::string   mdlname = m_session->GetSessionName() + ".mdl";
        std::ofstream mdlFile;

        if (m_energysteps)
        {
            mdlFile.open(mdlname.c_str());
        }

        //Time advance
        for(n = 0; n < nsteps; ++n)
        {
            // Advance velocity fields
            fields = m_integrationScheme[min(n,m_intSteps-1)]->TimeIntegrate(m_timestep, IntegrationSoln, m_integrationOps);
            
            m_time += m_timestep;
       		
            // Write out current time step
            if(m_session->GetComm()->GetRank() == 0
               && !((n+1)%m_infosteps))
            {
                cout << "Step: " << n+1 << "\t Time: " << m_time << endl;
            }

            // Write out energy data to file
            if(m_energysteps && !((n+1)%m_energysteps))
            {
                if(m_HomogeneousType != eNotHomogeneous)
                {
                    if(m_HomogeneousType == eHomogeneous1D)
                    {
                        Array<OneD, NekDouble> Energy(m_npointsZ,0.0);
                        Array<OneD, NekDouble> Energy_tmp(m_npointsZ,0.0);
			
                        for(i = 0; i < m_nConvectiveFields; ++i)
                        {
                            Energy_tmp = m_fields[i]->HomogeneousEnergy();
                            Vmath::Vadd(m_npointsZ,Energy_tmp,1,Energy,1,Energy,1);
                        }
                        mdlFile << m_time << "\t";
                        for(i = 0; i < m_npointsZ; i++)
                        {
                            mdlFile << Energy[i] << "\t";
                        }
                        mdlFile << endl;
                    }
                    else
                    {
                        ASSERTL0(false,"3D Homogeneous 2D energy dumping not implemented yet");
                    }
                }
                else
                {
                    NekDouble energy = 0.0;
                    for(i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_fields[i]->SetPhys(fields[i]);
                        m_fields[i]->SetPhysState(true);
                        NekDouble norm = L2Error(i, true);
                        energy += norm*norm;
                    }
                    mdlFile << m_time << "   " << 0.5*energy << endl;
                }
                
            }
            
            // dump data in m_fields->m_coeffs to file. 
            if(m_checksteps && n&&(!((n+1)%m_checksteps)))
            {
				if(integrate_in_wave_space)
				{
					for(i = 0; i< n_fields; i++)
					{
						m_fields[i]->SetWaveSpace(false);
						m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
						m_fields[i]->SetPhysState(true);
					}
					nchk++;
					Checkpoint_Output(nchk);
					for(i = 0; i< n_fields; i++)
					{
						m_fields[i]->SetWaveSpace(true);
						m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
						m_fields[i]->SetPhysState(false);
					}
				}
				else
				{
					for(i = 0; i < m_nConvectiveFields; ++i)
					{
						m_fields[i]->SetPhys(fields[i]);
						m_fields[i]->SetPhysState(true);
					}
					nchk++;
					Checkpoint_Output(nchk);
				}
            }

            
            if(m_steadyStateSteps && n && (!((n+1)%m_steadyStateSteps)))
            {
                if(CalcSteadyState() == true)
                {
                    cout << "Reached Steady State to tolerance " << m_steadyStateTol << endl; 
                    break;
                }
            }
        }
		
		if(integrate_in_wave_space)
		{
			for(i = 0; i< n_fields; i++)
			{
				m_fields[i]->SetWaveSpace(false);
				m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
				m_fields[i]->SetPhysState(true);
			}
		}
		else 
		{
			for(i = 0; i < m_nConvectiveFields; ++i)
			{
				m_fields[i]->SetPhys(fields[i]);
				m_fields[i]->SetPhysState(true);
			}
		}
		
        
        if (m_energysteps)
        {
            mdlFile.close();
        }
    }
    

    // Evaluation -N(V) for all fields except pressure using m_velocity
    void IncNavierStokes::EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                 Array<OneD, NekDouble> &wk)
    {
        int i;
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
                                 m_velocity,inarray,outarray,m_time,Deriv);
    }
	
    //time dependent boundary conditions updating
    
    void IncNavierStokes::SetBoundaryConditions(NekDouble time)
    {
        int  nvariables = m_fields.num_elements();
	
        for (int i = 0; i < nvariables; ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {	
                if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }
        }
    }
    
    // Decide if at a steady state if the discrerte L2 sum of the
    // coefficients is the same as the previous step to within the
    // tolerance m_steadyStateTol;
    bool IncNavierStokes::CalcSteadyState(void)
    {
        static NekDouble previousL2 = 0.0;
        bool returnval = false;

        NekDouble L2 = 0.0;
        
        // calculate L2 discrete summation 
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            L2 += Vmath::Dot(ncoeffs,m_fields[i]->GetCoeffs(),1,m_fields[i]->GetCoeffs(),1);
        }
        
        if(fabs(L2-previousL2) < ncoeffs*m_steadyStateTol)
        {
            returnval = true;
        }

        previousL2 = L2;

        return returnval;
    }

} //end of namespace

