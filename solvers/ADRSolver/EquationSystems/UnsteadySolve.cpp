///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadySolve.cpp
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
// Description: Generic information for Unseady 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <ADRSolver/EquationSystems/UnsteadySolve.h>

namespace Nektar
{

    UnsteadySolve::UnsteadySolve(SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession)
    {
        
        // load generic input parameters
        if(pSession->DefinesParameter("IO_InfoSteps"))
        {
            m_infosteps =  pSession->GetParameter("IO_InfoSteps");
        }

        if (pSession->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"))
        {
            int i;
            for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
            {
                bool match;
                pSession->MatchSolverInfo("TIMEINTEGRATIONMETHOD", LibUtilities::TimeIntegrationMethodMap[i], match, false);
                if (match)
                {
                    m_timeIntMethod = (LibUtilities::TimeIntegrationMethod) i;
                    break;
                }
            }
            ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
                                            "Invalid time integration type.");
        }

        SetInitialConditions();

    }

    UnsteadySolve::~UnsteadySolve()
    {

    }


    void UnsteadySolve::v_DoSolve() // Time integration routines in physical space
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->SetPhysState(false);
            fields[i]  = m_fields[i]->UpdatePhys();
        }

        // Declare an array of TimeIntegrationSchemes For multi-stage
        // methods, this array will have just one entry containing the
        // actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this
        //    is an initialization scheme)
        //  - the second scheme will used for the second timestep
        //    (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps
        //    (this will be the actual scheme)

        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(m_timeIntMethod)
        {
        case LibUtilities::eIMEXdirk_2_3_2:
        case LibUtilities::eIMEXdirk_3_4_3:
        case LibUtilities::eDIRKOrder2:
        case LibUtilities::eDIRKOrder3:
        case LibUtilities::eBackwardEuler:
        case LibUtilities::eForwardEuler:
        case LibUtilities::eClassicalRungeKutta4:
        case LibUtilities::eIMEXOrder1:
            {
                numMultiSteps = 1;
                
                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
                
                LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];
                
                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;
                
                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
                
                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);
                
                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
        case LibUtilities::eIMEXOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);
				
                // Used for all other time steps 
                LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod); 
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
            
        case LibUtilities::eIMEXOrder3:
            {
                numMultiSteps = 3;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);
		
                // Used in the first time step to initalize the scheme		
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eIMEXOrder1);
                LibUtilities::TimeIntegrationSchemeKey IntKey1(LibUtilities::eIMEXOrder2);
                
                // Used for all other time steps 
                LibUtilities::TimeIntegrationSchemeKey IntKey2(m_timeIntMethod); 
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                IntScheme[2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];
                
                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[2]->InitializeScheme(m_timestep,fields,m_time,m_ode);
                break;
            }
            default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }
        
        std::string outname = m_session->GetFilename() + ".his";
        std::ofstream hisFile (outname.c_str());

        for(n = 0; n < m_steps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,m_ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
            }
            
            m_time += m_timestep;

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%m_infosteps))
            {
                cout << "\rSteps: " << n+1 << "\t Time: " << m_time << "\t " << endl;
            }

            if(n&&(!((n+1)%m_checksteps)))
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->FwdTrans(fields[i],m_fields[i]->UpdateCoeffs());
                }
                Checkpoint_Output(nchk++);
                WriteHistoryData(hisFile);
            }
        }

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->UpdatePhys() = fields[i];
        }
    }
    

    void UnsteadySolve::v_PrintSummary(std::ostream &out)
    {
        ADRBase::Summary(out);
    }


}
