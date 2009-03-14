///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrection.cpp
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
// Description: Velocity Correction Scheme for the Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/VelocityCorrectionScheme.h>

namespace Nektar
{
    int nocase_cmp(const string & s1, const string& s2);


    /**
     * Basic construnctor
     */
    VelocityCorrectionScheme::VelocityCorrectionScheme(void):
        IncNavierStokes()
    {     
    }
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    VelocityCorrectionScheme::VelocityCorrectionScheme(string &fileNameString):
        IncNavierStokes(fileNameString)
    {

        m_advectionOps.DefineOdeRhs(&IncNavierStokes::EvaluateAdvection, this);
        

        // >>>> Need a method to define integration scheme order here. 
        // Get Integration scheme details
        LibUtilities::TimeIntegrationSchemeKey       IntKey(m_AdvectionintegrationType);
        m_AdvectionScheme = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

    }

    void VelocityCorrectionScheme::AdvanceInTime(int nsteps)
    {
        int nInitSteps;
        LibUtilities::TimeIntegrationSolutionSharedPtr 
            AdvectionSoln = m_IntScheme->InitializeScheme(m_timestep,
                                                          m_time,
                                                          nInitSteps,
                                                          fields,
                                                          m_advectionOps);


        for(n = nInitSteps; n < nsteps; ++n)
        {
            // Advance Velocity
            fields = m_advectionScheme->TimeIntegrate(m_timestep,
                                       AdvectionSoln,m_advectionOps);
            
            // Set Pressure Boundary Conditions
            ->Find element along edges.
                  -> Get offset and calculate vorticity and N.L. Terms
                  -> extrapolate. 
            
            // Pressure Forcing 

            
            // Solver Pressure Poisson Equation 
            

            // Viscous Term forcing
            
            
            // Solve Helmholtz system
            fields1 = m_diffusionScheme->TimeIntegrate(m_timestep,
                                                       Diffusion,
                                                       m_diffusionOps);
        }

    }
    
    
} //end of namespace

/**
* $Log: $
*
**/
