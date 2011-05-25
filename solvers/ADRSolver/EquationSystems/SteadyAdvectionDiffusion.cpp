///////////////////////////////////////////////////////////////////////////////
//
// File SteadAdvectionDiffusion.cpp
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
// Description: Steady Advection Diffusion solve routines 
//
///////////////////////////////////////////////////////////////////////////////
#include <ADRSolver/EquationSystems/SteadyAdvectionDiffusion.h>

namespace Nektar
{
    string SteadyAdvectionDiffusion::className = EquationSystemFactory::RegisterCreatorFunction("SteadyAdvectionDiffusion", SteadyAdvectionDiffusion::create);

    /**
     * @class SteadyAdvectionDiffusion
     * This is a solver class for solving the  problems.
     * - SteadyAdvectionDiffusion:
     *   \f$ c \cdot \nabla u -\nabla \cdot (\nabla u)  = f(x)\f$
     */

    SteadyAdvectionDiffusion::SteadyAdvectionDiffusion(SessionReaderSharedPtr& pSession)
        : EquationSystem(pSession),
          m_lambda(0.0)
    {
        // Define forcing function. 
        EquationSystem::InitialiseForce();    
        // Define Velocity fields     
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim); 
        EquationSystem::InitialiseBaseFlow(m_velocity);
    }
       
    SteadyAdvectionDiffusion::~SteadyAdvectionDiffusion()
    {

    }

    void SteadyAdvectionDiffusion::v_PrintSummary(std::ostream &out)
    {
        out << "\tLambda          : " << m_lambda << endl;
    }


    void SteadyAdvectionDiffusion::v_DoSolve()
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->LinearAdvectionDiffusionReactionSolve(m_velocity,
                                                               m_fields[i]->GetPhys(),
                                                               m_fields[i]->UpdateCoeffs(),
                                                               m_lambda);
            m_fields[i]->SetPhysState(false);
        }
    }

}
