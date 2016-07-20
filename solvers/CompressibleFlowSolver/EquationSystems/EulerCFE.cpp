///////////////////////////////////////////////////////////////////////////////
//
// File EulerCFE.cpp
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
// Description: Euler equations in consƒervative variables without artificial
// diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <CompressibleFlowSolver/EquationSystems/EulerCFE.h>

using namespace std;

namespace Nektar
{
    string EulerCFE::className = 
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "EulerCFE", EulerCFE::create, 
        "Euler equations in conservative variables without "
        "artificial diffusion.");
    
    EulerCFE::EulerCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : CompressibleFlowSystem(pSession)
    {
    }
    
    void EulerCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&EulerCFE::DoOdeRhs,        this);
            m_ode.DefineProjection (&EulerCFE::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }
    }
    
    /**
     * @brief Destructor for EulerCFE class.
     */
    EulerCFE::~EulerCFE()
    {
    }

    /**
     * @brief Compute the right-hand side.
     */
    void EulerCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advection->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time);

        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }
        
        // Add sponge layer if defined in the session file
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
        }
    }
    
    /**
     * @brief Compute the projection and call the method for imposing the 
     * boundary conditions in case of discontinuous projection.
     */
    void EulerCFE::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        
        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();
                
                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                }
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for Euler equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }

}
