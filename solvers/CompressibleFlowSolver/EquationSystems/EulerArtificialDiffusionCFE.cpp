///////////////////////////////////////////////////////////////////////////////
//
// File EulerArtificialDiffusionCFE.cpp
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
// Description: Euler equations in conservative variables with artificial diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/EulerArtificialDiffusionCFE.h>

namespace Nektar
{
    string EulerArtificialDiffusionCFE::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "EulerArtificialDiffusionCFE", EulerArtificialDiffusionCFE::create, 
            "Euler equations in conservative variables with "
            "artificial diffusion.");

    EulerArtificialDiffusionCFE::EulerArtificialDiffusionCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleFlowSystem(pSession)
    {
    }

    void EulerArtificialDiffusionCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        if(m_session->DefinesSolverInfo("PROBLEMTYPE"))
        {

            std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
            int i;
            for(i = 0; i < (int) SIZE_ProblemType; ++i)
            {
                if(NoCaseStringCompare(ProblemTypeMap[i],ProblemTypeStr) == 0)
                {
                    m_problemType = (ProblemType)i;
                    break;
                }
            }
        }
        else
        {
            m_problemType = (ProblemType)0;
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs    (&EulerArtificialDiffusionCFE::
                                    DoOdeRhs, this);
            m_ode.DefineProjection(&EulerArtificialDiffusionCFE::
                                    DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }

    }

    EulerArtificialDiffusionCFE::~EulerArtificialDiffusionCFE()
    {

    }

    void EulerArtificialDiffusionCFE::v_PrintSummary(std::ostream &out)
    {
        CompressibleFlowSystem::v_PrintSummary(out);
        out << "\tProblem Type    : " << ProblemTypeMap[m_problemType] << endl;
    }

    void EulerArtificialDiffusionCFE::v_SetInitialConditions(
        NekDouble initialtime, 
        bool      dumpInitialConditions)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);

        if(dumpInitialConditions)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }

    void EulerArtificialDiffusionCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
                
        Array<OneD, Array<OneD, NekDouble> > advVel;
                
        m_advection->Advect(nvariables, m_fields, advVel, inarray, outarray);
                
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }
    }

    void EulerArtificialDiffusionCFE::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
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
    
    void EulerArtificialDiffusionCFE::SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble                             time)
    {    
        int nvariables = m_fields.num_elements();
        int nq         = inarray[0].num_elements();
        int cnt        = 0;
    
        // loop over Boundary Regions
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Wall Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eWall)
            {
                WallBoundary(n, cnt, inarray);
            }
            
            // Wall Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eWallViscous)
            {
                ASSERTL0(false, "WallViscous is a wrong bc for the "
                         "Euler equations");
            }
    
            // Symmetric Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eSymmetry)
            {
                SymmetryBoundary(n, cnt, inarray);
            }
            
            // Inflow characteristic Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eInflowCFS)
            {
                InflowCFSBoundary(n, cnt, inarray);
            }
            
            // Outflow characteristic Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eOutflowCFS)
            {
                OutflowCFSBoundary(n, cnt, inarray);
            }
            
            // Extrapolation of the data at the boundaries
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eExtrapOrder0)
            {
                ExtrapOrder0Boundary(n, cnt, inarray);
            }
    
            // Time Dependent Boundary Condition (specified in meshfile)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() 
                == SpatialDomains::eTimeDependent)
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }
    
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }
}
