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
// Description: Euler equations in conservative variables with artificial
// diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/EulerADCFE.h>
#include <boost/algorithm/string.hpp>

namespace Nektar
{
    string EulerADCFE::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "EulerADCFE", EulerADCFE::create,
        "Euler equations in conservative variables with "
        "artificial diffusion.");
    
    EulerADCFE::EulerADCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : CompressibleFlowSystem(pSession)
    {
    }
    
    void EulerADCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();
        
        if (m_shockCaptureType == "Smooth")
        {
            ASSERTL0(m_fields.num_elements() == m_spacedim + 3,
                     "Not enough variables for smooth shock capturing; "
                     "make sure you have added eps to variable list.");
            m_smoothDiffusion = true;
        }
        
        m_diffusion->SetArtificialDiffusionVector(
            &EulerADCFE::GetArtificialDynamicViscosity, this);
        
        if(m_session->DefinesSolverInfo("PROBLEMTYPE"))
        {
            std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
            int i;
            for(i = 0; i < (int) SIZE_ProblemType; ++i)
            {
                if(boost::iequals(ProblemTypeMap[i], ProblemTypeStr))
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
            m_ode.DefineOdeRhs    (&EulerADCFE::
            DoOdeRhs, this);
            m_ode.DefineProjection(&EulerADCFE::
            DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }
    }
    
    EulerADCFE::~EulerADCFE()
    {
        
    }
    
    void EulerADCFE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(
            s, "Problem Type", ProblemTypeMap[m_problemType]);
    }
    
    void EulerADCFE::v_SetInitialConditions(
        NekDouble initialtime, 
        bool      dumpInitialConditions,
        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);
        
        if(dumpInitialConditions)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }
    
    void EulerADCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        
        Array<OneD, Array<OneD, NekDouble> > advVel;
        Array<OneD, Array<OneD, NekDouble> > outarrayAdv(nvariables);
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);

        for (i = 0; i < nvariables; ++i)
        {
            outarrayAdv[i] = Array<OneD, NekDouble>(npoints, 0.0);
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }
        
        m_advection->Advect(nvariables, m_fields, advVel, inarray,
                                        outarrayAdv, m_time);
        
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarrayAdv[i], 1);
        }
        
        m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff);
        
        if (m_shockCaptureType == "NonSmooth")
        {
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayAdv[i], 1,
                            outarrayDiff[i], 1,
                            outarray[i], 1);
            }
        }
        if(m_shockCaptureType == "Smooth")
        {
            const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();
            
            NekDouble pOrder = Vmath::Vmax(ExpOrder.num_elements(), ExpOrder, 1);
            
            Array <OneD, NekDouble > a_vel  (npoints, 0.0);
            Array <OneD, NekDouble > u_abs  (npoints, 0.0);
            Array <OneD, NekDouble > pres   (npoints, 0.0);
            Array <OneD, NekDouble > wave_sp(npoints, 0.0);
            
            GetPressure(inarray, pres);
            GetSoundSpeed(inarray, pres, a_vel);
            GetAbsoluteVelocity(inarray, u_abs);
            
            Vmath::Vadd(npoints, a_vel, 1, u_abs, 1, wave_sp, 1);
            
            NekDouble max_wave_sp = Vmath::Vmax(npoints, wave_sp, 1);
            
            Vmath::Smul(npoints,
                        m_C2,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);
            
            Vmath::Smul(npoints,
                        max_wave_sp,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);
            
            Vmath::Smul(npoints,
                        pOrder,
                        outarrayDiff[nvariables-1], 1,
                        outarrayDiff[nvariables-1], 1);
            
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayAdv[i], 1,
                            outarrayDiff[i], 1,
                            outarray[i], 1);
            }
            
            Array<OneD, Array<OneD, NekDouble> > outarrayForcing(nvariables);
            
            for (i = 0; i < nvariables; ++i)
            {
                outarrayForcing[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            GetForcingTerm(inarray, outarrayForcing);
            
            for (i = 0; i < nvariables; ++i)
            {
                // Add Forcing Term
                Vmath::Vadd(npoints,
                            outarray[i], 1,
                            outarrayForcing[i], 1,
                            outarray[i], 1);
            }
        }
    }
    
    void EulerADCFE::DoOdeProjection(
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
    
    void EulerADCFE::SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble                             time)
    {    
        std::string varName;
        int nvariables = m_fields.num_elements();
        int cnt        = 0;
        
        // loop over Boundary Regions
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Wall Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eWall)
            {
                WallBC(n, cnt, inarray);
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
                SymmetryBC(n, cnt, inarray);
            }
            
            // Riemann invariant characteristic Boundary Condition (CBC)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eRiemannInvariant)
            {
                RiemannInvariantBC(n, cnt, inarray);
            }
            
            // Extrapolation of the data at the boundaries
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == 
                SpatialDomains::eExtrapOrder0)
            {
                ExtrapOrder0BC(n, cnt, inarray);
            }
            
            // Time Dependent Boundary Condition (specified in meshfile)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() 
                == SpatialDomains::eTimeDependent)
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
            }
            
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }
}
