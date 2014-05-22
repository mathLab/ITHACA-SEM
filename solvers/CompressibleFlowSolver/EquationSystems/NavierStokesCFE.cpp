///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

namespace Nektar
{
    string NavierStokesCFE::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create, 
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleFlowSystem(pSession)
    {
    }

    void NavierStokesCFE::v_InitObject()
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
            m_ode.DefineOdeRhs     (&NavierStokesCFE::DoOdeRhs,        this);
            m_ode.DefineProjection (&NavierStokesCFE::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }

        /*
        m_checkpointFuncs["Sensor"] = boost::bind(&NavierStokesCFE::CPSensor, this, _1, _2);
        m_checkpointFuncs["SensorKappa"] = boost::bind(&NavierStokesCFE::CPSensorKappa, this, _1, _2);
        m_checkpointFuncs["SmoothVisc"] = boost::bind(&NavierStokesCFE::CPSmoothArtVisc, this, _1, _2);
        */
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    void NavierStokesCFE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type", ProblemTypeMap[m_problemType]);
    }

    void NavierStokesCFE::v_SetInitialConditions(
        NekDouble initialtime, 
        bool dumpInitialConditions,
        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);
        
        //insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);
        
        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements(); 
        
        if(Noise > 0.0)
        {
            for(int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot,Noise,noise,1,m_comm->GetColumnComm()->GetRank()+1);
                Vmath::Vadd(phystot,m_fields[i]->GetPhys(),1,noise,1,m_fields[i]->UpdatePhys(),1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
        }


        if (dumpInitialConditions)
        {
            // Dump initial conditions to file
            Checkpoint_Output(0);
        }
    }

    void NavierStokesCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        
        if(m_shockCaptureType == "Off")
        {
            Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > outarrayAdv(nvariables);
            Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
            
            Array<OneD, Array<OneD, NekDouble> > inarrayTemp(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            
            for (i = 0; i < nvariables; ++i)
            {
                outarrayAdv[i] = Array<OneD, NekDouble>(npoints, 0.0);
                outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            for (i = 0; i < nvariables-1; ++i)
            {
                inarrayTemp[i] = Array<OneD, NekDouble>(npoints, 0.0);
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            // Advection term in physical rhs form
            m_advection->Advect(nvariables, m_fields, advVel, inarray, outarrayAdv);

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(npoints, outarrayAdv[i], 1);
            }
            
            // Extract pressure and temperature
            Array<OneD, NekDouble > pressure   (npoints, 0.0);
            Array<OneD, NekDouble > temperature(npoints, 0.0);
            GetPressure(inarray, pressure);
            GetTemperature(inarray, pressure, temperature);
            
            // Extract velocities
            for (i = 1; i < nvariables-1; ++i)
            {
                Vmath::Vdiv(npoints,
                            inarray[i], 1,
                            inarray[0], 1,
                            inarrayTemp[i-1], 1);
            }
            
            // Copy velocities into new inarrayDiff
            for (i = 0; i < nvariables-2; ++i)
            {
                Vmath::Vcopy(npoints, inarrayTemp[i], 1, inarrayDiff[i], 1);
            }
            
            // Copy temperature into new inarrayDiffusion
            Vmath::Vcopy(npoints,
                         temperature, 1,
                         inarrayDiff[nvariables-2], 1);
            
            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff);
            
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
            Array<OneD, Array<OneD, NekDouble> > advVel;
            Array<OneD, Array<OneD, NekDouble> > outarrayAdv(nvariables);
            Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
            Array<OneD, Array<OneD, NekDouble> > outarrayForcing(nvariables);
            
            // contains u,v,w
            Array<OneD, Array<OneD, NekDouble> > inarrayTemp(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables+1);
            
            for (i = 0; i < nvariables; ++i)
            {
                outarrayAdv[i] = Array<OneD, NekDouble>(npoints, 0.0);
                outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
                outarrayForcing[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            for (i = 0; i < nvariables-1; ++i)
            {
                inarrayTemp[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            for (i = 0; i < nvariables+1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
            }
            
            // Advection term in physical rhs form
            m_advection->Advect(nvariables, m_fields, advVel, inarray, outarrayAdv);
            
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(npoints, outarrayAdv[i], 1);
            }
            
            // Extract pressure and temperature
            Array<OneD, NekDouble > pressure   (npoints, 0.0);
            Array<OneD, NekDouble > temperature(npoints, 0.0);
            Array<OneD, NekDouble > enthalpy   (npoints, 0.0);
            Array<OneD, NekDouble > energy     (npoints, 0.0);
            GetPressure(inarray, pressure);
            GetTemperature(inarray, pressure, temperature);
            GetEnthalpy(inarray, pressure, enthalpy);
            // Extract velocities
            for (i = 1; i < nvariables-2; ++i)
            {
                Vmath::Vdiv(npoints,
                            inarray[i], 1,
                            inarray[0], 1,
                            inarrayTemp[i-1], 1);
            }
            
            // Copy velocities into new inarrayDiff
            for (i = 0; i < nvariables-3; ++i)
            {
                Vmath::Vcopy(npoints, inarrayTemp[i], 1, inarrayDiff[i], 1);
            }
            
            // Copy temperature into new inarrayDiffusion
            Vmath::Vcopy(npoints,
                         temperature, 1,
                         inarrayDiff[nvariables-3], 1);
            
            // Copy density into new inarrayDiffusion
            Vmath::Vcopy(npoints,
                         inarray[0], 1,
                         inarrayDiff[nvariables-2], 1);
            
            // Copy artificial viscosity coefficient into new inarrayDiffusion
            Vmath::Vcopy(npoints,
                         inarray[nvariables-1], 1,
                         inarrayDiff[nvariables-1], 1);
            
            Vmath::Vdiv(npoints, inarray[nvariables-2], 1, inarray[0], 1, energy, 1);
            // Copy enthalpy into new inarrayDiffusion
            Vmath::Vcopy(npoints,
                         energy, 1,
                         inarrayDiff[nvariables], 1);
            
            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff);

            GetForcingTerm(inarray, outarrayForcing);
            
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayAdv[i], 1,
                            outarrayDiff[i], 1,
                            outarray[i], 1);
                
                // Add Forcing Term
                Vmath::Vadd(npoints,
                            outarray[i], 1,
                            outarrayForcing[i], 1,
                            outarray[i], 1);
            }
        }
        if (m_shockCaptureType == "NonSmooth")
        {
            ASSERTL0(false, "NS with non-smooth shock capturing not yet implemented");
        }
    }

    void NavierStokesCFE::DoOdeProjection(
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
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }
    
    void NavierStokesCFE::SetBoundaryConditions(
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
                ASSERTL0(false, "Wall is a wrong bc for the full "
                                "compressible Navier-Stokes equations");
            }
            
            // Wall Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eWallViscous)
            {
                WallViscousBC(n, cnt, inarray);
            }
            
            // artificial Condition
            if (m_fields[nvariables-1]->GetBndConditions()[n]->GetUserDefined() ==
                SpatialDomains::eArtificialViscosity)
            {
                ArtificialViscosityBC(n, cnt, inarray);
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
    
    void NavierStokesCFE::CPSensorKappa(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, NekDouble> &outarray)
    {
        const int npts = m_fields[0]->GetTotPoints();
        outarray = Array<OneD, NekDouble>(GetNcoeffs());
        Array<OneD, Array<OneD, NekDouble> > physfield(m_spacedim+2);
        
        for (int i = 0; i < m_spacedim+2; ++i)
        {
            physfield[i] = Array<OneD, NekDouble>(npts);
            m_fields[i]->BwdTrans(inarray[i], physfield[i]);
        }
        
        Array<OneD, NekDouble> sensor(npts,0.0);
        Array<OneD, NekDouble> SensorKappa(npts,0.0);
        GetSensor(physfield, sensor, SensorKappa);
        m_fields[0]->FwdTrans(SensorKappa, outarray);
    }
    
    void NavierStokesCFE::CPSensor(
                const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                      Array<OneD, NekDouble> &outarray)
    {
        const int npts = m_fields[0]->GetTotPoints();
        outarray = Array<OneD, NekDouble>(GetNcoeffs());
        Array<OneD, Array<OneD, NekDouble> > physfield(m_spacedim+2);
        
        for (int i = 0; i < m_spacedim+2; ++i)
        {
            physfield[i] = Array<OneD, NekDouble>(npts);
            m_fields[i]->BwdTrans(inarray[i], physfield[i]);
        }
        
        Array<OneD, NekDouble> sensor(npts,0.0);
        Array<OneD, NekDouble> SensorKappa(npts,0.0);
        GetSensor(physfield, sensor, SensorKappa);
        m_fields[0]->FwdTrans(sensor, outarray);
    }
    
    void NavierStokesCFE::CPSmoothArtVisc(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, NekDouble> &outarray)
    {
        const int npts = m_fields[0]->GetTotPoints();
        outarray = Array<OneD, NekDouble>(GetNcoeffs());
        Array<OneD, Array<OneD, NekDouble> > physfield(m_spacedim+3);
        
        for (int i = 0; i < m_spacedim+3; ++i)
        {
            physfield[i] = Array<OneD, NekDouble>(npts);
            m_fields[i]->BwdTrans(inarray[i], physfield[i]);
        }
        
        Array<OneD, NekDouble> eps_bar(npts, 0.0);
        GetSmoothArtificialViscosity(physfield, eps_bar);
        
        m_fields[0]->FwdTrans(eps_bar, outarray);
    }
}
