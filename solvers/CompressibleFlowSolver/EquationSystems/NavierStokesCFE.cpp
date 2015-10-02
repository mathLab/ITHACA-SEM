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
            for(int i = 0; i < (int) SIZE_ProblemType; ++i)
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
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    void NavierStokesCFE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        CompressibleFlowSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Problem Type",
                                    ProblemTypeMap[m_problemType]);
    }

    void NavierStokesCFE::v_SetInitialConditions(
        NekDouble initialtime,
        bool dumpInitialConditions,
        const int domain)
    {
        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements();

        if (Noise > 0.0)
        {
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1,
                                      m_comm->GetColumnComm()->GetRank()+1);
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1,
                            noise, 1, m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
            }
        }

        CompressibleFlowSystem::v_SetInitialConditions();

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
        m_advection->Advect(nvariables, m_fields, advVel, inarray,
                            outarrayAdv, time);

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
            Vmath::Vsub(npoints,
                        outarrayDiff[i], 1,
                        outarrayAdv[i], 1,
                        outarray[i], 1);
        }

        // Add sponge layer if defined in the session file
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
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
        int cnt        = 0;

        // loop over Boundary Regions
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            std::string type = m_fields[0]->GetBndConditions()[n]->GetUserDefined();
            SetCommonBC(type,n,time,cnt,inarray);
        }
    }
}
