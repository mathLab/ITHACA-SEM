///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSystem.cpp
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
// Description: Generic timestepping for Pulse Wave Solver
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <MultiRegions/ContField.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
#include <LibUtilities/BasicUtils/Timer.h>
using namespace std;

namespace Nektar
{

/**
 *  @class PulseWaveSystem
 *
 *  Initialises the arterial subdomains in m_vessels and sets up
 *  all domain-linking conditions (bifurcations, junctions,
 *  merging flows). Detects the network structure and assigns
 *  boundary conditons. Also provides the underlying timestepping
 *  framework for pulse wave solvers including the general
 *  timestepping routines.
 */

/**
 *  Processes SolverInfo parameters from the session file and sets
 *  up timestepping-specific code.
 *
 *  @param   m_Session        Session object to read parameters from.
 */
PulseWaveSystem::PulseWaveSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}

/**
 *  Destructor
 */
PulseWaveSystem::~PulseWaveSystem()
{
}

/**
 *  Initialisation routine for multidomain solver. Sets up the
 *  expansions for every arterial segment (m_vessels) and for one
 *  complete field m_outfield which is needed to write the
 *  postprocessing output. Also determines which upwind strategy
 *  is used (currently only upwinding scheme available) and reads
 *  blodd flow specific parameters from the inputfile
 *
 */
void PulseWaveSystem::v_InitObject()
{
    // Initialise base class
    UnsteadySystem::v_InitObject();

    // Read the geometry and the expansion information
    m_nDomains = m_graph->GetDomain().size();

    // Determine projectiontype
    ASSERTL0(m_session->MatchSolverInfo("Projection", "DisContinuous"),
             "Pulse solver only set up for Discontinuous projections");
    m_projectionType = MultiRegions::eDiscontinuous;
    ASSERTL0(m_graph->GetMeshDimension() == 1,
             "Pulse wave solver only set up for expansion dimension equal to 1");

    int i;
    m_nVariables = m_session->GetVariables().size();

    m_fields = Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVariables);
    m_vessels =
        Array<OneD, MultiRegions::ExpListSharedPtr>(m_nVariables * m_nDomains);

    /* If the PressureArea property is not specified, default to the Beta law;
     * it's the most well known and this way previous examples that did not
     * specify the tube law will still work.
     */
    if (m_session->DefinesSolverInfo("PressureArea"))
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance(
                m_session->GetSolverInfo("PressureArea"), m_vessels, m_session);
    }
    else
    {
        m_pressureArea = GetPressureAreaFactory().CreateInstance("Beta",
                                                          m_vessels, m_session);
    }

    m_pressure = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);

    m_PWV  = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_W1   = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_W2   = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);

    const std::vector<SpatialDomains::CompositeMap> domain = m_graph->GetDomain();

    SpatialDomains::BoundaryConditions Allbcs(m_session, m_graph);

    // Set up domains and put geometry to be only one space dimension.
    int cnt = 0;
    bool SetToOneSpaceDimension = true;

    if (m_session->DefinesCmdLineArgument("SetToOneSpaceDimension"))
    {
        std::string cmdline =
            m_session->GetCmdLineArgument<std::string>("SetToOneSpaceDimension");
        if (boost::to_upper_copy(cmdline) == "FALSE")
        {
            SetToOneSpaceDimension = false;
        }
    }

    for (i = 0; i < m_nDomains; ++i)
    {
        for (int j = 0; j < m_nVariables; ++j)
        {
            m_vessels[cnt++] =
                MemoryManager<MultiRegions::DisContField>::AllocateSharedPtr(
                    m_session, m_graph, domain[i], Allbcs,
                        m_session->GetVariable(j), SetToOneSpaceDimension);
        }
    }

    // Reset coeff and phys space to be continuous over all domains
    int totcoeffs = 0;
    int totphys   = 0;
    m_fieldPhysOffset = Array<OneD, int>(m_nDomains + 1, 0);
    for (i = 0; i < m_nDomains; ++i)
    {
        totcoeffs += m_vessels[i * m_nVariables]->GetNcoeffs();

        m_fieldPhysOffset[i] = totphys;
        totphys += m_vessels[i * m_nVariables]->GetTotPoints();
    }
    m_fieldPhysOffset[m_nDomains] = totphys;

    for (int n = 0; n < m_nVariables; ++n)
    {
        Array<OneD, NekDouble> coeffs(totcoeffs, 0.0);
        Array<OneD, NekDouble> phys(totphys, 0.0);
        Array<OneD, NekDouble> tmpcoeffs, tmpphys;

        m_vessels[n]->SetCoeffsArray(coeffs);
        m_vessels[n]->SetPhysArray(phys);

        int cnt  = m_vessels[n]->GetNcoeffs();
        int cnt1 = m_vessels[n]->GetTotPoints();

        for (i = 1; i < m_nDomains; ++i)
        {
            m_vessels[i * m_nVariables + n]->SetCoeffsArray(tmpcoeffs = coeffs + cnt);
            m_vessels[i * m_nVariables + n]->SetPhysArray(tmpphys = phys + cnt1);
            cnt  += m_vessels[i * m_nVariables + n]->GetNcoeffs();
            cnt1 += m_vessels[i * m_nVariables + n]->GetTotPoints();
        }
    }

    for (int i = 0; i < 2; ++i)
    {
        m_fields[i] = m_vessels[i];
    }

    // Zero all physical fields initially.
    ZeroPhysFields();

    // If Discontinuous Galerkin determine upwinding method to use
    for (int i = 0; i < (int)SIZE_UpwindTypePulse; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("UPWINDTYPEPULSE", UpwindTypeMapPulse[i],
                                                                  match, false);
        if (match)
        {
            m_upwindTypePulse = (UpwindTypePulse)i;
            break;
        }
    }

    // Load blood density and external pressure
    m_session->LoadParameter("rho", m_rho, 0.5);
    m_session->LoadParameter("pext", m_pext, 0.0);

    int nq = 0;
    /*
     *  Gets the Material Properties of each arterial segment
     *  specified in the inputfile from section MaterialProperties
     *  Also gets the Area at static equilibrium A_0 specified in the
     *  inputfile.
     *
     * Having found these points also extract the values at the
     * trace points and the normal direction consistent with the
     * left adjacent definition of Fwd and Bwd
     */
    m_beta             = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_beta_trace       = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_gamma            = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_alpha            = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_alpha_trace      = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_A_0              = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_A_0_trace        = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);
    m_trace_fwd_normal = Array<OneD, Array<OneD, NekDouble> >(m_nDomains);

    for (int omega = 0; omega < m_nDomains; omega++)
    {
        nq = m_vessels[2 * omega]->GetNpoints();
        m_fields[0] = m_vessels[2 * omega];

        m_beta[omega] = Array<OneD, NekDouble>(nq);
        GetFunction("MaterialProperties")
            ->Evaluate("beta", m_beta[omega], m_time, omega);

        // If the input file doesn't specify viscoelasticity, set it to zero.
        m_gamma[omega] = Array<OneD, NekDouble>(nq);

        if (m_session->DefinesFunction("Viscoelasticity"))
        {
            GetFunction("Viscoelasticity")
                ->Evaluate("gamma", m_gamma[omega], m_time, omega);
        }
        else
        {
            for (int j = 0; j < nq; ++j)
            {
                m_gamma[omega][j] = 0;
            }
        }

        /* If the input file doesn't specify strain-stiffening, set it to
         * 0.5 for elastic behaviour.
         */
        m_alpha[omega] = Array<OneD, NekDouble>(nq);

        if (m_session->DefinesFunction("StrainStiffening"))
        {
            GetFunction("StrainStiffening")
                ->Evaluate("alpha", m_alpha[omega], m_time, omega);
        }
        else
        {
            for (int j = 0; j < nq; ++j)
            {
                m_alpha[omega][j] = 0.5;
            }
        }

        m_A_0[omega] = Array<OneD, NekDouble>(nq);
        GetFunction("A_0")->Evaluate("A_0", m_A_0[omega], m_time, omega);

        int nqTrace = GetTraceTotPoints();

        m_beta_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_beta[omega], m_beta_trace[omega]);

        m_A_0_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_A_0[omega], m_A_0_trace[omega]);

        m_alpha_trace[omega] = Array<OneD, NekDouble>(nqTrace);
        m_fields[0]->ExtractTracePhys(m_alpha[omega], m_alpha_trace[omega]);

        if (SetToOneSpaceDimension)
        {
            m_trace_fwd_normal[omega] = Array<OneD, NekDouble>(nqTrace, 0.0);

            MultiRegions::ExpListSharedPtr trace = m_fields[0]->GetTrace();
            int nelmt_trace = trace->GetExpSize();

            Array<OneD, Array<OneD, NekDouble> > normals(nelmt_trace);

            for (int i = 0; i < nelmt_trace; ++i)
            {
                normals[i] = m_trace_fwd_normal[omega] + i;
            }

            // need to set to 1 for consistency since boundary
            // conditions may not have coordim = 1
            trace->GetExp(0)->GetGeom()->SetCoordim(1);

            trace->GetNormals(normals);
        }
    }

    SetUpDomainInterfaces();
}

void PulseWaveSystem::SetUpDomainInterfaces(void)
{
    map<int, std::vector<InterfacePointShPtr> > VidToDomain;

    /* Loop over domain and find out if we have any undefined
     * boundary conditions representing interfaces. If so make a
     * map based around vid and storing the domains that are
     * part of interfaces.
     */
    for (int omega = 0; omega < m_nDomains; ++omega)
    {
        int vesselID = omega * m_nVariables;

        for (int i = 0; i < 2; ++i)
        {
            if (m_vessels[vesselID]
                    ->GetBndConditions()[i]
                    ->GetBoundaryConditionType() == SpatialDomains::eNotDefined)
            {
                // Get Vid of interface
                int vid = m_vessels[vesselID]
                              ->UpdateBndCondExpansion(i)
                              ->GetExp(0)
                              ->GetGeom()
                              ->GetVid(0);
                cout << "Shared vertex ID: " << vid << endl;
                MultiRegions::ExpListSharedPtr trace = m_vessels[vesselID]->GetTrace();
                InterfacePointShPtr Ipt;

                bool finish = false;
                // find which elmt, the local vertex and the data offset of point
                for (int n = 0; n < m_vessels[vesselID]->GetExpSize(); ++n)
                {
                    for (int p = 0; p < 2; ++p)
                    {
                        if (m_vessels[vesselID]
                                ->GetTraceMap()
                                ->GetElmtToTrace()[n][p]
                                ->as<LocalRegions::Expansion>()
                                ->GetGeom()
                                ->GetVid(0) == vid)
                        {
                            int eid = m_vessels[vesselID]
                                          ->GetTraceMap()
                                          ->GetElmtToTrace()[n][p]
                                          ->GetElmtId();

                            int tid =
                                m_vessels[vesselID]->GetTrace()->GetCoeff_Offset(
                                    eid);
                            Ipt =
                                MemoryManager<InterfacePoint>::AllocateSharedPtr(
                                    vid, omega, n, p, tid, i);

                            cout << "Global VID of interface point: " << vid << endl;
                            cout << "Domain interface point belongs to: " << omega << endl;
                            cout << "Element ID of vertex: " << n << endl;
                            cout << "Vertex ID within local element: " << p << endl;
                            cout << "Element ID within the trace: " << tid << endl;
                            cout << "Position of boundary condition in region: " << i << endl;

                            finish = true;
                            break;
                        }
                    }
                    if (finish == true)
                    {
                        break;
                    }
                }

                VidToDomain[vid].push_back(Ipt);

                // finally reset boundary condition to Dirichlet
                m_vessels[vesselID]->GetBndConditions()[i]->SetBoundaryConditionType(
                    SpatialDomains::eDirichlet);

                m_vessels[vesselID + 1]
                    ->GetBndConditions()[i]
                    ->SetBoundaryConditionType(SpatialDomains::eDirichlet);
            }
        }
    }

    // loop over map and set up Interface information;
    for (auto &iter : VidToDomain)
    {
        if (iter.second.size() == 2) // Vessel jump interface
        {
            m_vesselJcts.push_back(iter.second);
        }
        else if (iter.second.size() == 3) // Bifurcation or Merging junction.
        {
            int nbeg = 0;
            int nend = 0;

            // Determine if bifurcation or merging junction
            // through number of element vertices that meet at
            // junction. Only one vertex using a m_elmtVert=1
            // indicates a bifurcation
            for (int i = 0; i < 3; ++i)
            {
                if (iter.second[i]->m_elmtVert == 0)
                {
                    nbeg += 1;
                }
                else
                {
                    nend += 1;
                }
            }

            // Set up Bifurcation information
            if (nbeg == 2)
            {
                // ensure first InterfacePoint is parent
                if (iter.second[0]->m_elmtVert ==
                    1) // m_elmtVert: Vertex id in local element
                {
                    m_bifurcations.push_back(iter.second);
                }
                else
                {
                    // order points according to Riemann solver convention
                    InterfacePointShPtr I;
                    // find merging vessel
                    if (iter.second[1]->m_elmtVert == 1)
                    {
                        I = iter.second[0];
                        iter.second[0] = iter.second[1];
                        iter.second[1] = I;
                    }
                    else if (iter.second[2]->m_elmtVert == 1)
                    {
                        I = iter.second[0];
                        iter.second[0] = iter.second[2];
                        iter.second[2] = I;
                    }
                    NEKERROR(ErrorUtil::ewarning,
                             "This routine has not been checked");
                }
            }
            else
            {
                // ensure last InterfacePoint is merged vessel
                if (iter.second[0]->m_elmtVert == 0)
                {
                    m_mergingJcts.push_back(iter.second);
                }
                else
                {
                    // order points according to Riemann solver convention
                    InterfacePointShPtr I;
                    // find merging vessel
                    if (iter.second[1]->m_elmtVert == 0)
                    {
                        I = iter.second[0];
                        iter.second[0] = iter.second[1];
                        iter.second[1] = I;
                    }
                    else if (iter.second[2]->m_elmtVert == 0)
                    {
                        I = iter.second[0];
                        iter.second[0] = iter.second[2];
                        iter.second[2] = I;
                    }
                    NEKERROR(ErrorUtil::ewarning,
                             "This routine has not been checked");
                }
            }

        }
        else
        {
            ASSERTL0(false, "Unknown junction type");
        }
    }
}

  /**
   *  Initialisation routine for multiple subdomain case. Sets the
   *  initial conditions for all arterial subdomains read from the
   *  inputfile. Sets the material properties and the A_0 area for
   *  all subdomains and fills the domain-linking boundary
   *  conditions with the initial values of their domain.
   */
void PulseWaveSystem::v_DoInitialise()
{
    if (m_session->GetComm()->GetRank() == 0)
    {
        cout << "Initial Conditions: " << endl;
    }

    /* Loop over all subdomains to initialize all with the Initial
     * Conditions read from the inputfile*/
    for (int omega = 0; omega < m_nDomains; omega++)
    {
        for (int i = 0; i < 2; ++i)
        {
            m_fields[i] = m_vessels[m_nVariables * omega + i];
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << "Subdomain: " << omega << endl;
        }

        SetInitialConditions(0.0, 0, omega);
    }

    // Reset to first definition
    for (int i = 0; i < 2; ++i)
    {
        m_fields[i] = m_vessels[i];
    }
}

/**
 * NEEDS Updating:
 *
 *  DoSolve routine for PulseWavePropagation with multiple
 *  subdomains taken from UnsteadySystem and modified for
 *  multidomain case. Initialises the time integration scheme (as
 *  specified in the session file), and perform the time
 *  integration. Within the timestepping loop the following is
 *  done: 1. Link all arterial segments according to the network
 *  structure, solve the Riemann problem between different
 *  arterial segments and assign the values to the boundary
 *  conditions (LinkSubdomains) 2. Every arterial segment is
 *  solved independently for this timestep. This is done by handing
 *  the solution vector \f$ \mathbf{u} \f$ and the right hand side
 *  m_ode, which is the PulseWavePropagation class in this example
 *  over to the time integration scheme
 */
void PulseWaveSystem::v_DoSolve()
{
    NekDouble IntegrationTime = 0.0;
    int i;
    int n;
    int nchk = 1;

    Array<OneD, Array<OneD, NekDouble> > fields(m_nVariables);

    for (int i = 0; i < m_nVariables; ++i)
    {
        fields[i] = m_vessels[i]->UpdatePhys();
        m_fields[i]->SetPhysState(false);
    }

    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    // Time loop
    for (n = 0; n < m_steps; ++n)
    {
        LibUtilities::Timer timer;
        timer.Start();
        fields = m_intScheme->TimeIntegrate(n, m_timestep, m_ode);
        m_time += m_timestep;
        timer.Stop();
        IntegrationTime += timer.TimePerTest(1);

        // Write out status information.
        if (m_session->GetComm()->GetRank() == 0 && !((n + 1) % m_infosteps))
        {
            cout << "Steps: " << n + 1 << "\t Time: " << m_time
                 << "\t Time-step: " << m_timestep << "\t" << endl;
        }

        // Transform data if needed
        if (!((n + 1) % m_checksteps))
        {
            for (i = 0; i < m_nVariables; ++i)
            {
                int cnt = 0;
                for (int omega = 0; omega < m_nDomains; omega++)
                {
                    m_vessels[omega * m_nVariables + i]->FwdTrans(
                        fields[i] + cnt,
                        m_vessels[omega * m_nVariables + i]->UpdateCoeffs());
                    cnt += m_vessels[omega * m_nVariables + i]->GetTotPoints();
                }
            }
            CheckPoint_Output(nchk++);
        }

    } // end of timeintegration

    // Copy Array To Vessel Phys Fields
    for (int i = 0; i < m_nVariables; ++i)
    {
        Vmath::Vcopy(fields[i].size(), fields[i], 1,
                     m_vessels[i]->UpdatePhys(), 1);
    }

    cout << "Time-integration timing: " << IntegrationTime << " s" << endl
         << endl;
}

void PulseWaveSystem::FillDataFromInterfacePoint(
    InterfacePointShPtr &I,
    const Array<OneD, const Array<OneD, NekDouble> > &fields, NekDouble &A,
    NekDouble &u, NekDouble &beta, NekDouble &A_0, NekDouble &alpha)
{
    int omega       = I->m_domain;
    int traceId     = I->m_traceId;
    int eid         = I->m_elmt;
    int vert        = I->m_elmtVert;
    int vesselID    = omega * m_nVariables;
    int phys_offset = m_vessels[vesselID]->GetPhys_Offset(eid);
    LocalRegions::ExpansionSharedPtr dumExp;
    Array<OneD, NekDouble> Tmp(1);

    m_vessels[vesselID]->GetExp(eid)->GetTracePhysVals(
        vert, dumExp, fields[0] + m_fieldPhysOffset[omega] + phys_offset, Tmp);
    A = Tmp[0];
    m_vessels[vesselID]->GetExp(eid)->GetTracePhysVals(
        vert, dumExp, fields[1] + m_fieldPhysOffset[omega] + phys_offset, Tmp);
    u = Tmp[0];

    beta = m_beta_trace[omega][traceId];
    A_0 = m_A_0_trace[omega][traceId];
    alpha = m_alpha_trace[omega][traceId];
}

void PulseWaveSystem::EnforceInterfaceConditions(
    const Array<OneD, const Array<OneD, NekDouble> > &fields)
{
    int dom, bcpos;
    Array<OneD, NekDouble> Au(3);
    Array<OneD, NekDouble> uu(3);
    Array<OneD, NekDouble> beta(3);
    Array<OneD, NekDouble> A_0(3);
    Array<OneD, NekDouble> alpha(3);

    // Enforce Bifurcations:
    for (int n = 0; n < m_bifurcations.size(); ++n)
    {
        for (int i = 0; i < 3; ++i)
        {
            FillDataFromInterfacePoint(m_bifurcations[n][i], fields, Au[i],
                                              uu[i], beta[i], A_0[i], alpha[i]);
        }
        // Solve the Riemann problem for a bifurcation
        BifurcationRiemann(Au, uu, beta, A_0, alpha);

        // Store the values into the right positions:
        for (int i = 0; i < 3; ++i)
        {
            dom   = m_bifurcations[n][i]->m_domain;
            bcpos = m_bifurcations[n][i]->m_bcPosition;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[i];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[i];
        }
    }

    // Enforce Bifurcations;
    for (int n = 0; n < m_mergingJcts.size(); ++n)
    {
        // Merged vessel
        for (int i = 0; i < 3; ++i)
        {
            FillDataFromInterfacePoint(m_mergingJcts[n][i], fields, Au[i],
                                              uu[i], beta[i], A_0[i], alpha[i]);
        }

        // Solve the Riemann problem for a merging vessel
        MergingRiemann(Au, uu, beta, A_0, alpha);

        // Store the values into the right positions:
        for (int i = 0; i < 3; ++i)
        {
            int dom   = m_mergingJcts[n][i]->m_domain;
            int bcpos = m_mergingJcts[n][i]->m_bcPosition;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[i];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[i];
        }
    }

    for (int n = 0; n < m_vesselJcts.size(); ++n)
    {
        for (int i = 0; i < 2; ++i)
        {
            FillDataFromInterfacePoint(m_vesselJcts[n][i], fields, Au[i], uu[i],
                                       beta[i], A_0[i], alpha[i]);
        }

        JunctionRiemann(Au, uu, beta, A_0, alpha);

        // Store the values into the right positions:
        for (int i = 0; i < 2; ++i)
        {
            int dom   = m_vesselJcts[n][i]->m_domain;
            int bcpos = m_vesselJcts[n][i]->m_bcPosition;
            m_vessels[dom * m_nVariables]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = Au[i];
            m_vessels[dom * m_nVariables + 1]
                ->UpdateBndCondExpansion(bcpos)
                ->UpdatePhys()[0] = uu[i];
        }
    }
}

/**
 *  Solves the Riemann problem at a bifurcation by assuming
 *  subsonic flow at both sides of the boundary and by
 *  applying conservation of mass and continuity of the total
 *  pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The
 *  other 3 missing equations come from the characteristic
 *  variables. For further information see "Pulse
 *  WavePropagation in the human vascular system" Section
 *  3.4.4
 */
void PulseWaveSystem::BifurcationRiemann(Array<OneD, NekDouble> &Au,
                                         Array<OneD, NekDouble> &uu,
                                         Array<OneD, NekDouble> &beta,
                                         Array<OneD, NekDouble> &A_0,
                                         Array<OneD, NekDouble> &alpha)
{
    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(3);
    Array<OneD, NekDouble> P_Au(3);
    Array<OneD, NekDouble> W_Au(3);
    NekMatrix<NekDouble> invJ(6, 6);
    NekVector<NekDouble> f(6);
    NekVector<NekDouble> dx(6);

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 100;

    // Forward characteristic
    m_pressureArea->GetW1(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);

    // Backward characteristics
    for (int i = 1; i < 3; ++i)
    {
        m_pressureArea->GetW2(W[i], uu[i], beta[i], Au[i], A_0[i], alpha[i]);
    }

    // Tolerances for the algorithm
    NekDouble Tol = 1.0E-10;

    // Newton Iteration
    while ((proceed) && (iter < MAX_ITER))
    {
        iter += 1;

        /*
         * We solve the six constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Forward characteristic:          W1(A_L,  U_L)  = W1(Au_L,  Uu_L)
         * 2. Backward characteristic 1:       W2(A_R1, U_R1) = W2(Au_R1, Uu_R1)
         * 3. Backward characteristic 2:       W2(A_R2, U_R2) = W2(Au_R2, Uu_R2)
         * 4. Conservation of mass:            Au_L * Uu_L    = Au_R1 * Uu_R1 +
         *                                                      Au_R2 * Uu_R2
         * 5. Continuity of total pressure 1:  rho * Uu_L  * Uu_L  / 2 + p(Au_L) =
         *                                     rho * Uu_R1 * Uu_R1 / 2 + p(Au_R1)
         * 6. Continuity of total pressure 2:  rho * Uu_L  * Uu_L  / 2 + p(Au_L) =
         *                                     rho * Uu_R2 * Uu_R2 / 2 + p(Au_R2)
         */
        for (int i = 0; i < 3; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                                                      alpha[i]);
        }

        m_pressureArea->GetW1(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW2(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);
        m_pressureArea->GetW2(W_Au[2], uu[2], beta[2], Au[2], A_0[2], alpha[2]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = W_Au[2] - W[2];
        f[3] = Au[0] * uu[0] - Au[1] * uu[1] - Au[2] * uu[2];
        f[4] = uu[0] * uu[0] + 2 * P_Au[0] / rho -
               uu[1] * uu[1] - 2 * P_Au[1] / rho;
        f[5] = uu[0] * uu[0] + 2 * P_Au[0] / rho -
               uu[2] * uu[2] - 2 * P_Au[2] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                                                 "Bifurcation");

        Multiply(dx, invJ, f);

        // Update the solution: x_new = x_old - dx
        for (int i = 0; i < 3; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 3];
        }

        // Check if the error of the solution is smaller than Tol
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }

        // Check if solver converges
        if (iter >= MAX_ITER)
        {
            ASSERTL0(false, "Riemann solver for Bifurcation did not converge.");
        }
    }
}

/**
 *  Solves the Riemann problem at an merging flow condition by
 *  assuming subsonic flow at both sides of the boundary and by
 *  applying conservation of mass and continuity of the total
 *  pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$ The other 3
 *  missing equations come from the characteristic variables. For
 *  further information see "Pulse WavePropagation in the human
 *  vascular system" Section 3.4.4
 */
void PulseWaveSystem::MergingRiemann(Array<OneD, NekDouble> &Au,
                                     Array<OneD, NekDouble> &uu,
                                     Array<OneD, NekDouble> &beta,
                                     Array<OneD, NekDouble> &A_0,
                                     Array<OneD, NekDouble> &alpha)
{
    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(3);
    Array<OneD, NekDouble> W_Au(3);
    Array<OneD, NekDouble> P_Au(3);
    NekMatrix<NekDouble> invJ(6, 6);
    NekVector<NekDouble> f(6);
    NekVector<NekDouble> dx(6);

    int proceed  = 1;
    int iter     = 0;
    int MAX_ITER = 15;

    // Backward characteristic
    m_pressureArea->GetW2(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);

    // Forward characteristics
    for (int i = 1; i < 3; ++i)
    {
        m_pressureArea->GetW1(W[i], uu[i], beta[i], Au[i], A_0[i], alpha[i]);
    }

    // Tolerances for the algorithm
    NekDouble Tol = 1.0E-10;

    // Newton Iteration
    while ((proceed) && (iter < MAX_ITER))
    {
        iter += 1;

        /*
         * We solve the six constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Backward characteristic:          W2(A_R,  U_R)  = W1(Au_R, Uu_R)
         * 2. Forward characteristic 1:         W1(A_L1, U_L1) = W1(Au_L1, Uu_R1)
         * 3. Forward characteristic 2:         W1(A_L2, U_L2) = W1(Au_L2, Uu_L2)
         * 4. Conservation of mass:             Au_R * Uu_R    = Au_L1 * Uu_L1 +
         *                                                       Au_L2 * Uu_L2
         * 5. Continuity of total pressure 1:  rho * Uu_R  * Uu_R  / 2 + p(Au_R) =
         *                                     rho * Uu_L1 * Uu_L1 / 2 + p(Au_L1)
         * 6. Continuity of total pressure 2:  rho * Uu_R  * Uu_R  / 2 + p(Au_R) =
         *                                     rho * Uu_L2 * Uu_L2 / 2 + p(Au_L2)
         */
        for (int i = 0; i < 3; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                                                      alpha[i]);
        }

        m_pressureArea->GetW2(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW1(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);
        m_pressureArea->GetW1(W_Au[2], uu[2], beta[2], Au[2], A_0[2], alpha[2]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = W_Au[2] - W[2];
        f[3] = Au[0] * uu[0] - Au[1] * uu[1] - Au[2] * uu[2];
        f[4] = uu[0] * uu[0] + 2 * P_Au[0] / rho -
               uu[1] * uu[1] - 2 * P_Au[1] / rho;
        f[5] = uu[0] * uu[0] + 2 * P_Au[0] / rho -
               uu[2] * uu[2] - 2 * P_Au[2] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                                                       "Merge");

        Multiply(dx, invJ, f);

        // Update the solution: x_new = x_old - dx
        for (int i = 0; i < 3; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 3];
        }

        // Check if the error of the solution is smaller than Tol
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }

        // Check if solver converges
        if (iter >= MAX_ITER)
        {
            ASSERTL0(false, "Riemann solver for Merging Flow did not converge");
        }
    }
}

/**
 *  Solves the Riemann problem at an interdomain junction by
 *  assuming subsonic flow at both sides of the boundary and
 *  by applying conservation of mass and continuity of the
 *  total pressure \f$ \frac{p}{rho} + \frac{u^{2}}{2}. \f$
 *  The other 2 missing equations come from the characteristic
 *  variables. For further information see "Pulse
 *  WavePropagation in the human vascular system" Section 3.4.
 */
void PulseWaveSystem::JunctionRiemann(Array<OneD, NekDouble> &Au,
                                      Array<OneD, NekDouble> &uu,
                                      Array<OneD, NekDouble> &beta,
                                      Array<OneD, NekDouble> &A_0,
                                      Array<OneD, NekDouble> &alpha)
{
    NekDouble rho = m_rho;
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> W_Au(2);
    Array<OneD, NekDouble> P_Au(2);
    NekMatrix<NekDouble> invJ(4, 4);
    NekVector<NekDouble> f(4);
    NekVector<NekDouble> dx(4);

    int proceed   = 1;
    int iter      = 0;
    int MAX_ITER  = 15;
    NekDouble Tol = 1.0E-10;

    // Forward and backward characteristics
    m_pressureArea->GetW1(W[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
    m_pressureArea->GetW2(W[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);

    while ((proceed) && (iter < MAX_ITER))
    {
        iter += 1;

        /*
         * We solve the four constraint equations via a multivariate Newton
         * iteration. Equations are:
         * 1. Forward characteristic:        W1(A_L, U_L) = W1(Au_L, Uu_L)
         * 2. Backward characteristic:       W2(A_R, U_R) = W2(Au_R, Uu_R)
         * 3. Conservation of mass:          Au_L * Uu_L = Au_R * Uu_R
         * 4. Continuity of total pressure:  rho * Uu_L * Uu_L / 2 + p(Au_L) =
         *                                   rho * Uu_R * Uu_R / 2 + p(Au_R)
         */
        for (int i = 0; i < 2; ++i)
        {
            m_pressureArea->GetPressure(P_Au[i], beta[i], Au[i], A_0[i], 0, 0,
                                                                      alpha[i]);
        }

        m_pressureArea->GetW1(W_Au[0], uu[0], beta[0], Au[0], A_0[0], alpha[0]);
        m_pressureArea->GetW2(W_Au[1], uu[1], beta[1], Au[1], A_0[1], alpha[1]);

        // Constraint equations set to zero
        f[0] = W_Au[0] - W[0];
        f[1] = W_Au[1] - W[1];
        f[2] = Au[0] * uu[0] - Au[1] * uu[1];
        f[3] = uu[0] * uu[0] + 2 * P_Au[0] / rho -
               uu[1] * uu[1] - 2 * P_Au[1] / rho;

        // Inverse Jacobian for x + dx = x - J^(-1) * f(x)
        m_pressureArea->GetJacobianInverse(invJ, Au, uu, beta, A_0, alpha,
                                                                    "Junction");

        Multiply(dx, invJ, f);

        // Update solution: x_new = x_old - dx
        for (int i = 0; i < 2; ++i)
        {
            uu[i] -= dx[i];
            Au[i] -= dx[i + 2];
        }

        // Check if the error of the solution is smaller than Tol.
        if (Dot(dx, dx) < Tol)
        {
            proceed = 0;
        }
    }

    if (iter >= MAX_ITER)
    {
        ASSERTL0(false, "Riemann solver for Junction did not converge");
    }
}

/**
 *  Writes the .fld file at the end of the simulation. Similar to the normal
 *  v_Output however the Multidomain output has to be prepared.
 */
void PulseWaveSystem::v_Output(void)
{
    /**
    * Write the field data to file. The file is named according to the session
    * name with the extension .fld appended.
    */
    std::string outname = m_sessionName + ".fld";

    WriteVessels(outname);
}

/**
 *  Writes the .fld file at the end of the simulation. Similar to the normal
 *  v_Output however the Multidomain output has to be prepared.
 */
void PulseWaveSystem::CheckPoint_Output(const int n)
{
    std::stringstream outname;
    outname << m_sessionName << "_" << n << ".chk";

    WriteVessels(outname.str());
}

/**
 * Writes the field data to a file with the given filename.
 * @param   outname         Filename to write to.
 */
void PulseWaveSystem::WriteVessels(const std::string &outname)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::string> variables = m_session->GetVariables();

    for (int n = 0; n < m_nDomains; ++n)
    {
        m_vessels[n * m_nVariables]->GetFieldDefinitions(FieldDef);
    }

    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    int nFieldDefPerDomain = FieldDef.size() / m_nDomains;
    int cnt;

    // Copy Data into FieldData and set variable
    for (int n = 0; n < m_nDomains; ++n)
    {
        // Outputs area and velocity
        for (int j = 0; j < m_nVariables; ++j)
        {
            for (int i = 0; i < nFieldDefPerDomain; ++i)
            {
                cnt = n * nFieldDefPerDomain + i;

                FieldDef[cnt]->m_fields.push_back(variables[j]);

                m_vessels[n * m_nVariables]->AppendFieldData(
                    FieldDef[cnt], FieldData[cnt],
                    m_vessels[n * m_nVariables + j]->UpdateCoeffs());
            }
        }

        // Outputs pressure
        Array<OneD, NekDouble> PFwd(m_vessels[n * m_nVariables]->GetNcoeffs());

        m_vessels[n * m_nVariables]->FwdTrans_IterPerExp(m_pressure[n], PFwd);

        for (int i = 0; i < nFieldDefPerDomain; ++i)
        {
            cnt = n * nFieldDefPerDomain + i;

            FieldDef[cnt]->m_fields.push_back("P");

            m_vessels[n * m_nVariables]->AppendFieldData(FieldDef[cnt],
                                                          FieldData[cnt], PFwd);
        }

        if (extraFields)
        {
            Array<OneD, NekDouble>  PWVFwd(m_vessels[n *
                                                   m_nVariables]->GetNcoeffs());
            Array<OneD, NekDouble>   W1Fwd(m_vessels[n *
                                                   m_nVariables]->GetNcoeffs());
            Array<OneD, NekDouble>   W2Fwd(m_vessels[n *
                                                   m_nVariables]->GetNcoeffs());

            m_vessels[n * m_nVariables]->FwdTrans_IterPerExp(m_PWV[n], PWVFwd);
            m_vessels[n * m_nVariables]->FwdTrans_IterPerExp(m_W1[n], W1Fwd);
            m_vessels[n * m_nVariables]->FwdTrans_IterPerExp(m_W2[n], W2Fwd);

            for (int i = 0; i < nFieldDefPerDomain; ++i)
            {
                cnt = n * nFieldDefPerDomain + i;

                FieldDef[cnt]->m_fields.push_back("c");
                FieldDef[cnt]->m_fields.push_back("W1");
                FieldDef[cnt]->m_fields.push_back("W2");

                m_vessels[n * m_nVariables]->AppendFieldData(FieldDef[cnt],
                                                        FieldData[cnt], PWVFwd);
                m_vessels[n * m_nVariables]->AppendFieldData(FieldDef[cnt],
                                                         FieldData[cnt], W1Fwd);
                m_vessels[n * m_nVariables]->AppendFieldData(FieldDef[cnt],
                                                         FieldData[cnt], W2Fwd);
            }
        }
    }

    // Update time in field info if required
    if (m_fieldMetaDataMap.find("Time") != m_fieldMetaDataMap.end())
    {
        m_fieldMetaDataMap["Time"] = boost::lexical_cast<std::string>(m_time);
    }

    LibUtilities::Write(outname, FieldDef, FieldData, m_fieldMetaDataMap);
}

/* Compute the error in the L2-norm
 * @param   field           The field to compare.
 * @param   exactsoln       The exact solution to compare with.
 * @param   Normalised      Normalise L2-error.
 * @returns                 Error in the L2-norm.
 */
NekDouble PulseWaveSystem::v_L2Error(unsigned int field,
                                     const Array<OneD, NekDouble> &exactsoln,
                                     bool Normalised)
{
    NekDouble L2error = 0.0;
    NekDouble L2error_dom;
    NekDouble Vol = 0.0;

    if (m_NumQuadPointsError == 0)
    {
        for (int omega = 0; omega < m_nDomains; omega++)
        {
            int vesselid = field + omega * m_nVariables;

            if (m_vessels[vesselid]->GetPhysState() == false)
            {
                m_vessels[vesselid]->BwdTrans(m_vessels[vesselid]->GetCoeffs(),
                                             m_vessels[vesselid]->UpdatePhys());
            }

            if (exactsoln.size())
            {
                L2error_dom = m_vessels[vesselid]->L2(
                    m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else if (m_session->DefinesFunction("ExactSolution"))
            {
                Array<OneD, NekDouble> exactsoln(
                 m_vessels[vesselid]->GetNpoints());

                LibUtilities::EquationSharedPtr vEqu =
                  m_session->GetFunction("ExactSolution", field, omega);
                GetFunction("ExactSolution")
                  ->Evaluate(m_session->GetVariable(field), exactsoln,
                             m_time);

                L2error_dom = m_vessels[vesselid]->L2(
                  m_vessels[vesselid]->GetPhys(), exactsoln);
            }
            else
            {
                L2error_dom =
                        m_vessels[vesselid]->L2(m_vessels[vesselid]->GetPhys());
            }

            L2error += L2error_dom * L2error_dom;

            if (Normalised == true)
            {
                Array<OneD, NekDouble> one(m_vessels[vesselid]->GetNpoints(),
                                                                           1.0);

                Vol += m_vessels[vesselid]->PhysIntegral(one);
            }
        }
    }
    else
    {
        ASSERTL0(false, "Not set up");
    }

    if (Normalised == true)
    {
        m_comm->AllReduce(Vol, LibUtilities::ReduceSum);

        L2error = sqrt(L2error / Vol);
    }
    else
    {
        L2error = sqrt(L2error);
    }

    return L2error;
}

/**
 * Compute the error in the L_inf-norm
 * @param   field           The field to compare.
 * @param   exactsoln       The exact solution to compare with.
 * @returns                 Error in the L_inft-norm.
 */
NekDouble PulseWaveSystem::v_LinfError(unsigned int field,
                             const Array<OneD, NekDouble> &exactsoln)
{
    NekDouble LinferrorDom, Linferror = -1.0;

    for (int omega = 0; omega < m_nDomains; ++omega)
    {
        int vesselid = field + omega * m_nVariables;

        if (m_NumQuadPointsError == 0)
        {
            if (m_vessels[vesselid]->GetPhysState() == false)
            {
                m_vessels[vesselid]->BwdTrans(m_vessels[vesselid]->GetCoeffs(),
                                            m_vessels[vesselid]->UpdatePhys());
            }

            if (exactsoln.size())
            {
                LinferrorDom =
                       m_vessels[vesselid]->Linf(m_vessels[vesselid]->GetPhys(),
                                                                     exactsoln);
            }
            else if (m_session->DefinesFunction("ExactSolution"))
            {
                Array<OneD, NekDouble>
                                   exactsoln(m_vessels[vesselid]->GetNpoints());

                GetFunction("ExactSolution")
                  ->Evaluate(m_session->GetVariable(field), exactsoln, m_time);

                LinferrorDom =
                       m_vessels[vesselid]->Linf(m_vessels[vesselid]->GetPhys(),
                                                                     exactsoln);
            }
            else
            {
                LinferrorDom = 0.0;
            }

            Linferror = (Linferror > LinferrorDom) ? Linferror : LinferrorDom;
        }
        else
        {
            ASSERTL0(false, "ErrorExtraPoints not allowed for this solver");
        }
    }
    return Linferror;
}

} // namespace Nektar
