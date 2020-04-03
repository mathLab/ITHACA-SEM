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

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    using namespace MultiRegions;

    string VelocityCorrectionScheme::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "VelocityCorrectionScheme",
            VelocityCorrectionScheme::create);

    /**
     * Constructor. Creates ...
     *
     * \param
     * \param
     */
    VelocityCorrectionScheme::VelocityCorrectionScheme(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          IncNavierStokes(pSession, pGraph),
          m_varCoeffLap(StdRegions::NullVarCoeffMap)
    {

    }

    void VelocityCorrectionScheme::v_InitObject()
    {
        int n;

        IncNavierStokes::v_InitObject();
        m_explicitDiffusion = false;

        // Set m_pressure to point to last field of m_fields;
        if (boost::iequals(m_session->GetVariable(m_fields.size()-1), "p"))
        {
            m_nConvectiveFields = m_fields.size()-1;
            m_pressure = m_fields[m_nConvectiveFields];
        }
        else
        {
            ASSERTL0(false,"Need to set up pressure field definition");
        }

        // Determine diffusion coefficients for each field
        m_diffCoeff = Array<OneD, NekDouble> (m_nConvectiveFields, m_kinvis);
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            std::string varName = m_session->GetVariable(n);
            if ( m_session->DefinesFunction("DiffusionCoefficient", varName))
            {
                LibUtilities::EquationSharedPtr ffunc
                    = m_session->GetFunction("DiffusionCoefficient", varName);
                m_diffCoeff[n] = ffunc->Evaluate();
            }
        }

        // Integrate only the convective fields
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            m_intVariables.push_back(n);
        }

        SetUpExtrapolation();
        SetUpSVV();

        m_session->MatchSolverInfo("SmoothAdvection", "True",
                                    m_SmoothAdvection, false);

        // set explicit time-intregration class operators
        m_ode.DefineOdeRhs(
            &VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs, this);

        // set implicit time-intregration class operators
        m_ode.DefineImplicitSolve(
            &VelocityCorrectionScheme::SolveUnsteadyStokesSystem, this);

        // Set up bits for flowrate.
        m_session->LoadParameter("Flowrate", m_flowrate, 0.0);
        m_session->LoadParameter("IO_FlowSteps", m_flowrateSteps, 0);
    }

    void VelocityCorrectionScheme::SetUpExtrapolation()
    {
        // creation of the extrapolation object
        if (m_equationType == eUnsteadyNavierStokes ||
            m_equationType == eUnsteadyStokes)
        {
            std::string vExtrapolation = v_GetExtrapolateStr();
            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = v_GetSubSteppingExtrapolateStr(
                    m_session->GetSolverInfo("Extrapolation"));
            }
            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation,
                m_session,
                m_fields,
                m_pressure,
                m_velocity,
                m_advObject);

            m_extrapolation->SubSteppingTimeIntegration(m_intScheme);
            m_extrapolation->GenerateHOPBCMap(m_session);
        }
    }

    /**
     * @brief Set up the Stokes solution used to impose constant flowrate
     * through a boundary.
     *
     * This routine solves a Stokes equation using a unit forcing direction,
     * specified by the user to be in the desired flow direction. This field can
     * then be used to correct the end of each timestep to impose a constant
     * volumetric flow rate through a user-defined boundary.
     *
     * There are three modes of operation:
     *
     * - Standard two-dimensional or three-dimensional simulations (e.g. pipes
     *   or channels)
     * - 3DH1D simulations where the forcing is not in the homogeneous
     *   direction (e.g. channel flow, where the y-direction of the 2D mesh
     *   is perpendicular to the wall);
     * - 3DH1D simulations where the forcing is in the homogeneous direction
     *   (e.g. pipe flow in the z-direction).
     *
     * In the first two cases, the user should define:
     * - the `Flowrate` parameter, which dictates the volumetric flux through
     *   the reference area
     * - tag a boundary region with the `Flowrate` user-defined type to define
     *   the reference area
     * - define a `FlowrateForce` function with components `ForceX`, `ForceY`
     *   and `ForceZ` that defines a unit forcing in the appropriate direction.
     *
     * In the latter case, the user should define only the `Flowrate`; the
     * reference area is taken to be the homogeneous plane and the force is
     * assumed to be the unit z-vector \f$ \hat{e}_z \f$.
     *
     * This routine solves a single timestep of the Stokes problem
     * (premultiplied by the backwards difference coefficient):
     *
     * \f[ \frac{\partial\mathbf{u}}{\partial t} = -\nabla p +
     * \nu\nabla^2\mathbf{u} + \mathbf{f} \f]
     *
     * with a zero initial condition to obtain a field \f$ \mathbf{u}_s \f$. The
     * flowrate is then corrected at each timestep \f$ n \f$ by adding the
     * correction \f$ \alpha\mathbf{u}_s \f$ where
     *
     * \f[ \alpha = \frac{\overline{Q} - Q(\mathbf{u^n})}{Q(\mathbf{u}_s)} \f]
     *
     * where \f$ Q(\cdot)\f$ is the volumetric flux through the appropriate
     * surface or line, which is implemented in
     * VelocityCorrectionScheme::MeasureFlowrate. For more details, see chapter
     * 3.2 of the thesis of D. Moxey (University of Warwick, 2011).
     */
    void VelocityCorrectionScheme::SetupFlowrate(NekDouble aii_dt)
    {
        m_flowrateBndID = -1;
        m_flowrateArea = 0.0;

        const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bcs =
            m_fields[0]->GetBndConditions();

        std::string forces[] = { "X", "Y", "Z" };
        Array<OneD, NekDouble> flowrateForce(m_spacedim, 0.0);

        // Set up flowrate forces.
        bool defined = true;
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string varName = std::string("Force") + forces[i];
            defined = m_session->DefinesFunction("FlowrateForce", varName);

            if (!defined && m_HomogeneousType == eHomogeneous1D)
            {
                break;
            }

            ASSERTL0(defined,
                     "A 'FlowrateForce' function must defined with components "
                     "[ForceX, ...] to define direction of flowrate forcing");

            LibUtilities::EquationSharedPtr ffunc
                = m_session->GetFunction("FlowrateForce", varName);
            flowrateForce[i] = ffunc->Evaluate();
        }

        // Define flag for case with homogeneous expansion and forcing not in the
        // z-direction
        m_homd1DFlowinPlane = false;
        if (defined && m_HomogeneousType == eHomogeneous1D)
        {
            m_homd1DFlowinPlane = true;
        }

        // For 3DH1D simulations, if force isn't defined then assume in
        // z-direction.
        if (!defined)
        {
            flowrateForce[2] = 1.0;
        }

        // Find the boundary condition that is tagged as the flowrate boundary.
        for (int i = 0; i < bcs.size(); ++i)
        {
            if (boost::iequals(bcs[i]->GetUserDefined(), "Flowrate"))
            {
                m_flowrateBndID = i;
                break;
            }
        }

        int tmpBr = m_flowrateBndID;
        m_comm->AllReduce(tmpBr, LibUtilities::ReduceMax);
        ASSERTL0(tmpBr >= 0 || m_HomogeneousType == eHomogeneous1D,
                 "One boundary region must be marked using the 'Flowrate' "
                 "user-defined type to monitor the volumetric flowrate.");

        // Extract an appropriate expansion list to represents the boundary.
        if (m_flowrateBndID >= 0)
        {
            // For a boundary, extract the boundary itself.
            m_flowrateBnd = m_fields[0]->GetBndCondExpansions()[m_flowrateBndID];
        }
        else if (m_HomogeneousType == eHomogeneous1D && !m_homd1DFlowinPlane)
        {
            // For 3DH1D simulations with no force specified, find the mean
            // (0th) plane.
            Array<OneD, unsigned int> zIDs = m_fields[0]->GetZIDs();
            int tmpId = -1;

            for (int i = 0; i < zIDs.size(); ++i)
            {
                if (zIDs[i] == 0)
                {
                    tmpId = i;
                    break;
                }
            }

            ASSERTL1(tmpId <= 0, "Should be either at location 0 or -1 if not "
                                 "found");

            if (tmpId != -1)
            {
                m_flowrateBnd = m_fields[0]->GetPlane(tmpId);
            }
        }

        // At this point, some processors may not have m_flowrateBnd set if they
        // don't contain the appropriate boundary. To calculate the area, we
        // integrate 1.0 over the boundary (which has been set up with the
        // appropriate subcommunicator to avoid deadlock), and then communicate
        // this to the other processors with an AllReduce.
        if (m_flowrateBnd)
        {
            Array<OneD, NekDouble> inArea(m_flowrateBnd->GetNpoints(), 1.0);
            m_flowrateArea = m_flowrateBnd->Integral(inArea);
        }
        m_comm->AllReduce(m_flowrateArea, LibUtilities::ReduceMax);

        // In homogeneous case with forcing not aligned to the z-direction,
        // redefine m_flowrateBnd so it is a 1D expansion
        if (m_HomogeneousType == eHomogeneous1D && m_homd1DFlowinPlane &&
            m_flowrateBnd)
        {
            // For 3DH1D simulations with no force specified, find the mean
            // (0th) plane.
            Array<OneD, unsigned int> zIDs = m_fields[0]->GetZIDs();
            m_planeID = -1;

            for (int i = 0; i < zIDs.size(); ++i)
            {
                if (zIDs[i] == 0)
                {
                    m_planeID = i;
                    break;
                }
            }

            ASSERTL1(m_planeID <= 0, "Should be either at location 0 or -1 if not "
                                 "found");

            if (m_planeID != -1)
            {
                m_flowrateBnd = m_fields[0]
                                    ->GetBndCondExpansions()[m_flowrateBndID]
                                    ->GetPlane(m_planeID);
            }
        }

        // Set up some storage for the Stokes solution (to be stored in
        // m_flowrateStokes) and its initial condition (inTmp), which holds the
        // unit forcing.
        int nqTot = m_fields[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble> > inTmp(m_spacedim);
        m_flowrateStokes = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

        for (int i = 0; i < m_spacedim; ++i)
        {
            inTmp[i] = Array<OneD, NekDouble>(
                nqTot, flowrateForce[i] * aii_dt);
            m_flowrateStokes[i] = Array<OneD, NekDouble>(nqTot, 0.0);

            if (m_HomogeneousType == eHomogeneous1D)
            {
                Array<OneD, NekDouble> inTmp2(nqTot);
                m_fields[i]->HomogeneousFwdTrans(inTmp[i], inTmp2);
                m_fields[i]->SetWaveSpace(true);
                inTmp[i] = inTmp2;
            }

            Vmath::Zero(
                m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(), 1);
        }

        // Create temporary extrapolation object to avoid issues with
        // m_extrapolation for HOPBCs using higher order timestepping schemes.
        ExtrapolateSharedPtr tmpExtrap = m_extrapolation;
        m_extrapolation = GetExtrapolateFactory().CreateInstance(
            "Standard", m_session, m_fields, m_pressure, m_velocity,
            m_advObject);

        // Finally, calculate the solution and the flux of the Stokes
        // solution. We set m_greenFlux to maximum numeric limit, which signals
        // to SolveUnsteadyStokesSystem that we don't need to apply a flowrate
        // force.
        m_greenFlux = numeric_limits<NekDouble>::max();
        m_flowrateAiidt = aii_dt;
        SolveUnsteadyStokesSystem(inTmp, m_flowrateStokes, 0.0, aii_dt);
        m_greenFlux = MeasureFlowrate(m_flowrateStokes);

        // If the user specified IO_FlowSteps, open a handle to store output.
        if (m_comm->GetRank() == 0 && m_flowrateSteps &&
            !m_flowrateStream.is_open())
        {
            std::string filename = m_session->GetSessionName();
            filename += ".prs";
            m_flowrateStream.open(filename.c_str());
            m_flowrateStream.setf(ios::scientific, ios::floatfield);
            m_flowrateStream << "# step      time            dP" << endl
                             << "# -------------------------------------------"
                             << endl;
        }

        m_extrapolation = tmpExtrap;
    }

    /**
     * @brief Measure the volumetric flow rate through the volumetric flow rate
     * reference surface.
     *
     * This routine computes the volumetric flow rate
     *
     * \f[
     * Q(\mathbf{u}) = \frac{1}{\mu(R)} \int_R \mathbf{u} \cdot d\mathbf{s}
     * \f]
     *
     * through the boundary region \f$ R \f$.
     */
    NekDouble VelocityCorrectionScheme::MeasureFlowrate(
        const Array<OneD, Array<OneD, NekDouble> > &inarray)
    {
        NekDouble flowrate = 0.0;

        if (m_flowrateBnd && m_flowrateBndID >= 0)
        {
            // If we're an actual boundary, calculate the vector flux through
            // the boundary.
            Array<OneD, Array<OneD, NekDouble> > boundary(m_spacedim);

            if(!m_homd1DFlowinPlane)
            {
                // General case
                for (int i = 0; i < m_spacedim; ++i)
                {
                    m_fields[i]->ExtractPhysToBnd(m_flowrateBndID, inarray[i],
                                                  boundary[i]);
                }
                flowrate = m_flowrateBnd->VectorFlux(boundary);
            }
            else if(m_planeID == 0)
            {
                //Homogeneous with forcing in plane. Calculate flux only on
                // the meanmode - calculateFlux necessary for hybrid
                // parallelisation.
                for (int i = 0; i < m_spacedim; ++i)
                {
                    m_fields[i]->GetPlane(m_planeID)->ExtractPhysToBnd(
                        m_flowrateBndID, inarray[i], boundary[i]);
                }

                // the flowrate is calculated on the mean mode so it needs to be
                // multiplied by LZ to be consistent with the general case.
                flowrate = m_flowrateBnd->VectorFlux(boundary) *
                           m_session->GetParameter("LZ");
            }
        }
        else if (m_flowrateBnd && !m_homd1DFlowinPlane)
        {
            // 3DH1D case with no Flowrate boundary defined: compute flux
            // through the zero-th (mean) plane.
            flowrate = m_flowrateBnd->Integral(inarray[2]);
        }

        // Communication to obtain the total flowrate
        if(!m_homd1DFlowinPlane && m_HomogeneousType == eHomogeneous1D)
        {
            m_comm->GetColumnComm()->AllReduce(flowrate, LibUtilities::ReduceSum);
        }
        else
        {
            m_comm->AllReduce(flowrate, LibUtilities::ReduceSum);
        }
        return flowrate / m_flowrateArea;
    }

    bool VelocityCorrectionScheme::v_PostIntegrate(int step)
    {
        if (m_flowrateSteps > 0)
        {
            if (m_comm->GetRank() == 0 && (step + 1) % m_flowrateSteps == 0)
            {
                m_flowrateStream << setw(8) << step << setw(16) << m_time
                                 << setw(16) << m_alpha << endl;
            }
        }

        return IncNavierStokes::v_PostIntegrate(step);
    }


    /**
     * Destructor
     */
    VelocityCorrectionScheme::~VelocityCorrectionScheme(void)
    {
    }

    /**
     *
     */
    void VelocityCorrectionScheme::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s,
                "Splitting Scheme", "Velocity correction (strong press. form)");

        if( m_extrapolation->GetSubStepName().size() )
        {
            SolverUtils::AddSummaryItem( s, "Substepping",
                                         m_extrapolation->GetSubStepName() );
        }

        string dealias = m_homogen_dealiasing ? "Homogeneous1D" : "";
        if (m_specHP_dealiasing)
        {
            dealias += (dealias == "" ? "" : " + ") + string("spectral/hp");
        }
        if (dealias != "")
        {
            SolverUtils::AddSummaryItem(s, "Dealiasing", dealias);
        }


        string smoothing = m_useSpecVanVisc ? "spectral/hp" : "";
        if (smoothing != "")
        {
            if(m_svvVarDiffCoeff == NullNekDouble1DArray)
            {
                SolverUtils::AddSummaryItem(
                   s, "Smoothing-SpecHP", "SVV (" + smoothing +
                   " Exp Kernel(cut-off = "
                   + boost::lexical_cast<string>(m_sVVCutoffRatio)
                   + ", diff coeff = "
                   + boost::lexical_cast<string>(m_sVVDiffCoeff)+"))");
            }
            else
            {
                if(m_IsSVVPowerKernel)
                {
                    SolverUtils::AddSummaryItem(
                       s, "Smoothing-SpecHP", "SVV (" + smoothing +
                       " Power Kernel (Power ratio ="
                       + boost::lexical_cast<string>(m_sVVCutoffRatio)
                       + ", diff coeff = "
                       + boost::lexical_cast<string>(m_sVVDiffCoeff)+"*Uh/p))");
                }
                else
                {
                    SolverUtils::AddSummaryItem(
                       s, "Smoothing-SpecHP", "SVV (" + smoothing +
                       " DG Kernel (diff coeff = "
                       + boost::lexical_cast<string>(m_sVVDiffCoeff)+"*Uh/p))");

                }
            }

        }

        if (m_useHomo1DSpecVanVisc && (m_HomogeneousType == eHomogeneous1D))
        {
            SolverUtils::AddSummaryItem(
                  s, "Smoothing-Homo1D", "SVV (Homogeneous1D - Exp Kernel(cut-off = "
                  + boost::lexical_cast<string>(m_sVVCutoffRatioHomo1D)
                  + ", diff coeff = "
                  + boost::lexical_cast<string>(m_sVVDiffCoeffHomo1D)+"))");
        }

    }

    /**
     *
     */
    void VelocityCorrectionScheme::v_DoInitialise(void)
    {
        m_F = Array<OneD, Array<OneD, NekDouble> > (m_nConvectiveFields);

        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_F[i] = Array< OneD, NekDouble> (m_fields[0]->GetTotPoints(), 0.0);
        }

        m_flowrateAiidt = 0.0;

        AdvectionSystem::v_DoInitialise();

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   =
                boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] =
                boost::lexical_cast<std::string>(m_timestep);

        // set boundary conditions here so that any normal component
        // correction are imposed before they are imposed on initial
        // field below
        SetBoundaryConditions(m_time);

	// Ensure the initial conditions have correct BCs  
        for(int i = 0; i < m_fields.size(); ++i)
        {
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
	    m_fields[i]->LocalToGlobal();
	    m_fields[i]->GlobalToLocal();
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
        }
    }


    /**
     *
     */
    void VelocityCorrectionScheme:: v_TransCoeffToPhys(void)
    {
        int nfields = m_fields.size() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Backward Transformation in physical space for time evolution
            m_fields[k]->BwdTrans_IterPerExp(m_fields[k]->GetCoeffs(),
                                             m_fields[k]->UpdatePhys());
        }
    }

    /**
     *
     */
    void VelocityCorrectionScheme:: v_TransPhysToCoeff(void)
    {

        int nfields = m_fields.size() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Forward Transformation in physical space for time evolution
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),
                                             m_fields[k]->UpdateCoeffs());
        }
    }

    /**
     *
     */
    Array<OneD, bool> VelocityCorrectionScheme::v_GetSystemSingularChecks()
    {
        int vVar = m_session->GetVariables().size();
        Array<OneD, bool> vChecks(vVar, false);
        vChecks[vVar-1] = true;
        return vChecks;
    }

    /**
     *
     */
    int VelocityCorrectionScheme::v_GetForceDimension()
    {
        return m_session->GetVariables().size() - 1;
    }

    /**
     * Explicit part of the method - Advection, Forcing + HOPBCs
     */
    void VelocityCorrectionScheme::v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        EvaluateAdvectionTerms(inarray, outarray);

        // Smooth advection
        if(m_SmoothAdvection)
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_pressure->SmoothField(outarray[i]);
            }
        }

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }

        // Calculate High-Order pressure boundary conditions
        m_extrapolation->EvaluatePressureBCs(inarray,outarray,m_kinvis);
    }

    /**
     * Implicit part of the method - Poisson + nConv*Helmholtz
     */
    void VelocityCorrectionScheme::SolveUnsteadyStokesSystem(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time,
        const NekDouble aii_Dt)
    {
        // Set up flowrate if we're starting for the first time or the value of
        // aii_Dt has changed.
        if (m_flowrate > 0.0 && (aii_Dt != m_flowrateAiidt))
        {
            SetupFlowrate(aii_Dt);
        }

        int physTot = m_fields[0]->GetTotPoints();

        // Substep the pressure boundary condition if using substepping
        m_extrapolation->SubStepSetPressureBCs(inarray,aii_Dt,m_kinvis);

        // Set up forcing term for pressure Poisson equation
        SetUpPressureForcing(inarray, m_F, aii_Dt);

        // Solve Pressure System
        SolvePressure (m_F[0]);

        // Set up forcing term for Helmholtz problems
        SetUpViscousForcing(inarray, m_F, aii_Dt);

        // Solve velocity system
        SolveViscous( m_F, outarray, aii_Dt);

        // Apply flowrate correction
        if (m_flowrate > 0.0 && m_greenFlux != numeric_limits<NekDouble>::max())
        {
            NekDouble currentFlux = MeasureFlowrate(outarray);
            m_alpha = (m_flowrate - currentFlux) / m_greenFlux;

            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Svtvp(physTot, m_alpha, m_flowrateStokes[i], 1,
                             outarray[i], 1, outarray[i], 1);
                //Enusre coeff space is updated for next time step
                m_fields[i]->FwdTrans_IterPerExp(outarray[i],
                                                 m_fields[i]->UpdateCoeffs());
                // Impsoe symmetry of flow on coeff space (good to enfore periodicity). 
                m_fields[i]->LocalToGlobal();
                m_fields[i]->GlobalToLocal();
            }
        }
    }

    /**
     * Forcing term for Poisson solver solver
     */
    void   VelocityCorrectionScheme::v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt)
    {
        int i;
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_velocity.size();

        m_fields[0]->PhysDeriv(eX,fields[0], Forcing[0]);

        for(i = 1; i < nvel; ++i)
        {
            // Use Forcing[1] as storage since it is not needed for the pressure
            m_fields[i]->PhysDeriv(DirCartesianMap[i],fields[i],Forcing[1]);
            Vmath::Vadd(physTot,Forcing[1],1,Forcing[0],1,Forcing[0],1);
        }

        Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);
    }

    /**
     * Forcing term for Helmholtz solver
     */
    void   VelocityCorrectionScheme::v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int phystot = m_fields[0]->GetTotPoints();

        // Grad p
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());

        int nvel = m_velocity.size();
        if(nvel == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(),
                                  Forcing[m_velocity[0]],
                                  Forcing[m_velocity[1]]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(),
                                  Forcing[m_velocity[0]],
                                  Forcing[m_velocity[1]],
                                  Forcing[m_velocity[2]]);
        }

        // zero convective fields.
        for(int i = nvel; i < m_nConvectiveFields; ++i)
        {
            Vmath::Zero(phystot,Forcing[i],1);
        }

        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            Blas::Daxpy(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            Blas::Dscal(phystot,1.0/m_diffCoeff[i],&(Forcing[i])[0],1);
        }
    }


    /**
     * Solve pressure system
     */
    void   VelocityCorrectionScheme::v_SolvePressure(
        const Array<OneD, NekDouble>  &Forcing)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficient for equation
        factors[StdRegions::eFactorLambda] = 0.0;

        // Solver Pressure Poisson Equation
        m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(), factors);

        // Add presure to outflow bc if using convective like BCs
        m_extrapolation->AddPressureToOutflowBCs(m_kinvis);
    }

    /**
     * Solve velocity system
     */
    void   VelocityCorrectionScheme::v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt)
    {
        StdRegions::ConstFactorMap factors;
        StdRegions::VarCoeffMap varCoeffMap = StdRegions::NullVarCoeffMap;
        MultiRegions::VarFactorsMap varFactorsMap =
            MultiRegions::NullVarFactorsMap;

        AppendSVVFactors(factors,varFactorsMap);

        // Solve Helmholtz system and put in Physical space
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            // Setup coefficients for equation
            factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_diffCoeff[i];
            m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(),
                                   factors, varCoeffMap,
                                   varFactorsMap);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }

    void  VelocityCorrectionScheme::SetUpSVV(void)
    {

        m_session->MatchSolverInfo("SpectralVanishingViscosity",
                                   "PowerKernel", m_useSpecVanVisc, false);

        if(m_useSpecVanVisc)
        {
            m_useHomo1DSpecVanVisc = true;
        }
        else
        {
            m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                       "PowerKernel", m_useSpecVanVisc, false);
        }

        if(m_useSpecVanVisc)
        {
            m_IsSVVPowerKernel = true;
        }
        else
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosity","DGKernel",
                                       m_useSpecVanVisc, false);
            if(m_useSpecVanVisc)
            {
                m_useHomo1DSpecVanVisc = true;
            }
            else
            {
                m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                           "DGKernel", m_useSpecVanVisc, false);
            }

            if(m_useSpecVanVisc)
            {
                m_IsSVVPowerKernel = false;
            }
        }

        //set up varcoeff kernel if PowerKernel or DG is specified
        if(m_useSpecVanVisc)
        {
            Array<OneD, Array<OneD, NekDouble> > SVVVelFields = NullNekDoubleArrayofArray;
            if(m_session->DefinesFunction("SVVVelocityMagnitude"))
            {
                if (m_comm->GetRank() == 0)
                {
                    cout << "Seting up SVV velocity from "
                        "SVVVelocityMagnitude section in session file" << endl;
                }
                int nvel = m_velocity.size();
                int phystot = m_fields[0]->GetTotPoints();
                SVVVelFields = Array<OneD, Array<OneD, NekDouble> >(nvel);
                vector<string> vars;
                for(int i = 0; i < nvel; ++i)
                {
                    SVVVelFields[i] = Array<OneD, NekDouble>(phystot);
                    vars.push_back(m_session->GetVariable(m_velocity[i]));
                }

                // Load up files into  m_fields;
                GetFunction("SVVVelocityMagnitude")
                    ->Evaluate(vars,SVVVelFields);
            }

            m_svvVarDiffCoeff = Array<OneD, NekDouble>(m_fields[0]->GetNumElmts());
            SVVVarDiffCoeff(1.0,m_svvVarDiffCoeff,SVVVelFields);
            m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  1.0);
        }
        else
        {
            m_svvVarDiffCoeff = NullNekDouble1DArray;
            m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  0.1);
        }

        // Load parameters for Spectral Vanishing Viscosity
        if(m_useSpecVanVisc == false)
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosity","True",
                                       m_useSpecVanVisc, false);
            if(m_useSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosity","ExpKernel",
                                           m_useSpecVanVisc, false);
            }
            m_useHomo1DSpecVanVisc = m_useSpecVanVisc;

            if(m_useSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP","True",
                                           m_useSpecVanVisc, false);
                if(m_useSpecVanVisc == false)
                {
                    m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP","ExpKernel",
                                               m_useSpecVanVisc, false);
                }
            }
        }


        // Case of only Homo1D kernel
        if(m_session->DefinesSolverInfo("SpectralVanishingViscosityHomo1D"))
        {
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "True", m_useHomo1DSpecVanVisc, false);
            if(m_useHomo1DSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                       "ExpKernel", m_useHomo1DSpecVanVisc, false);
            }
        }

        m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
        m_session->LoadParameter("SVVCutoffRatioHomo1D",m_sVVCutoffRatioHomo1D,m_sVVCutoffRatio);
        m_session->LoadParameter("SVVDiffCoeffHomo1D",  m_sVVDiffCoeffHomo1D,  m_sVVDiffCoeff);

        if(m_HomogeneousType == eHomogeneous1D)
        {
            ASSERTL0(m_nConvectiveFields > 2,
                "Expect to have three velocity fields with homogenous expansion");

            if(m_useHomo1DSpecVanVisc)
            {
                Array<OneD, unsigned int> planes;
                planes = m_fields[0]->GetZIDs();

                int num_planes = planes.size();
                Array<OneD, NekDouble> SVV(num_planes,0.0);
                NekDouble fac;
                int kmodes = m_fields[0]->GetHomogeneousBasis()->GetNumModes();
                int pstart;

                pstart = m_sVVCutoffRatioHomo1D*kmodes;

                for(int n = 0; n < num_planes; ++n)
                {
                    if(planes[n] > pstart)
                    {
                        fac = (NekDouble)((planes[n] - kmodes)*(planes[n] - kmodes))/
                            ((NekDouble)((planes[n] - pstart)*(planes[n] - pstart)));
                        SVV[n] = m_sVVDiffCoeffHomo1D*exp(-fac)/m_kinvis;
                    }
                }

                for(int i = 0; i < m_velocity.size(); ++i)
                {
                    m_fields[m_velocity[i]]->SetHomo1DSpecVanVisc(SVV);
                }
            }
        }

    }

    void VelocityCorrectionScheme::SVVVarDiffCoeff(
                     const NekDouble velmag,
                     Array<OneD, NekDouble> &diffcoeff,
                     const  Array<OneD, Array<OneD, NekDouble> >  &vel)
    {
        int phystot = m_fields[0]->GetTotPoints();
        int nel = m_fields[0]->GetNumElmts();
        int nvel,cnt;

        Array<OneD, NekDouble> tmp;

        Vmath::Fill(nel,velmag,diffcoeff,1);

        if(vel != NullNekDoubleArrayofArray)
        {
            Array<OneD, NekDouble> Velmag(phystot);
            nvel = vel.size();
            // calculate magnitude of v
            Vmath::Vmul(phystot,vel[0],1,vel[0],1,Velmag,1);
            for(int n = 1; n < nvel; ++n)
            {
                Vmath::Vvtvp(phystot,vel[n],1,vel[n],1,Velmag,1,
                             Velmag,1);
            }
            Vmath::Vsqrt(phystot,Velmag,1,Velmag,1);


            cnt = 0;
            Array<OneD, NekDouble> tmp;
            // calculate mean value of vel mag.
            for(int i = 0; i < nel; ++i)
            {
                int nq = m_fields[0]->GetExp(i)->GetTotPoints();
                tmp = Velmag + cnt;
                diffcoeff[i] = m_fields[0]->GetExp(i)->Integral(tmp);
                Vmath::Fill(nq,1.0,tmp,1);
                NekDouble area = m_fields[0]->GetExp(i)->Integral(tmp);
                diffcoeff[i] = diffcoeff[i]/area;
                cnt += nq;
            }
        }
        else
        {
            nvel = m_expdim;
        }

        if(m_expdim == 3)
        {
            LocalRegions::Expansion3DSharedPtr exp3D;
            for (int e = 0; e < nel; e++)
            {
                exp3D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                NekDouble h = 0;
                for(int i = 0; i < exp3D->GetNedges(); ++i)
                {

                    h = max(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->dist(
                             *(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }

                int p = 0;
                for(int i = 0; i < 3; ++i)
                {
                    p = max(p,exp3D->GetBasisNumModes(i)-1);
                }

                diffcoeff[e] *= h/p;
            }
        }
        else
        {
            LocalRegions::Expansion2DSharedPtr exp2D;
            for (int e = 0; e < nel; e++)
            {
                exp2D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                NekDouble h = 0;
                for(int i = 0; i < exp2D->GetNedges(); ++i)
                {

                   h = max(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->dist(
                             *(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }

                int p = 0;
                for(int i = 0; i < 2; ++i)
                {
                    p = max(p,exp2D->GetBasisNumModes(i)-1);
                }

                diffcoeff[e] *= h/p;
            }
        }
    }

    void VelocityCorrectionScheme::AppendSVVFactors(
                                 StdRegions::ConstFactorMap &factors,
                                 MultiRegions::VarFactorsMap &varFactorsMap)
    {

        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_kinvis;
            if(m_svvVarDiffCoeff != NullNekDouble1DArray)
            {
                if(m_IsSVVPowerKernel)
                {
                    varFactorsMap[StdRegions::eFactorSVVPowerKerDiffCoeff] =
                        m_svvVarDiffCoeff;
                }
                else
                {
                    varFactorsMap[StdRegions::eFactorSVVDGKerDiffCoeff] =
                        m_svvVarDiffCoeff;
                }
            }
        }

    }
} //end of namespace
