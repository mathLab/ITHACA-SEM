///////////////////////////////////////////////////////////////////////////////
//
// File: SubSteppingExtrapolate.cpp
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
// Description: Abstract base class for SubSteppingExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/SubSteppingExtrapolate.h>

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/TimeIntegration/TimeIntegrationSolution.h>

using namespace std;

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string SubSteppingExtrapolate::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "SubStepping",
        SubSteppingExtrapolate::create,
        "SubStepping");

    SubSteppingExtrapolate::SubSteppingExtrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr  pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : Extrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
        m_session->LoadParameter("SubStepCFL", m_cflSafetyFactor, 0.5);
        m_session->LoadParameter("MinSubSteps", m_minsubsteps,1);
        m_session->LoadParameter("MaxSubSteps", m_maxsubsteps,100);

        int dim = m_fields[0]->GetCoordim(0);
        m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(dim);

        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();
        for(int i = 0; i < dim; ++i)
        {
            m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
        }
        m_fields[0]->GetTrace()->GetNormals(m_traceNormals);

    }

    SubSteppingExtrapolate::~SubSteppingExtrapolate()
    {
    }

    void SubSteppingExtrapolate::v_EvaluatePressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N,   NekDouble kinvis)
    {
        ASSERTL0(false,"This method should not be called by Substepping routine");
    }

    void SubSteppingExtrapolate::v_SubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr & IntegrationScheme )
    {
        unsigned int order = IntegrationScheme->GetOrder();

        // Set to 1 for first step and it will then be increased in
        // time advance routines
        if( IntegrationScheme->GetName() == "BackwardEuler" ||
            (IntegrationScheme->GetName() == "BDFImplicit" &&
             (order == 1 || order == 2)) )
        {
            // Note RK first order SSP is just Forward Euler.
            std::string vSubStepIntScheme        = "RungeKutta";
            std::string vSubStepIntSchemeVariant = "SSP";
            int         vSubStepIntSchemeOrder   = order;

            if( m_session->DefinesSolverInfo( "SubStepIntScheme" ) )
            {
                vSubStepIntScheme =
		  m_session->GetSolverInfo( "SubStepIntScheme" );
                vSubStepIntSchemeVariant = "";
                vSubStepIntSchemeOrder = order;
            }

            m_subStepIntegrationScheme =
                LibUtilities::GetTimeIntegrationSchemeFactory().CreateInstance(
                    vSubStepIntScheme,
                    vSubStepIntSchemeVariant,
                    vSubStepIntSchemeOrder,
		    std::vector<NekDouble>() );

            int nvel = m_velocity.size();
            int ndim = order+1;

            // Fields for linear/quadratic interpolation
            m_previousVelFields = Array<OneD, Array<OneD, NekDouble> >(ndim*nvel);
            int ntotpts  = m_fields[0]->GetTotPoints();
            m_previousVelFields[0] = Array<OneD, NekDouble>(ndim*nvel*ntotpts);

            for( int i = 1; i < ndim*nvel; ++i )
            {
                m_previousVelFields[i] = m_previousVelFields[i-1] + ntotpts;
            }
        }
        else
        {
            ASSERTL0(0,"Integration method not suitable: Options include BackwardEuler or BDFImplicitOrder{1,2}");
        }

        m_intSteps = IntegrationScheme->GetNumIntegrationPhases();

        // set explicit time-integration class operators
        m_subStepIntegrationOps.DefineOdeRhs(&SubSteppingExtrapolate::SubStepAdvection, this);
        m_subStepIntegrationOps.DefineProjection(&SubSteppingExtrapolate::SubStepProjection, this);
    }

    /**
     * Explicit Advection terms used by SubStepAdvance time integration
     */
    void SubSteppingExtrapolate::SubStepAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD,       NekDouble> > &outarray,
        const NekDouble time)
    {
        int i;
        int nVariables     = inarray.size();
        int nQuadraturePts = inarray[0].size();

        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs();

        /// Define an auxiliary variable to compute the RHS
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for(i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }

        Array<OneD, Array<OneD, NekDouble> > Velfields(m_velocity.size());

        Velfields[0] = Array<OneD, NekDouble> (nQuadraturePts*m_velocity.size());

        for(i = 1; i < m_velocity.size(); ++i)
        {
            Velfields[i] = Velfields[i-1] + nQuadraturePts;
        }

        SubStepExtrapolateField(fmod(time,m_timestep), Velfields);

        m_advObject->Advect(m_velocity.size(), m_fields, Velfields, inarray, outarray, time);

        for(i = 0; i < nVariables; ++i)
        {
            m_fields[i]->IProductWRTBase(outarray[i],WeakAdv[i]);
            // negation requried due to sign of DoAdvection term to be consistent
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);
        }

        AddAdvectionPenaltyFlux(Velfields, inarray, WeakAdv);

        /// Operations to compute the RHS
        for(i = 0; i < nVariables; ++i)
        {
            // Negate the RHS
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);

            /// Multiply the flux by the inverse of the mass matrix
            m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);

            /// Store in outarray the physical values of the RHS
            m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
        }
    }

    /**
     * Projection used by SubStepAdvance time integration
     */
    void SubSteppingExtrapolate::SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time)
    {
        ASSERTL1(inarray.size() == outarray.size(),"Inarray and outarray of different sizes ");

        for(int i = 0; i < inarray.size(); ++i)
        {
            Vmath::Vcopy(inarray[i].size(),inarray[i],1,outarray[i],1);
        }
    }


    /**
     *
     */
    void SubSteppingExtrapolate::v_SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        const NekDouble Aii_Dt,
        NekDouble kinvis)
    {
        //int nConvectiveFields =m_fields.size()-1;
        Array<OneD, Array<OneD, NekDouble> > nullvelfields;

        m_pressureCalls++;

        // Calculate non-linear and viscous BCs at current level and
        // put in m_pressureHBCs[0]
        CalcNeumannPressureBCs(inarray,nullvelfields,kinvis);

        // Extrapolate to m_pressureHBCs to n+1
        ExtrapolateArray(m_pressureHBCs);

        // Add (phi,Du/Dt) term to m_presureHBC
        AddDuDt();

        // Copy m_pressureHBCs to m_PbndExp
        CopyPressureHBCsToPbndExp();

        // Evaluate High order outflow conditiosn if required.
        CalcOutflowBCs(inarray, kinvis);
    }


    /**
     *
     */
    void SubSteppingExtrapolate::v_SubStepSaveFields(const int nstep)
    {
        int i,n;
        int nvel = m_velocity.size();
        int npts = m_fields[0]->GetTotPoints();

        // rotate fields
        int nblocks = m_previousVelFields.size()/nvel;
        Array<OneD, NekDouble> save;

        // rotate storage space
        for(n = 0; n < nvel; ++n)
        {
            save = m_previousVelFields[(nblocks-1)*nvel+n];

            for(i = nblocks-1; i > 0; --i)
            {
                m_previousVelFields[i*nvel+n] = m_previousVelFields[(i-1)*nvel+n];
            }

            m_previousVelFields[n] = save;
        }

        // Put previous field
        for(i = 0; i < nvel; ++i)
        {
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),
                                              m_fields[m_velocity[i]]->UpdatePhys());
            Vmath::Vcopy(npts,m_fields[m_velocity[i]]->GetPhys(),1,
                         m_previousVelFields[i],1);
        }

        if(nstep == 0)// initialise all levels with first field
        {
            for(n = 0; n < nvel; ++n)
            {
                for(i = 1; i < nblocks; ++i)
                {
                    Vmath::Vcopy(npts,m_fields[m_velocity[n]]->GetPhys(),1,
                                 m_previousVelFields[i*nvel+n],1);

                }
            }
        }
    }

    /**
     *
     */
    void SubSteppingExtrapolate::v_SubStepAdvance(
        const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr &integrationSoln,
        int nstep,
        NekDouble time )
    {
        int n;
        int nsubsteps;

        NekDouble dt;

        Array<OneD, Array<OneD, NekDouble> > fields;

        static int ncalls = 1;
        int  nint         = min(ncalls++, m_intSteps);

        //this needs to change
        m_comm = m_fields[0]->GetComm()->GetRowComm();

        // Get the proper time step with CFL control
        dt = GetSubstepTimeStep();

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1;
        nsubsteps = max(m_minsubsteps, nsubsteps);

        ASSERTL0(nsubsteps < m_maxsubsteps,"Number of substeps has exceeded maximum");

        dt = m_timestep/nsubsteps;

        if (m_infosteps && !((nstep+1)%m_infosteps) && m_comm->GetRank() == 0)
        {
            cout << "Sub-integrating using "<< nsubsteps
                 << " steps over Dt = "     << m_timestep
                 << " (SubStep CFL="        << m_cflSafetyFactor << ")"<< endl;
        }

        for (int m = 0; m < nint; ++m)
        {
            // We need to update the fields held by the m_integrationSoln
            fields = integrationSoln->UpdateSolutionVector()[m];

            // Initialise NS solver which is set up to use a GLM method
            // with calls to EvaluateAdvection_SetPressureBCs and
            // SolveUnsteadyStokesSystem
            LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr
                SubIntegrationSoln = m_subStepIntegrationScheme->InitializeScheme( dt, fields, time, m_subStepIntegrationOps );

            for(n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(n, dt, SubIntegrationSoln,
                                                                   m_subStepIntegrationOps);
            }

            // set up HBC m_acceleration field for Pressure BCs
            IProductNormVelocityOnHBC(fields,m_iprodnormvel[m]);

            // Reset time integrated solution in m_integrationSoln
            integrationSoln->SetSolVector(m,fields);
        }
    }


    /**
     *
     */
    NekDouble SubSteppingExtrapolate::GetSubstepTimeStep()
    {
        int n_element      = m_fields[0]->GetExpSize();

        const Array<OneD, int> ExpOrder=m_fields[0]->EvalBasisNumModesMaxPerExp();

        const NekDouble cLambda = 0.2; // Spencer book pag. 317

        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields(m_velocity.size());

        for(int i = 0; i < m_velocity.size(); ++i)
        {
            velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
        }
        stdVelocity = GetMaxStdVelocity(velfields);

        for(int el = 0; el < n_element; ++el)
        {
            tstep[el] = m_cflSafetyFactor /
                (stdVelocity[el] * cLambda *
                 (ExpOrder[el]-1) * (ExpOrder[el]-1));
        }

        NekDouble TimeStep = Vmath::Vmin(n_element, tstep, 1);
        m_comm->AllReduce(TimeStep,LibUtilities::ReduceMin);

        return TimeStep;
    }




    void SubSteppingExtrapolate::AddAdvectionPenaltyFlux(
                                                         const Array<OneD, const Array<OneD, NekDouble> > &velfield,
                                                         const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                                                         Array<OneD, Array<OneD, NekDouble> > &Outarray)
    {
        ASSERTL1(
                 physfield.size() == Outarray.size(),
                 "Physfield and outarray are of different dimensions");

        int i;

        /// Number of trace points
        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();

        /// Number of spatial dimensions
        int nDimensions = m_bnd_dim;

        /// Forward state array
        Array<OneD, NekDouble> Fwd(3*nTracePts);

        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

        /// upwind numerical flux state array
        Array<OneD, NekDouble> numflux = Bwd + nTracePts;

        /// Normal velocity array
        Array<OneD, NekDouble> Vn (nTracePts, 0.0);

        // Extract velocity field along the trace space and multiply by trace normals
        for(i = 0; i < nDimensions; ++i)
        {
            m_fields[0]->ExtractTracePhys(velfield[i], Fwd);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        }

        for(i = 0; i < physfield.size(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            /// Note: Needs to have correct i value to get boundary conditions
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            /// Upwind between elements
            m_fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux);

            /// Construct difference between numflux and Fwd,Bwd
            Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
            Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

            /// Calculate the numerical fluxes multipling Fwd, Bwd and
            /// numflux by the normal advection velocity
            Vmath::Vmul(nTracePts, Fwd, 1, Vn, 1, Fwd, 1);
            Vmath::Vmul(nTracePts, Bwd, 1, Vn, 1, Bwd, 1);

            m_fields[0]->AddFwdBwdTraceIntegral(Fwd,Bwd,Outarray[i]);
        }
    }

    /**
     * Extrapolate field using equally time spaced field un,un-1,un-2, (at
     * dt intervals) to time n+t at order Ord
     */
    void SubSteppingExtrapolate::SubStepExtrapolateField(
                                                         NekDouble toff,
                                                         Array< OneD, Array<OneD, NekDouble> > &ExtVel)
    {
        int npts = m_fields[0]->GetTotPoints();
        int nvel = m_velocity.size();
        int i,j;
        Array<OneD, NekDouble> l(4);

        int ord = m_intSteps;

        // calculate Lagrange interpolants
        Vmath::Fill(4,1.0,l,1);

        for(i = 0; i <= ord; ++i)
        {
            for(j = 0; j <= ord; ++j)
            {
                if(i != j)
                {
                    l[i] *= (j*m_timestep+toff);
                    l[i] /= (j*m_timestep-i*m_timestep);
                }
            }
        }

        for(i = 0; i < nvel; ++i)
        {
            Vmath::Smul(npts,l[0],m_previousVelFields[i],1,ExtVel[i],1);

            for(j = 1; j <= ord; ++j)
            {
                Blas::Daxpy(npts,l[j],m_previousVelFields[j*nvel+i],1,
                            ExtVel[i],1);
            }
        }
    }

    /**
     *
     */
    void SubSteppingExtrapolate::v_MountHOPBCs(
                                               int HBCdata,
                                               NekDouble kinvis,
                                               Array<OneD, NekDouble> &Q,
                                               Array<OneD, const NekDouble> &Advection)
    {
        Vmath::Smul(HBCdata,-kinvis,Q,1,Q,1);
    }

    std::string SubSteppingExtrapolate::v_GetSubStepName(void)
    {
        return m_subStepIntegrationScheme->GetName();
    }

}
