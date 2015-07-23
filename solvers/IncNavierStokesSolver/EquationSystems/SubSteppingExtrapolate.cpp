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
// Description: Abstract base class for SubSteppingExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/SubSteppingExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

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
        m_session->LoadParameter("CFL", m_cflSafetyFactor, 0.5);
        m_session->LoadParameter("MinSubSteps", minsubsteps,0);
    }

    SubSteppingExtrapolate::~SubSteppingExtrapolate()
    {
    }
    
    void SubSteppingExtrapolate::v_SubSteppingTimeIntegration(
        int intMethod,
        const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)
    {
        int i;
        
        // Set to 1 for first step and it will then be increased in
        // time advance routines
        switch(intMethod)
        {
            case LibUtilities::eBackwardEuler:
            case LibUtilities::eBDFImplicitOrder1: 
            {
                m_subStepIntegrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance("ForwardEuler");
                
                // Fields for linear interpolation
                m_previousVelFields = Array<OneD, Array<OneD, NekDouble> >(2*m_fields.num_elements());                    
                int ntotpts  = m_fields[0]->GetTotPoints();
                m_previousVelFields[0] = Array<OneD, NekDouble>(2*m_fields.num_elements()*ntotpts);
                for(i = 1; i < 2*m_fields.num_elements(); ++i)
                {
                    m_previousVelFields[i] = m_previousVelFields[i-1] + ntotpts; 
                }
                
            }
            break;
            case LibUtilities::eBDFImplicitOrder2:
            {
                m_subStepIntegrationScheme = LibUtilities::GetTimeIntegrationWrapperFactory().CreateInstance("RungeKutta2_ImprovedEuler");
                    
                int nvel = m_velocity.num_elements();
                
                // Fields for quadratic interpolation
                m_previousVelFields = Array<OneD, Array<OneD, NekDouble> >(3*nvel);
                
                int ntotpts  = m_fields[0]->GetTotPoints();
                m_previousVelFields[0] = Array<OneD, NekDouble>(3*nvel*ntotpts);
                for(i = 1; i < 3*nvel; ++i)
                {
                    m_previousVelFields[i] = m_previousVelFields[i-1] + ntotpts; 
                }
                 
                }
            break;
            default:
                ASSERTL0(0,"Integration method not suitable: Options include BackwardEuler or BDFImplicitOrder1");
                break;
        }
	
        m_intSteps = IntegrationScheme->GetIntegrationSteps();

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
        int nVariables     = inarray.num_elements();
        int nQuadraturePts = inarray[0].num_elements();
        
        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        /// Define an auxiliary variable to compute the RHS 
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for(i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }
        
        Array<OneD, Array<OneD, NekDouble> > Velfields(m_velocity.num_elements());
        Array<OneD, int> VelIds(m_velocity.num_elements());
        
        Velfields[0] = Array<OneD, NekDouble> (nQuadraturePts*m_velocity.num_elements());
        
        for(i = 1; i < m_velocity.num_elements(); ++i)
        {
            Velfields[i] = Velfields[i-1] + nQuadraturePts; 
        }

        SubStepExtrapolateField(fmod(time,m_timestep), Velfields);
        
        m_advObject->Advect(m_velocity.num_elements(), m_fields, Velfields, inarray, outarray, time);
        
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
        
        //add the force
        /*if(m_session->DefinesFunction("BodyForce"))
        {
            if(m_SingleMode || m_HalfMode)
            {
                for(int i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_forces[i]->SetWaveSpace(true);    
                    m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),
                                          m_forces[i]->UpdatePhys());
                }
            }
            
            int nqtot = m_fields[0]->GetTotPoints();
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,outarray[i],1);
            }
            }*/
    }
        
    /** 
     * Projection used by SubStepAdvance time integration
     */
    void SubSteppingExtrapolate::SubStepProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,  
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        ASSERTL1(inarray.num_elements() == outarray.num_elements(),"Inarray and outarray of different sizes ");

        for(int i = 0; i < inarray.num_elements(); ++i)
        {
            Vmath::Vcopy(inarray[i].num_elements(),inarray[i],1,outarray[i],1);
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
        int nConvectiveFields =m_fields.num_elements()-1;
        Array<OneD, Array<OneD, NekDouble> > velfields(nConvectiveFields);
        
        for(int i = 0; i < nConvectiveFields; ++i)
        {
            velfields[i] = m_fields[m_velocity[i]]->GetPhys(); 
        }

        EvaluatePressureBCs(velfields,inarray,kinvis);
		
        AddDuDt(inarray,Aii_Dt);
    }

    /** 
     * 
     */
    void SubSteppingExtrapolate::v_SubStepSaveFields(
        const int nstep)
    {
        int i,n;
        int nvel = m_velocity.num_elements();
        int npts = m_fields[0]->GetTotPoints();

        // rotate fields 
        int nblocks = m_previousVelFields.num_elements()/nvel; 
        Array<OneD, NekDouble> save; 
        
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
        const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
        int nstep, 
        NekDouble time)
    {
        int n;
        int nsubsteps;
        
        NekDouble dt; 
        
        Array<OneD, Array<OneD, NekDouble> > fields, velfields;
        
        static int ncalls = 1;
        int  nint         = min(ncalls++, m_intSteps);
        
        Array<OneD, NekDouble> CFL(m_fields[0]->GetExpSize(), 
                                   m_cflSafetyFactor);
        //this needs to change
        m_comm = m_fields[0]->GetComm()->GetRowComm();

        // Get the proper time step with CFL control
        dt = GetSubstepTimeStep();

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1; 
        nsubsteps = max(minsubsteps, nsubsteps);

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
            LibUtilities::TimeIntegrationSolutionSharedPtr 
                SubIntegrationSoln = m_subStepIntegrationScheme->
                InitializeScheme(dt, fields, time, m_subStepIntegrationOps);
            
            for(n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(
                    n, dt, SubIntegrationSoln,
                    m_subStepIntegrationOps);
            }
            
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
        Array<OneD, int> ExpOrderList (n_element, ExpOrder);
        
        const NekDouble cLambda = 0.2; // Spencer book pag. 317
        
        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields(m_velocity.num_elements());
        
        for(int i = 0; i < m_velocity.num_elements(); ++i)
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
            physfield.num_elements() == Outarray.num_elements(),
            "Physfield and outarray are of different dimensions");
        
        int i;
        
        /// Number of trace points
        int nTracePts   = m_fields[0]->GetTrace()->GetNpoints();

        /// Number of spatial dimensions
        int nDimensions = m_curl_dim;

        Array<OneD, Array<OneD, NekDouble> > traceNormals(m_curl_dim);

        for(i = 0; i < nDimensions; ++i)
        {
            traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
        }

        //trace normals
        m_fields[0]->GetTrace()->GetNormals(traceNormals);

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
            m_fields[0]->ExtractTracePhys(m_fields[m_velocity[i]]->GetPhys(), Fwd);
            Vmath::Vvtvp(nTracePts, traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        }
        
        for(i = 0; i < physfield.num_elements(); ++i)
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
        int nvel = m_velocity.num_elements();
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
	
    /** 
     * 
     */
    void SubSteppingExtrapolate::AddDuDt(
        const Array<OneD, const Array<OneD, NekDouble> >  &N, 
        NekDouble Aii_Dt)
    {
        switch(m_velocity.num_elements())
        {
            case 1:
                ASSERTL0(
                    false,
                    "Velocity correction scheme not designed to have just one velocity component");
                break;
            case 2:
                AddDuDt2D(N,Aii_Dt);
                break;
            case 3:
                AddDuDt3D(N,Aii_Dt);
                break;
        }
    }
        
    /** 
     * 
     */    
    void SubSteppingExtrapolate::AddDuDt2D(
        const Array<OneD, const Array<OneD, NekDouble> >  &N, 
        NekDouble Aii_Dt)
    {
        int i,n;
        ASSERTL0(m_velocity.num_elements() == 2," Routine currently only set up for 2D");
        
        int pindex=2;

        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp,UBndExp,VBndExp;
        
        PBndConds = m_fields[pindex]->GetBndConditions();
        PBndExp   = m_fields[pindex]->GetBndCondExpansions();
        
        UBndExp   = m_fields[m_velocity[0]]->GetBndCondExpansions();
        VBndExp   = m_fields[m_velocity[1]]->GetBndCondExpansions();
        
        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansion1DSharedPtr Pbc;
        
        
        int cnt,elmtid,nq,offset, boundary,ncoeffs;
        
        Array<OneD, NekDouble> N1(m_pressureBCsMaxPts), N2(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> ubc(m_pressureBCsMaxPts),vbc(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> Pvals(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> Nu,Nv,Ptmp;
        
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {            
            std::string type = PBndConds[n]->GetUserDefined(); 
            
            if(boost::iequals(type,"H"))
            {
                for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and edge of this expansion. 
                    elmtid = m_pressureBCtoElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    offset = m_fields[0]->GetPhys_Offset(elmtid);
                    
                    Nu = N[0] + offset;
                    Nv = N[1] + offset; 
                    
                    Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[n]->GetExp(i));
                    nq       = Pbc->GetTotPoints();
                    ncoeffs  = Pbc->GetNcoeffs();
                    boundary = m_pressureBCtoTraceID[cnt];
                    
                    // Get velocity bc
                    UBndExp[n]->GetExp(i)->BwdTrans(UBndExp[n]->GetCoeffs() + UBndExp[n]->GetCoeff_Offset(i),ubc);
                    VBndExp[n]->GetExp(i)->BwdTrans(VBndExp[n]->GetCoeffs() + VBndExp[n]->GetCoeff_Offset(i),vbc);
                    
                    
                    // Get edge values and put into Nu,Nv
                    elmt->GetEdgePhysVals(boundary,Pbc,Nu,N1);
                    elmt->GetEdgePhysVals(boundary,Pbc,Nv,N2);
                    
                        
                    // Take different as Forward Euler but N1,N2
                    // actually contain the integration of the
                    // previous steps from the time integration
                    // scheme.
                    Vmath::Vsub(nq,ubc,1,N1,1,ubc,1);
                    Vmath::Vsub(nq,vbc,1,N2,1,vbc,1);
                    
                        
                    // Divide by aii_Dt to get correct Du/Dt.  This is
                    // because all coefficients in the integration
                    // scheme are normalised so u^{n+1} has unit
                    // coefficient and N is already multiplied by
                    // local coefficient when taken from integration
                    // scheme
                    Blas::Dscal(nq,1.0/Aii_Dt,&ubc[0],1);
                    Blas::Dscal(nq,1.0/Aii_Dt,&vbc[0],1);
                    
                    // subtrace off du/dt derivative 
                    Pbc->NormVectorIProductWRTBase(ubc,vbc,Pvals); 
                    
                    Vmath::Vsub(ncoeffs,Ptmp = PBndExp[n]->UpdateCoeffs()
                                +PBndExp[n]->GetCoeff_Offset(i),1, Pvals,1, 
                                Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),1);
                }
            }
            else  // No High order
            {
                cnt += PBndExp[n]->GetExpSize();
            }
        }        
    }
    
    /** 
     * 
     */   
    void SubSteppingExtrapolate::AddDuDt3D(
        const Array<OneD, const Array<OneD, NekDouble> >  &N, 
        NekDouble Aii_Dt)
    {
        int i,n;
        ASSERTL0(m_velocity.num_elements() == 3," Routine currently only set up for 3D");
        int pindex=3;                

        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp,UBndExp,VBndExp,WBndExp;

        PBndConds = m_fields[pindex]->GetBndConditions();
        PBndExp   = m_fields[pindex]->GetBndCondExpansions();
        
        UBndExp   = m_fields[m_velocity[0]]->GetBndCondExpansions();
        VBndExp   = m_fields[m_velocity[1]]->GetBndCondExpansions();
        WBndExp   = m_fields[m_velocity[2]]->GetBndCondExpansions();
        
        StdRegions::StdExpansionSharedPtr  elmt;
        StdRegions::StdExpansion2DSharedPtr Pbc;
            
        int cnt,elmtid,nq,offset, boundary,ncoeffs;        
        
        Array<OneD, NekDouble> N1(m_pressureBCsMaxPts), N2(m_pressureBCsMaxPts), 
            N3(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> ubc(m_pressureBCsMaxPts),vbc(m_pressureBCsMaxPts),
            wbc(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> Pvals(m_pressureBCsMaxPts);
        Array<OneD, NekDouble> Nu,Nv,Nw,Ptmp;
        
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {            
            std::string type = PBndConds[n]->GetUserDefined(); 
            
            if(boost::iequals(type,"H"))
            {
                for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and face of this expansion. 
                    elmtid = m_pressureBCtoElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    offset = m_fields[0]->GetPhys_Offset(elmtid);
                    
                    Nu = N[0] + offset;
                    Nv = N[1] + offset;
                    Nw = N[2] + offset;
                    
                    Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (PBndExp[n]->GetExp(i));
                    nq       = Pbc->GetTotPoints();
                    ncoeffs  = Pbc->GetNcoeffs();
                    boundary = m_pressureBCtoTraceID[cnt];
                    
                    // Get velocity bc
                    UBndExp[n]->GetExp(i)->BwdTrans(UBndExp[n]->GetCoeffs() + 
                                                    UBndExp[n]->GetCoeff_Offset(i),ubc);
                    VBndExp[n]->GetExp(i)->BwdTrans(VBndExp[n]->GetCoeffs() + 
                                                    VBndExp[n]->GetCoeff_Offset(i),vbc);
                    WBndExp[n]->GetExp(i)->BwdTrans(WBndExp[n]->GetCoeffs() + 
                                                    WBndExp[n]->GetCoeff_Offset(i),wbc);
                    
                    // Get edge values and put into N1,N2,N3
                    elmt->GetFacePhysVals(boundary,Pbc,Nu,N1);
                    elmt->GetFacePhysVals(boundary,Pbc,Nv,N2);
                    elmt->GetFacePhysVals(boundary,Pbc,Nw,N3);
                    
                    
                    // Take different as Forward Euler but N1,N2,N3
                    // actually contain the integration of the
                    // previous steps from the time integration
                    // scheme.
                    Vmath::Vsub(nq,ubc,1,N1,1,ubc,1);
                    Vmath::Vsub(nq,vbc,1,N2,1,vbc,1);
                    Vmath::Vsub(nq,wbc,1,N3,1,wbc,1);
                    
                    // Divide by aii_Dt to get correct Du/Dt.  This is
                    // because all coefficients in the integration
                    // scheme are normalised so u^{n+1} has unit
                    // coefficient and N is already multiplied by
                    // local coefficient when taken from integration
                    // scheme
                    Blas::Dscal(nq,1.0/Aii_Dt,&ubc[0],1);
                    Blas::Dscal(nq,1.0/Aii_Dt,&vbc[0],1);
                    Blas::Dscal(nq,1.0/Aii_Dt,&wbc[0],1);
                    
                    // subtrace off du/dt derivative 
                    Pbc->NormVectorIProductWRTBase(ubc,vbc,wbc,Pvals); 
                    
                    Vmath::Vsub(
                        ncoeffs,
                        Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),
                        1,Pvals,1, 
                        Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),
                        1);
                }
            }
            else  // No High order
            {
                cnt += PBndExp[n]->GetExpSize();
            }
        }        
    }    
}

