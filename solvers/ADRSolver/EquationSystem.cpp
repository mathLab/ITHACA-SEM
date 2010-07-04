#include <string>
using std::string;

#include <ADRSolver/EquationSystem.h>

namespace Nektar
{
    EquationSystem::EquationSystem(SessionReaderSharedPtr& pSession,
            LibUtilities::TimeIntegrationSchemeOperators& pOde)
        : ADRBase(pSession->getFilename(), true),
          mSession(pSession),
          mTimeIntMethod(LibUtilities::eClassicalRungeKutta4),
          mExplicitAdvection(true),
          mExplicitDiffusion(true),
          mExplicitReaction(true),
          mEpsilon(1.0)
    {
        ZeroPhysFields();

        mSession->matchSolverInfo("ADVECTIONADVANCEMENT", "Explicit", mExplicitAdvection, true);
        mSession->matchSolverInfo("DIFFUSIONADVANCEMENT", "Explicit", mExplicitDiffusion, true);
        mSession->matchSolverInfo("REACTIONADVANCEMENT",  "Explicit", mExplicitReaction,  true);

        if (mSession->definesSolverInfo("TIMEINTEGRATIONMETHOD"))
        {
            int i;
            for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
            {
                bool match;
                mSession->matchSolverInfo("TIMEINTEGRATIONMETHOD", LibUtilities::TimeIntegrationMethodMap[i], match, false);
                if (match)
                {
                    mTimeIntMethod = (LibUtilities::TimeIntegrationMethod)i;
                    break;
                }
            }
            ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
                                            "Invalid time integration type.");
        }

    }

    EquationSystem::~EquationSystem()
    {

    }

    void EquationSystem::doInitialise()
    {
        mSession->loadParameter("epsilon", mEpsilon, 0.0);

        if (!v_isSteady())
        {
            mSession->loadParameter("IO_InfoSteps",    mInfoSteps,  0);

            // check that any user defined boundary condition is indeed
            // implemented
            for(int n = 0;
                    n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {
                // Time Dependent Boundary Condition (if no use defined then
                // this is empty)
                if (m_fields[0]->GetBndConditions()[n]
                            ->GetUserDefined().GetEquation() != "")
                {
                    if (m_fields[0]->GetBndConditions()[n]
                            ->GetUserDefined().GetEquation() != "TimeDependent")
                    {
                        ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary "
                                       "condition");
                    }
                }
            }

        }
    }

    void EquationSystem::printSummary(std::ostream &out)
    {
        out << "=======================================================================" << endl;
        out << "\tEquation Type   : " << mSession->getSolverInfo("EQTYPE") << endl;
        ADRBase::SessionSummary(out);
        out << "\tAdvection       : " << (mExplicitAdvection ? "explicit" : "implicit") << endl;
        out << "\tDiffusion       : " << (mExplicitDiffusion ? "explicit" : "implicit") << endl;
        out << "\tReaction        : " << (mExplicitReaction  ? "explicit" : "implicit") << endl;
        v_printSummary(out);

        if (!isSteady())
        {
            out << "\tIntegration     : " << LibUtilities::TimeIntegrationMethodMap[mTimeIntMethod] << endl;
            ADRBase::TimeParamSummary(out);
        }
        out << "=======================================================================" << endl;
    }

    void EquationSystem::doOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        v_doOdeRhs(inarray, outarray, time);
    }

    void EquationSystem::doImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        v_doImplicitSolve(inarray, outarray, time, lambda);
    }

    void EquationSystem::setPhysForcingFunction()
    {
        SetPhysForcingFunctions(m_fields);
    }

    void EquationSystem::doSolveHelmholtz()
    {
        v_doSolveHelmholtz();
    }

    void EquationSystem::GeneralTimeIntegration(int nsteps,
                              LibUtilities::TimeIntegrationMethod IntMethod,
                              LibUtilities::TimeIntegrationSchemeOperators ode)
    {
        int i,n,nchk = 0;
        int ncoeffs = m_fields[0]->GetNcoeffs();
        int nvariables = m_fields.num_elements();

        // Set up wrapper to fields data storage.
        Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
        Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->SetPhysState(false);
            fields[i]  = m_fields[i]->UpdateCoeffs();
        }

        // Declare an array of TimeIntegrationSchemes
        // For multi-stage methods, this array will have just one entry containing
        // the actual multi-stage method...
        // For multi-steps method, this can have multiple entries
        //  - the first scheme will used for the first timestep (this is an initialization scheme)
        //  - the second scheme will used for the second timestep (this is an initialization scheme)
        //  - ...
        //  - the last scheme will be used for all other time-steps (this will be the actual scheme)
        Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
        LibUtilities::TimeIntegrationSolutionSharedPtr u;
        int numMultiSteps;

        switch(IntMethod)
        {
            case LibUtilities::eIMEXdirk_2_3_2:
            case LibUtilities::eIMEXdirk_3_4_3:
            case LibUtilities::eDIRKOrder2:
            case LibUtilities::eDIRKOrder3:
            case LibUtilities::eBackwardEuler:
            case LibUtilities::eForwardEuler:
            case LibUtilities::eClassicalRungeKutta4:
            {
                numMultiSteps = 1;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                LibUtilities::TimeIntegrationSchemeKey IntKey(IntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

                u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,ode);
                break;
            }
            case LibUtilities::eAdamsBashforthOrder2:
            {
                numMultiSteps = 2;

                IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

                // Used in the first time step to initalize the scheme
                LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eForwardEuler);

                // Used for all other time steps
                LibUtilities::TimeIntegrationSchemeKey IntKey1(IntMethod);
                IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

                // Initialise the scheme for the actual time integration scheme
                u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,ode);
                break;
            }
            default:
            {
                ASSERTL0(false,"populate switch statement for integration scheme");
            }
        }

        std::string outname = m_sessionName + ".his";
        std::ofstream hisFile (outname.c_str());

        for(n = 0; n < nsteps; ++n)
        {
            //----------------------------------------------
            // Perform time step integration
            //----------------------------------------------
            if( n < numMultiSteps-1)
            {
                // Use initialisation schemes
                fields = IntScheme[n]->TimeIntegrate(m_timestep,u,ode);
            }
            else
            {
                fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,ode);
            }

            m_time += m_timestep;

            if(m_projectionType==eGalerkin)
            {
//                SetBoundaryConditions(m_time);
            }

            //----------------------------------------------
            // Dump analyser information
            //----------------------------------------------
            if(!((n+1)%mInfoSteps))
            {
                cout << "\rSteps: " << n+1 << "\t Time: " << m_time << "\t " << flush;

            }

            if(n&&(!((n+1)%m_checksteps)))
            {
                for(i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->UpdateCoeffs() = fields[i];
                }
                Checkpoint_Output(nchk++);
                WriteHistoryData(hisFile);
            }
        }


        for(i = 0; i < nvariables; ++i)
        {
            m_fields[i]->UpdateCoeffs() = fields[i];
        }
    }

    bool EquationSystem::isSteady()
    {
        return v_isSteady();
    }

    LibUtilities::TimeIntegrationMethod EquationSystem::getTimeIntMethod()
    {
        return mTimeIntMethod;
    }


    void EquationSystem::evaluateFunction(Array<OneD, NekDouble>& pArray,
            SpatialDomains::ConstUserDefinedEqnShPtr pEqn) {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        if (pArray.num_elements() != nq)
        {
            pArray = Array<OneD, NekDouble>(nq);
        }
        for(int i = 0; i < nq; i++)
        {
            pArray[i] = pEqn->Evaluate(x0[i],x1[i],x2[i]);
        }

    }

    void EquationSystem::setBoundaryConditions(NekDouble time)
    {
      int nvariables = m_fields.num_elements();
      for (int i = 0; i < nvariables; ++i)
      {
          m_fields[i]->EvaluateBoundaryConditions(time);
      }
    }


    // Virtual functions
    void EquationSystem::v_printSummary(std::ostream &out)
    {

    }

    void EquationSystem::v_doOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        ASSERTL0(false, "ODE RHS not defined for this equation.");
    }

    void EquationSystem::v_doImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables = inarray.num_elements();
        int ncoeffs    = inarray[0].num_elements();
        int nq = m_fields[0]->GetNpoints();

        //Array<OneD, Array<OneD, NekDouble> > Forcing(1);
        //Forcing[0] = Array<OneD, NekDouble>(ncoeffs,0.0);

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs

        //MultiRegions::GlobalMatrixKey key(StdRegions::eMass);

        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(ncoeffs, -1.0/lambda/mEpsilon, inarray[i], 1, outarray[i], 1);

            // Update coeffs to m_fields
            m_fields[i]->UpdateCoeffs() = outarray[i];

            // Backward Transformation to nodal coefficients
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), m_fields[i]->UpdatePhys());
            // m_fields[i]->SetPhysState(true);

            NekDouble kappa = 1.0/lambda/mEpsilon;

            // Solve a system of equations with Helmholtz solver
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                            m_fields[i]->UpdateCoeffs(),kappa);
            m_fields[i]->SetPhysState(false);

            // The solution is Y[i]
            outarray[i] = m_fields[i]->GetCoeffs();
        }
    }

    void EquationSystem::v_doSolveHelmholtz()
    {
        ASSERTL0(false, "SolveHelmholtz not defined for this equation.");
    }

    bool EquationSystem::v_isSteady()
    {
        return false;
    }

    void EquationSystem::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL0(false, "v_GetFluxVector not defined.");
    }

    // Evaulate flux = m_fields*ivel for i th component of Vu for direction j
    void EquationSystem::v_GetFluxVector(const int i, const int j, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int k = 0; k < flux.num_elements(); ++k)
        {
            Vmath::Zero(GetNpoints(),flux[k],1);
        }
        Vmath::Vcopy(GetNpoints(),physfield[i],1,flux[j],1);
    }

    void EquationSystem::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        ASSERTL0(false, "Not implemented.");
    }

    void EquationSystem::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield,
                                 Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                                 Array<OneD, Array<OneD, NekDouble> > &numfluxY )
    {
        ASSERTL0(false, "Not implemented");
    }


    void EquationSystem::v_NumFluxforScalar(Array<OneD, Array<OneD, NekDouble> > &ufield,
                                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &uflux)
    {
        int i,j;
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = uflux.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);
        Array<OneD, NekDouble > fluxtemp (nTraceNumPoints,0.0);

        // Get the sign of (v \cdot n), v = an arbitrary vector

        //  Evaulate upwind flux of uflux = \hat{u} \phi \cdot u = u^{(+,-)} n
        for (j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables ; ++i)
            {
                //  Compute Forward and Backward value of ufield of i direction

                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uFwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uFwd

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uBwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uBwd

                m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],Fwd,Bwd,fluxtemp);

                // Imposing weak boundary condition with flux
                // if Vn >= 0, uflux = uBwd at Neumann, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick uflux = uBwd
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick uflux = uBwd

                // if Vn >= 0, uflux = uFwd at Neumann, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick uflux = uFwd
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick uflux = uFwd

                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforScalar(i,ufield[i],fluxtemp);
                }

                // if Vn >= 0, flux = uFwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uFwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uFwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                // else if Vn < 0, flux = uBwd*(tan_{\xi}^- \cdot \vec{n} ), i.e,
                // edge::eForward, uBwd \(\tan_{\xi}^Fwd \cdot \vec{n} )
                // edge::eBackward, uBwd \(\tan_{\xi}^Bwd \cdot \vec{n} )

                Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,fluxtemp,1,uflux[j][i],1);

            }
        }
    }

    void EquationSystem::v_NumFluxforVector(Array<OneD, Array<OneD, NekDouble> > &ufield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> >  &qflux)
    {
        int nTraceNumPoints = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = qfield.num_elements();

        NekDouble C11 = 1.0;
        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);

        Array<OneD, NekDouble > qFwd(nTraceNumPoints);
        Array<OneD, NekDouble > qBwd(nTraceNumPoints);
        Array<OneD, NekDouble > qfluxtemp(nTraceNumPoints,0.0);

        Array<OneD, NekDouble > uterm(nTraceNumPoints);

        // Evaulate upwind flux of qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
        for(int i = 0; i < nvariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTraceNumPoints,0.0);
            for(int j = 0; j < nqvar; ++j)
            {
                //  Compute Forward and Backward value of ufield of jth direction
                m_fields[i]->GetFwdBwdTracePhys(qfield[j][i],qFwd,qBwd);

                // if Vn >= 0, flux = uFwd, i.e.,
                //  edge::eForward, if V*n>=0 <=> V*n_F>=0, pick qflux = qBwd = q+
                //  edge::eBackward, if V*n>=0 <=> V*n_B<0, pick qflux = qBwd = q-

                // else if Vn < 0, flux = uBwd, i.e.,
                //  edge::eForward, if V*n<0 <=> V*n_F<0, pick qflux = qFwd = q-
                //  edge::eBackward, if V*n<0 <=> V*n_B>=0, pick qflux = qFwd =q+

                m_fields[i]->GetTrace()->Upwind(m_traceNormals[j],qBwd,qFwd,qfluxtemp);
                Vmath::Vmul(nTraceNumPoints,m_traceNormals[j],1,qfluxtemp,1,qfluxtemp,1);

                // Generate Stability term = - C11 ( u- - u+ )
                m_fields[i]->GetFwdBwdTracePhys(ufield[i],Fwd,Bwd);
                Vmath::Vsub(nTraceNumPoints,Fwd,1,Bwd,1,uterm,1);
                Vmath::Smul(nTraceNumPoints,-1.0*C11,uterm,1,uterm,1);

                //  Flux = {Fwd,Bwd}*(nx,ny,nz) + uterm*(nx,ny)
                Vmath::Vadd(nTraceNumPoints,uterm,1,qfluxtemp,1,qfluxtemp,1);

                // Imposing weak boundary condition with flux
                if(m_fields[0]->GetBndCondExpansions().num_elements())
                {
                    WeakPenaltyforVector(i,j,qfield[j][i],qfluxtemp,C11);
                }

                // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                // n_xi = n_x*tan_xi_x + n_y*tan_xi_y + n_z*tan_xi_z
                // n_xi = n_x*tan_eta_x + n_y*tan_eta_y + n_z*tan_eta_z
                Vmath::Vadd(nTraceNumPoints,qfluxtemp,1,qflux[i],1,qflux[i],1);
            }
        }
    }


    void EquationSystem::WeakPenaltyforScalar(const int var,
                                const Array<OneD, const NekDouble> &physfield,
                                Array<OneD, NekDouble> &penaltyflux,
                                NekDouble time)
    {
        int i, j, e, npoints, id1, id2;
        // Number of boundary regions
        int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
        int Nfps, numBDEdge;
        int nTraceNumPoints = GetTraceNpoints();
        int cnt = 0;

        Array<OneD, NekDouble > uplus(nTraceNumPoints);

        m_fields[var]->ExtractTracePhys(physfield,uplus);
        for(i = 0; i < nbnd; ++i)
        {
            // Number of boundary expansion related to that region
            numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
            // Evaluate boundary values g_D or g_N from input files
            SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(0);
            npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();
            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            for(j = 0; j < npoints; j++)
            {
                BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
            }

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                // Number of points on the expansion
                Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0) ;
                id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                // For Dirichlet boundary condition: uflux = g_D
                if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vcopy(Nfps,&BDphysics[id1],1,&penaltyflux[id2],1);
                }
                // For Neumann boundary condition: uflux = u+
                else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vcopy(Nfps,&uplus[id2],1,&penaltyflux[id2],1);
                }
            }
        }
    }

    // Diffusion: Imposing weak boundary condition for q with flux
    //  uflux = g_D  on Dirichlet boundary condition
    //  uflux = u_Fwd  on Neumann boundary condition
    void EquationSystem::WeakPenaltyforVector(const int var,
                              const int dir,
                              const Array<OneD, const NekDouble> &physfield,
                              Array<OneD, NekDouble> &penaltyflux,
                              NekDouble C11,
                              NekDouble time)
    {
        int i, j, e, npoints, id1, id2;
        int nbnd = m_fields[var]->GetBndCondExpansions().num_elements();
        int numBDEdge, Nfps;
        int nTraceNumPoints = GetTraceNpoints();
        Array<OneD, NekDouble > uterm(nTraceNumPoints);
        Array<OneD, NekDouble > qtemp(nTraceNumPoints);
        int cnt = 0;

        m_fields[var]->ExtractTracePhys(physfield,qtemp);

        for(i = 0; i < nbnd; ++i)
        {
            numBDEdge = m_fields[var]->GetBndCondExpansions()[i]->GetExpSize();
            // Evaluate boundary values g_D or g_N from input files
            SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(0);
            npoints = m_fields[var]->GetBndCondExpansions()[i]->GetNpoints();

            Array<OneD,NekDouble> BDphysics(npoints);
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[var]->GetBndCondExpansions()[i]->GetCoords(x0,x1,x2);
            for(j = 0; j < npoints; j++)
            {
                BDphysics[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
            }

            // Weakly impose boundary conditions by modifying flux values
            for (e = 0; e < numBDEdge ; ++e)
            {
                Nfps = m_fields[var]->GetBndCondExpansions()[i]->GetExp(e)->GetNumPoints(0);

                id1 = m_fields[var]->GetBndCondExpansions()[i]->GetPhys_Offset(e);
                id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

                // For Dirichlet boundary condition: qflux = q+ - C_11 (u+ - g_D) (nx, ny)
                if(m_fields[var]->GetBndConditions()[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&qtemp[id2],1,&penaltyflux[id2],1);
                }
                // For Neumann boundary condition: qflux = g_N
                else if((m_fields[var]->GetBndConditions()[i])->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                {
                    Vmath::Vmul(Nfps,&m_traceNormals[dir][id2],1,&BDphysics[id1],1,&penaltyflux[id2],1);
                }
            }
        }
    }
}
