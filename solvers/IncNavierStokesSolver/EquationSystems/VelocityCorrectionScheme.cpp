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
// Description: Velocity Correction Scheme for the Incompressible
// Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
    string VelocityCorrectionScheme::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("VelocityCorrectionScheme", VelocityCorrectionScheme::create);
    
     /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    VelocityCorrectionScheme::VelocityCorrectionScheme(
            const LibUtilities::SessionReaderSharedPtr& pSession):
        IncNavierStokes(pSession),
        m_showTimings(false)
    {
        
    }

    void VelocityCorrectionScheme::v_InitObject()
    {
        int cnt,n;

        UnsteadySystem::v_InitObject();
        IncNavierStokes::v_InitObject();
        // Set m_pressure to point to last field of m_fields; 
        if(NoCaseStringCompare(m_session->GetVariable(m_fields.num_elements()-1),"p") == 0)
        {
            m_nConvectiveFields = m_fields.num_elements()-1;
            m_pressure = m_fields[m_nConvectiveFields];
            m_pressureCalls = 1;
        }
        else
        {
            ASSERTL0(false,"Need to set up pressure field definition");
        }
        
        LibUtilities::TimeIntegrationMethod intMethod;
        std::string TimeIntStr = m_session->GetSolverInfo("TimeIntegrationMethod");
        int i;
        for(i = 0; i < (int) LibUtilities::SIZE_TimeIntegrationMethod; ++i)
        {
            if(NoCaseStringCompare(LibUtilities::TimeIntegrationMethodMap[i],TimeIntStr) == 0 )
            {
                intMethod = (LibUtilities::TimeIntegrationMethod)i; 
                break;
            }
        }
        
        ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod, "Invalid time integration type.");
        
        m_session->MatchSolverInfo("SpectralVanishingViscosity","True",m_useSpecVanVisc,false);
        m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
        m_session->LoadParameter("SVVDiffCoeff",m_sVVDiffCoeff,0.1);
        m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D","True",m_useHomo1DSpecVanVisc,false);
            
        if(m_HomogeneousType == eHomogeneous1D)
        {
            ASSERTL0(m_nConvectiveFields > 2,"Expect to have three velcoity fields with homogenous expansion");


            if(m_useHomo1DSpecVanVisc == false)
            {
                m_session->MatchSolverInfo("SpectralVanishingViscosity","True",m_useHomo1DSpecVanVisc,false);
            }

            if(m_useHomo1DSpecVanVisc)
            {
                
                Array<OneD, unsigned int> planes;
                planes = m_fields[0]->GetZIDs();

                int num_planes = planes.num_elements();
                Array<OneD, NekDouble> SVV(num_planes,0.0);
                NekDouble fac;
                int kmodes = m_fields[0]->GetHomogeneousBasis()->GetNumModes();
                int pstart;

                pstart = m_sVVCutoffRatio*kmodes;
                
                for(n = 0; n < num_planes; ++n)
                {
                    if(planes[n] > pstart)
                    {
                        fac = (NekDouble)((planes[n] - kmodes)*(planes[n] - kmodes))/
                            ((NekDouble)((planes[n] - pstart)*(planes[n] - pstart)));
                        SVV[n] = m_sVVDiffCoeff*exp(-fac)/m_kinvis;
                    }
                }

                for(i = 0; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->SetHomo1DSpecVanVisc(SVV);
                }
            }
            
        }

        m_session->MatchSolverInfo("SubSteppingScheme","True",m_subSteppingScheme,false);

        m_session->MatchSolverInfo("ShowTimings","True",m_showTimings,false);

        m_session->MatchSolverInfo("SmoothAdvection", "True",
                                   m_SmoothAdvection, false);

        if(m_subSteppingScheme)
        {
            
            ASSERTL0(m_projectionType == MultiRegions::eMixed_CG_Discontinuous,"Projection must be set to Mixed_CG_Discontinuous for substepping");
            
            m_session->LoadParameter("SubStepCFL", m_cflSafetyFactor, 0.5);
            
            // Set to 1 for first step and it will then be increased in
            // time advance routines
            switch(intMethod)
            {
            case LibUtilities::eBackwardEuler:
            case LibUtilities::eBDFImplicitOrder1: 
                {
                    m_intSteps = 1;
                    m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                    LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                    m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                        
                    LibUtilities::TimeIntegrationSchemeKey     SubIntKey(LibUtilities::eForwardEuler);
                    m_subStepIntegrationScheme = LibUtilities::TimeIntegrationSchemeManager()[SubIntKey];
                    
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
                    m_intSteps = 2;
                    m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                    
                    
                    LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eBackwardEuler);
                    m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
                    m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    
                    LibUtilities::TimeIntegrationSchemeKey     SubIntKey(LibUtilities::eRungeKutta2_ImprovedEuler);
                    
                    m_subStepIntegrationScheme = LibUtilities::TimeIntegrationSchemeManager()[SubIntKey];
                    
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
            
            // set explicit time-intregration class operators
                m_subStepIntegrationOps.DefineOdeRhs(&IncNavierStokes::SubStepAdvection, this);
                m_subStepIntegrationOps.DefineProjection(&IncNavierStokes::SubStepProjection, this);
                
        }
        else // Standard velocity correction scheme
        {
            
            // Set to 1 for first step and it will then be increased in
            // time advance routines
            switch(intMethod)
            {
                case LibUtilities::eIMEXOrder1: 
                {
                    m_intSteps = 1;
                    m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                    LibUtilities::TimeIntegrationSchemeKey       IntKey0(intMethod);
                    m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                }
                break;
                case LibUtilities::eIMEXOrder2: 
                {
                    m_intSteps = 2;
                    m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                    LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXOrder1);
                    m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    LibUtilities::TimeIntegrationSchemeKey       IntKey1(intMethod);
                    m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                }
                break;
                case LibUtilities::eIMEXOrder3: 
                {
                    m_intSteps = 3;
                    m_integrationScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> (m_intSteps);
                    LibUtilities::TimeIntegrationSchemeKey       IntKey0(LibUtilities::eIMEXdirk_3_4_3);
                    m_integrationScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
                    LibUtilities::TimeIntegrationSchemeKey       IntKey1(LibUtilities::eIMEXdirk_3_4_3);
                    m_integrationScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];
                    LibUtilities::TimeIntegrationSchemeKey       IntKey2(intMethod);
                    m_integrationScheme[2] = LibUtilities::TimeIntegrationSchemeManager()[IntKey2];
                }
                break;
                default:
                    ASSERTL0(0,"Integration method not suitable: Options include IMEXOrder1, IMEXOrder2 or IMEXOrder3");
                    break;
            }
            
            // set explicit time-intregration class operators
            m_integrationOps.DefineOdeRhs(&VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs, this);
            
        }
        // Count number of HBC conditions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds = m_pressure->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp = m_pressure->GetBndCondExpansions();
        
        // Set up mapping from pressure boundary condition to pressure
        // element details.
        m_pressure->GetBoundaryToElmtMap(m_pressureBCtoElmtID,m_pressureBCtoTraceID);
	
        // find the maximum values of points  for pressure BC evaluation
        m_pressureBCsMaxPts = 0; 
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            for(int i = 0; i < PBndExp[n]->GetExpSize(); ++i)
            {
                m_pressureBCsMaxPts = max(m_pressureBCsMaxPts, m_fields[0]->GetExp(m_pressureBCtoElmtID[cnt++])->GetTotPoints());
            }
        }

        // Storage array for high order pressure BCs
        m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
        m_acceleration = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);
        
        int HBCnumber = 0;
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                cnt += PBndExp[n]->GetNcoeffs();
                HBCnumber += PBndExp[n]->GetExpSize();
            }
        }
        
        if (HBCnumber > 0) 
        {
            m_acceleration[0] = Array<OneD, NekDouble>(cnt, 0.0);

            for(n = 0; n < m_intSteps; ++n)
            {
                m_pressureHBCs[n]   = Array<OneD, NekDouble>(cnt, 0.0);
                m_acceleration[n+1] = Array<OneD, NekDouble>(cnt, 0.0);
            }
        }

        // creating a Map to store the information regarding
        // High-Order pressure BCs to improve efficiency
        // currently only set up for Fourier expansions
        FillHOPBCMap(HBCnumber);
        
        // set implicit time-intregration class operators
        m_integrationOps.DefineImplicitSolve(&VelocityCorrectionScheme::SolveUnsteadyStokesSystem,this);
    }
        
    VelocityCorrectionScheme::~VelocityCorrectionScheme(void)
    {
        
    }
    
        
    void VelocityCorrectionScheme::v_PrintSummary(std::ostream &out)
    {
        if(m_subSteppingScheme)
        {
            cout <<  "\tSolver Type     : Velocity Correction with Substepping" <<endl;
        }
        else
        {
            cout <<  "\tSolver Type     : Velocity Correction" <<endl;
        }

        if(m_session->DefinesSolverInfo("EvolutionOperator"))
        {
            cout << "\tEvolutionOp     : " << m_session->GetSolverInfo("EvolutionOperator")<< endl;
            
        }
        else
        {
            cout << "\tEvolutionOp     : " << endl;
        }

        if(m_session->DefinesSolverInfo("Driver"))
        {
            cout << "\tDriver          : " << m_session->GetSolverInfo("Driver")<< endl;
            
        }
        else
        {
            cout << "\tDriver          : "<< endl;
        }
        
        TimeParamSummary(out);
        cout << "\tTime integ.     : " << LibUtilities::TimeIntegrationMethodMap[m_integrationScheme[m_intSteps-1]->GetIntegrationMethod()] << endl;
        
        if(m_subSteppingScheme)
        {
            cout << "\tSubstepping     : " << LibUtilities::TimeIntegrationMethodMap[m_subStepIntegrationScheme->GetIntegrationMethod()] << endl;
        }

        if(m_dealiasing)
        {
            cout << "\tDealiasing      : Homogeneous1D"  << endl;
        }
        
        if(m_specHP_dealiasing)
        {
            cout << "\tDealiasing      : Spectral/hp "  << endl;
        }

        if(m_useSpecVanVisc)
        {
            cout << "\tSmoothing       : Spectral vanishing viscosity (cut off ratio = " << m_sVVCutoffRatio << ", diff coeff = "<< m_sVVDiffCoeff << ")"<< endl;        
        }

        if(m_useHomo1DSpecVanVisc)
        {
            cout << "\tSmoothing       : Spectral vanishing viscosity (homogeneous1D, cut off ratio = " << m_sVVCutoffRatio << ", diff coeff = "<< m_sVVDiffCoeff << ")"<< endl;  
        }
    }

    void VelocityCorrectionScheme::v_DoInitialise(void)
    {

        UnsteadySystem::v_DoInitialise();

        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->LocalToGlobal();
            m_fields[i]->ImposeDirichletConditions(m_fields[i]->UpdateCoeffs());
            m_fields[i]->GlobalToLocal();
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
        }
        
        //insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);
	
        m_session->LoadParameter("Noise", Noise,0.0);
	
        if(Noise > 0.0)
        {
            for(int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot,Noise,noise,1,m_comm->GetColumnComm()->GetRank()+1);
                Vmath::Vadd(phystot,m_fields[i]->GetPhys(),1,noise,1,m_fields[i]->UpdatePhys(),1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
        }
    }
    
    void VelocityCorrectionScheme::v_DoSolve(void)
    {
        switch(m_equationType)
        {
        case eUnsteadyStokes: 
        case eUnsteadyNavierStokes:
        case eUnsteadyLinearisedNS:
            {  
                // Integrate from start time to end time
                AdvanceInTime(m_steps);
                break;
            }
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type for VelocityCorrectionScheme");
        }
    }

    
    void VelocityCorrectionScheme:: v_TransCoeffToPhys(void)
    {
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Backward Transformation in physical space for time evolution
            m_fields[k]->BwdTrans_IterPerExp(m_fields[k]->GetCoeffs(),
                                             m_fields[k]->UpdatePhys());
        }
    }
    
    void VelocityCorrectionScheme:: v_TransPhysToCoeff(void)
    {
        
        int nfields = m_fields.num_elements() - 1;
        for (int k=0 ; k < nfields; ++k)
        {
            //Forward Transformation in physical space for time evolution
            m_fields[k]->FwdTrans_IterPerExp(m_fields[k]->GetPhys(),m_fields[k]->UpdateCoeffs());
        }
    }
	
    Array<OneD, bool> VelocityCorrectionScheme::v_GetSystemSingularChecks()
    {
        int vVar = m_session->GetVariables().size();
        Array<OneD, bool> vChecks(vVar, false);
        vChecks[vVar-1] = true;
        return vChecks;
    }
    
    int VelocityCorrectionScheme::v_GetForceDimension()
    {
        return m_session->GetVariables().size() - 1;
    }

    void VelocityCorrectionScheme::EvaluateAdvection_SetPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                                    Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                                    const NekDouble time)
    {
        int nqtot        = m_fields[0]->GetTotPoints();
        
        Timer  timer;
        bool IsRoot = (m_comm->GetColumnComm()->GetRank())? false:true;

#if 0
        timer.Start();
        for(int k = 0; k < 1000; ++k)
        {
            m_fields[0]->IProductWRTBase(m_fields[0]->GetPhys(),
                                         m_fields[0]->UpdateCoeffs());

        }
        timer.Stop();
        cout << "\t 1000 Iprods   : "<< timer.TimePerTest(1) << endl;
#endif

#if 0
        timer.Start();
        Array<OneD, NekDouble> out (m_fields[0]->GetTotPoints());
        Array<OneD, NekDouble> out1(m_fields[0]->GetTotPoints());
        Array<OneD, NekDouble> out2(m_fields[0]->GetTotPoints());
        
        for(int k = 0; k < 10000; ++k)
        {
            m_fields[0]->PhysDeriv(out,out1,out2);

        }
        timer.Stop();
        cout << "\t 10000 Physderiv   : "<< timer.TimePerTest(1) << endl;
        exit(1);
#endif

        timer.Start();
        // evaluate convection terms
        m_advObject->DoAdvection(m_fields, m_nConvectiveFields, m_velocity,inarray,outarray,m_time);
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t Adv. eval Time   : "<< timer.TimePerTest(1) << endl;
        }
		
        if(m_SmoothAdvection && m_pressureCalls > 1)
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_pressure->SmoothField(outarray[i]);
            }
        }

        //add the force
        if(m_session->DefinesFunction("BodyForce"))
        {
            timer.Start();
            if(m_fields[0]->GetWaveSpace())
            {
                for(int i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_forces[i]->SetWaveSpace(true);					
                    m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),
                                          m_forces[i]->UpdatePhys());
                }
            }
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,
                            outarray[i],1);
            }
            timer.Stop();
            if(m_showTimings&&IsRoot)
            {
                cout << "\t Body ForceTime   : "<< timer.TimePerTest(1) << endl;
            }
        }
        
        if(m_pressureHBCs[0].num_elements() > 0)
        {
            // Set pressure BCs
            timer.Start();
            EvaluatePressureBCs(inarray, outarray); 
            timer.Stop();

            if(m_showTimings&&IsRoot)
            {
                cout << "\t Pressure BCs     : "<< timer.TimePerTest(1) << endl;
            }
        }
    }
    
    void VelocityCorrectionScheme::SolveUnsteadyStokesSystem(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                             Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                             const NekDouble time, 
                                                             const NekDouble aii_Dt)
    {
        int i,n;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, Array< OneD, NekDouble> > F(m_nConvectiveFields);
        StdRegions::ConstFactorMap factors;
        factors[StdRegions::eFactorLambda] = 0.0;
        Timer timer;
        bool IsRoot = (m_comm->GetColumnComm()->GetRank())? false:true;

        for(n = 0; n < m_nConvectiveFields; ++n)
        {
            F[n] = Array<OneD, NekDouble> (phystot);
        }
            
        SetBoundaryConditions(time);
        
        if(m_subSteppingScheme)
        {
            SubStepSetPressureBCs(inarray,aii_Dt);
        }
		
        // Pressure Forcing = Divergence Velocity; 
        timer.Start();
        SetUpPressureForcing(inarray, F, aii_Dt);
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t Pressure Forcing : "<< timer.TimePerTest(1) << endl;
	}

        // Solver Pressure Poisson Equation 
        timer.Start();
        m_pressure->HelmSolve(F[0], m_pressure->UpdateCoeffs(), NullFlagList, factors);
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t Pressure Solve   : "<< timer.TimePerTest(1) << endl;
        }

        // Viscous Term forcing
        timer.Start();
        SetUpViscousForcing(inarray, F, aii_Dt);
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t Viscous Forcing  : "<< timer.TimePerTest(1) << endl;
        }
        factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_kinvis;
        
        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_kinvis;
        }
        
        // Solve Helmholtz system and put in Physical space
        timer.Start();
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->HelmSolve(F[i], m_fields[i]->UpdateCoeffs(), NullFlagList, factors);            
        }
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t Visc. Slv.      : "<< timer.TimePerTest(1) << endl;
        }

        timer.Start();

        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
 
        timer.Stop();
        if(m_showTimings&&IsRoot)
        {
            cout << "\t BwdTrans : "<< timer.TimePerTest(1) << endl;
        
        }
    }
    
    /** 
     * 
     */
    void VelocityCorrectionScheme::SubStepSetPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &inarray, const NekDouble Aii_DT)
    {
        
        if(m_pressureHBCs[0].num_elements() > 0)
        {
            Array<OneD, Array<OneD, NekDouble> > velfields(m_velocity.num_elements());
            
            for(int i = 0; i < m_velocity.num_elements(); ++i)
            {
                velfields[i] = m_fields[m_velocity[i]]->GetPhys(); 
            }
            // Set pressure BCs
            EvaluatePressureBCs(velfields, inarray, Aii_DT);
        }
    }
    
    void VelocityCorrectionScheme::EvaluatePressureBCs(const Array<OneD, const Array<OneD, NekDouble> >  &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N, const NekDouble Aii_Dt)
    {		
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> accelerationTerm;
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        int  n,cnt;
        int  nint    = min(m_pressureCalls++,m_intSteps);
        int  nlevels = m_pressureHBCs.num_elements();
        
        PBndConds   = m_pressure->GetBndConditions();
        PBndExp     = m_pressure->GetBndCondExpansions();

        int acc_order = 0;
        accelerationTerm = Array<OneD, NekDouble>(m_acceleration[0].num_elements(), 0.0);

        // Rotate HOPBCs storage
        Rotate(m_pressureHBCs);

        // Rotate acceleration term
        Rotate(m_acceleration);

        // Calculate BCs at current level
        CalcPressureBCs(fields,N);
        
        // Copy High order values into storage array 
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                int nq = PBndExp[n]->GetNcoeffs();
                Vmath::Vcopy(nq,&(PBndExp[n]->GetCoeffs()[0]),1,&(m_pressureHBCs[0])[cnt],1);
                cnt += nq;
            }
        }

        //Calculate acceleration term at level n based on previous steps
        if (m_pressureCalls > 2)
        {
            acc_order = min(m_pressureCalls-2,m_intSteps);
            Vmath::Smul(cnt, StifflyStable_Gamma0_Coeffs[acc_order-1],
                             m_acceleration[0], 1,
                             accelerationTerm,  1);

            for(int i = 0; i < acc_order; i++)
            {
                Vmath::Svtvp(cnt, -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                                  m_acceleration[i+1], 1,
                                  accelerationTerm,    1,
                                  accelerationTerm,    1);
            }
        }
        
        // Adding acceleration term to HOPBCs
        Vmath::Svtvp(cnt, -1.0/m_timestep,
                          accelerationTerm,  1,
                          m_pressureHBCs[0], 1,
                          m_pressureHBCs[0], 1);

        // Extrapolate to n+1
        Vmath::Smul(cnt, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                         m_pressureHBCs[nint-1],    1,
                         m_pressureHBCs[nlevels-1], 1);

        for(n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(cnt,StifflyStable_Betaq_Coeffs[nint-1][n],
                         m_pressureHBCs[n],1,m_pressureHBCs[nlevels-1],1,
                         m_pressureHBCs[nlevels-1],1);
        }

        // Copy values of [dP/dn]^{n+1} in the pressure bcs storage.
        // m_pressureHBCS[nlevels-1] will be cancelled at next time step
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
            {
                int nq = PBndExp[n]->GetNcoeffs();
                Vmath::Vcopy(nq,&(m_pressureHBCs[nlevels-1])[cnt],1,&(PBndExp[n]->UpdateCoeffs()[0]),1);
                cnt += nq;
            }
        }
        
        if(m_subSteppingScheme)
        {
            AddDuDt(N,Aii_Dt);
        }
    }


    void VelocityCorrectionScheme::Rotate(
            Array<OneD, Array<OneD, NekDouble> > &input)
    {
        int  nlevels = input.num_elements();
        Array<OneD, NekDouble> tmp;

        tmp = input[nlevels-1];

        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }

        input[0] = tmp;
    }


    void VelocityCorrectionScheme::CalcPressureBCs(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        //  switch(m_nConvectiveFields)
        switch(m_velocity.num_elements())
        {
        case 1:
            ASSERTL0(false,"Velocity correction scheme not designed to have just one velocity component");
            break;
        case 2:
            CalcPressureBCs2D(fields,N);
            break;
        case 3:
            CalcPressureBCs3D(fields,N);
            break;
        }
    }
    
    void VelocityCorrectionScheme::CalcPressureBCs2D(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        int  i,n;
        
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansion1DSharedPtr Pbc;
        
        Array<OneD, NekDouble> Pvals;
        
        int cnt;
        int elmtid,nq,offset, boundary;
        
        Array<OneD, const NekDouble> U,V,Nu,Nv;
        // Elemental tempspace: call once adn set others as offset 
        Array<OneD, NekDouble> Uy(5*m_pressureBCsMaxPts); 
        Array<OneD, NekDouble> Vx = Uy + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Qx = Vx + m_pressureBCsMaxPts;
        Array<OneD, NekDouble> Qy = Qx + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Q  = Qy + m_pressureBCsMaxPts;
        
        for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
        {            
            SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined(); 
            
            if(type == SpatialDomains::eHigh)
            {
                for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    // find element and edge of this expansion. 
                    // calculate curl x curl v;
                    elmtid = m_pressureBCtoElmtID[cnt];
                    elmt   = m_fields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = m_fields[0]->GetPhys_Offset(elmtid);
                    
                    U = fields[m_velocity[0]] + offset;
                    V = fields[m_velocity[1]] + offset; 
                    
                    // Calculating vorticity Q = (dv/dx - du/dy)
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],V,Vx);
                    elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],U,Uy);  
                    
                    Vmath::Vsub(nq,Vx,1,Uy,1,Q,1);
                    
                    // Calculate  Curl(Q) = Qy i - Qx j 
                    elmt->PhysDeriv(Q,Qx,Qy);
                    
                    Nu = N[0] + offset;
                    Nv = N[1] + offset; 
                    
                    if(m_subSteppingScheme)
                    {
                        // Evaluate [- kinvis Curlx Curl V].n
                        // x-component (stored in Qy)
                        Vmath::Smul(nq,-m_kinvis,Qy,1,Qy,1);
                        
                        // y-component (stored in Qx )
                        Vmath::Smul(nq,m_kinvis,Qx,1,Qx,1);
                    }
                    else
                    {
                        // Evaluate [N - kinvis Curlx Curl V].n
                        // x-component (stored in Qy)
                        //Vmath::Zero(nq,Qy,1);
                        Vmath::Svtvp(nq,-m_kinvis,Qy,1,Nu,1,Qy,1);
                        
                        // y-component (stored in Qx )
                        //Vmath::Zero(nq,Qx,1);
                        Vmath::Svtvp(nq,m_kinvis,Qx,1,Nv,1,Qx,1);
                    }
                    Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[n]->GetExp(i));
                    
                    boundary = m_pressureBCtoTraceID[cnt];
                    
                    // Get edge values and put into Uy, Vx
                    elmt->GetEdgePhysVals(boundary,Pbc,Qy,Uy);
                    elmt->GetEdgePhysVals(boundary,Pbc,Qx,Vx);
                    
                    // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n) 
                    Pvals = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i);
                    Pbc->NormVectorIProductWRTBase(Uy,Vx,Pvals); 
                }
            }

            // setting if just standard BC not High order
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eTimeDependent) 
            {
                cnt += PBndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }
    }
    
    void VelocityCorrectionScheme::CalcPressureBCs3D(const Array<OneD, const Array<OneD, NekDouble> > &fields, const Array<OneD, const Array<OneD, NekDouble> >  &N)
    {
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp;
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
        int elmtid,nq,offset, boundary, n;
        
        Array<OneD, NekDouble> Pvals;
        Array<OneD, NekDouble> Uvals;
        Array<OneD, const NekDouble> U,V,W,Nu,Nv,Nw;

        Array<OneD, NekDouble> Uy(9*m_pressureBCsMaxPts);
        Array<OneD, NekDouble> Uz = Uy + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Vx = Uz + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Vz = Vx + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Wx = Vz + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Wy = Wx + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Qx = Wy + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Qy = Qx + m_pressureBCsMaxPts; 
        Array<OneD, NekDouble> Qz = Qy + m_pressureBCsMaxPts; 
            
        if(m_HomogeneousType == eHomogeneous1D)
        {
            Array<OneD, NekDouble> Wz(9*m_pressureBCsMaxPts);
            Array<OneD, NekDouble> Vxx = Wz  + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Vzz = Vxx + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Vxy = Vzz + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Uyy = Vxy + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Uzz = Uyy + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Uxy = Uzz + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Wxz = Uxy + m_pressureBCsMaxPts; 
            Array<OneD, NekDouble> Wyz = Wxz + m_pressureBCsMaxPts;

            Array<OneD, NekDouble> Ub(m_pressureBCsMaxPts);
            Array<OneD, NekDouble> Vb(m_pressureBCsMaxPts);

            StdRegions::StdExpansion1DSharedPtr Pbc;
            
            for(int j = 0 ; j < m_HBCdata.num_elements() ; j++)
            {
                Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (PBndExp[m_HBCdata[j].m_bndryElmtID]->GetExp(m_HBCdata[j].m_bndElmtOffset));
                
                // Picking up the element where the HOPBc is located
                m_elmt = m_fields[0]->GetExp(m_HBCdata[j].m_globalElmtID);
                
                // Using the physical offset to get the velocity (W is taken on the coniugated plane)
                U = fields[m_velocity[0]] + m_HBCdata[j].m_physOffset;
                V = fields[m_velocity[1]] + m_HBCdata[j].m_physOffset;
                W = fields[m_velocity[2]] + m_HBCdata[j].m_assPhysOffset;
                
                // Derivatives to build up the curl curl of the velocity
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],V,Vx);
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],U,Uy);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_wavenumber[j],W,1,Wz,1);
                
                // x-components of vorticity curl
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],Vx,Vxy);
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],Uy,Uyy);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],U,1,Uzz,1);
                    
                //x-component coming from the other plane
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],Wz,Wxz);
                
                // y-components of vorticity curl
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],Vx,Vxx);
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[0],Uy,Uxy);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt,m_negWavenumberSq[j],V,1,Vzz,1);
                
                //y-component coming from the other plane
                m_elmt->PhysDeriv(MultiRegions::DirCartesianMap[1],Wz,Wyz);
                
                // buinding up the curl of V adding the components
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Vxy,1,Uyy,1,Qx,1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Qx ,1,Uzz,1,Qx,1);
                Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Qx ,1,Wxz,1,Qx,1);
                
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Wyz,1,Vzz,1,Qy,1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt,Qy ,1,Vxx,1,Qy,1);
                Vmath::Vadd(m_HBCdata[j].m_ptsInElmt,Qy ,1,Uxy,1,Qy,1);
		
                // getting the advective term
                Nu = N[0] + m_HBCdata[j].m_physOffset;
                Nv = N[1] + m_HBCdata[j].m_physOffset;

                if(m_subSteppingScheme)
                {
                    // Evaluate [- kinvis Curlx Curl V]
                    // x-component (stored in Qx)
                    Vmath::Smul(m_HBCdata[j].m_ptsInElmt,-m_kinvis,Qx,1,Qx,1);
                    // y-component (stored in Qy)
                    Vmath::Smul(m_HBCdata[j].m_ptsInElmt,-m_kinvis,Qy,1,Qy,1);
                    // z-component (stored in Qz) not required for this approach
                    // the third component of the normal vector is always zero
                }
                else
                {
                    // Evaluate [N - kinvis Curlx Curl V]
                    // x-component (stored in Qx)
                    Vmath::Svtvp(m_HBCdata[j].m_ptsInElmt,-m_kinvis,Qx,1,Nu,1,Qx,1);
                    // y-component (stored in Qy)
                    Vmath::Svtvp(m_HBCdata[j].m_ptsInElmt,-m_kinvis,Qy,1,Nv,1,Qy,1);
                    // z-component (stored in Qz) not required for this approach
                    // the third component of the normal vector is always zero
                }
                
                // Get edge values and put into Uy, Vx
                m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Qy,Uy);
                m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,Qx,Vx);
                
                // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n) 
                Pvals = PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset);
                
                Pbc->NormVectorIProductWRTBase(Vx,Uy,Pvals);

                //update values for the acceleration term
                m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,U,Ub);
                m_elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID,Pbc,V,Vb);
                Uvals = (m_acceleration[0]) + m_HBCdata[j].m_coeffOffset;
                Pbc->NormVectorIProductWRTBase(Ub,Vb,Uvals);
            }
        }
        else if(m_HomogeneousType == eHomogeneous2D)
        {
            int phystot = m_fields[0]->GetTotPoints();
            Array< OneD, NekDouble>  Q(phystot);
            Array< OneD, NekDouble>  wx(phystot);
            Array< OneD, NekDouble>  vx(phystot);
            Array< OneD, NekDouble>  uy(phystot);
            Array< OneD, NekDouble>  uz(phystot);
            Array< OneD, NekDouble>  qx(phystot);
            Array< OneD, NekDouble>  qy(phystot);
            Array< OneD, NekDouble>  qz(phystot);
            
            //Vectors names are misleading, some of them have been re-used to save memory
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[0],fields[2],wx);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[0],fields[1],vx);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[1],fields[0],uy);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[2],fields[0],uz);
            Vmath::Vsub(phystot,uz,1,wx,1,qy,1);
            Vmath::Vsub(phystot,vx,1,uy,1,qz,1);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[1],qz,uy);
            m_pressure->PhysDeriv(MultiRegions::DirCartesianMap[2],qy,uz);
            Vmath::Vsub(phystot,uy,1,uz,1,qx,1);
            
            if(m_subSteppingScheme)
            {
                Vmath::Smul(phystot,-m_kinvis,qx,1,qx,1);
            }
            else
            {
                Vmath::Svtvp(phystot,-m_kinvis,qx,1,N[0],1,qx,1);
            }
		
            if(m_pressure->GetWaveSpace())
            {
                Q = qx;
            }
            else
            {
                m_pressure->HomogeneousFwdTrans(qx,Q);
            }
            
            for(int j = 0 ; j < m_HBCdata.num_elements(); j++)
            {				
                Qx = Q + m_HBCdata[j].m_physOffset;
                
                if(m_HBCdata[j].m_elmtTraceID == 0)
                {
                    (PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset))[0] = -1.0*Qx[0];
                }
                else if (m_HBCdata[j].m_elmtTraceID == 1)
                {
                    (PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()+PBndExp[m_HBCdata[j].m_bndryElmtID]->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset))[0] = Qx[m_HBCdata[j].m_ptsInElmt-1];
                }
                else 
                {
                    ASSERTL0(false,"In the 3D homogeneous 2D approach BCs edge ID can be just 0 or 1 ");
                }
            }
        }
        else if(m_HomogeneousType == eHomogeneous3D)
        {
            ASSERTL0(false,"High Order Pressure BC not required for this approach");
        }
        // Full 3D		
        else
        {
            int i, cnt;
            
            StdRegions::StdExpansionSharedPtr elmt;
            StdRegions::StdExpansion2DSharedPtr Pbc;
            Array<OneD,NekDouble> NoCalc; // Dummy array

            for(cnt = n = 0; n < PBndConds.num_elements(); ++n)
            {
                
                SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined();
                
                if(type == SpatialDomains::eHigh)
                {
                    for(i = 0; i < PBndExp[n]->GetExpSize(); ++i,cnt++)
                    {
                        // find element and face of this expansion. 
                        // calculate curl x curl v;
                        elmtid = m_pressureBCtoElmtID[cnt];
                        elmt   = m_fields[0]->GetExp(elmtid);
                        nq     = elmt->GetTotPoints();
                        offset = m_fields[0]->GetPhys_Offset(elmtid);
                        
                        U = fields[m_velocity[0]] + offset;
                        V = fields[m_velocity[1]] + offset; 
                        W = fields[m_velocity[2]] + offset;
                        
                        // Calculating vorticity Q = (dv/dx - du/dy)
                        elmt->PhysDeriv(U,NoCalc,Uy,Uz);
                        elmt->PhysDeriv(V,Vx,NoCalc,Vz);
                        elmt->PhysDeriv(W,Wx,Wy,NoCalc);

                        Vmath::Vsub(nq,Wy,1,Vz,1,Qx,1);
                        Vmath::Vsub(nq,Uz,1,Wx,1,Qy,1);
                        Vmath::Vsub(nq,Vx,1,Uy,1,Qz,1);
                        
                        // Calculate  NxQ = Curl(Q) = (Qzy-Qyz) i + (Qxz-Qzx) j 
                        //                 + (Qyx-Qxy) k
                        // NxQ = NxQ_x i + NxQ_y j + NxQ_z k
                        // Using the velocity derivatives memory space to
                        // store the vorticity derivatives.
                        // Qzy => Uy,  Qyz => U,  Qxz => Vx 
                        // Qzx => Vz,  Qyx => Wx, Qxy => Wy 

                        elmt->PhysDeriv(Qx,NoCalc,Wy,Vx);
                        elmt->PhysDeriv(Qy,Wx,NoCalc,Uz);
                        elmt->PhysDeriv(Qz,Vz,Uy,NoCalc);
                        
                        // Using the storage space associated with
                        // the 3 components of the vorticity to
                        // store the 3 components od the vorticity
                        // curl to save space Qx = Qzy-Qyz = Uy-Uz
                        // 
                        // Qy = Qxz-Qzx = Vx-Vz 
                        // Qz=  Qyx-Qxy = Wx-Wy
                        Vmath::Vsub(nq,Uy,1,Uz,1,Qx,1);
                        Vmath::Vsub(nq,Vx,1,Vz,1,Qy,1);
                        Vmath::Vsub(nq,Wx,1,Wy,1,Qz,1);
                        
                        Nu = N[0] + offset;
                        Nv = N[1] + offset;
                        Nw = N[2] + offset;
                        
                        if(m_subSteppingScheme)
                        {
                            // Evaluate [- kinvis Curlx Curl V]
                            // x-component (stored in Qx)
                            Vmath::Smul(nq,-m_kinvis,Qx,1,Qx,1);
                            // y-component (stored in Qy)
                            Vmath::Smul(nq,-m_kinvis,Qy,1,Qy,1);
                            // z-component (stored in Qz)
                            Vmath::Smul(nq,-m_kinvis,Qz,1,Qz,1);		
                        }
                        else
                        {
                            // Evaluate [N - kinvis Curlx Curl V]
                            // x-component (stored in Qx)
                            Vmath::Svtvp(nq,-m_kinvis,Qx,1,Nu,1,Qx,1);
                            // y-component (stored in Qy)
                            Vmath::Svtvp(nq,-m_kinvis,Qy,1,Nv,1,Qy,1);
                            // z-component (stored in Qz)
                            Vmath::Svtvp(nq,-m_kinvis,Qz,1,Nw,1,Qz,1);		
                        }
                        
                        Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (PBndExp[n]->GetExp(i));
                            
                        boundary = m_pressureBCtoTraceID[cnt];
                        // Get face values and put into Uy, Vx and Wx
                        elmt->GetFacePhysVals(boundary,Pbc,Qx,Uy);
                        elmt->GetFacePhysVals(boundary,Pbc,Qy,Vx);
                        elmt->GetFacePhysVals(boundary,Pbc,Qz,Wx);
                        
                        // calcuate (phi, dp/dn = [N-kinvis curl x curl v].n) 
                        Pvals = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i);
                        Pbc->NormVectorIProductWRTBase(Uy,Vx,Wx,Pvals); 

                    }
                }
                // setting if just standard BC no High order
                else if(type == SpatialDomains::eNoUserDefined || 
                        type == SpatialDomains::eTimeDependent)
                {
                    cnt += PBndExp[n]->GetExpSize();
                }
                else
                {
                    ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
                }
            }
        }
    }
    

    void VelocityCorrectionScheme::AddDuDt(const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt)
    {
        switch(m_velocity.num_elements())
        {
        case 1:
            ASSERTL0(false,"Velocity correction scheme not designed to have just one velocity component");
            break;
        case 2:
            AddDuDt2D(N,Aii_Dt);
            break;
        case 3:
            AddDuDt3D(N,Aii_Dt);
            break;
        }
    }
        
        
    void VelocityCorrectionScheme::AddDuDt2D(const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt)
    {
        int i,n;
        ASSERTL0(m_velocity.num_elements() == 2," Routine currently only set up for 2D");
        
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp,UBndExp,VBndExp;
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
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
            SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined(); 
            
            if(type == SpatialDomains::eHigh)
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
                    
                    Vmath::Vsub(ncoeffs,Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),1, Pvals,1, Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),1);
                }
            }
            // setting if just standard BC no High order
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eTimeDependent) 
            {
                cnt += PBndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }        
    }
    
        
    void VelocityCorrectionScheme::AddDuDt3D(const Array<OneD, const Array<OneD, NekDouble> >  &N, NekDouble Aii_Dt)
    {
        int i,n;
        ASSERTL0(m_velocity.num_elements() == 3," Routine currently only set up for 3D");
        
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp,UBndExp,VBndExp,WBndExp;
        
        PBndConds = m_pressure->GetBndConditions();
        PBndExp   = m_pressure->GetBndCondExpansions();
        
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
            SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined(); 
            
            if(type == SpatialDomains::eHigh)
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
                    UBndExp[n]->GetExp(i)->BwdTrans(UBndExp[n]->GetCoeffs() + UBndExp[n]->GetCoeff_Offset(i),ubc);
                    VBndExp[n]->GetExp(i)->BwdTrans(VBndExp[n]->GetCoeffs() + VBndExp[n]->GetCoeff_Offset(i),vbc);
                    WBndExp[n]->GetExp(i)->BwdTrans(WBndExp[n]->GetCoeffs() + WBndExp[n]->GetCoeff_Offset(i),wbc);
                    
                    // Get edge values and put into Nu,Nv
                    elmt->GetFacePhysVals(boundary,Pbc,Nu,N1);
                    elmt->GetFacePhysVals(boundary,Pbc,Nv,N2);
                    elmt->GetFacePhysVals(boundary,Pbc,Nw,N3);
                    
                    
                    // Take different as Forward Euler but N1,N2
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
                    
                    Vmath::Vsub(ncoeffs,Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),1, Pvals,1, Ptmp = PBndExp[n]->UpdateCoeffs()+PBndExp[n]->GetCoeff_Offset(i),1);
                }
            }
            // setting if just standard BC no High order
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eTimeDependent) 
            {
                cnt += PBndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }        
    }
    
    // Evaluate divergence of velocity field. 
    void   VelocityCorrectionScheme::SetUpPressureForcing(const Array<OneD, const Array<OneD, NekDouble> > &fields, Array<OneD, Array<OneD, NekDouble> > &Forcing, const NekDouble aii_Dt)
    {                
        int   i;
        int   physTot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> wk = Array<OneD, NekDouble>(physTot);
        int nvel = m_velocity.num_elements();
        
#if 1
        Vmath::Zero(physTot,Forcing[0],1);
        
        for(i = 0; i < nvel; ++i)
        {
            m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[i], wk);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
        }
#else
        if(nvel == 2)
        {
            m_fields[0]->PhysDeriv(fields[0],Forcing[0],NullNekDouble1DArray);
            m_fields[1]->PhysDeriv(fields[1],NullNekDouble1DArray,wk);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
        }
        else
        {
            m_fields[0]->PhysDeriv(fields[0],Forcing[0],NullNekDouble1DArray,NullNekDouble1DArray);
            m_fields[1]->PhysDeriv(fields[1],NullNekDouble1DArray,wk,NullNekDouble1DArray);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
            m_fields[2]->PhysDeriv(fields[2],NullNekDouble1DArray,NullNekDouble1DArray,wk);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
        }
#endif  
        Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);        
    }
    
    void   VelocityCorrectionScheme::SetUpViscousForcing(const Array<OneD, const Array<OneD, NekDouble> > &inarray, Array<OneD, Array<OneD, NekDouble> > &Forcing, const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int phystot = m_fields[0]->GetTotPoints();

        // Grad p
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
        
        int nvel = m_velocity.num_elements();
        if(nvel == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1],
                                  Forcing[2]);
        }
        
        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < nvel; ++i)
        {
            Blas::Daxpy(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            Blas::Dscal(phystot,1.0/m_kinvis,&(Forcing[i])[0],1);
        }
    }
	
    
    void VelocityCorrectionScheme::FillHOPBCMap(const int HOPBCnumber)
    {
            
        // Count number of HBC conditions
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > PBndConds = m_pressure->GetBndConditions();
        Array<OneD, MultiRegions::ExpListSharedPtr>  PBndExp = m_pressure->GetBndCondExpansions();
        
        ////////////////////////////////////////////////////////////////////////////
        if(m_HomogeneousType == eHomogeneous1D)
        {

            Array<OneD, unsigned int> planes;
            
            planes = m_fields[0]->GetZIDs();
            
            int num_planes = planes.num_elements();            
            int num_elm_per_plane = (m_fields[0]->GetExpSize())/num_planes;
            
            m_HBCdata = Array<OneD, HBCInfo>(HOPBCnumber);
            
            m_wavenumber      = Array<OneD, NekDouble>(HOPBCnumber);
            m_negWavenumberSq = Array<OneD, NekDouble>(HOPBCnumber);

            int coeff_count = 0;
            int exp_size, exp_size_per_plane;
            int j=0;
            int K;
            NekDouble sign = -1.0;
            int cnt = 0;

            for(int k = 0; k < num_planes; k++)
            {
                K = planes[k]/2;
                for(int n = 0 ; n < PBndConds.num_elements(); ++n)
                {
                    exp_size = PBndExp[n]->GetExpSize();
                    exp_size_per_plane = exp_size/num_planes;
                    if(PBndConds[n]->GetUserDefined() == SpatialDomains::eHigh)
                    {
                        for(int i = 0; i < exp_size_per_plane; ++i,cnt++)
                        {
                            m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
                            m_elmt      = m_fields[0]->GetExp(m_HBCdata[j].m_globalElmtID);
                            m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();         
                            m_HBCdata[j].m_physOffset = m_fields[0]->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
                            m_HBCdata[j].m_bndElmtOffset = i+k*exp_size_per_plane;       
                            m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
                            m_HBCdata[j].m_bndryElmtID = n;
                            m_HBCdata[j].m_coeffOffset = coeff_count;
                            coeff_count += m_elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
                            
                            if(m_SingleMode)
                            {
                                m_wavenumber[j]      = -2*M_PI/m_LhomZ;       
                                m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                            }
                            else if(m_HalfMode || m_MultipleModes)
                            {
                                m_wavenumber[j]      = 2*M_PI/m_LhomZ;       
                                m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                            }
                            else
                            {
                                m_wavenumber[j]     = 2*M_PI*sign*(NekDouble(K))/m_LhomZ; 
                                m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                            }
                            
                            int assElmtID;
                            if(k%2==0)
                            {
                                if(m_HalfMode)
                                {
                                    assElmtID = m_HBCdata[j].m_globalElmtID;
                                    
                                }
                                else
                                {
                                    assElmtID = m_HBCdata[j].m_globalElmtID + num_elm_per_plane;
                                }
                            }
                            else 
                            {
                                assElmtID = m_HBCdata[j].m_globalElmtID - num_elm_per_plane;
                            }
                            
                            m_HBCdata[j].m_assPhysOffset = m_fields[0]->GetPhys_Offset(assElmtID);
                            
                            j = j+1;
                        }
                    }
                    else // setting if just standard BC no High order
                    {
                        cnt += exp_size_per_plane;
                    }
                }
                sign = -1.0*sign;
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        else if(m_HomogeneousType == eHomogeneous2D)
        {
            m_HBCdata = Array<OneD, HBCInfo>(HOPBCnumber);
            
            int cnt = 0;
            int exp_size, exp_size_per_line;
            int j=0;
            
            for(int k1 = 0; k1 < m_npointsZ; k1++)
            {
                for(int k2 = 0; k2 < m_npointsY; k2++)
                {
                    for(int n = 0 ; n < PBndConds.num_elements(); ++n)
                    {
                        SpatialDomains::BndUserDefinedType type = PBndConds[n]->GetUserDefined();
                        
                        exp_size = PBndExp[n]->GetExpSize();
                        
                        exp_size_per_line = exp_size/(m_npointsZ*m_npointsY);
                        
                        if(type == SpatialDomains::eHigh)
                        {
                            for(int i = 0; i < exp_size_per_line; ++i,cnt++)
                            {
                                // find element and edge of this expansion. 
                                // calculate curl x curl v;
                                m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];
                                m_elmt      = m_fields[0]->GetExp(m_HBCdata[j].m_globalElmtID);
                                m_HBCdata[j].m_ptsInElmt = m_elmt->GetTotPoints();
                                m_HBCdata[j].m_physOffset = m_fields[0]->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
                                m_HBCdata[j].m_bndElmtOffset = i+(k1*m_npointsY+k2)*exp_size_per_line;
                                m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];                
                                m_HBCdata[j].m_bndryElmtID = n;
                            }
                        }
                        else
                        {
                            cnt += exp_size_per_line;
                        }
                    }
                }
            }
        }
        /////////////////////////////////////////////////////////////////////////
    }
    
} //end of namespace

/**
 * $Log: VelocityCorrectionScheme.cpp,v $
 * Revision 1.5  2010/01/28 15:17:59  abolis
 * Time-Dependent boundary conditions
 *
 * Revision 1.4  2009/12/09 13:16:58  abolis
 * Update related to regression test
 *
 * Revision 1.3  2009/09/10 10:42:49  pvos
 * Modification to bind object pointer rather than object itself to time-integration functions.
 * (Prevents unwanted copy-constructor calls)
 *
 * Revision 1.2  2009/09/06 22:31:16  sherwin
 * First working version of Navier-Stokes solver and input files
 *
 * Revision 1.1  2009/03/14 16:43:21  sherwin
 * First non-working version
 *
 *
 **/
