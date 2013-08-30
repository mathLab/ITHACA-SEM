///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.cpp
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
// Description: Incompressible Navier Stokes class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <cstdio>
#include <cstdlib>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Filters/Filter.h>
#include <iomanip>

namespace Nektar
{

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    IncNavierStokes::IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession):
        //EquationSystem(pSession),
        UnsteadySystem(pSession),
        m_subSteppingScheme(false),
        m_SmoothAdvection(false),
        m_steadyStateSteps(0)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        
        int i,j;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};
        
        // Set up Velocity field to point to the first m_expdim of m_fields; 
        m_velocity = Array<OneD,int>(m_spacedim);
        
        for(i = 0; i < m_spacedim; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(NoCaseStringCompare(velids[i],var) == 0)
                {
                    m_velocity[i] = j;
                    break;
                }
                
                if(j == numfields)
                {
                    std::string error = "Failed to find field: " + var; 
                    ASSERTL0(false,error.c_str());
                }
            }
        }
        
        // Set up equation type enum using kEquationTypeStr
        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("EQTYPE",kEquationTypeStr[i],match,false);
            if(match)
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }
        ASSERTL0(i != eEquationTypeSize,"EQTYPE not found in SOLVERINFO section");
        
        // This probably should to into specific implementations 
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eSteadyStokes: 
        case eSteadyOseen: 
        case eSteadyNavierStokes:
        case eSteadyLinearisedNS: 
            break;
        case eUnsteadyNavierStokes:
        case eUnsteadyStokes:
            {
                m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
                m_session->LoadParameter("IO_CFLSteps", m_cflsteps, 0);
                m_session->LoadParameter("SteadyStateSteps", m_steadyStateSteps, 0);
                m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 1e-6);
            
				
                // check to see if any user defined boundary condition is
                // indeed implemented
                
                for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
                {    
                    // Time Dependent Boundary Condition (if no user
                    // defined then this is empty)
                    ASSERTL0 (
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eNoUserDefined ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eWall_Forces ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eTimeDependent ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eRadiation ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eI,
                        "Unknown USERDEFINEDTYPE boundary condition");
                }
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
        
        m_session->LoadParameter("Kinvis", m_kinvis);
        
        if (m_equationType == eUnsteadyNavierStokes || m_equationType == eSteadyNavierStokes)
        {
            std::string vConvectiveType = "Convective";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                vConvectiveType = m_session->GetTag("AdvectiveType");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
    
        if (m_equationType == eUnsteadyLinearisedNS)// || m_equationType == eSteadyNavierStokes)
        {
            std::string vConvectiveType = "Linearised";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                //vConvectiveType = m_session->GetTag("Linearised");
                vConvectiveType = m_session->GetTag("AdvectiveType");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
    
        if(m_equationType == eUnsteadyStokes)
        {
            std::string vConvectiveType = "NoAdvection";
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
        
        // check to see if any Robin boundary conditions and if so set
        // up m_field to boundary condition maps;
        m_fieldsBCToElmtID  = Array<OneD, Array<OneD, int> >(m_fields.num_elements());
        m_fieldsBCToTraceID = Array<OneD, Array<OneD, int> >(m_fields.num_elements());
        m_fieldsRadiationFactor  = Array<OneD, Array<OneD, NekDouble> > (m_fields.num_elements());
        
        for (i = 0; i < m_fields.num_elements(); ++i)
        {
            bool Set = false;

            Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
            Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
            int radpts = 0;
            
            BndConds = m_fields[i]->GetBndConditions();
            BndExp   = m_fields[i]->GetBndCondExpansions();
            for(int n = 0; n < BndConds.num_elements(); ++n)
            {    
                if(BndConds[n]->GetUserDefined() == SpatialDomains::eRadiation)
                {
                    ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SpatialDomains::eRobin,
                             "Radiation boundary condition must be of type Robin <R>");
                    
                    if(Set == false)
                    {
                        m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);
                        Set = true;
                    }
                    radpts += BndExp[n]->GetTotPoints();
                }
            }

            m_fieldsRadiationFactor[i] = Array<OneD, NekDouble>(radpts);

            radpts = 0; // reset to use as a counter

            for(int n = 0; n < BndConds.num_elements(); ++n)
            {    
                if(BndConds[n]->GetUserDefined() == SpatialDomains::eRadiation)
                {
                    
                    int npoints    = BndExp[n]->GetNpoints();
                    Array<OneD, NekDouble> x0(npoints,0.0);
                    Array<OneD, NekDouble> x1(npoints,0.0);
                    Array<OneD, NekDouble> x2(npoints,0.0);
                    Array<OneD, NekDouble> tmpArray;

                    BndExp[n]->GetCoords(x0,x1,x2);
                    
                    LibUtilities::Equation coeff = 
                        boost::static_pointer_cast<
                    SpatialDomains::RobinBoundaryCondition
                        >(BndConds[n])->m_robinPrimitiveCoeff;
                    
                    coeff.Evaluate(x0,x1,x2,m_time, 
                                   tmpArray = m_fieldsRadiationFactor[i]+ radpts);
                    //Vmath::Neg(npoints,tmpArray = m_fieldsRadiationFactor[i]+ radpts,1);
                    radpts += npoints;
                }
            }
        }

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"] = boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<std::string>(m_timestep);
		
		// creation of the extrapolation object
		if(m_equationType == eUnsteadyNavierStokes)
		{
			std::string vExtrapolation = "StandardExtrapolate";
			if (m_session->DefinesTag("Extrapolation"))
            {
                //vConvectiveType = m_session->GetTag("Linearised");
                vExtrapolation = m_session->GetTag("Extrapolation");
            }
			
			m_extrapolation = GetExtrapolateFactory().CreateInstance(vExtrapolation, 
																	 m_session,
																	 m_fields,
																	 m_velocity);
		}
		
    }

	/**
	 * Distructor
	 */
    IncNavierStokes::~IncNavierStokes(void)
    {
    }
    
	/**
	 * Advance in time - time-stepping loop
	 */
    void IncNavierStokes::AdvanceInTime(int nsteps)
    {
        int i,n;
        int n_fields = m_fields.num_elements();
        static int nchk = 0;
        
        Timer timer;

        if(m_HomogeneousType == eHomogeneous1D)
        {
            for(i = 0; i < n_fields; ++i)
            {
                m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
                m_fields[i]->SetWaveSpace(true);
                m_fields[i]->SetPhysState(false);
            }
        }
    
        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >  fields(m_nConvectiveFields);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields[i]  = m_fields[i]->UpdatePhys();
        }
        
        // Initialise NS solver which is set up to use a GLM method
        // with calls to EvaluateAdvection_SetPressureBCs and
        // SolveUnsteadyStokesSystem
        m_integrationSoln = m_integrationScheme->InitializeScheme(
                                m_timestep, fields, m_time, m_integrationOps);

        std::vector<SolverUtils::FilterSharedPtr>::iterator x;
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Initialise(m_fields, m_time);
        }
        
        //Time advance
        for(n = 0; n < nsteps; ++n)
        {

            timer.Start();

			m_extrapolation->SubStepSaveFields(n);
			m_extrapolation->SubStepAdvance(n);
            
            fields = m_integrationScheme->TimeIntegrate(n, m_timestep, m_integrationSoln, m_integrationOps);
            
            m_time += m_timestep;
            
            timer.Stop();

            // Write out current time step
            if(m_infosteps && !((n+1)%m_infosteps) && m_comm->GetRank() == 0)
            {
                cout << "Step: " << n+1 << "  Time: " << 
                    m_time << " CPU-Time: " << timer.TimePerTest(1) << " s" << endl;
            }

            if(m_cflsteps && !((n+1)%m_cflsteps))
            {
                int elmtid;
                NekDouble cfl = GetCFLEstimate(elmtid);

                if(m_comm->GetRank() == 0)
                {
                    cout << "CFL (zero plane): "<< cfl << " (in elmt " << elmtid << ")" << endl;
                }
            }
            
            // dump data in m_fields->m_coeffs to file. 
            if(m_checksteps &&(!((n+1)%m_checksteps)))
            {
                //WriteCheckpoint_Ouptput();
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    for(i = 0; i< n_fields; i++)
                    {
                        m_fields[i]->SetWaveSpace(false);
                        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(true);
                    }
                    nchk++;
                    Checkpoint_Output(nchk);
                    for(i = 0; i< n_fields; i++)
                    {
                        m_fields[i]->SetWaveSpace(true);
                        m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(false);
                    }
                }
                else
                {
                    nchk++;
                    Checkpoint_Output(nchk);
                }
            }
            
            
            if(m_steadyStateSteps && n && (!((n+1)%m_steadyStateSteps)))
            {
                if(CalcSteadyState() == true)
                {
                    cout << "Reached Steady State to tolerance " << m_steadyStateTol << endl; 
                    break;
                }
            }

            // Transform data into coefficient space
            if (m_filters.size() > 0)
            {
                for (i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_fields[i]->SetPhys(fields[i]);
                    m_fields[i]->SetPhysState(true);
                }
            }
            
            std::vector<SolverUtils::FilterSharedPtr>::iterator x;
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Update(m_fields, m_time);
            }

        }
    
        if(m_HomogeneousType == eHomogeneous1D)
        {
            for(i = 0 ; i< n_fields ; i++)
            {
                m_fields[i]->SetWaveSpace(false);
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
            }
        }
        else 
        {
            for(i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[i]->SetPhys(fields[i]);
                m_fields[i]->SetPhysState(true);
            }
        }
            
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Finalise(m_fields, m_time);
        }
    }
    
	/**
	 * Get Flux
	 */
    void IncNavierStokes::v_GetFluxVector(const int i, 
                                          Array<OneD, Array<OneD, NekDouble> > &physfield,
                                            Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(), physfield[i], 1, m_fields[m_velocity[j]]->GetPhys(), 1, flux[j], 1);
        }
    }

    /**
	 * Calcualate numerical fluxes
	 */
    void IncNavierStokes::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                          Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        /// Counter variable
        int i;

        /// Number of trace points
        int nTracePts   = GetTraceNpoints();
        
        /// Number of spatial dimensions
        int nDimensions = m_spacedim;

        /// Forward state array
        Array<OneD, NekDouble> Fwd(2*nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;
        
        /// Normal velocity array
        Array<OneD, NekDouble> Vn (nTracePts, 0.0);
        
        // Extract velocity field along the trace space and multiply by trace normals
        for(i = 0; i < nDimensions; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_fields[m_velocity[i]]->GetPhys(), Fwd);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        }

        /// Compute the numerical fluxes at the trace points
        for(i = 0; i < numflux.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            /// Upwind between elements
            m_fields[i]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux[i]);

            /// Calculate the numerical fluxes multipling Fwd or Bwd
            /// by the normal advection velocity
            Vmath::Vmul(nTracePts, numflux[i], 1, Vn, 1, numflux[i], 1);
        }
    }

	/**
	 * Evaluation -N(V) for all fields except pressure using m_velocity
	 */
    void IncNavierStokes::EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                 Array<OneD, NekDouble> &wk)
    {
        int i;
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, NekDouble > Deriv;
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]]; 
        }
        // Set up Derivative work space; 
        if(wk.num_elements())
        {
            ASSERTL0(wk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }

        m_advObject->DoAdvection(m_fields,m_nConvectiveFields, 
                                 m_velocity,inarray,outarray,m_time,Deriv);
    }
    
    /**
	 * Time dependent boundary conditions updating
	 */
    void IncNavierStokes::SetBoundaryConditions(NekDouble time)
    {
        int  nvariables = m_fields.num_elements();
        
        for (int i = 0; i < nvariables; ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {    
                if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() ==
                   SpatialDomains::eTimeDependent)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }

            }

            SetRadiationBoundaryForcing(i); // Set Radiation conditions if required. 
        }
    }
    
    /**
	 * Probably should be pushed back into ContField? 
	 */
    void IncNavierStokes::SetRadiationBoundaryForcing(int fieldid)
    {
        int  i,n;
        
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
        
        
        BndConds = m_fields[fieldid]->GetBndConditions();
        BndExp   = m_fields[fieldid]->GetBndCondExpansions();
        
        StdRegions::StdExpansionSharedPtr   elmt;
        StdRegions::StdExpansionSharedPtr Bc;
        
        int cnt;
        int elmtid,nq,offset, boundary;
        Array<OneD, NekDouble> Bvals, U;
        int cnt1 = 0;


        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
        {            
            SpatialDomains::BndUserDefinedType type = BndConds[n]->GetUserDefined(); 
            
            if((BndConds[n]->GetBoundaryConditionType() == SpatialDomains::eRobin)&&(type == SpatialDomains::eRadiation))
            {
                for(i = 0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    elmtid = m_fieldsBCToElmtID[fieldid][cnt];
                    elmt   = m_fields[fieldid]->GetExp(elmtid);
                    offset = m_fields[fieldid]->GetPhys_Offset(elmtid);
                    
                    U = m_fields[fieldid]->UpdatePhys() + offset;
                    Bc = BndExp[n]->GetExp(i);
                    
                    boundary = m_fieldsBCToTraceID[fieldid][cnt];
                    
                    // Get edge values and put into ubc
                    nq = Bc->GetTotPoints();
                    Array<OneD, NekDouble> ubc(nq);
                    elmt->GetTracePhysVals(boundary,Bc,U,ubc);
                    
                    Vmath::Vmul(nq,&m_fieldsRadiationFactor[fieldid][cnt1 + 
                             BndExp[n]->GetPhys_Offset(i)],1,&ubc[0],1,&ubc[0],1);

                    Bvals = BndExp[n]->UpdateCoeffs()+BndExp[n]->GetCoeff_Offset(i);

                    Bc->IProductWRTBase(ubc,Bvals); 
                }
                cnt1 += BndExp[n]->GetTotPoints();
            }
            else if(type == SpatialDomains::eNoUserDefined || type == SpatialDomains::eWall_Forces || type == SpatialDomains::eTimeDependent || type == SpatialDomains::eHigh) 
            {
                cnt += BndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }
    }

    /** Decide if at a steady state if the discrerte L2 sum of the
	 * coefficients is the same as the previous step to within the
	 * tolerance m_steadyStateTol;
	 */
    bool IncNavierStokes::CalcSteadyState(void)
    {
        static NekDouble previousL2 = 0.0;
        bool returnval = false;
        
        NekDouble L2 = 0.0;
        
        // calculate L2 discrete summation 
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            L2 += Vmath::Dot(ncoeffs,m_fields[i]->GetCoeffs(),1,m_fields[i]->GetCoeffs(),1);
        }
        
        if(fabs(L2-previousL2) < ncoeffs*m_steadyStateTol)
        {
            returnval = true;
        }

        previousL2 = L2;

        return returnval;
    }
    
	/**
	 *
	 */
    Array<OneD, NekDouble> IncNavierStokes::GetElmtCFLVals(void)
    { 
        int n_element = m_fields[0]->GetExpSize(); 
        
        const Array<OneD, int> ExpOrder = GetNumExpModesPerExp();
        Array<OneD, int> ExpOrderList (n_element, ExpOrder);
        
        const NekDouble cLambda = 0.2; // Spencer book pag. 317
        
        Array<OneD, NekDouble> cfl        (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields; 
        
        if(m_HomogeneousType == eHomogeneous1D) // just do check on 2D info
        {
            velfields = Array<OneD, Array<OneD, NekDouble> >(2);

            for(int i = 0; i < 2; ++i)
            {
                velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
            }        
        }
        else
        {
            velfields = Array<OneD, Array<OneD, NekDouble> >(m_velocity.num_elements());

            for(int i = 0; i < m_velocity.num_elements(); ++i)
            {
                velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
            }        
        }
        stdVelocity = GetMaxStdVelocity(velfields);
        
        for(int el = 0; el < n_element; ++el)
        {
            cfl[el] =  m_timestep*(stdVelocity[el] * cLambda *
                                   (ExpOrder[el]-1) * (ExpOrder[el]-1));
        }
        
        return cfl;
    }
    
    /**
	 *
	 */
    NekDouble IncNavierStokes::GetCFLEstimate(int &elmtid)
    { 
        int n_element = m_fields[0]->GetExpSize(); 

        Array<OneD, NekDouble> cfl = GetElmtCFLVals();
        
        elmtid = Vmath::Imax(n_element,cfl,1);
        NekDouble CFL,CFL_loc;

        CFL = CFL_loc = cfl[elmtid];
        m_comm->AllReduce(CFL,LibUtilities::ReduceMax);

        // unshuffle elmt id if data is not stored in consecutive order. 
        elmtid = m_fields[0]->GetExp(elmtid)->GetGeom()->GetGlobalID();
        if(CFL != CFL_loc)
        {
            elmtid = -1;
        }

        m_comm->AllReduce(elmtid,LibUtilities::ReduceMax);
        
        if(m_HomogeneousType == eHomogeneous1D) // express element id with respect to plane
        {
            elmtid = elmtid%m_fields[0]->GetPlane(0)->GetExpSize();
        }
        return CFL;
    }
} //end of namespace

