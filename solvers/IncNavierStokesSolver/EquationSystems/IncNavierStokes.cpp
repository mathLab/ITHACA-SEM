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

#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Filters/Filter.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <complex>
#include <iostream>
#include <fstream> 
#include <sstream>


namespace Nektar
{

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    IncNavierStokes::IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession):
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
                if(boost::iequals(velids[i], var))
                {
                    m_velocity[i] = j;
                    break;
                }
                
                ASSERTL0(j != numfields, "Failed to find field: " + var);
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
                            SpatialDomains::eI ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eHighOutflow ||
                        m_fields[0]->GetBndConditions()[n]->GetUserDefined() ==
                            SpatialDomains::eWomersley,
                        "Unknown USERDEFINEDTYPE boundary condition");
                }
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
        
        m_session->LoadParameter("Kinvis", m_kinvis);
        
        // Default advection type per solver
        std::string vConvectiveType;
        switch(m_equationType)
        {
            case eUnsteadyStokes:
                vConvectiveType = "NoAdvection";
                break;
            case eUnsteadyNavierStokes:
            case eSteadyNavierStokes:
                vConvectiveType = "Convective";
                break;
            case eUnsteadyLinearisedNS:
                vConvectiveType = "Linearised";
                break;
            default:
                break;
        }

        // Check if advection type overridden
        if (m_session->DefinesTag("AdvectiveType") && m_equationType != eUnsteadyStokes)
        {
            vConvectiveType = m_session->GetTag("AdvectiveType");
        }

        // Initialise advection
        m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        
        // Forcing terms
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               v_GetForceDimension());

        // check to see if any Robin boundary conditions and if so set
        // up m_field to boundary condition maps;
        m_fieldsBCToElmtID  = Array<OneD, Array<OneD, int> >(numfields);
        m_fieldsBCToTraceID = Array<OneD, Array<OneD, int> >(numfields);
        m_fieldsRadiationFactor  = Array<OneD, Array<OneD, NekDouble> > (numfields);
        
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
		// Set up mapping for womersley BC - need to put this in a seperate loop for all the womersley info
		if(BndConds[n]->GetUserDefined() ==SpatialDomains::eWomersley){
			 m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);
		}
            }
        }

	// Set up maping for womersley BC - and load variabls
//        for (int i = 0; i < m_fields.num_elements(); ++i){
//            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n){
//		if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() 
//			==SpatialDomains::eWomersley){
//
//			NekDouble kt,T,R,n0,n1,n2;
//			int M;	
//		
//			m_session->LoadParameter("Period",T);
////			m_session->LoadParameter("Alpha",alpha);
//			m_session->LoadParameter("Radius",R);
//			m_session->LoadParameter("Modes",M);
//			m_session->LoadParameter("n0",n0);
//			m_session->LoadParameter("n1",n1);
//			m_session->LoadParameter("n2",n2);
//
//			NekDouble value;
//	        	Array<OneD, NekDouble> vel_r;
//	        	Array<OneD, NekDouble> vel_i;
//
//			std::ifstream file("fourier_womersley.txt");
//			std::string line;
//
//			int linenum=0;
//			while(std::getline(file,line)){
//      			        std::istringstream iss(line);
//				if (linenum == 0){
//			        	while(iss << value){
//					      vel_i.push_back (value);
//				        }
// 				}
//				else if (linenum == 1){
//					while(iss << value){
//					      vel_r.push_back (value);
//				        }
//
//				}
//				linenum += 1;
//				
//			}
//		}
//	    }
//	}

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"] = boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<std::string>(m_timestep);
		
        // creation of the extrapolation object
        if(m_equationType == eUnsteadyNavierStokes)
        {
            std::string vExtrapolation = "Standard";

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = m_session->GetSolverInfo("Extrapolation");
            }
                        
            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation, 
                m_session,
                m_fields,
		m_pressure,
                m_velocity,
                m_advObject);
        }
    }

    /**
     * Destructor
     */
    IncNavierStokes::~IncNavierStokes(void)
    {
    }

    
	/**
	 *
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
            ASSERTL0(wk.num_elements() >= nqtot*VelDim,
                     "Workspace is not sufficient");
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
        int i, n;
        std::string varName;
        int nvariables = m_fields.num_elements();
        
        for (i = 0; i < nvariables; ++i)
        {
            for(n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            { 
                if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() ==
                   SpatialDomains::eTimeDependent)
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
		 else if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() ==
                   SpatialDomains::eWomersley)  
                {
                    SetWomersleyBoundary(i,n);
                }

            }

            // Set Radiation conditions if required
            SetRadiationBoundaryForcing(i);
	    
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
            else if(type == SpatialDomains::eNoUserDefined ||
                    type == SpatialDomains::eWall_Forces ||
                    type == SpatialDomains::eTimeDependent ||
                    type == SpatialDomains::eHigh ||
                    type == SpatialDomains::eWomersley ||
                    type == SpatialDomains::eHighOutflow)
            {
                cnt += BndExp[n]->GetExpSize();
            }
            else
            {
                ASSERTL0(false,"Unknown USERDEFINEDTYPE in pressure boundary condition");
            }
        }
    }




    void IncNavierStokes::SetWomersleyBoundary(int fieldid,int bndid)
    {
	std::complex<NekDouble> za, zar, zJ0, zJ0r, zq, zvel, zJ0rJ0;
 	int  i,j,k;

	NekDouble kt,T,R,n0,n1,n2,x0,y0,z0;
	int M;	
		
	m_session->LoadParameter("Period",T);
//	m_session->LoadParameter("Alpha",alpha);
	m_session->LoadParameter("Radius",R);
	m_session->LoadParameter("Modes",M);
	m_session->LoadParameter("n0",n0);
	m_session->LoadParameter("n1",n1);
	m_session->LoadParameter("n2",n2);
	m_session->LoadParameter("x0",x0);
	m_session->LoadParameter("x1",y0);
	m_session->LoadParameter("x2",z0);
	
	// Womersley Number
	NekDouble alpha = R*sqrt(2*M_PI/T/m_kinvis);

	NekDouble normals[] = {n0,n1,n2};
	
	NekDouble vel_i[] = {.0000000000000000,	-0.0114652024383338,	0.0851029213457966,	-0.0063819332723328,	-0.0445916909877261,	0.0237233726269646,	0.0008464594812739,	0.0000176601880968,	-0.0026800533453590,	-0.0024003801647357,	0.0032485079266381,	-0.0001679460999985,	-0.0015374083917995,	-0.0007321330812402};
	NekDouble vel_r[] ={0.3789045336112559,	-0.1298052177588835,	-0.0005926433502513,	0.0322744583705603,	0.0051737128568664,	-0.0220534484187635,	0.0045170130867191,	0.0037315438343909,	0.0026464563011693,	-0.0050049685948089,	0.0002078141619409,	0.0021292792048836,	0.0003343912995364,	-0.0014819344644493};

	
	
	NekDouble r;

	std::complex<NekDouble> z1 (1.0,0.0);
	std::complex<NekDouble> zi (0.0,1.0);
	std::complex<NekDouble> z;

        
        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;
        
        
        BndConds = m_fields[fieldid]->GetBndConditions();
        BndExp   = m_fields[fieldid]->GetBndCondExpansions();
	
	StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansionSharedPtr bc;
        
        int cnt=0;
        int elmtid,offset, boundary,nfq;
	
	Array<OneD, NekDouble> Bvals,w;	

//	int npoints = BndExp[bndid]->GetNpoints();
//
//  	Array<OneD, NekDouble> x0(npoints,0.0);
//    	Array<OneD, NekDouble> x1(npoints,0.0);
//      Array<OneD, NekDouble> x2(npoints,0.0);
//	Array<OneD, NekDouble> zeros(npoints,0.0);
//	Array<OneD, NekDouble> w(npoints,0.0);
//
//        BndExp[bndid]->GetCoords(x0,x1,x2);

	//Loop over all expansions	
	for(i = 0; i < BndExp[bndid]->GetExpSize(); ++i,cnt++){
		// Get element id and offset
		elmtid = m_fieldsBCToElmtID[fieldid][cnt];
                elmt   = m_fields[fieldid]->GetExp(elmtid); 
                offset = m_fields[fieldid]->GetPhys_Offset(elmtid); 

		// Get Boundary and trace expansion
		bc = BndExp[bndid]->GetExp(i);
		boundary = m_fieldsBCToTraceID[fieldid][cnt];
		
		nfq=bc->GetTotPoints();
		w = m_fields[fieldid]->UpdatePhys() + offset;
		
	        Array<OneD, NekDouble> x(nfq,0.0);
    		Array<OneD, NekDouble> y(nfq,0.0);
        	Array<OneD, NekDouble> z(nfq,0.0);
		Array<OneD, NekDouble> wbc(nfq,0.0);
		bc->GetCoords(x,y,z);

		// Add edge values (trace) into the wbc 
		elmt->GetTracePhysVals(boundary,bc,w,wbc);
		//Compute womersley solution
	        for (j=0;j<nfq;j++){
			r = sqrt((x[j]-x0)*(x[j]-x0) + (y[j]-y0)*(y[j]-y0) + (z[j]-z0)*(z[j]-z0))/R;

			wbc[j] = vel_r[0]*(1. - r*r); // Compute Poiseulle Flow
			for (k=1; k<M; k++){
				kt = 2.0*M_PI*k*m_time/T;
				za = alpha*sqrt(k)/sqrt(2.0)*std::complex<NekDouble>(-1.0,1.0);
				zar = r*za;
				zJ0 = CompBessel(0,za);
				zJ0r = CompBessel(0,zar);
				zJ0rJ0 = zJ0r/zJ0;
				zq = std::exp(zi*kt)*std::complex<NekDouble>(vel_r[k],vel_i[k]);
//				zq = std::complex<double>(vel_r[k],vel_i[k])*std::complex<double>(cos(kt),sin(kt));
				zvel = zq*(z1-zJ0rJ0);
				wbc[j] = wbc[j]+zvel.real();
			}
		}
		// Multiply w by normal to get u,v,w component of velocity
		Vmath::Smul(nfq,normals[fieldid],wbc,1,wbc,1);

		Bvals = BndExp[bndid]->UpdateCoeffs()+BndExp[bndid]->GetCoeff_Offset(i);
		// Push back to Coeff space	
		bc->FwdTrans(wbc,Bvals);
	}

    }

/* Computes the Complex Bessel function of 1st kind integer order using series rep. - taken from numberical recipies in C.
*/
    std::complex<NekDouble> IncNavierStokes::CompBessel(int n, std::complex<NekDouble> y)
    {
	std::complex<NekDouble> z (1.0,0.0);
	std::complex<NekDouble> zbes (1.0,0.0);
	std::complex<NekDouble> zarg;
	NekDouble tol = 1e-15;
	int maxit = 10000;
	int i = 1;

	zarg = Rmul(-0.25,Cmul(y,y));

	while (Cabs(z) > tol && i <= maxit){
		z = Cmul(z,Rmul(1.0/i/(i+n),zarg));
		if  (Cabs(z) <= tol) break;
		zbes = zbes + z;
		i++;
	}

	zarg = Rmul(0.5,y);
	for (i=1;i<=n;i++){
		zbes = Cmul(zbes,zarg);
	}
	return zbes;

    }
/* Complex arithmatic - only used to ensure correct Wom solution - will remove and use c++/boost math*/
    std::complex<NekDouble> IncNavierStokes::Csub(std::complex<NekDouble> a, std::complex<NekDouble> b){
	std::complex<NekDouble> c;
	c.real() = a.real() - b.real();
	c.imag() = a.imag() - b.imag();
	return c;
    }
    std::complex<NekDouble> IncNavierStokes::Cmul(std::complex<NekDouble> a, std::complex<NekDouble> b){
	std::complex<NekDouble> c;
	c.real() = a.real()*b.real() - a.imag()*b.imag();
	c.imag() = a.imag()*b.real() + a.real()*b.imag();
	return c;
    }
    std::complex<NekDouble> IncNavierStokes::Cdiv(std::complex<NekDouble> a, std::complex<NekDouble> b){
	std::complex<NekDouble> c;
	NekDouble den,r;
	if (fabs(b.real())>=fabs(b.imag())){
	r = b.imag()/b.real();
	den = b.real() + r*b.imag();
	c.real() = (a.real() + r*a.imag())/den;
	c.imag() = (a.imag() - r*a.real())/den;
	}
	else{
	r = b.real()/b.imag();
	den = b.imag() + r*b.real();
	c.real() = (a.real()*r + a.imag())/den;
	c.imag() = (a.imag()*r - a.real())/den;

	}
	return c;
    }
    std::complex<NekDouble> IncNavierStokes::Rmul(NekDouble x, std::complex<NekDouble> a){
	std::complex<NekDouble> c;
	c.real() = x*a.real();
	c.imag() = x*a.imag();

	return c;
    }
    NekDouble IncNavierStokes::Cabs(std::complex<NekDouble> z){
	NekDouble x,y,ans,temp;
	
	x = fabs(z.real());
	y = fabs(z.imag());
	if (x==0.0){
	    ans=y;
    	}
	else if (y==0.0){
	    ans=x;
	}
	else if(x > y){	
	    temp = x/y;
	    ans = x*sqrt(1.0+temp*temp);
	}
	else{
	    temp = x/y;
	    ans = y*sqrt(1.0+temp*temp);
	}
	return ans;
     }

     /**
     * Add an additional forcing term programmatically.
     */
    void IncNavierStokes::AddForcing(const SolverUtils::ForcingSharedPtr& pForce)
    {
        m_forcing.push_back(pForce);
    }


    /**
     * Decide if at a steady state if the discrerte L2 sum of the
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
        int n_vel     = m_velocity.num_elements();
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
            velfields = Array<OneD, Array<OneD, NekDouble> >(n_vel);

            for(int i = 0; i < n_vel; ++i)
            {
                velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
            }        
        }

        stdVelocity = m_extrapolation->GetMaxStdVelocity(velfields);
        
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
        
        // express element id with respect to plane
        if(m_HomogeneousType == eHomogeneous1D)
        {
            elmtid = elmtid%m_fields[0]->GetPlane(0)->GetExpSize();
        }
        return CFL;
    }


    /**
     * Perform the extrapolation.
     */
    bool IncNavierStokes::v_PreIntegrate(int step)
    {
        m_extrapolation->SubStepSaveFields(step);
        m_extrapolation->SubStepAdvance(m_intSoln,step,m_time);
        return false;
    }


    /**
     * Estimate CFL and perform steady-state check
     */
    bool IncNavierStokes::v_PostIntegrate(int step)
    {
        if(m_cflsteps && !((step+1)%m_cflsteps))
        {
            int elmtid;
            NekDouble cfl = GetCFLEstimate(elmtid);

            if(m_comm->GetRank() == 0)
            {
                cout << "CFL (zero plane): "<< cfl << " (in elmt "
                     << elmtid << ")" << endl;
            }
        }

        if(m_steadyStateSteps && step && (!((step+1)%m_steadyStateSteps)))
        {
            if(CalcSteadyState() == true)
            {
                cout << "Reached Steady State to tolerance "
                     << m_steadyStateTol << endl;
                return true;
            }
        }

        return false;
    }
} //end of namespace

