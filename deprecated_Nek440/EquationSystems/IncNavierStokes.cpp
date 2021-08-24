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

#include "./IncNavierStokes.h"
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Filters/Filter.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>

#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>

#include <tinyxml.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

using namespace std;

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
        AdvectionSystem(pSession),
        m_SmoothAdvection(false),
        m_steadyStateSteps(0)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

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
        case eSteadyLinearisedNS:
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
        if (m_session->DefinesTag("AdvectiveType") && m_equationType != eUnsteadyStokes &&
            m_equationType != eSteadyLinearisedNS)
        {
            vConvectiveType = m_session->GetTag("AdvectiveType");
        }

        // Initialise advection
        m_advObject = SolverUtils::GetAdvectionFactory().CreateInstance(vConvectiveType, vConvectiveType);
        m_advObject->InitObject( m_session, m_fields);

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
                if(boost::iequals(BndConds[n]->GetUserDefined(),"Radiation"))
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
                if(boost::iequals(BndConds[n]->GetUserDefined(),"ZeroNormalComponent"))
                {
                    ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet,
                             "Zero Normal Component boundary condition option must be of type Dirichlet <D>");

                    if(Set == false)
                    {
                        m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);
                        Set = true;
                    }
                }
            }

            m_fieldsRadiationFactor[i] = Array<OneD, NekDouble>(radpts);

            radpts = 0; // reset to use as a counter

            for(int n = 0; n < BndConds.num_elements(); ++n)
            {
                if(boost::iequals(BndConds[n]->GetUserDefined(),"Radiation"))
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

        // Set up maping for womersley BC - and load variables
        for (int i = 0; i < m_fields.num_elements(); ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {
                if(boost::istarts_with(m_fields[i]->GetBndConditions()[n]->GetUserDefined(),"Womersley"))
                {

                    m_womersleyParams[n] = MemoryManager<WomersleyParams>::AllocateSharedPtr(m_spacedim);


#if 0
                    m_session->LoadParameter("Period",m_womersleyParams[n]->m_period);
                    m_session->LoadParameter("Radius",m_womersleyParams[n]->m_radius);

                    NekDouble n0,n1,n2;
                    m_session->LoadParameter("n0",n0);
                    m_session->LoadParameter("n1",n1);
                    m_session->LoadParameter("n2",n2);
                    m_womersleyParams[n]->m_axisnormal[0] = n0;
                    m_womersleyParams[n]->m_axisnormal[1] = n1;
                    m_womersleyParams[n]->m_axisnormal[2] = n2;

                    NekDouble x0,x1,x2;
                    m_session->LoadParameter("x0",x0);
                    m_session->LoadParameter("x1",x1);
                    m_session->LoadParameter("x2",x2);
                    m_womersleyParams[n]->m_axispoint[0] = x0;
                    m_womersleyParams[n]->m_axispoint[1] = x1;
                    m_womersleyParams[n]->m_axispoint[2] = x2;
#endif

                    // Read in fourier coeffs
                    SetUpWomersley(n,
                                   m_fields[i]->GetBndConditions()[n]->GetUserDefined());

                    m_fields[i]->GetBoundaryToElmtMap(m_fieldsBCToElmtID[i],m_fieldsBCToTraceID[i]);

                }
            }
        }

        // Set up Field Meta Data for output files
        m_fieldMetaDataMap["Kinvis"]   = boost::lexical_cast<std::string>(m_kinvis);
        m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<std::string>(m_timestep);
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
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        int i;
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);

        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]];
        }

        m_advObject->Advect(m_nConvectiveFields, m_fields,
                            velocity, inarray, outarray, m_time);
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
                if(m_fields[i]->GetBndConditions()[n]->IsTimeDependent())
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
                else if(boost::istarts_with(m_fields[i]->GetBndConditions()[n]->GetUserDefined(),"Womersley"))
                {
                    SetWomersleyBoundary(i,n);
                }
            }

            // Set Radiation conditions if required
            SetRadiationBoundaryForcing(i);
        }

        SetZeroNormalVelocity();
    }

    /**
     * Probably should be pushed back into ContField?
     */
    void IncNavierStokes::SetRadiationBoundaryForcing(int fieldid)
    {
        int  i,n;

        Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
        Array<OneD, MultiRegions::ExpListSharedPtr>                BndExp;

        BndConds = m_fields[fieldid]->GetBndConditions();
        BndExp   = m_fields[fieldid]->GetBndCondExpansions();

        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansionSharedPtr Bc;

        int cnt;
        int elmtid,nq,offset, boundary;
        Array<OneD, NekDouble> Bvals, U;
        int cnt1 = 0;

        for(cnt = n = 0; n < BndConds.num_elements(); ++n)
        {
            std::string type = BndConds[n]->GetUserDefined();

            if((BndConds[n]->GetBoundaryConditionType() == SpatialDomains::eRobin)&&(boost::iequals(type,"Radiation")))
            {
                for(i = 0; i < BndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    elmtid = m_fieldsBCToElmtID[m_velocity[fieldid]][cnt];
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
            else
            {
                cnt += BndExp[n]->GetExpSize();
            }
        }
    }


    void IncNavierStokes::SetZeroNormalVelocity()
    {
        // use static trip since cannot use UserDefinedTag for zero
        // velocity and have time dependent conditions
        static bool Setup  = false;

        if(Setup == true)
        {
            return;
        }
        Setup = true;

        int  i,n;

        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr > >
            BndConds(m_spacedim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
            BndExp(m_spacedim);


        for(i = 0; i < m_spacedim; ++i)
        {
            BndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
            BndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        StdRegions::StdExpansionSharedPtr elmt,Bc;

        int cnt;
        int elmtid,nq, boundary;

        Array<OneD, Array<OneD, NekDouble> > normals;
        Array<OneD, NekDouble> Bphys,Bcoeffs;

        int fldid = m_velocity[0];

        for(cnt = n = 0; n < BndConds[0].num_elements(); ++n)
        {
            if((BndConds[0][n]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)&& (boost::iequals(BndConds[0][n]->GetUserDefined(),"ZeroNormalComponent")))
            {
                for(i = 0; i < BndExp[0][n]->GetExpSize(); ++i,cnt++)
                {
                    elmtid   = m_fieldsBCToElmtID[fldid][cnt];
                    elmt     = m_fields[0]->GetExp(elmtid);
                    boundary = m_fieldsBCToTraceID[fldid][cnt];

                    normals = elmt->GetSurfaceNormal(boundary);

                    nq = BndExp[0][n]->GetExp(i)->GetTotPoints();
                    Array<OneD, NekDouble> normvel(nq,0.0);

                    for(int k = 0; k < m_spacedim; ++k)
                    {
                        Bphys  = BndExp[k][n]->UpdatePhys()+
                            BndExp[k][n]->GetPhys_Offset(i);
                        Bc  = BndExp[k][n]->GetExp(i);
                        Vmath::Vvtvp(nq,normals[k],1,Bphys,1,normvel,1,
                                     normvel,1);
                    }

                    // negate normvel for next step
                    Vmath::Neg(nq,normvel,1);

                    for(int k = 0; k < m_spacedim; ++k)
                    {
                        Bphys  = BndExp[k][n]->UpdatePhys()+
                            BndExp[k][n]->GetPhys_Offset(i);
                        Bcoeffs = BndExp[k][n]->UpdateCoeffs()+
                            BndExp[k][n]->GetCoeff_Offset(i);
                        Bc  = BndExp[k][n]->GetExp(i);
                        Vmath::Vvtvp(nq,normvel,1,normals[k],1,Bphys,1,
                                     Bphys,1);
                        Bc->FwdTrans_BndConstrained(Bphys,Bcoeffs);
                    }
                }
            }
            else
            {
                cnt += BndExp[0][n]->GetExpSize();
            }
        }
    }


    /**
     *  Womersley boundary condition defintion
     */
    void IncNavierStokes::SetWomersleyBoundary(const int fldid, const int bndid)
    {
        ASSERTL1(m_womersleyParams.count(bndid) == 1, "Womersley parameters for this boundary have not been set up");

        WomersleyParamsSharedPtr WomParam = m_womersleyParams[bndid];
        std::complex<NekDouble> za, zar, zJ0, zJ0r, zq, zvel, zJ0rJ0;
        int  i,j,k;

        int M = WomParam->m_wom_vel_r.size();

        NekDouble R = WomParam->m_radius;
        NekDouble T = WomParam->m_period;

        Array<OneD, NekDouble > normals = WomParam->m_axisnormal;
        Array<OneD, NekDouble > x0      = WomParam->m_axispoint;

        // Womersley Number
        NekDouble alpha = R*sqrt(2*M_PI/T/m_kinvis);
        NekDouble r,kt;

        std::complex<NekDouble> z1 (1.0,0.0);
        std::complex<NekDouble> zi (0.0,1.0);
        std::complex<NekDouble> comp_conj (-1.0,1.0); //complex conjugate

        Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

        BndExp   = m_fields[fldid]->GetBndCondExpansions();

        StdRegions::StdExpansionSharedPtr elmt;
        StdRegions::StdExpansionSharedPtr bc;
        int cnt=0;
        int elmtid,offset, boundary,nfq;

        Array<OneD, NekDouble> Bvals,w;

        //Loop over all expansions
        for(i = 0; i < BndExp[bndid]->GetExpSize(); ++i,cnt++)
        {
            // Get element id and offset
            elmtid = m_fieldsBCToElmtID[fldid][cnt];
            elmt   = m_fields[fldid]->GetExp(elmtid);
            offset = m_fields[fldid]->GetPhys_Offset(elmtid);

            // Get Boundary and trace expansion
            bc = BndExp[bndid]->GetExp(i);
            boundary = m_fieldsBCToTraceID[fldid][cnt];

            nfq=bc->GetTotPoints();
            w = m_fields[fldid]->UpdatePhys() + offset;

            Array<OneD, NekDouble> x(nfq,0.0);
            Array<OneD, NekDouble> y(nfq,0.0);
            Array<OneD, NekDouble> z(nfq,0.0);
            Array<OneD, NekDouble> wbc(nfq,0.0);
            bc->GetCoords(x,y,z);

            // Add edge values (trace) into the wbc
            elmt->GetTracePhysVals(boundary,bc,w,wbc);

            //Compute womersley solution
            for (j=0;j<nfq;j++)
            {
                //NOTE: only need to calculate these two once, could
                //be stored or precomputed?
                r = sqrt((x[j]-x0[0])*(x[j]-x0[0]) +
                         (y[j]-x0[1])*(y[j]-x0[1]) +
                         (z[j]-x0[2])*(z[j]-x0[2]))/R;

                wbc[j] = WomParam->m_wom_vel_r[0]*(1. - r*r); // Compute Poiseulle Flow

                for (k=1; k<M; k++)
                {
                    kt = 2.0 * M_PI * k * m_time / T;
                    za = alpha * sqrt((NekDouble)k/2.0) * comp_conj;
                    zar = r * za;
                    zJ0  = Polylib::ImagBesselComp(0,za);
                    zJ0r = Polylib::ImagBesselComp(0,zar);
                    zJ0rJ0 = zJ0r / zJ0;
                    zq = std::exp(zi * kt) * std::complex<NekDouble>(
                                                   WomParam->m_wom_vel_r[k],
                                                   WomParam->m_wom_vel_i[k]);
                    zvel = zq * (z1 - zJ0rJ0);
                    wbc[j] = wbc[j] + zvel.real();
                }
            }

            // Multiply w by normal to get u,v,w component of velocity
            Vmath::Smul(nfq,normals[fldid],wbc,1,wbc,1);

            Bvals = BndExp[bndid]->UpdateCoeffs()+
                    BndExp[bndid]->GetCoeff_Offset(i);
            // Push back to Coeff space
            bc->FwdTrans(wbc,Bvals);
        }
    }


    void IncNavierStokes::SetUpWomersley(const int bndid, std::string womStr)
    {
        std::string::size_type indxBeg = womStr.find_first_of(':') + 1;
        string filename = womStr.substr(indxBeg,string::npos);

        std::complex<NekDouble> coef;

#if 1
        TiXmlDocument doc(filename);

        bool loadOkay = doc.LoadFile();
        ASSERTL0(loadOkay,(std::string("Failed to load file: ") +
                           filename).c_str());

        TiXmlHandle docHandle(&doc);

        int err;    /// Error value returned by TinyXML.

        TiXmlElement *nektar = doc.FirstChildElement("NEKTAR");
        ASSERTL0(nektar, "Unable to find NEKTAR tag in file.");

        TiXmlElement *wombc = nektar->FirstChildElement("WOMERSLEYBC");
        ASSERTL0(wombc, "Unable to find WOMERSLEYBC tag in file.");

        // read womersley parameters
        TiXmlElement *womparam = wombc->FirstChildElement("WOMPARAMS");
        ASSERTL0(womparam, "Unable to find WOMPARAMS tag in file.");

        // Input coefficients
        TiXmlElement *params = womparam->FirstChildElement("W");
        map<std::string,std::string> Wparams;

        // read parameter list
        while (params)
        {

            std::string propstr;
            propstr = params->Attribute("PROPERTY");

            ASSERTL0(!propstr.empty(),"Failed to read PROPERTY value Womersley BC Parameter");


            std::string valstr;
            valstr = params->Attribute("VALUE");

            ASSERTL0(!valstr.empty(),"Failed to read VALUE value Womersley BC Parameter");

            std::transform(propstr.begin(),propstr.end(),propstr.begin(),
                           ::toupper);
            Wparams[propstr] = valstr;

            params = params->NextSiblingElement("W");
        }

        // Read parameters

        ASSERTL0(Wparams.count("RADIUS") == 1,
          "Failed to find Radius parameter in Womersley boundary conditions");
        std::vector<NekDouble> rad;
        ParseUtils::GenerateUnOrderedVector(
                                         Wparams["RADIUS"].c_str(),rad);
        m_womersleyParams[bndid]->m_radius = rad[0];

        ASSERTL0(Wparams.count("PERIOD") == 1,
          "Failed to find period parameter in Womersley boundary conditions");
        std::vector<NekDouble> period;
        ParseUtils::GenerateUnOrderedVector(
                                         Wparams["PERIOD"].c_str(),period);
        m_womersleyParams[bndid]->m_period = period[0];


        ASSERTL0(Wparams.count("AXISNORMAL") == 1,
          "Failed to find axisnormal parameter in Womersley boundary conditions");
        std::vector<NekDouble> anorm;
        ParseUtils::GenerateUnOrderedVector(
                                         Wparams["AXISNORMAL"].c_str(),anorm);
        m_womersleyParams[bndid]->m_axisnormal[0] = anorm[0];
        m_womersleyParams[bndid]->m_axisnormal[1] = anorm[1];
        m_womersleyParams[bndid]->m_axisnormal[2] = anorm[2];


        ASSERTL0(Wparams.count("AXISPOINT") == 1,
          "Failed to find axispoint parameter in Womersley boundary conditions");
        std::vector<NekDouble> apt;
        ParseUtils::GenerateUnOrderedVector(
                                         Wparams["AXISPOINT"].c_str(),apt);
        m_womersleyParams[bndid]->m_axispoint[0] = apt[0];
        m_womersleyParams[bndid]->m_axispoint[1] = apt[1];
        m_womersleyParams[bndid]->m_axispoint[2] = apt[2];

        // Read Temporal Foruier Coefficients.

        // Find the FourierCoeff tag
        TiXmlElement *coeff = wombc->FirstChildElement("FOURIERCOEFFS");

        // Input coefficients
        TiXmlElement *fval = coeff->FirstChildElement("F");

        int indx;
        int nextFourierCoeff = -1;

        while (fval)
        {
            nextFourierCoeff++;

            TiXmlAttribute *fvalAttr = fval->FirstAttribute();
            std::string attrName(fvalAttr->Name());

            ASSERTL0(attrName == "ID", (std::string("Unknown attribute name: ") + attrName).c_str());

            err = fvalAttr->QueryIntValue(&indx);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute ID.");

            std::string coeffStr = fval->FirstChild()->ToText()->ValueStr();
            vector<NekDouble> coeffvals;
            bool parseGood = ParseUtils::GenerateUnOrderedVector(coeffStr.c_str(),
                                                               coeffvals);
            ASSERTL0(parseGood,(std::string("Problem reading value of fourier coefficient, ID=") + boost::lexical_cast<string>(indx)).c_str());
            ASSERTL1(coeffvals.size() == 2,(std::string("Have not read two entries of Fourier coefficicent from ID="+ boost::lexical_cast<string>(indx)).c_str()));
            m_womersleyParams[bndid]->m_wom_vel_r.push_back(coeffvals[0]);
            m_womersleyParams[bndid]->m_wom_vel_i.push_back(coeffvals[1]);

            fval = fval->NextSiblingElement("F");
        }

#else
        std::ifstream file(filename);
        std::string line;

        ASSERTL1(file.is_open(),(std::string("Missing file ") + filename).c_str());
        int count = 0;
        while(std::getline(file,line))
        {
            std::stringstream stream(line);
            while(stream>>coef)
            {
                m_womersleyParams[bndid]->m_wom_vel_r.push_back(coef.real());
                m_womersleyParams[bndid]->m_wom_vel_i.push_back(coef.imag());
                count++;
            }
        }
#endif
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
        SetBoundaryConditions(m_time+m_timestep);
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

