///////////////////////////////////////////////////////////////////////////////
//
// File APE.cpp
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
// Description: Acoustic perturbation equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <APESolver/EquationSystems/APE.h>

namespace Nektar
{
string APE::className = GetEquationSystemFactory().RegisterCreatorFunction(
            "APE", APE::create,
            "Acoustic perturbation equations in conservative variables.");

APE::APE(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
{
}

void APE::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    // Load constant incompressible density
    m_session->LoadParameter("Rho0", m_Rho0, 1.204);

    // Load isentropic coefficient, Ratio of specific heats
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    // Determine TimeIntegrationMethod to use.
    ASSERTL0(m_session->DefinesSolverInfo("TIMEINTEGRATIONMETHOD"),
             "No TIMEINTEGRATIONMETHOD defined in session.");
    int i;
    for (i = 0; i < (int)LibUtilities::SIZE_TimeIntegrationMethod; ++i)
    {
        bool match;
        m_session->MatchSolverInfo("TIMEINTEGRATIONMETHOD",
                                   LibUtilities::TimeIntegrationMethodMap[i], match, false);
        if (match)
        {
            m_timeIntMethod = (LibUtilities::TimeIntegrationMethod) i;
            break;
        }
    }
    ASSERTL0(i != (int) LibUtilities::SIZE_TimeIntegrationMethod,
             "Invalid time integration type.");


    // if discontinuous  determine numerical flux to use
    if (m_projectionType == MultiRegions::eDiscontinuous)
    {
        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        int i;
        for (i = 0; i < (int)SIZE_UpwindType; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("UPWINDTYPE",
                                       UpwindTypeMap[i], match, false);
            if (match)
            {
                m_upwindType = (UpwindType) i;
                break;
            }
        }
        ASSERTL0(i != (int) SIZE_UpwindType,
                 "Invalid upwind type.");
    }
    else
    {
        m_upwindType = (UpwindType) 0;
    }

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs     (&APE::DoOdeRhs,        this);
        m_ode.DefineProjection (&APE::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit APE not set up.");
    }
}

APE::~APE()
{
    
}


/**
     * Initialises the time integration scheme (as specified in the session
     * file), and perform the time integration.
     */
void APE::v_DoSolve()
{
    int i,n,nchk = 0;
    int nq = m_fields[0]->GetTotPoints();
    int nvariables = m_fields.num_elements();

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble> >   fields(nvariables);
    Array<OneD, Array<OneD, NekDouble> >   tmp(nvariables);

    for(i = 0; i < nvariables; ++i)
    {
        m_fields[i]->SetPhysState(false);
        fields[i]  = m_fields[i]->UpdatePhys();
    }


    /* Declare an array of TimeIntegrationSchemes For multi-stage
        methods, this array will have just one entry containing the
        actual multi-stage method...
        For multi-steps method, this can have multiple entries
        - the first scheme will used for the first timestep (this
        is an initialization scheme)
        - the second scheme will used for the second timestep
        (this is an initialization scheme)
        - ...
        - the last scheme will be used for all other time-steps
        (this will be the actual scheme)*/

    Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr> IntScheme;
    LibUtilities::TimeIntegrationSolutionSharedPtr u;
    int numMultiSteps;

    switch(m_timeIntMethod)
    {
    case LibUtilities::eForwardEuler:
    case LibUtilities::eClassicalRungeKutta4:
    {
        numMultiSteps = 1;

        IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

        LibUtilities::TimeIntegrationSchemeKey IntKey(m_timeIntMethod);
        IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey];

        u = IntScheme[0]->InitializeScheme(m_timestep,fields,m_time,m_ode);
        break;
    }
    case LibUtilities::eAdamsBashforthOrder2:
    {
        numMultiSteps = 2;

        IntScheme = Array<OneD, LibUtilities::TimeIntegrationSchemeSharedPtr>(numMultiSteps);

        // Used in the first time step to initalize the scheme
        LibUtilities::TimeIntegrationSchemeKey IntKey0(LibUtilities::eClassicalRungeKutta4);

        // Used for all other time steps
        LibUtilities::TimeIntegrationSchemeKey IntKey1(m_timeIntMethod);
        IntScheme[0] = LibUtilities::TimeIntegrationSchemeManager()[IntKey0];
        IntScheme[1] = LibUtilities::TimeIntegrationSchemeManager()[IntKey1];

        // Initialise the scheme for the actual time integration scheme
        u = IntScheme[1]->InitializeScheme(m_timestep,fields,m_time,m_ode);
        break;
    }
    default:
    {
        ASSERTL0(false,"populate switch statement for integration scheme");
    }
    }

    std::string outname = m_session->GetFilename() + ".his";
    std::ofstream hisFile (outname.c_str());

    // Perform integration in time.
    for(n = 0; n < m_steps; ++n)
    {
        // Integrate over timestep.
        if( n < numMultiSteps-1)
        {
            // Use initialisation schemes if time step is less than the
            // number of steps in the scheme.
            fields = IntScheme[n]->TimeIntegrate(m_timestep,u,m_ode);
        }
        else
        {
            fields = IntScheme[numMultiSteps-1]->TimeIntegrate(m_timestep,u,m_ode);
        }

        // Increment time.
        m_time += m_timestep;

        // Write out status information.
        if(!((n+1)%m_infosteps))
        {
            cout << "Steps: " << n+1 << "\t Time: " << m_time << "\t " << endl;
        }

        // Write out checkpoint files.
        if(n&&(!((n+1)%m_checksteps)))
        {

            // update m_fields
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(nq,fields[i],1,m_fields[i]->UpdatePhys(),1);
            }

            // go to primitive variables
            v_ConservativeToPrimitive();

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
            Checkpoint_Output(nchk++);

            v_PrimitiveToConservative();
        }
    }

    // At the end of the time integration, store final solution.
    // update m_fields
    for(i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq,fields[i],1,m_fields[i]->UpdatePhys(),1);
    }

    // to to primitive variables
    v_ConservativeToPrimitive();

    for(i = 0; i < nvariables; ++i)
    {
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
    }

}

/**
     *
     */
void APE::v_DoInitialise()
{
    SetInitialConditions();

    v_PrimitiveToConservative();
}

/**
     *
     */
void APE::v_GenerateSummary(SolverUtils::SummaryList& s)
{
    UnsteadySystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Upwind Type", UpwindTypeMap[m_upwindType]);
    SolverUtils::AddSummaryItem(s, "Advection", (m_explicitAdvection ? "explicit" : "implicit"));
    SolverUtils::AddSummaryItem(s, "Integration Type", LibUtilities::TimeIntegrationMethodMap[m_timeIntMethod]);
    SolverUtils::AddSummaryItem(s, "Time Step", m_timestep);
    SolverUtils::AddSummaryItem(s, "No. of Steps", m_steps);
    SolverUtils::AddSummaryItem(s, "Checkpoints (steps)", m_checksteps);
}

/**
     *
     */
void APE::v_NumericalFlux(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &numflux)
{
    ASSERTL0(false, "This function is not implemented for this equation.");
}


/**
     *
     */
void APE::v_NumericalFlux(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &numfluxX,
        Array<OneD, Array<OneD, NekDouble> > &numfluxY )
{
    switch(m_expdim)
    {
    case 1:
        ASSERTL0(false,"1D not implemented for Acoustic perturbation equations");
        break;
    case 2:
        NumericalFlux2D(physfield,numfluxX,numfluxY);
        break;
    case 3:
        ASSERTL0(false,"3D not implemented for Acoustic perturbation equations");
        break;
    default:
        ASSERTL0(false,"Illegal dimension");
    }
}


/**
     *
     */
void APE::v_GetFluxVector(const int i, const int j,
                           Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
{
    v_GetFluxVector(i, physfield, flux);
}


void APE::v_GetFluxVector(const int i,
                             Array<OneD, Array<OneD, NekDouble> > &physfield,
                             Array<OneD, Array<OneD, NekDouble> > &flux)
{
    InitialiseBaseFlowAnalytical(basefield, m_time);

    ASSERTL1(flux.num_elements() == basefield.num_elements() - 1,
             "Dimension of flux array and velocity array do not match");

    int nq = physfield[0].num_elements();
    NekDouble tmp0 = 0.0;
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    if (i == 0)
    {
        // F_{adv,p',j} = \gamma p_0 u'_j + p' \bar{u}_j
        for (int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Zero(nq, flux[j], 1);

            // construct \gamma p_0 u'_j term
            Vmath::Smul(nq, m_gamma, basefield[basefield.num_elements()-1], 1, tmp1, 1);
            Vmath::Vmul(nq, tmp1, 1, physfield[j+1], 1, tmp1, 1);

            // construct p' \bar{u}_j term
            Vmath::Vmul(nq, physfield[0], 1, basefield[j], 1, tmp2, 1);

            // add both terms
            Vmath::Vadd(nq, tmp1, 1, tmp2, 1, flux[j], 1);
    }
    }
    else
    {
        // F_{adv,u'_i,j} = (p'/ rho + \bar{u}_k u'_k) \delta_{ij}
        for (int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Zero(nq, flux[j], 1);

            if (i-1 == j)
            {
                // contruct p'/ rho term
                tmp0 = 1 / m_Rho0;
                Vmath::Smul(nq, tmp0, physfield[0], 1, flux[j], 1);

                // construct \bar{u}_k u'_k term
                Vmath::Zero(nq, tmp1, 1);
                for (int k = 0; k < flux.num_elements(); ++k)
                {
                    Vmath::Vmul(nq, physfield[k+1], 1, basefield[k], 1, tmp2, 1);
                    Vmath::Vadd(nq, tmp1, 1, tmp2, 1, tmp1, 1);
                }

                // add terms
                Vmath::Vadd(nq, flux[j], 1, tmp1, 1, flux[j], 1);
            }
        }
    }

}


void APE::v_PrimitiveToConservative()
{

}

void APE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                         Array<OneD,       Array<OneD, NekDouble> >&outarray,
                   const NekDouble time)
{
    int i;
    int ndim       = m_spacedim;
    int nvariables = inarray.num_elements();
    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();
    
    switch(m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            //-------------------------------------------------------
            //inarray in physical space

            Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
            Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);

            for (i = 0; i < nvariables; ++i)
            {
                physarray[i] = Array<OneD, NekDouble>(nq);
                modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
            }

            ConservativeToPrimitive(inarray,physarray);

            // get the advection part
            // input: physical space
            // output: modal space

            // straighforward DG
            WeakDGAdvection(physarray, modarray, false, true);

            // negate the outarray since moving terms to the rhs
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,modarray[i],1);
            }

            // Add "source term"
            // input: physical space
            // output: modal space (JOSEF)
            AddSource(physarray,modarray);

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->MultiplyByElmtInvMass(modarray[i],modarray[i]);
                m_fields[i]->BwdTrans(modarray[i],outarray[i]);
            }
            break;
        }

        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            Array<OneD, Array<OneD, NekDouble> > physarray(nvariables);
            Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);

            for (i = 0; i < nvariables; ++i)
            {
                physarray[i] = Array<OneD, NekDouble>(nq);
                modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
            }

            ConservativeToPrimitive(inarray,physarray);

            Array<OneD, Array<OneD, NekDouble> > fluxvector(ndim);
            for(i = 0; i < ndim; ++i)
            {
                fluxvector[i]    = Array<OneD, NekDouble>(nq);
            }

            Array<OneD,NekDouble> tmp(nq);
            Array<OneD, NekDouble>tmp1(nq);

            for(i = 0; i < nvariables; ++i)
            {
                // Get the ith component of the  flux vector in (physical space)
                APE::GetFluxVector(i, physarray, fluxvector);

            Vmath::Zero(nq, outarray[i], 1);
            for (int j = 0; j < ndim; ++j)
            {
                // Get the ith component of the  flux vector in (physical space)
                m_fields[0]->PhysDeriv(j,fluxvector[j],tmp1);
                Vmath::Vadd(nq, outarray[i], 1, tmp1, 1, outarray[i], 1);
            }
                Vmath::Neg(nq,outarray[i],1);
            }

            // Add "source term"
            // input: physical space
            // output: modal space
            AddSource(physarray,outarray);
            break;
        }

        default:
            ASSERTL0(false,"Unknown projection scheme for the APE");
            break;
    }
}


void APE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                Array<OneD,       Array<OneD, NekDouble> >&outarray,
                          const NekDouble time)
{
    int i;
    int nvariables = inarray.num_elements();

    switch(m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            ConservativeToPrimitive(inarray,outarray);
            SetBoundaryConditions(outarray,time);
            PrimitiveToConservative(outarray,outarray);
            break;
        }

        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            ConservativeToPrimitive(inarray,outarray);
            UnsteadySystem::SetBoundaryConditions(time);
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(outarray[i],coeffs);
                m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
            }
            PrimitiveToConservative(outarray,outarray);
            break;
        }

        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
    }
}

void APE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray,
                                NekDouble time)
{
    std::string varName;
    int nvariables = m_fields.num_elements();
    int cnt = 0;

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    {
        // Wall Boundary Condition
        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWall)
        {
            if (m_expdim == 1)
            {
                WallBoundary1D(n,inarray);
            }
            else if (m_expdim == 2)
            {
                WallBoundary2D(n,cnt,inarray);
            }
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }
        cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

void APE::WallBoundary2D(int bcRegion, int cnt,
                         Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = physarray.num_elements();
    
    // get physical values of the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
    }
    
    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;
    
    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
    {
        npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
        id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));

        switch(m_expdim)
        {
            case 1:
            {
                // negate the forward flux
                Vmath::Neg(npts,&Fwd[1][id2],1);
                break;
            }

            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
                Vmath::Vvtvp(npts,&Fwd[2][id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);

                Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
                Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);

                // negate the normal flux
                Vmath::Neg(npts,tmp_n,1);

                // rotate back to Cartesian
                Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
                Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);

                Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
                Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
                break;
            }

            case 3:
                ASSERTL0(false,"3D not implemented for Acoustic perturbation equations");
                break;

            default:
                ASSERTL0(false,"Illegal expansion dimension");
                break;
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npts,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
        }
    }
}

void APE::WallBoundary1D(int bcRegion, Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    ASSERTL0(false,"1D not yet working for APE");
}



void APE::NumericalFlux1D(Array<OneD, Array<OneD, NekDouble> > &physfield,
                          Array<OneD, Array<OneD, NekDouble> > &numfluxX)
{
    ASSERTL0(false,"1D DG not yet working for APE");
}

//Evaluation of the upwinded DG fluxes
void APE::NumericalFlux2D(Array<OneD, Array<OneD, NekDouble> > &physfield,
                          Array<OneD, Array<OneD, NekDouble> > &numfluxX,
                          Array<OneD, Array<OneD, NekDouble> > &numfluxY)
{
    int i;
    //Number the points of the "shared" edges of the elements
    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables      = 3; //p', u', v'

    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > Bwd(nvariables);
    Array<OneD, Array<OneD, NekDouble> > rotbasefield(nvariables);
    Array<OneD, Array<OneD, NekDouble> > rotbasefieldBwd(nvariables);

    int nq = m_fields[0]->GetNpoints();
    Array<OneD,NekDouble> x0(nq);
    Array<OneD,NekDouble> x1(nq);
    Array<OneD,NekDouble> x2(nq);

    m_fields[0]->GetCoords(x0,x1,x2);

    for (i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        rotbasefield[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        rotbasefieldBwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }
    
    // get the physical values at the trace
    for (i = 0; i < nvariables; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
        m_fields[i]->GetFwdBwdTracePhys(basefield[i],rotbasefield[i],rotbasefieldBwd[i]);
    }

    // rotate the values to the normal direction
    NekDouble tmpX, tmpY;
    for (i = 0; i < nTraceNumPoints; ++i)
    {
        tmpX =  Fwd[1][i]*m_traceNormals[0][i]+Fwd[2][i]*m_traceNormals[1][i];
        tmpY = -Fwd[1][i]*m_traceNormals[1][i]+Fwd[2][i]*m_traceNormals[0][i];
        Fwd[1][i] = tmpX;
        Fwd[2][i] = tmpY;

        tmpX =  Bwd[1][i]*m_traceNormals[0][i]+Bwd[2][i]*m_traceNormals[1][i];
        tmpY = -Bwd[1][i]*m_traceNormals[1][i]+Bwd[2][i]*m_traceNormals[0][i];
        Bwd[1][i] = tmpX;
        Bwd[2][i] = tmpY;

        //rotates the baseflow
        tmpX =  rotbasefield[0][i]*m_traceNormals[0][i]+rotbasefield[1][i]*m_traceNormals[1][i];
        tmpY = -rotbasefield[0][i]*m_traceNormals[1][i]+rotbasefield[1][i]*m_traceNormals[0][i];
        rotbasefield[0][i] = rotbasefield[2][i];
        rotbasefield[1][i] = tmpX;
        rotbasefield[2][i] = tmpY;
    }

    // Solve the Riemann problem
    NekDouble pflux, uflux, vflux;

    for (i = 0; i < nTraceNumPoints; ++i)
    {
        //cout << "Tracepoint: "<<i<<" of "<<nTraceNumPoints<<endl;

        switch(m_upwindType)
        {
            case eUpwind:
            {
                RiemannSolverUpwind(Fwd[0][i],Fwd[1][i],Fwd[2][i],
                        Bwd[0][i],Bwd[1][i],Bwd[2][i],
                        rotbasefield[0][i],rotbasefield[1][i],rotbasefield[2][i],
                        pflux, uflux, vflux );
                break;
            }

            default:
                ASSERTL0(false,"populate switch statement for upwind flux");
                break;
        }

        // rotate back to Cartesian
        numfluxX[0][i]  = pflux*m_traceNormals[0][i];
        numfluxY[0][i]  = pflux*m_traceNormals[1][i];
        numfluxX[1][i] = (uflux*m_traceNormals[0][i] - vflux*m_traceNormals[1][i]) * m_traceNormals[0][i];
        numfluxY[1][i] = (uflux*m_traceNormals[0][i] - vflux*m_traceNormals[1][i]) * m_traceNormals[1][i];
        numfluxX[2][i] = (uflux*m_traceNormals[1][i] + vflux*m_traceNormals[0][i]) * m_traceNormals[0][i];
        numfluxY[2][i] = (uflux*m_traceNormals[1][i] + vflux*m_traceNormals[0][i]) * m_traceNormals[1][i];
    }

}

void APE::RiemannSolverUpwind(NekDouble pL,     NekDouble uL,    NekDouble vL,
                              NekDouble pR,     NekDouble uR,    NekDouble vR,
                              NekDouble P0,     NekDouble U0 ,   NekDouble V0,
                              NekDouble &pflux, NekDouble &uflux, NekDouble &vflux )
{
    int nvariables      = 2;
    Array<OneD, NekDouble> characteristic(4);
    Array<OneD, NekDouble> W(2);
    Array<OneD, NekDouble> lambda(nvariables);
    Array<OneD, NekDouble> upphysfield(3);

    // compute the wave speeds
    lambda[0]=U0 + sqrt(P0*m_gamma*m_Rho0)/m_Rho0;
    lambda[1]=U0 - sqrt(P0*m_gamma*m_Rho0)/m_Rho0;

    // calculate the caracteristic variables
    //left characteristics
    characteristic[0] = pL/2 + uL*sqrt(P0*m_gamma*m_Rho0)/2;
    characteristic[1] = pL/2 - uL*sqrt(P0*m_gamma*m_Rho0)/2;
    //right characteristics
    characteristic[2] = pR/2 + uR*sqrt(P0*m_gamma*m_Rho0)/2;
    characteristic[3] = pR/2 - uR*sqrt(P0*m_gamma*m_Rho0)/2;

    //take left or right value of characteristic variable
    for (int j=0; j<nvariables; j++)
    {
        if (lambda[j]>=0)
        {
            W[j]=characteristic[j];
        }
        if(lambda[j]<0)
        {
            W[j]=characteristic[j+2];
        }
    }

    //calculate conservative variables from characteristics
    upphysfield[0]= W[0]+W[1];
    upphysfield[1]= (W[0]-W[1])/sqrt(P0*m_gamma*m_Rho0);
    upphysfield[2]= vL;

    // compute the fluxes
    pflux = U0*upphysfield[0] + m_gamma*P0*upphysfield[1];
    uflux = U0*upphysfield[1]+V0*upphysfield[2] + upphysfield[0]/m_Rho0;
    vflux = 0.0;
}

void APE::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
                                        Array<OneD,       Array<OneD, NekDouble> >&physout)
{
    int i;
    int nq = GetTotPoints();
    int nvariables = physin.num_elements();
    
    if(physin.get() == physout.get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble> >tmp(nvariables);
        for(i = 0; i < nvariables; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
        }
        for(i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(nq,tmp[i],1,physout[i],1);
        }
    }
    else
    {
        for(i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(nq,physin[i],1,physout[i],1);
        }
    }
}


void APE::v_ConservativeToPrimitive( )
{
    
}

void APE::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
                                        Array<OneD,       Array<OneD, NekDouble> >&physout)
{  
    int i;
    int nq = GetTotPoints();
    int nvariables = physin.num_elements();
    
    if(physin.get() == physout.get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble> >tmp(nvariables);
        for(i = 0; i < nvariables; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
        }
        for(i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(nq,tmp[i],1,physout[i],1);
        }
    }
    else
    {
        for(i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(nq,physin[i],1,physout[i],1);
        }
    }
}

// Initialise baseflow from the inputfile
void APE::InitialiseBaseFlowAnalytical(Array<OneD, Array<OneD, NekDouble> > &base,
                                       const NekDouble time)
{
    base = Array<OneD, Array<OneD, NekDouble> >(m_spacedim+1);
    int nq = m_fields[0]->GetNpoints();
    std::string velStr[3] = {"U0","V0","P0"};

    for(int i = 0; i <= m_spacedim; ++i)
    {
        base[i] = Array<OneD, NekDouble> (nq,0.0);
        EvaluateFunction(velStr[i], base[i], "Baseflow", time);
    }
}


// Get sourceterm for p' equation from the inputfile
void APE::GetSource(Array<OneD, NekDouble> &source, const NekDouble time)
{
    EvaluateFunction("S", source, "Source", time);
}

// Add sourceterm for p' equation obtained from GetSource
void APE::AddSource(const Array< OneD, Array< OneD, NekDouble > > &inarray,
                          Array< OneD, Array< OneD, NekDouble > > &outarray)
{
    int ncoeffs = outarray[0].num_elements();
    int nq      = inarray[0].num_elements();
    
    Array<OneD, NekDouble> source(nq);
    
    switch(m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            GetSource(source,m_time);

            //Source term solely for the p' equation (outarray[0])
            m_fields[0]->IProductWRTBase(source,source);
            Vmath::Vadd(ncoeffs,source,1,outarray[0],1,outarray[0],1);
            break;
        }

        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            GetSource(source,m_time);

            //Source term solely for the p' equation (outarray[0])
            Vmath::Vadd(ncoeffs,source,1,outarray[0],1,outarray[0],1);
            break;
        }

        default:
            ASSERTL0(false,"Unknown projection scheme for the APE");
            break;
    }
}

} //end of namespace

