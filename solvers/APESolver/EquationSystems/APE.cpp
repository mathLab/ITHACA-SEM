///////////////////////////////////////////////////////////////////////////////
//
// File APE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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


/**
 * @brief Initialization object for the APE class.
 */
void APE::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    // TODO: We have a bug somewhere in the 1D boundary conditions. Therefore 1D
    // problems are currently disabled. This should get fixed in the future.
    ASSERTL0(m_spacedim > 1, "1D problems currently not supported by the APE class.");

    // Load constant incompressible density
    m_session->LoadParameter("Rho0", m_Rho0, 1.204);

    // Load isentropic coefficient, Ratio of specific heats
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    // Define Baseflow fields
    m_basefield = Array<OneD, Array<OneD, NekDouble> >(m_spacedim + 1);
    m_basefield_names.push_back("P0");
    m_basefield_names.push_back("U0");
    m_basefield_names.push_back("V0");
    m_basefield_names.push_back("W0");

    // Resize the advection velocities vector to dimension of the problem
    m_basefield_names.resize(m_spacedim + 1);

    // if discontinuous  determine numerical flux to use
    if (m_projectionType == MultiRegions::eDiscontinuous)
    {
        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Define the normal velocity fields
        if (m_fields[0]->GetTrace())
        {
            m_traceBasefield = Array<OneD, Array<OneD, NekDouble> > (m_spacedim+1);
            for (int i = 0; i < m_spacedim + 1; i++)
            {
                m_traceBasefield[i] = Array<OneD, NekDouble>(GetTraceNpoints());
            }
        }

        // Set up locations of velocity and base velocity vectors.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            // u', v', w'
            m_vecLocs[0][i] = 1 + i;
        }

        string riemName;
        m_session->LoadSolverInfo("UpwindType", riemName, "APEUpwind");
        m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(riemName);
        m_riemannSolver->SetVector("N",         &APE::GetNormals,   this);
        m_riemannSolver->SetVector("basefield", &APE::GetBasefield, this);
        m_riemannSolver->SetAuxVec("vecLocs",   &APE::GetVecLocs,   this);
        m_riemannSolver->SetParam ("Gamma",     &APE::GetGamma,     this);
        m_riemannSolver->SetParam ("Rho",       &APE::GetRho,       this);
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


/**
 * @brief Destructor for APE class.
 */
APE::~APE()
{
    
}


/**
 *
 */
void APE::v_DoInitialise()
{
    SetInitialConditions();
}


/**
 * @brief Compute the numerical flux through the element boundaries.
 *
 */
void APE::v_NumericalFlux(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &numflux)
{
    //Number the points of the shared edges of the elements
    int ntp  = GetTraceTotPoints();
    int nvar = physfield.num_elements();

    // temporary arrays
    Array<OneD, Array<OneD, NekDouble> >  Fwd(nvar);
    Array<OneD, Array<OneD, NekDouble> >  Bwd(nvar);

    for (int i = 0; i < nvar; ++i)
    {
        Fwd[i]  = Array<OneD, NekDouble>(ntp);
        Bwd[i]  = Array<OneD, NekDouble>(ntp);
    }

    // get the physical values at the trace
    for (int i = 0; i < nvar; ++i)
    {
        m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd[i],Bwd[i]);
    }

    // Solve the Riemann problem
    m_riemannSolver->Solve(m_spacedim, Fwd, Bwd, numflux);

}


/**
 *
 */
void APE::v_NumericalFlux(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &numfluxX,
        Array<OneD, Array<OneD, NekDouble> > &numfluxY )
{
    ASSERTL0(false, "This function is not implemented for this equation.");
}


/**
 * @brief Compute the flux vectors.
 */
void APE::v_GetFluxVector(const int i, const int j,
                           Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
{
    v_GetFluxVector(i, physfield, flux);
}


/**
 * @brief Compute the flux vectors.
 */
void APE::v_GetFluxVector(const int i,
                             Array<OneD, Array<OneD, NekDouble> > &physfield,
                             Array<OneD, Array<OneD, NekDouble> > &flux)
{
    UpdateBasefield();

    ASSERTL1(flux.num_elements() == m_basefield.num_elements() - 1,
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
            Vmath::Smul(nq, m_gamma, m_basefield[0], 1, tmp1, 1);
            Vmath::Vmul(nq, tmp1, 1, physfield[j+1], 1, tmp1, 1);

            // construct p' \bar{u}_j term
            Vmath::Vmul(nq, physfield[0], 1, m_basefield[j+1], 1, tmp2, 1);

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
                    Vmath::Vmul(nq, physfield[k+1], 1, m_basefield[k+1], 1, tmp2, 1);
                    Vmath::Vadd(nq, tmp1, 1, tmp2, 1, tmp1, 1);
                }

                // add terms
                Vmath::Vadd(nq, flux[j], 1, tmp1, 1, flux[j], 1);
            }
        }
    }

}


/**
 * @brief Compute the right-hand side.
 */
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

            Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);

            for (i = 0; i < nvariables; ++i)
            {
                modarray[i]  = Array<OneD, NekDouble>(ncoeffs);
            }

            // get the advection part
            // input: physical space
            // output: modal space

            // straighforward DG
            WeakDGAdvection(inarray, modarray, true, true);

            // negate the outarray since moving terms to the rhs
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,modarray[i],1);
            }

            // Add "source term"
            // input: physical space
            // output: modal space (JOSEF)
            AddSource(inarray, modarray);

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

            // deep copy
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(nq,inarray[i],1,physarray[i],1);
            }

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


/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void APE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                Array<OneD,       Array<OneD, NekDouble> >&outarray,
                          const NekDouble time)
{
    int i;
    int nvariables = inarray.num_elements();
    int nq = m_fields[0]->GetNpoints();

    // deep copy
    for(int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq,inarray[i],1,outarray[i],1);
    }

    switch(m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            SetBoundaryConditions(outarray,time);
            break;
        }

        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            UnsteadySystem::SetBoundaryConditions(time);
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for(i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(outarray[i],coeffs);
                m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
            }
            break;
        }

        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
    }
}


/**
 * @brief Apply the Boundary Conditions to the APE equations.
 */
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
            WallBC(n, cnt, inarray);
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


/**
 * @brief Wall boundary conditions for the APE equations.
 */
void APE::WallBC(int bcRegion, int cnt,
                 Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int nTracePts = GetTraceTotPoints();
    int nVariables = physarray.num_elements();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    // Get physical values of the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nVariables);
    for (int i = 0; i < nVariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
    }

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, nBCEdgePts;
    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    for (int e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

        // For 2D/3D, define: v* = v - 2(v.n)n
        Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);

        // Calculate (v.n)
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts,
                         &Fwd[1+i][id2], 1,
                         &m_traceNormals[i][id2], 1,
                         &tmp[0], 1,
                         &tmp[0], 1);
        }

        // Calculate 2.0(v.n)
        Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);

        // Calculate v* = v - 2.0(v.n)n
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts,
                         &tmp[0], 1,
                         &m_traceNormals[i][id2], 1,
                         &Fwd[1+i][id2], 1,
                         &Fwd[1+i][id2], 1);
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts,
                         &Fwd[i][id2], 1,
                         &(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1], 1);
        }
    }
}


/**
 * @brief sourceterm for p' equation obtained from GetSource
 */
void APE::AddSource(const Array< OneD, Array< OneD, NekDouble > > &inarray,
                    Array< OneD, Array< OneD, NekDouble > > &outarray)
{
    int ncoeffs = outarray[0].num_elements();
    int nq      = inarray[0].num_elements();
    Array<OneD, NekDouble> source(nq);

    EvaluateFunction("S", source, "Source", m_time);
    if ( m_projectionType == MultiRegions::eDiscontinuous )
    {
        m_fields[0]->IProductWRTBase(source,source);
    }
    Vmath::Vadd(ncoeffs,source,1,outarray[0],1,outarray[0],1);

}


void APE::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    std::vector<std::string>             &variables)
{
    UpdateBasefield();

    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    for (int i = 0; i < m_spacedim + 1; i++)
    {
        variables.push_back(m_basefield_names[i]);

        Array<OneD, NekDouble> tmpFwd(nCoeffs);
        m_fields[0]->FwdTrans(m_basefield[i], tmpFwd);
        fieldcoeffs.push_back(tmpFwd);
    }
}


/**
 * @brief Get the normal vectors.
 */
const Array<OneD, const Array<OneD, NekDouble> > &APE::GetNormals()
{
    return m_traceNormals;
}


/**
 * @brief Get the locations of the components of the directed fields within the fields array.
 */
const Array<OneD, const Array<OneD, NekDouble> > &APE::GetVecLocs()
{
    return m_vecLocs;
}


/**
 * @brief Get the baseflow field.
 */
const Array<OneD, const Array<OneD, NekDouble> > &APE::GetBasefield()
{
    for (int i = 0; i < m_spacedim +1; i++)
    {
        m_fields[0]->ExtractTracePhys(m_basefield[i], m_traceBasefield[i]);
    }
    return m_traceBasefield;
}


/**
 * @brief Get the heat capacity ratio.
 */
NekDouble APE::GetGamma()
{
    return m_gamma;
}


/**
 * @brief Get the density.
 */
NekDouble APE::GetRho()
{
    return m_Rho0;
}

void APE::UpdateBasefield()
{
    static NekDouble last_update = -1.0;

    if (m_time > last_update)
    {
        EvaluateFunction(m_basefield_names, m_basefield, "Baseflow");
        last_update = m_time;
    }
}


} //end of namespace

