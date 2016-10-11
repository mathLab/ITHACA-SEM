///////////////////////////////////////////////////////////////////////////////
//
// File APE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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
// Description: APE1/APE4 (Acoustic Perturbation Equations)
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <APESolver/EquationSystems/APE.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

using namespace std;

namespace Nektar
{
string APE::className = GetEquationSystemFactory().RegisterCreatorFunction(
            "APE", APE::create,
            "APE1/APE4 (Acoustic Perturbation Equations)");


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

    ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
             "Only Projection=DisContinuous supported by the APE class.");

    // Load isentropic coefficient, Ratio of specific heats
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    m_session->LoadParameter("IO_CFLSteps", m_cflsteps, 0);

    // Define Baseflow and source term fields
    switch (m_spacedim)
    {
        case 1:
        {
            for (int i = 0; i < m_spacedim + 2; ++i)
            {
                m_bfField = MemoryManager<MultiRegions::ContField1D>::
                    AllocateSharedPtr(m_session, m_graph);
            }
            break;
        }

        case 2:
        {
            for (int i = 0; i < m_spacedim + 2; ++i)
            {
                m_bfField = MemoryManager<MultiRegions::ContField2D>::
                    AllocateSharedPtr(m_session, m_graph);
            }
            break;
        }

        case 3:
        {
            for (int i = 0; i < m_spacedim + 2; ++i)
            {
                m_bfField = MemoryManager < MultiRegions::ContField3D >::
                    AllocateSharedPtr(m_session, m_graph);
            }
            break;
        }

        default:
        {

            ASSERTL0(false, "Expansion dimension not recognised");
            break;
        }
    }

    m_bfNames.push_back("p0");
    m_bfNames.push_back("rho0");
    m_bfNames.push_back("u0");
    m_bfNames.push_back("v0");
    m_bfNames.push_back("w0");

    // Resize the advection velocities vector to dimension of the problem
    m_bfNames.resize(m_spacedim + 2);

    // Initialize basefield
    m_bf = Array<OneD, Array<OneD, NekDouble> >(m_spacedim + 2);
    for (int i = 0; i < m_bf.num_elements(); ++i)
    {
        m_bf[i] = Array<OneD, NekDouble>(GetTotPoints());
    }
    EvaluateFunction(m_bfNames, m_bf, "Baseflow", m_time);

    m_forcing = SolverUtils::Forcing::Load(m_session, m_fields, m_spacedim + 1);

    // Do not forwards transform initial condition
    m_homoInitialFwd = false;

    // Define the normal velocity fields
    if (m_fields[0]->GetTrace())
    {
        m_traceBasefield = Array<OneD, Array<OneD, NekDouble> > (m_spacedim + 2);
        for (int i = 0; i < m_spacedim + 2; i++)
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
    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
                          riemName);
    m_riemannSolver->SetVector("N",         &APE::GetNormals,   this);
    m_riemannSolver->SetVector("basefield", &APE::GetBasefield, this);
    m_riemannSolver->SetAuxVec("vecLocs",   &APE::GetVecLocs,   this);
    m_riemannSolver->SetParam("Gamma",     &APE::GetGamma,     this);

    // Set up advection operator
    string advName;
    m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    m_advection = SolverUtils::GetAdvectionFactory()
                  .CreateInstance(advName, advName);
    m_advection->SetFluxVector(&APE::GetFluxVector, this);
    m_advection->SetRiemannSolver(m_riemannSolver);
    m_advection->InitObject(m_session, m_fields);

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&APE::DoOdeRhs,        this);
        m_ode.DefineProjection(&APE::DoOdeProjection, this);
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


NekDouble APE::GetCFLEstimate()
{
    int nElm = m_fields[0]->GetExpSize();
    const Array<OneD, int> expOrder = GetNumExpModesPerExp();

    Array<OneD, NekDouble> cfl(nElm, 0.0);
    Array<OneD, NekDouble> stdVelocity(nElm, 0.0);

    // Get standard velocity to compute the time-step limit
    GetStdVelocity(stdVelocity);

    // Factors to compute the time-step limit
    NekDouble alpha   = MaxTimeStepEstimator();
    NekDouble cLambda = 0.2; // Spencer book-317

    // Loop over elements to compute the time-step limit for each element
    for (int el = 0; el < nElm; ++el)
    {
        NekDouble lambdaMax = stdVelocity[el] * cLambda
            * (expOrder[el] - 1) * (expOrder[el] - 1);
        cfl[el] = m_timestep * lambdaMax / alpha;
    }

    // Get the minimum time-step limit and return the time-step
    NekDouble maxCFL = Vmath::Vmax(nElm, cfl, 1);
    m_comm->AllReduce(maxCFL, LibUtilities::ReduceMax);

    return maxCFL;
}


/**
 * @brief Return the flux vector for the APE equations.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux. flux[eq][dir][pt]
 */
void APE::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    int nq = physfield[0].num_elements();
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    ASSERTL1(flux[0].num_elements() == m_spacedim,
                 "Dimension of flux array and velocity array do not match");

    // F_{adv,p',j} = \gamma p_0 u'_j + p' \bar{u}_j
    for (int j = 0; j < m_spacedim; ++j)
    {
        Vmath::Zero(nq, flux[0][j], 1);

        // construct \gamma p_0 u'_j term
        Vmath::Smul(nq, m_gamma, m_bf[0], 1, tmp1, 1);
        Vmath::Vmul(nq, tmp1, 1, physfield[j+1], 1, tmp1, 1);

        // construct p' \bar{u}_j term
        Vmath::Vmul(nq, physfield[0], 1, m_bf[j+2], 1, tmp2, 1);

        // add both terms
        Vmath::Vadd(nq, tmp1, 1, tmp2, 1, flux[0][j], 1);
    }

    for (int i = 1; i < flux.num_elements(); ++i)
    {
        ASSERTL1(flux[i].num_elements() == m_spacedim,
                 "Dimension of flux array and velocity array do not match");

        // F_{adv,u'_i,j} = (p'/ \bar{rho} + \bar{u}_k u'_k) \delta_{ij}
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Zero(nq, flux[i][j], 1);

            if (i - 1 == j)
            {
                // contruct p'/ \bar{rho} term
                Vmath::Vdiv(nq, physfield[0], 1, m_bf[1], 1, flux[i][j], 1);

                // construct \bar{u}_k u'_k term
                Vmath::Zero(nq, tmp1, 1);
                for (int k = 0; k < m_spacedim; ++k)
                {
                    Vmath::Vvtvp(nq, physfield[k + 1], 1, m_bf[k + 2 ], 1, tmp1, 1, tmp1, 1);
                }

                // add terms
                Vmath::Vadd(nq, flux[i][j], 1, tmp1, 1, flux[i][j], 1);
            }
        }
    }
}


/**
 * @brief v_PostIntegrate
 */
bool APE::v_PreIntegrate(int step)
{
    EvaluateFunction(m_bfNames, m_bf, "Baseflow", m_time);

    Array<OneD, NekDouble> tmpC(GetNcoeffs());
    std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
    for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
    {
        for (int i = 0; i < (*x)->GetForces().num_elements(); ++i)
        {
            m_bfField->IProductWRTBase((*x)->GetForces()[i], tmpC);
            m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
            m_bfField->LocalToGlobal(tmpC, tmpC);
            m_bfField->GlobalToLocal(tmpC, tmpC);
            m_bfField->BwdTrans(tmpC, (*x)->UpdateForces()[i]);
        }
    }

    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        // ensure the field is C0-continuous
        m_bfField->IProductWRTBase(m_bf[i], tmpC);
        m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
        m_bfField->LocalToGlobal(tmpC, tmpC);
        m_bfField->GlobalToLocal(tmpC, tmpC);
        m_bfField->BwdTrans(tmpC, m_bf[i]);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}


/**
 * @brief v_PostIntegrate
 */
bool APE::v_PostIntegrate(int step)
{
    if (m_cflsteps && !((step + 1) % m_cflsteps))
    {
        NekDouble cfl = GetCFLEstimate();
        if (m_comm->GetRank() == 0)
        {
            cout << "CFL: " << cfl << endl;
        }
    }

    return UnsteadySystem::v_PostIntegrate(step);
}


/**
 * @brief Compute the right-hand side.
 */
void APE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                         Array<OneD,       Array<OneD, NekDouble> >&outarray,
                   const NekDouble time)
{
    int nVariables = inarray.num_elements();
    int nq = GetTotPoints();

    // WeakDG does not use advVel, so we only provide a dummy array
    Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);
    m_advection->Advect(nVariables, m_fields, advVel, inarray, outarray, time);

    // Negate the LHS terms
    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Neg(nq, outarray[i], 1);
    }

    std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
    for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
    {
        (*x)->Apply(m_fields, outarray, outarray, m_time);
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
    int nvariables = inarray.num_elements();
    int nq = m_fields[0]->GetNpoints();

    // deep copy
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
    }

    SetBoundaryConditions(outarray, time);
}


/**
 * @brief Apply the Boundary Conditions to the APE equations.
 */
void APE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray,
                                NekDouble time)
{
    std::string varName;
    int nvariables = m_fields.num_elements();
    int cnt        = 0;
    int nTracePts  = GetTraceTotPoints();

    // Extract trace for boundaries. Needs to be done on all processors to avoid
    // deadlock.
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // loop over Boundary Regions
    for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    {
        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"Wall"))
        {
            WallBC(n, cnt, Fwd, inarray);
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
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
                 Array<OneD, Array<OneD, NekDouble> > &Fwd,
                 Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int nVariables = physarray.num_elements();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

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
 * @brief Compute the advection velocity in the standard space
 * for each element of the expansion.
 *
 * @param stdV       Standard velocity field.
 */
void APE::GetStdVelocity(Array<OneD, NekDouble> &stdV)
{
    int nElm = m_fields[0]->GetExpSize();

    ASSERTL1(stdV.num_elements() ==  nElm,  "stdV malformed");

    Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
    Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim+1);
    LibUtilities::PointsKeyVector ptsKeys;

    // Zero output array
    Vmath::Zero(stdV.num_elements(), stdV, 1);

    int cnt = 0;

    for (int el = 0; el < nElm; ++el)
    {
        ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

        // Possible bug: not multiply by jacobian??
        const SpatialDomains::GeomFactorsSharedPtr metricInfo =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
        const Array<TwoD, const NekDouble> &gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
                ->GetDerivFactors(ptsKeys);

        int nq = m_fields[0]->GetExp(el)->GetTotPoints();

        for (int i = 0; i < m_spacedim; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(nq, 0.0);

            velocity[i] = Array<OneD, NekDouble>(nq, 0.0);
            for (int j = 0; j < nq; ++j)
            {
                // The total advection velocity is v+c, so we need to scale c by
                // adding it before we do the transformation.
                NekDouble c = sqrt(m_gamma * m_bf[0][cnt+j] / m_bf[1][cnt+j]);
                velocity[i][j] = m_bf[i+2][cnt+j] + c;
            }
        }

        // scale the velocity components
        if (metricInfo->GetGtype() == SpatialDomains::eDeformed)
        {
            // d xi/ dx = gmat = 1/J * d x/d xi
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, gmat[i], 1, velocity[0], 1, stdVelocity[i], 1);
                for (int j = 1; j < m_spacedim; ++j)
                {
                    Vmath::Vvtvp(nq, gmat[m_spacedim * j + i], 1, velocity[j],
                            1, stdVelocity[i], 1, stdVelocity[i], 1);
                }
            }
        }
        else
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Smul(nq, gmat[i][0], velocity[0], 1, stdVelocity[i], 1);
                for (int j = 1; j < m_spacedim; ++j)
                {
                    Vmath::Svtvp(nq, gmat[m_spacedim * j + i][0], velocity[j],
                            1, stdVelocity[i], 1, stdVelocity[i], 1);
                }
            }
        }

        // compute the max absolute velocity of the element
        for (int i = 0; i < nq; ++i)
        {
            NekDouble pntVelocity = 0.0;
            for (int j = 0; j < m_spacedim; ++j)
            {
                pntVelocity += stdVelocity[j][i] * stdVelocity[j][i];
            }
            pntVelocity = sqrt(pntVelocity);

            if (pntVelocity > stdV[el])
            {
                stdV[el] = pntVelocity;
            }
        }

        cnt += nq;
    }
}


void APE::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
    std::vector<std::string>             &variables)
{
    for (int i = 0; i < m_spacedim + 2; i++)
    {
        Array<OneD, NekDouble> tmpC(GetNcoeffs());

        // ensure the field is C0-continuous
        m_bfField->IProductWRTBase(m_bf[i], tmpC);
        m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
        m_bfField->LocalToGlobal(tmpC, tmpC);
        m_bfField->GlobalToLocal(tmpC, tmpC);

        variables.push_back(m_bfNames[i]);
        fieldcoeffs.push_back(tmpC);
    }

    int f = 0;
    std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
    for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
    {
        for (int i = 0; i < (*x)->GetForces().num_elements(); ++i)
        {
            Array<OneD, NekDouble> tmpC(GetNcoeffs());

            m_bfField->IProductWRTBase((*x)->GetForces()[i], tmpC);
            m_bfField->MultiplyByElmtInvMass(tmpC, tmpC);
            m_bfField->LocalToGlobal(tmpC, tmpC);
            m_bfField->GlobalToLocal(tmpC, tmpC);

            variables.push_back("F_" + boost::lexical_cast<string>(f) +
                                "_" + m_session->GetVariable(i));
            fieldcoeffs.push_back(tmpC);
        }
        f++;
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
    for (int i = 0; i < m_spacedim + 2; i++)
    {
        m_fields[0]->ExtractTracePhys(m_bf[i], m_traceBasefield[i]);
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

} //end of namespace

