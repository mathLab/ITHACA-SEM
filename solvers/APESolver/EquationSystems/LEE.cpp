///////////////////////////////////////////////////////////////////////////////
//
// File LEE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: Linearized Euler Equations
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <APESolver/EquationSystems/LEE.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/HexExp.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

using namespace std;

namespace Nektar
{
string LEE::className = GetEquationSystemFactory().RegisterCreatorFunction(
            "LEE", LEE::create,
            "Linearized Euler Equations");


LEE::LEE(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
{
}


/**
 * @brief Initialization object for the LEE class.
 */
void LEE::v_InitObject()
{
    UnsteadySystem::v_InitObject();

    // TODO: We have a bug somewhere in the 1D boundary conditions. Therefore 1D
    // problems are currently disabled. This should get fixed in the future.
//     ASSERTL0(m_spacedim > 1, "1D problems currently not supported by the LEE class.");

    ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
             "Only Projection=DisContinuous supported by the LEE class.");

    // Load isentropic coefficient, Ratio of specific heats
    m_session->LoadParameter("Gamma", m_gamma, 1.4);

    m_bfNames.push_back("p0");
    m_bfNames.push_back("rho0");
    m_bfNames.push_back("u0");
    m_bfNames.push_back("v0");
    m_bfNames.push_back("w0");

    // Resize the advection velocities vector to dimension of the problem
    m_bfNames.resize(m_spacedim + 2);

    // Initialize basefield
    m_bf = Array<OneD, Array<OneD, NekDouble> >(m_bfNames.size());
    for (int i = 0; i < m_bf.num_elements(); ++i)
    {
        m_bf[i] = Array<OneD, NekDouble>(GetTotPoints());
    }
    GetFunction("Baseflow", m_fields[0], true)->Evaluate(m_bfNames, m_bf, m_time);

    m_forcing = SolverUtils::Forcing::Load(m_session, m_fields, m_spacedim + 2);

    // Do not forwards transform initial condition
    m_homoInitialFwd = false;

    // Define the normal velocity fields
    m_bfFwdBwd = Array<OneD, Array<OneD, NekDouble> > (2*m_bfNames.size());
    for (int i = 0; i < m_bfFwdBwd.num_elements(); i++)
    {
        m_bfFwdBwd[i] = Array<OneD, NekDouble>(GetTraceNpoints(), 0.0);
    }

    // Set up locations of velocity and base velocity vectors.
    m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
    m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        // u', v', w'
        m_vecLocs[0][i] = 2 + i;
    }

    string riemName;
    m_session->LoadSolverInfo("UpwindType", riemName, "LEEUpwind");
    m_riemannSolver = SolverUtils::GetRiemannSolverFactory().CreateInstance(
                          riemName);
    m_riemannSolver->SetVector("N",         &LEE::GetNormals,   this);
    m_riemannSolver->SetVector("basefieldFwdBwd", &LEE::GetBasefieldFwdBwd, this);
    m_riemannSolver->SetAuxVec("vecLocs",   &LEE::GetVecLocs,   this);
    m_riemannSolver->SetParam("Gamma",      &LEE::GetGamma,     this);

    // Set up advection operator
    string advName;
    m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
    m_advection = SolverUtils::GetAdvectionFactory()
                  .CreateInstance(advName, advName);
    m_advection->SetFluxVector(&LEE::GetFluxVector, this);
    m_advection->SetRiemannSolver(m_riemannSolver);
    m_advection->InitObject(m_session, m_fields);

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&LEE::DoOdeRhs,        this);
        m_ode.DefineProjection(&LEE::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit LEE not set up.");
    }
}


/**
 * @brief Destructor for LEE class.
 */
LEE::~LEE()
{
    
}


/**
 * @brief Return the flux vector for the LEE equations.
 *
 * @param physfield   Fields.
 * @param flux        Resulting flux. flux[eq][dir][pt]
 */
void LEE::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    int nq = physfield[0].num_elements();
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    ASSERTL1(flux[0].num_elements() == m_spacedim,
                 "Dimension of flux array and velocity array do not match");

    Array<OneD, NekDouble> p0 = m_bf[0];
    Array<OneD, NekDouble> rho0 = m_bf[1];
    Array<OneD, Array<OneD, NekDouble> > u0(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        u0[i] = m_bf[2+i];
    }

    Array<OneD, NekDouble> p = physfield[0];
    Array<OneD, NekDouble> rho = physfield[1];
    Array<OneD, Array<OneD, NekDouble> > ru(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        ru[i] = physfield[2+i];
    }

    // F_{adv,p',j} = c0^2 * ru_j + u0_j * p
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 0;
        Vmath::Zero(nq, flux[i][j], 1);

        // c0^2 * ru_j
        Vmath::Smul(nq, m_gamma, p0, 1, tmp1, 1);
        Vmath::Vdiv(nq, tmp1, 1, rho0, 1, tmp1, 1);
        Vmath::Vmul(nq, tmp1, 1, ru[j], 1, tmp1, 1);

        // u0_j * p
        Vmath::Vmul(nq, u0[j], 1, p, 1, tmp2, 1);

        Vmath::Vadd(nq, tmp1, 1, tmp2, 1, flux[i][j], 1);
    }

    // F_{adv,rho',j} = u0_j * rho' + ru_j
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 1;
        Vmath::Zero(nq, flux[i][j], 1);

        // u0_j * rho'
        Vmath::Vmul(nq, u0[j], 1, rho, 1, tmp1, 1);

        Vmath::Vadd(nq, tmp1, 1, ru[j], 1, flux[i][j], 1);
    }

    for (int i = 2; i < flux.num_elements(); ++i)
    {
        ASSERTL1(flux[i].num_elements() == m_spacedim,
                 "Dimension of flux array and velocity array do not match");

        // F_{adv,u'_i,j} = ru_i * u0_j + delta_ij * p
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Zero(nq, flux[i][j], 1);

            // ru_i * u0_j
            Vmath::Vmul(nq, ru[i-2], 1, u0[j], 1, flux[i][j], 1);

            // kronecker delta
            if (i - 2 == j)
            {
                // delta_ij * p
                Vmath::Vadd(nq, p, 1, flux[i][j], 1, flux[i][j], 1);
            }
        }
    }
}



void LEE::GetLinTerm(const Array< OneD, const Array< OneD, NekDouble > > &inarray,
                     Array<OneD, Array<OneD, NekDouble> > &linTerm)
{
    int nq = GetTotPoints();

    Array<OneD, NekDouble> p0 = m_bf[0];
    Array<OneD, NekDouble> rho0 = m_bf[1];
    Array<OneD, Array<OneD, NekDouble> > u0(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        u0[i] = m_bf[2+i];
    }

    Array<OneD, NekDouble> p = inarray[0];
    Array<OneD, NekDouble> rho = inarray[1];
    Array<OneD, Array<OneD, NekDouble> > ru(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        ru[i] = inarray[2+i];
    }

    Array<OneD, NekDouble> grad(nq);
    Array<OneD, NekDouble> tmp1(nq);
    Array<OneD, NekDouble> tmp2(nq);

    // p
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 0;
        Vmath::Zero(nq, linTerm[i], 1);

        // (1-gamma) ( 1/rho0 * dp0/dx_j * ru_j - p * du0_j/dx_j )
        for (int j = 0; j < m_spacedim; ++j)
        {
            // 1/rho0 * dp0/dx_j * ru_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], p0, grad);
            Vmath::Vmul(nq, grad, 1, ru[j], 1, tmp1, 1);
            Vmath::Vdiv(nq, tmp1, 1, rho0, 1, tmp1, 1);

            // p * du0_j/dx_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], u0[j], grad);
            Vmath::Vmul(nq, grad, 1, p, 1, tmp2, 1);

            // (1-gamma) ( 1/rho0 * dp0/dx_j * ru_j - p * du0_j/dx_j )
            Vmath::Vsub(nq, tmp1, 1, tmp2, 1, tmp1, 1);
            Vmath::Smul(nq, (1-m_gamma), tmp1, 1, tmp1, 1);

            Vmath::Vadd(nq, tmp1, 1, linTerm[i], 1, linTerm[i], 1);
        }
    }

    // rho
    for (int j = 0; j < m_spacedim; ++j)
    {
        int i = 1;
        Vmath::Zero(nq, linTerm[i], 1);
    }


    // ru_i
    for (int i = 2; i < m_spacedim + 2; ++i)
    {
        Vmath::Zero(nq, linTerm[i], 1);

        // du0_i/dx_j * (u0_j * rho + ru_j)
        for (int j = 0; j < m_spacedim; ++j)
        {
            // d u0_i / d x_j
            m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[j], u0[i-2], grad);

            // u0_j * rho + ru_j
            Vmath::Vmul(nq, u0[j], 1, rho, 1, tmp1, 1);
            Vmath::Vadd(nq, ru[j], 1, tmp1, 1, tmp1, 1);

            Vmath::Vmul(nq, grad, 1, tmp1, 1, tmp1, 1);

            Vmath::Vadd(nq, tmp1, 1, linTerm[i], 1, linTerm[i], 1);
        }
    }

    for (int i = 0; i < m_spacedim + 2; ++i)
    {
        Array<OneD, NekDouble> tmpC(GetNcoeffs());

        m_fields[0]->FwdTrans(linTerm[i], tmpC);
        m_fields[0]->BwdTrans(tmpC, linTerm[i]);
    }
}


/**
 * @brief v_PostIntegrate
 */
bool LEE::v_PreIntegrate(int step)
{
    GetFunction("Baseflow", m_fields[0], true)->Evaluate(m_bfNames, m_bf, m_time);

    return UnsteadySystem::v_PreIntegrate(step);
}


/**
 * @brief v_PostIntegrate
 */
bool LEE::v_PostIntegrate(int step)
{
    return UnsteadySystem::v_PostIntegrate(step);
}


/**
 * @brief Compute the right-hand side.
 */
void LEE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                         Array<OneD,       Array<OneD, NekDouble> >&outarray,
                   const NekDouble time)
{
    int nVariables = inarray.num_elements();
    int nq = GetTotPoints();

    // WeakDG does not use advVel, so we only provide a dummy array
    Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);
    m_advection->Advect(nVariables, m_fields, advVel, inarray, outarray, time);

    Array<OneD, Array<OneD, NekDouble> > linT(nVariables);
    for (int i = 0; i < nVariables; ++i)
    {
        linT[i] = Array<OneD, NekDouble> (nq);
    }
    GetLinTerm(inarray, linT);
    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Vadd(nq, outarray[i], 1, linT[i], 1, outarray[i], 1);
    }

    // Negate the LHS terms
    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Neg(nq, outarray[i], 1);
    }

    std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
    for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
    {
        (*x)->Apply(m_fields, inarray, outarray, m_time);
    }
}


/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void LEE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
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

    UpdateBasefieldFwdBwd();

    SetBoundaryConditions(outarray, time);
}


/**
 * @brief Apply the Boundary Conditions to the LEE equations.
 */
void LEE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray,
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
    Array<OneD, Array<OneD, NekDouble> > bfFwd = GetBasefieldFwdBwd();

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
    {
        std::string userDefStr =
            m_fields[0]->GetBndConditions()[n]->GetUserDefined();

        if (!userDefStr.empty())
        {
            // Wall Boundary Condition
            if (boost::iequals(userDefStr, "Wall"))
            {
                WallBC(n, cnt, Fwd, inarray);
            }
            else if (boost::iequals(userDefStr, "WhiteNoise"))
            {
                WhiteNoiseBC(n, cnt, Fwd, bfFwd, inarray);
            }
            else if (boost::iequals(userDefStr, "RiemannInvariantBC"))
            {
                RiemannInvariantBC(n, cnt, Fwd, bfFwd, inarray);
            }
            else if (boost::iequals(userDefStr, "TimeDependent"))
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    varName = m_session->GetVariable(i);
                    m_fields[i]->EvaluateBoundaryConditions(time, varName);
                }
            }
            else
            {
                string errmsg = "Unrecognised boundary condition: ";
                errmsg += userDefStr;
                ASSERTL0(false, errmsg.c_str());
            }
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                varName = m_session->GetVariable(i);
                m_fields[i]->EvaluateBoundaryConditions(time, varName);
            }
        }

        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}


/**
 * @brief Wall boundary conditions for the LEE equations.
 */
void LEE::WallBC(int bcRegion, int cnt,
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
                         &Fwd[2+i][id2], 1,
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
                         &Fwd[2+i][id2], 1,
                         &Fwd[2+i][id2], 1);
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
 * @brief Outflow characteristic boundary conditions for compressible
 * flow problems.
 */
void LEE::RiemannInvariantBC(int bcRegion,
                             int cnt,
                             Array<OneD, Array<OneD, NekDouble> > &Fwd,
                             Array<OneD, Array<OneD, NekDouble> > &BfFwd,
                             Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int id1, id2, nBCEdgePts;
    int nVariables = physarray.num_elements();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    for (int e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]
                         ->GetBndCondExpansions()[bcRegion]
                         ->GetExp(e)
                         ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

        // Calculate (v.n)
        Array<OneD, NekDouble> Vn(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts,
                         &Fwd[1 + i][id2], 1,
                         &m_traceNormals[i][id2], 1,
                         &Vn[0], 1,
                         &Vn[0], 1);
        }

        // Calculate (v0.n)
        Array<OneD, NekDouble> Vn0(nBCEdgePts, 0.0);
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts,
                         &BfFwd[2 + i][id2], 1,
                         &m_traceNormals[i][id2], 1,
                         &Vn0[0], 1,
                         &Vn0[0], 1);
        }

        for (int i = 0; i < nBCEdgePts; ++i)
        {
            NekDouble c = sqrt(m_gamma * BfFwd[0][id2 + i] / BfFwd[1][id2 + i]);

            NekDouble l0 = Vn0[i] + c;
            NekDouble l1 = Vn0[i] - c;

            NekDouble h0, h1;

            // outgoing
            if (l0 > 0)
            {
                // p/2 + u*c*rho0/2
                h0 = Fwd[0][id2 + i] / 2 + Vn[i] * c * BfFwd[1][id2 + i] / 2;
            }
            // incoming
            else
            {
                h0 = 0.0;
            }

            // outgoing
            if (l1 > 0)
            {
                // p/2 - u*c*rho0/2
                h1 = Fwd[0][id2 + i] / 2 - Vn[i] * c * BfFwd[1][id2 + i] / 2;
            }
            // incoming
            else
            {
                h1 = 0.0;
            }

            // compute primitive variables
            // p = h0 + h1
            // u = ( h0 - h1) / (c*rho0)
            Fwd[0][id2 + i] = h0 + h1;
            NekDouble VnNew = (h0 - h1) / (c * BfFwd[1][id2 + i]);

            // adjust velocity pert. according to new value
            for (int j = 0; j < m_spacedim; ++j)
            {
                Fwd[1 + j][id2 + i] =
                    Fwd[1 + j][id2 + i] +
                    (VnNew - Vn[i]) * m_traceNormals[j][id2 + i];
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts,
                         &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1], 1);
        }
    }
}


/**
 * @brief Wall boundary conditions for the LEE equations.
 */
void LEE::WhiteNoiseBC(int bcRegion,
                       int cnt,
                       Array<OneD, Array<OneD, NekDouble> > &Fwd,
                       Array<OneD, Array<OneD, NekDouble> > &BfFwd,
                       Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int id1, id2, nBCEdgePts;
    int nVariables = physarray.num_elements();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    if (m_rng.count(bcRegion) == 0)
    {
        m_rng[bcRegion] = boost::mt19937(bcRegion);
    }

    ASSERTL0(
        m_fields[0]->GetBndConditions()[bcRegion]->GetBoundaryConditionType() ==
            SpatialDomains::eDirichlet,
        "WhiteNoise BCs must be Dirichlet type BCs");

    LibUtilities::Equation cond =
        boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(
            m_fields[0]->GetBndConditions()[bcRegion])
            ->m_dirichletCondition;
    NekDouble sigma = cond.Evaluate();

    ASSERTL0(sigma > NekConstants::kNekZeroTol,
             "sigma must be greater than zero");

    // random velocity perturbation
    if (m_whiteNoiseBC_lastUpdate < m_time)
    {
        m_whiteNoiseBC_lastUpdate = m_time;

        boost::normal_distribution<> dist(0, sigma);
        m_whiteNoiseBC_p = dist(m_rng[bcRegion]);
    }

    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    for (int e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]
                         ->GetBndCondExpansions()[bcRegion]
                         ->GetExp(e)
                         ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

        Array<OneD, Array<OneD, NekDouble> > tmp(nVariables);
        for (int i = 0; i < nVariables; ++i)
        {
            tmp[i] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        }

        // pressure perturbation
        Vmath::Fill(nBCEdgePts, m_whiteNoiseBC_p, &tmp[0][0], 1);

        // velocity perturbation
        for (int i = 0; i < nBCEdgePts; ++i)
        {
            NekDouble u = m_whiteNoiseBC_p /
                          sqrt(m_gamma * BfFwd[0][id2 + i] * BfFwd[1][id2 + i]);
            for (int j = 0; j < m_spacedim; ++j)
            {
                tmp[1 + j][i] = -1.0 * u * m_traceNormals[j][id2 + i];
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts,
                         &tmp[i][0], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1], 1);
        }
    }
}



void LEE::v_AuxFields(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<Array<OneD, NekDouble> > &fieldphys,
        std::vector<MultiRegions::ExpListSharedPtr> &expansions,
        std::vector<std::string> &variables)
{
    for (int i = 0; i < m_bfNames.size(); i++)
    {
        fieldphys.push_back(m_bf[i]);

        Array<OneD, NekDouble> tmpC(GetNcoeffs());
        m_fields[0]->FwdTrans(m_bf[i], tmpC);
        fieldcoeffs.push_back(tmpC);

        variables.push_back(m_bfNames[i]);

        expansions.push_back(m_fields[0]);
    }

    int f = 0;
    std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
    for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
    {
        for (int i = 0; i < (*x)->GetForces().num_elements(); ++i)
        {
            fieldphys.push_back((*x)->GetForces()[i]);

            Array<OneD, NekDouble> tmpC(GetNcoeffs());

            m_fields[0]->FwdTrans((*x)->GetForces()[i], tmpC);

            fieldcoeffs.push_back(tmpC);

            variables.push_back("F_" + boost::lexical_cast<string>(f) +
                                "_" + m_session->GetVariable(i));

            expansions.push_back(m_fields[0]);
        }
        f++;
    }

}


/**
 * @brief Get the normal vectors.
 */
const Array<OneD, const Array<OneD, NekDouble> > &LEE::GetNormals()
{
    return m_traceNormals;
}


/**
 * @brief Get the locations of the components of the directed fields within the fields array.
 */
const Array<OneD, const Array<OneD, NekDouble> > &LEE::GetVecLocs()
{
    return m_vecLocs;
}


/**
 * @brief Get the baseflow field.
 */
const Array<OneD, const Array<OneD, NekDouble> > &LEE::GetBasefieldFwdBwd()
{
    return m_bfFwdBwd;
}


void LEE::UpdateBasefieldFwdBwd()
{
    for (int i = 0; i < m_bfNames.size(); i++)
    {
        int j = (m_bfNames.size()) + i;
        m_fields[0]->GetFwdBwdTracePhys(m_bf[i], m_bfFwdBwd[i], m_bfFwdBwd[j]);
        CopyBoundaryTrace(m_bfFwdBwd[i], m_bfFwdBwd[j]);
    }
}


void LEE::CopyBoundaryTrace(const Array<OneD, NekDouble> &Fwd,
                                  Array<OneD, NekDouble> &Bwd)
{
    int cnt = 0;
    // loop over Boundary Regions
    for (int bcRegion = 0;
         bcRegion < m_fields[0]->GetBndConditions().num_elements();
         ++bcRegion)
    {
        // Copy the forward trace of the field to the backward trace
        int e, id2, npts;

        for (e = 0;
             e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
             ++e)
        {
            npts = m_fields[0]
                       ->GetBndCondExpansions()[bcRegion]
                       ->GetExp(e)
                       ->GetTotPoints();
            id2 = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(
                    cnt + e));

            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
        }

        cnt += m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    }
}


/**
 * @brief Get the heat capacity ratio.
 */
NekDouble LEE::GetGamma()
{
    return m_gamma;
}



} //end of namespace

