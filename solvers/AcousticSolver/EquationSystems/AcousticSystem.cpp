///////////////////////////////////////////////////////////////////////////////
//
// File AcousticSystem.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Kilian Lackhove
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: AcousticSystem
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

// Define variable to avoid deprecated warning in Boost 1.69.
#include <boost/version.hpp>
#if BOOST_VERSION >= 106900 && BOOST_VERSION < 107000
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif

#include <boost/core/ignore_unused.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

#include <AcousticSolver/EquationSystems/AcousticSystem.h>

using namespace std;

namespace Nektar
{

AcousticSystem::AcousticSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph), AdvectionSystem(pSession, pGraph),
      m_ip(-1), m_irho(-1), m_iu(1), m_conservative(false)
{
}

/**
 * @brief Initialization object for the AcousticSystem class.
 */
void AcousticSystem::v_InitObject()
{
    AdvectionSystem::v_InitObject();

    ASSERTL0(
        m_projectionType == MultiRegions::eDiscontinuous,
        "Only Projection=DisContinuous supported by the AcousticSystem class.");

    m_bfNames.push_back("c0sq");
    m_bfNames.push_back("rho0");
    m_bfNames.push_back("u0");
    m_bfNames.push_back("v0");
    m_bfNames.push_back("w0");

    // Resize the advection velocities vector to dimension of the problem
    m_bfNames.resize(m_spacedim + 2);

    m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                           m_fields, m_fields.size());

    // Do not forwards transform initial condition
    m_homoInitialFwd = false;

    // Set up locations of velocity and base velocity vectors.
    m_vecLocs    = Array<OneD, Array<OneD, NekDouble>>(1);
    m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        // u', v', w'
        m_vecLocs[0][i] = m_iu + i;
    }

    if (m_session->DefinesElement("Nektar/Coupling"))
    {
        TiXmlElement *vCoupling = m_session->GetElement("Nektar/Coupling");

        ASSERTL0(vCoupling->Attribute("TYPE"),
                 "Missing TYPE attribute in Coupling");
        string vType = vCoupling->Attribute("TYPE");
        ASSERTL0(!vType.empty(),
                 "TYPE attribute must be non-empty in Coupling");

        m_coupling = GetCouplingFactory().CreateInstance(vType, m_fields[0]);
    }

    m_whiteNoiseBC_lastUpdate = -1.0;
    m_whiteNoiseBC_p          = 0.0;
}

/**
 * @brief Destructor for AcousticSystem class.
 */
AcousticSystem::~AcousticSystem()
{
}

/**
 * @brief v_PreIntegrate
 */
bool AcousticSystem::v_PreIntegrate(int step)
{
    GetFunction("Baseflow", m_fields[0], true)
        ->Evaluate(m_bfNames, m_bf, m_time);

    if (m_coupling)
    {
        int numForceFields = 0;
        for (auto &x : m_forcing)
        {
            numForceFields += x->GetForces().size();
        }
        vector<string> varNames;
        Array<OneD, Array<OneD, NekDouble>> phys(
            m_fields.size() + m_bfNames.size() + numForceFields);
        for (int i = 0; i < m_fields.size(); ++i)
        {
            varNames.push_back(m_session->GetVariable(i));
            phys[i] = m_fields[i]->UpdatePhys();
        }
        for (int i = 0; i < m_bfNames.size(); ++i)
        {
            varNames.push_back(m_bfNames[i]);
            phys[m_fields.size() + i] = m_bf[i];
        }

        int f = 0;
        for (auto &x : m_forcing)
        {
            for (int i = 0; i < x->GetForces().size(); ++i)
            {
                phys[m_fields.size() + m_bfNames.size() + f + i] =
                    x->GetForces()[i];
                varNames.push_back("F_" + boost::lexical_cast<string>(f) + "_" +
                                   m_session->GetVariable(i));
            }
            f++;
        }

        m_coupling->Send(step, m_time, phys, varNames);
        m_coupling->Receive(step, m_time, phys, varNames);
    }

    return AdvectionSystem::v_PreIntegrate(step);
}

void AcousticSystem::v_Output()
{
    if (m_coupling)
    {
        m_coupling->Finalize();
    }

    AdvectionSystem::v_Output();
}

/**
 * @brief Compute the right-hand side.
 */
void AcousticSystem::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nVariables = inarray.size();
    int nq         = GetTotPoints();

    // WeakDG does not use advVel, so we only provide a dummy array
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
    m_advection->Advect(nVariables, m_fields, advVel, inarray, outarray, time);

    // Negate the LHS terms
    for (int i = 0; i < nVariables; ++i)
    {
        Vmath::Neg(nq, outarray[i], 1);
    }

    v_AddLinTerm(inarray, outarray);

    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, inarray, outarray, m_time);
    }
}

/**
 * @brief Compute the projection and call the method for imposing the
 * boundary conditions in case of discontinuous projection.
 */
void AcousticSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();
    int nq         = m_fields[0]->GetNpoints();

    // deep copy
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(nq, inarray[i], 1, outarray[i], 1);
    }

    UpdateBasefieldFwdBwd();

    SetBoundaryConditions(outarray, time);
}

/**
 * @brief Apply the Boundary Conditions to the AcousticSystem equations.
 */
void AcousticSystem::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &inarray, NekDouble time)
{
    std::string varName;
    int nvariables = m_fields.size();
    int cnt        = 0;
    int nTracePts  = GetTraceTotPoints();

    // Extract trace for boundaries. Needs to be done on all processors to avoid
    // deadlock.
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }
    Array<OneD, Array<OneD, NekDouble>> bfFwd = GetBasefieldFwdBwd();

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        std::string userDefStr =
            m_fields[0]->GetBndConditions()[n]->GetUserDefined();

        if (!userDefStr.empty())
        {
            // Wall Boundary Condition
            if (boost::iequals(userDefStr, "Wall"))
            {
                v_WallBC(n, cnt, Fwd, inarray);
            }
            else if (boost::iequals(userDefStr, "WhiteNoise"))
            {
                v_WhiteNoiseBC(n, cnt, Fwd, bfFwd, inarray);
            }
            else if (boost::iequals(userDefStr, "RiemannInvariantBC"))
            {
                v_RiemannInvariantBC(n, cnt, Fwd, bfFwd, inarray);
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
 * @brief Wall boundary conditions for the AcousticSystem equations.
 */
void AcousticSystem::v_WallBC(int bcRegion, int cnt,
                              Array<OneD, Array<OneD, NekDouble>> &Fwd,
                              Array<OneD, Array<OneD, NekDouble>> &physarray)
{
    int nVariables = physarray.size();

    const Array<OneD, const int> &traceBndMap = m_fields[0]->GetTraceBndMap();

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int id1, id2, nBCEdgePts;
    int eMax = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();

    for (int e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]
                         ->GetBndCondExpansions()[bcRegion]
                         ->GetExp(e)
                         ->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt + e]);

        // For 2D/3D, define: v* = v - 2(v.n)n
        Array<OneD, NekDouble> tmp(nBCEdgePts, 0.0);

        // Calculate (v.n)
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &Fwd[m_iu + i][id2], 1,
                         &m_traceNormals[i][id2], 1, &tmp[0], 1, &tmp[0], 1);
        }

        // Calculate 2.0(v.n)
        Vmath::Smul(nBCEdgePts, -2.0, &tmp[0], 1, &tmp[0], 1);

        // Calculate v* = v - 2.0(v.n)n
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(nBCEdgePts, &tmp[0], 1, &m_traceNormals[i][id2], 1,
                         &Fwd[m_iu + i][id2], 1, &Fwd[m_iu + i][id2], 1);
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

/**
 * @brief Wall boundary conditions for the AcousticSystem equations.
 */
void AcousticSystem::v_WhiteNoiseBC(
    int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd,
    Array<OneD, Array<OneD, NekDouble>> &BfFwd,
    Array<OneD, Array<OneD, NekDouble>> &physarray)
{
    boost::ignore_unused(Fwd);

    int id1, id2, nBCEdgePts;
    int nVariables = physarray.size();

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
        std::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(
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

        Array<OneD, Array<OneD, NekDouble>> tmp(nVariables);
        for (int i = 0; i < nVariables; ++i)
        {
            tmp[i] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        }

        // pressure perturbation
        Vmath::Fill(nBCEdgePts, m_whiteNoiseBC_p, &tmp[m_ip][0], 1);

        if (m_conservative)
        {
            for (int i = 0; i < nBCEdgePts; ++i)
            {
                // density perturbation
                tmp[m_irho][i] = m_whiteNoiseBC_p *
                                 BfFwd[m_spacedim + 2][id2 + i] /
                                 BfFwd[0][id2 + i];

                // velocity perturbation
                NekDouble ru = m_whiteNoiseBC_p / sqrt(BfFwd[0][id2 + i]);
                for (int j = 0; j < m_spacedim; ++j)
                {
                    tmp[m_iu + j][i] = -1.0 * ru * m_traceNormals[j][id2 + i];
                }
            }
        }
        else
        {
            for (int i = 0; i < nBCEdgePts; ++i)
            {
                // velocity perturbation
                NekDouble u = m_whiteNoiseBC_p /
                              (sqrt(BfFwd[0][id2 + i]) * BfFwd[1][id2 + i]);

                for (int j = 0; j < m_spacedim; ++j)
                {
                    tmp[m_iu + j][i] = -1.0 * u * m_traceNormals[j][id2 + i];
                }
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (int i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts, &tmp[i][0], 1,
                         &(m_fields[i]
                               ->GetBndCondExpansions()[bcRegion]
                               ->UpdatePhys())[id1],
                         1);
        }
    }
}

/**
 * @brief Compute the advection velocity in the standard space
 * for each element of the expansion.
 *
 * @return       Standard velocity field.
 */
Array<OneD, NekDouble> AcousticSystem::v_GetMaxStdVelocity(void)
{
    int nElm = m_fields[0]->GetExpSize();

    Array<OneD, NekDouble> stdV(nElm, 0.0);

    Array<OneD, Array<OneD, NekDouble>> stdVelocity(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> velocity(m_spacedim);
    LibUtilities::PointsKeyVector ptsKeys;

    int cnt = 0;

    for (int el = 0; el < nElm; ++el)
    {
        ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

        // Possible bug: not multiply by jacobian??
        const SpatialDomains::GeomFactorsSharedPtr metricInfo =
            m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
        const Array<TwoD, const NekDouble> &gmat =
            m_fields[0]
                ->GetExp(el)
                ->GetGeom()
                ->GetMetricInfo()
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
                NekDouble c    = sqrt(m_bf[0][cnt + j]);
                velocity[i][j] = m_bf[i + 2][cnt + j] + c;
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

    return stdV;
}

void AcousticSystem::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    for (int i = 0; i < m_bfNames.size(); i++)
    {
        variables.push_back(m_bfNames[i]);
        Array<OneD, NekDouble> tmpC(GetNcoeffs());
        m_fields[0]->FwdTrans(m_bf[i], tmpC);
        fieldcoeffs.push_back(tmpC);
    }

    int f = 0;
    for (auto &x : m_forcing)
    {
        for (int i = 0; i < x->GetForces().size(); ++i)
        {
            variables.push_back("F_" + boost::lexical_cast<string>(f) + "_" +
                                m_session->GetVariable(i));
            Array<OneD, NekDouble> tmpC(GetNcoeffs());
            m_fields[0]->FwdTrans(x->GetForces()[i], tmpC);
            fieldcoeffs.push_back(tmpC);
        }
        f++;
    }
}

/**
 * @brief Get the normal vectors.
 */
const Array<OneD, const Array<OneD, NekDouble>> &AcousticSystem::GetNormals()
{
    return m_traceNormals;
}

/**
 * @brief Get the locations of the components of the directed fields within the
 * fields array.
 */
const Array<OneD, const Array<OneD, NekDouble>> &AcousticSystem::GetVecLocs()
{
    return m_vecLocs;
}

/**
 * @brief Get the baseflow field.
 */
const Array<OneD, const Array<OneD, NekDouble>>
    &AcousticSystem::GetBasefieldFwdBwd()
{
    return m_bfFwdBwd;
}

void AcousticSystem::UpdateBasefieldFwdBwd()
{
    for (int i = 0; i < m_bfNames.size(); i++)
    {
        int j = m_bfNames.size() + i;
        m_fields[0]->GetFwdBwdTracePhys(m_bf[i], m_bfFwdBwd[i], m_bfFwdBwd[j]);
        CopyBoundaryTrace(m_bfFwdBwd[i], m_bfFwdBwd[j]);
    }
}

void AcousticSystem::CopyBoundaryTrace(const Array<OneD, NekDouble> &Fwd,
                                       Array<OneD, NekDouble> &Bwd)
{
    int cnt = 0;
    // loop over Boundary Regions
    for (int bcRegion = 0;
         bcRegion < m_fields[0]->GetBndConditions().size(); ++bcRegion)
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
                 m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(
                    cnt + e));

            Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
        }

        cnt += m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
    }
}

} // namespace Nektar
