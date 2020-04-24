///////////////////////////////////////////////////////////////////////////////
//
// File NonlinearPeregrine.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Nonlinear Boussinesq equations of Peregrine in
//              conservative variables (constant depth case)
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/NonlinearPeregrine.h>

using namespace std;

namespace Nektar
{

string NonlinearPeregrine::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
                "NonlinearPeregrine", NonlinearPeregrine::create,
                "Nonlinear Peregrine equations in conservative variables.");

NonlinearPeregrine::NonlinearPeregrine(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : ShallowWaterSystem(pSession, pGraph), m_factors()
{
    m_factors[StdRegions::eFactorLambda] = 0.0;
    m_factors[StdRegions::eFactorTau] = 1000000.0;
    // note: eFactorTau = 1.0 becomes unstable...
    // we need to investigate the behaviuor w.r.t. tau
}

void NonlinearPeregrine::v_InitObject()
{
    ShallowWaterSystem::v_InitObject();

    if (m_session->DefinesSolverInfo("PROBLEMTYPE"))
    {
        int i;
        std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
        for (i = 0; i < (int) SIZE_ProblemType; ++i)
        {
            if (boost::iequals(ProblemTypeMap[i], ProblemTypeStr))
            {
                m_problemType = (ProblemType) i;
                break;
            }
        }
    }
    else
    {
        m_problemType = (ProblemType) 0;
    }

    if (m_explicitAdvection)
    {
        m_ode.DefineOdeRhs(&NonlinearPeregrine::DoOdeRhs, this);
        m_ode.DefineProjection(&NonlinearPeregrine::DoOdeProjection, this);
    }
    else
    {
        ASSERTL0(false, "Implicit Peregrine not set up.");
    }

    // NB! At the moment only the constant depth case is
    // supported for the Peregrine eq.
    if (m_session->DefinesParameter("ConstDepth"))
    {
        m_const_depth = m_session->GetParameter("ConstDepth");
    }
    else
    {
        ASSERTL0(false, "Constant Depth not specified");
    }

    // Type of advection class to be used
    switch (m_projectionType)
    {
        // Continuous field
        case MultiRegions::eGalerkin:
        {
            ASSERTL0(false,
                    "Continuous projection type not supported for Peregrine.");
            break;
        }
        // Discontinuous field
        case MultiRegions::eDiscontinuous:
        {
            string advName;
            string diffName;
            string riemName;

            //---------------------------------------------------------------
            // Setting up advection and diffusion operators
            // NB: diffusion not set up for SWE at the moment
            //     but kept here for future use ...
            m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
            // m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
            m_advection = SolverUtils::GetAdvectionFactory().CreateInstance(
                    advName, advName);

            m_advection->SetFluxVector(&NonlinearPeregrine::GetFluxVector,
                    this);

            // Setting up Riemann solver for advection operator
            m_session->LoadSolverInfo("UpwindType", riemName, "NoSolver");

            m_riemannSolver =
                    SolverUtils::GetRiemannSolverFactory().CreateInstance(
                            riemName, m_session);

            // Setting up parameters for advection operator Riemann solver
            m_riemannSolver->SetParam("gravity",
                    &NonlinearPeregrine::GetGravity, this);
            m_riemannSolver->SetAuxVec("vecLocs",
                    &NonlinearPeregrine::GetVecLocs, this);
            m_riemannSolver->SetVector("N", &NonlinearPeregrine::GetNormals,
                    this);
            m_riemannSolver->SetScalar("depth", &NonlinearPeregrine::GetDepth,
                    this);

            // Concluding initialisation of advection / diffusion operators
            m_advection->SetRiemannSolver(m_riemannSolver);
            m_advection->InitObject(m_session, m_fields);
            break;
        }
        default:
        {
            ASSERTL0(false, "Unsupported projection type.");
            break;
        }
    }

}

NonlinearPeregrine::~NonlinearPeregrine()
{

}

// physarray contains the conservative variables
void NonlinearPeregrine::AddCoriolis(
        const Array<OneD, const Array<OneD, NekDouble> > &physarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray)
{

    int ncoeffs = GetNcoeffs();
    int nq = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // add to hu equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to hv equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            m_fields[0]->IProductWRTBase(tmp, mod);
            m_fields[0]->MultiplyByElmtInvMass(mod, mod);
            m_fields[0]->BwdTrans(mod, tmp);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            // add to hu equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[2], 1, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[1], 1, outarray[1], 1);

            // add to hv equation
            Vmath::Vmul(nq, m_coriolis, 1, physarray[1], 1, tmp, 1);
            Vmath::Neg(nq, tmp, 1);
            Vmath::Vadd(nq, tmp, 1, outarray[2], 1, outarray[2], 1);
            break;
        }
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }

}

// physarray contains the conservative variables
void NonlinearPeregrine::AddVariableDepth(
        const Array<OneD, const Array<OneD, NekDouble> > &physarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray)
{

    int ncoeffs = GetNcoeffs();
    int nq = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, m_bottomSlope[i], 1, physarray[0], 1, tmp, 1);
                Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);
                m_fields[0]->IProductWRTBase(tmp, mod);
                m_fields[0]->MultiplyByElmtInvMass(mod, mod);
                m_fields[0]->BwdTrans(mod, tmp);
                Vmath::Vadd(nq, tmp, 1, outarray[i + 1], 1, outarray[i + 1], 1);
            }
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            for (int i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vmul(nq, m_bottomSlope[i], 1, physarray[0], 1, tmp, 1);
                Vmath::Smul(nq, m_g, tmp, 1, tmp, 1);
                Vmath::Vadd(nq, tmp, 1, outarray[i + 1], 1, outarray[i + 1], 1);
            }
            break;
        }
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }

}

void NonlinearPeregrine::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
        Array<OneD, Array<OneD, NekDouble> >&outarray, const NekDouble time)
{
    int i;
    int nvariables = inarray.size();
    int ncoeffs = GetNcoeffs();
    int nq = GetTotPoints();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {

            //-------------------------------------------------------
            //inarray in physical space

            Array<OneD, Array<OneD, NekDouble> > modarray(nvariables);
            for (i = 0; i < nvariables; ++i)
            {
                modarray[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
            }
            //-------------------------------------------------------

            //-------------------------------------------------------
            // Compute the DG advection including the numerical flux
            // by using SolverUtils/Advection
            // Input and output in physical space
            Array<OneD, Array<OneD, NekDouble> > advVel;

            m_advection->Advect(nvariables - 1, m_fields, advVel, inarray,
                    outarray, time);
            //-------------------------------------------------------

            //-------------------------------------------------------
            // negate the outarray since moving terms to the rhs
            for (i = 0; i < nvariables - 1; ++i)
            {
                Vmath::Neg(nq, outarray[i], 1);
            }
            //-------------------------------------------------------

            //-------------------------------------------------
            // Add "source terms"
            // Input and output in physical space

            // Coriolis forcing
            if (m_coriolis.size() != 0)
            {
                AddCoriolis(inarray, outarray);
            }

            // Variable Depth
            if (m_constantDepth != true)
            {
                ASSERTL0(false,
                        "Variable depth not supported for the Peregrine "
                        "equations");
            }

            //-------------------------------------------------

            //---------------------------------------
            // As no more terms is required for the
            // continuity equation and we have aleady evaluated
            //  the values for h_t we are done for h
            //---------------------------------------

            //-------------------------------------------------
            // go to modal space
            m_fields[0]->IProductWRTBase(outarray[1], modarray[1]);
            m_fields[0]->IProductWRTBase(outarray[2], modarray[2]);

            // store f1 and f2 for later use (modal space)
            Array<OneD, NekDouble> f1(ncoeffs);
            Array<OneD, NekDouble> f2(ncoeffs);

            Vmath::Vcopy(ncoeffs, modarray[1], 1, f1, 1); // f1
            Vmath::Vcopy(ncoeffs, modarray[2], 1, f2, 1); // f2

            // Solve the remaining block-diagonal systems
            m_fields[0]->MultiplyByElmtInvMass(modarray[1], modarray[1]);
            m_fields[0]->MultiplyByElmtInvMass(modarray[2], modarray[2]);
            //---------------------------------------------

            //---------------------------------------------

            //-------------------------------------------------
            // create tmp fields to be used during
            // the dispersive section

            Array<OneD, Array<OneD, NekDouble> > coeffsfield(2);
            Array<OneD, Array<OneD, NekDouble> > physfield(2);

            for (i = 0; i < 2; ++i)
            {
                coeffsfield[i] = Array<OneD, NekDouble>(ncoeffs);
                physfield[i] = Array<OneD, NekDouble>(nq);
            }
            //---------------------------------------------

            //---------------------------------------------
            // Go from modal to physical space
            Vmath::Vcopy(nq, outarray[1], 1, physfield[0], 1);
            Vmath::Vcopy(nq, outarray[2], 1, physfield[1], 1);
            //---------------------------------------

            //---------------------------------------
            // Start for solve of mixed dispersive terms
            // using the 'WCE method'
            // (Eskilsson & Sherwin, JCP 2006)

            // constant depth case
            // \nabla \cdot (\nabla z) - invgamma z
            //    = - invgamma (\nabla \cdot {\bf f}_(2,3)

            NekDouble gamma = (m_const_depth * m_const_depth) * (1.0 / 3.0);
            NekDouble invgamma = 1.0 / gamma;

            int nTraceNumPoints = GetTraceTotPoints();
            Array<OneD, Array<OneD, NekDouble> > upwindX(1);
            Array<OneD, Array<OneD, NekDouble> > upwindY(1);
            upwindX[0] = Array<OneD, NekDouble>(nTraceNumPoints);
            upwindY[0] = Array<OneD, NekDouble>(nTraceNumPoints);
            //--------------------------------------------

            //--------------------------------------------
            // Compute the forcing function for the
            // wave continuity equation

            // Set boundary condidtions for z
            SetBoundaryConditionsForcing(physfield, time);

            // \nabla \phi \cdot f_{2,3}
            m_fields[0]->IProductWRTDerivBase(0, physfield[0], coeffsfield[0]);
            m_fields[0]->IProductWRTDerivBase(1, physfield[1], coeffsfield[1]);
            Vmath::Vadd(ncoeffs, coeffsfield[0], 1, coeffsfield[1], 1,
                                 coeffsfield[0], 1);
            Vmath::Neg(ncoeffs, coeffsfield[0], 1);

            // Evaluate  upwind numerical flux (physical space)
            NumericalFluxForcing(physfield, upwindX[0], upwindY[0]);

            m_fields[0]->AddTraceIntegral(upwindX[0], upwindY[0],
                                          coeffsfield[0]);
            m_fields[0]->MultiplyByElmtInvMass(coeffsfield[0], coeffsfield[0]);
            m_fields[0]->BwdTrans(coeffsfield[0], physfield[0]);

            Vmath::Smul(nq, -invgamma, physfield[0], 1, physfield[0], 1);

            // ok: forcing function for HelmSolve... done!
            //--------------------------------------

            //--------------------------------------
            // Solve the Helmhotz-type equation
            // for the wave continuity equation
            // (missing slope terms...)

            // note: this is just valid for the constant depth case:

            // as of now we need not to specify any
            // BC routine for the WCE: periodic
            // and zero Neumann (for walls)

            WCESolve(physfield[0], invgamma);

            Vmath::Vcopy(nq, physfield[0], 1, outarray[3], 1); // store z

            // ok: Wave Continuity Equation... done!
            //------------------------------------

            //------------------------------------
            // Return to the primary variables

            // (h {\bf u})_t = gamma \nabla z + {\bf f}_{2,3}

            Vmath::Smul(nq, gamma, physfield[0], 1, physfield[0], 1);

            // Set boundary conditions
            SetBoundaryConditionsContVariables(physfield[0], time);

            m_fields[0]->IProductWRTDerivBase(0, physfield[0], coeffsfield[0]);
            m_fields[1]->IProductWRTDerivBase(1, physfield[0], coeffsfield[1]);

            Vmath::Neg(ncoeffs, coeffsfield[0], 1);
            Vmath::Neg(ncoeffs, coeffsfield[1], 1);

            // Evaluate  upwind numerical flux (physical space)
            NumericalFluxConsVariables(physfield[0], upwindX[0], upwindY[0]);

            {
                Array<OneD, NekDouble> uptemp(nTraceNumPoints, 0.0);

                m_fields[0]->AddTraceIntegral(upwindX[0], uptemp,
                                              coeffsfield[0]);
                m_fields[0]->AddTraceIntegral(uptemp, upwindY[0],
                                              coeffsfield[1]);
            }

            Vmath::Vadd(ncoeffs, f1, 1, coeffsfield[0], 1, modarray[1], 1);
            Vmath::Vadd(ncoeffs, f2, 1, coeffsfield[1], 1, modarray[2], 1);

            m_fields[1]->MultiplyByElmtInvMass(modarray[1], modarray[1]);
            m_fields[2]->MultiplyByElmtInvMass(modarray[2], modarray[2]);

            m_fields[1]->BwdTrans(modarray[1], outarray[1]);
            m_fields[2]->BwdTrans(modarray[2], outarray[2]);

            // ok: returned to conservative variables... done!
            //---------------------

            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
            ASSERTL0(false, "Unknown projection scheme for the Peregrine");
            break;
        default:
            ASSERTL0(false, "Unknown projection scheme for the NonlinearSWE");
            break;
    }
}

void NonlinearPeregrine::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
        Array<OneD, Array<OneD, NekDouble> >&outarray,
        const NekDouble time)
{
    int i;
    int nvariables = inarray.size();

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {

            // Just copy over array
            int npoints = GetNpoints();

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }

            SetBoundaryConditions(outarray, time);
            break;
        }
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {

            EquationSystem::SetBoundaryConditions(time);
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(),0.0);

            for (i = 0; i < nvariables; ++i)
            {
                m_fields[i]->FwdTrans(inarray[i], coeffs);
                m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
            }
            break;
        }
        default:
            ASSERTL0(false, "Unknown projection scheme");
            break;
    }
}

//----------------------------------------------------
void NonlinearPeregrine::SetBoundaryConditions(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble time)
{

    int nvariables = m_fields.size();
    int cnt = 0;
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
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {

        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"Wall"))
        {
            WallBoundary2D(n, cnt, Fwd, inarray);
        }

        // Time Dependent Boundary Condition (specified in meshfile)
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            for (int i = 0; i < nvariables; ++i)
            {
                m_fields[i]->EvaluateBoundaryConditions(time);
            }
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

//----------------------------------------------------
/**
 * @brief Wall boundary condition.
 */
void NonlinearPeregrine::WallBoundary(int bcRegion, int cnt,
        Array<OneD, Array<OneD, NekDouble> > &Fwd,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    int i;
    int nvariables = physarray.size();

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
                                m_fields[0]->GetBndCondExpansions()[bcRegion];
    for (e = 0; e < bcexp->GetExpSize(); ++e)
    {
        npts = bcexp->GetExp(e)->GetTotPoints();
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
               m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        // For 2D/3D, define: v* = v - 2(v.n)n
        Array<OneD, NekDouble> tmp(npts, 0.0);

        // Calculate (v.n)
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, &Fwd[1 + i][id2], 1, &m_traceNormals[i][id2], 1,
                               &tmp[0],          1, &tmp[0],                 1);
        }

        // Calculate 2.0(v.n)
        Vmath::Smul(npts, -2.0, &tmp[0], 1, &tmp[0], 1);

        // Calculate v* = v - 2.0(v.n)n
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vvtvp(npts, &tmp[0],          1, &m_traceNormals[i][id2], 1,
                               &Fwd[1 + i][id2], 1, &Fwd[1 + i][id2],        1);
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            bcexp = m_fields[i]->GetBndCondExpansions()[bcRegion];
            Vmath::Vcopy(npts, &Fwd[i][id2], 1, &(bcexp->UpdatePhys())[id1], 1);
        }
    }
}

void NonlinearPeregrine::WallBoundary2D(
        int bcRegion,
        int cnt,
        Array<OneD, Array<OneD, NekDouble> > &Fwd,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
{
    boost::ignore_unused(physarray);

    int i;
    int nvariables = 3;

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
                                m_fields[0]->GetBndCondExpansions()[bcRegion];

    for (e = 0; e < bcexp->GetExpSize();
            ++e)
    {
        npts = bcexp->GetExp(e)->GetNumPoints(0);
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
               m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                // negate the forward flux
                Vmath::Neg(npts, &Fwd[1][id2], 1);
                break;
            }
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul (npts, &Fwd[1][id2], 1, &m_traceNormals[0][id2], 1,
                                   &tmp_n[0],    1);
                Vmath::Vvtvp(npts, &Fwd[2][id2], 1, &m_traceNormals[1][id2], 1,
                                   &tmp_n[0],    1, &tmp_n[0],               1);

                Vmath::Vmul (npts, &Fwd[1][id2], 1, &m_traceNormals[1][id2], 1,
                                   &tmp_t[0],    1);
                Vmath::Vvtvm(npts, &Fwd[2][id2], 1, &m_traceNormals[0][id2], 1,
                                   &tmp_t[0],    1, &tmp_t[0],               1);

                // negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                // rotate back to Cartesian
                Vmath::Vmul (npts, &tmp_t[0],    1, &m_traceNormals[1][id2], 1,
                                   &Fwd[1][id2], 1);
                Vmath::Vvtvm(npts, &tmp_n[0],    1, &m_traceNormals[0][id2], 1,
                                   &Fwd[1][id2], 1, &Fwd[1][id2],            1);

                Vmath::Vmul(npts,  &tmp_t[0],    1, &m_traceNormals[0][id2], 1,
                                   &Fwd[2][id2], 1);
                Vmath::Vvtvp(npts, &tmp_n[0],    1, &m_traceNormals[1][id2], 1,
                                   &Fwd[2][id2], 1, &Fwd[2][id2],            1);
                break;
            }
            case 3:
                ASSERTL0(false,
                        "3D not implemented for Shallow Water Equations");
                break;
            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nvariables; ++i)
        {
            bcexp = m_fields[i]->GetBndCondExpansions()[bcRegion];
            Vmath::Vcopy(npts, &Fwd[i][id2], 1, &(bcexp->UpdatePhys())[id1], 1);
        }
    }
}

// Physfield in conservative Form
void NonlinearPeregrine::GetFluxVector(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
{
    int i, j;
    int nq = m_fields[0]->GetTotPoints();

    NekDouble g = m_g;
    Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

    // Flux vector for the mass equation
    for (i = 0; i < m_spacedim; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, physfield[i + 1], 1, flux[0][i], 1);
    }

    GetVelocityVector(physfield, velocity);

    // Put (0.5 g h h) in tmp
    Array<OneD, NekDouble> tmp(nq);
    Vmath::Vmul(nq, physfield[0], 1, physfield[0], 1, tmp, 1);
    Vmath::Smul(nq, 0.5 * g, tmp, 1, tmp, 1);

    // Flux vector for the momentum equations
    for (i = 0; i < m_spacedim; ++i)
    {
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j],    1, physfield[i + 1], 1,
                            flux[i + 1][j], 1);
        }

        // Add (0.5 g h h) to appropriate field
        Vmath::Vadd(nq, flux[i + 1][i], 1, tmp, 1, flux[i + 1][i], 1);
    }

}

void NonlinearPeregrine::ConservativeToPrimitive(
        const Array<OneD, const Array<OneD, NekDouble> >&physin,
        Array<OneD, Array<OneD, NekDouble> >&physout)
{
    int nq = GetTotPoints();

    if (physin.get() == physout.get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble> > tmp(3);
        for (int i = 0; i < 3; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physin[i], 1, tmp[i], 1);
        }

        // \eta = h - d
        Vmath::Vsub(nq, tmp[0], 1, m_depth, 1, physout[0], 1);

        // u = hu/h
        Vmath::Vdiv(nq, tmp[1], 1, tmp[0], 1, physout[1], 1);

        // v = hv/ v
        Vmath::Vdiv(nq, tmp[2], 1, tmp[0], 1, physout[2], 1);
    }
    else
    {
        // \eta = h - d
        Vmath::Vsub(nq, physin[0], 1, m_depth, 1, physout[0], 1);

        // u = hu/h
        Vmath::Vdiv(nq, physin[1], 1, physin[0], 1, physout[1], 1);

        // v = hv/ v
        Vmath::Vdiv(nq, physin[2], 1, physin[0], 1, physout[2], 1);
    }
}

void NonlinearPeregrine::v_ConservativeToPrimitive()
{
    int nq = GetTotPoints();

    // u = hu/h
    Vmath::Vdiv(nq, m_fields[1]->GetPhys(),    1, m_fields[0]->GetPhys(), 1,
                    m_fields[1]->UpdatePhys(), 1);

    // v = hv/ v
    Vmath::Vdiv(nq, m_fields[2]->GetPhys(),    1, m_fields[0]->GetPhys(), 1,
                    m_fields[2]->UpdatePhys(), 1);

    // \eta = h - d
    Vmath::Vsub(nq, m_fields[0]->GetPhys(),    1, m_depth, 1,
                    m_fields[0]->UpdatePhys(), 1);
}

void NonlinearPeregrine::PrimitiveToConservative(
        const Array<OneD, const Array<OneD, NekDouble> >&physin,
        Array<OneD, Array<OneD, NekDouble> >&physout)
{

    int nq = GetTotPoints();

    if (physin.get() == physout.get())
    {
        // copy indata and work with tmp array
        Array<OneD, Array<OneD, NekDouble> > tmp(3);
        for (int i = 0; i < 3; ++i)
        {
            // deep copy
            tmp[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physin[i], 1, tmp[i], 1);
        }

        // h = \eta + d
        Vmath::Vadd(nq, tmp[0],     1, m_depth, 1, physout[0], 1);

        // hu = h * u
        Vmath::Vmul(nq, physout[0], 1, tmp[1],  1, physout[1], 1);

        // hv = h * v
        Vmath::Vmul(nq, physout[0], 1, tmp[2],  1, physout[2], 1);

    }
    else
    {
        // h = \eta + d
        Vmath::Vadd(nq, physin[0],  1, m_depth,   1, physout[0], 1);

        // hu = h * u
        Vmath::Vmul(nq, physout[0], 1, physin[1], 1, physout[1], 1);

        // hv = h * v
        Vmath::Vmul(nq, physout[0], 1, physin[2], 1, physout[2], 1);

    }

}

void NonlinearPeregrine::v_PrimitiveToConservative()
{
    int nq = GetTotPoints();

    // h = \eta + d
    Vmath::Vadd(nq, m_fields[0]->GetPhys(),    1, m_depth, 1,
                    m_fields[0]->UpdatePhys(), 1);

    // hu = h * u
    Vmath::Vmul(nq, m_fields[0]->GetPhys(),    1, m_fields[1]->GetPhys(), 1,
                    m_fields[1]->UpdatePhys(), 1);

    // hv = h * v
    Vmath::Vmul(nq, m_fields[0]->GetPhys(),    1, m_fields[2]->GetPhys(), 1,
                    m_fields[2]->UpdatePhys(), 1);
}

/**
 * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
 * \f$ h\mathbf{v} \f$.
 *
 * @param physfield  Momentum field.
 * @param velocity   Velocity field.
 */
void NonlinearPeregrine::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD, Array<OneD, NekDouble> > &velocity)
{
    const int npts = physfield[0].size();

    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Vdiv(npts, physfield[1 + i], 1, physfield[0], 1, velocity[i], 1);
    }
}

void NonlinearPeregrine::v_GenerateSummary(SolverUtils::SummaryList& s)
{
    ShallowWaterSystem::v_GenerateSummary(s);
    SolverUtils::AddSummaryItem(s, "Variables", "h  should be in field[0]");
    SolverUtils::AddSummaryItem(s, "", "hu should be in field[1]");
    SolverUtils::AddSummaryItem(s, "", "hv should be in field[2]");
    SolverUtils::AddSummaryItem(s, "", "z  should be in field[3]");
}

void NonlinearPeregrine::WCESolve(
        Array<OneD, NekDouble> &fce,
        NekDouble lambda)
{
    int nq = GetTotPoints();

    m_factors[StdRegions::eFactorLambda] = lambda;

    for (int j = 0; j < nq; j++)
    {
        (m_fields[3]->UpdatePhys())[j] = fce[j];
    }

    m_fields[3]->SetPhysState(true);

    m_fields[3]->HelmSolve(m_fields[3]->GetPhys(),
                           m_fields[3]->UpdateCoeffs(),
                           m_factors);

    m_fields[3]->BwdTrans(m_fields[3]->GetCoeffs(), m_fields[3]->UpdatePhys());

    m_fields[3]->SetPhysState(true);

    Vmath::Vcopy(nq, m_fields[3]->GetPhys(), 1, fce, 1);
}

void NonlinearPeregrine::NumericalFluxForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, NekDouble> &numfluxX,
        Array<OneD, NekDouble> &numfluxY)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();

    //-----------------------------------------------------
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(2);
    Array<OneD, Array<OneD, NekDouble> > Bwd(2);

    for (i = 0; i < 2; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        Bwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
    }
    //-----------------------------------------------------

    //-----------------------------------------------------
    // get the physical values at the trace
    // (any time-dependent BC previuosly put in fields[1] and [2]

    m_fields[1]->GetFwdBwdTracePhys(inarray[0], Fwd[0], Bwd[0]);
    m_fields[2]->GetFwdBwdTracePhys(inarray[1], Fwd[1], Bwd[1]);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // use centred fluxes for the numerical flux
    for (i = 0; i < nTraceNumPoints; ++i)
    {
        numfluxX[i] = 0.5 * (Fwd[0][i] + Bwd[0][i]);
        numfluxY[i] = 0.5 * (Fwd[1][i] + Bwd[1][i]);
    }
    //-----------------------------------------------------
}

void NonlinearPeregrine::SetBoundaryConditionsForcing(
        Array<OneD, Array<OneD, NekDouble> > &inarray,
        NekDouble time)
{
    boost::ignore_unused(time);

    int cnt = 0;

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // Use wall for all BC...
        // Wall Boundary Condition
        if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"Wall"))
        {
            WallBoundaryForcing(n, cnt, inarray);
        }

        //Timedependent Boundary Condition
        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            ASSERTL0(false, "time-dependent BC not implemented for Boussinesq");
        }
        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
    }
}

// fills up boundary expansion for field[1] and [2]
void NonlinearPeregrine::WallBoundaryForcing(
        int bcRegion,
        int cnt,
        Array<OneD, Array<OneD, NekDouble> >&inarray)
{

    //std::cout << " WallBoundaryForcing" << std::endl;

    int nTraceNumPoints = GetTraceTotPoints();
    int nvariables = 2;

    // get physical values of f1 and f2 for the forward trace
    Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
        m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
    }

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
                                m_fields[0]->GetBndCondExpansions()[bcRegion];
    for (e = 0; e < bcexp->GetExpSize(); ++e)
    {
        npts = bcexp->GetExp(e)->GetTotPoints();
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
               m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        switch (m_expdim)
        {
            case 1:
            {
                ASSERTL0(false, "1D not yet implemented for Boussinesq");
                break;
            }
            case 2:
            {
                Array<OneD, NekDouble> tmp_n(npts);
                Array<OneD, NekDouble> tmp_t(npts);

                Vmath::Vmul (npts, &Fwd[0][id2], 1, &m_traceNormals[0][id2], 1,
                                   &tmp_n[0],    1);
                Vmath::Vvtvp(npts, &Fwd[1][id2], 1, &m_traceNormals[1][id2], 1,
                                   &tmp_n[0],    1, &tmp_n[0],               1);

                Vmath::Vmul (npts, &Fwd[0][id2], 1, &m_traceNormals[1][id2], 1,
                                   &tmp_t[0],    1);
                Vmath::Vvtvm(npts, &Fwd[1][id2], 1, &m_traceNormals[0][id2], 1,
                                   &tmp_t[0],    1, &tmp_t[0],               1);

                // negate the normal flux
                Vmath::Neg(npts, tmp_n, 1);

                // rotate back to Cartesian
                Vmath::Vmul (npts, &tmp_t[0],    1, &m_traceNormals[1][id2], 1,
                                   &Fwd[0][id2], 1);
                Vmath::Vvtvm(npts, &tmp_n[0],    1, &m_traceNormals[0][id2], 1,
                                   &Fwd[0][id2], 1, &Fwd[0][id2],            1);

                Vmath::Vmul (npts, &tmp_t[0],    1, &m_traceNormals[0][id2], 1,
                                   &Fwd[1][id2], 1);
                Vmath::Vvtvp(npts, &tmp_n[0],    1, &m_traceNormals[1][id2], 1,
                                   &Fwd[1][id2], 1, &Fwd[1][id2],            1);
                break;
            }
            case 3:
                ASSERTL0(false, "3D not implemented for Boussinesq equations");
                break;
            default:
                ASSERTL0(false, "Illegal expansion dimension");
        }

        // copy boundary adjusted values into the boundary expansion
        bcexp = m_fields[1]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &Fwd[0][id2], 1, &(bcexp->UpdatePhys())[id1], 1);

        bcexp = m_fields[2]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &Fwd[1][id2], 1, &(bcexp->UpdatePhys())[id1], 1);
    }
}

void NonlinearPeregrine::SetBoundaryConditionsContVariables(
        Array<OneD, NekDouble> &inarray,
        NekDouble time)
{
    boost::ignore_unused(time);

    int cnt = 0;

    // loop over Boundary Regions
    for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
    {
        // Use wall for all
        // Wall Boundary Condition
        if(boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"Wall"))
        {
            WallBoundaryContVariables(n, cnt, inarray);
        }

        if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
        {
            WallBoundaryContVariables(n, cnt, inarray);
        }

        cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize() - 1;
    }
}

void NonlinearPeregrine::WallBoundaryContVariables(
        int bcRegion,
        int cnt,
        Array<OneD, NekDouble>&inarray)
{
    int nTraceNumPoints = GetTraceTotPoints();

    // get physical values of z for the forward trace
    Array<OneD, NekDouble> z(nTraceNumPoints);
    m_fields[0]->ExtractTracePhys(inarray, z);

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;
    MultiRegions::ExpListSharedPtr bcexp =
                                m_fields[0]->GetBndCondExpansions()[bcRegion];

    for (e = 0; e < bcexp->GetExpSize(); ++e)
    {
        npts = bcexp->GetExp(e)->GetTotPoints();
        id1  = bcexp->GetPhys_Offset(e);
        id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
               m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt + e));

        // copy boundary adjusted values into the boundary expansion
        // field[1] and field[2]
        bcexp = m_fields[1]->GetBndCondExpansions()[bcRegion];
        Vmath::Vcopy(npts, &z[id2], 1, &(bcexp->UpdatePhys())[id1], 1);

    }
}

void NonlinearPeregrine::NumericalFluxConsVariables(
        Array<OneD, NekDouble> &physfield, Array<OneD, NekDouble> &outX,
        Array<OneD, NekDouble> &outY)
{
    int i;
    int nTraceNumPoints = GetTraceTotPoints();

    //-----------------------------------------------------
    // get temporary arrays
    Array<OneD, Array<OneD, NekDouble> > Fwd(1);
    Array<OneD, Array<OneD, NekDouble> > Bwd(1);

    Fwd[0] = Array<OneD, NekDouble>(nTraceNumPoints);
    Bwd[0] = Array<OneD, NekDouble>(nTraceNumPoints);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // get the physical values at the trace
    // (we have put any time-dependent BC in field[1])

    m_fields[1]->GetFwdBwdTracePhys(physfield, Fwd[0], Bwd[0]);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // use centred fluxes for the numerical flux
    for (i = 0; i < nTraceNumPoints; ++i)
    {
        outX[i] = 0.5 * (Fwd[0][i] + Bwd[0][i]);
        outY[i] = 0.5 * (Fwd[0][i] + Bwd[0][i]);
    }
    //-----------------------------------------------------
}

// initial condition Laitone's first order solitary wave
void NonlinearPeregrine::LaitoneSolitaryWave(
        NekDouble amp,
        NekDouble d,
        NekDouble time,
        NekDouble x_offset)
{
    int nq = GetTotPoints();

    NekDouble A = 1.0;
    NekDouble C = sqrt(m_g * d) * (1.0 + 0.5 * (amp / d));

    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> x1(nq);
    Array<OneD, NekDouble> zeros(nq, 0.0);

    // get the coordinates (assuming all fields have the same
    // discretisation)
    m_fields[0]->GetCoords(x0, x1);

    for (int i = 0; i < nq; i++)
    {
        (m_fields[0]->UpdatePhys())[i] = amp * pow((1.0 / cosh(
                                sqrt(0.75 * (amp / (d * d * d))) *
                                (A * (x0[i] + x_offset) - C * time))), 2.0);
        (m_fields[1]->UpdatePhys())[i] = (amp / d) * pow((1.0 / cosh(
                                    sqrt(0.75 * (amp / (d * d * d))) *
                                    (A * (x0[i] + x_offset) - C * time)
                                )), 2.0) * sqrt(m_g * d);
    }

    Vmath::Sadd(nq, d, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
    Vmath::Vmul(nq,    m_fields[0]->GetPhys(), 1, m_fields[1]->GetPhys(),    1,
                                                  m_fields[1]->UpdatePhys(), 1);
    Vmath::Vcopy(nq, zeros, 1, m_fields[2]->UpdatePhys(), 1);
    Vmath::Vcopy(nq, zeros, 1, m_fields[3]->UpdatePhys(), 1);

    // Forward transform to fill the coefficient space
    for (int i = 0; i < 4; ++i)
    {
        m_fields[i]->SetPhysState(true);
        m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),
                              m_fields[i]->UpdateCoeffs());
    }

}

/**
 * @brief Set the initial conditions.
 */
void NonlinearPeregrine::v_SetInitialConditions(
        NekDouble initialtime,
        bool dumpInitialConditions,
        const int domain)
{
    boost::ignore_unused(domain);

    switch (m_problemType)
    {
        case eSolitaryWave:
        {
            LaitoneSolitaryWave(0.1, m_const_depth, 0.0, 0.0);
            break;
        }
        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime, false);
            break;
        }
    }

    if (dumpInitialConditions)
    {
        // Dump initial conditions to file
        Checkpoint_Output(0);
    }
}

} //end of namespace

