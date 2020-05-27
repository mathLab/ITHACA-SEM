///////////////////////////////////////////////////////////////////////////////
//
// File: ImageWarpingSystem.cpp
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
// Description: Image warping solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <ImageWarpingSolver/EquationSystems/ImageWarpingSystem.h>
#include <MultiRegions/ContField2D.h>

using namespace std;

namespace Nektar
{
    string ImageWarpingSystem::className =
        GetEquationSystemFactory().RegisterCreatorFunction(
            "ImageWarpingSystem",
            ImageWarpingSystem::create,
            "Image warping system.");

    ImageWarpingSystem::ImageWarpingSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
    }

    void ImageWarpingSystem::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

        // Define Velocity fields
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        int nq = m_fields[0]->GetNpoints();

        for(int i = 0; i < m_spacedim; ++i)
        {
            m_velocity[i] = Array<OneD, NekDouble> (nq,0.0);
        }

        // Bit of a hack: redefine u/v fields so they are continuous for
        // Helmholtz solve.
        MultiRegions::ContField2DSharedPtr fld =
            MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(2));
        m_fields[2] = fld;
        m_fields[3] = MemoryManager<MultiRegions::ContField2D>
            ::AllocateSharedPtr(*fld,m_graph,m_session->GetVariable(3));

        // Tell UnsteadySystem to only integrate first two fields (i.e. I and
        // phi).
        m_intVariables.push_back(0);
        m_intVariables.push_back(1);

        // Define the normal velocity fields
        if (m_fields[0]->GetTrace())
        {
            m_traceVn = Array<OneD, NekDouble>(GetTraceNpoints());
        }

        string advName;
        string riemName;
        m_session->LoadSolverInfo(
            "AdvectionType", advName, "WeakDG");
        m_advObject = SolverUtils::
            GetAdvectionFactory().CreateInstance(advName, advName);
        m_advObject->SetFluxVector(
            &ImageWarpingSystem::GetFluxVector, this);
        m_session->LoadSolverInfo(
            "UpwindType", riemName, "Upwind");
        m_riemannSolver = SolverUtils::
            GetRiemannSolverFactory().CreateInstance(riemName, m_session);
        m_riemannSolver->SetScalar(
            "Vn", &ImageWarpingSystem::GetNormalVelocity, this);

        m_advObject->SetRiemannSolver(m_riemannSolver);
        m_advObject->InitObject(m_session, m_fields);

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&ImageWarpingSystem::DoOdeRhs,        this);
            m_ode.DefineProjection (&ImageWarpingSystem::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    ImageWarpingSystem::~ImageWarpingSystem()
    {

    }

    void ImageWarpingSystem::DoOdeRhs(
        const Array<OneD, const Array<OneD,NekDouble> > &inarray,
              Array<OneD,       Array<OneD,NekDouble> > &outarray,
        const NekDouble time)
    {
        boost::ignore_unused(time);

        int npoints = GetNpoints();
        int ncoeffs = inarray[0].size();
        StdRegions::ConstFactorMap factors;

        // Load parameter alpha.
        m_session->LoadParameter("Alpha", m_alpha);

        ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
                 "CG not implemented yet.");

        // Set up storage arrays.
        Array<OneD, NekDouble> tmp(npoints);
        Array<OneD, NekDouble> alloc(3*npoints);
        Array<OneD, NekDouble> dIdx1(alloc);
        Array<OneD, NekDouble> dIdx2(alloc+npoints);
        Array<OneD, NekDouble> dIdx3(alloc+2*npoints);

        // Calculate grad I.
        m_fields[0]->PhysDeriv(inarray[0], dIdx1, dIdx2);

        // Set factors.
        // TODO: Check - should be -1?
        factors[StdRegions::eFactorLambda] = 1.0 / m_alpha / m_alpha;

        // Multiply by phi, and perform Helmholtz solve to calculate the
        // advection velocity field.
        for (int i = 0; i < 2; ++i)
        {
            Vmath::Vmul(npoints, &alloc[i*npoints], 1, inarray[1].get(), 1,
                        m_fields[i+2]->UpdatePhys().get(), 1);
            Vmath::Smul(npoints, 1/m_alpha/m_alpha, m_fields[i+2]->GetPhys().get(), 1,
                        m_fields[i+2]->UpdatePhys().get(), 1);
            m_fields[i+2]->HelmSolve(m_fields[i+2]->GetPhys(),
                                     m_fields[i+2]->UpdateCoeffs(), factors);
            m_fields[i+2]->BwdTrans(m_fields[i+2]->GetCoeffs(),
                                    m_velocity[i]);
        }

        // Calculate the weak advection operator for I and phi - result is put
        // in WeakAdv and is in physical space.
        m_advObject->Advect(2, m_fields, m_velocity, inarray,
                            outarray, 0.0);
        for(int i = 0; i < 2; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }

        // Calculate du/dx -> dIdx1, dv/dy -> dIdx2.
        m_fields[2]->PhysDeriv(m_velocity[0], dIdx1, dIdx3);
        m_fields[3]->PhysDeriv(m_velocity[1], dIdx3, dIdx2);

        // Calculate RHS = I*div(u) = I*du/dx + I*dv/dy -> dIdx1.
        Vmath::Vvtvvtp(npoints, dIdx1.get(), 1, inarray[0].get(), 1,
                       dIdx2.get(), 1, inarray[0].get(), 1,
                       dIdx1.get(), 1);

        // Take inner product to get to coefficient space.
        Array<OneD, NekDouble> tmp2(ncoeffs);
        m_fields[0]->IProductWRTBase(dIdx1, tmp2);

        // Multiply by elemental inverse mass matrix, backwards transform
        m_fields[0]->MultiplyByElmtInvMass(tmp2, tmp2);
        m_fields[0]->BwdTrans(tmp2,tmp);
        Vmath::Vadd(npoints, outarray[0], 1, tmp, 1, outarray[0], 1);
    }



    /**
     *
     */
    void ImageWarpingSystem::DoOdeProjection(const Array<OneD,
                                            const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                            const NekDouble time)
    {
        int nvariables = inarray.size();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(int i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
                }
            }
            break;
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

    /**
     * @brief Get the normal velocity
     */
    Array<OneD, NekDouble> &ImageWarpingSystem::GetNormalVelocity()
    {
        // Number of trace (interface) points
        int nTracePts = GetTraceNpoints();

        // Auxiliary variable to compute the normal velocity
        Array<OneD, NekDouble> tmp(nTracePts);

        // Reset the normal velocity
        Vmath::Zero(nTracePts, m_traceVn, 1);

        for (int i = 0; i < m_velocity.size(); ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], tmp);

            Vmath::Vvtvp(nTracePts,
                         m_traceNormals[i], 1,
                         tmp,               1,
                         m_traceVn,         1,
                         m_traceVn,         1);
        }

        return m_traceVn;
    }

    void ImageWarpingSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        ASSERTL1(flux[0].size() == m_velocity.size(),
                 "Dimension of flux array and velocity array do not match");

        int nq = physfield[0].size();

        for (int i = 0; i < flux.size(); ++i)
        {
            for (int j = 0; j < flux[0].size(); ++j)
            {
                Vmath::Vmul(nq, physfield[i], 1, m_velocity[j], 1,
                            flux[i][j], 1);
            }
        }
    }

    void ImageWarpingSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        AdvectionSystem::v_GenerateSummary(s);
    }
}
