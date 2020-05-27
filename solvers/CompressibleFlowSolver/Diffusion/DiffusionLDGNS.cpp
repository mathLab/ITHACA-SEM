///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDGNS.cpp
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
// Description: LDGNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include "DiffusionLDGNS.h"
#include <iostream>
#include <iomanip>

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
std::string DiffusionLDGNS::type = GetDiffusionFactory().
RegisterCreatorFunction("LDGNS", DiffusionLDGNS::create);

DiffusionLDGNS::DiffusionLDGNS()
{
}

void DiffusionLDGNS::v_InitObject(
    LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    m_session = pSession;
    m_session->LoadParameter ("Twall", m_Twall, 300.15);

    // Setting up the normals
    std::size_t nDim = pFields[0]->GetCoordim(0);
    std::size_t nTracePts = pFields[0]->GetTrace()->GetTotPoints();

    m_spaceDim = nDim;
    if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
    {
        m_spaceDim = 3;
    }

    m_diffDim = m_spaceDim - nDim;

    m_traceVel = Array<OneD, Array<OneD, NekDouble> >{m_spaceDim};
    m_traceNormals = Array<OneD, Array<OneD, NekDouble> >{m_spaceDim};
    for (std::size_t i = 0; i < m_spaceDim; ++i)
    {
        m_traceVel[i] = Array<OneD, NekDouble> {nTracePts, 0.0};
        m_traceNormals[i] = Array<OneD, NekDouble> {nTracePts};
    }
    pFields[0]->GetTrace()->GetNormals(m_traceNormals);

    // Create equation of state object
    std::string eosType;
    m_session->LoadSolverInfo("EquationOfState",
                              eosType, "IdealGas");
    m_eos = GetEquationOfStateFactory()
                            .CreateInstance(eosType, m_session);


    // Set up {h} reference on the trace for penalty term
    //
    // Note, this shold be replaced with something smarter when merging
    // LDG with IP

    // Get min h per element
    std::size_t nElements = pFields[0]->GetExpSize();
    Array<OneD, NekDouble> hEle{nElements, 1.0};
    for (std::size_t e = 0; e < nElements; e++)
    {
        NekDouble h{1.0e+10};
        std::size_t expDim = pFields[0]->GetShapeDimension();
        switch(expDim)
        {
            case 3:
            {
                LocalRegions::Expansion3DSharedPtr exp3D;
                exp3D = pFields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                for (std::size_t i = 0; i < exp3D->GetNedges(); ++i)
                {
                    h = std::min(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }

            case 2:
            {
                LocalRegions::Expansion2DSharedPtr exp2D;
                exp2D = pFields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                for (std::size_t i = 0; i < exp2D->GetNedges(); ++i)
                {
                    h = std::min(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }
            case 1:
            {
                LocalRegions::Expansion1DSharedPtr exp1D;
                exp1D = pFields[0]->GetExp(e)->as<LocalRegions::Expansion1D>();

                h = std::min(h, exp1D->GetGeom1D()->GetVertex(0)->
                    dist(*(exp1D->GetGeom1D()->GetVertex(1))));

            break;
            }
            default:
            {
                ASSERTL0(false,"Dimension out of bound.")
            }
        }

        // Store scaling
        hEle[e] = h;
    }
    // Expand h from elements to points
    std::size_t nPts = pFields[0]->GetTotPoints();
    Array<OneD, NekDouble> hElePts{nPts, 0.0};
    Array<OneD, NekDouble> tmp;
    for (std::size_t e = 0; e < pFields[0]->GetExpSize(); e++)
    {
        std::size_t nElmtPoints     = pFields[0]->GetExp(e)->GetTotPoints();
        std::size_t physOffset      = pFields[0]->GetPhys_Offset(e);
        Vmath::Fill(nElmtPoints, hEle[e],
            tmp = hElePts + physOffset, 1);
    }
    // Get Fwd and Bwd traces
    Array<OneD, NekDouble> Fwd{nTracePts, 0.0};
    Array<OneD, NekDouble> Bwd{nTracePts, 0.0};
    pFields[0]->GetFwdBwdTracePhys(hElePts, Fwd, Bwd);
    // Fix boundaries
    std::size_t cnt = 0;
    std::size_t nBndRegions = pFields[0]->GetBndCondExpansions().size();
    // Loop on the boundary regions
    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        // Number of boundary regions related to region 'i'
        std::size_t nBndEdges = pFields[0]->
        GetBndCondExpansions()[i]->GetExpSize();

        if (pFields[0]->GetBndConditions()[i]->GetBoundaryConditionType()
            == SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Get value from interior
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts = pFields[0]->
            GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

            std::size_t id2 = pFields[0]->GetTrace()->
            GetPhys_Offset(pFields[0]->GetTraceMap()->
                           GetBndCondIDToGlobalTraceID(cnt++));

            Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &Bwd[id2], 1);
        }
    }
    // Get average of traces
    Array<OneD, NekDouble> traceH{nTracePts, 1.0};
    m_traceOneOverH = Array<OneD, NekDouble> {nTracePts, 1.0};
    Vmath::Svtsvtp(nTracePts, 0.5, Fwd, 1, 0.5, Bwd, 1, traceH, 1);
    // Multiply by coefficient = - C11 / h
    m_session->LoadParameter ("LDGNSc11", m_C11, 1.0);
    Vmath::Sdiv(nTracePts, -m_C11, traceH, 1, m_traceOneOverH, 1);
}

/**
 * @brief Calculate weak DG Diffusion in the LDG form for the
 * Navier-Stokes (NS) equations:
 *
 * \f$ \langle\psi, \hat{u}\cdot n\rangle
 *   - \langle\nabla\psi \cdot u\rangle
 *     \langle\phi, \hat{q}\cdot n\rangle -
 *     (\nabla \phi \cdot q) \rangle \f$
 *
 * The equations that need a diffusion operator are those related
 * with the velocities and with the energy.
 *
 */
 void DiffusionLDGNS::v_Diffuse(
    const std::size_t                                  nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    std::size_t nCoeffs   = fields[0]->GetNcoeffs();
    Array<OneD, Array<OneD, NekDouble> > tmp2{nConvectiveFields};
    for (std::size_t i = 0; i < nConvectiveFields; ++i)
    {
        tmp2[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};
    }
    v_DiffuseCoeffs(nConvectiveFields, fields, inarray, tmp2, pFwd, pBwd);
    for (std::size_t i = 0; i < nConvectiveFields; ++i)
    {
        fields[i]->BwdTrans             (tmp2[i], outarray[i]);
    }
}

void DiffusionLDGNS::v_DiffuseCoeffs(
    const std::size_t                                  nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    std::size_t nDim      = fields[0]->GetCoordim(0);
    std::size_t nPts      = fields[0]->GetTotPoints();
    std::size_t nCoeffs   = fields[0]->GetNcoeffs();
    std::size_t nScalars  = inarray.size();
    std::size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();

    Array<OneD, NekDouble>               tmp1{nCoeffs};
    Array<OneD, Array<OneD, NekDouble> > tmp2{nConvectiveFields};

    TensorOfArray3D<NekDouble> derivativesO1{m_spaceDim};

    for (std::size_t j = 0; j < m_spaceDim; ++j)
    {
        derivativesO1[j]   = Array<OneD, Array<OneD, NekDouble> >{nScalars};

        for (std::size_t i = 0; i < nScalars; ++i)
        {
            derivativesO1[j][i]   = Array<OneD, NekDouble>{nPts, 0.0};
        }
    }

    DiffuseCalculateDerivative(fields, inarray, derivativesO1, pFwd, pBwd);

    // Initialisation viscous tensor
    m_viscTensor = TensorOfArray3D<NekDouble> {m_spaceDim};
    Array<OneD, Array<OneD, NekDouble> > viscousFlux{nConvectiveFields};

    for (std::size_t j = 0; j < m_spaceDim; ++j)
    {
        m_viscTensor[j] = Array<OneD, Array<OneD, NekDouble> >{nScalars+1};
        for (std::size_t i = 0; i < nScalars+1; ++i)
        {
            m_viscTensor[j][i] = Array<OneD, NekDouble>{nPts, 0.0};
        }
    }

    for (std::size_t i = 0; i < nConvectiveFields; ++i)
    {
        viscousFlux[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    DiffuseVolumeFlux(fields, inarray, derivativesO1, m_viscTensor);

    // Compute u from q_{\eta} and q_{\xi}
    // Obtain numerical fluxes
    DiffuseTraceFlux(fields, inarray, derivativesO1, m_viscTensor, viscousFlux,
                        pFwd, pBwd);

    for (std::size_t i = 0; i < nConvectiveFields; ++i)
    {
        tmp2[i] = Array<OneD, NekDouble>{nCoeffs, 0.0};

        for (std::size_t j = 0; j < nDim; ++j)
        {
            fields[i]->IProductWRTDerivBase(j, m_viscTensor[j][i], tmp1);
            Vmath::Vadd(nCoeffs, tmp1, 1, tmp2[i], 1, tmp2[i], 1);
        }

        // Evaulate  <\phi, \hat{F}\cdot n> - outarray[i]
        Vmath::Neg                      (nCoeffs, tmp2[i], 1);
        fields[i]->AddTraceIntegral     (viscousFlux[i], tmp2[i]);
        fields[i]->SetPhysState         (false);
        fields[i]->MultiplyByElmtInvMass(tmp2[i], outarray[i]);
    }
}

void DiffusionLDGNS::v_DiffuseCalculateDerivative(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    TensorOfArray3D<NekDouble>                        &qfields,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    std::size_t nDim      = fields[0]->GetCoordim(0);
    std::size_t nCoeffs   = fields[0]->GetNcoeffs();
    std::size_t nScalars  = inarray.size();
    std::size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nConvectiveFields = fields.size();

    Array<OneD, NekDouble>               tmp1{nCoeffs};
    Array<OneD, Array<OneD, NekDouble> > tmp2{nConvectiveFields};
    TensorOfArray3D<NekDouble> numericalFluxO1{m_spaceDim};

    for (std::size_t j = 0; j < m_spaceDim; ++j)
    {
        numericalFluxO1[j] = Array<OneD, Array<OneD, NekDouble> >{nScalars};

        for (std::size_t i = 0; i < nScalars; ++i)
        {
            numericalFluxO1[j][i] = Array<OneD, NekDouble>{nTracePts, 0.0};
        }
    }

    NumericalFluxO1(fields, inarray, numericalFluxO1, pFwd, pBwd);

    for (std::size_t j = 0; j < nDim; ++j)
    {
        for (std::size_t i = 0; i < nScalars; ++i)
        {
            fields[i]->IProductWRTDerivBase (j, inarray[i], tmp1);
            Vmath::Neg                      (nCoeffs, tmp1, 1);
            fields[i]->AddTraceIntegral     (numericalFluxO1[j][i], tmp1);
            fields[i]->SetPhysState         (false);
            fields[i]->MultiplyByElmtInvMass(tmp1, tmp1);
            fields[i]->BwdTrans             (tmp1, qfields[j][i]);
        }
    }
    // For 3D Homogeneous 1D only take derivatives in 3rd direction
    if (m_diffDim == 1)
    {
        for (std::size_t i = 0; i < nScalars; ++i)
        {
            qfields[2][i] = m_homoDerivs[i];
        }
    }
}

void DiffusionLDGNS::v_DiffuseVolumeFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
    const Array<OneD, Array<OneD, NekDouble>>           &inarray,
    TensorOfArray3D<NekDouble>                          &qfields,
    TensorOfArray3D<NekDouble>                          &VolumeFlux,
    Array< OneD, int >                                  &nonZeroIndex)
{

    boost::ignore_unused(fields, nonZeroIndex);
    m_fluxVectorNS(inarray, qfields, VolumeFlux);

}

void DiffusionLDGNS::v_DiffuseTraceFlux(
    const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
    const Array<OneD, Array<OneD, NekDouble>>           &inarray,
    TensorOfArray3D<NekDouble>                          &qfields,
    TensorOfArray3D<NekDouble>                          &VolumeFlux,
    Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
    const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
    const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
    Array< OneD, int >                                  &nonZeroIndex)
{
    boost::ignore_unused(inarray, qfields, nonZeroIndex);
    NumericalFluxO2(fields, VolumeFlux, TraceFlux, pFwd, pBwd);
}

/**
 * @brief Builds the numerical flux for the 1st order derivatives
 *
 */
void DiffusionLDGNS::NumericalFluxO1(
    const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
    const Array<OneD, Array<OneD, NekDouble> >               &inarray,
    TensorOfArray3D<NekDouble>                               &numericalFluxO1,
    const Array<OneD, Array<OneD, NekDouble> >               &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >               &pBwd)
{
    std::size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nScalars  = inarray.size();

    //Upwind
    Array<OneD, Array<OneD, NekDouble> > numflux{nScalars};
    for (std::size_t i = 0; i < nScalars; ++i)
    {
        numflux[i] = {pFwd[i]};
    }

    // Modify the values in case of boundary interfaces
    if (fields[0]->GetBndCondExpansions().size())
    {
        ApplyBCsO1(fields, inarray, pFwd, pBwd, numflux);
    }

    // Splitting the numerical flux into the dimensions
    for (std::size_t j = 0; j < m_spaceDim; ++j)
    {
        for (std::size_t i = 0; i < nScalars; ++i)
        {
            Vmath::Vmul(nTracePts, m_traceNormals[j], 1,numflux[i], 1,
                numericalFluxO1[j][i], 1);
        }
    }
}

/**
 * @brief Imposes appropriate bcs for the 1st order derivatives
 *
 */
void DiffusionLDGNS::ApplyBCsO1(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
          Array<OneD, Array<OneD, NekDouble> >        &fluxO1)
{
    boost::ignore_unused(pBwd);

    std::size_t nTracePts = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nScalars  = inarray.size();

    Array<OneD, NekDouble> tmp1{nTracePts, 0.0};
    Array<OneD, NekDouble> tmp2{nTracePts, 0.0};
    Array<OneD, NekDouble> Tw{nTracePts, m_Twall};

    Array< OneD, Array<OneD, NekDouble > > scalarVariables{nScalars};

    // Extract internal values of the scalar variables for Neumann bcs
    for (std::size_t i = 0; i < nScalars; ++i)
    {
        scalarVariables[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    // Compute boundary conditions for velocities
    for (std::size_t i = 0; i < nScalars-1; ++i)
    {
        // Note that cnt has to loop on nBndRegions and nBndEdges
        // and has to be reset to zero per each equation
        std::size_t cnt = 0;
        std::size_t nBndRegions = fields[i+1]->
            GetBndCondExpansions().size();
        for (std::size_t j = 0; j < nBndRegions; ++j)
        {
            if (fields[i+1]->GetBndConditions()[j]->
                GetBoundaryConditionType() == SpatialDomains::ePeriodic)
            {
                continue;
            }

            std::size_t nBndEdges = fields[i+1]->
                GetBndCondExpansions()[j]->GetExpSize();
            for (std::size_t e = 0; e < nBndEdges; ++e)
            {
                std::size_t nBndEdgePts = fields[i+1]->
                    GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                std::size_t id1 = fields[i+1]->
                    GetBndCondExpansions()[j]->GetPhys_Offset(e);

                std::size_t id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondIDToGlobalTraceID(cnt++));

                if (boost::iequals(fields[i]->GetBndConditions()[j]->
                    GetUserDefined(),"WallViscous") ||
                    boost::iequals(fields[i]->GetBndConditions()[j]->
                    GetUserDefined(),"WallAdiabatic"))
                {
                    // Reinforcing bcs for velocity in case of Wall bcs
                    Vmath::Zero(nBndEdgePts, &scalarVariables[i][id2], 1);

                }
                else if (
                    boost::iequals(fields[i]->GetBndConditions()[j]->
                    GetUserDefined(),"Wall") ||
                    boost::iequals(fields[i]->GetBndConditions()[j]->
                    GetUserDefined(),"Symmetry"))
                {
                    // Symmetry bc: normal velocity is zero
                    //    get all velocities at once because we need u.n
                    if (i==0)
                    {
                        //    tmp1 = -(u.n)
                        Vmath::Zero(nBndEdgePts, tmp1, 1);
                        for (std::size_t k = 0; k < nScalars-1; ++k)
                        {
                            Vmath::Vdiv(nBndEdgePts,
                              &(fields[k+1]->GetBndCondExpansions()[j]->
                                  UpdatePhys())[id1], 1,
                              &(fields[0]->GetBndCondExpansions()[j]->
                                  UpdatePhys())[id1], 1,
                              &scalarVariables[k][id2], 1);
                            Vmath::Vvtvp(nBndEdgePts,
                                    &m_traceNormals[k][id2], 1,
                                    &scalarVariables[k][id2], 1,
                                    &tmp1[0], 1,
                                    &tmp1[0], 1);
                        }
                        Vmath::Smul(nBndEdgePts, -1.0,
                                    &tmp1[0], 1,
                                    &tmp1[0], 1);

                        //    u_i - (u.n)n_i
                        for (std::size_t k = 0; k < nScalars-1; ++k)
                        {
                            Vmath::Vvtvp(nBndEdgePts,
                                    &tmp1[0], 1,
                                    &m_traceNormals[k][id2], 1,
                                    &scalarVariables[k][id2], 1,
                                    &scalarVariables[k][id2], 1);
                        }
                    }
                }
                else if (fields[i]->GetBndConditions()[j]->
                         GetBoundaryConditionType() ==
                         SpatialDomains::eDirichlet)
                {
                    // Imposing velocity bcs if not Wall
                    Vmath::Vdiv(nBndEdgePts,
                            &(fields[i+1]->GetBndCondExpansions()[j]->
                              UpdatePhys())[id1], 1,
                            &(fields[0]->GetBndCondExpansions()[j]->
                              UpdatePhys())[id1], 1,
                            &scalarVariables[i][id2], 1);
                }

                // For Dirichlet boundary condition: uflux = u_bcs
                if (fields[i]->GetBndConditions()[j]->
                    GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet)
                {
                    Vmath::Vcopy(nBndEdgePts,
                                 &scalarVariables[i][id2], 1,
                                 &fluxO1[i][id2], 1);
                }

                // For Neumann boundary condition: uflux = u_+
                else if ((fields[i]->GetBndConditions()[j])->
                         GetBoundaryConditionType() ==
                         SpatialDomains::eNeumann)
                {
                    Vmath::Vcopy(nBndEdgePts,
                                 &pFwd[i][id2], 1,
                                 &fluxO1[i][id2], 1);
                }

                // Building kinetic energy to be used for T bcs
                Vmath::Vmul(nBndEdgePts,
                            &scalarVariables[i][id2], 1,
                            &scalarVariables[i][id2], 1,
                            &tmp1[id2], 1);

                Vmath::Smul(nBndEdgePts, 0.5,
                            &tmp1[id2], 1,
                            &tmp1[id2], 1);

                Vmath::Vadd(nBndEdgePts,
                            &tmp2[id2], 1,
                            &tmp1[id2], 1,
                            &tmp2[id2], 1);
            }
        }
    }

    // Compute boundary conditions  for temperature
    std::size_t cnt = 0;
    std::size_t nBndRegions = fields[nScalars]->
    GetBndCondExpansions().size();
    for (std::size_t j = 0; j < nBndRegions; ++j)
    {
        if (fields[nScalars]->GetBndConditions()[j]->
            GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        std::size_t nBndEdges = fields[nScalars]->
        GetBndCondExpansions()[j]->GetExpSize();
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts = fields[nScalars]->
            GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

            std::size_t id1 = fields[nScalars]->
            GetBndCondExpansions()[j]->GetPhys_Offset(e);

            std::size_t id2 = fields[0]->GetTrace()->
            GetPhys_Offset(fields[0]->GetTraceMap()->
                           GetBndCondIDToGlobalTraceID(cnt++));

            // Imposing Temperature Twall at the wall
            if (boost::iequals(fields[nScalars]->GetBndConditions()[j]->
                GetUserDefined(),"WallViscous"))
            {
                Vmath::Vcopy(nBndEdgePts,
                             &Tw[0], 1,
                             &scalarVariables[nScalars-1][id2], 1);
            }
            // Imposing Temperature through condition on the Energy
            // for no wall boundaries (e.g. farfield)
            else if (fields[nScalars]->GetBndConditions()[j]->
                     GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                // Use equation of state to evaluate temperature
                NekDouble rho, ene;
                for (std::size_t n = 0; n < nBndEdgePts; ++n)
                {
                    rho = fields[0]->
                              GetBndCondExpansions()[j]->
                              GetPhys()[id1+n];
                    ene = fields[nScalars]->
                              GetBndCondExpansions()[j]->
                              GetPhys()[id1 +n] / rho - tmp2[id2+n];
                    scalarVariables[nScalars-1][id2+n] =
                            m_eos->GetTemperature(rho, ene);
                }
            }

            // For Dirichlet boundary condition: uflux = u_bcs
            if (fields[nScalars]->GetBndConditions()[j]->
                GetBoundaryConditionType() == SpatialDomains::eDirichlet &&
                !boost::iequals(fields[nScalars]->GetBndConditions()[j]
                ->GetUserDefined(), "WallAdiabatic"))
            {
                Vmath::Vcopy(nBndEdgePts,
                             &scalarVariables[nScalars-1][id2], 1,
                             &fluxO1[nScalars-1][id2], 1);

            }

            // For Neumann boundary condition: uflux = u_+
            else if (((fields[nScalars]->GetBndConditions()[j])->
                      GetBoundaryConditionType() == SpatialDomains::eNeumann) ||
                      boost::iequals(fields[nScalars]->GetBndConditions()[j]->
                                    GetUserDefined(), "WallAdiabatic"))
            {
                Vmath::Vcopy(nBndEdgePts,
                             &pFwd[nScalars-1][id2], 1,
                             &fluxO1[nScalars-1][id2], 1);

            }
        }
    }
}

/**
 * @brief Build the numerical flux for the 2nd order derivatives
 *
 */
void DiffusionLDGNS::NumericalFluxO2(
    const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
    TensorOfArray3D<NekDouble>                               &qfield,
    Array<OneD, Array<OneD, NekDouble> >                     &qflux,
    const Array<OneD, Array<OneD, NekDouble> >               &uFwd,
    const Array<OneD, Array<OneD, NekDouble> >               &uBwd)
{
    std::size_t nTracePts  = fields[0]->GetTrace()->GetTotPoints();
    std::size_t nVariables = fields.size();

    // Initialize penalty flux
    Array<OneD, Array<OneD, NekDouble> >  fluxPen{nVariables-1};
    for (std::size_t i = 0; i < nVariables-1; ++i)
    {
        fluxPen[i] = Array<OneD, NekDouble>{nTracePts, 0.0};
    }

    // Get penalty flux
    m_fluxPenaltyNS(uFwd, uBwd, fluxPen);

    // Evaluate Riemann flux
    // qflux = \hat{q} \cdot u = q \cdot n
    // Notice: i = 1 (first row of the viscous tensor is zero)

    Array<OneD, NekDouble > qFwd{nTracePts};
    Array<OneD, NekDouble > qBwd{nTracePts};
    Array<OneD, NekDouble > qtemp{nTracePts};
    Array<OneD, NekDouble > qfluxtemp{nTracePts, 0.0};
    std::size_t nDim = fields[0]->GetCoordim(0);
    for (std::size_t i = 1; i < nVariables; ++i)
    {
        qflux[i] = Array<OneD, NekDouble> {nTracePts, 0.0};
        for (std::size_t j = 0; j < nDim; ++j)
        {
            // Compute qFwd and qBwd value of qfield in position 'ji'
            fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);

            // Downwind
            Vmath::Vcopy(nTracePts, qBwd, 1, qfluxtemp, 1);

            // Multiply the Riemann flux by the trace normals
            Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                        qfluxtemp, 1);

            // Add penalty term
            Vmath::Vvtvp(nTracePts, m_traceOneOverH, 1, fluxPen[i-1], 1,
                qfluxtemp, 1, qfluxtemp, 1);

            // Impose weak boundary condition with flux
            if (fields[0]->GetBndCondExpansions().size())
            {
                ApplyBCsO2(fields, i, j, qfield[j][i], qFwd, qBwd, qfluxtemp);
            }

            // Store the final flux into qflux
            Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1, qflux[i], 1);
        }
    }
}


/**
 * @brief Imposes appropriate bcs for the 2nd order derivatives
 *
 */
void DiffusionLDGNS::ApplyBCsO2(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const std::size_t                                  var,
    const std::size_t                                  dir,
    const Array<OneD, const NekDouble>                &qfield,
    const Array<OneD, const NekDouble>                &qFwd,
    const Array<OneD, const NekDouble>                &qBwd,
          Array<OneD,       NekDouble>                &penaltyflux)
{
    boost::ignore_unused(qfield, qBwd);

    std::size_t cnt = 0;
    std::size_t nBndRegions = fields[var]->GetBndCondExpansions().size();
    // Loop on the boundary regions to apply appropriate bcs
    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        // Number of boundary regions related to region 'i'
        std::size_t nBndEdges = fields[var]->
        GetBndCondExpansions()[i]->GetExpSize();

        if (fields[var]->GetBndConditions()[i]->GetBoundaryConditionType()
            == SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Weakly impose bcs by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts = fields[var]->
            GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

            std::size_t id2 = fields[0]->GetTrace()->
            GetPhys_Offset(fields[0]->GetTraceMap()->
                           GetBndCondIDToGlobalTraceID(cnt++));

            // In case of Dirichlet bcs:
            // uflux = gD
            if (fields[var]->GetBndConditions()[i]->
               GetBoundaryConditionType() == SpatialDomains::eDirichlet
               && !boost::iequals(fields[var]->GetBndConditions()[i]->
                                  GetUserDefined(), "WallAdiabatic"))
            {
                Vmath::Vmul(nBndEdgePts,
                            &m_traceNormals[dir][id2], 1,
                            &qFwd[id2], 1,
                            &penaltyflux[id2], 1);
            }
            // 3.4) In case of Neumann bcs:
            // uflux = u+
            else if ((fields[var]->GetBndConditions()[i])->
                GetBoundaryConditionType() == SpatialDomains::eNeumann)
            {
                ASSERTL0(false,
                         "Neumann bcs not implemented for LDGNS");

                /*
                Vmath::Vmul(nBndEdgePts,
                            &m_traceNormals[dir][id2], 1,
                            &(fields[var]->
                              GetBndCondExpansions()[i]->
                              UpdatePhys())[id1], 1,
                            &penaltyflux[id2], 1);
                 */
            }
            else if (boost::iequals(fields[var]->GetBndConditions()[i]->
                                   GetUserDefined(), "WallAdiabatic"))
            {
                if ((var == m_spaceDim + 1))
                {
                    Vmath::Zero(nBndEdgePts, &penaltyflux[id2], 1);
                }
                else
                {

                    Vmath::Vmul(nBndEdgePts,
                                &m_traceNormals[dir][id2], 1,
                                &qFwd[id2], 1,
                                &penaltyflux[id2], 1);

                }
            }
        }
    }
}
}//end of namespace Nektar
