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
// Description: LDGNS diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include "DiffusionLDGNS.h"
#include <iostream>
#include <iomanip>

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
        int i;
        int nDim = pFields[0]->GetCoordim(0);
        int nTracePts = pFields[0]->GetTrace()->GetTotPoints();

        m_spaceDim = nDim;
        if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
        {
            m_spaceDim = 3;
        }

        m_diffDim = m_spaceDim - nDim;

        m_traceVel = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
        for(i = 0; i < m_spaceDim; ++i)
        {
            m_traceVel[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
            m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
        }
        pFields[0]->GetTrace()->GetNormals(m_traceNormals);

        // Create equation of state object
        std::string eosType;
        m_session->LoadSolverInfo("EquationOfState",
                                  eosType, "IdealGas");
        m_eos = GetEquationOfStateFactory()
                                .CreateInstance(eosType, m_session);
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
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
    {
        int nCoeffs   = fields[0]->GetNcoeffs();
        Array<OneD, Array<OneD, NekDouble> > tmp2(nConvectiveFields);
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
        }
        v_Diffuse_coeff(nConvectiveFields,fields,inarray,tmp2,pFwd,pBwd);
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            fields[i]->BwdTrans             (tmp2[i], outarray[i]);
        }
    }

    // TODO:: REPLACED
    void DiffusionLDGNS::v_Diffuse_coeff(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
              Array<OneD, Array<OneD, NekDouble> >        &outarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
    {
        int i, j;
        int nDim      = fields[0]->GetCoordim(0);
        int nScalars  = inarray.num_elements();
        int nPts      = fields[0]->GetTotPoints();
        int nCoeffs   = fields[0]->GetNcoeffs();
        int nTracePts = fields[0]->GetTrace()->GetTotPoints();

        Array<OneD, NekDouble>               tmp1(nCoeffs);
        Array<OneD, Array<OneD, NekDouble> > tmp2(nConvectiveFields);

        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                derivativesO1(m_spaceDim);

        for (j = 0; j < m_spaceDim; ++j)
        {
            derivativesO1[j]   = Array<OneD, Array<OneD, NekDouble> >(
                                                                nScalars);

            for (i = 0; i < nScalars; ++i)
            {
                derivativesO1[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
            }
        }

        DiffuseCalculateDerivative(nConvectiveFields,fields,inarray,derivativesO1,pFwd,pBwd);

        // Initialisation viscous tensor
        m_viscTensor = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                               (m_spaceDim);
        Array<OneD, Array<OneD, NekDouble> > viscousFlux(nConvectiveFields);

        for (j = 0; j < m_spaceDim; ++j)
        {
            m_viscTensor[j] = Array<OneD, Array<OneD, NekDouble> >(
                                                                nScalars+1);
            for (i = 0; i < nScalars+1; ++i)
            {
                m_viscTensor[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
            }
        }

        for (i = 0; i < nConvectiveFields; ++i)
        {
            viscousFlux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        DiffuseVolumeFlux(nConvectiveFields,fields,inarray,derivativesO1,m_viscTensor);

        // Compute u from q_{\eta} and q_{\xi}
        // Obtain numerical fluxes
        DiffuseTraceFlux(nConvectiveFields,fields,inarray,derivativesO1,m_viscTensor,viscousFlux,pFwd,pBwd);
        // v_NumericalFluxO2(fields, inarray, m_viscTensor, viscousFlux);

        for (i = 0; i < nConvectiveFields; ++i)
        {
            tmp2[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);

            for (j = 0; j < nDim; ++j)
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

     /**
     * @brief Calculate  LDG First Order derivatives
     * @param inarray 2D [u,v,T]
     * @param inarrayderivative 2D [du_dx,dv_dx,dT_dx],[du_dy,dv_dy,dT_dy]
     * The equations that need a diffusion operator are those related
     * with the velocities and with the energy.
     *
     */
    void DiffusionLDGNS::v_DiffuseCalculateDerivative(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
    {
        int nDim      = fields[0]->GetCoordim(0);
        int nScalars  = inarray.num_elements();
        int nPts      = fields[0]->GetTotPoints();
        int nCoeffs   = fields[0]->GetNcoeffs();
        int nTracePts = fields[0]->GetTrace()->GetTotPoints();

        Array<OneD, NekDouble>               tmp1(nCoeffs);
        Array<OneD, Array<OneD, NekDouble> > tmp2(nConvectiveFields);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > numericalFluxO1(m_spaceDim);

        for (int j = 0; j < m_spaceDim; ++j)
        {
            numericalFluxO1[j] = Array<OneD, Array<OneD, NekDouble> >(nScalars);

            for (int i = 0; i < nScalars; ++i)
            {
                numericalFluxO1[j][i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
        }

        // Compute the numerical fluxes for the first order derivatives
        v_NumericalFluxO1(fields, inarray, numericalFluxO1, pFwd, pBwd);

        for (int j = 0; j < nDim; ++j)
        {
            for (int i = 0; i < nScalars; ++i)
            {
                fields[i]->IProductWRTDerivBase (j, inarray[i], tmp1);
                Vmath::Neg                      (nCoeffs, tmp1, 1);
                fields[i]->AddTraceIntegral     (numericalFluxO1[j][i], tmp1);
                fields[i]->SetPhysState         (false);
                fields[i]->MultiplyByElmtInvMass(tmp1, tmp1);
                fields[i]->BwdTrans             (tmp1, inarrayderivative[j][i]);
            }
        }

        // For 3D Homogeneous 1D only take derivatives in 3rd direction
        if (m_diffDim == 1)
        {
            for (int i = 0; i < nScalars; ++i)
            {
                inarrayderivative[2][i] = m_homoDerivs[i];
            }
        }
    }

     /**
     * @brief Calculate LDG Diffusion Volume Flux
     * @param inarray 2D [u,v,T]
     * @param inarrayderivative 2D [du_dx,dv_dx,dT_dx],[du_dy,dv_dy,dT_dy]
     * @param VolumeFlux volumeflxu integrant
     * The equations that need a diffusion operator are those related
     * with the velocities and with the energy.
     *
     */
    void DiffusionLDGNS::v_DiffuseVolumeFlux(
        const int                                           nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        const Array<OneD, Array<OneD, NekDouble>>           &inarray,
        Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
        Array< OneD, int >                                  &nonZeroIndex)
    {
        //?Homogeneous? Need to change here to inarrayderivative
        //For 3D Homogeneous 1D only take derivatives in 3rd direction
        // if (m_diffDim == 1)
        // {
        //     for (i = 0; i < nScalars; ++i)
        //     {
        //         derivativesO1[2][i] = m_homoDerivs[i];
        //     }
        // }

        int nScalars  = inarray.num_elements();
        int nPts      = fields[0]->GetTotPoints();

        // Initialisation viscous tensor
        // ASSERTL0(VolumeFlux.num_elements()==m_spaceDim,'Array dimension of VolumeFlux is wrong');
        // ASSERTL0(VolumeFlux[0].num_elements()==(nScalars+1),'Primitive variables used should be one smaller than nConvectiveVariables');

        m_fluxVectorNS(inarray, inarrayderivative, VolumeFlux);

    }


     /**
     * @brief Calculate LDG Diffusion Trace Flux
     * @param inarray 2D [u,v,T]
     * @param inarrayderivative 2D [du_dx,dv_dx,dT_dx],[du_dy,dv_dy,dT_dy]
     * @param VolumeFlux volumeflxu integrant
     * The equations that need a diffusion operator are those related
     * with the velocities and with the energy.
     *
     */
    void DiffusionLDGNS::v_DiffuseTraceFlux(
        const int                                           nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
        const Array<OneD, Array<OneD, NekDouble>>           &inarray,
        Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
        Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
        const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
        const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
        Array< OneD, int >                                  &nonZeroIndex)
    {
        ////?Homogeneous? Need to change here to inarrayderivative
        // For 3D Homogeneous 1D only take derivatives in 3rd direction
        // if (m_diffDim == 1)
        // {
        //     for (i = 0; i < nScalars; ++i)
        //     {
        //         derivativesO1[2][i] = m_homoDerivs[i];
        //     }
        // }

        // Compute u from q_{\eta} and q_{\xi}
        // Obtain numerical fluxes
        v_NumericalFluxO2(fields, inarray, VolumeFlux, TraceFlux);
    }

    /**
     * @brief Builds the numerical flux for the 1st order derivatives
     *
     */
    void DiffusionLDGNS::v_NumericalFluxO1(
        const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
        const Array<OneD, Array<OneD, NekDouble> >               &inarray,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                        &numericalFluxO1,
        const Array<OneD, Array<OneD, NekDouble> >               &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >               &pBwd)
    {
        int i, j;
        int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
        int nScalars   = inarray.num_elements();
        int nDim       = fields[0]->GetCoordim(0);

        Array<OneD, NekDouble > Vn      (nTracePts, 0.0);

        // Get the normal velocity Vn
        for(i = 0; i < nDim; ++i)
        {
            Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1,
                         Vn, 1, Vn, 1);
        }

        // Store forwards/backwards space along trace space
        Array<OneD, NekDouble> Fwd;
        Array<OneD, NekDouble> Bwd;
        Array<OneD, Array<OneD, NekDouble> > numflux(nScalars);

        for (i = 0; i < nScalars; ++i)
        {
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                Fwd    = Array<OneD, NekDouble>(nTracePts);
                Bwd    = Array<OneD, NekDouble>(nTracePts);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd, Bwd);
            }
            else
            {
                Fwd    = pFwd[i];
                Bwd    = pBwd[i];
            }
            numflux[i] = Array<OneD, NekDouble>(nTracePts);
            fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux[i]);
        }

        // Extract internal values of the scalar variables for Neumann bcs
        Array< OneD, Array<OneD, NekDouble > > uplus(nScalars);

        for (i = 0; i < nScalars; ++i)
        {
            uplus[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            fields[i]->ExtractTracePhys(inarray[i], uplus[i]);
        }

        // Modify the values in case of boundary interfaces
        if (fields[0]->GetBndCondExpansions().num_elements())
        {
            v_WeakPenaltyO1(fields, inarray, uplus, numflux);
        }

        // Splitting the numerical flux into the dimensions
        for (j = 0; j < m_spaceDim; ++j)
        {
            for (i = 0; i < nScalars; ++i)
            {
                Vmath::Vmul(nTracePts, m_traceNormals[j], 1,
                            numflux[i], 1, numericalFluxO1[j][i], 1);
            }
        }
    }

    /**
     * @brief Imposes appropriate bcs for the 1st order derivatives
     *
     */
    void DiffusionLDGNS::v_WeakPenaltyO1(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &uplus,
              Array<OneD, Array<OneD, NekDouble> >        &penaltyfluxO1)
    {
        int cnt;
        int i, j, e;
        int id1, id2;

        int nBndEdgePts, nBndEdges, nBndRegions;

        int nTracePts = fields[0]->GetTrace()->GetTotPoints();
        int nScalars  = inarray.num_elements();

        Array<OneD, NekDouble> tmp1(nTracePts, 0.0);
        Array<OneD, NekDouble> tmp2(nTracePts, 0.0);
        Array<OneD, NekDouble> Tw(nTracePts, m_Twall);

        Array< OneD, Array<OneD, NekDouble > > scalarVariables(nScalars);

        // Extract internal values of the scalar variables for Neumann bcs
        for (i = 0; i < nScalars; ++i)
        {
            scalarVariables[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        // Compute boundary conditions for velocities
        for (i = 0; i < nScalars-1; ++i)
        {
            // Note that cnt has to loop on nBndRegions and nBndEdges
            // and has to be reset to zero per each equation
            cnt = 0;
            nBndRegions = fields[i+1]->
            GetBndCondExpansions().num_elements();
            for (j = 0; j < nBndRegions; ++j)
            {
                if (fields[i+1]->GetBndConditions()[j]->
                    GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {
                    continue;
                }

                nBndEdges = fields[i+1]->
                GetBndCondExpansions()[j]->GetExpSize();
                for (e = 0; e < nBndEdges; ++e)
                {
                    nBndEdgePts = fields[i+1]->
                    GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                    id1 = fields[i+1]->
                    GetBndCondExpansions()[j]->GetPhys_Offset(e);

                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));

                    if (boost::iequals(fields[i]->GetBndConditions()[j]->
                        GetUserDefined(),"WallViscous") ||
                        boost::iequals(fields[i]->GetBndConditions()[j]->
                        GetUserDefined(),"WallAdiabatic"))
                    {
                        // Reinforcing bcs for velocity in case of Wall bcs
                        Vmath::Zero(nBndEdgePts,
                                    &scalarVariables[i][id2], 1);

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
                            for (int k = 0; k < nScalars-1; ++k)
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
                            for (int k = 0; k < nScalars-1; ++k)
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
                                     &penaltyfluxO1[i][id2], 1);
                    }

                    // For Neumann boundary condition: uflux = u_+
                    else if ((fields[i]->GetBndConditions()[j])->
                             GetBoundaryConditionType() ==
                             SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts,
                                     &uplus[i][id2], 1,
                                     &penaltyfluxO1[i][id2], 1);
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
        cnt = 0;
        nBndRegions = fields[nScalars]->
        GetBndCondExpansions().num_elements();
        for (j = 0; j < nBndRegions; ++j)
        {
            if (fields[nScalars]->GetBndConditions()[j]->
                GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                continue;
            }

            nBndEdges = fields[nScalars]->
            GetBndCondExpansions()[j]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = fields[nScalars]->
                GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                id1 = fields[nScalars]->
                GetBndCondExpansions()[j]->GetPhys_Offset(e);

                id2 = fields[0]->GetTrace()->
                GetPhys_Offset(fields[0]->GetTraceMap()->
                               GetBndCondTraceToGlobalTraceMap(cnt++));

                // Imposing Temperature Twall at the wall
                // The original code is "fields[i]->", in which i depends on results of previous loop
                // never do that!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (boost::iequals(fields[nScalars]->GetBndConditions()[j]->
                    GetUserDefined(),"WallViscous"))
                {
                    Vmath::Vcopy(nBndEdgePts,
                                 &Tw[0], 1,
                                 &scalarVariables[nScalars-1][id2], 1);
                }
                // Imposing Temperature through condition on the Energy
                // for no wall boundaries (e.g. farfield)
                else if (fields[i]->GetBndConditions()[j]->
                         GetBoundaryConditionType() ==
                         SpatialDomains::eDirichlet)
                {
                    // Use equation of state to evaluate temperature
                    NekDouble rho, e;
                    for(int n = 0; n < nBndEdgePts; ++n)
                    {
                        rho = fields[0]->
                                  GetBndCondExpansions()[j]->
                                  GetPhys()[id1+n];
                        e   = fields[nScalars]->
                                  GetBndCondExpansions()[j]->
                                  GetPhys()[id1 +n] / rho - tmp2[id2+n];
                        scalarVariables[nScalars-1][id2+n] =
                                m_eos->GetTemperature(rho, e);
                    }
                }

                // For Dirichlet boundary condition: uflux = u_bcs
                if (fields[nScalars]->GetBndConditions()[j]->
                    GetBoundaryConditionType() ==
                    SpatialDomains::eDirichlet &&
                    !boost::iequals(
                        fields[nScalars]->GetBndConditions()[j]
                        ->GetUserDefined(), "WallAdiabatic"))
                {
                    Vmath::Vcopy(nBndEdgePts,
                                 &scalarVariables[nScalars-1][id2], 1,
                                 &penaltyfluxO1[nScalars-1][id2], 1);

                }

                // For Neumann boundary condition: uflux = u_+
                else if (((fields[nScalars]->GetBndConditions()[j])->
                          GetBoundaryConditionType() ==
                          SpatialDomains::eNeumann) ||
                         boost::iequals(fields[nScalars]->GetBndConditions()[j]->
                                        GetUserDefined(), "WallAdiabatic"))
                {
                    Vmath::Vcopy(nBndEdgePts,
                                 &uplus[nScalars-1][id2], 1,
                                 &penaltyfluxO1[nScalars-1][id2], 1);

                }
            }
        }
    }

    /**
     * @brief Build the numerical flux for the 2nd order derivatives
     *
     */
    void DiffusionLDGNS::v_NumericalFluxO2(
        const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
        const Array<OneD, Array<OneD, NekDouble> >               &ufield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
              Array<OneD, Array<OneD, NekDouble> >               &qflux)
    {
        int i, j;
        int nTracePts = fields[0]->GetTrace()->GetTotPoints();
        int nVariables   = fields.num_elements();
        int nDim         = fields[0]->GetCoordim(0);

        Array<OneD, NekDouble > Vn (nTracePts, 0.0);

        Array<OneD, NekDouble > qFwd     (nTracePts);
        Array<OneD, NekDouble > qBwd     (nTracePts);
        Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);

        // Get the normal velocity Vn
        for(i = 0; i < nDim; ++i)
        {
            Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1,
                         Vn, 1, Vn, 1);
        }

        Array<OneD, NekDouble > qtemp(nTracePts);

        // Evaulate Riemann flux
        // qflux = \hat{q} \cdot u = q \cdot n
        // Notice: i = 1 (first row of the viscous tensor is zero)
        for (i = 1; i < nVariables; ++i)
        {
            qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
            for (j = 0; j < nDim; ++j)
            {
                // Compute qFwd and qBwd value of qfield in position 'ji'
                fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);

                // Get Riemann flux of qflux --> LDG implies upwind
                fields[i]->GetTrace()->Upwind(Vn, qBwd, qFwd, qfluxtemp);

                // Multiply the Riemann flux by the trace normals
                Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1,
                            qfluxtemp, 1);

                // Extract the physical values of the solution at the boundaries
                fields[i]->ExtractTracePhys(qfield[j][i], qtemp);

                // Impose weak boundary condition with flux
                if (fields[0]->GetBndCondExpansions().num_elements())
                {
                    v_WeakPenaltyO2(fields, i, j, qfield[j][i], qtemp, qfluxtemp);
                }

                // Store the final flux into qflux
                Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux[i], 1,
                            qflux[i], 1);
            }
        }
    }


    /**
     * @brief Imposes appropriate bcs for the 2nd order derivatives
     *
     */
    void DiffusionLDGNS::v_WeakPenaltyO2(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const int                                          var,
        const int                                          dir,
        const Array<OneD, const NekDouble>                &qfield,
        const Array<OneD, const NekDouble>                &qtemp,
              Array<OneD,       NekDouble>                &penaltyflux)
    {
        int cnt = 0;
        int nBndEdges, nBndEdgePts;
        int i, e;
        int id2;

        int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();

        // Loop on the boundary regions to apply appropriate bcs
        for (i = 0; i < nBndRegions; ++i)
        {
            // Number of boundary regions related to region 'i'
            nBndEdges = fields[var]->
            GetBndCondExpansions()[i]->GetExpSize();

            if (fields[var]->GetBndConditions()[i]->GetBoundaryConditionType()
                == SpatialDomains::ePeriodic)
            {
                continue;
            }

            // Weakly impose bcs by modifying flux values
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = fields[var]->
                GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

                id2 = fields[0]->GetTrace()->
                GetPhys_Offset(fields[0]->GetTraceMap()->
                               GetBndCondTraceToGlobalTraceMap(cnt++));

                // In case of Dirichlet bcs:
                // uflux = gD
                if(fields[var]->GetBndConditions()[i]->
                   GetBoundaryConditionType() == SpatialDomains::eDirichlet
                   && !boost::iequals(fields[var]->GetBndConditions()[i]->
                                      GetUserDefined(), "WallAdiabatic"))
                {
                    Vmath::Vmul(nBndEdgePts,
                                &m_traceNormals[dir][id2], 1,
                                &qtemp[id2], 1,
                                &penaltyflux[id2], 1);
                }
                // 3.4) In case of Neumann bcs:
                // uflux = u+
                else if((fields[var]->GetBndConditions()[i])->
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
                else if(boost::iequals(fields[var]->GetBndConditions()[i]->
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
                                    &qtemp[id2], 1,
                                    &penaltyflux[id2], 1);

                    }
                }
            }
        }
    }

    /**
     * @brief Compute primary variables
     *
     */
    void DiffusionLDGNS::v_GetPrimVar(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>>         &inarray,
              Array<OneD, Array<OneD, NekDouble>>         &primVar)
    {
        int nDim = fields[0]->GetCoordim(0);
        int nPts = fields[0]->GetTotPoints();
        for(int i = 0; i < nDim; ++i)
            {
                primVar[i] = Array<OneD, NekDouble>(nPts, 0.0);
                primVar[i] = inarray[i];
            }
    }
}//end of namespace Nektar
