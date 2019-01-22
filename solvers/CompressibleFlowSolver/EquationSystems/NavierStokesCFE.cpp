///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

using namespace std;

namespace Nektar
{
    string NavierStokesCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          CompressibleFlowSystem(pSession, pGraph)
    {
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        // Get gas constant from session file and compute Cp
        NekDouble gasConstant;
        m_session->LoadParameter ("GasConstant",   gasConstant,   287.058);
        m_Cp      = m_gamma / (m_gamma - 1.0) * gasConstant;
        m_Cv      = m_Cp   /  m_gamma;

        m_session->LoadParameter ("Twall", m_Twall, 300.15);

        // Create equation of state object
        std::string eosType;
        m_session->LoadSolverInfo("EquationOfState",
                                  eosType, "IdealGas");
        m_eos = GetEquationOfStateFactory()
                                .CreateInstance(eosType, m_session);

        // Viscosity
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_mu,            1.78e-05);

        // Thermal conductivity or Prandtl
        if( m_session->DefinesParameter("thermalConductivity"))
        {
            ASSERTL0( !m_session->DefinesParameter("Pr"),
                 "Cannot define both Pr and thermalConductivity.");

            m_session->LoadParameter ("thermalConductivity",
                                        m_thermalConductivity);
            m_Prandtl = m_Cp * m_mu / m_thermalConductivity;
        }
        else
        {
            m_session->LoadParameter ("Pr", m_Prandtl, 0.72);
            m_thermalConductivity = m_Cp * m_mu / m_Prandtl;
        }

        string diffName, advName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);
        if("InteriorPenalty"==diffName)
        {
            SetBoundaryConditionsBwdWeight();
        }

        if (m_specHP_dealiasing)
        {
            m_diffusion->SetFluxVectorNS(
                &NavierStokesCFE::v_GetViscousFluxVectorDeAlias,
                this);
            m_diffusion->SetDiffusionFluxCons(
                &NavierStokesCFE::GetViscousFluxVectorConservVar, this);
            m_diffusion->SetFunctorDerivBndCond(
                &NavierStokesCFE::SetBoundaryConditionsDeriv, this);
            
            m_diffusion->SetSpecialBndTreat(
                &NavierStokesCFE::SpecialBndTreat, this);

            m_diffusion->SetDiffusionSymmFluxCons(
                &NavierStokesCFE::v_GetViscousSymmtrFluxConservVar, this);

            if (m_shockCaptureType != "Off")
            {
                m_diffusion->SetArtificialDiffusionVector(
                    &NavierStokesCFE::GetArtificialViscosity, this);
            }
        }
        else
        {
            m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                                          v_GetViscousFluxVector, this);
            m_diffusion->SetDiffusionFluxCons(
                &NavierStokesCFE::GetViscousFluxVectorConservVar, this);
            m_diffusion->SetFunctorDerivBndCond(
                &NavierStokesCFE::SetBoundaryConditionsDeriv, this);

            m_diffusion->SetSpecialBndTreat(
                &NavierStokesCFE::SpecialBndTreat, this);

            m_diffusion->SetDiffusionSymmFluxCons(
                &NavierStokesCFE::v_GetViscousSymmtrFluxConservVar, this);

            if (m_shockCaptureType != "Off")
            {
                m_diffusion->SetArtificialDiffusionVector(
                    &NavierStokesCFE::GetArtificialViscosity, this);
            }
        }

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);
        
        //Only NavierStokes equation and using weakDG,LDGNS can temparary use the codes
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
	    m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if(m_useUnifiedWeakIntegration)
        {
            if(advName=="WeakDG" && ((diffName=="LDGNS")||(diffName=="InteriorPenalty")))
            {
            }
            else
            {
                m_useUnifiedWeakIntegration=false;
                if(m_session->DefinesCmdLineArgument("verbose"))
                {
                    WARNINGL0(false, "useUnifiedWeakIntegration not coded for these parameters of Diffusion");
                }
            }
        }
    }

    void NavierStokesCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints,0.0);
        }

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if("InteriorPenalty"==diffName)
        {
            if(m_BndEvaluateTime<0.0)
            {
                ASSERTL0(false, "m_BndEvaluateTime not setup");
            }
            m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff,m_BndEvaluateTime,
                                pFwd, pBwd);
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);

            for (i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints);
                inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
                inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables-2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayofArray || 
                pBwd == NullNekDoubleArrayofArray)
            {
                inFwd = NullNekDoubleArrayofArray;
                inBwd = NullNekDoubleArrayofArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd,    inFwd[0]);
                m_varConv->GetPressure(pBwd,    inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables-2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables-2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff,
                                inFwd, inBwd);

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
            }
        }
    }


    void NavierStokesCFE::v_DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int ncoeffs    = GetNcoeffs();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(ncoeffs,0.0);
        }

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if("InteriorPenalty"==diffName)
        {
            if(m_BndEvaluateTime<0.0)
            {
                ASSERTL0(false, "m_BndEvaluateTime not setup");
            }
            m_diffusion->Diffuse_coeff(nvariables, m_fields, inarray, outarrayDiff,m_BndEvaluateTime,
                                pFwd, pBwd);
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(ncoeffs,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);

            for (i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints);
                inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
                inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables-2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayofArray || 
                pBwd == NullNekDoubleArrayofArray)
            {
                inFwd = NullNekDoubleArrayofArray;
                inBwd = NullNekDoubleArrayofArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd,    inFwd[0]);
                m_varConv->GetPressure(pBwd,    inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables-2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables-2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse_coeff(nvariables, m_fields, inarrayDiff, outarrayDiff,
                                inFwd, inBwd);

            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(ncoeffs,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }

            // m_diffusion->v_Diffuse_coeff(inarray, outarray, pFwd, pBwd);
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusion_coeff(inarray, outarray);
            }
        }
    }

    void NavierStokesCFE::v_DoDiffusionFlux(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>&VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>>             &TraceFlux,
        const Array<OneD, Array<OneD, NekDouble>>       &pFwd,
        const Array<OneD, Array<OneD, NekDouble>>       &pBwd)
    {
        int nDim= m_fields[0]->GetCoordim(0);
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints,0.0);
        }
        
        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");
        if("InteriorPenalty"==diffName)
        {
            Array<OneD,Array<OneD, Array<OneD, NekDouble>>> inarrayDiffderivative(nDim);
            for (int i = 0; i < nDim; i++)
            {
                inarrayDiffderivative[i]=Array<OneD, Array<OneD, NekDouble>> (nvariables);
                for(int j=0;j<nvariables;j++)
                {
                    inarrayDiffderivative[i][j]=Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            m_diffusion->DiffuseCalculateDerivative(nvariables,m_fields,inarray,inarrayDiffderivative,pFwd,pBwd);
            m_diffusion->DiffuseVolumeFlux(nvariables, m_fields, inarray,inarrayDiffderivative, VolumeFlux);
            m_diffusion->DiffuseTraceFlux(nvariables, m_fields, inarray,inarrayDiffderivative,VolumeFlux,TraceFlux,pFwd,pBwd);
        }
        else
        {
            Array<OneD, Array<OneD, NekDouble>> inarrayDiff;
            Array<OneD, Array<OneD, NekDouble>> inFwd;
            Array<OneD, Array<OneD, NekDouble>> inBwd;
            Array<OneD,Array<OneD, Array<OneD, NekDouble>>> inarrayDiffderivative(nDim);
            
            //Allocate memory
            inarrayDiff=Array<OneD, Array<OneD, NekDouble>> (nvariables - 1);
            inFwd=Array<OneD, Array<OneD, NekDouble>>(nvariables - 1);
            inBwd=Array<OneD, Array<OneD, NekDouble>> (nvariables - 1);
            inarrayDiffderivative=Array<OneD,Array<OneD, Array<OneD, NekDouble>>> (nDim);

            for (int i = 0; i < nvariables-1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>(npoints,0.0);
                inFwd[i]       = Array<OneD, NekDouble>(nTracePts,0.0);
                inBwd[i]       = Array<OneD, NekDouble>(nTracePts,0.0);
            }

            for (int i = 0; i < nDim; i++)
            {
                inarrayDiffderivative[i]=Array<OneD, Array<OneD, NekDouble>> (nvariables-1);
                for(int j=0;j<nvariables-1;j++)
                {
                    inarrayDiffderivative[i][j]=Array<OneD, NekDouble>(npoints,0.0);
                }
            }

            // Extract pressure
            //    (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables - 2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayofArray || pBwd == NullNekDoubleArrayofArray)
            {
                inFwd = NullNekDoubleArrayofArray;
                inBwd = NullNekDoubleArrayofArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd, inFwd[0]);
                m_varConv->GetPressure(pBwd, inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables - 2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables - 2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            // To notice, needs to firstly calculate volumeflux, traceflux uses it.
            m_diffusion->DiffuseCalculateDerivative(nvariables,m_fields,inarrayDiff,inarrayDiffderivative,inFwd,inBwd);
            m_diffusion->DiffuseVolumeFlux(nvariables, m_fields, inarrayDiff,inarrayDiffderivative, VolumeFlux);
            m_diffusion->DiffuseTraceFlux(nvariables, m_fields, inarrayDiff,inarrayDiffderivative,VolumeFlux,TraceFlux, inFwd, inBwd);

            //Artificial Diffusion need to implement
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusionFlux(inarray, VolumeFlux,TraceFlux);
            }
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        int nPts       = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivativesO1,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor)
    {
        int i, j;
        int nVariables = m_fields.num_elements();
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        int nPts_orig = physfield[0].num_elements();

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > thermalConductivity(nPts, 0.0);
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            m_varConv->GetDynamicViscosity(physfield[nVariables-2], mu);
            NekDouble tRa = m_Cp / m_Prandtl;
            Vmath::Smul(nPts, tRa, mu, 1, thermalConductivity, 1);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
            Vmath::Fill(nPts, m_thermalConductivity,
                        thermalConductivity, 1);
        }

        // Interpolate inputs and initialise interpolated output
        Array<OneD, Array<OneD, NekDouble> > vel_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             deriv_interp(m_spacedim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                             out_interp(m_spacedim);
        for (i = 0; i < m_spacedim; ++i)
        {
            // Interpolate velocity
            vel_interp[i]   = Array<OneD, NekDouble> (nPts);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], vel_interp[i]);

            // Interpolate derivatives
            deriv_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+1);
            for (j = 0; j < m_spacedim+1; ++j)
            {
                deriv_interp[i][j] = Array<OneD, NekDouble> (nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j], deriv_interp[i][j]);
            }

            // Output (start from j=1 since flux is zero for rho)
            out_interp[i] = Array<OneD,Array<OneD,NekDouble> > (m_spacedim+2);
            for (j = 1; j < m_spacedim+2; ++j)
            {
                out_interp[i][j] = Array<OneD, NekDouble> (nPts);
            }
        }

        // Velocity divergence
        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, deriv_interp[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0 (no need to dealias)
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts_orig, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, deriv_interp[i][j], 1,
                                  deriv_interp[j][i], 1,
                                  out_interp[i][j+1], 1);

                Vmath::Vmul(nPts, mu, 1,
                                  out_interp[i][j+1], 1,
                                  out_interp[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, out_interp[i][j+1], 1,
                                  divVel, 1,
                                  out_interp[i][j+1], 1);
                }
                else
                {
                    // Make symmetric
                    out_interp[j][i+1] = out_interp[i][j+1];
                }
            }
        }

        // Terms for the energy equation
        for (i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, out_interp[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, vel_interp[j], 1,
                               out_interp[i][j+1], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, thermalConductivity, 1,
                               deriv_interp[i][m_spacedim], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
        }

        // Project to original space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 1; j < m_spacedim+2; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    out_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
    }

    /**
     * @brief Return the flux vector for the IP diffusion problem.
     * \TODO:: RECODE USING VARIABLE CONVERTORS
     */
    void NavierStokesCFE::GetPrimDerivFromConsDeriv(
        const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &qfields,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &physderivatives)
    {
        //can be nPoints/nTracePoints
        int nPts                = inarray[0].num_elements();
        int nConvectiveFields   = inarray.num_elements();
        int nDim                = qfields.num_elements();
        Array<OneD,NekDouble> tmp(nPts,0.0);
        int nScalars=nConvectiveFields-1;
        Array<OneD,Array<OneD,NekDouble>> physfield(nScalars);
        for(int i=0;i<nScalars;i++)
        {
            physfield[i]=Array<OneD,NekDouble> (nPts,0.0);
        }

        //Transfer conservative variables to primal variables u,v,w
        //E
        Array<OneD,NekDouble> E(nPts,0.0);
        Vmath::Vdiv(nPts,&inarray[nConvectiveFields-1][0],1,&inarray[0][0],1,&E[0],1);
        //q2
        for(int i=0;i<nScalars-1;i++)
        {
            Vmath::Vdiv(nPts,&inarray[i+1][0],1,&inarray[0][0],1,&physfield[i][0],1);
            Vmath::Vvtvp(nPts,&physfield[i][0],1,&physfield[i][0],1,&tmp[0],1,&tmp[0],1);
        }
        //e=E-0.5q2
        Vmath::Smul(nPts,0.5,&tmp[0],1,&tmp[0],1);
        Vmath::Vsub(nPts,&E[0],1,&tmp[0],1,&physfield[nScalars-1][0],1);
        //T=e/Cv
        NekDouble oCv=1./m_Cv;
        Vmath::Smul(nPts,oCv,&physfield[nScalars-1][0],1,&physfield[nScalars-1][0],1);

        //Transfer conservative variable derivatives to primal variables du,dv,dw
        Array<OneD,NekDouble> orho1(nPts,0.0);
        Array<OneD,NekDouble> orho2(nPts,0.0);
        Vmath::Sdiv(nPts,1.0,&inarray[0][0],1,&orho1[0],1);
        Vmath::Vmul(nPts,&orho1[0],1,&orho1[0],1,&orho2[0],1);
        for(int i=0;i<nDim;i++)
        {
            for(int j=0;j<nScalars-1;j++)
            {
                Vmath::Vmul(nPts,&qfields[i][0][0],1,&physfield[j][0],1,&tmp[0],1);
                Vmath::Vsub(nPts,&qfields[i][j+1][0],1,&tmp[0],1,&physderivatives[i][j][0],1);
                Vmath::Vmul(nPts,&orho1[0],1,&physderivatives[i][j][0],1,&physderivatives[i][j][0],1);
            }
        }

        //Construct dT_dx,dT_dy,dT_dz
        for(int i=0;i<nDim;i++)
        {
            Vmath::Vmul(nPts,&qfields[i][0][0],1,&E[0],1,&tmp[0],1);
            Vmath::Vsub(nPts,&qfields[i][nConvectiveFields-1][0],1,&tmp[0],1,&physderivatives[i][nScalars-1][0],1);
            Vmath::Vmul(nPts,&orho1[0],1,&physderivatives[i][nScalars-1][0],1,&physderivatives[i][nScalars-1][0],1);

            for(int j=0;j<nScalars-1;j++)
            {
               Vmath::Vmul(nPts,&physfield[j][0],1,&physderivatives[i][j][0],1,&tmp[0],1);
               Vmath::Vsub(nPts,&physderivatives[i][nScalars-1][0],1,&tmp[0],1,&physderivatives[i][nScalars-1][0],1);
            }
            Vmath::Smul(nPts,oCv,&physderivatives[i][nScalars-1][0],1,&physderivatives[i][nScalars-1][0],1);
        }
    }

      /**
     * @brief Return the flux vector for the IP diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::GetViscousFluxVectorConservVar(
        const int                                                       nConvectiveFields,
        const int                                                       nDim,
        const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &qfields,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &outarray,
              Array< OneD, int >                                        &nonZeroIndex,    
        const Array<OneD, Array<OneD, NekDouble> >                      &normal,         
        const Array<OneD, NekDouble>                                    &ArtifDiffFactor)
    {
        int nPts=inarray[0].num_elements();
        int n_nonZero   =   nConvectiveFields-1;
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxVec;
        Array<OneD,Array<OneD,NekDouble>> outtmp(nConvectiveFields);

        for(int i=0; i<nConvectiveFields;i++)
        {
            outtmp[i]=Array<OneD,NekDouble>(nPts,0.0);
        }

        for(int i=0; i< outarray.num_elements(); i++)
        {
            for(int j=0; j< nConvectiveFields; j++)
            {
                Vmath::Zero(nPts,outarray[i][j],1);
            }
        }

        if(normal.num_elements())
        {
            fluxVec =   Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(nDim);

            for(int i=0; i< nDim; i++)
            {
                fluxVec[i]  =   Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);

                for(int j=0; j< nConvectiveFields; j++)
                {
                    fluxVec[i][j]  =   Array<OneD, NekDouble>(nPts,0.0);
                }
            }

            for(int nd=0;nd<nDim;nd++)
            {
                for(int nderiv=0;nderiv<nDim;nderiv++)
                {
                    GetViscousFluxBilinearForm(nConvectiveFields,nDim,nd,nderiv,inarray,qfields[nderiv],outtmp);

                    for(int j=0;j<nConvectiveFields;j++)
                    {
                        Vmath::Vadd(nPts,&outtmp[j][0],1,&fluxVec[nd][j][0],1,&fluxVec[nd][j][0],1);
                    }
                }

                for(int j=0;j<nConvectiveFields;j++)
                {
                    Vmath::Vvtvp(nPts,&normal[nd][0],1,&fluxVec[nd][j][0],1,&outarray[0][j][0],1,&outarray[0][j][0],1);
                }
            }
        }
        else
        {
            fluxVec = outarray;

            for(int nd=0;nd<nDim;nd++)
            {
                for(int nderiv=0;nderiv<nDim;nderiv++)
                {
                    GetViscousFluxBilinearForm(nConvectiveFields,nDim,nd,nderiv,inarray,qfields[nderiv],outtmp);

                    for(int j=0;j<nConvectiveFields;j++)
                    {
                        Vmath::Vadd(nPts,&outtmp[j][0],1,&fluxVec[nd][j][0],1,&fluxVec[nd][j][0],1);
                    }
                }
            }
        }

        if(ArtifDiffFactor.num_elements())
        {
            n_nonZero   =   nConvectiveFields;
            
            if(normal.num_elements())
            {
                Array<OneD, NekDouble> tmparray(nPts,0.0);
                for(int i=0; i< nDim; i++)
                {
                    Vmath::Vmul(nPts,ArtifDiffFactor,1,normal[i],1,tmparray,1);
                    for(int j=0; j< nConvectiveFields; j++)
                    {
                        Vmath::Vvtvp(nPts,tmparray,1,qfields[i][j],1,outarray[0][j],1,outarray[0][j],1);
                    }
                }
            }
            else
            {
                for(int i=0; i< nDim; i++)
                {
                    for(int j=0; j< nConvectiveFields; j++)
                    {
                        Vmath::Vvtvp(nPts,ArtifDiffFactor,1,qfields[i][j],1,outarray[i][j],1,outarray[i][j],1);
                    }
                }
            }
        }

        nonZeroIndex = Array< OneD, int > (n_nonZero,0);
        for(int i=1;i<n_nonZero+1; i++)
        {
            nonZeroIndex[n_nonZero-i] =   nConvectiveFields-i;
        }
    }


    /**
     * @brief For very special treatment. For general boundaries it should do nothing
     * TODO: check WallViscous Boundary and Twall setting and remove the special treatment here
     *
     */
    void NavierStokesCFE::SpecialBndTreat(
        const int                                           nConvectiveFields,
              Array<OneD,       Array<OneD, NekDouble> >    &consvar)
    {            
        int ndens       = 0;
        int nengy       = nConvectiveFields-1;
        int nvelst      = ndens + 1;
        int nveled      = nengy;
        
        int cnt;
        int j, e;
        int id2;

        int nBndEdgePts, nBndEdges, nBndRegions;

        NekDouble InternalEnergy  =   m_eos->GetInternalEnergy(m_Twall);
        
        Array<OneD, NekDouble> wallTotEngy;
        int nLengthArray    =0;

        // Compute boundary conditions  for Energy
        cnt = 0;
        nBndRegions = m_fields[nengy]->
        GetBndCondExpansions().num_elements();
        for (j = 0; j < nBndRegions; ++j)
        {
            if (m_fields[nengy]->GetBndConditions()[j]->
                GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                continue;
            }

            nBndEdges = m_fields[nengy]->
            GetBndCondExpansions()[j]->GetExpSize();
            for (e = 0; e < nBndEdges; ++e)
            {
                nBndEdgePts = m_fields[nengy]->
                GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                id2 = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                            GetBndCondTraceToGlobalTraceMap(cnt++));

                // Imposing Temperature Twall at the wall 
                if (boost::iequals(m_fields[nengy]->GetBndConditions()[j]->
                    GetUserDefined(),"WallViscous"))
                {
                    if(nBndEdgePts>nLengthArray)
                    {
                        wallTotEngy     =   Array<OneD, NekDouble> (nBndEdgePts,0.0);
                        nLengthArray    =   nBndEdgePts;
                    }
                    else
                    {
                        Vmath::Fill(nLengthArray,0.0,wallTotEngy,1);
                    }

                    for(int k=nvelst;k<nveled;k++)
                    {
                        Vmath::Vvtvp(nBndEdgePts,&consvar[k][id2],1,&consvar[k][id2],1,&wallTotEngy[0],1,&wallTotEngy[0],1);
                    }
                    Vmath::Vdiv(nBndEdgePts,&wallTotEngy[0],1,&consvar[ndens][id2],1,&wallTotEngy[0],1);
                    Vmath::Svtvp(nBndEdgePts,InternalEnergy,&consvar[ndens][id2],1,&wallTotEngy[0],1,&wallTotEngy[0],1);
                    

                    Vmath::Vcopy(nBndEdgePts, 
                                &wallTotEngy[0], 1, 
                                &consvar[nengy][id2], 1);
                }                    
            }
        }
    }

    void NavierStokesCFE::GetArtificialViscosity(
        const Array<OneD, Array<OneD, NekDouble> >  &inarray,
              Array<OneD,             NekDouble  >  &muav)
    {            
        m_artificialDiffusion->GetArtificialViscosity(inarray,muav);
    }

    void NavierStokesCFE::v_GetViscousSymmtrFluxConservVar(
        const int                                                       nConvectiveFields,
        const int                                                       nSpaceDim,
        const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
        const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
        Array< OneD, int >                                              &nonZeroIndex,    
        const Array<OneD, Array<OneD, NekDouble> >                      &normals)
    {
        int nPts                = inaverg[nConvectiveFields-1].num_elements();
        nonZeroIndex = Array<OneD, int>(nConvectiveFields-1,0);
        for(int i=0;i<nConvectiveFields-1;i++)
        {
            nonZeroIndex[i] =   i+1;
        }
        Array<OneD, Array<OneD, NekDouble> > outtmp( nConvectiveFields );
        for(int i=0;i<nConvectiveFields;i++)
        {
            outtmp[i] =   Array<OneD, NekDouble> (nPts, 0.0);
        }
        for(int nd = 0; nd< nSpaceDim; nd++)
        {
            for(int nderiv =0; nderiv<nSpaceDim;nderiv++)
            {
                GetViscousFluxBilinearForm(nConvectiveFields,nSpaceDim,nd,nderiv,inaverg,inarray,outtmp);

                for(int i=0;i<nConvectiveFields;i++)
                {
                    Vmath::Vvtvp(nPts,&outtmp[i][0],1,&normals[nderiv][0],1,&outarray[nd][i][0],1,&outarray[nd][i][0],1);
                }
            }
        }
    }
    void NavierStokesCFE::GetViscousFluxBilinearForm(
        const int                                                       nConvectiveFields,
        const int                                                       nSpaceDim,
        const int                                                       FluxDirection,
        const int                                                       DerivDirection,
        const Array<OneD, const Array<OneD, NekDouble> >                &inaverg,
        const Array<OneD, const Array<OneD, NekDouble> >                &injumpp,
              Array<OneD, Array<OneD, NekDouble> >                      &outarray)
    {
        int nPts                = inaverg[nConvectiveFields-1].num_elements();
        int nDim=m_spacedim;
        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inaverg,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
        }

        Array<OneD,Array<OneD, NekDouble>> outtmp = outarray;
        for(int i=0; i<nConvectiveFields;i++)
        {
            Vmath::Zero(nPts,&outarray[i][0],1);
            // outtmp[i]=Array<OneD,NekDouble> (nPts,0.0);
        }

        //TODO: to get primary variable outside.(even in the beginning of DoODERhs)
        Array<OneD,Array<OneD,NekDouble>> u(nDim);
        Array<OneD,Array<OneD,NekDouble>> u2(nDim);
        for(int i=0;i<nDim;i++)
        {
            u[i]=Array<OneD,NekDouble> (nPts,0.0);
            u2[i]=Array<OneD,NekDouble> (nPts,0.0);
        }
        Array<OneD,NekDouble> q2(nPts,0.0);
        Array<OneD,NekDouble> E_minus_q2(nPts,0.0);
        Array<OneD,NekDouble> orho(nPts,0.0);
        Array<OneD,NekDouble>tmp(nPts,0.0);
        Array<OneD,NekDouble>tmp1(nPts,0.0);

        //Constants
        int nDim_plus_one=nDim+1;
        int FluxDirection_plus_one=FluxDirection+1;
        NekDouble gamma=m_gamma;
        NekDouble Pr=m_Prandtl;
        NekDouble gammaoPr=gamma/Pr;
        NekDouble one_minus_gammaoPr=1.0-gammaoPr;
        const NekDouble OneThird=1./3.;
        const NekDouble TwoThird=2.*OneThird;
        const NekDouble FourThird=4.*OneThird;

        Vmath::Sdiv(nPts,1.0,&inaverg[0][0],1,&orho[0],1);
        for(int i=0;i<nDim;i++)
        {
            Vmath::Vmul(nPts,&inaverg[i+1][0],1,&orho[0],1,&u[i][0],1);
            Vmath::Vmul(nPts,&u[i][0],1,&u[i][0],1,&u2[i][0],1);
            Vmath::Vadd(nPts,&q2[0],1,&u2[i][0],1,&q2[0],1);
        }
        Vmath::Vmul(nPts,&inaverg[nDim_plus_one][0],1,&orho[0],1,&E_minus_q2[0],1);
        Vmath::Vsub(nPts,&E_minus_q2[0],1,&q2[0],1,&E_minus_q2[0],1);
        Vmath::Vmul(nPts,&mu[0],1,&orho[0],1,&tmp[0],1);

        int DerivDirection_plus_one=DerivDirection+1;
        if(DerivDirection==FluxDirection)
        {
            Vmath::Svtvp(nPts,OneThird,&u2[FluxDirection][0],1,&q2[0],1,&tmp1[0],1);
            Vmath::Svtvp(nPts,gammaoPr,&E_minus_q2[0],1,&tmp1[0],1,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp1[0],1,&injumpp[0][0],1,&tmp1[0],1);
            //orho is tmperary array
            Vmath::Svtvm(nPts,gammaoPr,&injumpp[nDim_plus_one][0],1,&tmp1[0],1,&orho[0],1);

            for(int i=0;i<nDim;i++)
            {
                int i_plus_one=i+1;
                //flux[rhou,rhov,rhow]
                Vmath::Vvtvm(nPts,&u[i][0],1,&injumpp[0][0],1,&injumpp[i_plus_one][0],1,&outtmp[i_plus_one][0],1);
                Vmath::Neg(nPts,&outtmp[i_plus_one][0],1);
                Vmath::Vmul(nPts,&tmp[0],1,&outtmp[i_plus_one][0],1,&outtmp[i_plus_one][0],1);
                //flux rhoE
                Vmath::Smul(nPts,one_minus_gammaoPr,&u[i][0],1,&tmp1[0],1);
                Vmath::Vvtvp(nPts,&tmp1[0],1,&injumpp[i_plus_one][0],1,&outtmp[nDim_plus_one][0],1,&outtmp[nDim_plus_one][0],1);

                if(i==FluxDirection)
                {
                    Vmath::Smul(nPts,FourThird,&outtmp[i_plus_one][0],1,&outtmp[i_plus_one][0],1);
                    Vmath::Smul(nPts,OneThird,&u[FluxDirection][0],1,&tmp1[0],1);
                    Vmath::Vvtvp(nPts,&tmp1[0],1,&injumpp[FluxDirection_plus_one][0],1,&outtmp[nDim_plus_one][0],1,&outtmp[nDim_plus_one][0],1);
                }
            }
            Vmath::Vadd(nPts,&orho[0],1,&outtmp[nDim_plus_one][0],1,&outtmp[nDim_plus_one][0],1);
            Vmath::Vmul(nPts,&tmp[0],1,&outtmp[nDim_plus_one][0],1,&outtmp[nDim_plus_one][0],1);

        }
        else
        {
            Vmath::Vvtvm(nPts,&u[DerivDirection][0],1,&injumpp[0][0],1,&injumpp[DerivDirection_plus_one][0],1,&tmp1[0],1);
            Vmath::Smul(nPts,TwoThird,&tmp1[0],1,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp[0],1,&tmp1[0],1,&outtmp[FluxDirection_plus_one][0],1);

            Vmath::Vvtvm(nPts,&u[FluxDirection][0],1,&injumpp[0][0],1,&injumpp[FluxDirection_plus_one][0],1,&tmp1[0],1);
            Vmath::Neg(nPts,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp[0],1,&tmp1[0],1,&outtmp[DerivDirection_plus_one][0],1);

            Vmath::Smul(nPts,OneThird,&u[FluxDirection][0],1,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp1[0],1,&u[DerivDirection][0],1,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp1[0],1,&injumpp[0][0],1,&tmp1[0],1);
            //previous orho as a tmperary memory because it is non-used any more
            Vmath::Smul(nPts,TwoThird,&u[FluxDirection][0],1,&orho[0],1);
            Vmath::Vmul(nPts,&orho[0],1,&injumpp[DerivDirection_plus_one][0],1,&orho[0],1);
            Vmath::Vadd(nPts,&tmp1[0],1,&orho[0],1,&tmp1[0],1);
            Vmath::Neg(nPts,&tmp1[0],1);
            Vmath::Vvtvp(nPts,&u[DerivDirection][0],1,&injumpp[FluxDirection_plus_one][0],1,&tmp1[0],1,&tmp1[0],1);
            Vmath::Vmul(nPts,&tmp[0],1,&tmp1[0],1,&outtmp[nDim_plus_one][0],1);
            
        }
    }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

    /**
     * @brief return part of viscous Jacobian: 
     * \todo flux derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhoE_dx] 
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesCFE::GetdFlux_dQx_2D( 
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U, 
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble u=U[1]/U[0];
        NekDouble v=U[2]/U[0];
        NekDouble E=U[3]/U[0];
        NekDouble q2=u*u+v*v;
        NekDouble e=E-0.5*q2;
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Pr= Cp *mu / kappa;
        //To notice, here is positive, which is consistent with 
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu/rho;

        // TODO: replace repeated divid by multiply
        (*OutputMatrix)(0,0)=tmp*(-4.0/3.0*u*nx-v*ny);
        (*OutputMatrix)(0,1)=tmp*(4.0/3.0*nx);
        (*OutputMatrix)(0,2)=tmp*ny;
        (*OutputMatrix)(0,3)=0.0;
        (*OutputMatrix)(1,0)=tmp*(-v*nx+2.0/3.0*u*ny);
        (*OutputMatrix)(1,1)=tmp*(-2.0/3.0*ny);
        (*OutputMatrix)(1,2)=tmp*nx;
        (*OutputMatrix)(1,3)=0.0;
        (*OutputMatrix)(2,0)=(4.0/3.0*u*u+v*v+gamma/Pr*(E-q2))*nx+1.0/3.0*u*v*ny;
        (*OutputMatrix)(2,0)=-tmp*(*OutputMatrix)(2,0);
        (*OutputMatrix)(2,1)=(4.0/3.0-gamma/Pr)*u*nx-2.0/3.0*v*ny;
        (*OutputMatrix)(2,1)=tmp*(*OutputMatrix)(2,1);
        (*OutputMatrix)(2,2)=(1-gamma/Pr)*v*nx+u*ny;
        (*OutputMatrix)(2,2)=tmp*(*OutputMatrix)(2,2);
        (*OutputMatrix)(2,3)=tmp*gamma/Pr*nx;
    }

     /**
     * @brief return part of viscous Jacobian: 
     * \todo flux derived with Qx=[drho_dy,drhou_dy,drhov_dy,drhoE_dy] 
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhoE]
     * Output: 2D 3*4 Matrix (flux with rho is zero)
     */
    void NavierStokesCFE::GetdFlux_dQy_2D( 
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U, 
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble rho=U[0];
        NekDouble u=U[1]/U[0];
        NekDouble v=U[2]/U[0];
        NekDouble E=U[3]/U[0];
        NekDouble q2=u*u+v*v;
        NekDouble e=E-0.5*q2;
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Pr= Cp *mu / kappa;
        //To notice, here is positive, which is consistent with 
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I Do like CFD"
        NekDouble tmp=mu/rho;
           
        (*OutputMatrix)(0,0)=tmp*(2.0/3.0*v*nx-u*ny);
        (*OutputMatrix)(0,1)=tmp*ny;
        (*OutputMatrix)(0,2)=tmp*(-2.0/3.0)*nx;
        (*OutputMatrix)(0,3)=0.0;
        (*OutputMatrix)(1,0)=tmp*(-u*nx-4.0/3.0*v*ny);
        (*OutputMatrix)(1,1)=tmp*nx;
        (*OutputMatrix)(1,2)=tmp*(4.0/3.0*ny);
        (*OutputMatrix)(1,3)=0.0;
        (*OutputMatrix)(2,0)=1.0/3.0*u*v*nx+(4.0/3.0*v*v+u*u+gamma/Pr*(E-q2))*ny;
        (*OutputMatrix)(2,0)=-tmp*(*OutputMatrix)(2,0);
        (*OutputMatrix)(2,1)=(1-gamma/Pr)*u*ny+v*nx;
        (*OutputMatrix)(2,1)=tmp*(*OutputMatrix)(2,1);
        (*OutputMatrix)(2,2)=(4.0/3.0-gamma/Pr)*v*ny-2.0/3.0*u*nx;
        (*OutputMatrix)(2,2)=tmp*(*OutputMatrix)(2,2);
        (*OutputMatrix)(2,3)=tmp*gamma/Pr*ny;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qx=[drho_dx,drhou_dx,drhov_dx,drhow_dx,drhoE_dx]
     * Output: 3D 4*5 Matrix (flux about rho is zero) 
     * OutputMatrix(dir=0)= dF_dQx;
     */
    void NavierStokesCFE::GetdFlux_dQx_3D( 
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U, 
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble u=U[1]/U[0];
        NekDouble v=U[2]/U[0];
        NekDouble w=U[3]/U[0];
        NekDouble E=U[4]/U[0];
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Pr= Cp *mu / kappa;
        //To notice, here is positive, which is consistent with 
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu/rho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;

        (*OutputMatrix)(0,0)=tmpx*(-4.0/3.0*u)+tmpy*(-v)+tmpz*(-w);
        (*OutputMatrix)(0,1)=tmpx*(4.0/3.0);
        (*OutputMatrix)(0,2)=tmpy;
        (*OutputMatrix)(0,3)=tmpz;
        (*OutputMatrix)(0,4)=0.0;
        (*OutputMatrix)(1,0)=tmpx*(-v)+tmpy*(2./3.*u);
        (*OutputMatrix)(1,1)=tmpy*(-2./3.);
        (*OutputMatrix)(1,2)=tmpx;
        (*OutputMatrix)(1,3)=0.0;
        (*OutputMatrix)(1,4)=0.0;
        (*OutputMatrix)(2,0)=tmpx*(-w)+tmpz*(2./3.*u);
        (*OutputMatrix)(2,1)=tmpz*(-2./3);
        (*OutputMatrix)(2,2)=0.0;
        (*OutputMatrix)(2,3)=tmpx;
        (*OutputMatrix)(2,4)=0.0;
        (*OutputMatrix)(3,0)=-tmpx*(4./3.*u*u+v*v+w*w+gamma/Pr*(E-q2))+tmpy*(-1./3.*u*v)+tmpz*(-1./3.*u*w);
        (*OutputMatrix)(3,1)=tmpx*(4./3.-gamma/Pr)*u+tmpy*(-2./3.*v)+tmpz*(-2./3.*w);
        (*OutputMatrix)(3,2)=tmpx*(1.0-gamma/Pr)*v+tmpy*u;
        (*OutputMatrix)(3,3)=tmpx*(1.0-gamma/Pr)*w+tmpz*u;
        (*OutputMatrix)(3,4)=tmpx*gamma/Pr;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qy=[drho_dy,drhou_dy,drhov_dy,drhow_dy,drhoE_dy]
     * Output: 3D 4*5 Matrix (flux about rho is zero) 
     * OutputMatrix(dir=1)= dF_dQy;
     */
    void NavierStokesCFE::GetdFlux_dQy_3D( 
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U, 
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble u=U[1]/U[0];
        NekDouble v=U[2]/U[0];
        NekDouble w=U[3]/U[0];
        NekDouble E=U[4]/U[0];
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Pr= Cp *mu / kappa;
        //To notice, here is positive, which is consistent with 
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu/rho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;

        (*OutputMatrix)(0,0)=tmpx*(2./3.*v)+tmpy*(-u);
        (*OutputMatrix)(0,1)=tmpy;
        (*OutputMatrix)(0,2)=tmpx*(-2./3.);
        (*OutputMatrix)(0,3)=0.0;
        (*OutputMatrix)(0,4)=0.0;
        (*OutputMatrix)(1,0)=tmpx*(-u)+tmpy*(-4./3.*v)+tmpz*(-w);
        (*OutputMatrix)(1,1)=tmpx;
        (*OutputMatrix)(1,2)=tmpy*(4./3.);
        (*OutputMatrix)(1,3)=tmpz;
        (*OutputMatrix)(1,4)=0.0;
        (*OutputMatrix)(2,0)=tmpy*(-w)+tmpz*(2./3.*v);
        (*OutputMatrix)(2,1)=0.0;
        (*OutputMatrix)(2,2)=tmpz*(-2./3.);
        (*OutputMatrix)(2,3)=tmpy;
        (*OutputMatrix)(2,4)=0.0;
        (*OutputMatrix)(3,0)=tmpx*(-1./3.*u*v)-tmpy*(u*u+4./3.*v*v+w*w+gamma/Pr*(E-q2))+tmpz*(-1./3.*v*w);
        (*OutputMatrix)(3,1)=tmpx*v+tmpy*(1-gamma/Pr)*u;
        (*OutputMatrix)(3,2)=tmpx*(-2./3.*u)+tmpy*(4./3.-gamma/Pr)*v+tmpz*(-2./3.*w);
        (*OutputMatrix)(3,3)=tmpy*(1-gamma/Pr)*w+tmpz*v;
        (*OutputMatrix)(3,4)=tmpy*gamma/Pr;
    }

    /**
     * @brief return part of viscous Jacobian derived with Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Input:
     * normals:Point normals
     * U=[rho,rhou,rhov,rhow,rhoE]
     * dir: means whether derive with
     * Qz=[drho_dz,drhou_dz,drhov_dz,drhow_dz,drhoE_dz]
     * Output: 3D 4*5 Matrix (flux about rho is zero) 
     * OutputMatrix(dir=2)= dF_dQz;
     */
    void NavierStokesCFE::GetdFlux_dQz_3D( 
        const Array<OneD, NekDouble>    &normals,
        const NekDouble                 &mu,
        const Array<OneD, NekDouble>    &U, 
              DNekMatSharedPtr          &OutputMatrix )
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble nz=normals[2];
        NekDouble rho=U[0];
        NekDouble u=U[1]/U[0];
        NekDouble v=U[2]/U[0];
        NekDouble w=U[3]/U[0];
        NekDouble E=U[4]/U[0];
        NekDouble q2=u*u+v*v+w*w;
        NekDouble e=E-0.5*q2;
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble T=e/Cv;
        //q_x=-kappa*dT_dx;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Pr= Cp *mu / kappa;
        //To notice, here is positive, which is consistent with 
        //"SYMMETRIC INTERIOR PENALTY DG METHODS FOR THE COMPRESSIBLE NAVIER-STOKES EQUATIONS"
        //But opposite to "I do like CFD"
        NekDouble tmp=mu/rho;
        NekDouble tmpx=tmp*nx;
        NekDouble tmpy=tmp*ny;
        NekDouble tmpz=tmp*nz;

        (*OutputMatrix)(0,0)=tmpx*(2./3.*w)+tmpz*(-u);
        (*OutputMatrix)(0,1)=tmpz;
        (*OutputMatrix)(0,2)=0.0;
        (*OutputMatrix)(0,3)=tmpx*(-2./3.);
        (*OutputMatrix)(0,4)=0.0;
        (*OutputMatrix)(1,0)=tmpy*(2./3.*w)+tmpz*(-v);
        (*OutputMatrix)(1,1)=0.0;
        (*OutputMatrix)(1,2)=tmpz;
        (*OutputMatrix)(1,3)=tmpy*(-2./3.);
        (*OutputMatrix)(1,4)=0.0;
        (*OutputMatrix)(2,0)=tmpx*(-u)+tmpy*(-v)+tmpz*(-4./3.*w);
        (*OutputMatrix)(2,1)=tmpx;
        (*OutputMatrix)(2,2)=tmpy;
        (*OutputMatrix)(2,3)=tmpz*(4./3.);
        (*OutputMatrix)(2,4)=0.0;
        (*OutputMatrix)(3,0)=tmpx*(-1./3.*u*w)+tmpy*(-1./3.*v*w)-tmpz*(u*u+v*v+4./3.*w*w+gamma/Pr*(E-q2));
        (*OutputMatrix)(3,1)=tmpx*w+tmpz*(1-gamma/Pr)*u;
        (*OutputMatrix)(3,2)=tmpy*w+tmpz*(1-gamma/Pr)*v;
        (*OutputMatrix)(3,3)=tmpx*(-2./3.*u)+tmpy*(-2./3.*v)+tmpz*(4./3.-gamma/Pr)*w;
        (*OutputMatrix)(3,4)=tmpz*gamma/Pr;
    }

    /**
     * @brief return part of viscous Jacobian 
     * Input:
     * normals:Point normals
     * mu: dynamicviscosity
     * dmu_dT: mu's derivative with T using Sutherland's law
     * U=[rho,rhou,rhov,rhoE]
     * Output: 3*4 Matrix (the flux about rho is zero)
     * OutputMatrix dFLux_dU,  the matrix sign is consistent with SIPG
     */
    void NavierStokesCFE::GetdFlux_dU_2D(
        const Array<OneD, NekDouble>                        &normals, 
        const NekDouble                                     mu, 
        const NekDouble                                     dmu_dT,
        const Array<OneD, NekDouble>                        &U,
        const Array<OneD, const Array<OneD, NekDouble> >    &qfield,
              DNekMatSharedPtr                              &OutputMatrix)
    {
        NekDouble nx=normals[0];
        NekDouble ny=normals[1];
        NekDouble U1=U[0];
        NekDouble U2=U[1];
        NekDouble U3=U[2];
        NekDouble U4=U[3];
        NekDouble dU1_dx=qfield[0][0];
        NekDouble dU2_dx=qfield[0][1];
        NekDouble dU3_dx=qfield[0][2];
        NekDouble dU4_dx=qfield[0][3];
        NekDouble dU1_dy=qfield[1][0];
        NekDouble dU2_dy=qfield[1][1];
        NekDouble dU3_dy=qfield[1][2];
        NekDouble dU4_dy=qfield[1][3];
        NekDouble R =m_varConv->GetGasconstant();
        NekDouble gamma=m_gamma;
        NekDouble kappa=m_thermalConductivity;
        NekDouble Cp=gamma / (gamma - 1.0) *R;
        NekDouble Cv=1.0/(gamma-1)*R;
        NekDouble Pr= Cp *mu / kappa;

        NekDouble orho1,orho2,orho3,orho4;
        NekDouble oCv=1./Cv;
        orho1=1.0/U1;
        orho2=1.0/(U1*U1);
        orho3=orho1*orho2;
        orho4=orho2*orho2;
        
        //Assume Fn=mu*Sn
        //Sn=Sx*nx+Sy*ny

        //Term1 mu's derivative with U: dmu_dU*Sn
        NekDouble u=U2/U1;
        NekDouble v=U3/U1;
        NekDouble du_dx=orho1*(dU2_dx-u*dU1_dx);
        NekDouble dv_dy=orho1*(dU3_dy-v*dU1_dy);
        NekDouble du_dy=orho1*(dU2_dy-u*dU1_dy);
        NekDouble dv_dx=orho1*(dU3_dx-v*dU1_dx);
        NekDouble s12=4./3.*du_dx-2./3.*dv_dy;
        NekDouble s13=du_dy+dv_dx;
        NekDouble s22=s13;
        NekDouble s23=4./3.*dv_dy-2./3.*du_dx;
        NekDouble snx=s12*nx+s22*ny;
        NekDouble sny=s13*nx+s23*ny;
        NekDouble snv=snx*u+sny*v;
        NekDouble qx=-gamma*mu/Pr*(orho1*dU4_dx-U[3]*orho2*dU1_dx-u*(orho1*dU2_dx-U[1]*orho2*dU1_dx)-v*(orho1*dU3_dx-U[2]*orho2*dU1_dx));
        NekDouble qy=-gamma*mu/Pr*(orho1*dU4_dy-U[3]*orho2*dU1_dy-u*(orho1*dU2_dy-U[1]*orho2*dU1_dy)-v*(orho1*dU3_dy-U[2]*orho2*dU1_dy));
        NekDouble qn=qx*nx+qy*ny;
        Array<OneD,NekDouble> tmp(3,0.0);
        tmp[0]=snx;
        tmp[1]=sny;
        tmp[2]=snv-qn/mu;
        Array<OneD,NekDouble> dT_dU (4,0.0);
        dT_dU[0]=oCv*(-orho2*U4+orho3*U2*U2+orho3*U3*U3);
        dT_dU[1]=-oCv*orho2*U2;   
        dT_dU[2]=-oCv*orho2*U3;
        dT_dU[3]=oCv*orho1;
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<4;j++)
            {
                (*OutputMatrix)(i,j)=dmu_dT*dT_dU[j]*tmp[i];
            }
        }
        
        //Term 2 +mu*dSn_dU
        NekDouble du_dx_dU1,du_dx_dU2;
        NekDouble du_dy_dU1,du_dy_dU2;
        NekDouble dv_dx_dU1,dv_dx_dU3;
        NekDouble dv_dy_dU1,dv_dy_dU3;
        NekDouble ds12_dU1,ds12_dU2,ds12_dU3;
        NekDouble ds13_dU1,ds13_dU2,ds13_dU3;
        NekDouble ds22_dU1,ds22_dU2,ds22_dU3;
        NekDouble ds23_dU1,ds23_dU2,ds23_dU3;
        NekDouble dsnx_dU1,dsnx_dU2,dsnx_dU3;
        NekDouble dsny_dU1,dsny_dU2,dsny_dU3;
        NekDouble dsnv_dU1,dsnv_dU2,dsnv_dU3;

        du_dx_dU1=-orho2*dU2_dx+2*orho3*U2*dU1_dx;
        du_dx_dU2=-orho2*dU1_dx;
        du_dy_dU1=-orho2*dU2_dy+2*orho3*U2*dU1_dy;
        du_dy_dU2=-orho2*dU1_dy;
        dv_dx_dU1=-orho2*dU3_dx+2*orho3*U3*dU1_dx;
        dv_dx_dU3=du_dx_dU2;
        dv_dy_dU1=-orho2*dU3_dy+2*orho3*U3*dU1_dy;
        dv_dy_dU3=du_dy_dU2;
        ds12_dU1=4./3.*du_dx_dU1-2./3.*dv_dy_dU1;
        ds12_dU2=4./3.*du_dx_dU2;
        ds12_dU3=-2./3.*dv_dy_dU3;
        ds13_dU1=du_dy_dU1+dv_dx_dU1;
        ds13_dU2=du_dy_dU2;
        ds13_dU3=dv_dx_dU3;
        ds22_dU1=ds13_dU1;
        ds22_dU2=ds13_dU2;
        ds22_dU3=ds13_dU3;
        ds23_dU1=4./3.*dv_dy_dU1-2./3.*du_dx_dU1;
        ds23_dU2=-2./3.*du_dx_dU2;
        ds23_dU3=4./3.*dv_dy_dU3;
        dsnx_dU1=ds12_dU1*nx+ds22_dU1*ny;
        dsnx_dU2=ds12_dU2*nx+ds22_dU2*ny;
        dsnx_dU3=ds12_dU3*nx+ds22_dU3*ny;
        dsny_dU1=ds13_dU1*nx+ds23_dU1*ny;
        dsny_dU2=ds13_dU2*nx+ds23_dU2*ny;
        dsny_dU3=ds13_dU3*nx+ds23_dU3*ny;
        dsnv_dU1=u*dsnx_dU1+v*dsny_dU1-orho2*U2*snx-orho2*U3*sny;
        dsnv_dU2=u*dsnx_dU2+v*dsny_dU2+orho1*snx;
        dsnv_dU3=u*dsnx_dU3+v*dsny_dU3+orho1*sny;
        (*OutputMatrix)(0,0)=(*OutputMatrix)(0,0)+mu*dsnx_dU1;
        (*OutputMatrix)(0,1)=(*OutputMatrix)(0,1)+mu*dsnx_dU2;
        (*OutputMatrix)(0,2)=(*OutputMatrix)(0,2)+mu*dsnx_dU3;
        (*OutputMatrix)(1,0)=(*OutputMatrix)(1,0)+mu*dsny_dU1;
        (*OutputMatrix)(1,1)=(*OutputMatrix)(1,1)+mu*dsny_dU2;
        (*OutputMatrix)(1,2)=(*OutputMatrix)(1,2)+mu*dsny_dU3;
        (*OutputMatrix)(2,0)=(*OutputMatrix)(2,0)+mu*dsnv_dU1;
        (*OutputMatrix)(2,1)=(*OutputMatrix)(2,1)+mu*dsnv_dU2;
        (*OutputMatrix)(2,2)=(*OutputMatrix)(2,2)+mu*dsnv_dU3;

        //Consider qn's effect (does not include mu's effect)
        NekDouble dqx_dU1,dqx_dU2,dqx_dU3,dqx_dU4;
        NekDouble dqy_dU1,dqy_dU2,dqy_dU3,dqy_dU4;
        NekDouble tmpx=-nx*mu*gamma/Pr;
        dqx_dU1=tmpx*(-orho2*dU4_dx+2*orho3*U4*dU1_dx+2*orho3*U2*dU2_dx-3*orho4*U2*U2*dU1_dx+2*orho3*U3*dU3_dx-3*orho4*U3*U3*dU1_dx);
        dqx_dU2=tmpx*(-orho2*dU2_dx+2*orho3*U2*dU1_dx);
        dqx_dU3=tmpx*(-orho2*dU3_dx+2*orho3*U3*dU1_dx);
        dqx_dU4=-tmpx*orho2*dU1_dx;
        NekDouble tmpy=-ny*mu*gamma/Pr;
        dqy_dU1=tmpy*(-orho2*dU4_dy+2*orho3*U4*dU1_dy+2*orho3*U2*dU2_dy-3*orho4*U2*U2*dU1_dy+2*orho3*U3*dU3_dy-3*orho4*U3*U3*dU1_dy);
        dqy_dU2=tmpy*(-orho2*dU2_dy+2*orho3*U2*dU1_dy);
        dqy_dU3=tmpy*(-orho2*dU3_dy+2*orho3*U3*dU1_dy);
        dqy_dU4=-tmpy*orho2*dU1_dy;
        (*OutputMatrix)(2,0)=(*OutputMatrix)(2,0)-dqx_dU1-dqy_dU1;
        (*OutputMatrix)(2,1)=(*OutputMatrix)(2,1)-dqx_dU2-dqy_dU2;
        (*OutputMatrix)(2,2)=(*OutputMatrix)(2,2)-dqx_dU3-dqy_dU3;
        (*OutputMatrix)(2,3)=(*OutputMatrix)(2,3)-dqx_dU4-dqy_dU4;
    }

    void NavierStokesCFE::v_MinusDiffusionFluxJacDirctn(
        const int                                                       nDirctn,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
              Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        int nPts            = m_fields[0]->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();  
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);
        Array<OneD, NekDouble > DmuDT              (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
            m_varConv->GetDmuDT(temperature,mu,DmuDT);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
        }

        NekDouble pointmu       = 0.0;
        NekDouble pointDmuDT    = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> locDmuDT;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > pointDerv(nSpaceDim);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > locDerv(nSpaceDim);
        for(int j = 0; j < nSpaceDim; j++)
        {   
            pointDerv[j] = Array<OneD, NekDouble>(nConvectiveFields,0.0);
            locDerv[j]   = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
        }

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
            int noffest             = GetPhys_Offset(nelmt);
                 
            for(int j = 0; j < nConvectiveFields; j++)
            {   
                locVars[j] = inarray[j]+noffest;
            }

            for(int j = 0; j < nSpaceDim; j++)
            {   
                for(int k = 0; k < nConvectiveFields; k++)
                {
                    locDerv[j][k] = qfields[j][k]+noffest;
                }
            }
            locmu       =   mu      + noffest;
            locDmuDT    =   DmuDT   + noffest;
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }
                for(int j = 0; j < nSpaceDim; j++)
                {   
                    for(int k = 0; k < nConvectiveFields; k++)
                    {
                        pointDerv[j][k] = locDerv[j][k][npnt];
                    }
                }

                pointmu     = locmu[npnt];
                pointDmuDT  = locDmuDT[npnt];

                GetDiffusionFluxJacPoint(nelmt,pointVar,pointDerv,pointmu,pointDmuDT,normals,PointFJac);
                for (int i =0; i < nConvectiveFields-1; i++)
                {
                    for (int j =0; j < nConvectiveFields; j++)
                    {
                        (*ElmtJac[nelmt][npnt])(i+1,j) -=  (*PointFJac)(i,j);
                    }
                }
            }
        }
    }

    // TODO: currently only for 2D
    void NavierStokesCFE::v_GetDiffusionFluxJacPoint(
            const int                                           nelmt,
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals, 
                 DNekMatSharedPtr                               &fluxJac)
    {
        GetdFlux_dU_2D(normals,mu,DmuDT,conservVar,conseDeriv,fluxJac);
    }

    void NavierStokesCFE::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr                            &explist,
        const Array<OneD, const Array<OneD, NekDouble> >                &normals,
        const int                                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
              Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    explist->GetExp();
        int ntotElmt            = (*expvect).size();
        int nPts                = explist->GetTotPoints();
        int nSpaceDim           = m_graph->GetSpaceDimension();  
        
        //Debug
        if(!ElmtJac.num_elements())
        {
            ElmtJac =   Array<OneD, Array<OneD, DNekMatSharedPtr> > (ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
                ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nElmtPnt);
                for(int npnt = 0; npnt < nElmtPnt; npnt++)
                {
                    ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
                }
            }
        }
        // Auxiliary variables
        Array<OneD, NekDouble > mu                 (nPts, 0.0);

        // Variable viscosity through the Sutherland's law
        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble > temperature        (nPts, 0.0);
            m_varConv->GetTemperature(inarray,temperature);
            m_varConv->GetDynamicViscosity(temperature, mu);
        }
        else
        {
            Vmath::Fill(nPts, m_mu, mu, 1);
        }

        NekDouble pointmu       = 0.0;
        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);
        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);
        Array<OneD, NekDouble> pointnormals(nSpaceDim,0.0);
        Array<OneD, Array<OneD, NekDouble> > locnormal(nSpaceDim);

        DNekMatSharedPtr PointFJac = MemoryManager<DNekMat>
                                ::AllocateSharedPtr(nConvectiveFields-1, nConvectiveFields);

        switch (nDervDir)
        {
            case 0:
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
                    int noffest             = explist->GetPhys_Offset(nelmt);
                        
                    for(int j = 0; j < nConvectiveFields; j++)
                    {   
                        locVars[j] = inarray[j]+noffest;
                    }

                    for(int j = 0; j < nSpaceDim; j++)
                    {   
                        locnormal[j] = normals[j]+noffest;
                    }

                    locmu       =   mu      + noffest;
                    for(int npnt = 0; npnt < nElmtPnt; npnt++)
                    {
                        for(int j = 0; j < nConvectiveFields; j++)
                        {
                            pointVar[j] = locVars[j][npnt];
                        }
                        for(int j = 0; j < nSpaceDim; j++)
                        {   
                            pointnormals[j] = locnormal[j][npnt];
                        }

                        pointmu     = locmu[npnt];

                        GetdFlux_dQx_2D(pointnormals,pointmu,pointVar,PointFJac);
                        for (int i =0; i < nConvectiveFields-1; i++)
                        {
                            for (int j =0; j < nConvectiveFields; j++)
                            {
                                (*ElmtJac[nelmt][npnt])(i+1,j) =  (*PointFJac)(i,j);
                            }
                        }
                        // (*ElmtJac[nelmt][npnt]) =   (*PointFJac);
                    }
                }
                break;
            }
            case 1:
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
                    int noffest             = explist->GetPhys_Offset(nelmt);
                        
                    for(int j = 0; j < nConvectiveFields; j++)
                    {   
                        locVars[j] = inarray[j]+noffest;
                    }

                    for(int j = 0; j < nSpaceDim; j++)
                    {   
                        locnormal[j] = normals[j]+noffest;
                    }

                    locmu       =   mu      + noffest;
                    for(int npnt = 0; npnt < nElmtPnt; npnt++)
                    {
                        for(int j = 0; j < nConvectiveFields; j++)
                        {
                            pointVar[j] = locVars[j][npnt];
                        }
                        for(int j = 0; j < nSpaceDim; j++)
                        {   
                            pointnormals[j] = locnormal[j][npnt];
                        }

                        pointmu     = locmu[npnt];
                        
                        GetdFlux_dQy_2D(pointnormals,pointmu,pointVar,PointFJac);
                        for (int i =0; i < nConvectiveFields-1; i++)
                        {
                            for (int j =0; j < nConvectiveFields; j++)
                            {
                                (*ElmtJac[nelmt][npnt])(i+1,j) =  (*PointFJac)(i,j);
                            }
                        }
                        // (*ElmtJac[nelmt][npnt]) =   (*PointFJac);
                    }
                }
                break;
            }
            default:
            {
                ASSERTL0(false, "only for 2D");
                break;
            }
        }
        return;
    }
    
    void NavierStokesCFE::v_CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
    {
        int nConvectiveFields = m_fields.num_elements();
        int npoints           = GetTotPoints();
        const Array<OneD, Array<OneD, NekDouble> >                  pFwd;
        const Array<OneD, Array<OneD, NekDouble> >                  pBwd;
        if(!qfield.num_elements())
        {
            qfield  =   Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >(m_spacedim);
            for(int i = 0; i< m_spacedim; i++)
            {
                qfield[i]   =   Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for(int j = 0; j< nConvectiveFields; j++)
                {
                    qfield[i][j]   =   Array<OneD, NekDouble>(npoints,0.0);
                }
            }
        }
        m_diffusion->DiffuseCalculateDerivative(nConvectiveFields,m_fields,inarray,qfield,pFwd,pBwd);
    }
#endif

}
