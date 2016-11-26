///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionSchemeWeakPressure.cpp
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
// Description: Velocity Correction Scheme for the Incompressible
// Navier Stokes equations using Weak Pressure formulation 
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionSchemeWeakPressure.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

using namespace std;

namespace Nektar
{
    string VCSWeakPressure::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "VCSWeakPressure", 
            VCSWeakPressure::create);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    VCSWeakPressure::VCSWeakPressure(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          VelocityCorrectionScheme(pSession)
    {
        
    }

    void VCSWeakPressure::v_InitObject()
    {
        int n;
        
        IncNavierStokes::v_InitObject();
        m_explicitDiffusion = false;

        // Set m_pressure to point to last field of m_fields;
        if (boost::iequals(m_session->GetVariable(m_fields.num_elements()-1), "p"))
        {
            m_nConvectiveFields = m_fields.num_elements()-1;
            m_pressure = m_fields[m_nConvectiveFields];
        }
        else
        {
            ASSERTL0(false,"Need to set up pressure field definition");
        }

        // Determine diffusion coefficients for each field
        m_diffCoeff = Array<OneD, NekDouble> (m_nConvectiveFields, m_kinvis);
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            std::string varName = m_session->GetVariable(n);
            if ( m_session->DefinesFunction("DiffusionCoefficient", varName))
            {
                LibUtilities::EquationSharedPtr ffunc
                    = m_session->GetFunction("DiffusionCoefficient", varName);
                m_diffCoeff[n] = ffunc->Evaluate();
            }
        }

        // creation of the extrapolation object
        if(m_equationType == eUnsteadyNavierStokes)
        {
            std::string vExtrapolation = "WeakPressure";

            if (m_session->DefinesSolverInfo("Extrapolation"))
            {
                vExtrapolation = m_session->GetSolverInfo("Extrapolation");

                // if substepping is specified redefine extrapolation to SubSteppingWeakPressure
                // so that SubStepSetPressureBCs is redfined to use AddVelBC instead of AddDuDtx
                if(boost::iequals(vExtrapolation,"SubStepping"))
                {
                    vExtrapolation = "SubSteppingWeakPressure";
                }
            }

            m_extrapolation = GetExtrapolateFactory().CreateInstance(
                vExtrapolation,
                m_session,
                m_fields,
                m_pressure,
                m_velocity,
                m_advObject);
        }

        // Integrate only the convective fields
        for (n = 0; n < m_nConvectiveFields; ++n)
        {
            m_intVariables.push_back(n);
        }
        
        m_saved_aii_Dt = Array<OneD, NekDouble>(m_nConvectiveFields,
                                                NekConstants::kNekUnsetDouble);

        // Load parameters for Spectral Vanishing Viscosity
        m_session->MatchSolverInfo("SpectralVanishingViscosity","True",
                                   m_useSpecVanVisc, false);
        m_useHomo1DSpecVanVisc = m_useSpecVanVisc;
        if(m_useSpecVanVisc == false)
        {
            m_session->MatchSolverInfo("SpectralVanishingViscositySpectralHP",
                                "True", m_useSpecVanVisc, false);
            m_session->MatchSolverInfo("SpectralVanishingViscosityHomo1D",
                                "True", m_useHomo1DSpecVanVisc, false);
        }
        m_session->LoadParameter("SVVCutoffRatio",m_sVVCutoffRatio,0.75);
        m_session->LoadParameter("SVVDiffCoeff",  m_sVVDiffCoeff,  0.1);

        m_session->MatchSolverInfo("SPECTRALHPDEALIASING","True",
                                   m_specHP_dealiasing,false);

        if(m_HomogeneousType == eHomogeneous1D)
        {
            ASSERTL0(m_nConvectiveFields > 2,"Expect to have three velocity fields with homogenous expansion");

            if(m_useHomo1DSpecVanVisc)
            {
                Array<OneD, unsigned int> planes;
                planes = m_fields[0]->GetZIDs();

                int num_planes = planes.num_elements();
                Array<OneD, NekDouble> SVV(num_planes,0.0);
                NekDouble fac;
                int kmodes = m_fields[0]->GetHomogeneousBasis()->GetNumModes();
                int pstart;

                pstart = m_sVVCutoffRatio*kmodes;
                
                for(n = 0; n < num_planes; ++n)
                {
                    if(planes[n] > pstart)
                    {
                        fac = (NekDouble)((planes[n] - kmodes)*(planes[n] - kmodes))/
                            ((NekDouble)((planes[n] - pstart)*(planes[n] - pstart)));
                        SVV[n] = m_sVVDiffCoeff*exp(-fac)/m_kinvis;
                    }
                }

                for(int i = 0; i < m_velocity.num_elements(); ++i)
                {
                    m_fields[m_velocity[i]]->SetHomo1DSpecVanVisc(SVV);
                }
            }
            
        }

        // set explicit time-intregration class operators
        m_ode.DefineOdeRhs(&VCSWeakPressure::EvaluateAdvection_SetPressureBCs, this);

        m_extrapolation->SubSteppingTimeIntegration(m_intScheme->GetIntegrationMethod(), m_intScheme);
        m_extrapolation->GenerateHOPBCMap(m_session);
        
        // set implicit time-intregration class operators
        m_ode.DefineImplicitSolve(&VCSWeakPressure::SolveUnsteadyStokesSystem,this);
    }
    
    /**
     * Destructor
     */
    VCSWeakPressure::~VCSWeakPressure(void)
    {        
    }
    
    /**
     * 
     */
    void VCSWeakPressure::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);

        SolverUtils::AddSummaryItem(s, "Splitting Scheme", "Velocity correction (weak press. form)");
        if (m_extrapolation->GetSubStepIntegrationMethod() !=
            LibUtilities::eNoTimeIntegrationMethod)
        {
            SolverUtils::AddSummaryItem(s, "Substepping", 
                             LibUtilities::TimeIntegrationMethodMap[
                              m_extrapolation->GetSubStepIntegrationMethod()]);
        }

        string dealias = m_homogen_dealiasing ? "Homogeneous1D" : "";
        if (m_specHP_dealiasing)
        {
            dealias += (dealias == "" ? "" : " + ") + string("spectral/hp");
        }
        if (dealias != "")
        {
            SolverUtils::AddSummaryItem(s, "Dealiasing", dealias);
        }

        string smoothing = m_useSpecVanVisc ? "spectral/hp" : "";
        if (m_useHomo1DSpecVanVisc && (m_HomogeneousType == eHomogeneous1D))
        {
            smoothing += (smoothing == "" ? "" : " + ") + string("Homogeneous1D");
        }
        if (smoothing != "")
        {
            SolverUtils::AddSummaryItem(
                s, "Smoothing", "SVV (" + smoothing + " SVV (cut-off = "
                + boost::lexical_cast<string>(m_sVVCutoffRatio)
                + ", diff coeff = "
                + boost::lexical_cast<string>(m_sVVDiffCoeff)+")");
        }

    }
        
    /**
     * Weak Forcing term for Poisson solver solver
     */ 
    void   VCSWeakPressure::v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {                
        int ncoeffs = m_fields[0]->GetNcoeffs();
        
        m_fields[0]->IProductWRTDerivBase(fields,Forcing[0]);

        // aii required since time integration scheme normalises against aii
        Vmath::Smul(ncoeffs,-1.0/aii_Dt,Forcing[0],1,Forcing[0],1);        
    }
    
    /**
     * Solve pressure system
     */
    void   VCSWeakPressure::v_SolvePressure(
        const Array<OneD, NekDouble>  &Forcing)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficient for equation
        factors[StdRegions::eFactorLambda] = 0.0;

        // Solver Pressure Poisson Equation
        m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(),
                              NullFlagList,
                              factors, StdRegions::NullVarCoeffMap,
                              NullNekDouble1DArray, false);

        // Add presure to outflow bc if using convective like BCs
        m_extrapolation->AddPressureToOutflowBCs(m_kinvis);
    }
    
} //end of namespace
