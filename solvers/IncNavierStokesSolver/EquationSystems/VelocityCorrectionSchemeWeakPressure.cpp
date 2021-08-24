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
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          VelocityCorrectionScheme(pSession, pGraph)
    {
        
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

        if( m_extrapolation->GetSubStepName().size() )
        {
            SolverUtils::AddSummaryItem( s, "Substepping", 
                                         m_extrapolation->GetSubStepName() );
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
                              factors, StdRegions::NullVarCoeffMap,
                              MultiRegions::NullVarFactorsMap,
                              NullNekDouble1DArray, false);

        // Add presure to outflow bc if using convective like BCs
        m_extrapolation->AddPressureToOutflowBCs(m_kinvis);
    }
    
} //end of namespace
