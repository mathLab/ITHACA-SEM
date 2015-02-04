///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrection.cpp
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
// Description: Velocity Correction Scheme w/ coordinate transformation 
// for the Incompressible Navier Stokes equations
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VCSMapping.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

namespace Nektar
{
    string VCSMapping::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "VCSMapping", 
            VCSMapping::create);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    VCSMapping::VCSMapping(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession),
          VelocityCorrectionScheme(pSession)  
    {

    }

    void VCSMapping::v_InitObject()
    {
        VelocityCorrectionScheme::v_InitObject();
        
        m_mapping = SolverUtils::Mapping::Load(m_session, m_fields); 
    }
    
    /**
     * Destructor
     */
    VCSMapping::~VCSMapping(void)
    {        
    }
    


    /**
     * Explicit part of the method - Advection, Forcing + HOPBCs
     */
    void VCSMapping::v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &outarray, 
        const NekDouble time)
    {
        EvaluateAdvectionTerms(inarray, outarray);

        // Smooth advection
        if(m_SmoothAdvection)
        {
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_pressure->SmoothField(outarray[i]);
            }
        }

        // Add forcing terms
        std::vector<SolverUtils::ForcingSharedPtr>::const_iterator x;
        for (x = m_forcing.begin(); x != m_forcing.end(); ++x)
        {
            (*x)->Apply(m_fields, inarray, outarray, time);
        }
        
        // Add mapping terms
        ApplyIncNSMappingForcing( outarray );
        
        // Calculate High-Order pressure boundary conditions
        m_extrapolation->EvaluatePressureBCs(inarray,outarray,m_kinvis);
    }
    
        
    /**
     * Forcing term for Poisson solver solver
     */ 
    void   VCSMapping::v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &fields, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {                
        int i;
        int physTot = m_fields[0]->GetTotPoints();
        int nvel = m_velocity.num_elements();
        Array<OneD, NekDouble> wk(physTot, 0.0);
        
        Vmath::Zero(physTot,Forcing[0],1);
        
        for(i = 0; i < nvel; ++i)
        {
            m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],fields[i], wk);
            Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
        }

        Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1);        
    }
    
    /**
     * Forcing term for Helmholtz solver
     */
    void   VCSMapping::v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {
        NekDouble aii_dtinv = 1.0/aii_Dt;
        int phystot = m_fields[0]->GetTotPoints();

        // Grad p
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),m_pressure->UpdatePhys());
        
        int nvel = m_velocity.num_elements();
        if(nvel == 2)
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1]);
        }
        else
        {
            m_pressure->PhysDeriv(m_pressure->GetPhys(), Forcing[0], Forcing[1],
                                  Forcing[2]);
        }
        
        // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
        // need to be updated for the convected fields.
        for(int i = 0; i < nvel; ++i)
        {
            Blas::Daxpy(phystot,-aii_dtinv,inarray[i],1,Forcing[i],1);
            Blas::Dscal(phystot,1.0/m_kinvis,&(Forcing[i])[0],1);
        }
    }
    
    /**
     * Solve pressure system
     */
    void   VCSMapping::v_SolvePressure(
        const Array<OneD, NekDouble>  &Forcing)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficient for equation
        factors[StdRegions::eFactorLambda] = 0.0;

        // Solver Pressure Poisson Equation
        m_pressure->HelmSolve(Forcing, m_pressure->UpdateCoeffs(), NullFlagList,
                              factors);
    }
    
    /**
     * Solve velocity system
     */
    void   VCSMapping::v_SolveViscous(
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficients for equation
        factors[StdRegions::eFactorLambda] = 1.0/aii_Dt/m_kinvis;
        if(m_useSpecVanVisc)
        {
            factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
            factors[StdRegions::eFactorSVVDiffCoeff]   = m_sVVDiffCoeff/m_kinvis;
        }

        // Solve Helmholtz system and put in Physical space
        for(int i = 0; i < m_nConvectiveFields; ++i)
        {
            m_fields[i]->HelmSolve(Forcing[i], m_fields[i]->UpdateCoeffs(),
                                   NullFlagList, factors);
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
        }
    }
    
    /**
     * Explicit terms of the mapping
     */
    void   VCSMapping::ApplyIncNSMappingForcing(
        Array<OneD, Array<OneD, NekDouble> >          &outarray)
    {
        int physTot = m_fields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> >          vel(m_nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> >          Forcing(m_nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> >          tmp(m_nConvectiveFields);
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            vel[i] = Array<OneD, NekDouble> (physTot, 0.0);
            Forcing[i] = Array<OneD, NekDouble> (physTot, 0.0);
            tmp[i] = Array<OneD, NekDouble> (physTot, 0.0);
        }
        Array<OneD, NekDouble>                        P (physTot, 0.0);
        
        if(m_fields[0]->GetWaveSpace())
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[0]->HomogeneousBwdTrans(m_fields[i]->GetPhys(),vel[i]);
            }
            m_fields[0]->HomogeneousBwdTrans(m_fields[m_nConvectiveFields]->GetPhys(),P);
        }
        else
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vcopy(physTot, m_fields[i]->GetPhys(), 1, vel[i], 1);
            }
            Vmath::Vcopy(physTot, m_fields[m_nConvectiveFields]->GetPhys(), 1, P, 1);
        }
        //Advection contribution
        m_mapping->IncNSAdvectionCorrection(vel, Forcing);

        // Pressure contribution
        m_mapping->IncNSPressureCorrection(P, tmp);
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
        }   
        // Viscous contribution
        m_mapping->IncNSViscousCorrection(vel, tmp);
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {            
            Vmath::Smul(physTot, m_kinvis, tmp[i], 1, tmp[i], 1);
            Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
        }

        // If necessary, transform to wavespace
        if(m_fields[0]->GetWaveSpace())
        {
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[0]->HomogeneousFwdTrans(Forcing[i],Forcing[i]);
            }
        }
        
        // Add to outarray 
        for (int i = 0; i < m_nConvectiveFields; ++i)
        {
            Vmath::Vadd(physTot, outarray[i], 1, Forcing[i], 1, outarray[i], 1);
        }          
    }
    
} //end of namespace
