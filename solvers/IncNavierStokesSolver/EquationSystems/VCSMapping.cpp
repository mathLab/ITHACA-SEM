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

       // Storage to extrapolate pressure forcing
        int physTot = m_fields[0]->GetTotPoints();
        int intSteps;
        int intMethod = m_intScheme->GetIntegrationMethod();
        switch(intMethod)
        {
            case LibUtilities::eIMEXOrder1:
            {
                intSteps = 1; 
            }
            break;
            case LibUtilities::eIMEXOrder2:
            {
                intSteps = 2;
            }
            break;
            case LibUtilities::eIMEXOrder3:
            {
                intSteps = 3;
            }
            break;
        }        
        m_presForcingCorrection = Array<OneD, Array<OneD, NekDouble> > (intSteps);
        for(int i = 0; i < m_presForcingCorrection.num_elements(); i++)
        {
            m_presForcingCorrection[i] = Array<OneD, NekDouble>(physTot,0.0);
        }
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
        if (m_mapping->HasConstantJacobian())
        {
            VelocityCorrectionScheme::v_SetUpPressureForcing(fields, Forcing, aii_Dt);
        }
        else
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_nConvectiveFields;
            Array<OneD, NekDouble> wk(physTot, 0.0);

            Array<OneD, NekDouble> Jac(physTot,0.0);
            m_mapping->GetJacobian(Jac);
            
            // Calculate div(-J*u/Dt)
            Vmath::Zero(physTot,Forcing[0],1);
            for(int i = 0; i < nvel; ++i)
            {
                if (m_fields[i]->GetWaveSpace())
                {
                    m_fields[i]->HomogeneousBwdTrans(fields[i],wk);
                }
                else
                {
                    Vmath::Vcopy(physTot, fields[i], 1, wk, 1);
                }
                Vmath::Vmul(physTot,wk,1,Jac,1,wk,1);
                if (m_fields[i]->GetWaveSpace())
                {
                    m_fields[i]->HomogeneousFwdTrans(wk,wk);
                }               
                m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[i],wk, wk);
                Vmath::Vadd(physTot,wk,1,Forcing[0],1,Forcing[0],1);
            }
            Vmath::Smul(physTot,1.0/aii_Dt,Forcing[0],1,Forcing[0],1); 

            //
            // If the mapping viscous terms are being treated explicitly 
            //        we need to apply a correction to the forcing
            if (!m_mapping->ImplicitViscous())
            {               
                bool wavespace = m_fields[0]->GetWaveSpace();
                m_fields[0]->SetWaveSpace(false);

                // 
                //  Part 1: div(J*grad(U/J . grad(J)))
                Array<OneD, Array<OneD, NekDouble> > tmp (nvel);
                Array<OneD, Array<OneD, NekDouble> > velocity (nvel);
                for(int i = 0; i < tmp.num_elements(); i++)
                {
                    tmp[i] = Array<OneD, NekDouble>(physTot,0.0);
                    velocity[i] = Array<OneD, NekDouble>(physTot,0.0);
                    if (wavespace)
                    {
                        m_fields[0]->HomogeneousBwdTrans(m_fields[i]->GetPhys(),
                                                         velocity[i]);
                    }
                    else
                    {
                        Vmath::Vcopy(physTot, m_fields[i]->GetPhys(), 1, 
                                                velocity[i], 1);
                    }
                }
                // Calculate wk = U.grad(J)
                m_mapping->DotGradJacobian(velocity, wk);
                // Calculate wk = (U.grad(J))/J
                Vmath::Vdiv(physTot, wk, 1, Jac, 1, wk, 1);
                // J*grad[(U.grad(J))/J]
                for(int i = 0; i < nvel; ++i)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], 
                                wk, tmp[i]);
                    Vmath::Vmul(physTot, Jac, 1, tmp[i], 1, tmp[i], 1);
                }          
                // div(J*grad[(U.grad(J))/J])
                Vmath::Zero(physTot, wk, 1);
                for(int i = 0; i < nvel; ++i)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], 
                                tmp[i], tmp[i]);
                    Vmath::Vadd(physTot, wk, 1, tmp[i], 1, wk, 1);
                }        

                // Part 2: grad(J) . curl(curl(U))
                m_mapping->CurlCurlField(velocity, tmp);
                m_mapping->DotGradJacobian(tmp, velocity[0]); // dont need velocity any more

                // Add two parts
                Vmath::Vadd(physTot, velocity[0], 1, wk, 1, wk, 1);

                // Multiply by kinvis
                Vmath::Smul(physTot, m_kinvis, wk, 1, wk, 1);

                // Extrapolate correction
                m_extrapolation->ExtrapolateArray(m_presForcingCorrection, wk, wk);

                // Put in wavespace
                if (wavespace)
                {
                    m_fields[0]->HomogeneousFwdTrans(wk,wk);
                }
                // Apply correction: Forcing = Forcing - correction
                Vmath::Vsub(physTot, Forcing[0], 1, wk, 1, Forcing[0], 1);

                m_fields[0]->SetWaveSpace(wavespace); 
            }
        }
    }
    
    /**
     * Forcing term for Helmholtz solver
     */
    void   VCSMapping::v_SetUpViscousForcing(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        Array<OneD, Array<OneD, NekDouble> > &Forcing, 
        const NekDouble aii_Dt)
    {
        if (m_mapping->HasConstantJacobian())
        {
            VelocityCorrectionScheme::v_SetUpViscousForcing(inarray, Forcing, aii_Dt);
        }
        else
        {
            NekDouble aii_dtinv = 1.0/aii_Dt;
            int physTot = m_fields[0]->GetTotPoints();

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
            
            // If the pressure terms are treated explicitly, we need to divide by J
            //    if they are implicit, we need to calculate G(p)
            Array<OneD, Array<OneD, NekDouble> > tmp (nvel);
            if (m_pressure->GetWaveSpace())
            {
                
                for (int i=0; i<nvel; i++)
                {                    
                    tmp[i] = Array<OneD, NekDouble>(physTot,0.0);
                    m_pressure->HomogeneousBwdTrans(Forcing[i],tmp[i]);
                }                
            }
            else
            {
                for (int i=0; i<nvel; i++)
                {                    
                    tmp[i] = Array<OneD, NekDouble>(physTot,0.0);
                    Vmath::Vcopy(physTot, Forcing[i], 1, tmp[i], 1);
                }
                
            }
            if (m_mapping->ImplicitPressure())
            {                
                m_mapping->RaiseIndex(tmp, Forcing);
            }
            else
            {
                Array<OneD, NekDouble> Jac(physTot,0.0);
                m_mapping->GetJacobian(Jac);
                for (int i=0; i<nvel; i++)
                {
                    Vmath::Vdiv(physTot, tmp[i], 1, Jac, 1, Forcing[i], 1);
                }               
            }
            // Transform back to wavespace
            if (m_pressure->GetWaveSpace())
            {
                for (int i=0; i<nvel; i++)
                {
                    m_pressure->HomogeneousFwdTrans(Forcing[i],Forcing[i]);
                }                
            }          
            
            // Subtract inarray/(aii_dt) and divide by kinvis. Kinvis will
            // need to be updated for the convected fields.
            for(int i = 0; i < nvel; ++i)
            {
                Blas::Daxpy(physTot,-aii_dtinv,inarray[i],1,Forcing[i],1);
                Blas::Dscal(physTot,1.0/m_kinvis,&(Forcing[i])[0],1);
            }            
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
        if (!m_mapping->ImplicitPressure())
        {
            m_mapping->IncNSPressureCorrection(P, tmp);
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
            }            
        }   
        // Viscous contribution
        if (!m_mapping->ImplicitViscous())
        {
            m_mapping->IncNSViscousCorrection(vel, tmp);
            for (int i = 0; i < m_nConvectiveFields; ++i)
            {            
                Vmath::Smul(physTot, m_kinvis, tmp[i], 1, tmp[i], 1);
                Vmath::Vadd(physTot, tmp[i], 1, Forcing[i], 1, Forcing[i], 1);
            }
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
