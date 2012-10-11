///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <boost/math/constants/constants.hpp>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

namespace Nektar
{
    string CompressibleFlowSystem::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "CompressibleFlowSystem", 
            CompressibleFlowSystem::create, 
            "Auxiliary functions for the compressible flow system.");
  
    CompressibleFlowSystem::CompressibleFlowSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }
    
    void CompressibleFlowSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        
	ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
		 "No UPWINDTYPE defined in session.");

        // Set up locations of velocity vector.
        m_velLoc = Array<OneD, NekDouble>(m_expdim);
        for (int i = 0; i < m_expdim; ++i)
        {
            m_velLoc[i] = i+1;
        }

        // Create Riemann solver instance, depending on the UPWINDTYPE specified
        // in the session file. Bind gamma, velLoc and the trace normals.
        m_riemannSolver = SolverUtils::GetRiemannSolverFactory().
            CreateInstance(m_session->GetSolverInfo("UPWINDTYPE"));
        m_riemannSolver->AddParam ("gamma",  
                                   &CompressibleFlowSystem::GetGamma,   this);
        m_riemannSolver->AddScalar("velLoc", 
                                   &CompressibleFlowSystem::GetVelLoc,  this);
        m_riemannSolver->AddVector("N",
                                   &CompressibleFlowSystem::GetNormals, this);
        
        // Create an advection object. For now, weak DG is hard-coded but
        // eventually this choice will be defined by the user, with a default
        // being set in UnsteadySystem. Bind flux vector and the Riemann solver.
        m_advection = SolverUtils::GetAdvectionFactory().
            CreateInstance("WeakDG", "WeakDG");
        m_advection->SetFluxVector(
            &CompressibleFlowSystem::GetFluxVector, this);
        m_advection->SetRiemannSolver(m_riemannSolver);
    }
    
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {
        
    }
    
    void CompressibleFlowSystem::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
    
    //----------------------------------------------------
    
    void CompressibleFlowSystem::WallBoundary(
        int                                   b,
        int                                   cnt, 
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    { 
        int i;
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }
        
        // Adjust the physical values of the trace to take 
        // user defined boundaries into account
        int e, id1, id2, npts;
        
        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[b]->GetExpSize();++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[b]->
                GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[b]->
                GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->
                GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            switch(m_expdim)
            {
                // Special case for 2D.
                case 2:
                {
                    Array<OneD, NekDouble> tmp_n(npts);
                    Array<OneD, NekDouble> tmp_t(npts);
                    
                    Vmath::Vmul (npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,
                                 &tmp_n[0],1);
                    Vmath::Vvtvp(npts,&Fwd[2][id2],1,&m_traceNormals[1][id2],1,
                                 &tmp_n[0],1,&tmp_n[0],1);
	      
                    Vmath::Vmul (npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,
                                 &tmp_t[0],1);
                    Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,
                                 &tmp_t[0],1,&tmp_t[0],1);
                    
                    // negate the normal flux
                    Vmath::Neg  (npts,tmp_n,1);		      
                    
                    // rotate back to Cartesian
                    Vmath::Vmul (npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,
                                 &Fwd[1][id2],1);
                    Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,
                                 &Fwd[1][id2],1,&Fwd[1][id2],1);
	      
                    Vmath::Vmul (npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,
                                 &Fwd[2][id2],1);
                    Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,
                                 &Fwd[2][id2],1,&Fwd[2][id2],1);
                    break;
                }
                
                // For 1D/3D, define: v* = v - (v.n)n so that v*.n = 0
                case 1:
                case 3:
                {
                    Array<OneD,NekDouble> tmp(npts);
                    
                    // Calculate (v.n)
                    for (i = 0; i < m_expdim; ++i)
                    {
                        Vmath::Vvtvp(npts,
                                     &Fwd[1+i][id2],1,&m_traceNormals[i][id2],1,
                                     &tmp[0],1,&tmp[0],1);
                    }
                    
                    for (i = 0; i < m_expdim; ++i)
                    {
                        Vmath::Vvtvm(npts,&tmp[0],1,&m_traceNormals[i][id2],1,
                                     &Fwd[1+i][id2],1,&Fwd[1+i][id2],1);
                        Vmath::Neg  (npts,&Fwd[1+i][id2],1);
                    }

                    break;
                }
                default:
                    ASSERTL0(false,"Illegal expansion dimension");
            }
            
            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npts,&Fwd[i][id2],1,&(m_fields[i]->
                    GetBndCondExpansions()[b]->UpdatePhys())[id1],1);
            }
        }
    }
  
    void CompressibleFlowSystem::SymmetryBoundary(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
    {  
        int i;
        int nTraceNumPoints = GetTraceTotPoints();
        int nvariables      = physarray.num_elements();
        
        // get physical values of the forward trace (from exp to phys)
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }
        
        int e, id1, id2, npts;
        
        for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            
            switch(m_expdim)
            {
                case 1:
                {
                    ASSERTL0(false,"1D not yet implemented for the Compressible Flow Equations");
                }
                break;
                case 2:
                {
                    Array<OneD, NekDouble> tmp_t(npts);
                    
                    Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
                    Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);
                    
                    Array<OneD, NekDouble> tmp_n(npts,0.0);
                    
                    // rotate back to Cartesian
                    Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
                    Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);
                    
                    Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
                    Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
                }
                break;
                case 3:
                    ASSERTL0(false,"3D not implemented for the Compressible Flow Equations");
                    break;
                default:
                    ASSERTL0(false,"Illegal expansion dimension");
            }
            
            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npts,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);	
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     * 
     * @param i           Component of the flux vector to calculate.
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const int                                   i,
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        int j, nq = m_fields[0]->GetTotPoints();
        
        if (i == 0)
        {
            // Flux vector for the rho equation.
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vcopy(nq, physfield[j+1], 1, flux[j], 1);
            }
        } 
        else if (i >= 1 && i <= m_expdim)
        {
            // Flux vector for the velocity fields.
            Array<OneD, NekDouble> pressure(nq);
            GetPressure(physfield,pressure);
            
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vmul(nq, physfield[j+1], 1, physfield[i], 1, flux[j], 1);
                Vmath::Vdiv(nq, flux[j], 1, physfield[0], 1, flux[j], 1);
            }
            
            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i-1], 1, pressure, 1, flux[i-1], 1);
        }
        else if (i == m_expdim+1) 
        {
            // Flux vector for the total energy field.
            Array<OneD, NekDouble> pressure(nq);
            GetPressure(physfield, pressure);
            Vmath::Vadd(nq, physfield[m_expdim+1], 1, pressure, 1, pressure, 1);
            
            for (j = 0; j < m_expdim; ++j)
            {
                Vmath::Vdiv(nq, physfield[j+1], 1, physfield[0], 1, flux[j], 1);
                Vmath::Vmul(nq, flux[j], 1, pressure, 1, flux[j], 1);
            }
        }
        else
        {
            ASSERTL0(false, "Invalid vector index.");
        }
    }
    
    /**
     * @brief Calculate the pressure field \f$ p =
     * (\gamma-1)(E-\frac{1}{2}\rho\| \mathbf{v} \|^2) \f$ assuming an ideal gas
     * law.
     * 
     * @param physfield  Input momentum.
     * @param pressure   Computed pressure field.
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD,                   NekDouble>   &pressure)
    {
        NekDouble gamma = m_gamma;
        int       npts  = m_fields[0]->GetTotPoints();
        double    alpha = -0.5;
        
        // Calculate ||rho v||^2.
        Vmath::Zero(npts, pressure, 1);
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vvtvp(npts, physfield[1+i], 1, physfield[1+i], 1, 
                               pressure,       1, pressure,       1);
        }
        // Divide by rho to get rho*||v||^2.
        Vmath::Vdiv (npts, pressure, 1, physfield[0], 1, pressure, 1);
        // pressure <- E - 0.5*pressure
        Vmath::Svtvp(npts,     alpha, 
                     pressure, 1, physfield[m_expdim+1], 1, pressure, 1);
        // Multiply by (gamma-1).
        Vmath::Smul (npts, m_gamma-1, pressure, 1, pressure, 1);
    }
    
    /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ \rho\mathbf{v} \f$.
     * 
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void CompressibleFlowSystem::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int npts = m_fields[0]->GetTotPoints();
        
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vdiv(npts, physfield[1+i], 1, physfield[0], 1, 
                              velocity[i],    1);
        }
    }
  
    /**
     * @brief Compute the temperature \f$ T = p/\rho R \f$.
     * 
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param temperature  The resulting temperature \f$ T \f$.
     */
    void CompressibleFlowSystem::GetTemperature(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &pressure,
        Array<OneD,             NekDouble  > &temperature)
    {
        for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
        {
            temperature[i] = pressure[i]/(physfield[0][i]*m_GasConstant);
        }
    }
    
    /**
     * @brief Compute the sound speed \f$ c = sqrt(\gamma p/\rho) \f$.
     * 
     * @param physfield    Input physical field.
     * @param pressure     The pressure field corresponding to physfield.
     * @param soundspeed   The resulting sound speed \f$ c \f$.
     */
    void CompressibleFlowSystem::GetSoundSpeed(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD,             NekDouble  > &pressure,
              Array<OneD,             NekDouble  > &soundspeed)
    {
        for (int i = 0; i < m_fields[0]->GetTotPoints(); ++i)
        {
            soundspeed[i] = sqrt(m_gamma*pressure[i]/physfield[0][i]);
        }
    }
    
    /**
     * @brief Compute the mach number \f$ M = \| \mathbf{v} \|^2 / c \f$.
     * 
     * @param physfield    Input physical field.
     * @param soundfield   The speed of sound corresponding to physfield.
     * @param mach         The resulting mach number \f$ M \f$.
     */
    void CompressibleFlowSystem::GetMach(
        Array<OneD, Array<OneD, NekDouble> > &physfield,
        Array<OneD,             NekDouble  > &soundspeed,
        Array<OneD,             NekDouble  > &mach)
    {
        const int npts = m_fields[0]->GetTotPoints();

        Vmath::Zero(npts, mach, 1);
        
        for (int i = 0; i < m_expdim; ++i)
        {
            Vmath::Vvtvp(npts, physfield[1+i], 1, physfield[1+i], 1, 
                               mach,           1, mach,           1);
        }
        Vmath::Vdiv(npts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(npts, mach, 1, physfield[0], 1, mach, 1);
        Vmath::Vdiv(npts, mach, 1, soundspeed,   1, mach, 1);
    }

    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     * 
     * @param physarray  Physical field.
     * @param ExpOrder   Polynomial order of each element in the mesh.
     * @param CFL        CFL number for each element in the mesh.
     * @param timeCFL    ?
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
        const Array<OneD, Array<OneD,NekDouble> > physarray, 
        const Array<OneD, int>                    ExpOrder, 
        const Array<OneD, NekDouble>              CFL,
        NekDouble                                 timeCFL)
    { 
        int nvariables     = m_fields.num_elements();
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize(); 
        
        Array<OneD, NekDouble> tstep      (n_element,0.0);
        Array<OneD, NekDouble> stdVelocity(n_element);
        GetStdVelocity(physarray, stdVelocity);
        
        // TODO: This should be implemented as a virtual function inside
        // StdExpansion.
        map<StdRegions::ExpansionType,double> minLengths;
        
        minLengths[StdRegions::eTriangle     ] = 0.5 / sqrt(2.0);
        minLengths[StdRegions::eQuadrilateral] = 0.5;
        minLengths[StdRegions::eHexahedron   ] = 0.25;
        
        for(int el = 0; el < n_element; ++el)
        {
            //tstep[el] = CFL[el]/stdVelocity[el];
            tstep[el] = CFL[el]*minLengths[
                m_fields[0]->GetExp(el)->DetExpansionType()]/stdVelocity[el];
        }
        
        double minDt = Vmath::Vmin(n_element, tstep, 1);
        m_comm->AllReduce(minDt, LibUtilities::ReduceMin);
        return minDt;
    }
    
    void CompressibleFlowSystem::GetStdVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,                   NekDouble>   &stdV)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int npts           = 0;

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_expdim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_expdim);
        Array<OneD, NekDouble>               pressure   (nTotQuadPoints);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        
        // Zero output array.
        Vmath::Zero(stdV.num_elements(), stdV, 1);
        
        for (int i = 0; i < m_expdim; ++i)
        {
            velocity   [i] = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }
        GetVelocityVector(inarray,velocity);
        GetPressure      (inarray,pressure);
        GetSoundSpeed    (inarray,pressure,soundspeed);
        
        for(int el = 0; el < n_element; ++el)
        { 
            const Array<OneD, const NekDouble> &jac  = 
                m_fields[0]->GetExp(el)->GetGeom()->GetJac();
            const Array<TwoD, const NekDouble> &gmat = 
                m_fields[0]->GetExp(el)->GetGeom()->GetGmat();
            
            int nqtot = m_fields[0]->GetExp(el)->GetTotPoints();
            
            if(m_fields[0]->GetExp(el)->GetGeom()->GetGtype() == 
                   SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < m_expdim; ++i)
                {
                    Vmath::Zero(nqtot, stdVelocity[i], 1);
                    for (int j = 0; j < m_expdim; ++j)
                    {
                        Vmath::Vvtvp(nqtot, gmat[m_expdim*j+i], 1, velocity[j], 
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < m_expdim; ++i)
                {
                    Vmath::Zero(nqtot, stdVelocity[i], 1);
                    for (int j = 0; j < m_expdim; ++j)
                    {
                        Vmath::Svtvp(nqtot, gmat[m_expdim*j+i][0], velocity[j], 
                                     1, stdVelocity[i], 1, stdVelocity[i], 1);
                    }
                }
            }
            
            for(int i = 0; i < nqtot; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < m_expdim; ++j)
                {
                    pntVelocity += stdVelocity[j][i]*stdVelocity[j][i];
                }
                pntVelocity = sqrt(pntVelocity) + soundspeed[npts];
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
                npts++;
            }
        }
    }
}
