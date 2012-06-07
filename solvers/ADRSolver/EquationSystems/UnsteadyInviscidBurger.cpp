///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyInviscidBurger.cpp
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
// Description: Unsteady inviscid Burger solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyInviscidBurger.h>

namespace Nektar
{
    string UnsteadyInviscidBurger::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyInviscidBurger", UnsteadyInviscidBurger::create, "Inviscid Burger equation");
    
    UnsteadyInviscidBurger::UnsteadyInviscidBurger(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }
    
    void UnsteadyInviscidBurger::v_InitObject()
    {
        UnsteadySystem::v_InitObject();
        
        // Useless parameter
        //m_session->LoadParameter("wavefreq",   m_waveFreq, 0.0);
        
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyInviscidBurger::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyInviscidBurger::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }
    
    UnsteadyInviscidBurger::~UnsteadyInviscidBurger()
    {
        
    }
    
    
    
    void UnsteadyInviscidBurger::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                                Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                          const NekDouble time)
    {
        /// Counter variables
        int i, j;
        
        /// Always equal to one because there is only the field u
        int nVariables = inarray.num_elements();
        
        /// Number of quadrature points
        int nQuadraturePts  = GetNpoints();
        
        /// Number of elements
        int nElements       = m_fields[0]->GetExpSize();
        
        /// Switch on the projection type (Discontinuous or Continuous)
        switch (m_projectionType)
        {
            case MultiRegions::eDiscontinuousGalerkin:
            {
                /// Discontinuous Galerkin approach standard
                if(m_discontinuousApproach == "StandardDG")
                {
                    /// Get the number of coefficients
                    int ncoeffs    = GetNcoeffs();
                    
                    /// Define an auxiliary variable to compute the RHS 
                    Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
                    WeakAdv[0] = Array<OneD, NekDouble>(ncoeffs*nVariables);
                    for(i = 1; i < nVariables; ++i)
                    {
                        WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
                    }
                    
                    /// Call the method to compute the weak flux
                    WeakDGAdvection(inarray, WeakAdv, true, true);
                    
                    /// Operations to compute the RHS
                    for(i = 0; i < nVariables; ++i)
                    {
                        /// Multiply the flux by the inverse of the mass matrix
                        m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
                        
                        /// Store in outarray the physical values of the RHS
                        m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
                        
                        /// Negate the RHS
                        Vmath::Neg(nQuadraturePts, outarray[i], 1);
                    }
                    
                }
                
                /// Flux reconstruction approach
                else if(m_discontinuousApproach == "FR-DG" || m_discontinuousApproach == "FR-SD" || m_discontinuousApproach == "FR-HU")
                {        
                    
                    /// Call the method to compute the strong divergence of the flux
                    StrongFRAdvection(inarray, outarray, true);
                    
                    /// Array to store the Jacobian and its inverse
                    Array<OneD, const NekDouble>jac(nElements);
                    Array<OneD, NekDouble>      jacobian(nElements);
                    Array<OneD, NekDouble>      tmparray;
                    
                    
                    /// Evaluation of the jacobian of each element
                    for(i = 0; i < nElements; i++)
                    {
                        jac         = m_fields[0]->GetExp(i)->GetGeom1D()->GetJac();
                        jacobian[i] = jac[0];
                    }
                    
                    /// Operations to compute the RHS
                    for(i = 0; i < nVariables; ++i)
                    {
                        for(j = 0; j < nElements; j++)
                        {
                            Vmath::Smul(nQuadraturePts/nElements, 1/jacobian[j], 
                                        tmparray = outarray[i] + j*nQuadraturePts/nElements, 1.0, 
                                        tmparray = outarray[i] + j*nQuadraturePts/nElements, 1.0);
                        }
                        
                        Vmath::Neg(nQuadraturePts, outarray[i], 1);
                    }
                }
                break;
            }
                
            /// Continuous approach
            case MultiRegions::eGalerkin:
            {
                /// Calculate -V \cdot Grad(u);
                for(i = 0; i < nVariables; ++i)
                {
                    AdvectionNonConservativeForm(m_velocity, inarray[i], outarray[i]);
                    Vmath::Neg(nQuadraturePts, outarray[i], 1);
                }
                break;
            }
        }
    }
    
    
    
    void UnsteadyInviscidBurger::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                       Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                                 const NekDouble time)
    {
        /// Counter variable
        int i;
        
        /// Always equal to one because there is only the field u
        int nVariables = inarray.num_elements();
        
        /// Set the boundary conditions (1D periodic discontinuous still doesn't work)
        SetBoundaryConditions(time);
        
        /// Switch on the projection type (Discontinuous or Continuous)
        switch(m_projectionType)
        {
            /// Discontinuous projection
            case MultiRegions::eDiscontinuousGalerkin:
            {
                /// Number of quadrature points
                int nQuadraturePts = GetNpoints();
                
                /// Just copy over array                
                for(i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(nQuadraturePts, inarray[i], 1, outarray[i], 1);
                }
                break;
            }
                
            /// Continuous projection
            case MultiRegions::eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());
                
                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i], coeffs, false);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs, outarray[i]);
                }
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }
    
    
    
    void UnsteadyInviscidBurger::v_GetFluxVector(const int i, 
                                                 Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                 Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(), physfield[i], 1, physfield[i], 1, flux[j], 1);
            Vmath::Smul(GetNpoints(), 0.5, flux[j], 1, flux[j], 1);
        }
    }
    
    
    
    /// Evaulate flux = m_fields * ivel for ith component of Vu for direction j
    void UnsteadyInviscidBurger::v_GetFluxVector(const int i, 
                                                 const int j, 
                                                 Array<OneD, Array<OneD, NekDouble> > &physfield,
                                                 Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL0(false, "should never arrive here ...");
    }
    
    
    
    void UnsteadyInviscidBurger::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                                 Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        /// Counter variable
        int i;
        
        /// Number of trace points
        int nTracePts   = GetTraceNpoints();
        
        /// Number of spatial dimensions
        int nDimensions = m_spacedim;
        
        /// Number of elements
        int nElements = m_fields[0]->GetExpSize();
        
        /// Forward state array
        Array<OneD, NekDouble > Fwd(nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble > Bwd(nTracePts);
        
        /// Normal velocity array
        Array<OneD, NekDouble > Vn (nTracePts, 0.0);
        
        // Extract velocity field along the trace space and multiply by trace normals
        m_fields[0]->ExtractTracePhys(physfield[0], Fwd);
        for(i = 0; i < nDimensions; ++i)
        {
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        } 
        
        /// Compute the numerical fluxes at the trace points
        for(i = 0; i < numflux.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
            
            /// Upwind between elements
            m_fields[i]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux[i]);
            
            /// Calculate the numerical fluxes multipling Fwd or Bwd by the normal advection velocity
            Vmath::Vmul(nTracePts, numflux[i], 1, Vn, 1, numflux[i], 1);
            Vmath::Smul(nTracePts, 0.5, numflux[i], 1, numflux[i], 1);

        }
    }
}

