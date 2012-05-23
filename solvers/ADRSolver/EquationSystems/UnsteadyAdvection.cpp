///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvection.cpp
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
// Description: Unsteady advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ADRSolver/EquationSystems/UnsteadyAdvection.h>

namespace Nektar
{
    string UnsteadyAdvection::className = GetEquationSystemFactory().RegisterCreatorFunction("UnsteadyAdvection", UnsteadyAdvection::create, "Unsteady Advection equation.");

    UnsteadyAdvection::UnsteadyAdvection(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    void UnsteadyAdvection::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        /// Define a global variable to store the advection velocities
        m_velocity = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
        
        /// Read the advection velocities from session file 
        std::vector<std::string> vel;
        vel.push_back("Vx");
        vel.push_back("Vy");
        vel.push_back("Vz");
        
        /// Resize the advection velocitites vector to the dimension of the problem
        vel.resize(m_spacedim);

        /// Store in the global variable m_velocity the advection velocities
        EvaluateFunction(vel, m_velocity, "AdvectionVelocity");

        /// If explicit it computes the RHS and the PROJECTION for the time integration
        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&UnsteadyAdvection::DoOdeRhs,        this);
            m_ode.DefineProjection (&UnsteadyAdvection::DoOdeProjection, this);
        }
        /// Otherwise it gives an error because there is no implicit integration at the moment
        else
        {
            ASSERTL0(false, "Implicit unsteady Advection not set up.");
        }
    }

    UnsteadyAdvection::~UnsteadyAdvection()
    {

    }

    void UnsteadyAdvection::DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                                           Array<OneD,        Array<OneD, NekDouble> >&outarray,
                                     const NekDouble time)
    {
        /// Counter variables
        int i, j;
        
        /// Always equal to one because there is only the field u for the linear advection equation 
        int nVariables      = inarray.num_elements();
        
        /// Number of quadrature points
        int nQuadraturePts  = GetNpoints();
        
        /// Number of elements
        int nElements       = m_fields[0]->GetExpSize();

        /// Switch on the projection type (Discontinuous or Continuous)
        switch (m_projectionType)
        {
            /// Discontinuous approach
            case MultiRegions::eDiscontinuousGalerkin:
            {
                /// Discontinuous Galerkin approach standard
                if(m_discontinuousApproach == "StandardDG")
                {
                    /// Get the number of coefficients
                    int ncoeffs    = inarray[0].num_elements();
                    
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
                /// Calculate - ( V \cdot Grad(u) );
                for(i = 0; i < nVariables; ++i)
                {
                    AdvectionNonConservativeForm(m_velocity, inarray[i], outarray[i]);
                    Vmath::Neg(nQuadraturePts, outarray[i], 1);
                }
                break;
            }
        }
    }



    /**
     *
     */
    void UnsteadyAdvection::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                  Array<OneD,       Array<OneD, NekDouble> >&outarray,
                                            const NekDouble time)
    {
        int i;
        int nVariables = inarray.num_elements();
        SetBoundaryConditions(time);

        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuousGalerkin:
            {
                // Just copy over array
                int nQuadraturePts = GetNpoints();

                for(i = 0; i < nVariables; ++i)
                {
                    Vmath::Vcopy(nQuadraturePts,inarray[i],1,outarray[i],1);
                }
            }
            break;
        case MultiRegions::eGalerkin:
            {
                Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

                for(i = 0; i < nVariables; ++i)
                {
                    m_fields[i]->FwdTrans(inarray[i],coeffs,false);
                    m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
                }
                break;
            }
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }


    void UnsteadyAdvection::v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield,
                           Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(),physfield[i],1,
                m_velocity[j],1,flux[j],1);
        }
    }

    void UnsteadyAdvection::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        int i;

        int nTraceNumPoints = GetTraceNpoints();
        int nvel = m_spacedim; //m_velocity.num_elements();

        Array<OneD, NekDouble > Fwd(nTraceNumPoints);
        Array<OneD, NekDouble > Bwd(nTraceNumPoints);
        Array<OneD, NekDouble > Vn (nTraceNumPoints,0.0);		
        
        //Get Edge Velocity - Could be stored if time independent
        for(i = 0; i < nvel; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_velocity[i], Fwd);
            Vmath::Vvtvp(nTraceNumPoints,m_traceNormals[i],1,Fwd,1,Vn,1,Vn,1);
        }
        
        for(i = 0; i < numflux.num_elements(); ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(physfield[i],Fwd,Bwd);
            
            //evaulate upwinded m_fields[i]
            if (m_expdim == 1)
            {
                m_fields[i]->GetTrace1D()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            }
            else if (m_expdim == 2)
            {
                m_fields[i]->GetTrace()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            }
            else if (m_expdim == 3)
            {
                m_fields[i]->GetTrace3D()->Upwind(Vn,Fwd,Bwd,numflux[i]);
            }
            // calculate m_fields[i]*Vn
            Vmath::Vmul(nTraceNumPoints,numflux[i],1,Vn,1,numflux[i],1);
        }
    }
    
    void UnsteadyAdvection::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
    }
}
