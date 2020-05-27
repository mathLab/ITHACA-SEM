///////////////////////////////////////////////////////////////////////////////
//
// File VCSMapping.h
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
// Description: Velocity Correction Scheme with coordinate transformation header 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VCSMAPPING_H
#define NEKTAR_SOLVERS_VCSMAPPING_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <GlobalMapping/Mapping.h>

namespace Nektar
{
    class VCSMapping: public VelocityCorrectionScheme
    {
    public:

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<VCSMapping>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;


        /// Constructor.
        VCSMapping(const LibUtilities::SessionReaderSharedPtr& pSession,
                   const SpatialDomains::MeshGraphSharedPtr &pGraph);

        // 
        void ApplyIncNSMappingForcing (
                const Array<OneD, const Array<OneD, NekDouble> >  &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray);

        virtual ~VCSMapping();

        virtual void v_InitObject();

    protected:
        // Mapping object
        GlobalMapping::MappingSharedPtr             m_mapping;
        
        bool                                        m_verbose;
        
        // Flags defining how pressure and viscous mapping terms 
        //should be treated
        bool                                        m_implicitPressure;
        bool                                        m_implicitViscous;
        bool                                        m_neglectViscous;
        // Tolerance and relaxation parameters for pressure and viscous
        //       systems (when solved iteratively)
        NekDouble                                   m_pressureTolerance;
        NekDouble                                   m_viscousTolerance;
        NekDouble                                   m_pressureRelaxation;
        NekDouble                                   m_viscousRelaxation;
        
        // Pressure gradient (to avoid duplicate calculations)
        Array<OneD, Array<OneD, NekDouble> >        m_gradP;

        // Virtual functions     
        virtual void v_DoInitialise(void);
        
        
        virtual void v_SetUpPressureForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        
        virtual void v_SetUpViscousForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt);
        
        virtual void v_SolvePressure( const Array<OneD, NekDouble>  &Forcing);
        
        virtual void v_SolveViscous( 
                    const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble aii_Dt);
        
        virtual void v_EvaluateAdvection_SetPressureBCs(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time);
    
    private:        
        Array<OneD, Array<OneD, NekDouble> >    m_presForcingCorrection;
        
        // Correction needed for convective terms = N(u) - ( -(u \nabla) u)
        //     velPhys is the velocity field (transformed for physical space)
        void MappingAdvectionCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        // Correction needed for time-derivative terms 
        //     = U_coord^j u^i_,j - u^j U_coord^i_,j
        //     vel     is the velocity field (can be in wavespace)
        //     velPhys is the velocity field (transformed for physical space)
        void MappingAccelerationCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &vel,
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        // Correction needed for pressure terms   
        //     = -g^(ij)p_j + (\nabla p)/J for variable Jacobian
        //     = -g^(ij)p_j + (\nabla p)   for constant Jacobian
        //         the pressure field can be in wavespace
        void MappingPressureCorrection(
            Array<OneD, Array<OneD, NekDouble> >              &outarray);

        // Correction needed for viscous terms = g^jk u^i_{,jk}-(\nabla^2 u)
        //     vel     is the velocity field (can be in wavespace)
        //     velPhys is the velocity field (transformed for physical space)
        void MappingViscousCorrection(
            const Array<OneD, const Array<OneD, NekDouble> >  &velPhys,
            Array<OneD, Array<OneD, NekDouble> >              &outarray);
        
    };

    typedef std::shared_ptr<VCSMapping>
                VCSMappingSharedPtr;

} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H
