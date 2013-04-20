///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.h
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

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/Diffusion/Diffusion.h>


#define EPSILON 0.000001

#define CROSS(dest, v1, v2){                 \
          dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
          dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
          dest[2] = v1[0] * v2[1] - v1[1] * v2[0];}

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2){       \
          dest[0] = v1[0] - v2[0]; \
          dest[1] = v1[1] - v2[1]; \
          dest[2] = v1[2] - v2[2];}

namespace Nektar
{  
    /**
     * 
     */
    class CompressibleFlowSystem: public SolverUtils::UnsteadySystem
    {
    public:

        friend class MemoryManager<CompressibleFlowSystem>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            return MemoryManager<CompressibleFlowSystem>::AllocateSharedPtr(pSession);
        }
        /// Name of class
        static std::string className;
      
        virtual ~CompressibleFlowSystem();
        
        /// Function to calculate the stability limit for DG/CG.
        NekDouble GetStabilityLimit(int n);
        
        /// Function to calculate the stability limit for DG/CG 
        /// (a vector of them).
        Array<OneD, NekDouble> GetStabilityLimitVector(
            const Array<OneD,int> &ExpOrder);
      
    protected:
        SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
        SolverUtils::RiemannSolverSharedPtr m_riemannSolverLDG;
        SolverUtils::AdvectionSharedPtr     m_advection;
        SolverUtils::DiffusionSharedPtr     m_diffusion;
        Array<OneD, NekDouble>              m_velLoc;
        NekDouble                           m_gamma;
        NekDouble                           m_pInf;
        NekDouble                           m_rhoInf;
        NekDouble                           m_uInf;
        NekDouble                           m_vInf;
        NekDouble                           m_wInf;
        NekDouble                           m_gasConstant;
        NekDouble                           m_Twall;
        std::string                         m_ViscosityType;
        NekDouble                           m_mu;
        NekDouble                           m_thermalConductivity;
        NekDouble                           m_Cp;
        NekDouble                           m_Prandtl;


        CompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        virtual void v_InitObject();
      
        /// Print a summary of time stepping parameters.
        virtual void v_PrintSummary(std::ostream &out);

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        void GetViscousFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);
        void WallBoundary(
            int                                               bcRegion,
            int                                               cnt,
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void WallBoundaryViscous(
            int                                               bcRegion,
            int                                               cnt,
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void SymmetryBoundary(
            int                                               bcRegion,
            int                                               cnt,
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void InflowCFSBoundary(
            int                                               bcRegion, 
            int                                               cnt, 
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void OutflowCFSBoundary(
            int                                               bcRegion, 
            int                                               cnt, 
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void ExtrapOrder0Boundary(
            int                                               bcRegion, 
            int                                               cnt, 
            Array<OneD, Array<OneD, NekDouble> >             &physarray);
        void GetVelocityVector(
            const Array<OneD,       Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,       Array<OneD,       NekDouble> >&velocity);
        void GetSoundSpeed(
            const Array<OneD,       Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,                         NekDouble>  &pressure,
                  Array<OneD,                         NekDouble>  &soundspeed);
        void GetMach(
                  Array<OneD,       Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,                         NekDouble>  &soundspeed,
                  Array<OneD,                         NekDouble>  &mach);
        void GetTemperature(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,                         NekDouble>  &pressure,
                  Array<OneD,                         NekDouble>  &temperature);
        void GetPressure(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,                         NekDouble>  &pressure);
        void GetPressure(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
            const Array<OneD, const Array<OneD,       NekDouble> >&velocity,
                  Array<OneD,                         NekDouble>  &pressure);
        void GetDynamicViscosity(
            const Array<OneD, const Array<OneD,       NekDouble> >&physfield,
                  Array<OneD,                         NekDouble  >&mu);
        void GetStdVelocity(
            const Array<OneD, const Array<OneD,       NekDouble> >&inarray,
                  Array<OneD,                         NekDouble>  &stdV);
      
        virtual NekDouble v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray);

        virtual void v_SetInitialConditions(
            NekDouble initialtime = 0.0,
            bool dumpInitialConditions = true)
        {
        }

        NekDouble GetGasConstant()
        {
            return m_gasConstant;
        }
        
        NekDouble GetGamma()
        {
            return m_gamma;
        }
      
        const Array<OneD, NekDouble> &GetVelLoc()
        {
            return m_velLoc;
        }
      
        const Array<OneD, const Array<OneD, NekDouble> > &GetNormals()
        {
            return m_traceNormals;
        }
      
    private:

    };
}
#endif
