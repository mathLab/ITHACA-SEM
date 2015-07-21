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
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdHexExp.h>

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
            return MemoryManager<CompressibleFlowSystem>::
                AllocateSharedPtr(pSession);
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
        Array<OneD, Array<OneD, NekDouble> >m_vecLocs;
        NekDouble                           m_gamma;
        NekDouble                           m_pInf;
        NekDouble                           m_rhoInf;
        NekDouble                           m_uInf;
        NekDouble                           m_vInf;
        NekDouble                           m_wInf;
        NekDouble                           m_UInf;
        NekDouble                           m_gasConstant;
        NekDouble                           m_Twall;
        std::string                         m_ViscosityType;
        std::string                         m_shockCaptureType;
        std::string                         m_EqTypeStr;
        NekDouble                           m_mu;
        NekDouble                           m_Skappa;
        NekDouble                           m_Kappa;
        NekDouble                           m_mu0;
        NekDouble                           m_FacL;
        NekDouble                           m_FacH;
        NekDouble                           m_eps_max;
        NekDouble                           m_thermalConductivity;
        NekDouble                           m_Cp;
        NekDouble                           m_C1;
        NekDouble                           m_C2;
        NekDouble                           m_hFactor;
        NekDouble                           m_Prandtl;
        NekDouble                           m_amplitude;
        NekDouble                           m_omega;

        // L2 error file
        std::ofstream m_errFile;

        // Check for steady state at step interval
        int m_steadyStateSteps;

        // Tolerance to which steady state should be evaluated at
        NekDouble m_steadyStateTol;

        // Forcing term
        std::vector<SolverUtils::ForcingSharedPtr> m_forcing;
        StdRegions::StdQuadExpSharedPtr            m_OrthoQuadExp;
        StdRegions::StdHexExpSharedPtr             m_OrthoHexExp;
        bool                                       m_smoothDiffusion;


        // Pressure storage for PressureOutflowFileBC
        Array<OneD, NekDouble> m_pressureStorage;

        // Field storage for PressureInflowFileBC
        Array<OneD, Array<OneD, NekDouble> > m_fieldStorage;

        // Storage for L2 norm error
        Array<OneD, Array<OneD, NekDouble> > m_un;

        CompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        virtual void v_InitObject();

        /// Print a summary of time stepping parameters.
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        void GetViscousFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);
        void GetFluxVectorPDESC(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        void GetViscousFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &derivatives,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &viscousTensor);

        void SetCommonBC(const std::string &userDefStr,
                         const int n,
                         const NekDouble time,
                         int &cnt,
                         Array<OneD, Array<OneD, NekDouble> > &inarray);
        
        void WallBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void WallViscousBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void SymmetryBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void RiemannInvariantBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void PressureOutflowNonReflectiveBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void PressureOutflowBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void PressureOutflowFileBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void PressureInflowFileBC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void ExtrapOrder0BC(
            int                                                 bcRegion,
            int                                                 cnt,
            Array<OneD, Array<OneD, NekDouble> >               &physarray);
        void GetVelocityVector(
            const Array<OneD,       Array<OneD, NekDouble> > &physfield,
                  Array<OneD,       Array<OneD, NekDouble> > &velocity);
        void GetSoundSpeed(
            const Array<OneD,       Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &soundspeed);
        void GetMach(
                  Array<OneD,       Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &soundspeed,
                  Array<OneD,                   NekDouble>   &mach);
        void GetTemperature(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &temperature);
        void GetPressure(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure);
        void GetPressure(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
            const Array<OneD, const Array<OneD, NekDouble> > &velocity,
                  Array<OneD,                   NekDouble>   &pressure);
        void GetEnthalpy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &pressure,
                  Array<OneD,                   NekDouble>   &enthalpy);
        void GetEntropy(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
            const Array<OneD, const             NekDouble>   &pressure,
            const Array<OneD, const             NekDouble>   &temperature,
                  Array<OneD,                   NekDouble>   &entropy);
        void GetSmoothArtificialViscosity(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD,             NekDouble  > &eps_bar);
        void GetDynamicViscosity(
            const Array<OneD, const NekDouble> &temperature,
                  Array<OneD,       NekDouble> &mu);
        void GetStdVelocity(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,                   NekDouble>   &stdV);

        virtual bool v_PostIntegrate(int step);
        bool CalcSteadyState(bool output);

        void GetSensor(
            const Array<OneD, const Array<OneD, NekDouble> > &physarray,
                  Array<OneD,                   NekDouble>   &Sensor,
                  Array<OneD,                   NekDouble>   &SensorKappa);
        void GetElementDimensions(
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
                  Array<OneD,                   NekDouble >  &hmin);
        void GetAbsoluteVelocity(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,                   NekDouble>   &Vtot);
        void GetArtificialDynamicViscosity(
            const Array<OneD,  Array<OneD, NekDouble> > &physfield,
                  Array<OneD,              NekDouble>   &mu_var);
        void SetVarPOrderElmt(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD,                   NekDouble>   &PolyOrder);
        void GetForcingTerm(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > outarrayForcing);
        virtual NekDouble v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray);
        virtual void v_SetInitialConditions(
            NekDouble initialtime           = 0.0,
            bool      dumpInitialConditions = true,
            const int domain                = 0);

        NekDouble GetGasConstant()
        {
            return m_gasConstant;
        }

        NekDouble GetGamma()
        {
            return m_gamma;
        }

        const Array<OneD, const Array<OneD, NekDouble> > &GetVecLocs()
        {
            return m_vecLocs;
        }

        const Array<OneD, const Array<OneD, NekDouble> > &GetNormals()
        {
            return m_traceNormals;
        }

        virtual void v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables);
    };
}
#endif
