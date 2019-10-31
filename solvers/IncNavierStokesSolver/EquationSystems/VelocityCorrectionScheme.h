///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionScheme.h
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
// Description: Velocity Correction Scheme header 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H
#define NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{
    class VelocityCorrectionScheme: public IncNavierStokes
    {
    public:

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p =
                MemoryManager<VelocityCorrectionScheme>::AllocateSharedPtr(
                    pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

        /// Constructor.
        VelocityCorrectionScheme(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph);

        virtual ~VelocityCorrectionScheme();

        virtual void v_InitObject();
        
        void SetUpPressureForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &fields,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt)
        {
            v_SetUpPressureForcing( fields, Forcing, aii_Dt);
        }

        void SetUpViscousForcing(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &Forcing,
                    const NekDouble aii_Dt)
        {
            v_SetUpViscousForcing( inarray, Forcing, aii_Dt);
        }
        
        void SolvePressure( const Array<OneD, NekDouble>  &Forcing)
        {
            v_SolvePressure( Forcing);
        }
        
        void SolveViscous(
                    const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble aii_Dt)
        {
            v_SolveViscous( Forcing, outarray, aii_Dt);
        }

        void SolveUnsteadyStokesSystem(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time,
                    const NekDouble a_iixDt);

        void EvaluateAdvection_SetPressureBCs(
                    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                    const NekDouble time)
        {
            v_EvaluateAdvection_SetPressureBCs( inarray, outarray, time);
        }

    protected:
        /// bool to identify if spectral vanishing viscosity is active.
        bool m_useHomo1DSpecVanVisc;
        /// bool to identify if spectral vanishing viscosity is active.
        bool m_useSpecVanVisc;
        /// cutt off ratio from which to start decayhing modes
        NekDouble m_sVVCutoffRatio;
        /// Diffusion coefficient of SVV modes
        NekDouble m_sVVDiffCoeff;
        NekDouble m_sVVCutoffRatioHomo1D;
        /// Diffusion coefficient of SVV modes in homogeneous 1D Direction
        NekDouble m_sVVDiffCoeffHomo1D;
        /// Array of coefficient if power kernel is used in SVV
        Array<OneD, NekDouble> m_svvVarDiffCoeff;
        /// Identifier for Power Kernel otherwise DG kernel
        bool m_IsSVVPowerKernel;
        /// Diffusion coefficients (will be kinvis for velocities)
        Array<OneD, NekDouble> m_diffCoeff;

        /// Variable Coefficient map for the Laplacian which can be activated as part of SVV or otherwise
        StdRegions::VarCoeffMap m_varCoeffLap;

        /// Desired volumetric flowrate
        NekDouble m_flowrate;
        /// Area of the boundary through which we are measuring the flowrate
        NekDouble m_flowrateArea;
        // Bool to identify 3D1HD with forcing explicitly defined
        bool m_homd1DFlowinPlane;
        /// Flux of the Stokes function solution
        NekDouble m_greenFlux;
        /// Current flowrate correction
        NekDouble m_alpha;
        /// Boundary ID of the flowrate reference surface
        int m_flowrateBndID;
        /// Plane ID for cases with homogeneous expansion
        int m_planeID;
        /// Flowrate reference surface
        MultiRegions::ExpListSharedPtr m_flowrateBnd;
        /// Stokes solution used to impose flowrate
        Array<OneD, Array<OneD, NekDouble> > m_flowrateStokes;
        /// Output stream to record flowrate
        std::ofstream m_flowrateStream;
        /// Interval at which to record flowrate data
        int m_flowrateSteps;
        /// Value of aii_dt used to compute Stokes flowrate solution.
        NekDouble m_flowrateAiidt;

        void SetupFlowrate(NekDouble aii_dt);
        NekDouble MeasureFlowrate(
            const Array<OneD, Array<OneD, NekDouble> > &inarray);

        // Virtual functions
        virtual bool v_PostIntegrate(int step);

        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

        virtual void v_TransCoeffToPhys(void);

        virtual void v_TransPhysToCoeff(void);

        virtual void v_DoInitialise(void);

        virtual Array<OneD, bool> v_GetSystemSingularChecks();

        virtual int v_GetForceDimension();
        
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

        virtual bool v_RequireFwdTrans()
        {
            return false;
        }

        virtual std::string v_GetExtrapolateStr(void)
        {
            return "Standard";
        }
        
        virtual std::string v_GetSubSteppingExtrapolateStr(
                                               const std::string &instr)
        {
            return instr;
        }
        
        Array<OneD, Array< OneD, NekDouble> > m_F;

        void SetUpSVV(void);
        void SetUpExtrapolation(void);
        
        void SVVVarDiffCoeff(const NekDouble velmag, 
                             Array<OneD, NekDouble> &diffcoeff,
                             const Array<OneD, Array<OneD, NekDouble> >
                             &vel = NullNekDoubleArrayofArray);
        void AppendSVVFactors(
                              StdRegions::ConstFactorMap &factors,
                              MultiRegions::VarFactorsMap &varFactorsMap);
    private:
        
    };

    typedef std::shared_ptr<VelocityCorrectionScheme>
                VelocityCorrectionSchemeSharedPtr;

} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H
