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

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.h>
#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <CompressibleFlowSolver/BoundaryConditions/CFSBndCond.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>


namespace Nektar
{
    /**
     *
     */
    class CompressibleFlowSystem: public SolverUtils::AdvectionSystem,
                                  public SolverUtils::FluidInterface
    {
    public:

        friend class MemoryManager<CompressibleFlowSystem>;

        virtual ~CompressibleFlowSystem();

        /// Function to calculate the stability limit for DG/CG.
        NekDouble GetStabilityLimit(int n);

        /// Function to calculate the stability limit for DG/CG
        /// (a vector of them).
        Array<OneD, NekDouble> GetStabilityLimitVector(
            const Array<OneD,int> &ExpOrder);

        /// Function to get estimate of min h/p factor per element
        Array<OneD, NekDouble>  GetElmtMinHP(void);

        virtual void GetPressure(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &pressure);

        virtual void GetDensity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &density);

        virtual bool HasConstantDensity()
        {
            return false;
        }

        virtual void GetVelocity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> >       &velocity);

    protected:
        SolverUtils::DiffusionSharedPtr     m_diffusion;
        ArtificialDiffusionSharedPtr        m_artificialDiffusion;
        Array<OneD, Array<OneD, NekDouble> >m_vecLocs;
        NekDouble                           m_gamma;
        std::string                         m_shockCaptureType;

        // Parameters for exponential filtering
        NekDouble                           m_filterAlpha;
        NekDouble                           m_filterExponent;
        NekDouble                           m_filterCutoff;
        bool                                m_useFiltering;

        // Parameters for local time-stepping
        bool                                m_useLocalTimeStep;

        // Auxiliary object to convert variables
        VariableConverterSharedPtr          m_varConv;

        // User defined boundary conditions
        std::vector<CFSBndCondSharedPtr>    m_bndConds;

        // Forcing term
        std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

        NekDouble                           m_BndEvaluateTime;

        CompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph);

        virtual void v_InitObject();

        void InitialiseParameters();

        void InitAdvection();

        void DoOdeRhs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time);
        void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time);

        void DoAdvection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        void DoDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);
        void DoDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> >             &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >       &physfield,
            TensorOfArray3D<NekDouble>                       &flux);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >       &physfield,
            TensorOfArray3D<NekDouble>                       &flux);

        void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time);
        
        void SetBoundaryConditionsBwdWeight();

        void GetElmtTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD, NekDouble> &tstep);

        virtual NekDouble v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray);
        virtual void v_SetInitialConditions(
            NekDouble initialtime           = 0.0,
            bool      dumpInitialConditions = true,
            const int domain                = 0);

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

        virtual void v_DoDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        virtual void v_DoDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity();

    };
}
#endif
