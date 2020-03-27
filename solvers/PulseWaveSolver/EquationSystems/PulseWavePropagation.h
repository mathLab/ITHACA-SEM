///////////////////////////////////////////////////////////////////////////////
//
// File PulseWavePropagation.h
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
// Description: Pulse Wave Propagation solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PULSEWAVESOLVER_EQUATIONSYSTEMS_PULSEWAVEPROPAGATION_H
#define NEKTAR_SOLVERS_PULSEWAVESOLVER_EQUATIONSYSTEMS_PULSEWAVEPROPAGATION_H

#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>
#include <PulseWaveSolver/EquationSystems/PulseWaveBoundary.h>
#include <PulseWaveSolver/EquationSystems/PulseWavePressureArea.h>
#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
    class PulseWavePropagation : public PulseWaveSystem
    {
    public:
        friend class MemoryManager<PulseWavePropagation>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph)
        {
            EquationSystemSharedPtr p = MemoryManager<PulseWavePropagation>
                ::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }
	
        /// Name of class
        static std::string className;
        
        virtual ~PulseWavePropagation();
	
        // Functions for Riemann solver
        Array<OneD, NekDouble> &GetA0();
        Array<OneD, NekDouble> &GetBeta();
        Array<OneD, NekDouble> &GetN();
        NekDouble               GetRho();
        NekDouble               GetPext();
    protected:
        PulseWavePropagation(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph);

        void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                      Array<OneD,  Array<OneD, NekDouble> > &outarray,
                      const NekDouble time);
        
        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                             Array<OneD,  Array<OneD, NekDouble> > &outarray,
                             const NekDouble time);
        
        void SetPulseWaveBoundaryConditions(const Array<OneD,const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD, Array<OneD, NekDouble> >&outarray, 
                                            const NekDouble time);
        virtual void v_InitObject();

        /// DG Pulse Wave Propagation routines:
        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        SolverUtils::RiemannSolverSharedPtr     m_riemannSolver;
        SolverUtils::AdvectionSharedPtr         m_advObject;

        Array<OneD, PulseWaveBoundarySharedPtr> m_Boundary;

        PulseWavePressureAreaSharedPtr m_pressureArea;
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);
    };
}

#endif
