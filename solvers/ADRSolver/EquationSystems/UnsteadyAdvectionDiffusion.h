///////////////////////////////////////////////////////////////////////////////
//
// File UnsteadyAdvectionDiffusion.h
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
// Description: Unsteady advection-diffusion solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_UNSTEADYADVECTIONDIFFUSION_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>

namespace Nektar
{
    class UnsteadyAdvectionDiffusion : public SolverUtils::AdvectionSystem
    {
    public:
        friend class MemoryManager<UnsteadyAdvectionDiffusion>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession) {
            SolverUtils::EquationSystemSharedPtr p
                = MemoryManager<UnsteadyAdvectionDiffusion>::
                                        AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        /// Destructor
        virtual ~UnsteadyAdvectionDiffusion();

    protected:
        bool m_useSpecVanVisc;
        NekDouble m_sVVCutoffRatio;   // cut off ratio from which to start decayhing modes
        NekDouble m_sVVDiffCoeff;     // Diffusion coefficient of SVV modes
        SolverUtils::RiemannSolverSharedPtr     m_riemannSolver;
        SolverUtils::DiffusionSharedPtr         m_diffusion;
        Array<OneD, Array<OneD, NekDouble> >    m_velocity;
        Array<OneD, NekDouble>                  m_traceVn;

        // Plane (used only for Discontinous projection
        //        with 3DHomogenoeus1D expansion)
        int                                     m_planeNumber;

        /// Session reader
        UnsteadyAdvectionDiffusion(
            const LibUtilities::SessionReaderSharedPtr& pSession);

        /// Evaluate the flux at each solution point for the advection part
        void GetFluxVectorAdv(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        /// Evaluate the flux at each solution point for the diffusion part
        void GetFluxVectorDiff(
            const int i, 
            const int j,
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> > &derivatives,
                  Array<OneD, Array<OneD, NekDouble> > &flux);

        /// Compute the RHS
        virtual void DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time);

        /// Perform the projection
        void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

        /// Solve implicitly the diffusion term
        virtual void DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD,       Array<OneD, NekDouble> >&outarray,
            NekDouble time,
            NekDouble lambda);

        /// Get the normal velocity
        Array<OneD, NekDouble> &GetNormalVelocity();

        /// Initialise the object
        virtual void v_InitObject();

        /// Print Summary
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

    private:
        NekDouble m_waveFreq;
        NekDouble m_epsilon;
    };
}

#endif
