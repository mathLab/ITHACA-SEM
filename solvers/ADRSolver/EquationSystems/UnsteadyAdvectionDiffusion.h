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

using namespace Nektar::SolverUtils;

namespace Nektar
{
    class UnsteadyAdvectionDiffusion : public UnsteadySystem
    {
    public:
        friend class MemoryManager<UnsteadyAdvectionDiffusion>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession) {
            EquationSystemSharedPtr p = MemoryManager<UnsteadyAdvectionDiffusion>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        virtual ~UnsteadyAdvectionDiffusion();

    protected:
        UnsteadyAdvectionDiffusion(
                const LibUtilities::SessionReaderSharedPtr& pSession);

        virtual void DoOdeRhs(const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time);

        virtual void DoImplicitSolve(const Array<OneD, const Array<OneD,      NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                      NekDouble time,
                      NekDouble lambda);

        virtual void v_InitObject();

        virtual void v_GetFluxVector(const int i, Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &flux);

        virtual void v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, Array<OneD, Array<OneD, NekDouble> > &numflux);

        // Print Summary
        virtual void v_PrintSummary(std::ostream &out);

    private:
        NekDouble m_waveFreq;
        NekDouble m_epsilon;
        Array<OneD, Array<OneD, NekDouble> > m_velocity;
    };
}

#endif
