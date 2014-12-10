///////////////////////////////////////////////////////////////////////////////
//
// File APE.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2014 Kilian Lackhove
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
// Description: Acoustic perturbation equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_H
#define NEKTAR_SOLVERS_APESOLVER_EQUATIONSYSTEMS_APE_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{     

class APE : public UnsteadySystem
{
    public:

        friend class MemoryManager<APE>;

        /// Creates an instance of this class
        static EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession)
        {
            EquationSystemSharedPtr p = MemoryManager<APE>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        /// Destructor
        virtual ~APE();


    protected:

        SolverUtils::AdvectionSharedPtr                 m_advection;
        SolverUtils::RiemannSolverSharedPtr             m_riemannSolver;
        Array<OneD, Array<OneD, NekDouble> >            m_traceBasefield;
        Array<OneD, Array<OneD, NekDouble> >            m_vecLocs;
        /// Constant incompressible density (APE)
        NekDouble                                       m_Rho0;
        /// Isentropic coefficient, Ratio of specific heats (APE)
        NekDouble                                       m_gamma;
        Array<OneD, Array<OneD, NekDouble> >            m_basefield;
        std::vector<std::string>                        m_basefield_names;

        /// Initialises UnsteadySystem class members.
        APE(const LibUtilities::SessionReaderSharedPtr& pSession);

        virtual void v_InitObject();

        /// Sets up initial conditions.
        virtual void v_DoInitialise();

        void DoOdeRhs(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                            Array<OneD,  Array<OneD, NekDouble> > &outarray,
                      const NekDouble time);

        void DoOdeProjection(const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                                   Array<OneD,  Array<OneD, NekDouble> > &outarray,
                             const NekDouble time);

        void GetFluxVector(
                const Array<OneD, Array<OneD, NekDouble> > &physfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        void AddSource(Array< OneD, Array< OneD, NekDouble > >& outarray);

        virtual void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
                                      std::vector<std::string>             &variables);

        const Array<OneD, const Array<OneD, NekDouble> > &GetNormals();

        const Array<OneD, const Array<OneD, NekDouble> > &GetVecLocs();

        const Array<OneD, const Array<OneD, NekDouble> > &GetBasefield();

        NekDouble GetGamma();

        NekDouble GetRho();

        void UpdateBasefield();

    private:

        void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &physarray, NekDouble time);

        void WallBC(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray);
};
}

#endif 

