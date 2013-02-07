///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.h
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_ADVECTION
#define NEKTAR_SOLVERUTILS_ADVECTION

#include <string>
#include <boost/function.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        typedef boost::function<void (
            const int, 
            const Array<OneD, Array<OneD, NekDouble> >&,
            Array<OneD, Array<OneD, NekDouble> >&)> AdvectionFluxVecCB;
        
        class Advection
        {
        public:
            SOLVER_UTILS_EXPORT void InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            SOLVER_UTILS_EXPORT void Advect(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray);

            template<typename FuncPointerT, typename ObjectPointerT> 
            void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxVector = boost::bind(func, obj, _1, _2, _3);
            }
            
            inline void SetRiemannSolver(RiemannSolverSharedPtr riemann)
            {
                m_riemann = riemann;
            }
            
        protected:
            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
            {
                
            };
                        
            virtual void v_Advect(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &advVel,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray)=0;
            
            AdvectionFluxVecCB     m_fluxVector;
            RiemannSolverSharedPtr m_riemann;
        }; 
        
        /// A shared pointer to an EquationSystem object
        typedef boost::shared_ptr<Advection> AdvectionSharedPtr;
        
        /// Datatype of the NekFactory used to instantiate classes derived
        /// from the Advection class.
        typedef LibUtilities::NekFactory<std::string, Advection, std::string> AdvectionFactory;
        SOLVER_UTILS_EXPORT AdvectionFactory& GetAdvectionFactory();
    }
}

#endif
