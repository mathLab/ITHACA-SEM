///////////////////////////////////////////////////////////////////////////////
//
// File: Diffusion.h
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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSION
#define NEKTAR_SOLVERUTILS_DIFFUSION

#include <string>
#include <functional>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> > &,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&)>
                                            DiffusionFluxVecCB;

        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&)>
                                            DiffusionFluxVecCBNS;

        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD, Array<OneD, NekDouble> >&)>
                                            DiffusionFluxPenaltyNS;

        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD,             NekDouble  >&)>
                                            DiffusionArtificialDiffusion;

        class Diffusion
        {
        public:

            SOLVER_UTILS_EXPORT virtual ~Diffusion()
            {};

            SOLVER_UTILS_EXPORT void InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            SOLVER_UTILS_EXPORT void Diffuse(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void FluxVec(
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                                &fluxvector);

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxVector = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3);
            }

            void SetFluxVector(DiffusionFluxVecCB fluxVector)
            {
                m_fluxVector = fluxVector;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFluxVectorNS(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxVectorNS = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3);
            }

            void SetFluxVectorNS(DiffusionFluxVecCBNS fluxVector)
            {
                m_fluxVectorNS = fluxVector;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFluxPenaltyNS(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxPenaltyNS = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3);
            }

            void SetFluxPenaltyNS(DiffusionFluxPenaltyNS flux)
            {
                m_fluxPenaltyNS = flux;
            }

            inline void SetHomoDerivs(Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                v_SetHomoDerivs(deriv);
            }

            virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &GetFluxTensor()
            {
                return v_GetFluxTensor();
            }

        protected:
            DiffusionFluxVecCB              m_fluxVector;
            DiffusionFluxVecCBNS            m_fluxVectorNS;
            DiffusionFluxPenaltyNS          m_fluxPenaltyNS;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
            {
                boost::ignore_unused(pSession, pFields);
            };

            virtual void v_Diffuse(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray)=0;

            virtual void v_SetHomoDerivs(
                Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                boost::ignore_unused(deriv);
            }

            virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &v_GetFluxTensor()
            {
                static Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmp;
                return tmp;
            }
        };

        /// A shared pointer to an EquationSystem object
        typedef std::shared_ptr<Diffusion> DiffusionSharedPtr;

        /// Datatype of the NekFactory used to instantiate classes derived
        /// from the Diffusion class.
        typedef LibUtilities::NekFactory<std::string, Diffusion, std::string> DiffusionFactory;
        SOLVER_UTILS_EXPORT DiffusionFactory& GetDiffusionFactory();
    }
}

#endif
