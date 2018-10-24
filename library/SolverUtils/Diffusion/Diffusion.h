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
// Description: Abstract base class for diffusion.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DIFFUSION
#define NEKTAR_SOLVERUTILS_DIFFUSION

#include <string>
#include <functional>

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
            const int, 
            const int, 
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD, Array<OneD, NekDouble> >&)> DiffusionFluxVecCB;
        
        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >&)> 
                                            DiffusionFluxVecCBNS;
        
        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
                  Array<OneD,             NekDouble  >&)>
                                            DiffusionArtificialDiffusion;

        /**
         * Parameter list meaning:
         *  1rd: field conservative variables
         *  2th: Devrivatives of field conservative varialbes
         *  3rd: the current time for time-dependent boundary
         *  4th: Fwd of field conservative variables        optional
         *  5th: Fwd of Devrivatives(2nd)                   optional 
         * 
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const Array<OneD, const Array<OneD, NekDouble> >                    &,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &,
            NekDouble                                                            ,
            const Array<OneD, const Array<OneD, NekDouble> >                    &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &)> FunctorDerivBndCond;
        
        /**
         * Parameter list meaning:
         *  1st: nvariables
         *  2nd: nspaceDimension
         *  3rd: field conservative variables
         *  4th: Devrivatives of field conservative varialbes
         *  5th: nonzero flux index array,          optional
         *  6th: normal vectors                     optional 
         *  7th: aritificial diffusion facotres     optional
         * 
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const int                                                       ,
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array< OneD, int >                                        &,    
            const Array<OneD, Array<OneD, NekDouble> >                      &,           
            const Array<OneD, NekDouble>                                    &)> DiffusionFluxVec;

        class Diffusion
        {
        public:
            SOLVER_UTILS_EXPORT void InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);
                        
            SOLVER_UTILS_EXPORT void Diffuse(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse_coeff(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd= NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd= NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                NekDouble                                           time,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse_coeff(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >        &outarray,
                NekDouble                                           time,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd= NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd= NullNekDoubleArrayofArray);
            
            SOLVER_UTILS_EXPORT void FluxVec(
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
                                                                &fluxvector);
            
            template<typename FuncPointerT, typename ObjectPointerT> 
            void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxVector = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, std::placeholders::_4,
                    std::placeholders::_5);
            }
            
            void SetFluxVectorVec(DiffusionFluxVecCB fluxVector)
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

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetArtificialDiffusionVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_ArtificialDiffusionVector = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2);
            }

            void SetFluxVectorNS(DiffusionFluxVecCBNS fluxVector)
            {
                m_fluxVectorNS = fluxVector;
            }

            void SetDiffusionFluxVec(DiffusionFluxVec flux)
            {
                m_Diffusionflux =   flux;
            }

            void SetFunctorDerivBndCond(FunctorDerivBndCond DerivBndCond)
            {
                m_FunctorDerivBndCond = DerivBndCond;
            }

            void DiffusionFlux(
                const int                                                       nConvectiveFields,
                const int                                                       nDim,
                const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &qfields,
                    Array<OneD, Array<OneD, NekDouble> >                        &outarray,
                    Array< OneD, int >                                          &nonZeroIndex       = NullInt1DArray,
                const Array<OneD, Array<OneD, NekDouble> >                      &normal             = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >                      &ArtifDiffFactor    = NullNekDoubleArrayofArray);

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
            DiffusionArtificialDiffusion    m_ArtificialDiffusionVector;
            DiffusionFluxVec                m_Diffusionflux;
            FunctorDerivBndCond             m_FunctorDerivBndCond;

            NekDouble                       m_time=0.0;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
            {
                
            };
                        
            virtual void v_Diffuse(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray)=0;

            virtual void v_Diffuse_coeff(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            virtual void v_SetHomoDerivs(
                Array<OneD, Array<OneD, NekDouble> > &deriv)
            {

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
