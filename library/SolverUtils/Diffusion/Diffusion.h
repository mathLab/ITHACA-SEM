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

        /**
         * Parameter list meaning:
         *  1th: Devrivatives of field conservative varialbes
         *  2rd: the current time for time-dependent boundary
         *  3th: Fwd of field conservative variables        optional
         *  4th: Fwd of Devrivatives(2nd)                   optional 
         * 
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const Array<OneD, const Array<OneD, NekDouble> >                    &,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &,
            NekDouble                                                            ,
            const Array<OneD, const Array<OneD, NekDouble> >                    &,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &)> FunctorDerivBndCond;
        
        /**
         * Parameter list meaning:
         *  1nd: nspaceDimension
         *  2rd: field conservative variables
         *  3th: Devrivatives of field conservative varialbes
         *  4th: nonzero flux index array,          optional
         *  5th: normal vectors                     optional 
         *  6th: aritificial diffusion facotres     optional
         * 
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array< OneD, int >                                        &,    
            const Array<OneD, Array<OneD, NekDouble> >                      &,           
            const Array<OneD, NekDouble>                                    &)> DiffusionFluxCons;

        /**
         * Parameter list meaning:
         *  1nd: nspaceDimension
         *  2rd: trace conservative variables for Diffusion Flux Jacobian
         *  3th: trace conservative variables( usually the jump of trace value)
         *  4th: trace symmetric flux
         *  5th: nonzero flux index array,          optional
         *  6th: normal vectors                     optional 
         * 
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, NekDouble > >                     &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array< OneD, int >                                        &,    
            const Array<OneD, Array<OneD, NekDouble> >                      &)> DiffusionSymmFluxCons;

        /**
         * Parameter list meaning:
         *  1rd: trace conservative variables
         */
        typedef std::function<void (
                  Array<OneD,       Array<OneD, NekDouble> >                &)> SpecialBndTreat;

        class Diffusion
        {
        public:

            SOLVER_UTILS_EXPORT virtual ~Diffusion()
            {};

            SOLVER_UTILS_EXPORT void InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields);

            SOLVER_UTILS_EXPORT void Diffuse(
                const std::size_t                                 nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse_coeff(
                const std::size_t                                 nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse(
                const std::size_t                                 nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                NekDouble                                           time,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse_coeff(
                const std::size_t                                 nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                    Array<OneD, Array<OneD, NekDouble> >        &outarray,
                NekDouble                                           time,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd = NullNekDoubleArrayofArray);
            
            // Diffusion Calculate the physical derivatives
            SOLVER_UTILS_EXPORT void DiffuseCalculateDerivative(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
                const Array<OneD, Array<OneD, NekDouble>>         &pFwd =NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>         &pBwd =NullNekDoubleArrayofArray)
            {
                v_DiffuseCalculateDerivative(fields, inarray,inarrayderivative,pFwd, pBwd);
            }
            
            // Diffusion Volume FLux
            SOLVER_UTILS_EXPORT void DiffuseVolumeFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array< OneD, int >                                  &nonZeroIndex       =   NullInt1DArray)    
            {
                v_DiffuseVolumeFlux(fields, inarray,inarrayderivative,VolumeFlux,nonZeroIndex);
            }
            
            // Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT void DiffuseTraceFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                    Array<OneD, Array<OneD, NekDouble> >            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd           =   NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd           =   NullNekDoubleArrayofArray,
                Array< OneD, int >                                  &nonZeroIndex   =   NullInt1DArray)    
            {
                v_DiffuseTraceFlux(fields, inarray,inarrayderivative,VolumeFlux,TraceFlux,pFwd, pBwd,nonZeroIndex);
            }

            SOLVER_UTILS_EXPORT void AddDiffusionSymmFluxToCoeff(
                const std::size_t                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                v_AddDiffusionSymmFluxToCoeff(nConvectiveFields,fields,inarray,qfield,VolumeFlux,outarray,pFwd,pBwd);
            }

            SOLVER_UTILS_EXPORT void AddDiffusionSymmFluxToPhys(
                const std::size_t                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                v_AddDiffusionSymmFluxToPhys(nConvectiveFields,fields,inarray,qfield,VolumeFlux,outarray,pFwd,pBwd);
            }

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

            SOLVER_UTILS_EXPORT void SetFluxVector(DiffusionFluxVecCB fluxVector)
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

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetDiffusionFluxCons(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorDiffusionfluxCons = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                               std::placeholders::_5, std::placeholders::_6,
                               std::placeholders::_7);
            }
            
            void SetDiffusionFluxCons(DiffusionFluxCons flux)
            {
                m_FunctorDiffusionfluxCons =   flux;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFunctorDerivBndCond(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorDerivBndCond = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                                std::placeholders::_5);
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetArtificialDiffusionVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_ArtificialDiffusionVector = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2);
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetSpecialBndTreat(FuncPointerT func, ObjectPointerT obj)
            {
                m_SpecialBndTreat = std::bind(
                    func, obj, std::placeholders::_1);
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetDiffusionSymmFluxCons(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorSymmetricfluxCons = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                               std::placeholders::_5, std::placeholders::_6);
            }

            SOLVER_UTILS_EXPORT void SetFunctorDerivBndCond(FunctorDerivBndCond DerivBndCond)
            {
                m_FunctorDerivBndCond = DerivBndCond;
            }

            SOLVER_UTILS_EXPORT inline void SetHomoDerivs(Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                v_SetHomoDerivs(deriv);
            }

            SOLVER_UTILS_EXPORT virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &GetFluxTensor()
            {
                return v_GetFluxTensor();
            }
            SOLVER_UTILS_EXPORT void GetAVmu(
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                      Array<OneD, NekDouble >                               &muvar,
                      Array<OneD, NekDouble >                               &MuVarTrace);

            SOLVER_UTILS_EXPORT void ConsVarAveJump(
                const std::size_t                                   nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump)
            {
                v_ConsVarAveJump(nConvectiveFields,npnts,vFwd,vBwd,aver,jump);
            }

            SOLVER_UTILS_EXPORT const Array<OneD, const Array<OneD, NekDouble> > &GetTraceNormal()
            {
                return v_GetTraceNormal();
            }
            
        protected:
            DiffusionFluxVecCB              m_fluxVector;
            DiffusionFluxVecCBNS            m_fluxVectorNS;
            DiffusionFluxPenaltyNS          m_fluxPenaltyNS;
            DiffusionArtificialDiffusion    m_ArtificialDiffusionVector;
            DiffusionFluxCons               m_FunctorDiffusionfluxCons;
            FunctorDerivBndCond             m_FunctorDerivBndCond;
            SpecialBndTreat                 m_SpecialBndTreat;
            DiffusionSymmFluxCons           m_FunctorSymmetricfluxCons;

            NekDouble                       m_time=0.0;

            SOLVER_UTILS_EXPORT virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
            {
                boost::ignore_unused(pSession, pFields);
            };

            SOLVER_UTILS_EXPORT virtual void v_Diffuse(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd);
            
            SOLVER_UTILS_EXPORT virtual void v_Diffuse_coeff(
                const std::size_t                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> > &pFwd,
                const Array<OneD, Array<OneD, NekDouble> > &pBwd);
        
            SOLVER_UTILS_EXPORT virtual void v_ConsVarAveJump(
                const std::size_t                                   nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump);

            // Diffusion Flux, calculate the physical derivatives
            SOLVER_UTILS_EXPORT virtual void v_DiffuseCalculateDerivative(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            // Diffusion Volume Flux
            SOLVER_UTILS_EXPORT virtual void v_DiffuseVolumeFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array< OneD, int >                                  &nonZeroIndex) ;
            
             // Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT virtual void v_DiffuseTraceFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
                Array< OneD, int >                                  &nonZeroIndex);

            SOLVER_UTILS_EXPORT virtual void v_AddDiffusionSymmFluxToCoeff(
                const std::size_t                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                boost::ignore_unused(nConvectiveFields,fields,inarray,qfield,VolumeFlux,outarray,pFwd,pBwd);

            }
            SOLVER_UTILS_EXPORT virtual void v_AddDiffusionSymmFluxToPhys(
                const std::size_t                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                boost::ignore_unused(nConvectiveFields,fields,inarray,qfield,VolumeFlux,outarray,pFwd,pBwd);

            }
            
            SOLVER_UTILS_EXPORT virtual void v_SetHomoDerivs(
                Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                boost::ignore_unused(deriv);
            }

            SOLVER_UTILS_EXPORT virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &v_GetFluxTensor()
            {
                static Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmp;
                return tmp;
            }

            SOLVER_UTILS_EXPORT virtual const Array<OneD, const Array<OneD, NekDouble> > &v_GetTraceNormal();
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
