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
            const TensorOfArray3D<NekDouble>           &,
                  TensorOfArray3D<NekDouble>           &)>
                                            DiffusionFluxVecCB;

        typedef std::function<void (
            const Array<OneD, Array<OneD, NekDouble> >&,
                  TensorOfArray3D<NekDouble>          &,
                  TensorOfArray3D<NekDouble>          &)>
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
         *  1st: field conservative variables
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
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &)> FunctorDerivBndCond;

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
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const TensorOfArray3D<NekDouble>                                &,
                  TensorOfArray3D<NekDouble>                                &,
                  Array< OneD, int >                                        &,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, NekDouble>                                    &)>
        DiffusionFluxCons;

        /**
         * Parameter list meaning:
         *  1st: nvariables
         *  2nd: nspaceDimension
         *  3rd: trace conservative variables for Diffusion Flux Jacobian
         *  4th: trace conservative variables( usually the jump of trace value)
         *  5th: trace symmetric flux
         *  6th: nonzero flux index array,          optional
         *  7th: normal vectors                     optional
         *
         * a null pointer need to be passed for optional parameters
         */
        typedef std::function<void (
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, NekDouble > >                     &,
                  TensorOfArray3D<NekDouble>                                &,
                  Array< OneD, int >                                        &,
            const Array<OneD, Array<OneD, NekDouble> >                      &)>
                DiffusionSymmFluxCons;

        /**
         * Parameter list meaning:
         *  1rd: trace conservative variables
         */
        typedef std::function<void (
                  Array<OneD,       Array<OneD, NekDouble> >                &)>
                  SpecialBndTreat;

        /**
         * Parameter list meaning:
         *  1st: trace conservative variables
         *  2rd: dynamic viscosity
         */
        typedef std::function<void (
            const Array<OneD, const Array<OneD, NekDouble> >                &,
            Array<OneD, NekDouble>                                          &)>
        CalcViscosity;

        class Diffusion
        {
        public:
            ///Params for Ducros sensor
            Array<OneD, NekDouble> m_divVel;
            Array<OneD, NekDouble> m_divVelSquare;
            Array<OneD, NekDouble> m_curlVelSquare;

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

            SOLVER_UTILS_EXPORT void DiffuseCoeffs(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pBwd = NullNekDoubleArrayofArray);

            SOLVER_UTILS_EXPORT void Diffuse(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                NekDouble                                         time,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pBwd = NullNekDoubleArrayofArray);
            
            SOLVER_UTILS_EXPORT void DiffuseCoeffs(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
                TensorOfArray3D<NekDouble>                        &qfield,
                Array< OneD, int >                                &nonZeroIndex)
            {
                v_DiffuseCoeffs(nConvectiveFields, fields, inarray, outarray,
                                pFwd, pBwd,qfield,nonZeroIndex);
            }

            SOLVER_UTILS_EXPORT void DiffuseCoeffs(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                NekDouble                                         time,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >
                    &pBwd = NullNekDoubleArrayofArray);
            
            // Diffusion Calculate the physical derivatives
            SOLVER_UTILS_EXPORT void DiffuseCalcDerivative(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                TensorOfArray3D<NekDouble>                        &qfields, 
                const Array<OneD, Array<OneD, NekDouble>>         &pFwd
                      =NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>         &pBwd
                      =NullNekDoubleArrayofArray);

            /// Diffusion Volume FLux

            SOLVER_UTILS_EXPORT void DiffuseVolumeFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array< OneD, int >
                            &nonZeroIndex       =   NullInt1DArray)
            {
                v_DiffuseVolumeFlux(fields, inarray, qfields,
                                    VolumeFlux, nonZeroIndex);
            }

            /// Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT void DiffuseTraceFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>
                    &pFwd           =   NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>
                    &pBwd           =   NullNekDoubleArrayofArray,
                Array< OneD, int >
                    &nonZeroIndex   =   NullInt1DArray)
            {
                v_DiffuseTraceFlux(fields, inarray, qfields,
                    VolumeFlux, TraceFlux, pFwd, pBwd, nonZeroIndex);
            }

            SOLVER_UTILS_EXPORT void DiffuseTraceSymmFlux(
                const int                                 nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &fields,
                const Array<OneD, Array<OneD, NekDouble>>          &inarray,
                const TensorOfArray3D<NekDouble>                   &qfield,
                const TensorOfArray3D<NekDouble>                   &VolumeFlux,
                TensorOfArray3D<NekDouble>                         &SymmFlux,
                const Array<OneD, Array<OneD, NekDouble>>          &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>          &pBwd,
                Array< OneD, int >                                 &nonZeroIndex,
                Array<OneD, Array<OneD, NekDouble> >               &solution_Aver = NullNekDoubleArrayofArray,
                Array<OneD, Array<OneD, NekDouble> >               &solution_jump = NullNekDoubleArrayofArray)
            {
                v_DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,
                                       VolumeFlux,SymmFlux,pFwd,pBwd,
                                       nonZeroIndex,solution_Aver,solution_jump);
            }

            /// Add symmetric flux to field in coeff space
            SOLVER_UTILS_EXPORT void AddDiffusionSymmFluxToCoeff(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                TensorOfArray3D<NekDouble>                          &qfield,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                v_AddDiffusionSymmFluxToCoeff(nConvectiveFields, fields,
                    inarray, qfield, VolumeFlux, outarray, pFwd, pBwd);
            }

            /// Add symmetric flux to field in coeff physical space
            SOLVER_UTILS_EXPORT void AddDiffusionSymmFluxToPhys(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                TensorOfArray3D<NekDouble>                          &qfield,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                v_AddDiffusionSymmFluxToPhys(nConvectiveFields, fields, inarray,
                    qfield, VolumeFlux, outarray, pFwd, pBwd);
            }

            SOLVER_UTILS_EXPORT void AddSymmFluxIntegralToOffDiag(
                const int                                      nConvectiveFields,
                const int                                           nDim,
                const int                                           nPts,
                const int                                           nTracePts,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, const int >                       &nonZeroIndex,
                TensorOfArray3D<NekDouble>                          &Fwdflux,
                TensorOfArray3D<NekDouble>                          &Bwdflux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray);

            SOLVER_UTILS_EXPORT void FluxVec(
                    TensorOfArray3D<NekDouble> &fluxvector);

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
            {
                m_fluxVector = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3);
            }

            SOLVER_UTILS_EXPORT void SetFluxVector(
                DiffusionFluxVecCB fluxVector)
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
            void SetDiffusionFluxConsTrace(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorDiffusionfluxConsTrace = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                               std::placeholders::_5, std::placeholders::_6,
                               std::placeholders::_7);
            }

            void SetDiffusionFluxConsTrace(DiffusionFluxCons flux)
            {
                m_FunctorDiffusionfluxConsTrace =   flux;
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetArtificialDiffusionVector(
                FuncPointerT func, ObjectPointerT obj)
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
            void SetCalcViscosity(FuncPointerT func, ObjectPointerT obj)
            {
                m_CalcViscosity = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2);
            }

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetDiffusionSymmFluxCons(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorSymmetricfluxCons = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                               std::placeholders::_5, std::placeholders::_6);
            }

            inline void SetHomoDerivs(Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                v_SetHomoDerivs(deriv);
            }

            virtual TensorOfArray3D<NekDouble> &GetFluxTensor()
            {
                return v_GetFluxTensor();
            }
            /// Get the mu of artifical viscosity(AV)
            SOLVER_UTILS_EXPORT void GetAVmu(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, NekDouble >                           &muvar,
                Array<OneD, NekDouble >                           &MuVarTrace);

            /// Get the average and jump value of conservative variables on trace
            SOLVER_UTILS_EXPORT void ConsVarAveJump(
                const std::size_t                             nConvectiveFields,
                const size_t                                        npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump)
            {
                v_ConsVarAveJump(nConvectiveFields, npnts, vFwd, vBwd, aver,
                    jump);
            }

            /// Get trace normal
            SOLVER_UTILS_EXPORT const Array<OneD, const Array<OneD, NekDouble> >
                &GetTraceNormal()
            {
                return v_GetTraceNormal();
            }

        protected:
            DiffusionFluxVecCB              m_fluxVector;
            DiffusionFluxVecCBNS            m_fluxVectorNS;
            DiffusionFluxPenaltyNS          m_fluxPenaltyNS;
            DiffusionArtificialDiffusion    m_ArtificialDiffusionVector;
            DiffusionFluxCons               m_FunctorDiffusionfluxCons;
            DiffusionFluxCons               m_FunctorDiffusionfluxConsTrace;
            SpecialBndTreat                 m_SpecialBndTreat;
            CalcViscosity                   m_CalcViscosity;
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
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            SOLVER_UTILS_EXPORT virtual void v_DiffuseCoeffs(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            SOLVER_UTILS_EXPORT virtual void v_DiffuseCoeffs(
                const std::size_t                            nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                const Array<OneD, Array<OneD, NekDouble> >        &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &vBwd,
                TensorOfArray3D<NekDouble>                        &qfield,
                Array< OneD, int >                               &nonZeroIndex);
        

            SOLVER_UTILS_EXPORT virtual void v_ConsVarAveJump(
                const std::size_t                             nConvectiveFields,
                const size_t                                        npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump);

            /// Diffusion Flux, calculate the physical derivatives
            SOLVER_UTILS_EXPORT virtual void v_DiffuseCalcDerivative(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                TensorOfArray3D<NekDouble>                        &qfields,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd);

            /// Diffusion Volume Flux
            SOLVER_UTILS_EXPORT virtual void v_DiffuseVolumeFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array< OneD, int >
                    &nonZeroIndex) ;

             /// Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT virtual void v_DiffuseTraceFlux(
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                TensorOfArray3D<NekDouble>                          &Qfields,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
                Array< OneD, int >
                    &nonZeroIndex);

            SOLVER_UTILS_EXPORT virtual void v_DiffuseTraceSymmFlux(
                const int                                    nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>  &fields,
                const Array<OneD, Array<OneD, NekDouble>>          &inarray,
                const TensorOfArray3D<NekDouble>                   &qfield,
                const TensorOfArray3D<NekDouble>                   &VolumeFlux,
                TensorOfArray3D<NekDouble>                         &SymmFlux,
                const Array<OneD, Array<OneD, NekDouble>>          &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>          &pBwd,
                Array< OneD, int >                                 &nonZeroIndex,
                Array<OneD, Array<OneD, NekDouble> >               &solution_Aver,
                Array<OneD, Array<OneD, NekDouble> >               &solution_jump)
            {
                boost::ignore_unused(nConvectiveFields, fields, inarray, qfield,
                                     VolumeFlux, SymmFlux, pFwd, pBwd,
                                     nonZeroIndex,solution_Aver,solution_jump);
                ASSERTL0(false, "Not defined function v_DiffuseTraceSymmFlux.");
            }

            SOLVER_UTILS_EXPORT virtual void v_AddDiffusionSymmFluxToCoeff(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                TensorOfArray3D<NekDouble>                          &qfield,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                boost::ignore_unused(nConvectiveFields, fields, inarray, qfield,
                                        VolumeFlux, outarray, pFwd, pBwd);

            }
            SOLVER_UTILS_EXPORT virtual void v_AddDiffusionSymmFluxToPhys(
                const std::size_t                             nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                TensorOfArray3D<NekDouble>                          &qfield,
                TensorOfArray3D<NekDouble>                          &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {
                boost::ignore_unused(nConvectiveFields, fields, inarray, qfield,
                                        VolumeFlux, outarray, pFwd, pBwd);

            }

            virtual void v_SetHomoDerivs(
                Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                boost::ignore_unused(deriv);
            }

            virtual TensorOfArray3D<NekDouble> &v_GetFluxTensor()
            {
                static TensorOfArray3D<NekDouble> tmp;
                return tmp;
            }

            SOLVER_UTILS_EXPORT virtual
                    const Array<OneD, const Array<OneD, NekDouble> >
                    &v_GetTraceNormal();

            /// Compute primary derivatives
            SOLVER_UTILS_EXPORT virtual void v_GetPrimVar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                  Array<OneD, Array<OneD, NekDouble>>         &primVar);

            /// Compute divergence and curl squared
            SOLVER_UTILS_EXPORT void GetDivCurl(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &pVarDer);
            
        };

    /// A shared pointer to an EquationSystem object
    typedef std::shared_ptr<SolverUtils::Diffusion> DiffusionSharedPtr;
    
    /// Datatype of the NekFactory used to instantiate classes derived
    /// from the Diffusion class.
    typedef LibUtilities::NekFactory<std::string, Diffusion, std::string>
    DiffusionFactory;
    SOLVER_UTILS_EXPORT DiffusionFactory& GetDiffusionFactory();
}
}

#endif
