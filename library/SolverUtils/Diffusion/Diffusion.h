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
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array< OneD, int >                                        &,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, NekDouble>                                    &)> DiffusionFluxCons;

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
            const int                                                       ,
            const Array<OneD, Array<OneD, NekDouble> >                      &,
            const Array<OneD, Array<OneD, NekDouble > >                     &,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >        &,
                  Array< OneD, int >                                        &,
            const Array<OneD, Array<OneD, NekDouble> >                      &)> DiffusionSymmFluxCons;

        /**
         * Parameter list meaning:
         *  1st: nvariables
         *  2rd: trace conservative variables
         */
        typedef std::function<void (
            const int                                                       ,
                  Array<OneD,       Array<OneD, NekDouble> >                &)> SpecialBndTreat;

        /**
         * Parameter list meaning:
         *  1st: trace conservative variables
         *  2rd: dynamic viscosity
         */
        typedef std::function<void (
            const Array<OneD, const Array<OneD, NekDouble> >                &,
            Array<OneD, NekDouble>                                          &)> CalcViscosity;

        class Diffusion
        {
        public:
            ///Params for Ducros sensor
            Array<OneD, NekDouble> m_divVel;
            Array<OneD, NekDouble> m_divVelSquare;
            Array<OneD, NekDouble> m_curlVelSquare;

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
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
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
                Array<OneD, Array<OneD, NekDouble> >              &outarray,
                NekDouble                                           time,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd= NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd= NullNekDoubleArrayofArray);
            
            SOLVER_UTILS_EXPORT void Diffuse_coeff(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                Array<OneD, Array<OneD, NekDouble> >                        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >                  &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >                  &vBwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &qfield,
                Array< OneD, int >                                          &nonZeroIndex)
            {
                v_Diffuse_coeff(nConvectiveFields, fields, inarray, outarray, vFwd, vBwd,qfield,nonZeroIndex);
            }
            
            // Diffusion Calculate the physical derivatives
            SOLVER_UTILS_EXPORT void DiffuseCalculateDerivative(
                const int                                         nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
                const Array<OneD, Array<OneD, NekDouble>>         &pFwd =NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>         &pBwd =NullNekDoubleArrayofArray);

            // Diffusion Volume FLux
            SOLVER_UTILS_EXPORT void DiffuseVolumeFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array< OneD, int >                                  &nonZeroIndex       =   NullInt1DArray)
            {
                v_DiffuseVolumeFlux(nConvectiveFields, fields, inarray,inarrayderivative,VolumeFlux,nonZeroIndex);
            }

             // Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT void DiffuseTraceFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                    Array<OneD, Array<OneD, NekDouble> >            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd           =   NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd           =   NullNekDoubleArrayofArray,
                Array< OneD, int >                                  &nonZeroIndex   =   NullInt1DArray)
            {
                v_DiffuseTraceFlux(nConvectiveFields, fields, inarray,inarrayderivative,VolumeFlux,TraceFlux,pFwd, pBwd,nonZeroIndex);
            }

            // Diffusion term Trace Flux
            SOLVER_UTILS_EXPORT void DiffuseTraceFlux(
                const int                                                       nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
                const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qFwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qBwd,
                const Array<OneD, NekDouble>                                    &MuAVTrace,
                Array< OneD, int >                                              &nonZeroIndex   =   NullInt1DArray,
                const Array<OneD, Array<OneD, NekDouble>>                       &Aver           =   NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble>>                       &Jump           =   NullNekDoubleArrayofArray)
            {
                v_DiffuseTraceFlux(nConvectiveFields, fields, inarray,qfield,TraceFlux,pFwd, pBwd,qFwd,qBwd,MuAVTrace,nonZeroIndex,Aver,Jump);
            }

            SOLVER_UTILS_EXPORT void DiffuseTraceSymmFlux(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble>>                   &inarray,
                const Array<OneD,Array<OneD, Array<OneD, NekDouble> > >     &qfield,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &VolumeFlux,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &SymmFlux,
                const Array<OneD, Array<OneD, NekDouble>>                   &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                   &pBwd,
                Array< OneD, int >                                          &nonZeroIndex,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_Aver = NullNekDoubleArrayofArray,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_jump = NullNekDoubleArrayofArray)
            {
                v_DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,
                            VolumeFlux,SymmFlux,pFwd,pBwd,nonZeroIndex,solution_Aver,solution_jump);
            }

            SOLVER_UTILS_EXPORT void AddDiffusionSymmFluxToCoeff(
                const int                                           nConvectiveFields,
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
                const int                                           nConvectiveFields,
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

            SOLVER_UTILS_EXPORT void AddSymmFluxIntegralToOffDiag(
                const int                                                           nConvectiveFields,
                const int                                                           nDim,
                const int                                                           nPts,
                const int                                                           nTracePts,
                const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
                const Array<OneD, const int >                                       &nonZeroIndex,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                  &Fwdflux,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                  &Bwdflux,
                Array<OneD, Array<OneD, NekDouble> >                                &outarray);

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

            template<typename FuncPointerT, typename ObjectPointerT>
            void SetDiffusionFluxCons(FuncPointerT func, ObjectPointerT obj)
            {
                m_FunctorDiffusionfluxCons = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2,
                               std::placeholders::_3, std::placeholders::_4,
                               std::placeholders::_5, std::placeholders::_6,
                               std::placeholders::_7, std::placeholders::_8);
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
            void SetSpecialBndTreat(FuncPointerT func, ObjectPointerT obj)
            {
                m_SpecialBndTreat = std::bind(
                    func, obj, std::placeholders::_1, std::placeholders::_2);
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
                               std::placeholders::_5, std::placeholders::_6,
                               std::placeholders::_7);
            }

            void SetFunctorDerivBndCond(FunctorDerivBndCond DerivBndCond)
            {
                m_FunctorDerivBndCond = DerivBndCond;
            }

            inline void SetHomoDerivs(Array<OneD, Array<OneD, NekDouble> > &deriv)
            {
                v_SetHomoDerivs(deriv);
            }

            virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &GetFluxTensor()
            {
                return v_GetFluxTensor();
            }
            void GetAVmu(
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                      Array<OneD, NekDouble >                               &muvar,
                      Array<OneD, NekDouble >                               &MuVarTrace);

            void ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump)
            {
                v_ConsVarAveJump(nConvectiveFields,npnts,vFwd,vBwd,aver,jump);
            }

            const Array<OneD, const Array<OneD, NekDouble> > &GetTraceNormal()
            {
                return v_GetTraceNormal();
            }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
            void MinusVolumDerivJacToMat(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD,
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir,
                Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >              &gmtxarray)
            {
                v_MinusVolumDerivJacToMat(nConvectiveFields,pFields,ElmtJacArray,nDervDir,gmtxarray);
            }
            void MinusVolumDerivJacToMat( 
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD, 
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir, 
                Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >              &gmtxarray)
            {
                v_MinusVolumDerivJacToMat(nConvectiveFields,pFields,ElmtJacArray,nDervDir,gmtxarray);
            }
#endif

        protected:
            DiffusionFluxVecCB              m_fluxVector;
            DiffusionFluxVecCBNS            m_fluxVectorNS;
            DiffusionArtificialDiffusion    m_ArtificialDiffusionVector;
            DiffusionFluxCons               m_FunctorDiffusionfluxCons;
            FunctorDerivBndCond             m_FunctorDerivBndCond;
            SpecialBndTreat                 m_SpecialBndTreat;
            CalcViscosity                   m_CalcViscosity;
            DiffusionSymmFluxCons           m_FunctorSymmetricfluxCons;

            NekDouble                       m_time=0.0;

            virtual void v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
            {

            };

            virtual void v_Diffuse(
                const int                                         nConvectiveFields,
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

            virtual void v_Diffuse_coeff(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                Array<OneD, Array<OneD, NekDouble> >                        &outarray,
                const Array<OneD, Array<OneD, NekDouble> >                  &vFwd,
                const Array<OneD, Array<OneD, NekDouble> >                  &vBwd,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &qfield,
                Array< OneD, int >                                          &nonZeroIndex);
        
            virtual void v_ConsVarAveJump(
                const int                                           nConvectiveFields,
                const int                                           npnts,
                const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
                const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                      Array<OneD,       Array<OneD, NekDouble> >    &aver,
                      Array<OneD,       Array<OneD, NekDouble> >    &jump);

            // Diffusion Flux, calculate the physical derivatives
            virtual void v_DiffuseCalculateDerivative(
                const int nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > > &inarrayderivative,
                const Array<OneD, Array<OneD, NekDouble> >        &pFwd = NullNekDoubleArrayofArray,
                const Array<OneD, Array<OneD, NekDouble> >        &pBwd = NullNekDoubleArrayofArray);

            // Diffusion Volume Flux
            virtual void v_DiffuseVolumeFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array< OneD, int >                                  &nonZeroIndex) ;

             // Diffusion term Trace Flux
            virtual void v_DiffuseTraceFlux(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble>>           &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &inarrayderivative,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
                Array< OneD, int >                                  &nonZeroIndex);

            virtual void v_DiffuseTraceFlux(
                const int                                                       nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
                const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
                Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
                const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qFwd,
                const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qBwd,
                const Array<OneD, NekDouble>                                    &MuAVTrace,
                Array< OneD, int >                                              &nonZeroIndex  ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Aver          ,
                const Array<OneD, Array<OneD, NekDouble>>                       &Jump          );

            virtual void v_DiffuseTraceSymmFlux(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
                const Array<OneD, Array<OneD, NekDouble>>                   &inarray,
                const Array<OneD,Array<OneD, Array<OneD, NekDouble> > >     &qfield,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &VolumeFlux,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >          &SymmFlux,
                const Array<OneD, Array<OneD, NekDouble>>                   &pFwd,
                const Array<OneD, Array<OneD, NekDouble>>                   &pBwd,
                Array< OneD, int >                                          &nonZeroIndex,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_Aver,
                Array<OneD, Array<OneD, NekDouble> >                        &solution_jump)
            {
                ASSERTL0(false, "Not defined function v_DiffuseTraceSymmFlux.");
            }

            virtual void v_AddDiffusionSymmFluxToCoeff(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {

            }
            virtual void v_AddDiffusionSymmFluxToPhys(
                const int                                           nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                const Array<OneD, Array<OneD, NekDouble> >          &inarray,
                Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
                Array<OneD, Array<OneD, NekDouble> >                &outarray,
                const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
                const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
            {

            }

            virtual void v_SetHomoDerivs(
                Array<OneD, Array<OneD, NekDouble> > &deriv)
            {

            }

            virtual Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &v_GetFluxTensor()
            {
                static Array<OneD, Array<OneD, Array<OneD, NekDouble> > > tmp;
                return tmp;
            }

            virtual const Array<OneD, const Array<OneD, NekDouble> > &v_GetTraceNormal();

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
            virtual void v_MinusVolumDerivJacToMat(
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD,
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir,
                Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >              &gmtxarray);
            virtual void v_MinusVolumDerivJacToMat( 
                const int                                                   nConvectiveFields,
                const Array<OneD, MultiRegions::ExpListSharedPtr>           &pFields,
                const Array<OneD, const Array<OneD,  Array<OneD, 
                    Array<OneD,  Array<OneD,  NekDouble> > > > >            &ElmtJacArray,
                const int                                                   nDervDir, 
                Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >              &gmtxarray);
#endif


            /// Compute primary derivatives
            virtual void v_GetPrimVar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble>>         &inarray,
                  Array<OneD, Array<OneD, NekDouble>>         &primVar);

            /// Compute divergence and curl squared
            void GetDivCurl(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &pVarDer);

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
