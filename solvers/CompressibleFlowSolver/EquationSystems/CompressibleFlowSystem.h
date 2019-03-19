///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.h
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
// Description: Auxiliary functions for the compressible flow system
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_EQUATIONSYSTEMS_COMPRESSIBLEFLOWSYSTEM_H

#include <CompressibleFlowSolver/ArtificialDiffusion/ArtificialDiffusion.h>
#include <CompressibleFlowSolver/Misc/VariableConverter.h>
#include <CompressibleFlowSolver/BoundaryConditions/CFSBndCond.h>
#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    /**
     *
     */
    class CompressibleFlowSystem: public SolverUtils::AdvectionSystem
    {
    public:

        friend class MemoryManager<CompressibleFlowSystem>;

        virtual ~CompressibleFlowSystem();

        /// Function to calculate the stability limit for DG/CG.
        NekDouble GetStabilityLimit(int n);

        /// Function to calculate the stability limit for DG/CG
        /// (a vector of them).
        Array<OneD, NekDouble> GetStabilityLimitVector(
            const Array<OneD,int> &ExpOrder);

    protected:
        SolverUtils::DiffusionSharedPtr     m_diffusion;
        ArtificialDiffusionSharedPtr        m_artificialDiffusion;
        Array<OneD, Array<OneD, NekDouble> >m_vecLocs;
        NekDouble                           m_gamma;
        std::string                         m_shockCaptureType;

        // Parameters for exponential filtering
        NekDouble                           m_filterAlpha;
        NekDouble                           m_filterExponent;
        NekDouble                           m_filterCutoff;
        bool                                m_useFiltering;

        // Parameters for local time-stepping
        bool                                m_useLocalTimeStep;

        // Auxiliary object to convert variables
        VariableConverterSharedPtr          m_varConv;

        // User defined boundary conditions
        std::vector<CFSBndCondSharedPtr>    m_bndConds;

        // Forcing term
        std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

        enum PreconditionerType
        {
            eNull,    ///< No Solution type specified
            eDiagonal,
            eSparse,
        };

        PreconditionerType                  m_PrecMatStorage;

        CompressibleFlowSystem(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr& pGraph);

        virtual void v_InitObject();

        void InitialiseParameters();

        void InitAdvection();

        void DoOdeRhs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time);
        void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time);

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble> &inarray,
            Array<OneD, NekDouble >&outarray);

        void preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray);
        
        void MinusOffDiag2Rhs(
            const int nvariables,
            const int nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
                  Array<OneD,       Array<OneD, NekDouble> >    &outarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &tmpTrace);

        void AddMatNSBlkDiag_volume(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                  &gmtxarray);

        void AddMatNSBlkDiag_boundary(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                  &gmtxarray,
            Array<OneD, DNekBlkMatSharedPtr >                               &TraceJac,
            Array<OneD, DNekBlkMatSharedPtr >                               &TraceJacDeriv,
            Array<OneD, Array<OneD, NekDouble> >                            &TraceJacDerivSign);

        void ElmtVarInvMtrx(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);
        
        void GetTraceJac(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, DNekBlkMatSharedPtr >                               &TraceJac,
            Array<OneD, DNekBlkMatSharedPtr >                               &TraceJacDeriv,
            Array<OneD, Array<OneD, NekDouble> >                            &TraceJacDerivSign);

        void NumCalRiemFluxJac(
            const int                                                       nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
            const Array<OneD, Array<OneD, NekDouble> >                      &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            const NekDouble                                                 &time,
            const Array<OneD, Array<OneD, NekDouble> >                      &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >                      &Bwd,
            DNekBlkMatSharedPtr                                             &FJac,
            DNekBlkMatSharedPtr                                             &BJac);

        void PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw);

        void CoutBlkMat(
            DNekBlkMatSharedPtr &gmtx, 
            const unsigned int nwidthcolm=12);

        void CoutStandardMat(
            DNekMatSharedPtr &loc_matNvar,
            const unsigned int nwidthcolm=12);

        void Cout1DArrayBlkMat(
            Array<OneD, DNekBlkMatSharedPtr> &gmtxarray,
            const unsigned int nwidthcolm=12);
        
        void Cout2DArrayBlkMat(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
            const unsigned int nwidthcolm=12);

        void Cout2DArrayStdMat(
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &gmtxarray,
            const unsigned int nwidthcolm=12);

        void Fill2DArrayOfBlkDiagonalMat(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
            const NekDouble valu);

        void DoImplicitSolve_phy2coeff(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                Array<OneD,       Array<OneD, NekDouble> >&out,
            const NekDouble time,
            const NekDouble lambda);

        void DoImplicitSolve_coeff(
            const Array<OneD, const Array<OneD, NekDouble> >&inpnts,
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                Array<OneD,       Array<OneD, NekDouble> >&out,
            const NekDouble time,
            const NekDouble lambda);

        void AllocatePrecondBlkDiag_coeff(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

         void GetpreconditionerNSBlkDiag_coeff(
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
            Array<OneD, DNekBlkMatSharedPtr >                   &TraceJac,
            Array<OneD, DNekBlkMatSharedPtr >                   &TraceJacDeriv,
            Array<OneD, Array<OneD, NekDouble> >                &TraceJacDerivSign);

        void MatrixMultiply_MatrixFree_coeff(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);

        void MatrixMultiply_MatrixFree_coeff_dualtimestep(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);

        void DebugNumCalJac_coeff(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);
            
        void DebugNumCalElmtJac_coeff(
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars,
            const int nelmt);
        
        void NonlinSysEvaluator_coeff(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

        void NonlinSysEvaluator_coeff_out(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

        // void NonlinSysEvaluator_coeff_bnd(
        //         Array<OneD, Array<OneD, NekDouble> > &inarray,
        //         Array<OneD, Array<OneD, NekDouble> > &out);

        void DoOdeRhs_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time);

        
        void DoAdvection_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        void MultiplyElmtInvMass_PlusSource(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda);

        void GetFluxVectorJacDirctn(
            const int                                           nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >   &ElmtJac);
        void GetFluxVectorJacPoint(
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
                 DNekMatSharedPtr                       &fluxJac);
        
        void CalTraceNumericalFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
            const NekDouble                                                     time,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
                  Array<OneD, int >                                             &nonZeroIndex,
                  Array<OneD, Array<OneD, NekDouble> >                          &traceflux);

        void CalVisFluxDerivJac(
            const int                                                       nConvectiveFields,
            const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
            const Array<OneD, Array<OneD, NekDouble> >                      &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >                      &Bwd,
            DNekBlkMatSharedPtr                                             &BJac);

        void MinusDiffusionFluxJacDirctn(
            const int                                                       nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
        {
            v_MinusDiffusionFluxJacDirctn(nDirctn,inarray, qfields,ElmtJac);
        }

        void GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
        {
            v_GetFluxDerivJacDirctn(explist,normals,nDervDir,inarray,ElmtJac);
        }
        // void GetFluxDerivJacDirctn(
        //     const MultiRegions::ExpListSharedPtr                            &explist,
        //     const int                                                       nFluxDir,
        //     const int                                                       nDervDir,
        //     const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        //           Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
        // {
        //     v_GetFluxDerivJacDirctn(explist,nFluxDir,nDervDir,inarray,ElmtJac);
        // }
        void CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
        {
            v_CalphysDeriv(inarray, qfield);
        }
        void GetDiffusionFluxJacDirctn(
            const int nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);
        void GetDiffusionFluxJacPoint(
            const int                                           nelmt,
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals, 
                  DNekMatSharedPtr                              &fluxJac)
        {
            v_GetDiffusionFluxJacPoint(nelmt,conservVar,conseDeriv,mu,DmuDT,normals,fluxJac);
        }
#endif

        void DoAdvection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble                                   time,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);

        void DoDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);
        void DoDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd);

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >               &physfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time);
        
        void SetBoundaryConditionsBwdWeight();

        void SetBoundaryConditionsDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                    &physarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &dervarray,
            NekDouble                                                           time,
            const Array<OneD, const Array<OneD, NekDouble> >                    &pFwd       = NullNekDoubleArrayofArray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &pDervFwd   = NullNekDoubleArrayofArrayofArray);

        void GetElmtTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD, NekDouble> &tstep);

        void GetViscousSymmtrFluxConservVar(
            const int                                                       nConvectiveFields,
            const int                                                       nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
            Array< OneD, int >                                              &nonZeroIndex,    
            const Array<OneD, Array<OneD, NekDouble> >                      &normals)
        {
            v_GetViscousSymmtrFluxConservVar(nConvectiveFields,nSpaceDim,inaverg,inarray,outarray,nonZeroIndex,normals);
        }

        virtual NekDouble v_GetTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray);
        virtual void v_SetInitialConditions(
            NekDouble initialtime           = 0.0,
            bool      dumpInitialConditions = true,
            const int domain                = 0);

        NekDouble GetGamma()
        {
            return m_gamma;
        }

        const Array<OneD, const Array<OneD, NekDouble> > &GetVecLocs()
        {
            return m_vecLocs;
        }

        const Array<OneD, const Array<OneD, NekDouble> > &GetNormals()
        {
            return m_traceNormals;
        }

        virtual void v_ExtraFldOutput(
            std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
            std::vector<std::string>             &variables);

        virtual void v_DoDiffusion(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
        {
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
            }
        }

        virtual void v_DoDiffusion_coeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
        {
            // Do nothing by default
        }

        virtual void v_DoDiffusionFlux(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble>>              &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
        {
            //Artificial Diffusion need to implement
            if (m_shockCaptureType != "Off")
            {
                m_artificialDiffusion->DoArtificialDiffusionFlux(inarray, VolumeFlux,TraceFlux);
            }
        }

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity();

        virtual void v_GetViscousSymmtrFluxConservVar(
            const int                                                       nConvectiveFields,
            const int                                                       nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
            Array< OneD, int >                                              &nonZeroIndex,    
            const Array<OneD, Array<OneD, NekDouble> >                      &normals);
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        virtual void v_GetDiffusionFluxJacPoint(
            const int                                           nelmt,
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals, 
                  DNekMatSharedPtr                              &fluxJac);
        virtual void v_CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
        {}

        virtual void v_MinusDiffusionFluxJacDirctn(
            const int                                                       nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);
        virtual void v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);

        // virtual void v_GetFluxDerivJacDirctn(
        //     const MultiRegions::ExpListSharedPtr                            &explist,
        //     const int                                                       nFluxDir,
        //     const int                                                       nDervDir,
        //     const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
        //           Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);
#endif
    };
}
#endif
