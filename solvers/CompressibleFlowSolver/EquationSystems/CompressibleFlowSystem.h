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

#include <boost/core/ignore_unused.hpp>

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
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>
#include <LibUtilities/LinearAlgebra/NekNonlinSys.h>

#define DEMO_IMPLICITSOLVER_JFNK_COEFF
#define CFS_DEBUGMODE
namespace Nektar
{
    /**
     *
     */
    class CompressibleFlowSystem: public SolverUtils::AdvectionSystem,
                                  public SolverUtils::FluidInterface
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

        /// Function to get estimate of min h/p factor per element
        Array<OneD, NekDouble>  GetElmtMinHP(void);

        void InitialiseNonlinSysSolver();

        void NonlinSysEvaluator1D(
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const bool                          &flag);

        void NonlinSysEvaluatorCoeff(
            TensorOfArray2D<NekDouble>  &inarray,
            TensorOfArray2D<NekDouble>  &out);
        void DoOdeRhsCoeff(
            const TensorOfArray2D<NekDouble>    &inarray,
            TensorOfArray2D<NekDouble>          &outarray,
            const NekDouble                     time);
        
        void DoAdvectionCoeff(
            const TensorOfArray2D<NekDouble>    &inarray,
            TensorOfArray2D<NekDouble>          &outarray,
            const NekDouble                     time,
            const TensorOfArray2D<NekDouble>    &pFwd,
            const TensorOfArray2D<NekDouble>    &pBwd);
        void DoImplicitSolvePhysToCoeff(
            const TensorOfArray2D<NekDouble>    &inpnts,
            TensorOfArray2D<NekDouble>          &outpnt,
            const NekDouble                     time,
            const NekDouble                     lambda);
        void DoImplicitSolveCoeff(
            const TensorOfArray2D<NekDouble>    &inpnts,
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const NekDouble                     time,
            const NekDouble                     lambda);
        bool UpdatePrecMatCheck(
            TensorOfArray1D<NekDouble>  &res);
        void CalPrecMat(
            const TensorOfArray2D<NekDouble>    &inpnts,
            const NekDouble                     time,
            const NekDouble                     lambda);
         void CalcRefValues(
            const TensorOfArray1D<NekDouble>    &inarray);
    
        virtual void GetPressure(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &pressure);

        virtual void GetDensity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, NekDouble>                     &density);

        virtual bool HasConstantDensity()
        {
            return false;
        }

        virtual void GetVelocity(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
                  Array<OneD, Array<OneD, NekDouble> >       &velocity);

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

        /// Store physical artificial viscosity
        Array<OneD, NekDouble>              m_muav;

        /// Store physical artificial viscosity
        Array<OneD, NekDouble>              m_muavTrace;

        // Parameters for local time-stepping
        bool                                m_useLocalTimeStep;

        bool                                m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT;
        bool                                m_DEBUG_VISCOUS_JAC_MAT;
        bool                                m_DEBUG_ADVECTION_JAC_MAT;
        bool                                m_updatePrecMatFlag;

        // 1: Adv; 2: Dif; Default: all
        int                                 m_DebugAdvDiffSwitch; 
       // 1: Vol; 2: Trace; Default: all
        int                                 m_DebugVolTraceSwitch; 
       // 1: Con; 2: Deriv; Default: all
        int                                 m_DebugConsDerivSwitch; 

        int                                 m_DebugInvMassSwitch   ; 
        int                                 m_DebugPlusSourceSwitch; 

        int                                 m_DebugIPSymmFluxJacSwitch; 

        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>
            m_StdDMatDataDBB;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble>>>>> m_StdDMatDataDBDB;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekSingle>>>> 
            m_StdSMatDataDBB;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekSingle>>>>> m_StdSMatDataDBDB;

        int                                m_nPadding = 1;

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
        NekDouble                           m_BndEvaluateTime;
        NekDouble                           m_TimeIntegLambda;
        TensorOfArray2D<NekDouble>          m_solutionPhys;

        NekDouble                           m_JacobiFreeEps;
        NekDouble                           m_inArrayNorm = -1.0;
        NekDouble                           m_NewtonAbsoluteIteTol;
        int                                 m_TotNewtonIts = 0;
        int                                 m_TotImpStages = 0;
        int                                 m_StagesPerStep = 0;
        Array<OneD, NekDouble>              m_magnitdEstimat;

        NekDouble                                   m_NewtonRelativeIteTol;
        /// cfl number for local time step(notice only for jfnk other see m_cflSafetyFactor)
        NekDouble                                   m_cflLocTimestep;
        /// In Jacobi iteration the SOR relaxation parameter
        NekDouble                                   m_SORRelaxParam;
        /// two strategies: time accurate or not.
        int                                        m_JFNKTimeAccurate;
        int                                         m_PrcdMatFreezNumb;
        /// preconditioning steps
        int                                         m_JFNKPrecondStep;
        int                                         m_MaxNonlinIte;

        LibUtilities::NekNonlinSysSharedPtr         m_nonlinsol;
        LibUtilities::NekSysOperators         m_NekSysOp;
        
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

        void preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);
 
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble>     &inarray,
            Array<OneD, NekDouble >          &outarray,
            const TypeNekBlkMatSharedPtr     &PrecMatVars,
            const DataType                   &tmpDataType);

        void preconditionerBlkSORCoeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &flag);

        template<typename DataType>
        void MinusOffDiag2Rhs(
            const int                                         nvariables,
            const int                                         nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble> >  &inarray,
            Array<OneD,       Array<OneD, NekDouble> >        &outarray,
            bool                                          flagUpdateDervFlux,
            Array<OneD,       Array<OneD, NekDouble> >        &FwdFluxDeriv,
            Array<OneD,       Array<OneD, NekDouble> >        &BwdFluxDeriv,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>  &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>  &wspTrace,
            Array<OneD, Array<OneD, DataType > >          &wspTraceDataType,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                 &TraceJacArray,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                 &TraceJacDerivArray,
            const Array<OneD, Array<OneD, DataType> >      &TraceJacDerivSign,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,
                      Array<OneD,DataType >>>>>     &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_volume(
            const Array<OneD, const Array<OneD, NekDouble> >   &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >                  &qfield,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >  &gmtxarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                 &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >   &StdMatDataDBDB);

        template<typename DataType>
        void CalcVolJacStdMat(
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                   &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >     &StdMatDataDBDB);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_boundary(
            const Array<OneD, const Array<OneD, NekDouble> >     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >    &gmtxarray,
            Array<OneD, TypeNekBlkMatSharedPtr >                 &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >            &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >             &TraceJacDerivSign,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>   &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void ElmtVarInvMtrx(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            TypeNekBlkMatSharedPtr                            &gmtVar,
            const DataType                                    &tmpDatatype);
    
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetTraceJac(
            const Array<OneD, const Array<OneD, NekDouble> >                        &inarray,
            TensorOfArray3D<NekDouble>          &qfield,
            Array<OneD, TypeNekBlkMatSharedPtr >                                    &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >                                    &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >                                     &TraceJacDerivSign,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>   &TraceIPSymJacArray);
       
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void NumCalRiemFluxJac(
            const int                                                               nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                       &fields,
            const Array<OneD, Array<OneD, NekDouble> >                              &AdvVel,
            const Array<OneD, Array<OneD, NekDouble> >                              &inarray,
            TensorOfArray3D<NekDouble>          &qfield,
            const NekDouble                                                         &time,
            const Array<OneD, Array<OneD, NekDouble> >                              &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >                              &Bwd,
            TypeNekBlkMatSharedPtr                                                  &FJac,
            TypeNekBlkMatSharedPtr                                                  &BJac,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>   &TraceIPSymJacArray);

        void PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TranSamesizeBlkDiagMatIntoArray(
            const TypeNekBlkMatSharedPtr                        &BlkMat,
            Array<OneD,Array<OneD,Array<OneD,DataType >>>       &MatArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TransTraceJacMatToArray(
            const Array<OneD, TypeNekBlkMatSharedPtr >          &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr >          &TraceJacDeriv,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>
                              &TraceJacArray,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>
                              &TraceJacDerivArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill2DArrayOfBlkDiagonalMat(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >   &gmtxarray,
            const DataType                                      valu);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill1DArrayOfBlkDiagonalMat( 
            Array<OneD, TypeNekBlkMatSharedPtr >    &gmtxarray,
            const DataType                          valu);

        void DoImplicitSolve_phy2Coeff(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                Array<OneD,       Array<OneD, NekDouble> >&out,
            const NekDouble time,
            const NekDouble lambda);

        void DoImplicitSolveCoeff(
            const Array<OneD, const Array<OneD, NekDouble> >&inpnts,
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                Array<OneD,       Array<OneD, NekDouble> >&out,
            const NekDouble time,
            const NekDouble lambda);

        template<typename TypeNekBlkMatSharedPtr>
        void AllocatePrecondBlkDiagCoeff(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const int                                          &nscale=1 );

        inline void AllocateNekBlkMatDig(
            SNekBlkMatSharedPtr                         &mat,
            const Array<OneD, unsigned int >            nrow,
            const Array<OneD, unsigned int >            ncol)
        {
            mat = MemoryManager<SNekBlkMat>
                ::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
            SNekMatSharedPtr loc_matNvar;
            for(int nelm = 0; nelm < nrow.size(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<SNekMat>::
                    AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                mat->SetBlock(nelm,nelm,loc_matNvar);
            }
        }

        inline void AllocateNekBlkMatDig(
            DNekBlkMatSharedPtr                         &mat,
            const Array<OneD, unsigned int >            nrow,
            const Array<OneD, unsigned int >            ncol)
        {
            mat = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrow, ncol, eDIAGONAL);
            DNekMatSharedPtr loc_matNvar;
            for(int nelm = 0; nelm < nrow.size(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<DNekMat>::
                    AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                mat->SetBlock(nelm,nelm,loc_matNvar);
            }
        }


        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetpreconditionerNSBlkDiagCoeff(
            const Array<OneD, const Array<OneD, NekDouble> >                                &inarray,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >                               &gmtxarray,
            TypeNekBlkMatSharedPtr                                                          &gmtVar,
            Array<OneD, TypeNekBlkMatSharedPtr >                                            &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >                                            &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >                                             &TraceJacDerivSign,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                       &TraceJacArray,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                       &TraceJacDerivArray,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>           &TraceIPSymJacArray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                 &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >   &StdMatDataDBDB);

        void MatrixMultiplyMatrixFreeCoeffCentral(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);

        void MatrixMultiplyMatrixFreeCoeffDualtimestep(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out,
            const  bool                   &controlFlag);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void MultiplyElmtInvMass_PlusSource(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const NekDouble                                    dtlamda,
            const DataType                                     tmpDataType);

        void GetFluxVectorJacDirctnElmt(
            const int                                           nConvectiveFields,
            const int                                           nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >          &locVars,
            const Array<OneD, NekDouble>                        &normals,
            DNekMatSharedPtr                                    &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                &PntJacArray);

        void GetFluxVectorJacPoint(
            const int                                   nConvectiveFields,
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
            DNekMatSharedPtr                            &fluxJac);
        
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
            TensorOfArray3D<NekDouble>      &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qFwd,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
                  Array<OneD, int >                                             &nonZeroIndex,
                  Array<OneD, Array<OneD, NekDouble> >                          &traceflux);
        
        void CalTraceIPSymmFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
            const NekDouble                                                     time,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
            Array<OneD, int >                                                   &nonZeroIndex,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                  &traceflux);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void CalVisFluxDerivJac(
            const int                                                       nConvectiveFields,
            const Array<OneD, Array<OneD, NekDouble> >                      &inarray,
            const Array<OneD, Array<OneD, NekDouble> >                      &Fwd,
            const Array<OneD, Array<OneD, NekDouble> >                      &Bwd,
            TypeNekBlkMatSharedPtr                                          &BJac,
            DataType                                                        &tmpDataType);

        void MinusDiffusionFluxJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, Array<OneD,  Array<OneD, NekDouble> > >       &locDerv,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, NekDouble>                                    &locDmuDT,
            const Array<OneD, NekDouble>                                    &normals,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
        {
            v_MinusDiffusionFluxJacDirctnElmt(nConvectiveFields,nElmtPnt,locVars,locDerv,locmu,locDmuDT,normals,wspMat,PntJacArray);
        }

        void GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray,
            const int                                                       nfluxDir)
        {
            v_GetFluxDerivJacDirctn(explist,normals,nDervDir,inarray,ElmtJacArray,nfluxDir);
        }

        void GetFluxDerivJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const int                                                       nDervDir,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, Array<OneD, NekDouble> >                      &locnormal,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
        {
            v_GetFluxDerivJacDirctnElmt(nConvectiveFields,nElmtPnt,nDervDir,locVars,locmu,locnormal,wspMat,PntJacArray);
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

        void CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
        {
            v_CalphysDeriv(inarray, qfield);
        }
        void GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac)
        {
            v_GetDiffusionFluxJacPoint(conservVar,conseDeriv,mu,DmuDT,normals,fluxJac);
        }

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
        void DoDiffusionCoeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> >             &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd);
        void MatrixMultiplyMatrixFreeCoeff(
            const  TensorOfArray1D<NekDouble>   &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const bool                          &flag = false);

        void GetFluxVector(
            const Array<OneD, Array<OneD, NekDouble> >       &physfield,
            TensorOfArray3D<NekDouble>                       &flux);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >       &physfield,
            TensorOfArray3D<NekDouble>                       &flux);

        void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time);

        void SetBoundaryConditionsBwdWeight();

        void GetElmtTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD, NekDouble> &tstep);

        void CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
            Array<OneD, NekDouble>                              &mu,
            Array<OneD, NekDouble>                              &DmuDT)
        {
            v_CalcMuDmuDT(inarray,mu,DmuDT);
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
            boost::ignore_unused(inarray, outarray, pFwd, pBwd);
            // Do nothing by default
        }
        
        virtual void v_DoDiffusionCoeff(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
        {
            boost::ignore_unused(inarray, outarray, pFwd, pBwd);
        }


        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity(
            const NekDouble SpeedSoundFactor);

        virtual void v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2);
        virtual void v_CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, NekDouble>                           &mu,
            Array<OneD, NekDouble>                           &DmuDT)
        {
            boost::ignore_unused(inarray, mu, DmuDT);
        }
                
        virtual void v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac);

        virtual void v_CalphysDeriv(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD,       Array<OneD, Array<OneD, NekDouble> > >  &qfield)
        {
            boost::ignore_unused(inarray, qfield);
        }

        virtual void v_MinusDiffusionFluxJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, Array<OneD,  Array<OneD, NekDouble> > >       &locDerv,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, NekDouble>                                    &locDmuDT,
            const Array<OneD, NekDouble>                                    &normals,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray);

        virtual void v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray,
            const int                                                       nfluxDir);

        virtual void v_GetFluxDerivJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const int                                                       nDervDir,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, Array<OneD, NekDouble> >                      &locnormal,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray);

        virtual void v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac);

    };
}
#endif
