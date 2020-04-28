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
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/Expansion2D.h>

#define CFS_DEBUGMODE
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

        /// Function to get estimate of min h/p factor per element
        Array<OneD, NekDouble>  GetElmtMinHP(void);
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

        NekDouble                           m_JFEps;
        
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

#ifdef CFS_DEBUGMODE
       // 1: Adv; 2: Dif; Default: all
        int                                 m_DebugAdvDiffSwitch; 
       // 1: Vol; 2: Trace; Default: all
        int                                 m_DebugVolTraceSwitch; 
       // 1: Con; 2: Deriv; Default: all
        int                                 m_DebugConsDerivSwitch; 


        int                                 m_DebugNumJacMatSwitch;
        int                                 m_DebugOutputJacMatSwitch;

        int                                 m_DebugInvMassSwitch   ; 
        int                                 m_DebugPlusSourceSwitch; 

        int                                 m_DebugIPSymmFluxJacSwitch; 
        int                                 m_DebugNumJacBSOR;
#endif

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >                m_StdDMatDataDBB;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > >  m_StdDMatDataDBDB;

        Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekSingle> > > >                m_StdSMatDataDBB;
        Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekSingle> > > > >  m_StdSMatDataDBDB;
        int                                 m_nPadding = 1;
#endif

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
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble>                                &inarray,
            Array<OneD, NekDouble >                                     &outarray,
            const Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >     &PrecMatVars,
            const DataType                                              &tmpDataType);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble>                                &inarray,
            Array<OneD, NekDouble >                                     &outarray,
            const TypeNekBlkMatSharedPtr                                &PrecMatVars,
            const DataType                                              &tmpDataType);

        void preconditioner_NumJac(
            const Array<OneD, NekDouble>                                                &inarray,
            Array<OneD, NekDouble >                                                     &outarray,
            const Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                        &PrecMatVars,
            const Array<OneD, Array<OneD, NekDouble > >                                 &PrecMatVarsOffDiag);
        void MinusOffDiag2RhsNumJac(
            const int                                                                   nvariables,
            const int                                                                   nCoeffs,
            const Array<OneD, NekDouble>                                                &inarray,
            Array<OneD, NekDouble>                                                      &outarray,
            const Array<OneD, Array<OneD, NekDouble > >                                 &PrecMatVarsOffDiag);
            
        void preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &flag);

        // void MinusOffDiag2Rhs(
        //     const int nvariables,
        //     const int nCoeffs,
        //     const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
        //           Array<OneD,       Array<OneD, NekDouble> >    &outarray,
        //     bool                                                flagUpdateDervFlux,
        //           Array<OneD,       Array<OneD, NekDouble> >    &FwdFluxDeriv,
        //           Array<OneD,       Array<OneD, NekDouble> >    &BwdFluxDeriv,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &qfield,
        //     Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &tmpTrace);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void MinusOffDiag2Rhs(
            const int                                               nvariables,
            const int                                               nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble> >        &inarray,
            Array<OneD,       Array<OneD, NekDouble> >              &outarray,
            bool                                                    flagUpdateDervFlux,
            Array<OneD,       Array<OneD, NekDouble> >              &FwdFluxDeriv,
            Array<OneD,       Array<OneD, NekDouble> >              &BwdFluxDeriv,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >      &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >      &tmpTrace,
            const Array<OneD, TypeNekBlkMatSharedPtr >              &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr >              &TraceJacDeriv,
            const Array<OneD, Array<OneD, DataType> >               &TraceJacDerivSign);

        template<typename DataType>
        void MinusOffDiag2Rhs(
            const int                                                                       nvariables,
            const int                                                                       nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble> >                                &inarray,
            Array<OneD,       Array<OneD, NekDouble> >                                      &outarray,
            bool                                                                            flagUpdateDervFlux,
            Array<OneD,       Array<OneD, NekDouble> >                                      &FwdFluxDeriv,
            Array<OneD,       Array<OneD, NekDouble> >                                      &BwdFluxDeriv,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                              &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                              &wspTrace,
            Array<OneD, Array<OneD, DataType > >                                            &wspTraceDataType,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                 &TraceJacArray,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>                 &TraceJacDerivArray,
            const Array<OneD, Array<OneD, DataType> >                                       &TraceJacDerivSign,
            const Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>     &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_volume(
            const Array<OneD, const Array<OneD, NekDouble> >                                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >                  &qfield,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >                               &gmtxarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                 &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >   &StdMatDataDBDB);

        template<typename DataType>
        void CalcVolJacStdMat(
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                   &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >     &StdMatDataDBDB);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void AddMatNSBlkDiag_boundary(
            const Array<OneD, const Array<OneD, NekDouble> >                        &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >                      &qfield,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >                       &gmtxarray,
            Array<OneD, TypeNekBlkMatSharedPtr >                                    &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >                                    &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >                                     &TraceJacDerivSign,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>   &TraceIPSymJacArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void ElmtVarInvMtrx(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            TypeNekBlkMatSharedPtr                            &gmtVar,
            const DataType                                    &tmpDatatype);
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void ElmtVarInvMtrx(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const DataType                                    &tmpDatatype);
        
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetTraceJac(
            const Array<OneD, const Array<OneD, NekDouble> >                        &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >          &qfield,
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
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >          &qfield,
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

#ifdef CFS_DEBUGMODE

        void NumJacElemental(
            DNekMatSharedPtr    &NumericalJacobianMatrix,
            const int                 RowElementID,
            const int                 ColElementID);
        
        void CalOffDiagJacByMinusOffDiagElemental(
            DNekMatSharedPtr    &MinusoffJacobianMatrix,
            const int                 RowElementID,
            const int                 ColElementID);
        void DebugCheckJac(
            const int                 RowElementID,
            const int                 ColElementID);
        
        void DebugNumCalJac_coeff(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                      &gmtxarray,
            Array<OneD, Array<OneD, NekDouble > >                               &JacOffDiagArray =NullNekDoubleArrayofArray);
            
        void DebugNumCalElmtJac_coeff(
            Array<OneD, Array<OneD, DNekMatSharedPtr> >                         &ElmtPrecMatVars ,
            const int                                                           nelmt,
            Array<OneD, Array<OneD, NekDouble > >                               &JacOffDiagArray);
        void DebugNumCalElmtJac_coeff(
            Array<OneD, Array<OneD, DNekMatSharedPtr> >                         &ElmtPrecMatVars ,
            const int                                                           nelmt);

        void NonlinSysEvaluator_coeff_out(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

        template<typename TypeNekBlkMatSharedPtr>
        void CoutBlkMat(
            TypeNekBlkMatSharedPtr &gmtx, 
            const unsigned int nwidthcolm=12);

        template<typename TypeNekBlkMatSharedPtr>
        void CoutStandardMat(
            TypeNekBlkMatSharedPtr &loc_matNvar,
            const unsigned int nwidthcolm=12);

        template<typename TypeNekBlkMatSharedPtr>
        void Cout1DArrayBlkMat(
            Array<OneD, TypeNekBlkMatSharedPtr> &gmtxarray,
            const unsigned int nwidthcolm=12);
        
        template<typename TypeNekBlkMatSharedPtr>
        void Cout2DArrayBlkMat(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const unsigned int nwidthcolm=12);

        template<typename TypeNekBlkMatSharedPtr>
        void Cout2DArrayStdMat(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const unsigned int nwidthcolm=12);
#endif
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TranSamesizeBlkDiagMatIntoArray(
            const TypeNekBlkMatSharedPtr                        &BlkMat,
            Array<OneD,Array<OneD,Array<OneD,DataType >>>       &MatArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void TransTraceJacMatToArray(
            const Array<OneD, TypeNekBlkMatSharedPtr >                      &TraceJac,
            const Array<OneD, TypeNekBlkMatSharedPtr >                      &TraceJacDeriv,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>       &TraceJacArray,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>       &TraceJacDerivArray);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill2DArrayOfBlkDiagonalMat(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >   &gmtxarray,
            const DataType                                      valu);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void Fill1DArrayOfBlkDiagonalMat( 
            Array<OneD, TypeNekBlkMatSharedPtr >    &gmtxarray,
            const DataType                          valu);

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

        template<typename TypeNekBlkMatSharedPtr>
        void AllocatePrecondBlkDiag_coeff(
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
            for(int nelm = 0; nelm < nrow.num_elements(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<SNekMat>::AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
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
            for(int nelm = 0; nelm < nrow.num_elements(); ++nelm)
            {
                int nrowsVars = nrow[nelm];
                int ncolsVars = ncol[nelm];
                
                loc_matNvar = MemoryManager<DNekMat>::AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                mat->SetBlock(nelm,nelm,loc_matNvar);
            }
        }


        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetpreconditionerNSBlkDiag_coeff(
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
        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void GetpreconditionerNSBlkDiag_coeff(
            const Array<OneD, const Array<OneD, NekDouble> >                                &inarray,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >                               &gmtxarray,
            Array<OneD, TypeNekBlkMatSharedPtr >                                            &TraceJac,
            Array<OneD, TypeNekBlkMatSharedPtr >                                            &TraceJacDeriv,
            Array<OneD, Array<OneD, DataType> >                                             &TraceJacDerivSign,
            Array<OneD,Array<OneD,Array<OneD,Array<OneD,Array<OneD,DataType >>>>>           &TraceIPSymJacArray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > >                 &StdMatDataDBB,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, DataType> > > > >   &StdMatDataDBDB);

        void MatrixMultiply_MatrixFree_coeff(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);
        void MatrixMultiply_MatrixFree_coeff_central(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);
        void MatrixMultiply_MatrixFree_coeff_FourthCentral(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out);

        void MatrixMultiply_MatrixFree_coeff_dualtimestep(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out,
            const  bool                   &controlFlag);

        void NonlinSysEvaluator_coeff(
                Array<OneD, Array<OneD, NekDouble> > &inarray,
                Array<OneD, Array<OneD, NekDouble> > &out);

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

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void MultiplyElmtInvMass_PlusSource(
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> > &gmtxarray,
            const NekDouble                                    dtlamda,
            const DataType                                     tmpDataType);

        void GetFluxVectorJacDirctn(
            const int                                           nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >    &inarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray);

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
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
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

        void MinusDiffusionFluxJacDirctn(
            const int                                                       nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray)
        {
            v_MinusDiffusionFluxJacDirctn(nDirctn,inarray, qfields,ElmtJacArray);
        }

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
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac)
        {
            v_GetDiffusionFluxJacPoint(conservVar,conseDeriv,mu,DmuDT,normals,fluxJac);
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

        void CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, NekDouble>                                          &mu,
            Array<OneD, NekDouble>                                          &DmuDT)
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
            // Do nothing by default
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

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity(const NekDouble SpeedSoundFactor=1.0);

        virtual void v_GetViscousSymmtrFluxConservVar(
            const int                                                       nConvectiveFields,
            const int                                                       nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
            Array< OneD, int >                                              &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> >                      &normals);
        
        virtual void v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2);
        virtual void v_CalcMuDmuDT(
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, NekDouble>                                          &mu,
            Array<OneD, NekDouble>                                          &DmuDT)
        {
        }
                
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
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
        {}

        virtual void v_MinusDiffusionFluxJacDirctn(
            const int                                                       nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray);
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
