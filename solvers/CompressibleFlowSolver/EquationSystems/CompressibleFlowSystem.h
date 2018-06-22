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

#ifdef DEMO_IMPLICITSOLVER_JFNK 
        void DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&forc,
                  Array<OneD,       Array<OneD, NekDouble> >&sol,
            const NekDouble time,
            const NekDouble lambda);

        void NonlinSysEvaluator(
                  Array<OneD,       Array<OneD, NekDouble> >&inarray,
                  Array<OneD,       Array<OneD, NekDouble> >&out);

        void MatrixMultiply(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);

        void MatrixMultiply_MatrixFree(
            const  Array<OneD, NekDouble> &inarray,
                   Array<OneD, NekDouble >&out);


        void AllocatePrecondBlkDiag(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

        void GetpreconditionerNSBlkDiag(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

        void MultiplyElmtBwdInvMass(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda);

        void MultiplyElmtBwdInvMassFwd(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda);
        
        void DebugNumCalElmtJac(
            Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars,
            const int nelmt);
        void DebugNumCalJac(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);
        
#endif
        
        
        void preconditioner(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&out);
        void preconditioner_BlkDiag(
            const Array<OneD, NekDouble> &inarray,
            Array<OneD, NekDouble >&outarray);

        void AddMatNSBlkDiag_volume(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

        void AddMatNSBlkDiag_boundary(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

        void ElmtVarInvMtrx(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);
        
        Array<OneD, DNekBlkMatSharedPtr> GetTraceJac(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray);

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

        void Cout2DArrayBlkMat(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
            const unsigned int nwidthcolm=12);

        void Fill2DArrayOfBlkDiagonalMat(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
            const NekDouble valu);

#define DEMO_IMPLICITSOLVER_JFNK_COEFF
#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

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
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray);

        void MatrixMultiply_MatrixFree_coeff(
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

        void NonlinSysEvaluator_coeff_bnd(
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

        void MultiplyElmtInvMass_PlusSource(
            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda);
        

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
        
        Array<OneD, Array<OneD, DNekMatSharedPtr> > 
            GetFluxVectorJacDirctn(const int nDirctn,
                         const Array<OneD, const Array<OneD, NekDouble> >&inarray);
        void GetFluxVectorJacPoint(
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
                 DNekMatSharedPtr                       &fluxJac);
        void GetFluxVectorDeAlias(
            const Array<OneD, Array<OneD, NekDouble> >         &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

        void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time);

        void GetElmtTimeStep(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD, NekDouble> &tstep);

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

        virtual Array<OneD, NekDouble> v_GetMaxStdVelocity();
    };
}
#endif
