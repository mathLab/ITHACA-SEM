///////////////////////////////////////////////////////////////////////////////
//
// File PreconCfsBRJ.h
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
// Description: PreconCfsBRJ header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSBRJ
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFSBRJ

#include <CompressibleFlowSolver/Preconditioner/PreconCfsOp.h>

namespace Nektar
{
    /**
     * Block Relaxed(weighted) Jacobi iterative (BRJ) Preconditioner for CFS  
     * 
     */
    class PreconCfsBRJ : public PreconCfsOp
    {
    public:

        friend class MemoryManager<PreconCfsBRJ>;

        /// Creates an instance of this class
        static PreconCfsOpSharedPtr create(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::CommSharedPtr &vComm)
        {
            PreconCfsOpSharedPtr p = MemoryManager<PreconCfsBRJ>::
                            AllocateSharedPtr(pFields, pSession, vComm);
            return p;
        }

        ///Name of the class
        static std::string className;

        PreconCfsBRJ(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::CommSharedPtr &vComm);
        ~PreconCfsBRJ() {};

        virtual bool UpdatePreconMatCheck(
            const Array<OneD, const NekDouble>  &res,
            const NekDouble                     dtLambda);

    protected:

        int         m_PreconItsStep;
        int         m_BRJRelaxParam;

        Array<OneD, Array<OneD, DNekBlkMatSharedPtr>>   m_PreconMatVars;
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>>   m_PreconMatVarsSingle;
    
        Array<OneD, Array<OneD, NekDouble>>     m_PreconMatVarsOffDiag;
        DNekBlkMatSharedPtr                     m_PreconMat;
        SNekBlkMatSharedPtr                     m_PreconMatSingle;
        Array<OneD, DNekBlkMatSharedPtr >       m_TraceJac;
        TensorOfArray4D<NekDouble>              m_TraceJacArray;
        Array<OneD, DNekBlkMatSharedPtr >       m_TraceJacDeriv;
        TensorOfArray4D<NekDouble>              m_TraceJacDerivArray;
        Array<OneD, Array<OneD, NekDouble>>     m_TraceJacDerivSign;
        Array<OneD, SNekBlkMatSharedPtr >       m_TraceJacSingle;
        TensorOfArray4D<NekSingle>              m_TraceJacArraySingle;
        Array<OneD, SNekBlkMatSharedPtr >       m_TraceJacDerivSingle;
        TensorOfArray4D<NekSingle>              m_TraceJacDerivArraySingle;
        Array<OneD, Array<OneD, NekSingle> >    m_TraceJacDerivSignSingle;
        TensorOfArray5D<NekSingle>              m_TraceIPSymJacArraySingle;
        TensorOfArray5D<NekDouble>              m_TraceIPSymJacArray;

        TensorOfArray4D<NekDouble>          m_StdDMatDataDBB;
        TensorOfArray5D<NekDouble>          m_StdDMatDataDBDB;
        TensorOfArray4D<NekSingle>          m_StdSMatDataDBB;
        TensorOfArray5D<NekSingle>          m_StdSMatDataDBDB;

        PrecType                            m_PreconMatStorage;
    
        virtual void v_InitObject();

    private:

        virtual void v_DoPreconCfs(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag);

        virtual void v_BuildPreconCfs(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, const Array<OneD, NekDouble>>   &intmp,
            const NekDouble                                   time,
            const NekDouble                                   lambda);

        template<typename DataType, typename TypeNekBlkMatSharedPtr>
        void PreconBlkDiag(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, NekDouble>    &inarray,
                  Array<OneD, NekDouble>    &outarray,
            const TypeNekBlkMatSharedPtr    &PreconMatVars,
            const DataType                  &tmpDataType);

        template<typename DataType>
        void MinusOffDiag2Rhs(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const int                                      nvariables,
            const int                                      nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble>>&inarray,
            Array<OneD, Array<OneD, NekDouble>>            &outarray,
            bool                                           flagUpdateDervFlux,
            Array<OneD, Array<OneD, NekDouble>>            &FwdFluxDeriv,
            Array<OneD, Array<OneD, NekDouble>>            &BwdFluxDeriv,
            TensorOfArray3D<NekDouble>                     &qfield,
            TensorOfArray3D<NekDouble>                     &wspTrace,
            Array<OneD, Array<OneD, DataType>>             &wspTraceDataType,
            const TensorOfArray4D<DataType>                &TraceJacArray,
            const TensorOfArray4D<DataType>                &TraceJacDerivArray,
            const Array<OneD, const Array<OneD, DataType>> &TraceJacDerivSign,
            const TensorOfArray5D<DataType>                &TraceIPSymJacArray);

        template<typename TypeNekBlkMatSharedPtr>
        void AllocatePreconBlkDiagCoeff(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>>  &gmtxarray,
            const int                                         &nscale = 1);
        
        inline void AllocateNekBlkMatDig(
                  SNekBlkMatSharedPtr       &mat,
            const Array<OneD, unsigned int> nrow,
            const Array<OneD, unsigned int> ncol)
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

    };
}

#endif
