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
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::CommSharedPtr &vComm)
        {
            PreconCfsOpSharedPtr p = MemoryManager<PreconCfsBRJ>::
                            AllocateSharedPtr(pSession, vComm);
            return p;
        }

        ///Name of the class
        static std::string className;

        PreconCfsBRJ(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::CommSharedPtr &vComm);
        ~PreconCfsBRJ() {};

    protected:

        int         m_PreconItsStep;
        int         m_BRJRelaxParam;

        virtual void v_InitObject();

    private:

        virtual void v_DoPreconCfs(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag);

        virtual void v_BuildPreconCfs();
    };
}

#endif
