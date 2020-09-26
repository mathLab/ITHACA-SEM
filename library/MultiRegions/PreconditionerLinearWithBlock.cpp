///////////////////////////////////////////////////////////////////////////////
//
// File Preconditioner.cpp
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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerLinearWithBlock.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <LocalRegions/MatrixKey.h>
#include <cmath>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */

        string PreconditionerLinearWithBlock::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "FullLinearSpaceWithBlock",
                    PreconditionerLinearWithBlock::create,
                    "Full Linear space and block preconditioning");

       /**
         * @class PreconditionerLinearWithBlock
         *
         * This class implements preconditioning for the conjugate
	 * gradient matrix solver.
	 */

        PreconditionerLinearWithBlock::PreconditionerLinearWithBlock(
            const std::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }

        /**
         *
         */
        void PreconditionerLinearWithBlock::v_InitObject()
        {
            m_linSpacePrecon = GetPreconFactory().CreateInstance("FullLinearSpace",m_linsys.lock(),m_locToGloMap.lock());
            m_blockPrecon = GetPreconFactory().CreateInstance("Block",m_linsys.lock(),m_locToGloMap.lock());
        }

        /**
         *
         */
        void PreconditionerLinearWithBlock::v_BuildPreconditioner()
        {
            m_linSpacePrecon->BuildPreconditioner();
            m_blockPrecon->BuildPreconditioner();
	}


        /**
         *
         */
        void PreconditionerLinearWithBlock::v_DoPreconditioner(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            int nGlobal = pInput.size();

            Array<OneD, NekDouble> OutputBlock(nGlobal, 0.0);
            Array<OneD, NekDouble> OutputLinear(nGlobal, 0.0);
            Array<OneD, NekDouble> InputLinear(nGlobal, 0.0);

            //Apply Low Energy preconditioner
            m_blockPrecon->DoPreconditioner(pInput, OutputBlock);

            //Apply linear space preconditioner
            m_linSpacePrecon->DoPreconditionerWithNonVertOutput(pInput, pOutput,OutputBlock);
        }

    }
}






