///////////////////////////////////////////////////////////////////////////////
//
// File PreconditionerLinearWithDiag.cpp
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
// Description: Preconditioner definition for Diagonal with Linear Subspace
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerLinearWithDiag.h>
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

        string PreconditionerLinearWithDiag::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "FullLinearSpaceWithDiagonal",
                    PreconditionerLinearWithDiag::create,
                    "Full linear space and diagonal preconditioning");

       /**
         * @class PreconditionerLinearWithDiag
         *
         * This class implements preconditioning for the conjugate
	 * gradient matrix solver.
	 */

        PreconditionerLinearWithDiag::PreconditionerLinearWithDiag(
            const std::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }

        /**
         *
         */
        void PreconditionerLinearWithDiag::v_InitObject()
        {
            m_linSpacePrecon = GetPreconFactory().CreateInstance("FullLinearSpace",m_linsys.lock(),m_locToGloMap.lock());
            m_diagonalPrecon = GetPreconFactory().CreateInstance("Diagonal",m_linsys.lock(),m_locToGloMap.lock());
        }

        /**
         *
         */
        void PreconditionerLinearWithDiag::v_BuildPreconditioner()
        {
            m_linSpacePrecon->BuildPreconditioner();
            m_diagonalPrecon->BuildPreconditioner();
	}


        /**
         *
         */
        void PreconditionerLinearWithDiag::v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {

            Array<OneD, NekDouble> OutputDiag(pOutput.size());
            m_diagonalPrecon->DoPreconditioner(pInput, OutputDiag);

            // Since linear preconditioner just copies other entries
            // this will only modify the linear space degrees of
            // freedom
            m_linSpacePrecon->DoPreconditionerWithNonVertOutput(pInput, pOutput,OutputDiag);
        }

    }
}






