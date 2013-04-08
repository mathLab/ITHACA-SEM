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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VDmathArray.hpp>
#include <MultiRegions/PreconditionerLLE.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <LocalRegions/MatrixKey.h>
#include <math.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * Registers the class with the Factory.
         */

        string PreconditionerLLE::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "FullLinearSpaceWithLowEnergyBlock",
                    PreconditionerLLE::create,
                    "Full Linear space and low energy block preconditioning");
 
       /**
         * @class PreconditionerLLE
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */
        
        PreconditionerLLE::PreconditionerLLE(
            const boost::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }
       
        /**
         *
         */ 
        void PreconditionerLLE::v_InitObject()
        {
            m_linSpacePrecon = GetPreconFactory().CreateInstance("FullLinearSpace",m_linsys.lock(),m_locToGloMap);
            m_lowEnergyPrecon = GetPreconFactory().CreateInstance("LowEnergyBlock",m_linsys.lock(),m_locToGloMap);
        }

        /**
         *
         */
        void PreconditionerLLE::v_DoTransformToLowEnergy(
            Array<OneD, NekDouble>& pInOut,
            int offset)
        {
            m_lowEnergyPrecon->DoTransformToLowEnergy(pInOut,offset);
        }

        /**
         *
         */
        void PreconditionerLLE::v_DoTransformFromLowEnergy(
            Array<OneD, NekDouble>& pInput)
        {
            m_lowEnergyPrecon->DoTransformFromLowEnergy(pInput);
        }


        DNekScalBlkMatSharedPtr PreconditionerLLE::
        v_TransformedSchurCompl(int offset, const boost::shared_ptr<DNekScalBlkMat > &loc_mat)
	{
            DNekScalBlkMatSharedPtr returnval;
            returnval=m_lowEnergyPrecon->TransformedSchurCompl(offset,loc_mat);
            return returnval;
        }

        /**
         *
         */
        void PreconditionerLLE::v_BuildPreconditioner()
        {
            m_lowEnergyPrecon->BuildPreconditioner();
            m_linSpacePrecon->BuildPreconditioner();
	}


        /**
         *
         */
        void PreconditionerLLE::v_DoPreconditioner(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            Array<OneD, NekDouble> OutputLowEnergy(pOutput.num_elements());
            Array<OneD, NekDouble> OutputLinear(pOutput.num_elements());
            Array<OneD, NekDouble> OutputInverseMultiply(pOutput.num_elements());

            m_lowEnergyPrecon->DoPreconditioner(pInput, OutputLowEnergy);
            m_lowEnergyPrecon->DoMultiplybyInverseTransformationMatrix(OutputLowEnergy, OutputInverseMultiply);
            //m_linSpacePrecon->DoPreconditionerWithNonVertOutput(pInput, OutputLinear,OutputInverseMultiply);
            m_lowEnergyPrecon->DoTransformToLowEnergy(OutputInverseMultiply,pOutput);
        }

    }
}






