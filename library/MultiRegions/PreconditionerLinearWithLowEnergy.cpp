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
#include <MultiRegions/PreconditionerLinearWithLowEnergy.h>
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

        string PreconditionerLinearWithLowEnergy::className
                = GetPreconFactory().RegisterCreatorFunction(
                    "FullLinearSpaceWithLowEnergyBlock",
                    PreconditionerLinearWithLowEnergy::create,
                    "Full Linear space and low energy block preconditioning");
 
       /**
         * @class PreconditionerLinearWithLowEnergy
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */
        
        PreconditionerLinearWithLowEnergy::PreconditionerLinearWithLowEnergy(
            const boost::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }
       
        /**
         *
         */ 
        void PreconditionerLinearWithLowEnergy::v_InitObject()
        {
            m_linSpacePrecon = GetPreconFactory().CreateInstance("FullLinearSpace",m_linsys.lock(),m_locToGloMap);
            m_lowEnergyPrecon = GetPreconFactory().CreateInstance("LowEnergyBlock",m_linsys.lock(),m_locToGloMap);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformToLowEnergy(
            Array<OneD, NekDouble>& pInOut,
            int offset)
        {
            m_lowEnergyPrecon->DoTransformToLowEnergy(pInOut,offset);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformFromLowEnergy(
            Array<OneD, NekDouble>& pInput)
        {
            m_lowEnergyPrecon->DoTransformFromLowEnergy(pInput);
        }


        DNekScalMatSharedPtr PreconditionerLinearWithLowEnergy::
        v_TransformedSchurCompl(int offset, const boost::shared_ptr<DNekScalMat > &loc_mat)
	{
            DNekScalMatSharedPtr returnval;
            returnval=m_lowEnergyPrecon->TransformedSchurCompl(offset,loc_mat);
            return returnval;
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_BuildPreconditioner()
        {
            m_linSpacePrecon->BuildPreconditioner();
            m_lowEnergyPrecon->BuildPreconditioner();
	}


        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoPreconditioner(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            int nGlobal = pInput.num_elements();

            Array<OneD, NekDouble> OutputLowEnergy(nGlobal, 0.0);
            Array<OneD, NekDouble> tmp(nGlobal, 0.0);
            Array<OneD, NekDouble> OutputLinear(nGlobal, 0.0);
            Array<OneD, NekDouble> InputLinear(nGlobal, 0.0);

            //Apply Low Energy preconditioner
            m_lowEnergyPrecon->DoPreconditioner(pInput, OutputLowEnergy);

            //Transform input from low energy to original basis
            m_lowEnergyPrecon->DoMultiplybyInverseTransformationMatrix(pInput, InputLinear);

            //Apply linear space preconditioner
            m_linSpacePrecon->DoPreconditionerWithNonVertOutput(InputLinear, OutputLinear, tmp);

            m_lowEnergyPrecon->DoMultiplybyInverseTransposedTransformationMatrix(OutputLinear,pOutput);

            Vmath::Vadd(nGlobal,pOutput,1,OutputLowEnergy,1,pOutput,1);
        }

    }
}






