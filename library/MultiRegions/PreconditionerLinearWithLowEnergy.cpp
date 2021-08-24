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
#include <MultiRegions/PreconditionerLowEnergy.h>
#include <MultiRegions/PreconditionerLinearWithLowEnergy.h>
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
            const std::shared_ptr<GlobalLinSys> &plinsys,
            const AssemblyMapSharedPtr &pLocToGloMap)
            : Preconditioner(plinsys, pLocToGloMap)
        {
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_InitObject()
        {
            m_linSpacePrecon = GetPreconFactory().CreateInstance("FullLinearSpace",m_linsys.lock(),m_locToGloMap.lock());
            m_lowEnergyPrecon = GetPreconFactory().CreateInstance("LowEnergyBlock",m_linsys.lock(),m_locToGloMap.lock());

            //Set up multiplicity array for inverse transposed transformation matrix
            int nDirBnd     = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();
            int nGlobHomBnd = m_locToGloMap.lock()->GetNumGlobalBndCoeffs() - nDirBnd;
            int nLocBnd     = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            m_invMultiplicity = Array<OneD,NekDouble>(nGlobHomBnd,1.0);
            Array<OneD,NekDouble> loc(nLocBnd);

            // need to scatter from global array to handle sign changes
            m_locToGloMap.lock()->GlobalToLocalBnd(m_invMultiplicity, loc, nDirBnd);

            // Now assemble values back together to get multiplicity
            m_locToGloMap.lock()->AssembleBnd(loc,m_invMultiplicity, nDirBnd);
            Vmath::Sdiv(nGlobHomBnd,1.0,m_invMultiplicity,1,m_invMultiplicity,1);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformBasisToLowEnergy(
             Array<OneD, NekDouble>& pInOut)
        {
            m_lowEnergyPrecon->DoTransformBasisToLowEnergy(pInOut);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformCoeffsFromLowEnergy(
            Array<OneD, NekDouble>& pInput)
        {
            m_lowEnergyPrecon->DoTransformCoeffsFromLowEnergy(pInput);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformCoeffsToLowEnergy(
               const Array<OneD, NekDouble>& pInput,
               Array<OneD, NekDouble>& pOutput)
        {
            m_lowEnergyPrecon->DoTransformCoeffsToLowEnergy(pInput,pOutput);
        }

        /**
         *
         */
        void PreconditionerLinearWithLowEnergy::v_DoTransformBasisFromLowEnergy(
               const Array<OneD, NekDouble>& pInput,
               Array<OneD, NekDouble>& pOutput)
        {
            m_lowEnergyPrecon->DoTransformBasisFromLowEnergy(pInput,pOutput);
        }
        

        DNekScalMatSharedPtr PreconditionerLinearWithLowEnergy::
        v_TransformedSchurCompl(int n, int offset,
                                const std::shared_ptr<DNekScalMat > &loc_mat)
	{
            DNekScalMatSharedPtr returnval;
            returnval=m_lowEnergyPrecon->TransformedSchurCompl(n,offset, loc_mat);
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
            int nDirBndDofs     = m_locToGloMap.lock()->GetNumGlobalDirBndCoeffs();
            int nGlobHomBndDofs = m_locToGloMap.lock()->GetNumGlobalBndCoeffs() - nDirBndDofs;
            int nLocBndDofs     = m_locToGloMap.lock()->GetNumLocalBndCoeffs();

            Array<OneD, NekDouble> tmp(nGlobHomBndDofs, 0.0);

            Array<OneD, NekDouble> InputLinear(nGlobHomBndDofs);
            Array<OneD, NekDouble> OutputLowEnergy(nGlobHomBndDofs);
            Array<OneD, NekDouble> OutputLinear(nGlobHomBndDofs);
            Array<OneD, NekDouble> local(nLocBndDofs);


            //Transform input from low energy to original basis
            Vmath::Vmul(nGlobHomBndDofs,m_invMultiplicity,1,
                        pInput,1,OutputLinear,1);
            m_locToGloMap.lock()->GlobalToLocalBnd(OutputLinear,local,nDirBndDofs);
            m_lowEnergyPrecon->DoTransformBasisFromLowEnergy(local,local);
            m_locToGloMap.lock()->AssembleBnd(local,InputLinear,nDirBndDofs);

            //Apply linear space preconditioner
            m_linSpacePrecon->DoPreconditionerWithNonVertOutput
                (InputLinear, OutputLinear, tmp);

            // transform coefficients back to low energy space
            m_locToGloMap.lock()->GlobalToLocalBnd(OutputLinear,local,nDirBndDofs);
            m_lowEnergyPrecon->DoTransformCoeffsToLowEnergy(local,local);
            m_locToGloMap.lock()->LocalBndToGlobal(local,pOutput,nDirBndDofs,false);

            //Apply Low Energy preconditioner
            m_lowEnergyPrecon->DoPreconditioner(pInput, OutputLowEnergy);

            ASSERTL1(pOutput.size() >= nGlobHomBndDofs, "Output array is not correct");
            Vmath::Vadd(nGlobHomBndDofs,pOutput,1,OutputLowEnergy,1,pOutput,1);
        }

    }
}






