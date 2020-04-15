///////////////////////////////////////////////////////////////////////////////
//
// File Preconditioner.h
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
// Description: Preconditioner header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLINEARWITHLOWENERGY_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLINEARWITHLOWENERGY_H
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class PreconditionerLinearWithLowEnergy;
        typedef std::shared_ptr<PreconditionerLinearWithLowEnergy>  PreconditionerLinearWithLowEnergySharedPtr;

        class PreconditionerLinearWithLowEnergy: public Preconditioner
	{
        public:
            /// Creates an instance of this class
            static PreconditionerSharedPtr create(
                        const std::shared_ptr<GlobalLinSys> &plinsys,
                        const std::shared_ptr<AssemblyMap> &pLocToGloMap)
            {
	        PreconditionerSharedPtr p = MemoryManager<PreconditionerLinearWithLowEnergy>::AllocateSharedPtr(plinsys,pLocToGloMap);
	        p->InitObject();
	        return p;
            }

            /// Name of class
            static std::string className;

            MULTI_REGIONS_EXPORT PreconditionerLinearWithLowEnergy(
                         const std::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~PreconditionerLinearWithLowEnergy() {}

	protected:
            PreconditionerSharedPtr m_linSpacePrecon;
            PreconditionerSharedPtr m_lowEnergyPrecon;

            Array<OneD, NekDouble>  m_invMultiplicity;
	private:

            virtual void v_InitObject();

            virtual void v_DoTransformBasisToLowEnergy(
                Array<OneD, NekDouble>& pInOut);

            virtual void v_DoTransformCoeffsFromLowEnergy(
                Array<OneD, NekDouble>& pInOut);

            virtual void v_DoTransformCoeffsToLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

            virtual void v_DoTransformBasisFromLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);


            virtual DNekScalMatSharedPtr
                v_TransformedSchurCompl(int n, int offset,
                const std::shared_ptr<DNekScalMat > &loc_mat);

            virtual void v_DoPreconditioner(                
                      const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

            virtual void v_BuildPreconditioner();

        };
    }
}

#endif
