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

#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONER_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONER_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Communication/GsLib.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

#include <memory>

namespace Nektar
{
    namespace MultiRegions
    {
        class AssemblyMap;
        typedef std::shared_ptr<AssemblyMap> AssemblyMapSharedPtr;

        class Preconditioner;
        typedef std::shared_ptr<Preconditioner>  PreconditionerSharedPtr;

        static PreconditionerSharedPtr NullPreconditionerSharedPtr;

        typedef LibUtilities::NekFactory< std::string, Preconditioner, 
            const std::shared_ptr<GlobalLinSys>&,
            const std::shared_ptr<AssemblyMap>& > PreconFactory;
        PreconFactory& GetPreconFactory();

        class Preconditioner
        {
        public:
            MULTI_REGIONS_EXPORT Preconditioner(
                         const std::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~Preconditioner() {}

	    inline void DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

            inline void DoPreconditionerWithNonVertOutput(
                         const Array<OneD, NekDouble>& pInput,
                         Array<OneD, NekDouble>& pOutput,
                         const Array<OneD, NekDouble>& pNonVertOutput,
                         Array<OneD, NekDouble>& pVertForce = NullNekDouble1DArray);


	    inline void DoTransformBasisToLowEnergy(
                                     Array<OneD, NekDouble>& pInOut);

	    inline void DoTransformCoeffsFromLowEnergy(
                                     Array<OneD, NekDouble>& pInOut);

	    inline void DoTransformCoeffsToLowEnergy(
                                     const Array<OneD, NekDouble>& pInput,
                                     Array<OneD, NekDouble>& pOutput);

	    inline void DoTransformBasisFromLowEnergy(
                                     const Array<OneD, NekDouble>& pInput,
                                     Array<OneD, NekDouble>& pOutput);
            
	    inline void BuildPreconditioner();

   	    inline void InitObject();

            Array<OneD, NekDouble> AssembleStaticCondGlobalDiagonals();

             inline const DNekScalBlkMatSharedPtr&
                GetBlockTransformedSchurCompl() const;
            
            inline const DNekScalBlkMatSharedPtr&
                GetBlockCMatrix() const;
            
            inline const DNekScalBlkMatSharedPtr&
                GetBlockInvDMatrix() const;
            
            inline const DNekScalBlkMatSharedPtr&
                GetBlockSchurCompl() const;
        
            inline const DNekScalBlkMatSharedPtr&
                GetBlockTransformationMatrix() const;
            
            inline const DNekScalBlkMatSharedPtr&
                GetBlockTransposedTransformationMatrix() const;

            inline DNekScalMatSharedPtr TransformedSchurCompl(
                   int offset, int bndoffset, 
                   const std::shared_ptr<DNekScalMat > &loc_mat);

	protected:
            const std::weak_ptr<GlobalLinSys>   m_linsys;
            PreconditionerType                  m_preconType;
            DNekMatSharedPtr                    m_preconditioner;
            std::weak_ptr<AssemblyMap>          m_locToGloMap;
            LibUtilities::CommSharedPtr         m_comm;

            virtual DNekScalMatSharedPtr v_TransformedSchurCompl(
                        int offset, int bndoffset,
                        const std::shared_ptr<DNekScalMat > &loc_mat);


	private:

            void NullPreconditioner(void);

	    virtual void v_InitObject();

	    virtual void v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

            virtual void v_DoPreconditionerWithNonVertOutput(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput,
                const Array<OneD, NekDouble>& pNonVertOutput,
                Array<OneD, NekDouble>& pVertForce);

            
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

	    virtual void v_BuildPreconditioner();

            static std::string lookupIds[];
            static std::string def;
	};
        typedef std::shared_ptr<Preconditioner>  PreconditionerSharedPtr;

        /**
         *
         */
        inline void Preconditioner::InitObject()
        {
            v_InitObject();
        }

        /**
         *
         */ 
        inline DNekScalMatSharedPtr Preconditioner::TransformedSchurCompl(
                            int offset, int bndoffset,
                            const std::shared_ptr<DNekScalMat > &loc_mat)
        {
            return v_TransformedSchurCompl(offset,bndoffset,loc_mat);
        }

        /**
         *
         */
        inline void Preconditioner::DoPreconditioner(
            const Array<OneD, NekDouble> &pInput,
                  Array<OneD, NekDouble> &pOutput)
        {
	    v_DoPreconditioner(pInput,pOutput);
        }
        

        /**
         *
         */
        inline void Preconditioner::DoPreconditionerWithNonVertOutput(
            const Array<OneD, NekDouble>& pInput,
                  Array<OneD, NekDouble>& pOutput,
            const Array<OneD, NekDouble>& pNonVertOutput,
                  Array<OneD, NekDouble>& pVertForce)
        {
            v_DoPreconditionerWithNonVertOutput(pInput,pOutput,pNonVertOutput,
                                                pVertForce);
        }

        /**
         *
         */
        inline void Preconditioner::DoTransformBasisToLowEnergy(
            Array<OneD, NekDouble>& pInOut)
        {
	    v_DoTransformBasisToLowEnergy(pInOut);
        }

        /**
         *
         */
        inline void Preconditioner::DoTransformCoeffsFromLowEnergy(
            Array<OneD, NekDouble>& pInput)
        {
	    v_DoTransformCoeffsFromLowEnergy(pInput);
        }

        /**
         *
         */
        inline void Preconditioner::DoTransformCoeffsToLowEnergy(
            const Array<OneD, NekDouble>& pInput,
            Array<OneD, NekDouble>& pOutput)
        {
            v_DoTransformCoeffsToLowEnergy(pInput,pOutput);
        }
           
        /**
         *
         */
        inline void Preconditioner::DoTransformBasisFromLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                      Array<OneD, NekDouble>& pOutput)
        {
            v_DoTransformBasisFromLowEnergy(pInput,pOutput);
        }

        /**
         *
         */
        inline void Preconditioner::BuildPreconditioner()
        {
	    v_BuildPreconditioner();
        }
    }
}

#endif
