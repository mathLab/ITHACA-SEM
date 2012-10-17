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
// Description: Preconditioner header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONER_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONER_H
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>


namespace Nektar
{
    namespace MultiRegions
    {
        class Preconditioner;
        typedef boost::shared_ptr<Preconditioner>  PreconditionerSharedPtr;

        typedef LibUtilities::NekFactory< std::string, Preconditioner, 
            const boost::shared_ptr<GlobalLinSys>&,
            const boost::shared_ptr<AssemblyMap>& > PreconFactory;
        PreconFactory& GetPreconFactory();

        class Preconditioner
	{
        public:
            MULTI_REGIONS_EXPORT Preconditioner(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~Preconditioner() {}

	    inline void DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

   	    inline void InitObject();

            Array<OneD, NekDouble> AssembleStaticCondGlobalDiagonals();

            const inline DNekMatSharedPtr &GetTransformationMatrix(void) const;

            const inline DNekMatSharedPtr &GetTransposedTransformationMatrix(void) const;

	protected:

            const boost::weak_ptr<GlobalLinSys>         m_linsys;

            PreconditionerType                          m_preconType;

            DNekMatSharedPtr                            m_preconditioner;

            boost::shared_ptr<AssemblyMap>              m_locToGloMap;

	private:

            void NullPreconditioner(void);

	    virtual void v_InitObject();

	    virtual void v_DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

            virtual const DNekMatSharedPtr& v_GetTransformationMatrix(void) const;

            virtual const DNekMatSharedPtr& v_GetTransposedTransformationMatrix(void) const;

            static std::string lookupIds[];
            static std::string def;
	};
        typedef boost::shared_ptr<Preconditioner>  PreconditionerSharedPtr;

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
        inline const DNekMatSharedPtr& Preconditioner::GetTransformationMatrix(void) const
        {
	  return v_GetTransformationMatrix();
        }

        /**
         *
         */
        inline const DNekMatSharedPtr& Preconditioner::GetTransposedTransformationMatrix(void) const
        {
	  return v_GetTransposedTransformationMatrix();
        }

        /**
         *
         */
        inline void Preconditioner::DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput)
        {
	    v_DoPreconditioner(pInput,pOutput);
        }

    }
}

#endif
