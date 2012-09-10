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
        class Preconditioner
	{
        public:
            MULTI_REGIONS_EXPORT Preconditioner(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~Preconditioner() {}

	    void DoPreconditioner(
                const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

	protected:

            const boost::weak_ptr<GlobalLinSys>         m_linsys;

            PreconditionerType                          m_preconType;
	    StdRegions::StdExpansionSharedPtr           vExp;

            DNekMatSharedPtr                            m_preconditioner;
	    DNekScalBlkMatSharedPtr                     GloBlkMat;

            DNekScalMatSharedPtr                        bndry_mat;

	    DNekMatSharedPtr                            m_vertexedgefacetransformmatrix;
            DNekMatSharedPtr                            m_vertexedgefacecoupling;
	    DNekMatSharedPtr                            m_edgefacecoupling;
	    DNekMatSharedPtr                            m_transformationmatrix;
	    DNekMatSharedPtr                            m_transposedtransformationmatrix;
	    DNekMatSharedPtr                            m_efedgefacecoupling;
	    DNekMatSharedPtr                            m_effacefacecoupling;
	    DNekMatSharedPtr                            m_edgefacetransformmatrix;

            boost::shared_ptr<AssemblyMap>              m_locToGloMap;

            Array<OneD, int>                            vertModeLocation;
            Array<OneD, Array<OneD, unsigned int> >     edgeModeLocation;
            Array<OneD, Array<OneD, unsigned int> >     faceModeLocation;

            Array<OneD, Array<OneD, unsigned int> >     MatEdgeLocation;
            Array<OneD, Array<OneD, unsigned int> >     MatFaceLocation;

	private:

            void NullPreconditioner(void);

            void DiagonalPreconditionerSum(void);

	    void StaticCondDiagonalPreconditionerSum(void);

	    void InverseLinearSpacePreconditioner(void);

	    void StaticCondInverseLinearSpacePreconditioner(void);

	    void SetUpLowEnergyBasis(void);

            void CreateLinearFiniteElmentSpace(void);

            void CreateReferenceGeometryAndMatrix(void);

	    void SetLowEnergyModes_Rv(void);

	    void SetLowEnergyModes_Ref(void);

	    void LowEnergyPreconditioner(void);

	    void BlockPreconditioner(void);

	    void VertexEdgeFaceMatrix(void);

            Array<OneD, NekDouble> AssembleStaticCondGlobalDiagonals();

            static std::string lookupIds[];
            static std::string def;
	};
        typedef boost::shared_ptr<Preconditioner>  PreconditionerSharedPtr;
    }
}

#endif
