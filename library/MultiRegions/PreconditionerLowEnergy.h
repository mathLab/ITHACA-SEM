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
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOWENERGY_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOWENERGY_H
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>


namespace Nektar
{
    namespace MultiRegions
    {
        class PreconditionerLowEnergy;
        typedef boost::shared_ptr<PreconditionerLowEnergy>  PreconditionerLowEnergySharedPtr;

        class PreconditionerLowEnergy: public Preconditioner
	{
        public:
            /// Creates an instance of this class
            static PreconditionerSharedPtr create(
                        const boost::shared_ptr<GlobalLinSys> &plinsys,
                        const boost::shared_ptr<AssemblyMap>
                                                               &pLocToGloMap)
            {
	        PreconditionerSharedPtr p = MemoryManager<PreconditionerLowEnergy>::AllocateSharedPtr(plinsys,pLocToGloMap);
	        p->InitObject();
	        return p;
            }

            /// Name of class
            static std::string className1;
            static std::string className2;
            static std::string className3;

            MULTI_REGIONS_EXPORT PreconditionerLowEnergy(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~PreconditionerLowEnergy() {}

	    /*const DNekMatSharedPtr& GetTransformationMatrix() const;

	    const DNekMatSharedPtr& GetTransposedTransformationMatrix() const;

	    const DNekMatSharedPtr& GetInverseTransformationMatrix() const;

	    const DNekMatSharedPtr& GetInverseTransposedTransformationMatrix() const;*/
	    
	protected:

            const boost::weak_ptr<GlobalLinSys>         m_linsys;

            PreconditionerType                          m_preconType;
	    StdRegions::StdExpansionSharedPtr           vExp;

            DNekMatSharedPtr                            m_preconditioner;
	    DNekScalBlkMatSharedPtr                     GloBlkMat;

            DNekScalMatSharedPtr                        bnd_mat;

	    DNekMatSharedPtr                            m_vertexedgefacetransformmatrix;
            DNekMatSharedPtr                            m_vertexedgefacecoupling;
	    DNekMatSharedPtr                            m_edgefacecoupling;
	    DNekMatSharedPtr                            m_transformationmatrix;
	    DNekMatSharedPtr                            m_inversetransformationmatrix;
	    DNekMatSharedPtr                            m_transposedtransformationmatrix;
	    DNekMatSharedPtr                            m_inversetransposedtransformationmatrix;
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

	    void InverseLinearSpacePreconditioner(void);

	    void StaticCondInverseLinearSpacePreconditioner(void);

	    void SetUpLowEnergyBasis(void);

            void CreateLinearFiniteElmentSpace(void);

            void CreateReferenceGeometryAndMatrix(void);

	    void SetLowEnergyModes_Rv(void);

	    void SetLowEnergyModes_Ref(void);

	    void SetUpInverseTransformationMatrix(void);

	    void LowEnergyPreconditioner(void);

	    void BlockPreconditioner(void);

	    void VertexEdgeFaceMatrix(void);

            virtual void v_InitObject();

            virtual void v_DoPreconditioner(                
                      const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput);

	    virtual const DNekMatSharedPtr& v_GetTransformationMatrix() const;

	    virtual const DNekMatSharedPtr& v_GetTransposedTransformationMatrix() const;

	    virtual const DNekMatSharedPtr& v_GetInverseTransformationMatrix() const;

	    virtual const DNekMatSharedPtr& v_GetInverseTransposedTransformationMatrix() const;

	};
    }
}

#endif
