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
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/HexExp.h>


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
            static std::string className;

            MULTI_REGIONS_EXPORT PreconditionerLowEnergy(
                const boost::shared_ptr<GlobalLinSys> &plinsys,
                const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
                virtual ~PreconditionerLowEnergy() {}

	protected:

            const boost::weak_ptr<GlobalLinSys> m_linsys;
            boost::shared_ptr<AssemblyMap> m_locToGloMap;

	    DNekBlkMatSharedPtr m_BlkMat;
            DNekBlkMatSharedPtr m_RBlk;

            DNekBlkMatSharedPtr m_InvRBlk;

            int m_nummodesmax;

            std::map<LibUtilities::ShapeType, DNekScalMatSharedPtr> m_maxRmat;
            std::map<LibUtilities::ShapeType, Array<OneD, unsigned int> > m_vertMapMaxR; 
            std::map<LibUtilities::ShapeType, Array<OneD, Array<OneD, unsigned int> > > m_edgeMapMaxR; 
            std::map<LibUtilities::ShapeType, LocalRegions::ExpansionSharedPtr > m_maxElmt; 
            
            Array<OneD, NekDouble>  m_locToGloSignMult;
            Array<OneD, NekDouble>  m_multiplicity;
            Array<OneD, int>        m_map;

            bool m_signChange;
            
            // store how many consecutive similar blocks there are in R and Rinv
            std::vector<std::pair<int,int> >  m_sameBlock;  
            
	private:

            void SetUpReferenceElements(void);
            
            void CreateMultiplicityMap(void);

            void SetupBlockTransformationMatrix(void);

            DNekMatSharedPtr ExtractLocMat(StdRegions::StdExpansionSharedPtr  &locExp);


        void ModifyPrismTransformationMatrix(
                LocalRegions::TetExpSharedPtr TetExp,
                LocalRegions::PrismExpSharedPtr PrismExp,
                DNekMatSharedPtr Rmodprism,
                DNekMatSharedPtr RTmodprism);

            void LocalTransformToLowEnergy(
                DNekScalMatSharedPtr RTmat,
                LocalRegions::HexExpSharedPtr maxTetExp);

            SpatialDomains::TetGeomSharedPtr CreateRefTetGeom(void);
            SpatialDomains::PrismGeomSharedPtr CreateRefPrismGeom(void);
            SpatialDomains::HexGeomSharedPtr CreateRefHexGeom(void);

            virtual void v_InitObject();

            virtual void v_DoPreconditioner(                
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

            virtual void v_DoTransformToLowEnergy(
                Array<OneD, NekDouble>& pInOut,
                int offset);

            virtual void v_DoTransformToLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);
            
            virtual void v_DoTransformFromLowEnergy(
                Array<OneD, NekDouble>& pInOut);

            virtual void v_DoMultiplybyInverseTransformationMatrix(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

            virtual void v_DoMultiplybyInverseTransposedTransformationMatrix(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);
            
            virtual void v_BuildPreconditioner();
            
            virtual DNekScalMatSharedPtr
                v_TransformedSchurCompl(int n, int offset, 
                                        const boost::shared_ptr<DNekScalMat > &loc_mat);
        };
    }
}

#endif
