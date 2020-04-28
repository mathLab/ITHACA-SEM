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
#ifndef NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOWENERGY_H
#define NEKTAR_LIB_MULTIREGIONS_PRECONDITIONERLOWENERGY_H
#include <MultiRegions/GlobalLinSys.h>
#include <MultiRegions/Preconditioner.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>
#include <LocalRegions/TetExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/HexExp.h>


namespace Nektar
{
    namespace MultiRegions
    {
        class PreconditionerLowEnergy;
        typedef std::shared_ptr<PreconditionerLowEnergy>  PreconditionerLowEnergySharedPtr;

        class PreconditionerLowEnergy: public Preconditioner
	{
        public:
            /// Creates an instance of this class
            static PreconditionerSharedPtr create(
                const std::shared_ptr<GlobalLinSys> &plinsys,
                const std::shared_ptr<AssemblyMap> &pLocToGloMap)
            {
	        PreconditionerSharedPtr p = MemoryManager<PreconditionerLowEnergy>::AllocateSharedPtr(plinsys,pLocToGloMap);
	        p->InitObject();
	        return p;
            }

            /// Name of class
            static std::string className;

            MULTI_REGIONS_EXPORT PreconditionerLowEnergy(
                const std::shared_ptr<GlobalLinSys> &plinsys,
                const AssemblyMapSharedPtr &pLocToGloMap);

            MULTI_REGIONS_EXPORT
                virtual ~PreconditionerLowEnergy() {}

	protected:

	    DNekBlkMatSharedPtr m_BlkMat;
            DNekBlkMatSharedPtr m_RBlk;
            DNekBlkMatSharedPtr m_InvRBlk;

            
            Array<OneD, NekDouble>  m_variablePmask;

            bool m_signChange;
            
            // store how many consecutive similar blocks there
            // are in R and Rinv
            std::vector<std::pair<int,int> >  m_sameBlock;  
            
            virtual void v_DoTransformBasisToLowEnergy(
                Array<OneD, NekDouble>& pInOut);

            virtual void v_DoTransformCoeffsFromLowEnergy(
                Array<OneD, NekDouble>& pInOut);

            virtual void v_DoTransformBasisFromLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

            virtual void v_DoTransformCoeffsToLowEnergy(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);
            
	private:

            void SetupBlockTransformationMatrix(void);

            typedef std::map<LibUtilities::ShapeType, DNekScalMatSharedPtr>
                ShapeToDNekMap;
            typedef std::map<LibUtilities::ShapeType,
                LocalRegions::ExpansionSharedPtr > ShapeToExpMap;
            typedef std::map<LibUtilities::ShapeType,
                Array<OneD, unsigned int> > ShapeToIntArrayMap;
            typedef std::map<LibUtilities::ShapeType,
                Array<OneD, Array<OneD, unsigned int> > >
                ShapeToIntArrayArrayMap;
            
            void SetUpReferenceElements(ShapeToDNekMap &maxRmat,
                                        ShapeToExpMap &maxElmt,
                                        ShapeToIntArrayMap      &vertMapMaxR,
                                        ShapeToIntArrayArrayMap &edgeMapMaxR);
            
            void SetUpPyrMaxRMat(int nummodesmax,
                                 LocalRegions::PyrExpSharedPtr &PyrExp,
                                 ShapeToDNekMap          &maxRmat,
                                 ShapeToIntArrayMap      &vertMapMaxR,
                                 ShapeToIntArrayArrayMap &edgeMapMaxR,
                                 ShapeToIntArrayArrayMap &faceMapMaxR);

            void ReSetTetMaxRMat(int nummodesmax,
                                 LocalRegions::TetExpSharedPtr &TetExp,
                                 ShapeToDNekMap          &maxRmat,
                                 ShapeToIntArrayMap      &vertMapMaxR,
                                 ShapeToIntArrayArrayMap &edgeMapMaxR,
                                 ShapeToIntArrayArrayMap &faceMapMaxR);


            void ReSetPrismMaxRMat(int nummodesmax,
                                   LocalRegions::PrismExpSharedPtr &PirsmExp,
                                   ShapeToDNekMap          &maxRmat,
                                   ShapeToIntArrayMap      &vertMapMaxR,
                                   ShapeToIntArrayArrayMap &edgeMapMaxR,
                                   ShapeToIntArrayArrayMap &faceMapMaxR,
                                   bool UseTetOnly);
            
            DNekMatSharedPtr ExtractLocMat(
                          StdRegions::StdExpansionSharedPtr &locExp,
                          DNekScalMatSharedPtr              &maxRmat,
                          LocalRegions::ExpansionSharedPtr  &expMax,
                          Array<OneD, unsigned int>         &vertMapMaxR,
                          Array<OneD, Array<OneD, unsigned int> > &edgeMapMaxR);
            
            void CreateVariablePMask(void);
            
            SpatialDomains::TetGeomSharedPtr   CreateRefTetGeom(void);
            SpatialDomains::PyrGeomSharedPtr   CreateRefPyrGeom(void);
            SpatialDomains::PrismGeomSharedPtr CreateRefPrismGeom(void);
            SpatialDomains::HexGeomSharedPtr   CreateRefHexGeom(void);

            virtual void v_InitObject();

            virtual void v_DoPreconditioner(                
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput);

            virtual void v_BuildPreconditioner();
            
            virtual DNekScalMatSharedPtr
                v_TransformedSchurCompl(int n, int offset, 
                             const std::shared_ptr<DNekScalMat > &loc_mat);
        };
    }
}

#endif
