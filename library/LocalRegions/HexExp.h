///////////////////////////////////////////////////////////////////////////////
//
// File HexExp.h
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
// Description: Header routines for Hex expansion 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef HEXEXP_H
#define HEXEXP_H

#include <LibUtilities/BasicUtils/NekManager.hpp>  // for NekManager

#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion3D.h>
#include <StdRegions/StdHexExp.h>
#include <SpatialDomains/HexGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class HexGeom;
        typedef std::shared_ptr<HexGeom> HexGeomSharedPtr;
    }

    namespace LocalRegions 
    {
        class HexExp;

        typedef std::shared_ptr<HexExp> HexExpSharedPtr;
        typedef std::vector< HexExpSharedPtr > HexExpVector;
 
        class HexExp: virtual public StdRegions::StdHexExp, virtual public Expansion3D
        {
        public:
    
            LOCAL_REGIONS_EXPORT HexExp( const LibUtilities::BasisKey &Ba, 
                            const LibUtilities::BasisKey &Bb,
                            const LibUtilities::BasisKey &Bc, 
                            const SpatialDomains::HexGeomSharedPtr &geom);

            LOCAL_REGIONS_EXPORT HexExp( const HexExp &T );

            LOCAL_REGIONS_EXPORT ~HexExp();

        protected:
            //------------------------------
            //    Integration Method
            //------------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_Integral( 
                const Array<OneD, const NekDouble> &inarray );


            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD, NekDouble> &out_d0,
                      Array<OneD, NekDouble> &out_d1,
                      Array<OneD, NekDouble> &out_d2);  
            
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(const int dir,
                const Array<OneD, const NekDouble>& inarray,
                Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void v_PhysDirectionalDeriv(
                const Array<OneD, const NekDouble> &inarray,
                const Array<OneD, const NekDouble> &direction,
                      Array<OneD,       NekDouble> &out);


            //---------------------------------------
            // Transforms
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans(
                const Array<OneD, const NekDouble> & inarray, 
                Array<OneD,NekDouble> &outarray);


            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble> &inarray,
                Array<OneD,       NekDouble> &outarray,
                bool multiplybyweights = true);

            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void IProductWRTDerivBase_SumFac(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void IProductWRTDerivBase_MatOp(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void v_IProductWRTDirectionalDerivBase(
                const Array<OneD, const NekDouble>& direction,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray)
            {
                IProductWRTDirectionalDerivBase_SumFac(direction, inarray,
                                                       outarray);
            }

            LOCAL_REGIONS_EXPORT void IProductWRTDirectionalDerivBase_SumFac(
                const Array<OneD, const NekDouble>& direction,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);


            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_StdPhysEvaluate(
                const Array<OneD, const NekDouble> &Lcoord,
                const Array<OneD, const NekDouble> &physvals);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coords, 
                const Array<OneD, const NekDouble> & physvals);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(
                const Array<OneD, const NekDouble> &Lcoords, 
                      Array<OneD,NekDouble> &coords);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(
                      Array<OneD,       NekDouble> &coords_1,
                      Array<OneD,       NekDouble> &coords_2,
                      Array<OneD,       NekDouble> &coords_3);

      
            //---------------------------------------
            // Helper functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual 
                LibUtilities::ShapeType v_DetShapeType() const;
    
            LOCAL_REGIONS_EXPORT virtual
                StdRegions::StdExpansionSharedPtr v_GetStdExp(void) const;

            LOCAL_REGIONS_EXPORT virtual
                StdRegions::StdExpansionSharedPtr v_GetLinStdExp(void) const;
            
            LOCAL_REGIONS_EXPORT virtual void v_ExtractDataToCoeffs(
                const NekDouble *data,
                const std::vector<unsigned int > &nummodes,  
                const int mode_offset,
                NekDouble * coeffs,
                std::vector<LibUtilities::BasisType> &fromType);

            LOCAL_REGIONS_EXPORT virtual 
                bool v_GetFaceDGForwards(const int i) const;

            LOCAL_REGIONS_EXPORT virtual void v_GetFacePhysMap(
                const int                face,
                Array<OneD, int>        &outarray);

            LOCAL_REGIONS_EXPORT void v_ComputeFaceNormal(const int face);      

            //---------------------------------------
            // Operator creation functions
            //---------------------------------------
	    LOCAL_REGIONS_EXPORT virtual void v_MassMatrixOp(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

	    LOCAL_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

	    LOCAL_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                const int k1, 
                const int k2, 
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

	    LOCAL_REGIONS_EXPORT virtual void v_WeakDerivMatrixOp(const int i,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);
            
	    LOCAL_REGIONS_EXPORT virtual void v_WeakDirectionalDerivMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);
            
	    LOCAL_REGIONS_EXPORT virtual void v_MassLevelCurvatureMatrixOp(
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

	    LOCAL_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT void v_GeneralMatrixOp_MatOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey); 
            
            LOCAL_REGIONS_EXPORT virtual void v_ReduceOrderCoeffs(
                int                                 numMin,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            
            LOCAL_REGIONS_EXPORT virtual void v_SVVLaplacianFilter(
                    Array<OneD, NekDouble> &array,
                    const StdRegions::StdMatrixKey &mkey);

            //---------------------------------------
            // Matrix creation functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr  CreateMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekScalMatSharedPtr v_GetLocMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual 
                DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT void v_DropLocStaticCondMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual void v_ComputeLaplacianMetric();


        private:
            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            HexExp();

            virtual void v_LaplacianMatrixOp_MatFree_Kernel(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                      Array<OneD,       NekDouble> &wsp);
        };

    } //end of namespace
} //end of namespace

#endif //HEX_EXP_H
