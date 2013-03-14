///////////////////////////////////////////////////////////////////////////////
//
// File TetExp.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TETEXP_H
#define TETEXP_H

#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/TetGeom.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions
    {

        class TetExp: virtual public StdRegions::StdTetExp, virtual public Expansion3D
        {

        public:
            LOCAL_REGIONS_EXPORT TetExp( const LibUtilities::BasisKey &Ba,
                            const LibUtilities::BasisKey &Bb,
                            const LibUtilities::BasisKey &Bc,
                            const SpatialDomains::TetGeomSharedPtr &geom);

            LOCAL_REGIONS_EXPORT TetExp(const TetExp &T);

            LOCAL_REGIONS_EXPORT ~TetExp();

        protected:
            //-----------------------------
            // Integration Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual
                NekDouble v_Integral(const Array<OneD, const NekDouble> &inarray);


            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD, NekDouble> &out_d0,
                      Array<OneD, NekDouble> &out_d1,
                      Array<OneD, NekDouble> &out_d2);

            //-----------------------------
            // Transforms
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans(
                const Array<OneD, const NekDouble> & inarray,
                      Array<OneD,NekDouble> &outarray);

            //-----------------------------
            // Inner product functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);


            //-----------------------------
            // Evaluation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coords);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble> &coords,
                const Array<OneD, const NekDouble> & physvals);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(
                Array<OneD,NekDouble> &coords_0,
                Array<OneD,NekDouble> &coords_1,
                Array<OneD,NekDouble> &coords_2);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(
                const Array<OneD, const NekDouble> &Lcoords,
                      Array<OneD,NekDouble> &coords);

            //-----------------------------
            // Helper functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_WriteToFile( 
                std::ofstream &outfile,
                OutputFormat format,
                const bool dumpVar = true,
                std::string var = "v");

            LOCAL_REGIONS_EXPORT virtual 
                LibUtilities::ShapeType v_DetShapeType() const;

            LOCAL_REGIONS_EXPORT virtual const 
                SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const;

            LOCAL_REGIONS_EXPORT virtual const 
                SpatialDomains::GeometrySharedPtr v_GetGeom() const;

            LOCAL_REGIONS_EXPORT virtual const 
                SpatialDomains::Geometry3DSharedPtr& v_GetGeom3D() const;

            LOCAL_REGIONS_EXPORT virtual int v_GetCoordim();

            LOCAL_REGIONS_EXPORT virtual 
                StdRegions::Orientation v_GetFaceOrient(int face);

            LOCAL_REGIONS_EXPORT virtual void v_GetFacePhysVals(
                const int                                face,
                const StdRegions::StdExpansionSharedPtr &FaceExp,
                const Array<OneD, const NekDouble>      &inarray,
                      Array<OneD,       NekDouble>      &outarray,
                StdRegions::Orientation                  orient);

            LOCAL_REGIONS_EXPORT void v_ComputeFaceNormal(const int face);      

            LOCAL_REGIONS_EXPORT virtual 
                NekDouble v_Linf(const Array<OneD, const NekDouble> &sol);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_Linf();

            LOCAL_REGIONS_EXPORT virtual 
                NekDouble v_L2(const Array<OneD, const NekDouble> &sol);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_L2();


            //-----------------------------
            // Operator creation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp_MatFree(
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

            LOCAL_REGIONS_EXPORT virtual void v_LaplacianMatrixOp_MatFree(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,NekDouble> &outarray,
                const StdRegions::StdMatrixKey &mkey);


            //-----------------------------
            // Matrix creation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr  CreateMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekScalMatSharedPtr v_GetLocMatrix(
                const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual 
                DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(
                const MatrixKey &mkey);


        private:
            SpatialDomains::Geometry3DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            TetExp();

            LOCAL_REGIONS_EXPORT void GeneralMatrixOp_MatOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const StdRegions::StdMatrixKey      &mkey);
            LOCAL_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            LOCAL_REGIONS_EXPORT void LaplacianMatrixOp_MatFree_Kernel(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                      Array<OneD,       NekDouble> &wsp);

	    LOCAL_REGIONS_EXPORT SpatialDomains::TetGeomSharedPtr CreateEquilateralTetGeom();
        };

        // type defines for use of TetExp in a boost vector
        typedef boost::shared_ptr<TetExp> TetExpSharedPtr;
        typedef std::vector< TetExpSharedPtr > TetExpVector;
        typedef std::vector< TetExpSharedPtr >::iterator TetExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // TETEXP_H
