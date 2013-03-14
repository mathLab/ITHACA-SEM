///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.h
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
// Description: Header file for SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SEGEXP_H
#define SEGEXP_H

#include <StdRegions/StdSegExp.h>
#include <SpatialDomains/Geometry1D.h>
#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

//#include <fstream>

namespace Nektar
{
    namespace LocalRegions
    {

        class SegExp: virtual public StdRegions::StdSegExp, virtual public Expansion1D
        {

        public:
            LOCAL_REGIONS_EXPORT SegExp(
                    const LibUtilities::BasisKey &Ba,
                    const SpatialDomains::Geometry1DSharedPtr &geom);

            LOCAL_REGIONS_EXPORT SegExp(const SegExp &S);

            LOCAL_REGIONS_EXPORT ~SegExp();

        protected:
            //----------------------------
            // Integration Methods
            //----------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_Integral(
                    const Array<OneD, const NekDouble>& inarray);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,NekDouble> &out_d0,
                          Array<OneD,NekDouble> &out_d1 = NullNekDouble1DArray,
                          Array<OneD,NekDouble> &out_d2 = NullNekDouble1DArray);

            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv_s(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &out_ds);

            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv_n(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble>& out_dn);

            //-----------------------------
            // Transforms
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans_BndConstrained(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            //-----------------------------
            // Inner product functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& base,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          int coll_check);

            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                    const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> & outarray);

            LOCAL_REGIONS_EXPORT virtual void v_NormVectorIProductWRTBase(
                    const Array<OneD, const NekDouble> &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                          Array<OneD, NekDouble> &outarray);

            //-----------------------------
            // Evaluation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coord);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& coord,
                    const Array<OneD, const NekDouble> & physvals);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(
                    Array<OneD,NekDouble> &coords_1,
                    Array<OneD,NekDouble> &coords_2 = NullNekDouble1DArray,
                    Array<OneD,NekDouble> &coords_3 = NullNekDouble1DArray);

            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(
                    const Array<OneD, const NekDouble>& Lcoords,
                          Array<OneD,NekDouble> &coords);

            //-----------------------------
            // Helper functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_WriteToFile(
                          std::ofstream &outfile,
                          OutputFormat format,
                    const bool dumpVar = true,
                          std::string var = "v");

            LOCAL_REGIONS_EXPORT virtual int v_GetCoordim();

            LOCAL_REGIONS_EXPORT virtual const
                    SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const;

            LOCAL_REGIONS_EXPORT virtual const
                    SpatialDomains::GeometrySharedPtr v_GetGeom() const;

            LOCAL_REGIONS_EXPORT virtual const
                    SpatialDomains::Geometry1DSharedPtr& v_GetGeom1D() const;

            LOCAL_REGIONS_EXPORT virtual void 
                    v_SetCoeffsToOrientation(StdRegions::Orientation dir);

            LOCAL_REGIONS_EXPORT virtual void v_SetCoeffsToOrientation(
                    StdRegions::Orientation dir,
                    Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual int v_GetNumPoints(const int dir) const;

            LOCAL_REGIONS_EXPORT virtual int v_GetNcoeffs(void) const;

            LOCAL_REGIONS_EXPORT virtual const
                    LibUtilities::BasisSharedPtr&  v_GetBasis(int dir) const;

            LOCAL_REGIONS_EXPORT virtual int v_NumBndryCoeffs() const;

            LOCAL_REGIONS_EXPORT virtual int v_NumDGBndryCoeffs() const;

            LOCAL_REGIONS_EXPORT virtual void v_ComputeVertexNormal(const int vertex);

            LOCAL_REGIONS_EXPORT virtual StdRegions::Orientation v_GetPorient(int point);

            LOCAL_REGIONS_EXPORT virtual SpatialDomains::GeomType  v_MetricInfoType();

            LOCAL_REGIONS_EXPORT virtual void v_ExtractDataToCoeffs(
                    const NekDouble *data,
                    const std::vector<unsigned int > &nummodes,
                    const int mode_offset,
                    NekDouble *coeffs);


            LOCAL_REGIONS_EXPORT virtual void v_SetUpPhysTangents(
                    const StdRegions::StdExpansionSharedPtr &exp2D,
                    const int edge);

            LOCAL_REGIONS_EXPORT virtual const
                    Array<OneD, const NekDouble>&  v_GetPhysNormals(void);

            LOCAL_REGIONS_EXPORT virtual NekDouble
                    v_Linf(const Array<OneD, const NekDouble>& sol);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_Linf();

            LOCAL_REGIONS_EXPORT virtual NekDouble
                    v_L2(const Array<OneD, const NekDouble>& sol);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_L2();

            //-----------------------------
            // Operator creation functions
            //-----------------------------
            LOCAL_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::StdMatrixKey     &mkey);

            LOCAL_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                    const StdRegions::StdMatrixKey     &mkey);

            //-----------------------------
            // Matrix creation functions
            //-----------------------------

            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr
                    v_GenMatrix(const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalMatSharedPtr CreateMatrix(const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekMatSharedPtr
                    v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);

            LOCAL_REGIONS_EXPORT DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(
                    const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekScalMatSharedPtr
                    v_GetLocMatrix(const MatrixKey &mkey);

            LOCAL_REGIONS_EXPORT virtual DNekScalBlkMatSharedPtr
                    v_GetLocStaticCondMatrix(const MatrixKey &mkey);


        private:
            SpatialDomains::Geometry1DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess>
                    m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess>
                    m_staticCondMatrixManager;

            SegExp();

            LOCAL_REGIONS_EXPORT void ReverseCoeffsAndSign(
                    const Array<OneD,NekDouble> &inarray,
                          Array<OneD,NekDouble> &outarray);


            /// \todo Same method exists in ExpList and everyone references
            ///       ExpList::MultiplyByElmtInvMass. Remove this one?
            LOCAL_REGIONS_EXPORT void MultiplyByElmtInvMass(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,NekDouble> &outarray);

        };

        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<SegExp>      SegExpSharedPtr;
        typedef std::vector< SegExpSharedPtr > SegExpVector;
        typedef std::vector< SegExpSharedPtr >::iterator SegExpVectorIter;


    } //end of namespace
} //end of namespace

#endif // SEGEXP_H
