///////////////////////////////////////////////////////////////////////////////
//
// File StdTetExp.h
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDTETEXP_H
#define NEKTAR_LIB_STDREGIONS_STDTETEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {
        class StdMatrixKey;

        namespace StdTetData
        {
            /**
             * Adds up the number of cells in a truncated Nc by Nc by Nc
             * pyramid, where the longest Na rows and longest Nb columns are
             * kept. Example: (Na, Nb, Nc) = (3, 4, 5); The number of
             * coefficients is the sum of the elements of the following
             * matrix:
             *
             * |5  4  3  2  0|
             * |4  3  2  0   |
             * |3  2  0      |
             * |0  0         |
             * |0            |
             * 
             * Sum = 28 = number of tet coefficients.
             */
            inline int getNumberOfCoefficients(int Na, int Nb, int Nc)
            {
                int nCoef = 0;
                for (int a = 0; a < Na; ++a)
                {
                    for (int b = 0; b < Nb - a; ++b)
                    {
                        for (int c = 0; c < Nc - a - b; ++c)
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

        class StdTetExp : virtual public StdExpansion3D
        {

        public:

            STD_REGIONS_EXPORT StdTetExp();
            STD_REGIONS_EXPORT StdTetExp(
                const LibUtilities::BasisKey &Ba,
                const LibUtilities::BasisKey &Bb, 
                const LibUtilities::BasisKey &Bc);
            STD_REGIONS_EXPORT StdTetExp(
                const LibUtilities::BasisKey &Ba,
                const LibUtilities::BasisKey &Bb,
                const LibUtilities::BasisKey &Bc, 
                NekDouble *coeffs,
                NekDouble *phys);
            STD_REGIONS_EXPORT StdTetExp(const StdTetExp &T);
            STD_REGIONS_EXPORT ~StdTetExp();

            ExpansionType DetExpansionType() const
            {
                return eTetrahedron;
            }





            /** \brief Single Point Evaluation */
            STD_REGIONS_EXPORT NekDouble PhysEvaluate3D(const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT NekDouble PhysEvaluate3D(const Array<OneD, const NekDouble>& coords,  const Array<OneD, const NekDouble> & physvals);


        protected:
            //-------------------------------
            // Integration Methods
            //-------------------------------
            STD_REGIONS_EXPORT void TripleTensorProduct(
                const Array<OneD, const NekDouble>& fx,
                const Array<OneD, const NekDouble>& gy,
                const Array<OneD, const NekDouble>& hz,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT NekDouble TripleInnerProduct(
                const Array<OneD, const NekDouble>& fxyz,
                const Array<OneD, const NekDouble>& wx,
                const Array<OneD, const NekDouble>& wy,
                const Array<OneD, const NekDouble>& wz);
            STD_REGIONS_EXPORT NekDouble Integral3D(
                const Array<OneD, const NekDouble>& inarray,
                const Array<OneD, const NekDouble>& w0,
                const Array<OneD, const NekDouble>& w1,
                const Array<OneD, const NekDouble>& w2);
            STD_REGIONS_EXPORT NekDouble v_Integral(
                const Array<OneD, const NekDouble>& inarray);

            //----------------------------
            // Differentiation Methods
            //----------------------------
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& out_dx,
                      Array<OneD,       NekDouble>& out_dy,
                      Array<OneD,       NekDouble>& out_dz);
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const int dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& out_d0,
                      Array<OneD,       NekDouble>& out_d1,
                      Array<OneD,       NekDouble>& out_d2);

            //---------------------------------------
            // Transforms
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT void BwdTrans_SumFacKernel(
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray,
                      Array<OneD,       NekDouble>& wsp,
                      bool                          doCheckCollDir0,
                      bool                          doCheckCollDir1,
                      bool                          doCheckCollDir2);
            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray );

            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_MatOp(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT void IProductWRTBase_SumFacKernel (
                const Array<OneD, const NekDouble>& base0,
                const Array<OneD, const NekDouble>& base1,
                const Array<OneD, const NekDouble>& base2,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray,
                      Array<OneD,       NekDouble>& wsp,
                      bool                          doCheckCollDir0,
                      bool                          doCheckCollDir1,
                      bool                          doCheckCollDir2);
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_MatOp(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_SumFac(
                const int                           dir,
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);


            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coords);
            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& coords,
                const Array<OneD, const NekDouble>& physvals);
            STD_REGIONS_EXPORT virtual void v_GetCoords(
                Array<OneD, NekDouble> &coords_x,
                Array<OneD, NekDouble> &coords_y,
                Array<OneD, NekDouble> &coords_z);
            STD_REGIONS_EXPORT virtual void v_FillMode(
                const int                     mode, 
                      Array<OneD, NekDouble>& outarray);
            
            
            //---------------------------
            // Helper functions
            //---------------------------
            STD_REGIONS_EXPORT virtual int  v_GetNverts() const;
            STD_REGIONS_EXPORT virtual int  v_GetNedges() const;
            STD_REGIONS_EXPORT virtual int  v_GetNfaces() const;
            STD_REGIONS_EXPORT virtual ExpansionType v_DetExpansionType() const;
            STD_REGIONS_EXPORT virtual int  v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int  v_NumDGBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int  v_GetEdgeNcoeffs(const int i) const;
 	    STD_REGIONS_EXPORT virtual int  v_GetTotalEdgeIntNcoeffs() const;
            STD_REGIONS_EXPORT virtual int  v_GetFaceNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int  v_GetFaceIntNcoeffs(
                    const int i) const;
 	    STD_REGIONS_EXPORT virtual int  v_GetTotalFaceIntNcoeffs() const;
            STD_REGIONS_EXPORT virtual int  v_GetFaceNumPoints(const int i) const;
            STD_REGIONS_EXPORT virtual LibUtilities::PointsKey v_GetFacePointsKey(
                    const int i, const int j) const;
            STD_REGIONS_EXPORT virtual int  v_CalcNumberOfCoefficients(
                    const std::vector<unsigned int> &nummodes,
                          int                       &modes_offset);
            STD_REGIONS_EXPORT virtual const LibUtilities::BasisKey 
                    v_DetFaceBasisKey(const int i, const int k) const;
            STD_REGIONS_EXPORT virtual LibUtilities::BasisType v_GetEdgeBasisType(
                    const int i) const;
            STD_REGIONS_EXPORT virtual void v_WriteToFile(
                    std::ofstream &outfile,
                    OutputFormat format,
                    const bool dumpVar = true,
                    std::string var = "v");
            STD_REGIONS_EXPORT virtual void v_WriteCoeffsToFile(
                    std::ofstream &outfile);
            STD_REGIONS_EXPORT virtual bool v_IsBoundaryInteriorExpansion();


            //--------------------------
            // Mappings
            //--------------------------
            STD_REGIONS_EXPORT virtual void v_GetFaceToElementMap(
                    const int                  fid,
                    const Orientation      faceOrient,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD,          int> &signarray,
                    int                        nummodesA = -1,
                    int                        nummodesB = -1);
            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(int localVertexId);
            STD_REGIONS_EXPORT virtual void v_GetEdgeInteriorMap(const int eid,
                    const Orientation edgeOrient,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD, int> &signarray);
            STD_REGIONS_EXPORT virtual void v_GetFaceInteriorMap(
                    const int fid,
                    const Orientation faceOrient,
                    Array<OneD, unsigned int> &maparray,
                    Array<OneD, int>& signarray);
            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                    Array<OneD, unsigned int>& outarray);
            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                    Array<OneD, unsigned int>& outarray);


            //---------------------------------------
            // Wrapper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD,       NekDouble>& outarray);

        private:
            //---------------------------------------
            // Private helper functions
            //---------------------------------------
            STD_REGIONS_EXPORT int  GetMode(const int i, const int j, const int k);
        };

        typedef boost::shared_ptr<StdTetExp> StdTetExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDTETEXP_H
