///////////////////////////////////////////////////////////////////////////////
//
// File StdSegExp.h
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
// Description: Header file for Standard Segment Expansions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H
#define NEKTAR_LIBS_STDREGIONS_STDSEGEXP_H

#include <StdRegions/StdExpansion1D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

        class StdSegExp;
        typedef std::shared_ptr<StdSegExp> StdSegExpSharedPtr;

        /// Class representing a segment element in reference space

        /// All interface of this class sits in StdExpansion class

        class StdSegExp: virtual public StdExpansion1D
        {
        public:

            STD_REGIONS_EXPORT StdSegExp();

            STD_REGIONS_EXPORT StdSegExp(const LibUtilities::BasisKey &Ba);

            STD_REGIONS_EXPORT StdSegExp(const StdSegExp &T);

            STD_REGIONS_EXPORT ~StdSegExp();


    protected:


            //----------------------------
            // Integration Methods
            //----------------------------

            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
                    const Array<OneD, const NekDouble>& inarray);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------

            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &out_d0,
                    Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                    Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                    const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &out_d0,
                    Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
                    Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray);

            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv (
                    const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            //-----------------------------
            // Transforms
            //-----------------------------

            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_FwdTrans_BndConstrained(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            //----------------------------
            // Inner product functions
            //----------------------------

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& base,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray,
                    int coll_check);

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> &outarray,
                    bool multiplybyweights = true);

            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase (
                    const int dir,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);



            //----------------------------
            // Evaluations Methods
            //---------------------------

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                    const Array<OneD, const NekDouble>& Lcoords,
                    const Array<OneD, const NekDouble>& physvals);

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,NekDouble> &outarray,
                    const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,NekDouble> &outarray,
                    const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_SVVLaplacianFilter(
                    Array<OneD, NekDouble> &array,
                    const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_ExponentialFilter(
                          Array<OneD, NekDouble> &array,
                    const NekDouble        alpha,
                    const NekDouble        exponent,
                    const NekDouble        cutoff);

            STD_REGIONS_EXPORT virtual void v_MultiplyByStdQuadratureMetric(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_FillMode(
                    const int mode,
                    Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_GetCoords(
                    Array<OneD, NekDouble> &coords_0,
                    Array<OneD, NekDouble> &coords_1,
                    Array<OneD, NekDouble> &coords_2);


            //----------------------------
            // Public Mappings
            //---------------------------

            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                    Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                    Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual int v_GetVertexMap(int localVertexId,
                                                          bool useCoeffPacking = false);

            //----------------------------
            // Helper functions
            //---------------------------

            STD_REGIONS_EXPORT virtual int v_GetNverts() const;
            STD_REGIONS_EXPORT virtual int v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int v_NumDGBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual bool v_IsBoundaryInteriorExpansion();
            STD_REGIONS_EXPORT virtual int v_CalcNumberOfCoefficients(
                    const std::vector<unsigned int> &nummodes,
                    int &modes_offset);
            STD_REGIONS_EXPORT virtual LibUtilities::ShapeType
                v_DetShapeType() const;

            //----------------------------
            // Wrapper functions
            //---------------------------

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(
                    const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                    const StdMatrixKey &mkey);

            //---------------------------------------
            // Output interpolation functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_GetSimplexEquiSpacedConnectivity(
                    Array<OneD, int> &conn,
                    bool standard = true);

            // Operator evaluation functions
            //---------------------------------------

	    STD_REGIONS_EXPORT virtual void v_ReduceOrderCoeffs(
                    int numMin,
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD, NekDouble> &outarray);

        private:

        };


    } //end of namespace
} //end of namespace

#endif //STDSEGEXP_H

