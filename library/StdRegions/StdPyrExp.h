///////////////////////////////////////////////////////////////////////////////
//
// File StdPyrExp.h
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
// Description: Header field for pyramidic routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGIONS_STDPYREXP_H
#define NEKTAR_LIBS_STDREGIONS_STDPYREXP_H

#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>
#include <boost/tuple/tuple.hpp>

namespace Nektar
{
    namespace StdRegions
    {
        class StdPyrExp : virtual public StdExpansion3D
        {
        public:
            typedef boost::tuple<unsigned char, unsigned char, unsigned char, bool> triple;
    
            STD_REGIONS_EXPORT StdPyrExp();

            STD_REGIONS_EXPORT StdPyrExp(const LibUtilities::BasisKey &Ba, 
                                         const LibUtilities::BasisKey &Bb, 
                                         const LibUtilities::BasisKey &Bc);

            STD_REGIONS_EXPORT StdPyrExp(const LibUtilities::BasisKey &Ba, 
                                         const LibUtilities::BasisKey &Bb, 
                                         const LibUtilities::BasisKey &Bc,
                                         NekDouble *coeffs, 
                                         NekDouble *phys);

            STD_REGIONS_EXPORT StdPyrExp(const StdPyrExp &T);

            STD_REGIONS_EXPORT ~StdPyrExp();

            //---------------------------------------
            // Integration/public 3D methods
            //---------------------------------------
            STD_REGIONS_EXPORT void TripleTensorProduct(
                const Array<OneD, const NekDouble>& fx, 
                const Array<OneD, const NekDouble>& gy, 
                const Array<OneD, const NekDouble>& hz, 
                const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> & outarray);

            STD_REGIONS_EXPORT NekDouble TripleInnerProduct(
                const Array<OneD, const NekDouble>& fxyz,
                const Array<OneD, const NekDouble>& wx,
                const Array<OneD, const NekDouble>& wy,
                const Array<OneD, const NekDouble>& wz);

            STD_REGIONS_EXPORT NekDouble Integral3D(
                const Array<OneD, const NekDouble>& inarray,
                const Array<OneD, const NekDouble>& wx,
                const Array<OneD, const NekDouble>& wy,
                const Array<OneD, const NekDouble>& wz);

            STD_REGIONS_EXPORT std::vector<triple> &GetMap()
            {
                return m_map;
            }

            STD_REGIONS_EXPORT vector<int> &GetRMap()
            {
                return m_rmap;
            }
            STD_REGIONS_EXPORT int GetTetMode(int I, int J, int K);
            STD_REGIONS_EXPORT int GetModTetMode(int I, int J, int K);

        protected:
            //---------------------------------------
            // Differentiation/integration Methods
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &out_d0,
                      Array<OneD,       NekDouble> &out_d1,
                      Array<OneD,       NekDouble> &out_d2);
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &out_d0,
                      Array<OneD,       NekDouble> &out_d1,
                      Array<OneD,       NekDouble> &out_d2);
            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
                const Array<OneD, const NekDouble>& inarray);
            STD_REGIONS_EXPORT void v_MultiplyByStdQuadratureMetric(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            //---------------------------------------
            // Transforms
            //---------------------------------------

            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2);

            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);


            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2);

            /*
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray, 
                      Array<OneD,       NekDouble> &outarray); 
            */

            //---------------------------------------
            // Evaluation functions
            //---------------------------------------

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(
                const Array<OneD, const NekDouble>& xi,
                const Array<OneD, const NekDouble>& physvals);

            STD_REGIONS_EXPORT virtual void v_GetCoords(
                Array<OneD, NekDouble> & xi_x, 
                Array<OneD, NekDouble> & xi_y, 
                Array<OneD, NekDouble> & xi_z);

            STD_REGIONS_EXPORT virtual void v_FillMode(
                const int mode, 
                Array<OneD, NekDouble> &outarray);  

            //---------------------------------------
            // Helper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual int v_GetNverts() const;
            STD_REGIONS_EXPORT virtual int v_GetNedges() const;
            STD_REGIONS_EXPORT virtual int v_GetNfaces() const;
            STD_REGIONS_EXPORT virtual LibUtilities::ShapeType v_DetShapeType() const;
            STD_REGIONS_EXPORT virtual int v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int v_GetEdgeNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetFaceNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_GetFaceIntNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int v_CalcNumberOfCoefficients(
                const std::vector<unsigned int> &nummodes, 
                int &modes_offset);
            STD_REGIONS_EXPORT virtual LibUtilities::BasisType v_GetEdgeBasisType(
                const int i) const;

            //---------------------------------------
            // Mappings
            //---------------------------------------
            STD_REGIONS_EXPORT virtual void v_GetFaceToElementMap(
                const int                  fid,
                const Orientation      faceOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD,          int> &signarray,
                int                        nummodesA=-1,
                int                        nummodesB=-1);
            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(int localVertexId,
                                                           bool useCoeffPacking = false);
            /*
            STD_REGIONS_EXPORT virtual void v_GetEdgeInteriorMap(
                const int eid,
                const Orientation edgeOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD, int> &signarray);
            STD_REGIONS_EXPORT virtual void v_GetFaceInteriorMap(
                const int fid,
                const Orientation faceOrient,
                Array<OneD, unsigned int> &maparray,
                Array<OneD, int>& signarray);
            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(
                Array<OneD, unsigned int> &outarray);
            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(
                Array<OneD, unsigned int>& outarray);
            */

            //---------------------------------------
            // Wrapper functions
            //---------------------------------------
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey);

        private:
            //---------------------------------------
            // Private helper functions
            //---------------------------------------
            vector<triple> m_map;
            vector<int> m_rmap;
            LibUtilities::BasisSharedPtr m_base_A;
            LibUtilities::BasisSharedPtr m_base_C;
        };    
        typedef boost::shared_ptr<StdPyrExp> StdPyrExpSharedPtr;
    } //end of namespace
} //end of namespace

#endif //STDPYREXP_H
