///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriExp.h
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
// Description: Header for NodalTriExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIEXP_H
#define NODALTRIEXP_H

#include <StdRegions/StdNodalTriExp.h>
#include <SpatialDomains/TriGeom.h>

#include <LocalRegions/MatrixKey.h>
#include <LocalRegions/LocalRegionsDeclspec.h>
#include <LocalRegions/Expansion2D.h>

namespace Nektar
{
    namespace LocalRegions 
    {  
        
        class NodalTriExp: virtual public StdRegions::StdNodalTriExp,
                           virtual public Expansion2D
        {
        public:
            /** \brief Constructor using BasisKey class for quadrature
                points and order definition */
            LOCAL_REGIONS_EXPORT NodalTriExp(
                const LibUtilities::BasisKey           &Ba,
                const LibUtilities::BasisKey           &Bb,
                const LibUtilities::PointsType          Ntype,
                const SpatialDomains::TriGeomSharedPtr &geom);
            
            /// Copy Constructor
            LOCAL_REGIONS_EXPORT NodalTriExp(const NodalTriExp &T); 
            
            /// Destructor
            LOCAL_REGIONS_EXPORT ~NodalTriExp();

        protected:
            //----------------------------
            // Integration Methods
            //----------------------------
            virtual NekDouble v_Integral(
                const Array<OneD, const NekDouble> &inarray);
            
            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            virtual void v_PhysDeriv(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &out_d0,
                      Array<OneD,       NekDouble> &out_d1,
                      Array<OneD,       NekDouble> &out_d2
                                                      = NullNekDouble1DArray);

            //---------------------------------------
            // Transforms
            //---------------------------------------
            virtual void v_FwdTrans(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            //---------------------------------------
            // Inner product functions
            //---------------------------------------
            virtual void v_IProductWRTBase(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_IProductWRTDerivBase(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_IProductWRTBase_SumFac(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_IProductWRTBase_MatOp(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_IProductWRTDerivBase_SumFac(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);
            virtual void v_IProductWRTDerivBase_MatOp(
                const int                           dir,
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray);

            //---------------------------------------
            // Evaluation functions
            //---------------------------------------
            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(
                            const Array<OneD, const NekDouble>& Lcoords,
                                  Array<OneD,NekDouble> &coords);
            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(
                                  Array<OneD,       NekDouble> &coords_1,
                                  Array<OneD,       NekDouble> &coords_2,
                                  Array<OneD,       NekDouble> &coords_3);
            
            //---------------------------------------
            // Matrix creation functions
            //---------------------------------------
            virtual DNekMatSharedPtr v_GenMatrix(
                const StdRegions::StdMatrixKey &mkey);

        private:           
            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            virtual DNekMatSharedPtr v_GenNBasisTransMatrix()
            {
                return StdNodalTriExp::GenNBasisTransMatrix();
            }
            
            virtual void v_GetNodalPoints(Array<OneD, const NekDouble> &x, 
                                          Array<OneD, const NekDouble> &y)
            {
                return StdNodalTriExp::GetNodalPoints(x,y);
            }
            
            
            virtual DNekScalMatSharedPtr v_GetLocMatrix(const MatrixKey &mkey)
            {
                return m_matrixManager[mkey];
            }
            
            virtual DNekScalBlkMatSharedPtr v_GetLocStaticCondMatrix(const MatrixKey &mkey)
            {
                return m_staticCondMatrixManager[mkey];
            }

            virtual DNekMatSharedPtr CreateStdMatrix(
                const StdRegions::StdMatrixKey &mkey);
            virtual DNekScalMatSharedPtr CreateMatrix(
                const MatrixKey &mkey);
            virtual DNekScalBlkMatSharedPtr CreateStaticCondMatrix(
                const MatrixKey &mkey);

            void v_ComputeEdgeNormal(const int edge);
        };
    
        // type defines for use of TriExp in a boost vector
        typedef boost::shared_ptr<NodalTriExp> NodalTriExpSharedPtr;
        typedef std::vector< NodalTriExpSharedPtr > NodalTriExpVector;
        typedef std::vector< NodalTriExpSharedPtr >::iterator NodalTriExpVectorIter;
    
    } //end of namespace
} //end of namespace

#endif // NODALTRIEXP_H

