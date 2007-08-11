///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrixKeys.h
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
// Description: Headers for StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDMATRIXKEY_H
#define STDMATRIXKEY_H

#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion;

        class StdMatrixKey
        {
        public:
            StdMatrixKey( const StdRegions::MatrixType matrixType, 
                          const StdRegions::ShapeType shapeType, 
                          const StdRegions::StdExpansion &stdExpansion,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);
            
            StdMatrixKey(const MatrixType matrixType, 
                         const ShapeType shapeType,
                         const ConstArray<OneD,LibUtilities::BasisSharedPtr> &base,
                         const int ncoeffs,
                         LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);
            
            StdMatrixKey(const StdMatrixKey& rhs);

            virtual ~StdMatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                bool operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;
            };

            /// Used for finding value given the key in NekManager.
            friend bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            friend bool opLess::operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;

            MatrixType GetMatrixType() const
            {
                return m_matrixType;
            }

            ShapeType GetShapeType() const
            {
                return m_shapeType;
            }

            LibUtilities::PointsType GetNodalPointsType() const
            {
                return m_nodalPointsType;
            }
           
            int GetNcoeffs() const
            {
                return m_ncoeffs;
            }

            inline const ConstArray<OneD, LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return m_base;
            }


            inline const LibUtilities::BasisSharedPtr GetBasis(int dir) const
            {
                return(m_base[dir]);
            }

        protected:
            ShapeType m_shapeType;
            ConstArray<OneD, LibUtilities::BasisSharedPtr> m_base;

            unsigned int m_ncoeffs;
            MatrixType m_matrixType;
            LibUtilities::PointsType m_nodalPointsType;

        private:
            StdMatrixKey();
        };

        std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs);

        typedef  boost::shared_ptr<StdMatrixKey> StdMatrixKeySharedPtr;

    } // end of namespace
} // end of namespace

#endif //STDMATRIXKEY_H

/**
* $Log: StdMatrixKey.h,v $
* Revision 1.16  2007/07/26 02:39:21  bnelson
* Fixed Visual C++ compiler errors when compiling in release mode.
*
* Revision 1.15  2007/07/22 23:04:26  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.14  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.13  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.12  2007/05/25 17:48:58  jfrazier
* Changed the class so the bases are not stored as an Array, but rather a ConstArray.
*
* Revision 1.11  2007/05/22 02:01:50  bnelson
* Changed Array::size to Array::num_elements.
*
* Fixed some compiler errors in assertions.
*
* Revision 1.10  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.9  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.8  2007/04/08 03:36:58  jfrazier
* Updated to use SharedArray consistently and minor reformatting.
*
* Revision 1.7  2007/04/04 21:49:25  sherwin
* Update for SharedArray
*
* Revision 1.6  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.5  2007/03/21 20:56:43  sherwin
* Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv
*
* Revision 1.4  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.3  2007/03/05 08:07:14  sherwin
* Modified so that StdMatrixKey has const calling arguments in its constructor.
*
* Revision 1.2  2007/03/02 12:01:53  sherwin
* Update for working version of LocalRegions/Project1D
*
* Revision 1.1  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
*
***/
