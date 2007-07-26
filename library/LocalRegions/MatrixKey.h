///////////////////////////////////////////////////////////////////////////////
//
// File MatrixKeys.h
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
// Description: Headers for MatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MATRIXKEY_H
#define MATRIXKEY_H

#include <LocalRegions/LocalRegions.hpp>
#include <StdRegions/StdMatrixKey.h>
#include <SpatialDomains/GeomFactors.h>

namespace Nektar
{
    namespace LocalRegions
    {

        class MatrixKey
        {
        public:
            MatrixKey( StdRegions::MatrixType matrixType, 
                       StdRegions::ShapeType shapeType, 
                       StdRegions::StdExpansion &stdExpansion,
                       double scalefactor = 1.0,
                       LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            virtual ~MatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                bool operator()(const MatrixKey &lhs, const MatrixKey &rhs) const;
            };

            /// Used for finding value given the key in NekManager.
            friend bool operator<(const MatrixKey &lhs, const MatrixKey &rhs);
            friend bool opLess::operator()(const MatrixKey &lhs, 
                const MatrixKey &rhs) const;

            StdRegions::MatrixType GetMatrixType() const
            {
                return m_stdMatKey->GetMatrixType(); 
            }

            StdRegions::ShapeType GetShapeType() const
            {
                return m_stdMatKey->GetShapeType();
            }

            int GetNcoeffs() const
            {
                return m_stdMatKey->GetNcoeffs();
            }

            const ConstArray<OneD, LibUtilities::BasisSharedPtr> &GetBase()
            {
                return m_stdMatKey->GetBase();
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(int dir) const
            {
                return m_stdMatKey->GetBasis(dir);
            }

            const StdRegions::StdMatrixKeySharedPtr &GetStdMatKey() const 
            {
                return m_stdMatKey;
            }

            SpatialDomains::GeomFactorsSharedPtr GetMetricInfo() const
            {
                return m_metricinfo;
            }

            const double GetScaleFactor() const 
            {
                return m_scalefactor;
            }

        protected:
            MatrixKey();

            StdRegions::StdMatrixKeySharedPtr     m_stdMatKey;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo; 
            NekDouble                             m_scalefactor;

        private:
        };

        std::ostream& operator<<(std::ostream& os, const MatrixKey& rhs);

    } // end of namespace
} // end of namespace

#endif //STDMATRIXKEY_H

/**
* $Log: MatrixKey.h,v $
* Revision 1.11  2007/07/12 12:53:00  sherwin
* Updated to have a helmholtz matrix
*
* Revision 1.10  2007/07/10 17:17:22  sherwin
* Introduced Scaled Matrices into the MatrixManager
*
* Revision 1.9  2007/05/28 08:35:25  sherwin
* Updated for localregions up to Project1D
*
* Revision 1.8  2007/05/27 16:10:28  bnelson
* Update to new Array type.
*
* Revision 1.7  2007/04/08 03:33:30  jfrazier
* Minor reformatting and fixing SharedArray usage.
*
* Revision 1.6  2007/04/04 21:49:23  sherwin
* Update for SharedArray
*
* Revision 1.5  2007/03/29 19:02:05  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.4  2007/03/25 15:48:21  sherwin
* UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
* Revision 1.3  2007/03/20 09:13:37  kirby
* new geomfactor routines; update for metricinfo; update style
*
* Revision 1.2  2007/03/14 21:24:07  sherwin
* Update for working version of MultiRegions up to ExpList1D
*
* Revision 1.1  2007/03/04 20:42:14  sherwin
* Keys for matrix managers
*
* Revision 1.1  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
*
***/
