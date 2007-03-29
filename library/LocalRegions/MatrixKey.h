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
		       StdRegions::StdExpansion &stdExpansion);

            virtual ~MatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                bool operator()(const MatrixKey &lhs, const MatrixKey &rhs);
            };

            /// Used for finding value given the key in NekManager.
            friend bool operator<(const MatrixKey &lhs, const MatrixKey &rhs);
            friend bool opLess::operator()(const MatrixKey &lhs, 
					   const MatrixKey &rhs);

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

            SharedArray<LibUtilities::BasisSharedPtr> GetBase() const
            {
                return m_stdMatKey->GetBase();
            }

            StdRegions::StdMatrixKeySharedPtr GetStdMatKey() const 
		{
                return m_stdMatKey;
            }

            SpatialDomains::GeomFactorsSharedPtr GetMetricInfo() const
            {
                return m_metricinfo;
            }

        protected:
            MatrixKey();

            StdRegions::StdMatrixKeySharedPtr  m_stdMatKey;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo; 

        private:
        };

        std::ostream& operator<<(std::ostream& os, const MatrixKey& rhs);

    } // end of namespace
} // end of namespace

#endif //STDMATRIXKEY_H

/**
* $Log: MatrixKey.h,v $
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
