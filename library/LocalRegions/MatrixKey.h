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

#include <StdRegions/StdMatrixKey.h>
#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

namespace Nektar
{
    namespace LocalRegions
    {

        class MatrixKey : public StdRegions::StdMatrixKey
        {
        public:
            LOCAL_REGIONS_EXPORT MatrixKey(const StdRegions::MatrixType matrixType,
                      const LibUtilities::ShapeType shapeType,
                      const StdRegions::StdExpansion &stdExpansion,
                      const StdRegions::ConstFactorMap &factorMap = StdRegions::NullConstFactorMap,
                      const StdRegions::VarCoeffMap &varCoeffMap = StdRegions::NullVarCoeffMap,
                      LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            LOCAL_REGIONS_EXPORT MatrixKey(const MatrixKey& mkey,
                          const StdRegions::MatrixType matrixType);
            
            LOCAL_REGIONS_EXPORT MatrixKey(const StdRegions::StdMatrixKey &mkey);

            virtual ~MatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                LOCAL_REGIONS_EXPORT bool operator()(const MatrixKey &lhs, const MatrixKey &rhs) const;
            };

            /// Used for finding value given the key in NekManager.
            LOCAL_REGIONS_EXPORT friend bool operator<(const MatrixKey &lhs, const MatrixKey &rhs);
            LOCAL_REGIONS_EXPORT friend bool opLess::operator()(const MatrixKey &lhs, 
                const MatrixKey &rhs) const;

            SpatialDomains::GeomFactorsSharedPtr GetMetricInfo() const
            {
                return m_metricinfo;
            }

        protected:
            MatrixKey();

            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo; 

        private:
        };

    } // end of namespace
} // end of namespace

#endif //STDMATRIXKEY_H

