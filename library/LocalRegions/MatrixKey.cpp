///////////////////////////////////////////////////////////////////////////////
//
// File MatrixKey.cpp
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
// Description: Definition of MatrixKey based on StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/Expansion.h>
#include <LocalRegions/MatrixKey.h>

namespace Nektar
{
    namespace LocalRegions
    {
        MatrixKey::MatrixKey(const StdRegions::MatrixType matrixType,
                  const LibUtilities::ShapeType shapeType,
                  const StdRegions::StdExpansion &stdExpansion,
                  const StdRegions::ConstFactorMap &factorMap,
                  const StdRegions::VarCoeffMap &varCoeffMap,
                  LibUtilities::PointsType nodalType) :
            StdMatrixKey(matrixType, shapeType, stdExpansion, factorMap, varCoeffMap, nodalType),
            m_metricinfo( ( dynamic_cast<const Expansion&>( stdExpansion ) ).GetMetricInfo() )
        {
        }

        MatrixKey::MatrixKey(const MatrixKey& mkey,
                      const StdRegions::MatrixType matrixType) :
            StdRegions::StdMatrixKey(mkey, matrixType),
            m_metricinfo(mkey.m_metricinfo)
        {
        }

        MatrixKey::MatrixKey(const StdRegions::StdMatrixKey &mkey) :
            StdRegions::StdMatrixKey(mkey)
        {
        }

        bool MatrixKey::opLess::operator()(const MatrixKey &lhs, const MatrixKey &rhs) const
        {        
            {
                return (lhs.GetMatrixType() < rhs.GetMatrixType());
            }
        }

        bool operator<(const MatrixKey &lhs, const MatrixKey &rhs)
        {
            if(lhs.m_metricinfo.get() < rhs.m_metricinfo.get())
            {
                return true;
            }


            if(lhs.m_metricinfo.get() > rhs.m_metricinfo.get())
            {
                return false;
            }    
            
            return (*dynamic_cast<const StdRegions::StdMatrixKey*>(&lhs)
                    < *dynamic_cast<const StdRegions::StdMatrixKey*>(&rhs));
        }

    }
}
