///////////////////////////////////////////////////////////////////////////////
//
// File: NekConstantSizedVector.hpp
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
// Description: Constant sized and variable sized vectors have the same interface
//              and performt he same algorithms - the only real difference is how
//              the data is stored.  These methods allow use to use the same 
//              algorithm regardless of data storage mechanism.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP

#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorTypeTraits.hpp>

namespace Nektar
{
    template<typename DataType>
    std::vector<DataType> fromString(const std::string& str)
    {
        std::vector<DataType> result;

        try
        {
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep("(<,>) ");
            tokenizer tokens(str, sep);
            for( tokenizer::iterator strIter = tokens.begin(); strIter != tokens.end(); ++strIter)
            {
                result.push_back(boost::lexical_cast<DataType>(*strIter));
            }
        }

        return result;
    }

    /// \todo Do the Norms with Blas where applicable.
    template<typename DataType, unsigned int dim, unsigned int space>
    DataType L1Norm(const NekVector<DataType, dim, space>& v)
    {
        typedef NekVector<DataType, dim, space> VectorType;

        DataType result(0);
        for(VectorType::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        {
            result += NekVectorTypeTraits<DataType>::abs(*iter);
        }

        return result;
    }
    
    template<typename DataType, unsigned int dim, unsigned int space>
    DataType L2Norm(const NekVector<DataType, dim, space>& v)
    {
        typedef NekVector<DataType, dim, space> VectorType;

        DataType result(0);
        for(VectorType::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        {
            DataType v = NekVectorTypeTraits<DataType>::abs(*iter);
            result += v*v;
        }
        return sqrt(result);
    }
    
    template<typename DataType, unsigned int dim, unsigned int space>
    DataType InfinityNorm(const NekVector<DataType, dim, space>& v) 
    {
        DataType result = NekVectorTypeTraits<DataType>::abs(v[0]);
        for(unsigned int i = 1; i < v.GetDimension(); ++i)
        {
            result = std::max(NekVectorTypeTraits<DataType>::abs(v[i]), result);
        }
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space> Negate(const NekVector<DataType, dim, space>& v) 
    {
        NekVector<DataType, dim, space> temp(v);
        for(unsigned int i=0; i < temp.GetDimension(); ++i)
        {
            temp(i) = -temp(i);
        }
        return temp;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    void PlusEqual(NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetDimension() == rhs.GetDimension(), "Two vectors must have the same size in PlusEqual.")
        
        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
        {
            lhs[i] += rhs[i];
        }
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    void MinusEqual(NekVector<DataType, dim, space>& lhs, const NekVector<DataType, dim, space>& rhs)
    {
        ASSERTL1(lhs.GetDimension() == rhs.GetDimension(), "Two vectors must have the same size in MinusEqual.")
        
        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
        {
            lhs[i] -= rhs[i];
        }
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    void TimesEqual(NekVector<DataType, dim, space>& lhs, const DataType& rhs)
    {
        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
        {
            lhs[i] *= rhs;
        }
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    void DivideEqual(NekVector<DataType, dim, space>& lhs, const DataType& rhs)
    {
        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
        {
            lhs[i] /= rhs;
        }
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    DataType Magnitude(const NekVector<DataType, dim, space>& v) 
    {
        DataType result = DataType(0);

        for(unsigned int i = 0; i < v.GetDimension(); ++i)
        {
            result += v[i]*v[i];
        }
        return sqrt(result);
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    DataType Dot(const NekVector<DataType, dim, space>& lhs, 
                 const NekVector<DataType, dim, space>& rhs) 
    {
        ASSERTL1( lhs.GetDimension() == rhs.GetDimension(), "Dot, dimension of the two operands must be identical.");

        DataType result = DataType(0);
        for(unsigned int i = 0; i < lhs.GetDimension(); ++i)
        {
            result += lhs[i]*rhs[i];
        }

        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    void Normalize(NekVector<DataType, dim, space>& v)
    {
        DataType m = v.Magnitude();
        if( m > DataType(0) )
        {
            v /= m;
        }
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    NekVector<DataType, dim, space> Cross(const NekVector<DataType, dim, space>& lhs,
                                          const NekVector<DataType, dim, space>& rhs)
    {
        BOOST_STATIC_ASSERT(dim==3 || dim == 0);
        ASSERTL1(lhs.GetDimension() == 3 && rhs.GetDimension() == 3, "Cross is only valid for 3D vectors.");

        DataType first = lhs.y()*rhs.z() - lhs.z()*rhs.y();
        DataType second = lhs.z()*rhs.x() - lhs.x()*rhs.z();
        DataType third = lhs.x()*rhs.y() - lhs.y()*rhs.x();

        NekVector<DataType, dim, space> result(first, second, third);
        return result;
    }

    template<typename DataType, unsigned int dim, unsigned int space>
    std::string AsString(const NekVector<DataType, dim, space>& v)
    {
        unsigned int d = v.GetRows();
        std::string result = "(";
        for(unsigned int i = 0; i < d; ++i)
        {
            result += boost::lexical_cast<std::string>(v[i]);
            if( i < dim-1 )
            {
                result += ", ";
            }
        }
        result += ")";
        return result;
    }

}

#endif NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP

