/////////////////////////////////////////////////////////////////////////////////
////
//// File: NekConstantSizedVector.hpp
////
//// For more information, please see: http://www.nektar.info
////
//// The MIT License
////
//// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//// Department of Aeronautics, Imperial College London (UK), and Scientific
//// Computing and Imaging Institute, University of Utah (USA).
////
//// License for the specific language governing rights and limitations under
//// Permission is hereby granted, free of charge, to any person obtaining a
//// copy of this software and associated documentation files (the "Software"),
//// to deal in the Software without restriction, including without limitation
//// the rights to use, copy, modify, merge, publish, distribute, sublicense,
//// and/or sell copies of the Software, and to permit persons to whom the
//// Software is furnished to do so, subject to the following conditions:
////
//// The above copyright notice and this permission notice shall be included
//// in all copies or substantial portions of the Software.
////
//// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//// DEALINGS IN THE SOFTWARE.
////
//// Description: Constant sized and variable sized vectors have the same interface
////              and performt he same algorithms - the only real difference is how
////              the data is stored.  These methods allow use to use the same
////              algorithm regardless of data storage mechanism.
////
/////////////////////////////////////////////////////////////////////////////////

//#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP
//#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP

//#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
//#include <LibUtilities/LinearAlgebra/NekVectorTypeTraits.hpp>
//#include <ExpressionTemplates/ExpressionTemplates.hpp>
//#include <boost/tokenizer.hpp>

//namespace Nektar
//{

    


////    template<typename DataType>
////    void AddEqual(NekVector<DataType>& lhs, const NekVector<DataType>& rhs)
////    {
////        ASSERTL1(lhs.GetDimension() == rhs.GetDimension(), "Two vectors must have the same size in PlusEqual.")
////        DataType* lhs_buf = lhs.GetRawPtr();
////        const DataType* rhs_buf = rhs.GetRawPtr();
////        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
////        {
////            lhs_buf[i] += rhs_buf[i];
////        }
////    }



////    template<typename DataType>
////    void MultiplyEqual(NekVector<DataType>& lhs, const DataType& rhs)
////    {
////        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
////        {
////            lhs[i] *= rhs;
////        }
////    }

////    template<typename DataType>
////    void DivideEqual(NekVector<DataType>& lhs, const DataType& rhs)
////    {
////        for(unsigned int i=0; i < lhs.GetDimension(); ++i)
////        {
////            lhs[i] /= rhs;
////        }
////    }



//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES


//#endif




//}

//#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_VECTOR_COMMON_HPP

