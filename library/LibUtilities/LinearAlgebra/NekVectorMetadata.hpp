///////////////////////////////////////////////////////////////////////////////
//
// File: NekVectorMetadata.hpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_VECTOR_METADATA_HPP
#define NEKTAR_LIB_UTILITIES_NEK_VECTOR_METADATA_HPP

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>

namespace Nektar
{
    // Interface for expression templates.
    class NekVectorMetadata
    {
        public:
            template<typename VectorType>
            explicit NekVectorMetadata(const VectorType& vec) :
                Rows(vec.GetDimension())
            {
            }

            NekVectorMetadata(const NekVectorMetadata& rhs) :
                Rows(rhs.Rows)
            {
            }

            static NekVectorMetadata CreateForNegation(const NekVectorMetadata& rhs)
            {
                return NekVectorMetadata(rhs);
            }

            static NekVectorMetadata CreateForAddition(const NekVectorMetadata& lhs, const NekVectorMetadata& rhs)
            {
                ASSERTL1(lhs.Rows == rhs.Rows, "Vector dimensions must agree in operator+");
                return NekVectorMetadata(lhs);
            }

            static NekVectorMetadata CreateForMultiplication(const NekMatrixMetadata& lhs, const NekVectorMetadata& rhs)
            {
                ASSERTL1(lhs.Columns == rhs.Rows, "Matrix dimensions must agree in operator*");
                NekVectorMetadata result;
                result.Rows = lhs.Rows;
                return result;
            }


            NekVectorMetadata& operator=(const NekVectorMetadata& rhs)
            {
                Rows = rhs.Rows;
                return *this;
            }

            unsigned int Rows;

        private:
            NekVectorMetadata() :
                Rows(0)
            {
            }
    };
}

#endif //NEKTAR_LIB_UTILITIES_NEK_VECTOR_METADATA_HPP

/**
    $Log: NekVectorMetadata.hpp,v $
    Revision 1.1  2006/09/14 02:06:17  bnelson
    Fixed gcc compiler errors.

 **/


