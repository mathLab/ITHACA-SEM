///////////////////////////////////////////////////////////////////////////////
//
// File: NekMatrix.hpp
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
// Description: Generic Matrix
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

#include <loki/Factory.h>
#include <loki/Singleton.h>
#include <loki/Typelist.h>

#include <boost/shared_ptr.hpp>
#include <boost/call_traits.hpp>
#include <boost/concept_check.hpp>

#include <LibUtilities/NekMemoryManager.hpp>
#include <LibUtilities/ErrorUtil.hpp>
#include <LibUtilities/NekVector.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        enum NekMatrixForm
        {
            eFull,
            eDiagonal
//             eZero,
//             eDiagonal,
//             eSquareSymmetric,
//             eSquareSymmetricPositiveDefinite,
//             eSymmetricPositiveDefiniteBanded,
//             eSquareGeneral,
//             eSquareGeneralBanded
        };

        class OutOfBoundsError
        {
        };

        /// NekMatrix Class
        /// \param width the width of the matrix.  Setting this to 0 allows for dynamically sized matrices.
        /// \param height The height of the matrix.  Setting this to 0 allows for dynamically size matrices.
        ///
        /// Some design decisions about the class.
        /// Only one factory object.  All impl objects need to implement the
        /// Initialize methods.  I could have made many constructors for the
        /// impl object constructors, but then I would have had a lot of messy
        /// Factory objects to created, one for each constructor.
        template<typename DataType, unsigned int space = 0>
        class NekMatrix
        {

            public:
                NekMatrix(unsigned int rows, unsigned int columns, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(rows),
                    m_columns(columns),
                    m_data()
                {
                    switch(form)
                    {
                        case eFull:
                        {
                            m_data = new DataType[rows*columns];//boost::shared_array<DataType>(new DataType[rows*columns]);
                        }
                        break;

                        case eDiagonal:
                        {
                            m_data = new DataType[rows]; //boost::shared_array<DataType>(new DataType[rows]);
                        }
                        break;
                    }
                }

                NekMatrix(unsigned int rows, unsigned int columns, const DataType* const ptr,
                        NekMatrixForm form = eFull) :
                        m_form(form),
                        m_rows(rows),
                        m_columns(columns),
                        m_data()
                {
                    switch(form)
                    {
                        case eFull:
                        {
                            m_data = new DataType[rows*columns]; //boost::shared_array<DataType>(new DataType[rows*columns]);
                            std::copy(ptr, ptr+rows*columns, begin());
                        }
                        break;

                        case eDiagonal:
                        {
                            m_data = new DataType[rows]; //boost::shared_array<DataType>(new DataType[rows]);
                            std::copy(ptr, ptr+rows, begin());
                        }
                        break;
                    }
                }

                NekMatrix(unsigned int rows, unsigned int columns, typename boost::call_traits<DataType>::param_type d,
                        NekMatrixForm form = eFull) :
                        m_form(form),
                        m_rows(rows),
                        m_columns(columns),
                        m_data()
                {
                    switch(form)
                    {
                        case eFull:
                        {
                            m_data = new DataType[rows*columns]; //boost::shared_array<DataType>(new DataType[rows*columns]);
                        }
                        break;

                        case eDiagonal:
                        {
                            m_data = new DataType[rows]; //boost::shared_array<DataType>(new DataType[rows]);
                        }
                        break;
                    }

                    std::fill(begin(), end(), d);
                }

                unsigned int rows() const { return m_rows; }
                unsigned int columns() const { return m_columns; }

                inline typename boost::call_traits<DataType>::reference operator()
                    (unsigned int rowNumber, unsigned int colNumber)
                {
                    if( rowNumber >= m_rows || colNumber >= m_columns )
                    {
                        throw OutOfBoundsError();
                    }

                    switch(m_form)
                    {
                        case eFull:
                        {
                            return m_data[rowNumber*m_columns + colNumber];
                        }
                        break;

                        case eDiagonal:
                        {
                            return m_data[rowNumber];
                        }
                        break;

                        default:
                        {
                            static DataType result;
                            return result;
                        }
                    }
                }

                inline typename boost::call_traits<DataType>::const_reference operator()
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    if( rowNumber >= m_rows || colNumber >= m_columns )
                    {
                        throw OutOfBoundsError();
                    }

                    switch(m_form)
                    {
                        case eFull:
                        {
                            return m_data[rowNumber*m_columns + colNumber];
                        }
                        break;

                        case eDiagonal:
                        {
                            return m_data[rowNumber];
                        }
                        break;

                        default:
                        {
                            static DataType result;
                            return result;
                        }
                    }
                }


                typename boost::call_traits<DataType>::reference getValue
                        (unsigned int rowNumber, unsigned int colNumber)
                {
                    return (*this)(rowNumber, colNumber);
                }

                typename boost::call_traits<DataType>::const_reference getValue
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    return (*this)(rowNumber, colNumber);
                }

                void setValue(unsigned int rowNumber, unsigned int colNumber, typename boost::call_traits<DataType>::param_type rhs)
                {
                    (*this)(rowNumber, colNumber) = rhs;
                }

                DataType* getPtr(unsigned int rowNumber, unsigned int colNumber)
                {
                    switch(m_form)
                    {
                        case eFull:
                        {
                            return &m_data[rowNumber*m_columns + colNumber];
                        }
                        break;

                        case eDiagonal:
                        {
                            return &m_data[rowNumber];
                        }
                        break;

                        default:
                        {
                            return NULL;
                        }
                    }
                }

                DataType* begin()
                {
                    return &m_data[0];
                }

                DataType* end()
                {
                    switch(m_form)
                    {
                        case eFull:
                            return &m_data[m_rows*m_columns];
                            break;

                        case eDiagonal:
                            return &m_data[m_rows];
                            break;

                        default:
                            return NULL;
                    }
                }

                const DataType* begin() const
                {
                    return &m_data[0];
                }

                const DataType* end() const
                {
                    switch(m_form)
                    {
                        case eFull:
                            return &m_data[m_rows*m_columns];
                            break;

                        case eDiagonal:
                            return &m_data[m_rows];
                            break;

                        default:
                            return NULL;
                    }
                }

                NekMatrix<DataType, space> operator+=(const NekMatrix<DataType, space>& rhs)
                {
                    ASSERTL0(rows() == rhs.rows() && columns() == rhs.columns(), "Matrix dimensions must agree in operator+");
//                     for(unsigned int i = 0; i < rows(); ++i)
//                     {
//                         for(unsigned int j = 0; j < columns(); ++j)
//                         {
//                             (*this)(i,j) += rhs(i,j);
//                         }
//                     }
                    DataType* lhs_data = begin();
                    const DataType* rhs_data = rhs.begin();

                    for(unsigned int i = 0; i < m_rows*m_columns; ++i)
                    {
                        lhs_data[i] += rhs_data[i];
                    }
                    return *this;
                }

            private:
                NekMatrixForm m_form;
                unsigned int m_rows;
                unsigned int m_columns;
                //boost::shared_array<DataType> m_data;
                DataType* m_data;
        };


    template<typename DataType, unsigned int space>
    NekMatrix<DataType, space> operator+(
            const NekMatrix<DataType, space>& lhs,
            const NekMatrix<DataType, space>& rhs)
        {
            NekMatrix<DataType, space> result(lhs);
            result += rhs;
            return result;
        }


        template<typename DataType, unsigned int space>
        NekMatrix<DataType, space> operator*(
        const NekMatrix<DataType, space>& lhs,
        const NekMatrix<DataType, space>& rhs)
        {
            ASSERTL0(lhs.columns() == rhs.rows(), "Invalid matrix dimensions in operator*");

            NekMatrix<DataType, space> result(lhs.rows(), rhs.columns());

            for(unsigned int i = 0; i < result.rows(); ++i)
            {
                for(unsigned int j = 0; j < result.columns(); ++j)
                {
                DataType t = DataType(0);

                // Set the result(i,j) element.
                for(unsigned int k = 0; k < lhs.columns(); ++k)
                {
                t += lhs(i,k)*rhs(k,j);
                }
                result(i,j) = t;
                }
            }

            return result;
        }

        template<typename DataType, unsigned int space>
        NekVector<DataType, 0, space> operator*(
        const NekMatrix<DataType, space>& lhs,
        const NekVector<DataType, 0, space>& rhs)
        {
            ASSERTL0(lhs.columns() == rhs.dimension(), "Invalid matrix dimensions in operator*");

            NekVector<DataType, 0, space> result(rhs.dimension(), DataType(0));

            for(unsigned int i = 0; i < lhs.columns(); ++i)
            {
                DataType t = DataType(0);
                for(unsigned int j = 0; j < rhs.dimension(); ++j)
                {
                    t += lhs(i,j)*rhs(j);
                }
                result(i) = t;
            }

            return result;
        }

    }
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekMatrix.hpp,v $
    Revision 1.7  2006/05/29 03:40:12  bnelson
    Changed the implementation from individual NekMatrixImpl objects for each type of matrix to an enumeration to store the type.

    Revision 1.6  2006/05/25 03:02:40  bnelson
    Added Matrix/Vector multiplication.

    Revision 1.5  2006/05/18 04:21:06  bnelson
    Removed the ability to specify the rows and columns in the template parameter list.  If this behavior is desired we'll need to create a fixed size array class.

    Added multiplication to the arrays.

    Revision 1.4  2006/05/15 05:06:55  bnelson
    Removed use of MemoryManager pending review of some memory problems.

    Revision 1.3  2006/05/15 04:13:36  bnelson
    no message

    Revision 1.2  2006/05/14 21:32:03  bnelson
 *** empty log message ***

    Revision 1.1  2006/05/04 18:57:43  kirby
 *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


 **/


