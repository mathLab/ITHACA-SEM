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
            eFull//,
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

        // Used to force the NekMatrix factory to register.
        template<typename MatrixType>
        class NekMatrixFunctionRegistration
        {
            public:
                void constraints()
                {
                    MatrixType* m;
                    m->RegisterCreateFuncs();
                }
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
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, rows, columns))
                {
                    m_impl->Initialize();
                }

                NekMatrix(const DataType* const ptr, unsigned int rows, unsigned int columns,
                          NekMatrixForm form = eFull) :
                    m_form(form),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, rows, columns))
                {
                    m_impl->Initialize(ptr);
                }

                NekMatrix(typename boost::call_traits<DataType>::param_type d,
                          unsigned int rows, unsigned int columns, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, rows, columns))
                {
                    m_impl->Initialize(d);
                }

                unsigned int rows() const { return m_impl->rows(); }
                unsigned int columns() const { return m_impl->columns(); }

                typename boost::call_traits<DataType>::reference operator()
                        (unsigned int rowNumber, unsigned int colNumber)
                {
                    if( rowNumber >= rows() || colNumber >= columns() )
                    {
                        throw OutOfBoundsError();
                    }

                    return (*m_impl)(rowNumber, colNumber);
                }

                typename boost::call_traits<DataType>::const_reference operator()
                        (unsigned int rowNumber, unsigned int colNumber) const
                {
                    assert(fullRegistered);
                    if( rowNumber >= rows() || colNumber >= columns() )
                    {
                        throw OutOfBoundsError();
                    }

                    return (*m_impl)(rowNumber, colNumber);
                }

                bool RegisterCreateFuncs()
                {
                    return fullRegistered;
                }

                typedef NekMatrix<DataType, space> MyType;
                BOOST_CLASS_REQUIRE(MyType, LibUtilities, NekMatrixFunctionRegistration);

                NekMatrix<DataType, space> operator+=(const NekMatrix<DataType, space>& rhs)
                {
                    ASSERTL0(rows() == rhs.rows() && columns() == rhs.columns(), "Matrix dimensions must agree in operator+");
                    for(unsigned int i = 0; i < rows(); ++i)
                    {
                        for(unsigned int j = 0; j < columns(); ++j)
                        {
                            (*this)(i,j) += rhs(i,j);
                        }
                    }

                    return *this;
                }

            private:
                // The implementation classes prevent the need for enumerated types
                // and large nested switch statements.
                class NekMatrixImpl
                {
                    public:
                        NekMatrixImpl(unsigned int rows, unsigned int columns) :
                            m_rows(rows),
                            m_columns(columns)
                        {
                        }

                        NekMatrixImpl(const NekMatrixImpl& rhs) :
                            m_rows(rhs.m_rows),
                            m_columns(rhs.m_columns)
                        {
                        }

                        virtual ~NekMatrixImpl()
                        {
                        };

                        NekMatrixImpl& operator=(const NekMatrixImpl& rhs)
                        {
                            m_rows = rhs.m_rows;
                            m_columns = rhs.m_columns;
                            return *this;
                        }

                        virtual void Initialize() = 0;
                        virtual void Initialize(const DataType* const p) = 0;
                        virtual void Initialize(typename boost::call_traits<DataType>::param_type d) = 0;

                        unsigned int rows() const { return m_rows; }
                        unsigned int columns() const { return m_columns; }

                        virtual typename boost::call_traits<DataType>::reference operator()
                                (unsigned int rowNumber, unsigned int colNumber) = 0;


                        virtual typename boost::call_traits<DataType>::const_reference operator()
                                (unsigned int rowNumber, unsigned int colNumber) const = 0;

                    protected:


                    private:
                        unsigned int m_rows;
                        unsigned int m_columns;
                };

                // Factory to generate implementations based upon the MatrixForm.
                typedef Loki::SingletonHolder< Loki::Factory<NekMatrixImpl, NekMatrixForm,
                    Loki::TL::MakeTypelist<unsigned int, unsigned int>::Result > >
                    NekMatrixImplFactory;

                class NekFullMatrixImpl : public NekMatrixImpl
                {
                    public:
                        NekFullMatrixImpl(unsigned int rows, unsigned int columns) :
                            NekMatrixImpl(rows, columns),
                            m_data()
                        {
                        }

                        NekFullMatrixImpl(const NekFullMatrixImpl& rhs) :
                            NekMatrixImpl(rhs),
                            m_data()
                        {
                            CreateSharedArray();
                            std::copy(rhs.begin(), rhs.end(), begin());
                        }

                        NekFullMatrixImpl& operator=(const NekFullMatrixImpl& rhs)
                        {
                            if( this != &rhs )
                            {
                                NekMatrixImpl::operator=(rhs);
                                CreateSharedArray();
                                std::copy(rhs.begin(), rhs.end(), begin());
                            }

                            return *this;
                        }

                        virtual ~NekFullMatrixImpl() {}

                        virtual void Initialize()
                        {
                            CreateSharedArray();

                            // Don't do anything else, it is up to the user
                            // to initialize.
                        }

                        virtual void Initialize(const DataType* const p)
                        {
                            CreateSharedArray();
                            std::copy(p, p+(rows()*columns()), begin());
                        }

                        virtual void Initialize(typename boost::call_traits<DataType>::param_type d)
                        {
                            CreateSharedArray();
                            std::fill(begin(), end(), d);
                        }

                        DataType* begin()
                        {
                            return &m_data[0];
                        }

                        DataType* end()
                        {
                            return &m_data[rows()*columns()];
                        }

                        typename boost::call_traits<DataType>::reference operator()
                            (unsigned int rowNumber, unsigned int colNumber)
                        {
                            return m_data[rowNumber*columns() + colNumber];
                        }

                        typename boost::call_traits<DataType>::const_reference operator()
                                (unsigned int rowNumber, unsigned int colNumber) const
                        {
                            return m_data[rowNumber*columns() + colNumber];
                        }

                        static const MemoryManager::MemoryPoolEnabler MemoryPoolEnabled = MemoryManager::eEnabled;

                    protected:
                        void CreateSharedArray()
                        {
                            m_data = MemoryManager::AllocateSharedArray<DataType>(rows()*columns());
                            //m_data = new DataType[rows()*columns()];
                            //m_data = boost::shared_array<DataType>(new DataType[rows()*columns()]);
                        }

                    private:
                        NekFullMatrixImpl();
                        boost::shared_array<DataType> m_data;
                        //DataType* m_data;
                };


                static NekFullMatrixImpl* createFullMatrix(unsigned int rows, unsigned int cols)
                {
                    //NekFullMatrixImpl* result = MemoryManager::Allocate<NekFullMatrixImpl>(rows, cols);
                    NekFullMatrixImpl* result = new NekFullMatrixImpl(rows, cols);
                    return result;
                }

                static const bool fullRegistered;

                NekMatrixForm m_form;

                // Change back to the shared pointer after finishing the
                // NekFactory.
                //boost::shared_ptr<NekMatrixImpl> m_impl;
                NekMatrixImpl* m_impl;
        };

        template<typename DataType, unsigned int space>
        const bool NekMatrix<DataType, space>::fullRegistered =
                NekMatrix<DataType, space>::NekMatrixImplFactory::
                Instance().Register(eFull, NekMatrix<DataType, space>::createFullMatrix);


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

