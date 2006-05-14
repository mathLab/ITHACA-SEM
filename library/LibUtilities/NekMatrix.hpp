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

#include <LibUtilities/NekMemoryManager.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        enum NekMatrixForm
        {
            eFull,
            eZero,
            eDiagonal,
            eSquareSymmetric,
            eSquareSymmetricPositiveDefinite,
            eSymmetricPositiveDefiniteBanded,
            eSquareGeneral,
            eSquareGeneralBanded
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
        template<typename DataType, unsigned int NumRows = 0, unsigned int NumColumns = 0, unsigned int space = 0>
        class NekMatrix
        {

            public:
                NekMatrix(NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(NumRows),
                    m_columns(NumColumns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
                {
                    m_impl->Initialize();
                }

                NekMatrix(unsigned int rows, unsigned int columns, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(rows),
                    m_columns(columns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
                {
                    m_impl->Initialize();
                }

                explicit NekMatrix(const DataType* const ptr, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(NumRows),
                    m_columns(NumColumns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
                {
                    m_impl->Initialize(ptr);
                }

                NekMatrix(const DataType* const ptr, unsigned int rows, unsigned int columns,
                          NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(rows),
                    m_columns(columns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
                {
                    m_impl->Initialize(ptr);
                }

                NekMatrix(typename boost::call_traits<DataType>::param_type d, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(NumRows),
                    m_columns(NumColumns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
                {
                    m_impl->Initialize(d);
                }

                NekMatrix(typename boost::call_traits<DataType>::param_type d,
                          unsigned int rows, unsigned int columns, NekMatrixForm form = eFull) :
                    m_form(form),
                    m_rows(rows),
                    m_columns(columns),
                    m_impl(NekMatrixImplFactory::Instance().CreateObject(form, m_rows, m_columns))
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
                    if( rowNumber >= rows() || colNumber >= columns() )
                    {
                        throw OutOfBoundsError();
                    }

                    return (*m_impl)(rowNumber, colNumber);
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

                        virtual ~NekMatrixImpl() = 0;

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

                    protected:
                        void CreateSharedArray()
                        {
                            m_data = MemoryManager::AllocateSharedArray<DataType>(rows()*columns());
                        }

                    private:
                        NekFullMatrixImpl();
                        boost::shared_array<DataType> m_data;
                };


                NekFullMatrixImpl* createFullMatrix()
                {
                    return MemoryManager::Allocate<NekFullMatrixImpl>();
                }

                static const bool fullRegistered;


                NekMatrixForm m_form;
                unsigned int m_rows;
                unsigned int m_columns;
                boost::shared_ptr<NekMatrixImpl> m_impl;
        };

        template<typename DataType, unsigned int width, unsigned int height, unsigned int space>
        const bool NekMatrix<DataType, width, height, space>::fullRegistered =
                NekMatrix<DataType, width, height, space>::NekMatrixImplFactory::
                Instance().Register(eFull, NekMatrix<DataType, width, height, space>::createFullMatrix);
    }
}

#endif //NEKTAR_LIB_UTILITIES_NEK_MATRIX_HPP

/**
    $Log: NekMatrix.hpp,v $
    Revision 1.1  2006/05/04 18:57:43  kirby
    *** empty log message ***

    Revision 1.1  2006/04/11 02:00:43  bnelson
    Initial Revision


**/

