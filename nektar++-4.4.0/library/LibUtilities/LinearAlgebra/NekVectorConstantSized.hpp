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
//// Description: Generic N-Dimensional Vector.
////
/////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP
//#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP
//
//#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
//#include <LibUtilities/LinearAlgebra/NekVectorTypeTraits.hpp>
//#include <LibUtilities/LinearAlgebra/NekVectorCommon.hpp>
//#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
//
//#include <ExpressionTemplates/ExpressionTemplates.hpp>
//#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
//
//
//#include <algorithm>
//#include <boost/static_assert.hpp>
//
//
//namespace Nektar
//{
//
//
//    template<typename dataType, typename dim, typename space>
//    class NekVector<dataType, dim, space>
//    {
//        public:
//            typedef dataType DataType;
//
//            NekVector() :
//                m_impl()
//            {
//                std::fill_n(m_impl, dim::Value, DataType());
//            }
//            
//            explicit NekVector(typename boost::call_traits<DataType>::const_reference a) :
//                m_impl()
//            {
//                std::fill_n(m_impl, dim::Value, a);
//            }
//
//            /// \brief Constructs a vector with initial value a.
//            ///
//            /// \param size The size of the vector.  Since the size of the vector is specified in the template, this parameter is ignored.
//            /// \param a The value with which to initialize the vector.
//            ///
//            /// This constructor is provided for generic methods which don't know if the vector they 
//            /// are creating is constant sized or variable sized.  
//            NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a) :
//                m_impl()
//            {
//                ASSERTL1(size == dim::Value, std::string("Attempting to construct a constant sized vector of size ") +
//                                             boost::lexical_cast<std::string>(dim::Value) + 
//                                             std::string(" with a size ") +
//                                             boost::lexical_cast<std::string>(size));
//                std::fill_n(m_impl, dim::Value, a);
//            }
//            
//            explicit NekVector(const std::string& vectorValues) :
//                m_impl()
//            {
//                std::vector<DataType> values = FromString<DataType>(vectorValues);
//
//                ASSERTL0(values.size() == dim::Value, "Error converting string values to vector");
//
//                std::copy(values.begin(), values.end(), &m_impl[0]);
//            }
//
//            NekVector(typename boost::call_traits<DataType>::const_reference x,
//                      typename boost::call_traits<DataType>::const_reference y,
//                      typename boost::call_traits<DataType>::const_reference z) :
//                m_impl()
//            {
//                BOOST_STATIC_ASSERT(dim::Value == 3);
//                m_impl[0] = x;
//                m_impl[1] = y;
//                m_impl[2] = z;
//            }
//
//            NekVector(typename boost::call_traits<DataType>::const_reference x,
//                      typename boost::call_traits<DataType>::const_reference y,
//                      typename boost::call_traits<DataType>::const_reference z,
//                      typename boost::call_traits<DataType>::const_reference w) :
//                m_impl()
//            {
//                BOOST_STATIC_ASSERT(dim::Value == 4);
//                m_impl[0] = x;
//                m_impl[1] = y;
//                m_impl[2] = z;
//                m_impl[3] = w;
//            }
//
//            NekVector(const NekVector<DataType, dim, space>& rhs) :
//                m_impl()
//            {
//                std::copy(rhs.m_impl, rhs.m_impl+dim::Value, m_impl);
//            }
//
//           explicit NekVector(const DataType* const ptr) :
//                m_impl()
//            {
//                std::copy(ptr, ptr+dim::Value, &m_impl[0]);
//            }
//
//            ~NekVector()
//            {
//#ifdef _DEBUG
//                std::fill_n(m_impl, dim::Value, DataType());
//#endif //_DEBUG
//            }
//            
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//            template<typename L, typename Op, typename R>
//            NekVector(const expt::Node<L, Op, R>& rhs) :
//                m_impl()
//            {
//				#ifdef _DEBUG
//                boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
//                    MatrixSize<Node<L, Op, R>, typename Node<L, Op, R>::Indices, 0>::GetRequiredSize(rhs.GetData());
//				ASSERTL0(sizes.get<0>() == dim::Value, "Data sizes are not equal.");
//				#endif
//
//                /// TODO Make sure this works correctly with eWrapper
//                //BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, VariableSizedVector, space> > ));
//                expt::ExpressionEvaluator::Evaluate(rhs, *this);
//            }
//#endif
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//            template<typename L, typename Op, typename R>
//            NekVector<DataType, dim, space>& operator=(const expt::Node<L, Op, R>& rhs)
//            {
//				rhs.Evaluate(*this);
//				return *this;
//            }
//#endif
//            
//            /// \brief Returns the number of dimensions for the point.
//            unsigned int GetDimension() const
//            {
//                return dim::Value;
//            }
//
//            /// \brief Treating the vector as a column vector, how many rows it has.
//            unsigned int GetRows() const
//            {
//                return dim::Value;
//            }
//
//            const DataType* GetRawPtr() const
//            {
//                return &m_impl[0];
//            }
//
//            typedef const DataType* const_iterator;
//            const_iterator begin() const { return GetRawPtr(); }
//            const_iterator end() const { return GetRawPtr() + GetDimension(); }
//
//            typename boost::call_traits<DataType>::const_reference operator()(unsigned int i) const
//            {
//                ASSERTL1(( i >= 0) && (i < GetDimension()), "Invalid access to m_data via parenthesis operator");
//                return m_impl[i];
//            }
//
//            typename boost::call_traits<DataType>::const_reference operator[](unsigned int i) const
//            {
//                return m_impl[i];
//            }
//
//            typename boost::call_traits<DataType>::const_reference x() const
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 1);
//                return m_impl[0];
//            }
//
//            typename boost::call_traits<DataType>::const_reference y() const
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 2);
//                return m_impl[1];
//            }
//
//            typename boost::call_traits<DataType>::const_reference z() const
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 3);
//                return m_impl[2];
//            }
//            
//            typename boost::call_traits<DataType>::const_reference w() const
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 4);
//                return m_impl[3];
//            }
//
//            
////            typename boost::call_traits<DataType>::reference x() 
////            {
////                BOOST_STATIC_ASSERT(dim::Value >= 1);
////                return m_impl[0];
////            }
////
////            typename boost::call_traits<DataType>::reference y() 
////            {
////                BOOST_STATIC_ASSERT(dim::Value >= 2);
////                return m_impl[1];
////            }
////
////            typename boost::call_traits<DataType>::reference z() 
////            {
////                BOOST_STATIC_ASSERT(dim::Value >= 3);
////                return m_impl[2];
////            }
////
////
////            typename boost::call_traits<DataType>::reference w() 
////            {
////                BOOST_STATIC_ASSERT(dim::Value >= 4);
////                return m_impl[3];
////            }
//
//            // Unitary operators
//            NekVector<DataType, dim, space> operator-() const { return Negate(*this); }
//
//            DataType Magnitude() const { return Nektar::Magnitude(*this); }
//            
//            DataType Dot(const NekVector<DataType, dim, space>& rhs) const { return Nektar::Dot(*this, rhs); }
//            
//            NekVector<DataType, dim, space> Cross(const NekVector<DataType, dim, space>& rhs) const
//            {
//                return Nektar::Cross(*this, rhs);
//            }
//
//            std::string AsString() const { return Nektar::AsString(*this); }
//
//            // Norms
//            DataType L1Norm() const { return Nektar::L1Norm(*this); }
//            DataType L2Norm() const { return Nektar::L2Norm(*this); }
//            DataType InfinityNorm() const { return Nektar::InfinityNorm(*this); }
//            
//            PointerWrapper GetWrapperType() const { return eCopy; }
//            
//        protected:
//            DataType* GetImpl() { return &m_impl[0]; }
//            
//            NekVector<DataType, dim, space>& operator=(const NekVector<DataType, dim, space>& rhs)
//            {
//                std::copy(rhs.m_impl, rhs.m_impl+dim::Value, m_impl);
//                return *this;
//            }
//            
//            NekVector<DataType, dim, space>& operator=(const NekVector<DataType, VariableSizedVector, space>& rhs)
//            {
//                ASSERTL0(GetDimension() == rhs.GetDimension(), "Assignment to a constant sized vector must have the same number of elements.");
//                std::copy(rhs.GetRawPtr(), rhs.GetRawPtr()+dim::Value, m_impl);
//                return *this;
//            }
//            
//        private:
//            DataType m_impl[dim::Value];
//
//    };
//    
//    // \param DataType The type of data held by each element of the vector.
//    // \param dim The number of elements in the vector.  If set to 0, the vector
//    //            will have a variable number of elements.
//    // \param space The space of the vector.
//    template<typename DataType, typename dim, typename space>
//    class NekVector : public NekVector<DataType, dim, space>
//    {
//        public:
//            typedef NekVector<DataType, dim, space> BaseType;
//            
//            /// \brief Creates a constant sized vector with each element initialized to
//            ///        the default value for DataType.
//            NekVector() : BaseType() {}
//
//            /// \brief Creates a constant sized vector with each element set to a.
//            /// \param a The value to assign to each element of the vector.
//            explicit NekVector(typename boost::call_traits<DataType>::const_reference a) :
//                BaseType(a) {}
//                
//            NekVector(unsigned int size, typename boost::call_traits<DataType>::const_reference a) :
//                BaseType(size, a) {}
//
//            /// \brief Creates a vector from the elements in a delimited string.
//            /// \param vectorValues A string the the vector values.
//            /// 
//            /// 
//            explicit NekVector(const std::string& vectorValues) : BaseType(vectorValues) {}
//
//            NekVector(typename boost::call_traits<DataType>::const_reference x,
//                      typename boost::call_traits<DataType>::const_reference y,
//                      typename boost::call_traits<DataType>::const_reference z) :
//                    BaseType(x, y, z) {}
//
//            NekVector(typename boost::call_traits<DataType>::const_reference x,
//                      typename boost::call_traits<DataType>::const_reference y,
//                      typename boost::call_traits<DataType>::const_reference z,
//                      typename boost::call_traits<DataType>::const_reference w) :
//                    BaseType(x,y,z,w) {}
//
//            NekVector(const NekVector<DataType, dim, space>& rhs) :
//                BaseType(rhs) {}
//
//            explicit NekVector(const DataType* const ptr) :
//                BaseType(ptr) {}
//            
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//            template<typename L, typename Op, typename R>
//            NekVector(const expt::Node<L, Op, R>& rhs) :
//                BaseType()
//            {
//				#ifdef _DEBUG
//                boost::tuple<unsigned int, unsigned int, unsigned int> sizes = 
//                    MatrixSize<expt::Node<L, Op, R>, typename expt::Node<L, Op, R>::Indices, 0>::GetRequiredSize(rhs.GetData());
//				ASSERTL0(sizes.get<0>() == dim::Value, "Data sizes are not equal.");
//				#endif
//
//                /// TODO Make sure this works correctly with eWrapper
//                //BOOST_MPL_ASSERT(( boost::is_same<typename Expression<ExpressionPolicyType>::ResultType, NekVector<DataType, VariableSizedVector, space> > ));
//                expt::ExpressionEvaluator::Evaluate(rhs, *this);
//            }
//#endif
//
//#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES
//            template<typename L, typename Op, typename R>
//            NekVector<DataType, dim, space>& operator=(const expt::Node<L, Op, R>& rhs)
//            {
//				rhs.Evaluate(*this);
//				return *this;
//            }
//#endif
//            
//
//
//            NekVector<DataType, dim, space>& operator=(const NekVector<DataType, dim, space>& rhs)
//            {
//                BaseType::operator=(rhs);
//                return *this;
//            }
//
//            NekVector<DataType, dim, space>& operator=(const NekVector<DataType, VariableSizedVector, space>& rhs)
//            {
//                BaseType::operator=(rhs);
//                return *this;
//            }
//            
//			using BaseType::GetRawPtr;
//            DataType* GetRawPtr()
//            {
//                return this->GetImpl();
//            }
//            
//                        
//            typedef DataType* iterator;
//            
//            
//            iterator begin() { return GetRawPtr(); }
//            iterator end() { return GetRawPtr() + this->GetDimension(); }
//            
//                        
//            /// \brief Returns i^{th} element.
//            /// \param i The element to return.
//            /// \pre i < dim
//            /// \return A reference to the i^{th} element.
//            ///
//            /// Retrieves the i^{th} element.  Since it returns a reference you may
//            /// assign a new value (i.e., p(2) = 3.2;)
//            ///
//            /// This operator performs range checking.
//            typename boost::call_traits<DataType>::reference operator()(unsigned int i)
//            {
//                ASSERTL1((i >= 0) && (i < this->GetDimension()), "Invalid access to m_data via parenthesis operator");
//                return this->GetImpl()[i];
//            }
//
//            typename boost::call_traits<DataType>::reference operator[](unsigned int i)
//            {
//                return this->GetImpl()[i];
//            }
//
//            using BaseType::operator();
//            using BaseType::operator[];
//            
//            using BaseType::x;
//            typename boost::call_traits<DataType>::reference x() 
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 1);
//                return this->GetImpl()[0];
//            }
//
//            using BaseType::y;
//            typename boost::call_traits<DataType>::reference y() 
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 2);
//                return this->GetImpl()[1];
//            }
//
//            using BaseType::z;
//            typename boost::call_traits<DataType>::reference z() 
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 3);
//                return this->GetImpl()[2];
//            }
//
//            using BaseType::w;
//            typename boost::call_traits<DataType>::reference w() 
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 4);
//                return this->GetImpl()[3];
//            }
//
//            void SetX(typename boost::call_traits<DataType>::const_reference val)
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 1);
//                this->GetImpl()[0] = val;
//            }
//
//            void SetY(typename boost::call_traits<DataType>::const_reference val)
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 2);
//                this->GetImpl()[1] = val;
//            }
//
//            void SetZ(typename boost::call_traits<DataType>::const_reference val)
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 3);
//                this->GetImpl()[2] = val;
//            }
//
//            void SetW(typename boost::call_traits<DataType>::const_reference val)
//            {
//                BOOST_STATIC_ASSERT(dim::Value >= 4);
//                this->GetImpl()[3] = val;
//            }
//
//            /// Arithmetic Routines
//
//            NekVector<DataType, dim, space>& operator+=(const NekVector<DataType, dim, space>& rhs)
//            {
//                AddEqual(*this, rhs);
//                return *this;
//            }
//
//            NekVector<DataType, dim, space>& operator-=(const NekVector<DataType, dim, space>& rhs)
//            {
//                SubtractEqual(*this, rhs);
//                return *this;
//            }
//
//            NekVector<DataType, dim, space>& operator*=(typename boost::call_traits<DataType>::const_reference rhs)
//            {
//                MultiplyEqual(*this, rhs);
//                return *this;
//            }
//            
//            NekVector<DataType, dim, space>& operator/=(typename boost::call_traits<DataType>::const_reference rhs)
//            {
//                DivideEqual(*this, rhs);
//                return *this;
//            }
//
//            
//            
//            void Normalize() { return Nektar::Normalize(*this); }
//            
//        
//
//    };    
//
//    
//    
//}
//
//#endif // NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_NEK_CONSTANT_SIZED_VECTOR_HPP
