///////////////////////////////////////////////////////////////////////////////
//
// File: Accumulator.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ACCUMULATOR_HPP
#define NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ACCUMULATOR_HPP

#include <boost/call_traits.hpp>

namespace Nektar
{

    /// \brief Wraps the data passed through the expression.
    ///
    /// Instead of using the result type of an expression as the acumulator,
    /// the data is wrapped in this accumulator object then passed through 
    /// the expression.  The main advantage to this approach is that we can 
    /// track if the accumulator has been initialized yet or not.
    template<typename DataType>
    class Accumulator
    {
        public:
            explicit Accumulator(typename boost::call_traits<DataType>::reference data) :
                m_isInitialized(false),
                m_data(data)
            {
            }
            
            /// \brief Returns true if the accumulator contains valid data, false otherwise.
            ///
            /// When starting an expression evaluation, the accumulator starts out in an uninitialized
            /// state and this method will return false.  At some point the accumulator will obtain
            /// valid values, and after this point this method will return true.
            inline bool IsInitialized() const 
            {
                return m_isInitialized; 
            }
            
            /// \brief Returns the data held by the accumulator.
            ///
            /// Without a little more infrastructure, it is impossible to determine 
            /// if the data has been modified after this call is made.  Therefore, 
            /// pessimistically, we will assume that it will be modified so all 
            /// subsequent calls to IsInitialized will return true.
            typename boost::call_traits<DataType>::reference operator*() const
            {
                m_isInitialized = true;
                return m_data;
            }

        private:
            /// By definition an accumulator stores a particular value and updates 
            /// it as the expression is being evaluated.  So a copy constructor and 
            /// assignment operator don't really make much sense.
            Accumulator(const Accumulator<DataType>& rhs);
            Accumulator<DataType>& operator=(const Accumulator<DataType>& rhs);
            
            mutable bool m_isInitialized;
            typename boost::call_traits<DataType>::reference m_data;
    };
}

#endif //NEKTAR_LIB_UTILITIES_EXPRESSION_TEMPLATES_ACCUMULATOR_HPP


/**
    $Log: Accumulator.hpp,v $
    Revision 1.8  2007/10/13 03:33:05  bnelson
    *** empty log message ***

    Revision 1.7  2007/10/04 03:48:53  bnelson
    *** empty log message ***

    Revision 1.6  2007/10/03 02:58:03  bnelson
    *** empty log message ***

    Revision 1.5  2007/08/16 02:14:21  bnelson
    Moved expression templates to the Nektar namespace.

    Revision 1.4  2007/02/15 06:55:07  bnelson
    Added comments and removed some unused code.

    Revision 1.3  2007/01/30 23:37:15  bnelson
    *** empty log message ***

    Revision 1.2  2007/01/16 17:37:55  bnelson
    Wrapped everything with #ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

    Revision 1.1  2007/01/16 05:29:49  bnelson
    Major improvements for expression templates.

 **/

