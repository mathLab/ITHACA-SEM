///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP
#define NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP

#include <boost/utility/enable_if.hpp>
#include "CommutativeTransform.hpp"
#include <boost/typeof/typeof.hpp>
#include <boost/mpl/and.hpp>
#include "Operators.hpp"



namespace Nektar
{

    template<typename NodeType>
    struct AssociativeTransform;

    template<typename T1, typename RootOp,
             typename R1, typename ROp, typename R2>
    struct AssociativeTransform<Node<T1, RootOp, Node<R1, ROp, R2> > >
    {
        typedef Node<Node<T1, RootOp, R1>, ROp, R2> TransformedNodeType;
    };

    

    
}



#endif //NEKTAR_EXPRESSION_TEMPLATES_ASSOCIATIVE_TRANSFORM_HPP

