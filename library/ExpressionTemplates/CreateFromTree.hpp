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

#ifndef EXPRESSION_TEMPLATES_CREATE_FROM_TREE_HPP
#define EXPRESSION_TEMPLATES_CREATE_FROM_TREE_HPP

#ifdef NEKTAR_USE_EXPRESSION_TEMPLATES

namespace expt
{
    /// \brief Creates an accumulator to be used to evaluate the given tree.
    ///
    /// During expression evaluation, there are times when the creation of a 
    /// temporary is unavoidable.  In these cases, we need to create an object 
    /// of the appropriate data type to be used during evaluation.
    ///
    /// In some cases, the accumulator needs information from the tree.  An example
    /// is a tree using matrices, in which case the accumulator needs to know the 
    /// expected size of the result to size memory appropriately.  When this occurs, 
    /// the user is responsible for creating a specialized version of CreateFromTree
    /// to create the accumulator appropriately.
    template<typename DataType, typename NodeType, typename Indices, unsigned int StartIndex>
    struct CreateFromTree
    {
        template<typename ArgVectorType>
        static DataType Apply(const ArgVectorType& tree)
        {
            return DataType();
        }
    };

}

#endif //NEKTAR_USE_EXPRESSION_TEMPLATES
#endif //EXPRESSION_TEMPLATES_CREATE_FROM_TREE_HPP