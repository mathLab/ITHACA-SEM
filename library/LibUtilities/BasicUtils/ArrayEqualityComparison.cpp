///////////////////////////////////////////////////////////////////////////////
//
// File ArrayEqualityComparison.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    bool operator==(const Array<OneD, NekDouble>& lhs, const Array<OneD, NekDouble>& rhs) 
    {     
        return IsEqual(lhs,rhs);
    }


    bool IsEqual(const Array<OneD, const NekDouble>& lhs,
                 const Array<OneD, const NekDouble>& rhs,
                 NekDouble tol)
    {
        if( lhs.num_elements() != rhs.num_elements() )
        {
            return false;
        }
        
        if( lhs.data() == rhs.data() )
        {
            return true;
        }
        
        for(unsigned int i = 0; i < lhs.num_elements(); ++i)
        {
            if( fabs(lhs[i]-rhs[i]) > tol )
            {
                return false;
            }
        }
        
        return true;
    }


    bool operator==(const Array<TwoD, NekDouble>& lhs, const Array<TwoD, NekDouble>& rhs) 
    {
        return IsEqual(lhs,rhs);
    }

    bool IsEqual(const Array<TwoD, const NekDouble>& lhs,
                 const Array<TwoD, const NekDouble>& rhs,
                 NekDouble tol)
    {
        if( (lhs.GetRows() != rhs.GetRows()) || 
            (lhs.GetColumns() != rhs.GetColumns()) )
        {
            return false;
        }
        
        if( lhs.data() == rhs.data() )
        {
            return true;
        }

        for(unsigned int i = 0; i < lhs.GetRows(); ++i)
        {
            for(unsigned int j = 0; j < lhs.GetColumns(); ++j)
            {
                if( fabs(lhs[i][j]-rhs[i][j]) > tol )
                {
                    return false;
                }
            }
        }
        
        return true;
    }
}
