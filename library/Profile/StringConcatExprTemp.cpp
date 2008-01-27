///////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////

#include "StringConcatExprTemp.h"
#include <LibUtilities/ExpressionTemplates/ExpressionTemplates.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryExpression.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpression.hpp>

template<typename ExpressionType>
void DoAssign(std::string& result, const ExpressionType& exp)
{
    std::string r(exp.GetMetadata().GetSize());
    exp.Apply(result);
}

void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2)
{
    Assign(result, Nektar::MakeExpr(str1) + Nektar::MakeExpr(str2));
}
                                     
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3)
{
    Assign(result, Nektar::MakeExpr(str1) + Nektar::MakeExpr(str2) + Nektar::MakeExpr(str3));
}

                                    
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4)
{
    Assign(result, Nektar::MakeExpr(str1) + Nektar::MakeExpr(str2) + Nektar::MakeExpr(str3) + Nektar::MakeExpr(str4));
}
                                     
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5)
{
    Assign(result, Nektar::MakeExpr(str1) + Nektar::MakeExpr(str2) + Nektar::MakeExpr(str3) +
            Nektar::MakeExpr(str4) + Nektar::MakeExpr(str5));
}
                                                                          
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6)
{
    Assign(result, Nektar::MakeExpr(str1) + Nektar::MakeExpr(str2) + Nektar::MakeExpr(str3) + 
            Nektar::MakeExpr(str4) + Nektar::MakeExpr(str5) + Nektar::MakeExpr(str6));
}
                                                                                                               
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6,
                const std::string& str7)
{
    Assign(result, Nektar::MakeExpr(str1) + 
            Nektar::MakeExpr(str2) + Nektar::MakeExpr(str3) + Nektar::MakeExpr(str4) + 
            Nektar::MakeExpr(str5) + Nektar::MakeExpr(str6) + Nektar::MakeExpr(str7));
}

namespace Nektar
{
    void NekAdd(std::string& result, const std::string& lhs, const std::string& rhs)
    {
        result = lhs;
        result += rhs;
    }
    
    void NekAddEqual(std::string& result, const std::string& rhs)
    {
        result += rhs;
    }
    
    std::string NekAdd(const std::string& lhs, const std::string& rhs)
    {
        std::string result;
        NekAdd(result, lhs, rhs);
        return result;
    }
}
