///////////////////////////////////////////////////////////////////////////////
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

#ifndef NEKTAR_PROFILE_STRING_CONCAT_EXPRE_TEMP_H
#define NEKTAR_PROFILE_STRING_CONCAT_EXPRE_TEMP_H

#define NEKTAR_USE_EXPRESSION_TEMPLATES
#include <LibUtilities/ExpressionTemplates/Expression.hpp>
#include <LibUtilities/ExpressionTemplates/ExpressionMetadata.hpp>
#include <LibUtilities/ExpressionTemplates/ConstantExpressionTraits.hpp>
#include <LibUtilities/ExpressionTemplates/BinaryOperators.hpp>

#include <string>

namespace Nektar
{
    class StringMetadata
    {
        public:
            explicit StringMetadata(const std::string& value) : m_size(value.size()) {}
            StringMetadata(const StringMetadata& rhs) : m_size(rhs.m_size) {}
            StringMetadata(const StringMetadata& lhs, const StringMetadata& rhs) : m_size(lhs.m_size + rhs.m_size) {}
            StringMetadata() : m_size(0) {}
            StringMetadata& operator=(const StringMetadata& rhs)
            {
                m_size = rhs.m_size; 
                return *this;
            }
            
            unsigned int GetSize() const { return m_size; }
        private:
            unsigned int m_size;
    };

    template<>
    class ConstantExpressionTraits<std::string>
    {
        public:
            typedef std::string result_type;
            typedef StringMetadata MetadataType;
    };
    
    template<>
    class BinaryExpressionMetadataTraits<std::string, std::string, Nektar::AddOp>
    {
        public:
            typedef StringMetadata MetadataType;
    };
    
    void NekAdd(std::string& result, const std::string& lhs, const std::string& rhs);
    
    void NekAddEqual(std::string& result, const std::string& rhs);
    
    std::string NekAdd(const std::string& lhs, const std::string& rhs);
}

void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2);
                                     
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3);
                                     
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4);
                                     
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5);
                                                                          
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6);
                                                                                                               
void AddStringsExprTemp(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6,
                const std::string& str7);
#endif //NEKTAR_PROFILE_STRING_CONCAT_EXPRE_TEMP_H

