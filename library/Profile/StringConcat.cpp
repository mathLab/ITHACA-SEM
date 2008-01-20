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

#include "StringConcat.h"

void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2)
{
    std::string r = str1 + str2;
}
                                     
void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3)
{
    result = str1 + str2 + str3;
}

                                    
void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4)
{
    result = str1 + str2 + str3 + str4;
}
                                     
void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5)
{
    result = str1 + str2 + str3 + str4 + str5;
}
                                                                          
void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6)
{
    result = str1 + str2 + str3 + str4 + str5 + str6 ;
}
                                                                                                               
void AddStrings(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6,
                const std::string& str7)
{
    result = str1 + str2 + str3 + str4 + str5 + str6 + str7;
}


void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2)
{
    std::string r = str1;
    r += str2;
}
                                     
void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3)
{
    result = str1;
    result += str2;
    result += str3;
}

                                    
void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4)
{
    result = str1;
    result += str2;
    result += str3;
    result += str4;
}
                                     
void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5)
{
    result = str1;
    result += str2;
    result += str3;
    result += str4;
    result += str5;
}
                                                                          
void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6)
{
    result = str1;
    result += str2;
    result += str3;
    result += str4;
    result += str5;
    result += str6;
}
                                                                                                               
void AddStringsAccum(std::string& result, const std::string& str1,
                const std::string& str2,
                const std::string& str3,
                const std::string& str4,
                const std::string& str5,
                const std::string& str6,
                const std::string& str7)
{
    result = str1;
    result += str2;
    result += str3;
    result += str4;
    result += str5;
    result += str6;
    result += str7;
}

