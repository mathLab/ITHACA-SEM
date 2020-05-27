///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Operator top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Operator.h>
#include <Collections/Collection.h>

namespace Nektar {
namespace Collections {

/**
 *
 */
bool operator< (OperatorKey const &p1, OperatorKey const &p2)
{
    if (std::get<0>(p1) < std::get<0>(p2))
    {
        return true;
    }
    if (std::get<0>(p1) > std::get<0>(p2))
    {
        return false;
    }
    if (std::get<1>(p1) < std::get<1>(p2))
    {
        return true;
    }
    if (std::get<1>(p1) > std::get<1>(p2))
    {
        return false;
    }
    if (std::get<2>(p1) < std::get<2>(p2))
    {
        return true;
    }
    if (std::get<2>(p1) > std::get<2>(p2))
    {
        return false;
    }

    if (std::get<3>(p1) < std::get<3>(p2))
    {
        return true;
    }
    if (std::get<3>(p1) > std::get<3>(p2))
    {
        return false;
    }

    return false;
}


/**
 *
 */
std::ostream &operator<<(std::ostream &os, OperatorKey const &p)
{
    os << LibUtilities::ShapeTypeMap[std::get<0>(p)] << ", "
       << OperatorTypeMap           [std::get<1>(p)] << ", "
       << ImplementationTypeMap     [std::get<2>(p)] << ", "
       << (std::get<3>(p) ? "Nodal" : "Modal");
    return os;
}


/**
 *
 */
Operator::~Operator()
{
}


/**
 *
 */
OperatorFactory& GetOperatorFactory()
{
    static OperatorFactory instance;
    return instance;
}



// simple operator Map evaluation
OperatorImpMap SetFixedImpType(ImplementationType defaultType)
{
    OperatorImpMap opMap;

    for(int i = 0; i < SIZE_OperatorType; ++i)
    {
        opMap[(OperatorType)i] = defaultType;
    }

    return opMap;
}

}
}
