///////////////////////////////////////////////////////////////////////////////
//
// File: Collection.cpp
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
// Description: Collection top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <Collections/Collection.h>
#include <sstream>

using namespace std;

namespace Nektar {
namespace Collections {

/**
 *
 */
Collection::Collection(
        vector<StdRegions::StdExpansionSharedPtr>    pCollExp,
        OperatorImpMap                              &impTypes):
    m_collExp(pCollExp), 
    m_impTypes(impTypes)
{
    // Initialise geometry data.
    m_geomData = MemoryManager<CoalescedGeomData>::AllocateSharedPtr();

#if 0 
    // Loop over all operator types.
    for (int i = 0; i < SIZE_OperatorType; ++i)
    {
        OperatorType opType = (OperatorType)i;
        ImplementationType impType;

        auto it = impTypes.find(opType);
        if (it != impTypes.end())
        {
            impType = it->second;
            OperatorKey opKey(pCollExp[0]->DetShapeType(), opType, impType,
                              pCollExp[0]->IsNodalNonTensorialExp());

            stringstream ss;
            ss << opKey;
            ASSERTL0(GetOperatorFactory().ModuleExists(opKey),
                 "Requested unknown operator "+ss.str());

            m_ops[opType] = GetOperatorFactory().CreateInstance(
                                                opKey, pCollExp, m_geomData);
        }
    }
#endif    
}

void Collection::Initialise(const OperatorType opType)
{
    if(!HasOperator(opType))
    {
        auto it = m_impTypes.find(opType);
        
        if (it != m_impTypes.end())
        {
            ImplementationType impType = it->second;
            OperatorKey opKey(m_collExp[0]->DetShapeType(), opType, impType,
                              m_collExp[0]->IsNodalNonTensorialExp());
            
            stringstream ss;
            ss << opKey;
            ASSERTL0(GetOperatorFactory().ModuleExists(opKey),
                     "Requested unknown operator "+ss.str());
            
            m_ops[opType] = GetOperatorFactory().CreateInstance(
                                          opKey, m_collExp, m_geomData);
        }
        else
        {
            NEKERROR(ErrorUtil::ewarning,
                     "Failed to determine implmentation to initialise "
                     "collection operator: " +
                     std::string(Collections::OperatorTypeMap[opType]));
        }
    }
}
}
}

