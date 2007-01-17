///////////////////////////////////////////////////////////////////////////////
//
// File: NekManager.hpp
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

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_MANAGER_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_MANAGER_HPP

#include <algorithm>
#include <map>

#include <boost/function.hpp>
#include <boost/call_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/concept_check.hpp>

#include <LibUtilities/Foundations/Points.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

using namespace std;

namespace Nektar
{
    namespace LibUtilities
    {
        template <typename KeyType, typename ValueT>
        class NekManager
        {
            public:
                BOOST_CLASS_REQUIRE(KeyType, boost, LessThanComparableConcept);

                typedef boost::shared_ptr<ValueT> ValueType;
                typedef boost::function<ValueType (const KeyType& key)> CreateFuncType;
                typedef std::map<KeyType, ValueType> ValueContainer;
                typedef std::map<KeyType, CreateFuncType, opLess> CreateFuncContainer;

                NekManager() :
                    m_values(), 
                    m_globalCreateFunc(),
                    m_keySpecificCreateFuncs()
                {
                };

                explicit NekManager(CreateFuncType f) :
                    m_values(),
                    m_globalCreateFunc(f),
                    m_keySpecificCreateFuncs()
                {
                }
                
                ~NekManager() {}
        	
                void RegisterCreator(typename boost::call_traits<KeyType>::const_reference key, 
                                     const CreateFuncType& createFunc)
                {
                    m_keySpecificCreateFuncs[key] = createFunc;
                }

                ValueType operator[](typename boost::call_traits<KeyType>::const_reference key)
                {
                    typename ValueContainer::iterator found = m_values.find(key);
                    if( found != m_values.end() )
                    {
                        return (*found).second;
                    }
                    else
                    {
                        static ValueType result;
                    
                        // No object, create a new one.
                        CreateFuncType f = m_globalCreateFunc;
                        typename CreateFuncContainer::iterator keyFound = m_keySpecificCreateFuncs.find(key);
                        if( keyFound != m_keySpecificCreateFuncs.end() )
                        {
                            f = (*keyFound).second;
                        }

                        if( f )
                        {
                            result = f(key);
                            m_values[key] = result;
                            return m_values[key];
                        }
                        else
                        {
                            std::string keyAsString = boost::lexical_cast<std::string>(key);
                            std::string message = std::string("No create func found for key ") + keyAsString;
                            NEKERROR(ErrorUtil::efatal, message.c_str());
                        }
                        
                        return result;
                    }
                }
                
            private:
                NekManager(const NekManager<KeyType, ValueType>& rhs);
                NekManager<KeyType, ValueType>& operator=(const NekManager<KeyType, ValueType>& rhs);

                ValueContainer m_values;
                CreateFuncType m_globalCreateFunc;
                CreateFuncContainer m_keySpecificCreateFuncs;
        };
    }
}

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_MANAGER_HPP
