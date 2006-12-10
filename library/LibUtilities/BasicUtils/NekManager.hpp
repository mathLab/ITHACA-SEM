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

#include <utility>
#include <vector>
#include <algorithm>
#include <map>

#include <functional> // for bind1st
#include <iostream>

#include <boost/mem_fn.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/call_traits.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

using namespace std;

template <typename keyT, typename valueT,
    typename valueContainer = std::map<keyT, valueT>,
    typename creatorContainer = std::map<keyT, boost::function<valueT (const keyT& key)> > >
class NekManager
{
    public:
        typedef boost::function<valueT (const keyT& key)> CreateFunc;

        NekManager() :
            m_Values(), 
            m_CreateFuncs()
        {
        };
        
        ~NekManager() {}
	
        void RegisterCreator(typename boost::call_traits<keyT>::const_reference key, 
                             const CreateFunc& createFunc)
        {
            m_CreateFuncs[key] = createFunc;
        }

        void AddValue(typename boost::call_traits<keyT>::const_reference key, 
                      typename boost::call_traits<valueT>::const_reference value)
        {
            typename valueContainer::iterator found = m_Values.find(key);
            if( found == m_Values.end() )
            {
                m_Values[key] = value;
            }
            else
            {
                std::string errorMessage = "Adding an object for key " +
                        boost::lexical_cast<std::string>(key) + 
                        " to a NekManager, but a value for this key already exists.";

                throw std::runtime_error(errorMessage);
            }
        }

        typename boost::call_traits<valueT>::reference 
        GetValue(typename boost::call_traits<keyT>::const_reference key)
        {
            typename valueContainer::iterator found = m_Values.find(key);
            if( found == m_Values.end() )
            {
                // No object, create a new one.
                typename creatorContainer::iterator keyFound = m_CreateFuncs.find(key);
                if( keyFound != m_CreateFuncs.end() )
                {
                    valueT newObj = (*keyFound).second(key);
                    m_Values[key] = newObj;
                    return m_Values[key];
                }
                else
                {
                    std::string errorMessage = "Attempting to obtain an object for key " +
                            boost::lexical_cast<std::string>(key) + 
                            " in a NekManager but no create function has been registered";
                    throw std::runtime_error(errorMessage);
                }
            }
            else
            {
                return (*found).second;
            }
        }

    private:
        NekManager(const NekManager<keyT, valueT, valueContainer, creatorContainer>& rhs);
        NekManager<keyT, valueT, valueContainer, creatorContainer>& operator=(const NekManager<keyT, valueT, valueContainer, creatorContainer>& rhs);

        valueContainer m_Values;
        creatorContainer m_CreateFuncs;
};

#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_MANAGER_HPP
