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

#include <map>
#include <sstream>
#include <memory>
#include <functional>

#ifdef NEKTAR_USE_THREAD_SAFETY
#include <boost/thread/shared_mutex.hpp>
#include <boost/thread/locks.hpp>
#endif

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
#ifdef NEKTAR_USE_THREAD_SAFETY
        typedef boost::unique_lock<boost::shared_mutex> WriteLock;
        typedef boost::shared_lock<boost::shared_mutex> ReadLock;
#endif

        template <typename KeyType>
        struct defOpLessCreator
        {
            bool operator()(const KeyType &lhs, const KeyType &rhs) const
            {
                return lhs < rhs;
            }
        };

        template <typename KeyType, typename ValueT, typename opLessCreator = defOpLessCreator<KeyType> >
        class NekManager
        {
            public:
                typedef std::shared_ptr<ValueT> ValueType;
                typedef std::function<ValueType (const KeyType& key)> CreateFuncType;
                typedef std::map<KeyType, ValueType> ValueContainer;
                typedef std::shared_ptr<ValueContainer> ValueContainerShPtr;
                typedef std::map<KeyType, CreateFuncType, opLessCreator> CreateFuncContainer;
                typedef std::map<std::string, std::shared_ptr<ValueContainer> > ValueContainerPool;
                typedef std::shared_ptr<bool> BoolSharedPtr;
                typedef std::map<std::string, BoolSharedPtr> FlagContainerPool;

                NekManager(std::string whichPool="") :
                    m_values(),
                    m_globalCreateFunc(),
                    m_keySpecificCreateFuncs()

                {
                    if (!whichPool.empty())
                    {
                        auto iter = m_ValueContainerPool.find(whichPool);
                        if (iter != m_ValueContainerPool.end())
                        {
                            m_values = iter->second;
                            m_managementEnabled = m_managementEnabledContainerPool[whichPool];
                        }
                        else
                        {
                            m_values = ValueContainerShPtr(new ValueContainer);
                            m_ValueContainerPool[whichPool] = m_values;
                            if (m_managementEnabledContainerPool.find(whichPool) == m_managementEnabledContainerPool.end())
                            {
                                m_managementEnabledContainerPool[whichPool] = BoolSharedPtr(new bool(true));
                            }
                            m_managementEnabled = m_managementEnabledContainerPool[whichPool];
                        }
                    }
                    else
                    {
                        m_values = ValueContainerShPtr(new ValueContainer);
                        m_managementEnabled = BoolSharedPtr(new bool(true));
                    }
                };


                explicit NekManager(CreateFuncType f, std::string whichPool="") :
                    m_values(),
                    m_globalCreateFunc(f),
                    m_keySpecificCreateFuncs()
                {
                    if (!whichPool.empty())
                    {
#ifdef NEKTAR_USE_THREAD_SAFETY
                        ReadLock v_rlock(m_mutex); // reading static members
#endif
                        auto iter = m_ValueContainerPool.find(whichPool);
                        if (iter != m_ValueContainerPool.end())
                        {
                            m_values = iter->second;
                            m_managementEnabled = m_managementEnabledContainerPool[whichPool];
                        }
                        else
                        {
#ifdef NEKTAR_USE_THREAD_SAFETY
                            v_rlock.unlock();
                            // now writing static members.  Apparently
                            // upgrade_lock has less desirable properties than
                            // just dropping read lock, grabbing write lock.
                            // write will block until all reads are done, but
                            // reads cannot be acquired if write lock is
                            // blocking.  In this context writes are supposed to
                            // be rare.
                            WriteLock v_wlock(m_mutex);
#endif
                            m_values = ValueContainerShPtr(new ValueContainer);
                            m_ValueContainerPool[whichPool] = m_values;
                            if (m_managementEnabledContainerPool.find(whichPool) == m_managementEnabledContainerPool.end())
                            {
                                m_managementEnabledContainerPool[whichPool] = BoolSharedPtr(new bool(true));
                            }
                            m_managementEnabled = m_managementEnabledContainerPool[whichPool];
                        }

                    }
                    else
                    {
                        m_values = ValueContainerShPtr(new ValueContainer);
                        m_managementEnabled = BoolSharedPtr(new bool(true));
                    }
                }

                ~NekManager()
                {
                }

                /// Register the given function and associate it with the key.
                /// The return value is just to facilitate calling statically.
                bool RegisterCreator(const KeyType& key,
                                     const CreateFuncType& createFunc)
                {
                    m_keySpecificCreateFuncs[key] = createFunc;

                    return true;
                }

                /// Register the Global Create Function.
                /// The return value is just to facilitate calling statically.
                bool RegisterGlobalCreator(const CreateFuncType& createFunc)
                {
                    m_globalCreateFunc = createFunc;

                    return true;
                }

                bool AlreadyCreated(const KeyType &key)
                {
                    bool value = false;
                    auto found = m_values->find(key);
                    if( found != m_values->end() )
                    {
                        value = true;
                    }

                    return value;
                }

                ValueType operator[](const KeyType &key)
                {
                    auto found = m_values->find(key);

                    if( found != m_values->end() )
                    {
                        return (*found).second;
                    }
                    else
                    {
                        // No object, create a new one.
                        CreateFuncType f = m_globalCreateFunc;
                        auto keyFound = m_keySpecificCreateFuncs.find(key);
                        if( keyFound != m_keySpecificCreateFuncs.end() )
                        {
                            f = (*keyFound).second;
                        }

                        if( f )
                        {
                            ValueType v = f(key);
                            if (*m_managementEnabled)
                            {
                                (*m_values)[key] = v;
                            }
                            return v;
                        }
                        else
                        {
                            std::stringstream ss;
                            ss << key;
                            std::string message = "No create func found for key " + ss.str();
                            NEKERROR(ErrorUtil::efatal, message.c_str());
                            static ValueType result;
                            return result;
                        }
                    }
                }

                void DeleteObject(const KeyType &key)
                {
                    auto found = m_values->find(key);

                    if( found != m_values->end() )
                    {
                        m_values->erase(found);
                    }
                }

                static void ClearManager(std::string whichPool = "")
                {
                    if (!whichPool.empty())
                    {
#ifdef NEKTAR_USE_THREAD_SAFETY
                        WriteLock v_wlock(m_mutex);
#endif
                        auto x = m_ValueContainerPool.find(whichPool);
                        ASSERTL1(x != m_ValueContainerPool.end(),
                                "Could not find pool " + whichPool);
                        x->second->clear();
                    }
                    else
                    {
#ifdef NEKTAR_USE_THREAD_SAFETY
                        WriteLock v_wlock(m_mutex);
#endif

                        for (auto &x : m_ValueContainerPool)
                        {
                            x.second->clear();
                        }
                    }
                }

                static bool PoolCreated(std::string whichPool)
                {
                    bool value = false;
                    auto x = m_ValueContainerPool.find(whichPool);
                    if (x != m_ValueContainerPool.end())
                    {
                        value = true;
                    }
                    return value;
                }

                static void EnableManagement(std::string whichPool = "")
                {
                    if (!whichPool.empty())
                    {
#ifdef NEKTAR_USE_THREAD_SAFETY
                        WriteLock v_wlock(m_mutex);
#endif

                        auto x = m_managementEnabledContainerPool.find(whichPool);
                        if (x != m_managementEnabledContainerPool.end())
                        {
                            (*x->second) = true;
                        }
                        else
                        {
                            m_managementEnabledContainerPool[whichPool] = BoolSharedPtr(new bool(true));
                        }
                    }
                }

                static void DisableManagement(std::string whichPool = "")
                {
                    if (!whichPool.empty())
                    {
#ifdef NEKTAR_USE_THREAD_SAFETY
                        WriteLock v_wlock(m_mutex);
#endif
                        auto x = m_managementEnabledContainerPool.find(whichPool);
                        if (x != m_managementEnabledContainerPool.end())
                        {
                            (*x->second) = false;
                        }
                        else
                        {
                            m_managementEnabledContainerPool[whichPool] = BoolSharedPtr(new bool(false));
                        }
                    }
                }

            private:
                NekManager<KeyType, ValueType, opLessCreator>& operator=(const NekManager<KeyType, ValueType, opLessCreator>& rhs);
                NekManager(const NekManager<KeyType, ValueType, opLessCreator>& rhs);

                ValueContainerShPtr m_values;
                BoolSharedPtr m_managementEnabled;
                static ValueContainerPool m_ValueContainerPool;
                static FlagContainerPool m_managementEnabledContainerPool;
                CreateFuncType m_globalCreateFunc;
                CreateFuncContainer m_keySpecificCreateFuncs;
#ifdef NEKTAR_USE_THREAD_SAFETY
                static boost::shared_mutex m_mutex;
#endif
        };
        template <typename KeyType, typename ValueT, typename opLessCreator> typename NekManager<KeyType, ValueT, opLessCreator>::ValueContainerPool NekManager<KeyType, ValueT, opLessCreator>::m_ValueContainerPool;
        template <typename KeyType, typename ValueT, typename opLessCreator> typename NekManager<KeyType, ValueT, opLessCreator>::FlagContainerPool NekManager<KeyType, ValueT, opLessCreator>::m_managementEnabledContainerPool;
#ifdef NEKTAR_USE_THREAD_SAFETY
        template <typename KeyType, typename ValueT, typename opLessCreator>
            typename boost::shared_mutex NekManager<KeyType, ValueT, opLessCreator>::m_mutex;
#endif
    }
}


#endif //NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_MANAGER_HPP
