///////////////////////////////////////////////////////////////////////////////
//
// File: NekFactory.hpp
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
// Description: Factory pattern class for Nektar
//
///////////////////////////////////////////////////////////////////////////////

// Primary definition and generator for specialised object factories.
#ifndef BOOST_PP_IS_ITERATING

    #ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_FACTORY_HPP
    #define NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_FACTORY_HPP

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/arithmetic/sub.hpp>
    #include <boost/preprocessor/punctuation/comma_if.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>

    #include <boost/shared_ptr.hpp>

    #include <iostream>
    #include <map>
    #include <string>

    #include <LibUtilities/BasicUtils/ErrorUtil.hpp>

    #ifndef MAX_PARAM
    #define MAX_PARAM 5  // default maximum number of parameters to support
    #endif

namespace Nektar
{
    namespace LibUtilities
    {
        // For unused template parameters.
        struct none {};

        // Generate parameter typenames with default type of 'none'
        #define FACTORY_print(z, n, data) BOOST_PP_CAT(data, n) = none

        /**
         * @class NekFactory
         * \brief Provides a generic Factory class.
         *
         * Implements a generic object factory. Class-types which use a
         * potentially arbitrary number of parameters may be used with
         * specialised forms of the NekFactory. An upper
         * limit on the number of parameters is imposed by the MAX_PARAM
         * preprocessor definition in the NekFactory.hpp file. The
         * specialisations are generated at compile type using Boost
         * preprocessor by through repeated inclusion of the NekFactory.hpp
         * file.
         *
         * To allow a class to be instantiated by the factory, the
         * following are required in each class definition (in the case of
         * a single parameter):
         * \code
         *   static [baseclass]* create([paramtype1] &P) {
         *      return new [derivedclass](P);
         *   }
         *   static std::string className;
         * \endcode
         * and outside the class definition in the implementation:
         * \code
         *   std::string [derivedclass]::className
         *      = Factory<std::string,[baseclass],[paramtype1]>::
         *          RegisterCreatorFunction("[derivedclass]",
         *                              [derivedclass]::create,"Description");
         * \endcode
         * The assignment of the static variable className is done through the
         * call to RegisterCreatorFunction, which registers the class with the
         * factory prior to the start of the main() routine.
         *
         * To create an instance of a derived class, for instance:
         * \code
         *   [baseclass]* var_name =
         *      Factory<std::string,[baseclass],[paramtype1]>
         *              ::CreateInstance("[derivedclass]",Param1);
         * \endcode
         */
        template <typename tKey,	       	// reference tag (e.g. string, int)
                  typename tBase,	       	// base class
                  BOOST_PP_ENUM(MAX_PARAM, FACTORY_print, typename tParam) >
        class NekFactory
        {
            public:
                /// Description datatype
                typedef std::string tDescription;
                /// Comparison predicator of key
                typedef std::less<tKey> tPredicator;
                /// Shared pointer to an object of baseclass type.
                typedef boost::shared_ptr<tBase> tBaseSharedPtr;
                /// CreatorFunction type which takes parameter and returns base
                /// class shared pointer.
                typedef tBaseSharedPtr (*CreatorFunction) (BOOST_PP_ENUM_PARAMS(MAX_PARAM, tParam));

                /// Define a struct to hold the information about a module.
                struct ModuleEntry
                {
                    ModuleEntry(CreatorFunction pFunc, const tDescription pDesc)
                        : m_func(pFunc),
                          m_desc(pDesc)
                    {
                    }

                    /// Function used to create instance of class.
                    CreatorFunction m_func;
                    /// Description of class for use in listing available classes.
                    tDescription m_desc;
                };

                /// Factory map between key and module data.
                typedef std::map<tKey, ModuleEntry, tPredicator> TMapFactory;
                /// Iterator for factory map
                typedef typename TMapFactory::iterator TMapFactoryIterator;


            public:
                NekFactory() {}

                /**
                 * @brief Create an instance of the class referred to by \c idKey.
                 *
                 * Searches the factory's map for the given key and returns a shared
                 * base class pointer to a new instance of the associated class.
                 * @param   idKey           Key of class to create.
                 * @param   x               Parameter to pass to class constructor.
                 * @returns                 Base class pointer to new instance.
                 */
                tBaseSharedPtr CreateInstance(tKey idKey BOOST_PP_COMMA_IF(MAX_PARAM)
                            BOOST_PP_ENUM_BINARY_PARAMS(MAX_PARAM, tParam, x))
                {
                    // Now try and find the key in the map.
                    TMapFactoryIterator it = getMapFactory()->find(idKey);

                    // If successful, check the CreatorFunction is defined and
                    // create a new instance of the class.
                    if (it != getMapFactory()->end())
                    {
                        if (it->second.m_func)
                        {
                            try
                            {
                                return it->second.m_func(BOOST_PP_ENUM_PARAMS(MAX_PARAM, x));
                            }
                            catch (const std::string& s)
                            {
                                std::stringstream errstr;
                                errstr << "Unable to create module: " << idKey << "\n";
                                errstr << s;
                                ASSERTL0(false, errstr.str());
                            }
                        }
                    }

                    // If we get this far, the key doesn't exist, so throw an error.
                    std::stringstream errstr;
                    errstr << "No such module: " << idKey << std::endl;
                    PrintAvailableClasses(errstr);
                    ASSERTL0(false, errstr.str());
                    return tBaseSharedPtr();
                }


                /**
                 * @brief Register a class with the factory.
                 *
                 * This function is called by each class in a static context (prior
                 * to the execution of main()) and creates an entry for the class
                 * in the factory's map.
                 * @param   idKey           Key used to reference the class.
                 * @param   classCreator    Function to call to create an instance
                 *                          of this class.
                 * @param   pDesc           Optional description of class.
                 * @returns                 The given key \c idKey.
                 */
                tKey RegisterCreatorFunction(tKey idKey, CreatorFunction classCreator,
                                             tDescription pDesc = "") 
                {
                    ModuleEntry e(classCreator, pDesc);
                    getMapFactory()->insert(std::pair<tKey,ModuleEntry>(idKey, e));
                    return idKey;
                }


                /**
                 * @brief Checks if a particular module is available.
                 */
                bool ModuleExists(tKey idKey)
                {
                    // Now try and find the key in the map.
                    TMapFactoryIterator it = getMapFactory()->find(idKey);

                    if (it != getMapFactory()->end())
                    {
                        return true;
                    }
                    return false;
                }


                /**
                 * @brief Prints the available classes to stdout.
                 */
                void PrintAvailableClasses(std::ostream& pOut = std::cout)
                {
                    pOut << std::endl << "Available classes: " << std::endl;
                    TMapFactoryIterator it;
                    for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
                    {
                        pOut << std::endl << "Available classes: " << std::endl;
                        TMapFactoryIterator it;
                        for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
                        {
                            pOut << "  " << it->first;
                            if (it->second.m_desc != "")
                            {
                                pOut << ":" << std::endl << "    "
                                          << it->second.m_desc << std::endl;
                            }
                            else
                            {
                                pOut << std::endl;
                            }
                        }
                    }
                }


                /**
                 * @brief Retrieves a key, given a description
                 */
                tKey GetKey(tDescription pDesc)
                {
                    TMapFactoryIterator it;
                    for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
                    {
                        if (it->second.m_desc == pDesc)
                        {
                            return it->first;
                        }
                    }
                    std::string errstr = "Module '"
                            + boost::lexical_cast<std::string>(pDesc)
                            + "' is not known.";
                    ASSERTL0(false, errstr);
                }


                /**
                 * @brief Returns the description of a class
                 */
                std::string GetClassDescription(tKey idKey)
                {
                    // Now try and find the key in the map.
                    TMapFactoryIterator it = getMapFactory()->find(idKey);

                    std::stringstream errstr;
                    errstr << "No such module: " << idKey << std::endl;
                    ASSERTL0 (it != getMapFactory()->end(), errstr.str());
                    return it->second.m_desc;
                }

            protected:
                /**
                 * @brief Ensure the factory's map is created.
                 * @returns                 The factory's map.
                 */
                TMapFactory* getMapFactory() 
                {
                    return &mMapFactory;
                }

            private:
                NekFactory(const NekFactory& rhs);
                NekFactory& operator=(const NekFactory& rhs);

                TMapFactory mMapFactory;

        };

        #undef FACTORY_print

        #define BOOST_PP_ITERATION_LIMITS (0, MAX_PARAM-1)
        #define BOOST_PP_FILENAME_1 "LibUtilities/BasicUtils/NekFactory.hpp"
        #include BOOST_PP_ITERATE()

    }
}

    #endif // end NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_FACTORY_HPP
    #undef MAX_PARAM

// Specialisations for the different numbers of parameters.
#else
    // Define the number of parameters
    #define n BOOST_PP_ITERATION()
    // Define macro for printing the non-required template parameters
    #define FACTORY_print(z, n, data) data

    template < typename tKey,
               typename tBase BOOST_PP_COMMA_IF(n)
               BOOST_PP_ENUM_PARAMS(n, typename tParam) >
    class NekFactory< tKey, tBase,
                BOOST_PP_ENUM_PARAMS(n, tParam)BOOST_PP_COMMA_IF(n)
                BOOST_PP_ENUM(BOOST_PP_SUB(MAX_PARAM,n), FACTORY_print, none) >
    {
    public:
        typedef std::string tDescription;
        typedef std::less<tKey> tPredicator;
        typedef boost::shared_ptr<tBase> tBaseSharedPtr;
        typedef tBaseSharedPtr (*CreatorFunction) (BOOST_PP_ENUM_PARAMS(n, tParam));

        struct ModuleEntry
        {
            ModuleEntry(CreatorFunction pFunc, const tDescription pDesc)
                : m_func(pFunc),
                  m_desc(pDesc)
            {
            }
            CreatorFunction m_func;
            tDescription m_desc;
        };
        typedef std::map<tKey, ModuleEntry, tPredicator> TMapFactory;
        typedef typename TMapFactory::iterator TMapFactoryIterator;

        NekFactory() {}
        
        tBaseSharedPtr CreateInstance(tKey idKey BOOST_PP_COMMA_IF(n)
                BOOST_PP_ENUM_BINARY_PARAMS(n, tParam, x))
        {
            TMapFactoryIterator it = getMapFactory()->find(idKey);
            if (it != getMapFactory()->end())
            {
                if (it->second.m_func)
                {
                    try
                    {
                        return it->second.m_func(BOOST_PP_ENUM_PARAMS(n, x));
                    }
                    catch (const std::string& s)
                    {
                        std::stringstream errstr;
                        errstr << "Unable to create module: " << idKey << "\n";
                        errstr << s;
                        ASSERTL0(false, errstr.str());
                    }
                }
            }
            std::stringstream errstr;
            errstr << "No such module: " << idKey << std::endl;
            PrintAvailableClasses(errstr);
            ASSERTL0(false, errstr.str());
            return tBaseSharedPtr();
        }

        tKey RegisterCreatorFunction(tKey idKey,
                                            CreatorFunction classCreator,
                                            tDescription pDesc = "") {
            ModuleEntry e(classCreator, pDesc);
            getMapFactory()->insert(std::pair<tKey,ModuleEntry>(idKey, e));
            return idKey;
        }

        bool ModuleExists(tKey idKey)
        {
            // Now try and find the key in the map.
            TMapFactoryIterator it = getMapFactory()->find(idKey);

            if (it != getMapFactory()->end())
            {
                return true;
            }
            return false;
        }

        void PrintAvailableClasses(std::ostream& pOut = std::cout)
        {
            pOut << std::endl << "Available classes: " << std::endl;
            TMapFactoryIterator it;
            for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
            {
                pOut << "  " << it->first;
                if (it->second.m_desc != "")
                {
                    pOut << ":" << std::endl << "    "
                              << it->second.m_desc << std::endl;
                }
                else
                {
                    pOut << std::endl;
                }
            }
        }

        tKey GetKey(tDescription pDesc)
        {
            TMapFactoryIterator it;
            for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
            {
                if (it->second.m_desc == pDesc)
                {
                    return it->first;
                }
            }
            std::string errstr = "Module '"
                    + boost::lexical_cast<std::string>(pDesc)
                    + "' is not known.";
            ASSERTL0(false, errstr);
        }

        std::string GetClassDescription(tKey idKey)
        {
            // Now try and find the key in the map.
            TMapFactoryIterator it = getMapFactory()->find(idKey);

            std::stringstream errstr;
            errstr << "No such module: " << idKey << std::endl;
            ASSERTL0 (it != getMapFactory()->end(), errstr.str());
            return it->second.m_desc;
        }

    protected:
        TMapFactory * getMapFactory() {
            return &mMapFactory;
        }

    private:
        NekFactory(const NekFactory& rhs);
        NekFactory& operator=(const NekFactory& rhs);

        TMapFactory mMapFactory;
    };
    #undef n
    #undef FACTORY_print

#endif

